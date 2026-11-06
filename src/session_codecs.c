/*
 * The codecs of a .cellucid-session file (session format 4) that reading one
 * needs (R/session-file.R), exactly as the web app writes them
 * (app/session/history-file.js, codecs/cell-set.js) and cellucid-python
 * reads them:
 *
 * - inflating one gzip (or zlib) stream to exactly its declared size, with no
 *   bytes after it (a truncated, longer, shorter or concatenated stream is an
 *   error, never a hang: base R's memDecompress() does not end on two
 *   concatenated gzip members);
 * - an array chunk stored as a cell set (a strictly increasing Uint32Array:
 *   sorted deltas, a bitset or runs).
 *
 * Every unsigned LEB128 varint must be canonical (no redundant byte), and
 * every count is checked against the bytes that remain before anything is
 * allocated from it.
 */

#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>

#define MAX_SAFE_INTEGER 9007199254740991.0

/* --- Inflate ---------------------------------------------------------------- */

/*
 * One gzip (`gzip` TRUE) or zlib stream, inflated. `expected` is the exact
 * size of the output, or negative (or NA) when unknown: the output then grows.
 *
 * Every size is held in size_t and bounded by R_XLEN_T_MAX, the longest raw
 * vector this R can return: on 32-bit R (the i386 build of Windows R 4.1)
 * R_xlen_t is int and size_t 32 bits, so neither a bound nor a product may
 * be formed in R_xlen_t. zlib counts each buffer in uInt, so input and output
 * longer than UINT_MAX are passed in pieces: nothing is refused for its size
 * but an output R cannot hold.
 */
SEXP cellucid_inflate(SEXP bytes, SEXP gzip, SEXP expected) {
  if (TYPEOF(bytes) != RAWSXP) error("inflate: bytes must be a raw vector");
  double wanted = asReal(expected);
  int is_gzip = asLogical(gzip) == TRUE;
  const char *kind = is_gzip ? "gzip" : "zlib";
  const size_t most = (size_t)R_XLEN_T_MAX;
  int known = wanted >= 0;
  if (known && wanted != floor(wanted)) error("inflate: the declared size must be a whole number of bytes");
  if (known && wanted > (double)R_XLEN_T_MAX) {
    error("inflate: %.0f bytes are more than one R raw vector holds", wanted);
  }
  const unsigned char *input = RAW(bytes);
  const size_t input_length = (size_t)XLENGTH(bytes);
  size_t capacity;
  if (known) {
    /* One byte more than declared: a stream holding more is told apart. */
    capacity = (size_t)wanted + 1;
  } else {
    double guess = (double)input_length * 4 + 1024;
    capacity = guess > (double)most ? most + 1 : (size_t)guess;
  }
  unsigned char *out = (unsigned char *)R_alloc(capacity, 1);

  z_stream stream;
  memset(&stream, 0, sizeof(stream));
  if (inflateInit2(&stream, is_gzip ? 16 + MAX_WBITS : MAX_WBITS) != Z_OK) {
    error("inflate: zlib could not start");
  }
  size_t fed = 0;
  size_t produced = 0;
  int status = Z_OK;
  for (;;) {
    if (produced == capacity) {
      /* Declared: one byte past it. Unknown: one byte past what R holds. */
      if (known || capacity > most) break;
      size_t grown = capacity > most / 2 ? most + 1 : capacity * 2;
      unsigned char *larger = (unsigned char *)R_alloc(grown, 1);
      memcpy(larger, out, produced);
      out = larger;
      capacity = grown;
    }
    if (stream.avail_in == 0 && fed < input_length) {
      size_t piece = input_length - fed;
      if (piece > UINT_MAX) piece = UINT_MAX;
      stream.next_in = (Bytef *)(input + fed);
      stream.avail_in = (uInt)piece;
      fed += piece;
    }
    size_t room = capacity - produced;
    if (room > UINT_MAX) room = UINT_MAX;
    stream.next_out = out + produced;
    stream.avail_out = (uInt)room;
    status = inflate(&stream, Z_NO_FLUSH);
    produced += room - stream.avail_out;
    if (status == Z_STREAM_END) break;
    /* No input is left to give: the stream ended mid-way. */
    if (status == Z_BUF_ERROR && stream.avail_in == 0 && fed == input_length) break;
    if (status != Z_OK && status != Z_BUF_ERROR) {
      const char *reason = stream.msg != NULL ? stream.msg : "corrupt data";
      inflateEnd(&stream);
      error("Invalid %s stream (%s)", kind, reason);
    }
  }
  int ended = status == Z_STREAM_END;
  size_t left = (size_t)stream.avail_in + (input_length - fed);
  inflateEnd(&stream);
  if (known && produced > (size_t)wanted) {
    error("Session chunk uncompressedBytes is smaller than the gzip payload");
  }
  if (!known && produced > most) error("inflate: the %s stream holds more bytes than one R raw vector holds", kind);
  if (!ended) error("Invalid %s stream (truncated stream)", kind);
  if (left != 0) error("Invalid %s stream (trailing or concatenated stream)", kind);
  if (known && produced != (size_t)wanted) {
    error("Session chunk uncompressedBytes does not match the decoded payload (%.0f declared, %.0f decoded)",
          wanted, (double)produced);
  }
  SEXP result = PROTECT(allocVector(RAWSXP, (R_xlen_t)produced));
  if (produced > 0) memcpy(RAW(result), out, produced);
  UNPROTECT(1);
  return result;
}

/* --- Varints ------------------------------------------------------------------ */

/* One canonical varint at *offset (advanced past it), up to 8 bytes. */
static double read_uvarint(const unsigned char *data, R_xlen_t length, R_xlen_t *offset, const char *label) {
  if (*offset < 0 || *offset >= length) error("%s: offset is outside the payload", label);
  double value = 0;
  double scale = 1;
  int count = 0;
  for (;;) {
    if (*offset >= length) error("%s: truncated", label);
    unsigned char byte = data[(*offset)++];
    count++;
    value += (double)(byte & 0x7F) * scale;
    if ((byte & 0x80) == 0) {
      if (value > MAX_SAFE_INTEGER) error("%s: value exceeds the exact integer range", label);
      int minimum = 1;
      double reach = 128;
      while (value >= reach) {
        minimum++;
        reach *= 128;
      }
      if (count != minimum) error("%s: noncanonical encoding", label);
      return value;
    }
    scale *= 128;
    if (count >= 8) error("%s: varint too long", label);
  }
}

/*
 * Check that data[start, length) holds exactly `count` canonical varints of
 * at most 7 bytes (the array rule of session_codecs.py), with the reasons in
 * its order: trailing bytes or too few varints first, then a varint too long,
 * then a redundant byte. Nothing is allocated, so a count read from the
 * payload is checked against the bytes before any array is made from it.
 */
static void check_uvarint_array(const unsigned char *data, R_xlen_t start, R_xlen_t length, double count,
                                const char *label) {
  R_xlen_t size = length - start;
  if (count == 0) {
    if (size != 0) error("%s: trailing bytes", label);
    return;
  }
  double ends = 0;
  for (R_xlen_t position = start; position < length; position++) {
    if ((data[position] & 0x80) == 0) ends++;
  }
  if (ends > count || (ends == count && (data[length - 1] & 0x80) != 0)) error("%s: trailing bytes", label);
  if (ends < count) error("%s: truncated (expected %.0f varints, found %.0f)", label, count, ends);
  int bytes = 0;
  for (R_xlen_t position = start; position < length; position++) {
    bytes++;
    if ((data[position] & 0x80) == 0) {
      if (bytes > 7) error("%s: varint too long", label);
      bytes = 0;
    }
  }
  bytes = 0;
  for (R_xlen_t position = start; position < length; position++) {
    bytes++;
    if ((data[position] & 0x80) == 0) {
      if (bytes > 1 && data[position] == 0) error("%s: noncanonical varint encoding", label);
      bytes = 0;
    }
  }
}

/* The `count` varints check_uvarint_array() accepted, into `values`. */
static void read_uvarint_array(const unsigned char *data, R_xlen_t start, R_xlen_t count, double *values) {
  R_xlen_t position = start;
  for (R_xlen_t index = 0; index < count; index++) {
    double value = 0;
    double scale = 1;
    for (;;) {
      unsigned char byte = data[position++];
      value += (double)(byte & 0x7F) * scale;
      if ((byte & 0x80) == 0) break;
      scale *= 128;
    }
    values[index] = value;
  }
}

/* --- Cell sets -------------------------------------------------------------- */

/* Members of a little-endian bitset of `count` bits, 0-based, as doubles; no
 * bit set past it. */
static SEXP bitset_members(const unsigned char *bits, double count, const char *label) {
  R_xlen_t bytes = (R_xlen_t)ceil(count / 8);
  R_xlen_t members = 0;
  for (R_xlen_t byte = 0; byte < bytes; byte++) {
    unsigned value = bits[byte];
    while (value != 0) {
      value &= value - 1;
      members++;
    }
  }
  R_xlen_t tail = (R_xlen_t)fmod(count, 8);
  if (tail != 0 && bytes > 0 && (bits[bytes - 1] >> tail) != 0) {
    error("%s: bits are set past its length", label);
  }
  SEXP out = PROTECT(allocVector(REALSXP, members));
  double *values = REAL(out);
  R_xlen_t filled = 0;
  for (R_xlen_t byte = 0; byte < bytes; byte++) {
    unsigned value = bits[byte];
    for (int bit = 0; value != 0 && bit < 8; bit++) {
      if (value & (1u << bit)) values[filled++] = (double)byte * 8 + bit;
    }
  }
  UNPROTECT(1);
  return out;
}

/*
 * One cell set (codecs/cell-set.js decodeCellSet), as 0-based indices below
 * `cell_count` (at most 2^32), strictly increasing, as doubles. Byte 0 is
 * the encoding: 0 sorted-delta, 1 bitset, 2 runs; any other is refused.
 */
SEXP cellucid_decode_cell_set(SEXP payload, SEXP cell_count) {
  if (TYPEOF(payload) != RAWSXP) error("decode_cell_set: payload must be a raw vector");
  double cells = asReal(cell_count);
  if (!(cells >= 0) || cells > 4294967296.0 || cells != floor(cells)) {
    error("decode_cell_set: cell_count must be a whole number from 0 to 2^32");
  }
  const unsigned char *data = RAW(payload);
  R_xlen_t length = XLENGTH(payload);
  if (length == 0) error("Cell set must be a non-empty Uint8Array.");
  int encoding = data[0];

  if (encoding == 1) {
    if ((double)length != 1 + ceil(cells / 8)) {
      error("decode_cell_set: bitset has the wrong length for the dataset cell count");
    }
    return bitset_members(data + 1, cells, "decode_cell_set bitset");
  }

  if (encoding == 0) {
    R_xlen_t offset = 1;
    double count_double = read_uvarint(data, length, &offset, "decode_uvarint");
    if (count_double > cells) {
      error("decode_cell_set: %.0f cells exceed the dataset's %.0f", count_double, cells);
    }
    if (count_double > (double)R_XLEN_T_MAX) error("decode_cell_set: %.0f cells are more than one R vector holds", count_double);
    R_xlen_t count = (R_xlen_t)count_double;
    check_uvarint_array(data, offset, length, count_double, "decode_cell_set");
    SEXP out = PROTECT(allocVector(REALSXP, count));
    double *values = REAL(out);
    read_uvarint_array(data, offset, count, values);
    double running = 0;
    for (R_xlen_t index = 0; index < count; index++) {
      if (index > 0 && values[index] == 0) error("decode_cell_set: indices must be strictly increasing");
      running += values[index];
      if (running >= cells) error("decode_cell_set: index exceeds the dataset cell count");
      values[index] = running;
    }
    UNPROTECT(1);
    return out;
  }

  if (encoding == 2) {
    R_xlen_t offset = 1;
    double runs_double = read_uvarint(data, length, &offset, "decode_uvarint");
    if (2 * runs_double > (double)(length - offset)) error("decode_cell_set runs: truncated");
    check_uvarint_array(data, offset, length, 2 * runs_double, "decode_cell_set runs");
    R_xlen_t runs = (R_xlen_t)runs_double;
    double *pairs = (double *)R_alloc(runs > 0 ? 2 * runs : 1, sizeof(double));
    read_uvarint_array(data, offset, 2 * runs, pairs);
    double total = 0;
    double end = 0;
    for (R_xlen_t run = 0; run < runs; run++) {
      double gap = pairs[2 * run];
      double size = pairs[2 * run + 1];
      if (size < 1 || (run > 0 && gap < 1)) error("decode_cell_set: runs must be non-empty and separated");
      end += gap + size;
      if (end > cells) error("decode_cell_set: run exceeds the dataset cell count");
      total += size;
    }
    if (total > (double)R_XLEN_T_MAX) error("decode_cell_set: %.0f cells are more than one R vector holds", total);
    SEXP out = PROTECT(allocVector(REALSXP, (R_xlen_t)total));
    double *values = REAL(out);
    R_xlen_t filled = 0;
    double position = 0;
    for (R_xlen_t run = 0; run < runs; run++) {
      position += pairs[2 * run];
      R_xlen_t size = (R_xlen_t)pairs[2 * run + 1];
      for (R_xlen_t step = 0; step < size; step++) values[filled++] = position + (double)step;
      position += (double)size;
    }
    UNPROTECT(1);
    return out;
  }

  error("decode_cell_set: encoding %d is not one this version reads", encoding);
  return R_NilValue;
}
