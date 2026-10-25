/*
 * The codecs of a .cellucid-session file (session format 2) that applying a
 * session to an object needs (R/session-apply.R), exactly as the web app
 * writes them and cellucid-python reads them (session_codecs.py):
 *
 * - inflating one gzip (or zlib) stream to exactly its declared size, with no
 *   bytes after it (a truncated, longer, shorter or concatenated stream is an
 *   error, never a hang: base R's memDecompress() does not end on two
 *   concatenated gzip members);
 * - a highlight's cell set (codecs/cell-set.js): sorted deltas, a bitset,
 *   runs, or a mask over a category's cells, with that category's digest;
 * - a user-defined categorical field's codes
 *   (contributors/user-defined-codes.js): raw, runs, or bit-packed.
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

static uint32_t mix(uint32_t x) {
  x = (x ^ (x >> 16)) * 0x45D9F3Bu;
  x = (x ^ (x >> 16)) * 0x45D9F3Bu;
  return x ^ (x >> 16);
}

/* The web app's digestCellSet: the sum, mod 2^32, of one mixed word per index. */
static uint32_t digest_indices(const int *indices, R_xlen_t count) {
  uint32_t sum = 0;
  for (R_xlen_t index = 0; index < count; index++) sum += mix((uint32_t)indices[index]);
  return sum;
}

SEXP cellucid_digest_cell_set(SEXP indices) {
  if (TYPEOF(indices) == INTSXP) return ScalarReal((double)digest_indices(INTEGER(indices), XLENGTH(indices)));
  if (TYPEOF(indices) != REALSXP) error("digest: indices must be numbers");
  uint32_t sum = 0;
  const double *values = REAL(indices);
  for (R_xlen_t index = 0; index < XLENGTH(indices); index++) {
    double value = values[index];
    if (!(value >= 0 && value <= 4294967295.0) || value != (double)(uint32_t)value) {
      error("digest: indices must be whole numbers from 0 to 2^32 - 1");
    }
    sum += mix((uint32_t)value);
  }
  return ScalarReal((double)sum);
}

/* Members of a little-endian bitset of `count` bits, 0-based; no bit set past it. */
static SEXP bitset_members(const unsigned char *bits, R_xlen_t count, const char *label) {
  R_xlen_t bytes = (count + 7) / 8;
  R_xlen_t members = 0;
  for (R_xlen_t index = 0; index < count; index++) {
    if (bits[index >> 3] & (1u << (index & 7))) members++;
  }
  if (count % 8 != 0 && bytes > 0 && (bits[bytes - 1] >> (count % 8)) != 0) {
    error("%s: bits are set past its length", label);
  }
  SEXP out = PROTECT(allocVector(INTSXP, members));
  int *values = INTEGER(out);
  R_xlen_t filled = 0;
  for (R_xlen_t index = 0; index < count; index++) {
    if (bits[index >> 3] & (1u << (index & 7))) values[filled++] = (int)index;
  }
  UNPROTECT(1);
  return out;
}

/*
 * One highlight membership (highlights/cells/<group>), as 0-based cell
 * indices, strictly increasing. Byte 0 is the encoding: 0 sorted-delta, 1
 * bitset, 2 runs, 3 reference-mask, 4 reference-all. Encodings 3 and 4 need
 * `reference`: the 0-based cells of the group's category, ascending.
 */
SEXP cellucid_decode_cell_set(SEXP payload, SEXP cell_count, SEXP reference) {
  if (TYPEOF(payload) != RAWSXP) error("decode_cell_set: payload must be a raw vector");
  double cells_double = asReal(cell_count);
  if (!(cells_double >= 0) || cells_double > INT_MAX) error("decode_cell_set: cell_count must be a non-negative integer");
  R_xlen_t cells = (R_xlen_t)cells_double;
  const unsigned char *data = RAW(payload);
  R_xlen_t length = XLENGTH(payload);
  if (length == 0) error("decode_cell_set: empty payload");
  int encoding = data[0];

  if (encoding == 3 || encoding == 4) {
    if (reference == R_NilValue) {
      error("decode_cell_set: this cell set is stored relative to a category; its cells are required");
    }
    if (TYPEOF(reference) != INTSXP) error("decode_cell_set: reference must be integer cell indices");
    R_xlen_t offset = 1;
    double reference_count = read_uvarint(data, length, &offset, "decode_uvarint");
    if (offset + 4 > length) error("decode_cell_set: truncated");
    uint32_t digest = (uint32_t)data[offset] | ((uint32_t)data[offset + 1] << 8) |
                      ((uint32_t)data[offset + 2] << 16) | ((uint32_t)data[offset + 3] << 24);
    R_xlen_t given = XLENGTH(reference);
    if (reference_count != (double)given || digest != digest_indices(INTEGER(reference), given)) {
      error("The category this selection was saved relative to no longer holds the same cells.");
    }
    R_xlen_t body = offset + 4;
    if (encoding == 4) {
      if (body != length) error("decode_cell_set: trailing bytes");
      return duplicate(reference);
    }
    if (length - body != (given + 7) / 8) error("decode_cell_set: reference mask has the wrong length");
    SEXP members = PROTECT(bitset_members(data + body, given, "decode_cell_set mask"));
    int *positions = INTEGER(members);
    const int *cells_of = INTEGER(reference);
    for (R_xlen_t index = 0; index < XLENGTH(members); index++) positions[index] = cells_of[positions[index]];
    UNPROTECT(1);
    return members;
  }

  if (encoding == 1) {
    if (length != 1 + (cells + 7) / 8) {
      error("decode_cell_set: bitset has the wrong length for the dataset cell count");
    }
    return bitset_members(data + 1, cells, "decode_cell_set bitset");
  }

  if (encoding == 0) {
    R_xlen_t offset = 1;
    double count_double = read_uvarint(data, length, &offset, "decode_uvarint");
    if (count_double > (double)cells) {
      error("decode_cell_set: %.0f cells exceed the dataset's %.0f", count_double, (double)cells);
    }
    R_xlen_t count = (R_xlen_t)count_double;
    check_uvarint_array(data, offset, length, count_double, "decode_cell_set");
    double *deltas = (double *)R_alloc(count > 0 ? count : 1, sizeof(double));
    read_uvarint_array(data, offset, count, deltas);
    SEXP out = PROTECT(allocVector(INTSXP, count));
    int *values = INTEGER(out);
    for (R_xlen_t index = 1; index < count; index++) {
      if (deltas[index] == 0) error("decode_cell_set: indices must be strictly increasing");
    }
    double running = 0;
    for (R_xlen_t index = 0; index < count; index++) {
      running += deltas[index];
      if (running >= (double)cells) error("decode_cell_set: index exceeds the dataset cell count");
      values[index] = (int)running;
    }
    UNPROTECT(1);
    return out;
  }

  if (encoding == 2) {
    R_xlen_t offset = 1;
    double runs_double = read_uvarint(data, length, &offset, "decode_uvarint");
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
      total += size;
    }
    if (runs > 0 && end > (double)cells) error("decode_cell_set: run exceeds the dataset cell count");
    SEXP out = PROTECT(allocVector(INTSXP, (R_xlen_t)total));
    int *values = INTEGER(out);
    R_xlen_t filled = 0;
    double position = 0;
    for (R_xlen_t run = 0; run < runs; run++) {
      position += pairs[2 * run];
      R_xlen_t size = (R_xlen_t)pairs[2 * run + 1];
      for (R_xlen_t step = 0; step < size; step++) values[filled++] = (int)(position + step);
      position += size;
    }
    UNPROTECT(1);
    return out;
  }

  error("decode_cell_set: encoding %d is not one this version reads", encoding);
  return R_NilValue;
}

/* --- User-defined codes ------------------------------------------------------- */

/*
 * One user-defined categorical field's codes (user-defined/codes/<id>), as
 * doubles (a Uint32Array code can exceed R's integers): `u8 encoding`,
 * `uvarint length`, then 0 raw little-endian codes of the declared width, 1
 * runs (`uvarint runCount`, then runCount x (uvarint value, uvarint run)),
 * 2 packed (`u8 bits`, 1, 2 or 4, then the codes from each byte's least
 * significant bit).
 */
SEXP cellucid_decode_user_codes(SEXP payload, SEXP expected_length, SEXP codes_type) {
  if (TYPEOF(payload) != RAWSXP) error("decode_user_defined_codes: payload must be a raw vector");
  const char *type = CHAR(STRING_ELT(codes_type, 0));
  int width;
  double maximum;
  if (strcmp(type, "Uint8Array") == 0) {
    width = 1;
    maximum = 255;
  } else if (strcmp(type, "Uint16Array") == 0) {
    width = 2;
    maximum = 65535;
  } else if (strcmp(type, "Uint32Array") == 0) {
    width = 4;
    maximum = 4294967295.0;
  } else {
    error("decode_user_defined_codes: expected_codes_type must be Uint8Array, Uint16Array or Uint32Array");
  }
  double expected = asReal(expected_length);
  if (!(expected >= 0) || expected > 4294967295.0) {
    error("decode_user_defined_codes: expected_length exceeds Uint32 capacity");
  }
  const unsigned char *data = RAW(payload);
  R_xlen_t length = XLENGTH(payload);
  if (length < 2) error("decode_user_defined_codes: payload too short");
  int encoding = data[0];
  R_xlen_t offset = 1;
  double declared = read_uvarint(data, length, &offset, "decode_uvarint");
  if (declared != expected) {
    error("decode_user_defined_codes: declared length %.0f does not match expected length %.0f", declared, expected);
  }
  /* A Uint32 length can exceed R_xlen_t on 32-bit R, where it is an int. */
  if (declared > (double)R_XLEN_T_MAX) {
    error("decode_user_defined_codes: %.0f codes are more than one R vector holds", declared);
  }
  R_xlen_t count = (R_xlen_t)declared;

  if (encoding == 0) {
    if ((double)(length - offset) != declared * width) {
      error("decode_user_defined_codes: raw payload has the wrong length or trailing bytes");
    }
    SEXP out = PROTECT(allocVector(REALSXP, count));
    double *values = REAL(out);
    for (R_xlen_t index = 0; index < count; index++) {
      const unsigned char *at = data + offset + index * width;
      uint32_t code = at[0];
      if (width >= 2) code |= (uint32_t)at[1] << 8;
      if (width == 4) code |= ((uint32_t)at[2] << 16) | ((uint32_t)at[3] << 24);
      values[index] = (double)code;
    }
    UNPROTECT(1);
    return out;
  }

  if (encoding == 2) {
    if (offset >= length) error("decode_user_defined_codes: packed payload is truncated");
    int bits = data[offset];
    if (bits != 1 && bits != 2 && bits != 4) error("decode_user_defined_codes: packed codes must use 1, 2 or 4 bits");
    R_xlen_t body = offset + 1;
    double needed = floor((declared * bits + 7) / 8);
    if ((double)(length - body) != needed) {
      error("decode_user_defined_codes: packed payload has the wrong length or trailing bytes");
    }
    /* Bit positions reach eight times the payload's length: 64-bit, as
       R_xlen_t is 32-bit on 32-bit R. */
    uint64_t used_bits = (uint64_t)count * (uint64_t)bits;
    if (used_bits % 8 != 0 && (data[length - 1] >> (used_bits % 8)) != 0) {
      error("decode_user_defined_codes: packed codes set bits past their length");
    }
    SEXP out = PROTECT(allocVector(REALSXP, count));
    double *values = REAL(out);
    unsigned mask = (1u << bits) - 1u;
    for (R_xlen_t index = 0; index < count; index++) {
      uint64_t bit = (uint64_t)index * (uint64_t)bits;
      values[index] = (double)((data[body + (R_xlen_t)(bit >> 3)] >> (bit & 7)) & mask);
    }
    UNPROTECT(1);
    return out;
  }

  if (encoding != 1) error("decode_user_defined_codes: unknown encoding type %d", encoding);
  double runs_double = read_uvarint(data, length, &offset, "decode_uvarint");
  if (declared == 0) {
    if (runs_double != 0) error("decode_user_defined_codes: empty runs payload must contain zero runs");
  } else if (runs_double == 0 || runs_double > declared) {
    error("decode_user_defined_codes: run count cannot fill the declared length");
  }
  check_uvarint_array(data, offset, length, 2 * runs_double, "decode_user_defined_codes runs");
  R_xlen_t runs = (R_xlen_t)runs_double;
  double *pairs = (double *)R_alloc(runs > 0 ? 2 * runs : 1, sizeof(double));
  read_uvarint_array(data, offset, 2 * runs, pairs);
  /* Each rule over every run, in session_codecs.py's order. */
  for (R_xlen_t run = 0; run < runs; run++) {
    if (pairs[2 * run] > maximum) error("decode_user_defined_codes: run value exceeds its declared code type");
  }
  for (R_xlen_t run = 0; run < runs; run++) {
    if (pairs[2 * run + 1] == 0) error("decode_user_defined_codes: run length must be positive");
  }
  for (R_xlen_t run = 1; run < runs; run++) {
    if (pairs[2 * run] == pairs[2 * (run - 1)]) error("decode_user_defined_codes: adjacent runs must have distinct values");
  }
  double total = 0;
  for (R_xlen_t run = 0; run < runs; run++) total += pairs[2 * run + 1];
  if (total != declared) error("decode_user_defined_codes: runs do not fill the declared length");
  SEXP out = PROTECT(allocVector(REALSXP, count));
  double *values = REAL(out);
  R_xlen_t filled = 0;
  for (R_xlen_t run = 0; run < runs; run++) {
    R_xlen_t size = (R_xlen_t)pairs[2 * run + 1];
    for (R_xlen_t step = 0; step < size; step++) values[filled++] = pairs[2 * run];
  }
  UNPROTECT(1);
  return out;
}
