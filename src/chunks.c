/*
 * The bytes of one Zarr chunk, made in one pass from an R vector: the values
 * first .. first + count - 1 written little-endian in the store's dtype, then
 * zeros up to the chunk's full length (every fill value this package writes
 * for a fixed-width dtype is zero). One allocation, the raw vector returned,
 * and no intermediate copy: base R's writeBin() of a slice makes the slice,
 * then the bytes, and a NaN needs a second pass over them.
 *
 * And the gene-major arrays of a genes x cells dgCMatrix, one block of
 * positions at a time, made without transposing the whole matrix
 * (cellucid_gene_major_block below; R/expression-input.R).
 */

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

enum dtype { F4, F8, I4, U4, I2, I1, B1, I8 };

static enum dtype parse_dtype(SEXP dtype) {
  if (TYPEOF(dtype) != STRSXP || XLENGTH(dtype) != 1) Rf_error("A chunk's dtype must be one string.");
  const char *name = CHAR(STRING_ELT(dtype, 0));
  if (!strcmp(name, "<f4")) return F4;
  if (!strcmp(name, "<f8")) return F8;
  if (!strcmp(name, "<i4")) return I4;
  if (!strcmp(name, "<u4")) return U4;
  if (!strcmp(name, "<i2")) return I2;
  if (!strcmp(name, "|i1")) return I1;
  if (!strcmp(name, "|b1")) return B1;
  if (!strcmp(name, "<i8")) return I8;
  Rf_error("Unknown dtype %s.", name);
  return F4;
}

static int dtype_size(enum dtype type) {
  switch (type) {
    case F4: case I4: case U4: return 4;
    case F8: case I8: return 8;
    case I2: return 2;
    default: return 1;
  }
}

/*
 * Sizes are bounded by R_XLEN_T_MAX, the longest vector this R allocates, and
 * checked in double before anything is cast or multiplied in R_xlen_t: on
 * 32-bit R (the i386 build of Windows R 4.1) R_xlen_t is an int, where a
 * larger double has no value to cast to and a product wraps to a shorter
 * vector than is then written.
 */

/* A count or position R passes as a double: a whole number R_xlen_t holds. */
static R_xlen_t xlength_of(double value, const char *what) {
  if (!(value >= 0) || value > (double)R_XLEN_T_MAX || value != floor(value)) {
    Rf_error("%s must be a whole number from 0 to %.0f.", what, (double)R_XLEN_T_MAX);
  }
  return (R_xlen_t)value;
}

/* The length of a raw vector of `count` items of `size` bytes. */
static R_xlen_t raw_length(double count, double size, const char *what) {
  double bytes = count * size;
  if (!(bytes >= 0) || bytes > (double)R_XLEN_T_MAX) {
    Rf_error("%s would take %.0f bytes, more than one R raw vector holds.", what, bytes);
  }
  return (R_xlen_t)bytes;
}

/* Little-endian stores that do not depend on the host's byte order. */
static inline void put32(unsigned char *at, uint32_t v) {
  at[0] = (unsigned char)v;
  at[1] = (unsigned char)(v >> 8);
  at[2] = (unsigned char)(v >> 16);
  at[3] = (unsigned char)(v >> 24);
}

static inline void put64(unsigned char *at, uint64_t v) {
  put32(at, (uint32_t)v);
  put32(at + 4, (uint32_t)(v >> 32));
}

/*
 * R's NA and NaN are NaNs with platform-dependent payloads and sign: every one
 * is written as numpy's canonical quiet NaN (0x7FC00000, 0x7FF8000000000000),
 * so the bytes do not depend on the platform.
 */
static inline void put_f4(unsigned char *at, double v) {
  if (isnan(v)) {
    put32(at, 0x7FC00000u);
    return;
  }
  float f = (float)v;
  uint32_t bits;
  memcpy(&bits, &f, 4);
  put32(at, bits);
}

static inline void put_f8(unsigned char *at, double v) {
  if (isnan(v)) {
    put64(at, 0x7FF8000000000000ull);
    return;
  }
  uint64_t bits;
  memcpy(&bits, &v, 8);
  put64(at, bits);
}

/*
 * A whole number held as a double, as the integer dtypes store it: the value
 * modulo 2^width in two's complement (as.integer() then writeBin() writes an
 * unsigned value of 2^31 or more as its negative twin), and NaN as R's
 * NA_integer_ (0x80000000).
 */
static inline int64_t whole(double v) {
  if (isnan(v)) return INT32_MIN;
  return (int64_t)v;
}

/* One value of `source` (an R vector) at `k`, as a double. */
static inline double value_at(SEXP source, int type, R_xlen_t k) {
  if (type == REALSXP) return REAL(source)[k];
  int v = type == INTSXP ? INTEGER(source)[k] : LOGICAL(source)[k];
  return v == NA_INTEGER ? NAN : (double)v;
}

SEXP cellucid_encode_values(SEXP values, SEXP first, SEXP count, SEXP chunk_length, SEXP dtype) {
  int type = TYPEOF(values);
  if (type != REALSXP && type != INTSXP && type != LGLSXP) {
    Rf_error("A chunk's values must be numbers or logicals.");
  }
  enum dtype out_type = parse_dtype(dtype);
  R_xlen_t n = xlength_of(Rf_asReal(count), "A chunk's value count");
  R_xlen_t length = xlength_of(Rf_asReal(chunk_length), "A chunk's length");
  int size = dtype_size(out_type);
  R_xlen_t byte_length = raw_length((double)length, size, "A chunk");
  double first_value = Rf_asReal(first);
  if (!(first_value >= 1) || n > length || first_value - 1 + (double)n > (double)XLENGTH(values)) {
    Rf_error("A chunk's values lie outside the vector.");
  }
  R_xlen_t start = (R_xlen_t)first_value - 1;
  SEXP out = PROTECT(Rf_allocVector(RAWSXP, byte_length));
  unsigned char *bytes = RAW(out);
  memset(bytes + n * size, 0, (size_t)((length - n) * size));
  /* The common cases at memory speed: int32 to int32, doubles to floats. */
  if (out_type == I4 && type == INTSXP) {
    const int *x = INTEGER(values) + start;
    for (R_xlen_t k = 0; k < n; k++) put32(bytes + 4 * k, (uint32_t)x[k]);
  } else if (out_type == F4 && type == REALSXP) {
    const double *x = REAL(values) + start;
    for (R_xlen_t k = 0; k < n; k++) put_f4(bytes + 4 * k, x[k]);
  } else if (out_type == F8 && type == REALSXP) {
    const double *x = REAL(values) + start;
    for (R_xlen_t k = 0; k < n; k++) put_f8(bytes + 8 * k, x[k]);
  } else {
    for (R_xlen_t k = 0; k < n; k++) {
      double v = value_at(values, type, start + k);
      unsigned char *at = bytes + size * k;
      switch (out_type) {
        case F4: put_f4(at, v); break;
        case F8: put_f8(at, v); break;
        case I4: case U4: put32(at, (uint32_t)whole(v)); break;
        case I8: put64(at, (uint64_t)whole(v)); break;
        case I2: {
          uint16_t w = (uint16_t)whole(v);
          at[0] = (unsigned char)w;
          at[1] = (unsigned char)(w >> 8);
          break;
        }
        case I1: at[0] = (unsigned char)whole(v); break;
        case B1: at[0] = (unsigned char)(!isnan(v) && v != 0); break;
      }
    }
  }
  UNPROTECT(1);
  return out;
}

/*
 * Positions first .. last - 1 (0-based, last exclusive) of the gene-major
 * arrays of a genes x cells dgCMatrix (`i`, `p`, `x`) over the genes it
 * keeps: the transpose's row indices (the cells) as little-endian integers of
 * `index_size` bytes, and its values in `data_dtype` ("<f4" or "<f8").
 * `map` gives each row of the matrix its output gene (0-based) or -1 when it
 * is left out, `origin` each output gene its row; `pointers` are the
 * transpose's column pointers over the kept genes (cellucid_row_pointers),
 * n_genes + 1 doubles. Keeping genes (gene_identifiers) copies nothing: the
 * map is read where a copy of the kept rows used to be made.
 *
 * The block holds the output genes g_low .. g_high whose values reach into
 * it; their rows lie between the smallest and the largest of their origins.
 * A dgCMatrix keeps each column's row indices strictly ascending, so every
 * cell's values in that range are one run of its column, found by a binary
 * search; rows outside the block's genes are passed over. The cells are
 * visited in order, so each gene's values come out in its cells' order,
 * which is the transpose's. When the genes are kept in the matrix's own
 * order (every gene, or a subset in order) the range holds only the block's
 * rows and the rows left out between them.
 */
SEXP cellucid_gene_major_block(SEXP i, SEXP p, SEXP x, SEXP map, SEXP origin, SEXP pointers,
                               SEXP first, SEXP last, SEXP index_size, SEXP data_dtype) {
  if (TYPEOF(i) != INTSXP || TYPEOF(p) != INTSXP || TYPEOF(x) != REALSXP || TYPEOF(pointers) != REALSXP ||
      TYPEOF(map) != INTSXP || TYPEOF(origin) != INTSXP) {
    Rf_error("The gene-major block needs a dgCMatrix's i, p and x, the row map and origins, and the pointers as doubles.");
  }
  enum dtype out_type = parse_dtype(data_dtype);
  if (out_type != F4 && out_type != F8) Rf_error("Expression values are stored as <f4 or <f8.");
  int value_size = dtype_size(out_type);
  int cell_size = Rf_asInteger(index_size);
  if (cell_size != 4 && cell_size != 8) Rf_error("Cell indices are 4 or 8 bytes.");
  const int *rows = INTEGER(i);
  const int *columns = INTEGER(p);
  const double *values = REAL(x);
  const int *target = INTEGER(map);
  const int *source = INTEGER(origin);
  const double *pointer = REAL(pointers);
  R_xlen_t n_genes = XLENGTH(pointers) - 1;
  R_xlen_t n_cells = XLENGTH(p) - 1;
  if (XLENGTH(origin) != n_genes) Rf_error("The origins name %d genes, the pointers %d.", (int)XLENGTH(origin), (int)n_genes);
  int64_t low = (int64_t)Rf_asReal(first);
  int64_t high = (int64_t)Rf_asReal(last);
  if (n_genes < 1 || low < 0 || high < low || high > (int64_t)pointer[n_genes]) {
    Rf_error("The block lies outside the matrix's values.");
  }
  R_xlen_t n = raw_length((double)(high - low), 1, "A gene-major block");
  SEXP cells_out = PROTECT(Rf_allocVector(RAWSXP, raw_length((double)n, cell_size, "A block's cell indices")));
  SEXP values_out = PROTECT(Rf_allocVector(RAWSXP, raw_length((double)n, value_size, "A block's values")));
  if (n > 0) {
    /* The first gene with a value at `low`, and the last at `high - 1`. */
    R_xlen_t g_low = 0, g_high;
    {
      R_xlen_t a = 0, b = n_genes; /* pointer[a] <= low < pointer[b] */
      while (b - a > 1) {
        R_xlen_t m = a + (b - a) / 2;
        if (pointer[m] <= (double)low) a = m; else b = m;
      }
      g_low = a;
      a = g_low;
      b = n_genes;
      while (b - a > 1) {
        R_xlen_t m = a + (b - a) / 2;
        if (pointer[m] <= (double)(high - 1)) a = m; else b = m;
      }
      g_high = a;
    }
    R_xlen_t n_block_genes = g_high - g_low + 1;
    int r_low = source[g_low], r_high = source[g_low];
    int64_t *next = (int64_t *)R_alloc((size_t)n_block_genes, sizeof(int64_t));
    for (R_xlen_t g = 0; g < n_block_genes; g++) {
      next[g] = (int64_t)pointer[g_low + g];
      int row = source[g_low + g];
      if (row < r_low) r_low = row;
      if (row > r_high) r_high = row;
    }
    unsigned char *cell_bytes = RAW(cells_out);
    unsigned char *value_bytes = RAW(values_out);
    for (R_xlen_t c = 0; c < n_cells; c++) {
      int k = columns[c], end = columns[c + 1];
      if (k == end || rows[end - 1] < r_low || rows[k] > r_high) continue;
      /* The column's first value at a row >= r_low. */
      int a = k, b = end;
      while (a < b) {
        int m = a + (b - a) / 2;
        if (rows[m] < r_low) a = m + 1; else b = m;
      }
      for (k = a; k < end && rows[k] <= r_high; k++) {
        int gene = target[rows[k]];
        if (gene < g_low || gene > g_high) continue;
        int64_t position = next[gene - g_low]++;
        if (position < low || position >= high) continue;
        R_xlen_t at = (R_xlen_t)(position - low);
        if (cell_size == 4) put32(cell_bytes + 4 * at, (uint32_t)c);
        else put64(cell_bytes + 8 * at, (uint64_t)c);
        if (value_size == 4) put_f4(value_bytes + 4 * at, values[k]);
        else put_f8(value_bytes + 8 * at, values[k]);
      }
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, cells_out);
  SET_VECTOR_ELT(out, 1, values_out);
  UNPROTECT(3);
  return out;
}

/*
 * Bytes offset .. offset + count - 1 (0-based) of a raw vector, then zeros up
 * to `total` bytes: one chunk cut from a gene-major block (one copy, where
 * R's own subsetting would first build an index vector as long as the chunk).
 */
SEXP cellucid_raw_slice(SEXP bytes, SEXP offset, SEXP count, SEXP total) {
  if (TYPEOF(bytes) != RAWSXP) Rf_error("A slice is cut from a raw vector.");
  R_xlen_t from = xlength_of(Rf_asReal(offset), "A slice's offset");
  R_xlen_t n = xlength_of(Rf_asReal(count), "A slice's length");
  R_xlen_t size = xlength_of(Rf_asReal(total), "A slice's padded length");
  if (n > size || (double)from + (double)n > (double)XLENGTH(bytes)) Rf_error("The slice lies outside the bytes.");
  SEXP out = PROTECT(Rf_allocVector(RAWSXP, size));
  if (n > 0) memcpy(RAW(out), RAW(bytes) + from, (size_t)n);
  if (size > n) memset(RAW(out) + n, 0, (size_t)(size - n));
  UNPROTECT(1);
  return out;
}

/*
 * One binary command frame of the viewer channel (R/viewer-protocol.R): the
 * header's length (uint32, little-endian), the header (UTF-8 JSON, padded by
 * the caller to a multiple of 4 bytes), then `cells` as little-endian uint32.
 * The cells are checked in the same pass that writes them: whole numbers from
 * 0 to 2^32 - 1, strictly ascending, at least one (or none, `allow_empty`,
 * for a setVisibility of every cell). NULL when they are not, for the caller
 * to say so; one allocation, the frame, either way.
 */
SEXP cellucid_cells_frame(SEXP header, SEXP cells, SEXP allow_empty) {
  if (TYPEOF(header) != RAWSXP) Rf_error("The frame's header must be bytes.");
  int type = TYPEOF(cells);
  if (type != INTSXP && type != REALSXP && !(Rf_isNull(cells) && Rf_asLogical(allow_empty) == TRUE)) {
    return R_NilValue;
  }
  R_xlen_t n = Rf_isNull(cells) ? 0 : XLENGTH(cells);
  if (n == 0 && Rf_asLogical(allow_empty) != TRUE) return R_NilValue;
  R_xlen_t header_length = XLENGTH(header);
  /* The header's length is a uint32 of the frame (compared in 64 bits). */
  if ((uint64_t)header_length > UINT32_MAX) Rf_error("The frame's header is too long.");
  SEXP out = PROTECT(Rf_allocVector(RAWSXP, raw_length(4 + (double)header_length + 4 * (double)n, 1, "A command frame")));
  unsigned char *bytes = RAW(out);
  put32(bytes, (uint32_t)header_length);
  memcpy(bytes + 4, RAW(header), (size_t)header_length);
  unsigned char *body = bytes + 4 + header_length;
  double previous = -1;
  if (type == INTSXP) {
    const int *x = INTEGER(cells);
    for (R_xlen_t k = 0; k < n; k++) {
      int v = x[k];
      if (v == NA_INTEGER || v < 0 || (double)v <= previous) {
        UNPROTECT(1);
        return R_NilValue;
      }
      previous = v;
      put32(body + 4 * k, (uint32_t)v);
    }
  } else if (type == REALSXP) {
    const double *x = REAL(cells);
    for (R_xlen_t k = 0; k < n; k++) {
      double v = x[k];
      if (!(v >= 0) || v > 4294967295.0 || v != floor(v) || v <= previous) {
        UNPROTECT(1);
        return R_NilValue;
      }
      previous = v;
      put32(body + 4 * k, (uint32_t)v);
    }
  }
  UNPROTECT(1);
  return out;
}
