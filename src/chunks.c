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
 * which is the transpose's; with `cells`, the cells are visited in that
 * order and output row `c` is the `c`-th of them (a store in priority order:
 * each gene's rows ascend by priority, made from the matrix as it is, never
 * from a reordered copy). When the genes are kept in the matrix's own
 * order (every gene, or a subset in order) the range holds only the block's
 * rows and the rows left out between them.
 */
SEXP cellucid_gene_major_block(SEXP i, SEXP p, SEXP x, SEXP map, SEXP origin, SEXP pointers, SEXP cells,
                               SEXP first, SEXP last, SEXP index_size, SEXP data_dtype) {
  if (TYPEOF(i) != INTSXP || TYPEOF(p) != INTSXP || TYPEOF(x) != REALSXP || TYPEOF(pointers) != REALSXP ||
      TYPEOF(map) != INTSXP || TYPEOF(origin) != INTSXP || TYPEOF(cells) != INTSXP) {
    Rf_error("The gene-major block needs a dgCMatrix's i, p and x, the row map and origins, the cell order and the pointers as doubles.");
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
  /* The cells in the order the output rows take them (0-based columns): the
   * priority order of a store in priority order, none for the matrix's own. */
  const int *order = XLENGTH(cells) > 0 ? INTEGER(cells) : NULL;
  if (order != NULL && XLENGTH(cells) != n_cells) Rf_error("The cell order names %d cells, the matrix %d.", (int)XLENGTH(cells), (int)n_cells);
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
      R_xlen_t column = order == NULL ? c : order[c];
      int k = columns[column], end = columns[column + 1];
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


/* One (row, value) pair of a column being put in row order. */
struct renumbered { int64_t row; double value; };

/*
 * The `n` pairs in ascending row (rows are distinct and below 2^`bits`), into
 * `pairs`; `spare` is as long. A short column (a graph's row) by insertion;
 * a long one (a gene) by a least-significant-digit radix sort, a byte of the
 * row a pass: about ten times faster than qsort() on a gene of 10^4 values,
 * which made priority order cost as much as the rest of writing X.
 */
static void sort_by_row(struct renumbered *pairs, struct renumbered *spare, R_xlen_t n, int bits) {
  if (n < 64) {
    for (R_xlen_t k = 1; k < n; k++) {
      struct renumbered item = pairs[k];
      R_xlen_t m = k;
      while (m > 0 && pairs[m - 1].row > item.row) {
        pairs[m] = pairs[m - 1];
        m--;
      }
      pairs[m] = item;
    }
    return;
  }
  struct renumbered *from = pairs, *to = spare;
  for (int shift = 0; shift < bits; shift += 8) {
    R_xlen_t count[257] = {0};
    for (R_xlen_t k = 0; k < n; k++) count[((from[k].row >> shift) & 0xFF) + 1]++;
    for (int d = 0; d < 256; d++) count[d + 1] += count[d];
    for (R_xlen_t k = 0; k < n; k++) to[count[(from[k].row >> shift) & 0xFF]++] = from[k];
    struct renumbered *swap = from;
    from = to;
    to = swap;
  }
  if (from != pairs) memcpy(pairs, from, (size_t)n * sizeof(struct renumbered));
}

/*
 * Positions first .. last - 1 (0-based) of the column-compressed arrays of a
 * dgCMatrix (`i`, `p`, `x`) whose rows are renumbered -- row `r` becomes
 * `new_row[r]` (0-based) -- and whose columns are taken in the order
 * `column_order` gives (0-based; empty for the matrix's own): output column
 * `g` is column `column_order[g]`, its values ascending by new row, and
 * `pointers` (integers, one more than the columns) are the output's column
 * pointers (`p` itself when the columns keep their order). Two uses, both
 * made from the matrix as it is, never from a reordered copy:
 *   - a cells x genes expression matrix (a column per gene) in priority
 *     order: its cells renumbered, its genes in their order;
 *   - a symmetric neighbour graph (obsp) in priority order: both axes
 *     permuted, row r of the output being column order[r] of the graph (its
 *     column-compressed arrays are its row-compressed ones), each row's
 *     columns renumbered and sorted.
 * Each output column reaching into the block is read whole and sorted
 * (R/expression-input.R sizes blocks so that a column is sorted about once).
 * `part` 0 returns list(row index bytes, value bytes); 1 or 2 returns only
 * the row indices or only the values, as one chunk: zeros after them up to
 * `padded` items (a graph's chunk, made when it is read and let go: one
 * allocation, the chunk itself).
 */
SEXP cellucid_renumbered_block(SEXP i, SEXP p, SEXP x, SEXP new_row, SEXP column_order, SEXP pointers,
                               SEXP first, SEXP last, SEXP index_size, SEXP data_dtype, SEXP part, SEXP padded) {
  if (TYPEOF(i) != INTSXP || TYPEOF(p) != INTSXP || TYPEOF(x) != REALSXP || TYPEOF(new_row) != INTSXP ||
      TYPEOF(column_order) != INTSXP || TYPEOF(pointers) != INTSXP) {
    Rf_error("The renumbered block needs a dgCMatrix's i, p and x, the new row of each row, the column order and the output pointers.");
  }
  enum dtype out_type = parse_dtype(data_dtype);
  if (out_type != F4 && out_type != F8) Rf_error("Stored values are <f4 or <f8.");
  int value_size = dtype_size(out_type);
  int cell_size = Rf_asInteger(index_size);
  if (cell_size != 4 && cell_size != 8) Rf_error("Row indices are 4 or 8 bytes.");
  const int *rows = INTEGER(i);
  const int *columns = INTEGER(p);
  const double *values = REAL(x);
  const int *renumber = INTEGER(new_row);
  R_xlen_t n_rows = XLENGTH(new_row);
  const int *pointer = INTEGER(pointers);
  R_xlen_t n_columns = XLENGTH(p) - 1;
  const int *source = XLENGTH(column_order) > 0 ? INTEGER(column_order) : NULL;
  if (source != NULL && XLENGTH(column_order) != n_columns) {
    Rf_error("The column order names %d columns, the matrix %d.", (int)XLENGTH(column_order), (int)n_columns);
  }
  if (XLENGTH(pointers) != n_columns + 1) Rf_error("The output pointers name %d columns, the matrix %d.", (int)XLENGTH(pointers) - 1, (int)n_columns);
  int64_t low = (int64_t)Rf_asReal(first);
  int64_t high = (int64_t)Rf_asReal(last);
  if (low < 0 || high < low || n_columns < 1 || high > (int64_t)pointer[n_columns] ||
      pointer[n_columns] != columns[n_columns]) {
    Rf_error("The block lies outside the matrix's values.");
  }
  R_xlen_t n = raw_length((double)(high - low), 1, "A block");
  int which = Rf_asInteger(part);
  if (which < 0 || which > 2) Rf_error("A block's part is 0 (both), 1 (row indices) or 2 (values).");
  R_xlen_t length = which == 0 ? n : xlength_of(Rf_asReal(padded), "A chunk's length");
  if (length < n) Rf_error("A chunk of %.0f items cannot hold %.0f.", (double)length, (double)n);
  SEXP cells_out = PROTECT(which == 2 ? R_NilValue : Rf_allocVector(RAWSXP, raw_length((double)length, cell_size, "A block's row indices")));
  SEXP values_out = PROTECT(which == 1 ? R_NilValue : Rf_allocVector(RAWSXP, raw_length((double)length, value_size, "A block's values")));
  unsigned char *cell_bytes = which == 2 ? NULL : RAW(cells_out);
  unsigned char *value_bytes = which == 1 ? NULL : RAW(values_out);
  if (cell_bytes != NULL && length > n) memset(cell_bytes + n * cell_size, 0, (size_t)((length - n) * cell_size));
  if (value_bytes != NULL && length > n) memset(value_bytes + n * value_size, 0, (size_t)((length - n) * value_size));
  if (n > 0) {
    R_xlen_t a = 0, b = n_columns; /* pointer[a] <= low < pointer[b] */
    while (b - a > 1) {
      R_xlen_t m = a + (b - a) / 2;
      if ((int64_t)pointer[m] <= low) a = m; else b = m;
    }
    /* One buffer, as long as the longest column reaching into the block, used
       by each column in turn (R_alloc'd per column, every column's pairs would
       be held until the call returns: 16 bytes a value of the whole block). */
    int longest = 0;
    for (R_xlen_t g = a; g < n_columns && (int64_t)pointer[g] < high; g++) {
      R_xlen_t column = source == NULL ? g : source[g];
      if (column < 0 || column >= n_columns || columns[column + 1] - columns[column] != pointer[g + 1] - pointer[g]) {
        Rf_error("The output pointers do not match the matrix's columns in that order.");
      }
      if (columns[column + 1] - columns[column] > longest) longest = columns[column + 1] - columns[column];
    }
    struct renumbered *pairs = longest > 0 ? (struct renumbered *)R_alloc((size_t)longest, sizeof(struct renumbered)) : NULL;
    struct renumbered *spare = longest >= 64 ? (struct renumbered *)R_alloc((size_t)longest, sizeof(struct renumbered)) : NULL;
    /* The bits a new row takes: rows are below n_rows. */
    int bits = 0;
    while (bits < 63 && ((int64_t)1 << bits) < (int64_t)n_rows) bits++;
    for (R_xlen_t g = a; g < n_columns && (int64_t)pointer[g] < high; g++) {
      R_xlen_t column = source == NULL ? g : source[g];
      int start = columns[column], end = columns[column + 1];
      if (end == start) continue;
      for (int k = start; k < end; k++) {
        if (rows[k] < 0 || rows[k] >= n_rows) Rf_error("A row index lies outside the %d rows renumbered.", (int)n_rows);
        if (renumber[rows[k]] < 0 || renumber[rows[k]] >= n_rows) Rf_error("A new row lies outside the %d rows.", (int)n_rows);
        pairs[k - start].row = renumber[rows[k]];
        pairs[k - start].value = values[k];
      }
      sort_by_row(pairs, spare, end - start, bits);
      int64_t base = pointer[g];
      for (int k = 0; k < end - start; k++) {
        int64_t position = base + k;
        if (position < low || position >= high) continue;
        R_xlen_t at = (R_xlen_t)(position - low);
        if (cell_bytes != NULL) {
          if (cell_size == 4) put32(cell_bytes + 4 * at, (uint32_t)pairs[k].row);
          else put64(cell_bytes + 8 * at, (uint64_t)pairs[k].row);
        }
        if (value_bytes != NULL) {
          if (value_size == 4) put_f4(value_bytes + 4 * at, pairs[k].value);
          else put_f8(value_bytes + 8 * at, pairs[k].value);
        }
      }
    }
  }
  if (which != 0) {
    UNPROTECT(2);
    return which == 1 ? cells_out : values_out;
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, cells_out);
  SET_VECTOR_ELT(out, 1, values_out);
  UNPROTECT(3);
  return out;
}
