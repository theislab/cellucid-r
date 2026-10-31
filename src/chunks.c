/*
 * The bytes of one Zarr chunk, made in one pass from an R vector: the values
 * first .. first + count - 1 written little-endian in the store's dtype, then
 * zeros up to the chunk's full length (every fill value this package writes
 * for a fixed-width dtype is zero). One allocation, the raw vector returned,
 * and no intermediate copy: base R's writeBin() of a slice makes the slice,
 * then the bytes, and a NaN needs a second pass over them.
 *
 * And the gene-major arrays of a genes x cells dgCMatrix, one block of
 * positions at a time, made without transposing the whole matrix, on every
 * core (cellucid_gene_major_block below; R/expression-input.R), and laid out
 * in the background for a served store (cellucid_layouts_start).
 */

#include "gene_major.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef _WIN32
#include <fcntl.h>
#include <pthread.h>
#include <sys/mman.h>
#include <time.h>
#include <unistd.h>
#endif

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
 * rows and the rows left out between them. On `threads` threads the cells
 * are cut into pieces: each piece's values of each gene are counted, so each
 * piece knows where its values of a gene begin, then placed, every piece on
 * its own thread; the bytes are the same on any number of threads.
 */
/* What a gene-major pass reads: a genes x cells dgCMatrix and how its rows map to output genes. */
struct gene_major_input {
  const int *rows, *columns;
  const double *values;
  const int *target, *source, *order;
  const double *pointer;
  R_xlen_t n_genes, n_cells;
  int cell_size, value_size;
};

/*
 * The positions low .. high - 1 a pass makes: the output genes g_low ..
 * g_high whose values reach into them, and the matrix rows r_low .. r_high
 * those genes lie between.
 */
struct gene_block {
  int64_t low, high;
  R_xlen_t g_low, g_high;
  int r_low, r_high;
};

static struct gene_block gene_block_of(const struct gene_major_input *in, int64_t low, int64_t high) {
  const double *pointer = in->pointer;
  struct gene_block b;
  b.low = low;
  b.high = high;
  /* The first gene with a value at `low`, and the last at `high - 1`. */
  R_xlen_t a = 0, z = in->n_genes; /* pointer[a] <= low < pointer[z] */
  while (z - a > 1) {
    R_xlen_t m = a + (z - a) / 2;
    if (pointer[m] <= (double)low) a = m; else z = m;
  }
  b.g_low = a;
  z = in->n_genes;
  while (z - a > 1) {
    R_xlen_t m = a + (z - a) / 2;
    if (pointer[m] <= (double)(high - 1)) a = m; else z = m;
  }
  b.g_high = a;
  b.r_low = b.r_high = in->source[b.g_low];
  for (R_xlen_t g = b.g_low; g <= b.g_high; g++) {
    int row = in->source[g];
    if (row < b.r_low) b.r_low = row;
    if (row > b.r_high) b.r_high = row;
  }
  return b;
}

/* Cells visited between two looks at a stop, a pause and the progress made. */
#define LAYOUT_STEP 256

struct layout_segment;
static int segment_step(struct layout_segment *s, int64_t c);

/*
 * The block's values held by cells from .. to - 1 (in output order): with
 * `cell_bytes` NULL, counted into `next` (one count a block gene); else each
 * placed at its gene's cursor in `next`, and written when its position lies
 * in the block (at cell_bytes and value_bytes, position b->low at their
 * start). A dgCMatrix keeps each column's row indices strictly ascending, so
 * a cell's values in the block's rows are one run of its column, found by a
 * binary search. `s` (may be NULL), a layout's segment, is told every
 * LAYOUT_STEP cells how far the pass is, and may stop or pause it there
 * (segment_step); returns 1 when stopped.
 */
static int gene_major_cells(const struct gene_major_input *in, const struct gene_block *b, int64_t from, int64_t to,
                            int64_t *next, unsigned char *cell_bytes, unsigned char *value_bytes,
                            struct layout_segment *s) {
  const int *rows = in->rows;
  const int *columns = in->columns;
  const R_xlen_t g_low = b->g_low, g_high = b->g_high;
  const int r_low = b->r_low, r_high = b->r_high;
  for (int64_t c = from; c < to; c++) {
    if (s != NULL && ((c - from) & (LAYOUT_STEP - 1)) == 0 && segment_step(s, c)) return 1;
    R_xlen_t column = in->order == NULL ? (R_xlen_t)c : in->order[c];
    int k = columns[column], end = columns[column + 1];
    if (k == end || rows[end - 1] < r_low || rows[k] > r_high) continue;
    /* The column's first value at a row >= r_low. */
    int a = k, z = end;
    while (a < z) {
      int m = a + (z - a) / 2;
      if (rows[m] < r_low) a = m + 1; else z = m;
    }
    if (cell_bytes == NULL) {
      for (k = a; k < end && rows[k] <= r_high; k++) {
        int gene = in->target[rows[k]];
        if (gene >= g_low && gene <= g_high) next[gene - g_low]++;
      }
      continue;
    }
    for (k = a; k < end && rows[k] <= r_high; k++) {
      int gene = in->target[rows[k]];
      if (gene < g_low || gene > g_high) continue;
      int64_t position = next[gene - g_low]++;
      if (position < b->low || position >= b->high) continue;
      R_xlen_t at = (R_xlen_t)(position - b->low);
      if (in->cell_size == 4) put32(cell_bytes + 4 * at, (uint32_t)c);
      else put64(cell_bytes + 8 * at, (uint64_t)c);
      if (in->value_size == 4) put_f4(value_bytes + 4 * at, in->values[k]);
      else put_f8(value_bytes + 8 * at, in->values[k]);
    }
  }
  return 0;
}

/*
 * Work shared among threads: task(argument + k * size) for k = 0 .. n - 1,
 * each on a thread of its own (the first on the caller's), joined before it
 * returns. A thread that could not start has its task run on the caller's
 * thread after. No task touches R.
 */
typedef void (*cellucid_task)(void *);

struct task_start {
  cellucid_task task;
  void *argument;
};

#ifdef _WIN32
static DWORD WINAPI task_thread(LPVOID argument) {
  struct task_start *start = (struct task_start *)argument;
  start->task(start->argument);
  return 0;
}
#else
static void *task_thread(void *argument) {
  struct task_start *start = (struct task_start *)argument;
  start->task(start->argument);
  return NULL;
}
#endif

static void run_tasks(cellucid_task task, void *arguments, size_t size, int n) {
  if (n <= 0) return;
  struct task_start *starts = (struct task_start *)malloc((size_t)n * sizeof(struct task_start));
#ifdef _WIN32
  HANDLE *threads = (HANDLE *)malloc((size_t)n * sizeof(HANDLE));
#else
  pthread_t *threads = (pthread_t *)malloc((size_t)n * sizeof(pthread_t));
#endif
  int *running = (int *)calloc((size_t)n, sizeof(int));
  if (starts == NULL || threads == NULL || running == NULL) {
    /* No room to start threads: every task here, one after another. */
    for (int k = 0; k < n; k++) task((char *)arguments + (size_t)k * size);
  } else {
    for (int k = 1; k < n; k++) {
      starts[k].task = task;
      starts[k].argument = (char *)arguments + (size_t)k * size;
#ifdef _WIN32
      threads[k] = CreateThread(NULL, 0, task_thread, &starts[k], 0, NULL);
      running[k] = threads[k] != NULL;
#else
      running[k] = pthread_create(&threads[k], NULL, task_thread, &starts[k]) == 0;
#endif
    }
    task(arguments);
    for (int k = 1; k < n; k++) {
      if (!running[k]) {
        task((char *)arguments + (size_t)k * size);
        continue;
      }
#ifdef _WIN32
      WaitForSingleObject(threads[k], INFINITE);
      CloseHandle(threads[k]);
#else
      pthread_join(threads[k], NULL);
#endif
    }
  }
  free(starts);
  free(threads);
  free(running);
}

/*
 * One piece of a block made on threads: the cells first .. last - 1, the
 * cursors of the block's genes there (or, while counting, their counts).
 */
struct gene_piece {
  const struct gene_major_input *in;
  const struct gene_block *b;
  int64_t first, last;
  int64_t *next;
  unsigned char *cell_bytes, *value_bytes;
};

static void gene_piece_task(void *argument) {
  struct gene_piece *piece = (struct gene_piece *)argument;
  gene_major_cells(piece->in, piece->b, piece->first, piece->last, piece->next, piece->cell_bytes, piece->value_bytes, NULL);
}

/*
 * The block's positions made from `n` pieces of cells, one a thread: each
 * piece's cursors given (`counted` 1), or (0) its values of each block gene
 * counted first, each piece's cursors then following the previous piece's
 * from the genes' pointers. The bytes do not depend on the pieces: each
 * gene's values still come out in its cells' order. 0, or -1 when the counts
 * disagree with the pointers.
 */
static int gene_major_pieces(struct gene_piece *pieces, int n, int counted) {
  if (n <= 0) return 0;
  const struct gene_major_input *in = pieces[0].in;
  const struct gene_block *b = pieces[0].b;
  R_xlen_t n_block_genes = b->g_high - b->g_low + 1;
  if (!counted) {
    /* Every piece writes the same block: counted with no bytes, then placed into it. */
    unsigned char *cells[2] = {pieces[0].cell_bytes, pieces[0].value_bytes};
    for (int k = 0; k < n; k++) {
      memset(pieces[k].next, 0, (size_t)n_block_genes * sizeof(int64_t));
      pieces[k].cell_bytes = NULL;
    }
    run_tasks(gene_piece_task, pieces, sizeof(struct gene_piece), n);
    for (R_xlen_t g = 0; g < n_block_genes; g++) {
      int64_t cursor = (int64_t)in->pointer[b->g_low + g];
      for (int k = 0; k < n; k++) {
        int64_t count = pieces[k].next[g];
        pieces[k].next[g] = cursor;
        cursor += count;
      }
      if (cursor != (int64_t)in->pointer[b->g_low + g + 1]) return -1;
    }
    for (int k = 0; k < n; k++) {
      pieces[k].cell_bytes = cells[0];
      pieces[k].value_bytes = cells[1];
    }
  }
  run_tasks(gene_piece_task, pieces, sizeof(struct gene_piece), n);
  return 0;
}

/*
 * Positions low .. high - 1 made directly, the cells cut into `threads`
 * pieces (each counted, then placed, on a thread of its own), at cell_bytes
 * and value_bytes. 0, -1 when the counts disagree with the pointers, -2 out
 * of memory. `visited` (may be NULL) adds the cells passed over.
 */
static int gene_major_direct(const struct gene_major_input *in, int64_t low, int64_t high, int threads,
                             unsigned char *cell_bytes, unsigned char *value_bytes, double *visited) {
  if (high <= low) return 0;
  struct gene_block b = gene_block_of(in, low, high);
  R_xlen_t n_block_genes = b.g_high - b.g_low + 1;
  int n = threads < 1 ? 1 : threads;
  if ((R_xlen_t)n > in->n_cells) n = in->n_cells > 0 ? (int)in->n_cells : 1;
  struct gene_piece *pieces = (struct gene_piece *)malloc((size_t)n * sizeof(struct gene_piece));
  int64_t *cursors = (int64_t *)malloc((size_t)n * (size_t)n_block_genes * sizeof(int64_t));
  if (pieces == NULL || cursors == NULL) {
    free(pieces);
    free(cursors);
    return -2;
  }
  for (int k = 0; k < n; k++) {
    pieces[k].in = in;
    pieces[k].b = &b;
    pieces[k].first = (int64_t)in->n_cells * k / n;
    pieces[k].last = (int64_t)in->n_cells * (k + 1) / n;
    pieces[k].next = cursors + (size_t)k * (size_t)n_block_genes;
    pieces[k].cell_bytes = cell_bytes;
    pieces[k].value_bytes = value_bytes;
  }
  int failed = n == 1 ? 0 : gene_major_pieces(pieces, n, 0);
  if (n == 1) {
    /* One piece: its cursors are the genes' pointers, nothing to count. */
    for (R_xlen_t g = 0; g < n_block_genes; g++) cursors[g] = (int64_t)in->pointer[b.g_low + g];
    gene_major_cells(in, &b, 0, in->n_cells, cursors, cell_bytes, value_bytes, NULL);
  }
  if (visited != NULL) *visited += (double)in->n_cells * (n == 1 ? 1 : 2);
  free(pieces);
  free(cursors);
  return failed;
}

/* The arguments of a gene-major pass, checked on R's thread. */
static struct gene_major_input gene_major_input_of(SEXP i, SEXP p, SEXP x, SEXP map, SEXP origin, SEXP pointers,
                                                   SEXP cells, SEXP index_size, SEXP data_dtype) {
  if (TYPEOF(i) != INTSXP || TYPEOF(p) != INTSXP || TYPEOF(x) != REALSXP || TYPEOF(pointers) != REALSXP ||
      TYPEOF(map) != INTSXP || TYPEOF(origin) != INTSXP || TYPEOF(cells) != INTSXP) {
    Rf_error("The gene-major block needs a dgCMatrix's i, p and x, the row map and origins, the cell order and the pointers as doubles.");
  }
  enum dtype out_type = parse_dtype(data_dtype);
  if (out_type != F4 && out_type != F8) Rf_error("Expression values are stored as <f4 or <f8.");
  struct gene_major_input in;
  in.value_size = dtype_size(out_type);
  in.cell_size = Rf_asInteger(index_size);
  if (in.cell_size != 4 && in.cell_size != 8) Rf_error("Cell indices are 4 or 8 bytes.");
  in.rows = INTEGER(i);
  in.columns = INTEGER(p);
  in.values = REAL(x);
  in.target = INTEGER(map);
  in.source = INTEGER(origin);
  in.pointer = REAL(pointers);
  in.n_genes = XLENGTH(pointers) - 1;
  in.n_cells = XLENGTH(p) - 1;
  /* The cells in the order the output rows take them (0-based columns): the
   * priority order of a store in priority order, none for the matrix's own. */
  in.order = XLENGTH(cells) > 0 ? INTEGER(cells) : NULL;
  if (in.order != NULL && XLENGTH(cells) != in.n_cells) Rf_error("The cell order names %d cells, the matrix %d.", (int)XLENGTH(cells), (int)in.n_cells);
  if (XLENGTH(origin) != in.n_genes) Rf_error("The origins name %d genes, the pointers %d.", (int)XLENGTH(origin), (int)in.n_genes);
  if (in.n_genes < 1) Rf_error("The block lies outside the matrix's values.");
  return in;
}

/* Threads R asks for: at least one. */
static int threads_of(SEXP threads) {
  int n = Rf_asInteger(threads);
  return n == NA_INTEGER || n < 1 ? 1 : n;
}

SEXP cellucid_gene_major_block(SEXP i, SEXP p, SEXP x, SEXP map, SEXP origin, SEXP pointers, SEXP cells,
                               SEXP first, SEXP last, SEXP index_size, SEXP data_dtype, SEXP threads) {
  struct gene_major_input in = gene_major_input_of(i, p, x, map, origin, pointers, cells, index_size, data_dtype);
  int64_t low = (int64_t)Rf_asReal(first);
  int64_t high = (int64_t)Rf_asReal(last);
  if (low < 0 || high < low || high > (int64_t)in.pointer[in.n_genes]) {
    Rf_error("The block lies outside the matrix's values.");
  }
  R_xlen_t n = raw_length((double)(high - low), 1, "A gene-major block");
  SEXP cells_out = PROTECT(Rf_allocVector(RAWSXP, raw_length((double)n, in.cell_size, "A block's cell indices")));
  SEXP values_out = PROTECT(Rf_allocVector(RAWSXP, raw_length((double)n, in.value_size, "A block's values")));
  int failed = gene_major_direct(&in, low, high, threads_of(threads), RAW(cells_out), RAW(values_out), NULL);
  if (failed == -2) Rf_error("Out of memory for a gene-major block.");
  if (failed != 0) Rf_error("The gene-major pointers do not match the matrix's values.");
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


/*
 * The `n` pairs in ascending row (rows are distinct and below 2^`bits`), into
 * `pairs`; `spare` is as long. A short column (a graph's row) by insertion;
 * a long one (a gene) by a least-significant-digit radix sort, a byte of the
 * row a pass: about ten times faster than qsort() on a gene of 10^4 values,
 * which made priority order cost as much as the rest of writing X.
 */
void cellucid_sort_by_row(struct renumbered *pairs, struct renumbered *spare, R_xlen_t n, int bits) {
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
/* What a renumbered pass reads: a dgCMatrix whose rows are renumbered, its columns in an order. */
struct renumbered_input {
  const int *rows, *columns, *renumber, *source, *pointer;
  const double *values;
  R_xlen_t n_rows, n_columns;
  int cell_size, value_size, bits;
};

/*
 * Positions low .. high - 1 of the renumbered column-compressed arrays
 * (struct renumbered_input) at cell_bytes and value_bytes (either may be
 * NULL: that part is not written), position low at their start. `pairs`
 * and `spare` hold the longest column reaching into the block. Returns 0,
 * or -1 when a row or a new row lies outside the rows renumbered.
 */
static int renumbered_fill(const struct renumbered_input *in, int64_t low, int64_t high,
                           unsigned char *cell_bytes, unsigned char *value_bytes,
                           struct renumbered *pairs, struct renumbered *spare) {
  if (high <= low) return 0;
  const int *pointer = in->pointer;
  R_xlen_t a = 0, b = in->n_columns; /* pointer[a] <= low < pointer[b] */
  while (b - a > 1) {
    R_xlen_t m = a + (b - a) / 2;
    if ((int64_t)pointer[m] <= low) a = m; else b = m;
  }
  for (R_xlen_t g = a; g < in->n_columns && (int64_t)pointer[g] < high; g++) {
    R_xlen_t column = in->source == NULL ? g : in->source[g];
    int start = in->columns[column], end = in->columns[column + 1];
    if (end == start) continue;
    for (int k = start; k < end; k++) {
      int row = in->rows[k];
      if (row < 0 || row >= in->n_rows) return -1;
      int renumbered = in->renumber[row];
      if (renumbered < 0 || renumbered >= in->n_rows) return -1;
      pairs[k - start].row = renumbered;
      pairs[k - start].value = in->values[k];
    }
    cellucid_sort_by_row(pairs, spare, end - start, in->bits);
    int64_t base = pointer[g];
    for (int k = 0; k < end - start; k++) {
      int64_t position = base + k;
      if (position < low || position >= high) continue;
      R_xlen_t at = (R_xlen_t)(position - low);
      if (cell_bytes != NULL) {
        if (in->cell_size == 4) put32(cell_bytes + 4 * at, (uint32_t)pairs[k].row);
        else put64(cell_bytes + 8 * at, (uint64_t)pairs[k].row);
      }
      if (value_bytes != NULL) {
        if (in->value_size == 4) put_f4(value_bytes + 4 * at, pairs[k].value);
        else put_f8(value_bytes + 8 * at, pairs[k].value);
      }
    }
  }
  return 0;
}

/* The arguments of a renumbered pass, checked on R's thread. */
static struct renumbered_input renumbered_input_of(SEXP i, SEXP p, SEXP x, SEXP new_row, SEXP column_order, SEXP pointers,
                                                   SEXP index_size, SEXP data_dtype) {
  if (TYPEOF(i) != INTSXP || TYPEOF(p) != INTSXP || TYPEOF(x) != REALSXP || TYPEOF(new_row) != INTSXP ||
      TYPEOF(column_order) != INTSXP || TYPEOF(pointers) != INTSXP) {
    Rf_error("The renumbered block needs a dgCMatrix's i, p and x, the new row of each row, the column order and the output pointers.");
  }
  enum dtype out_type = parse_dtype(data_dtype);
  if (out_type != F4 && out_type != F8) Rf_error("Stored values are <f4 or <f8.");
  struct renumbered_input in;
  in.value_size = dtype_size(out_type);
  in.cell_size = Rf_asInteger(index_size);
  if (in.cell_size != 4 && in.cell_size != 8) Rf_error("Row indices are 4 or 8 bytes.");
  in.rows = INTEGER(i);
  in.columns = INTEGER(p);
  in.values = REAL(x);
  in.renumber = INTEGER(new_row);
  in.n_rows = XLENGTH(new_row);
  in.pointer = INTEGER(pointers);
  in.n_columns = XLENGTH(p) - 1;
  in.source = XLENGTH(column_order) > 0 ? INTEGER(column_order) : NULL;
  if (in.source != NULL && XLENGTH(column_order) != in.n_columns) {
    Rf_error("The column order names %d columns, the matrix %d.", (int)XLENGTH(column_order), (int)in.n_columns);
  }
  if (XLENGTH(pointers) != in.n_columns + 1) Rf_error("The output pointers name %d columns, the matrix %d.", (int)XLENGTH(pointers) - 1, (int)in.n_columns);
  if (in.n_columns < 1 || in.pointer[in.n_columns] != in.columns[in.n_columns]) Rf_error("The block lies outside the matrix's values.");
  for (R_xlen_t g = 0; g < in.n_columns; g++) {
    R_xlen_t column = in.source == NULL ? g : in.source[g];
    if (column < 0 || column >= in.n_columns || in.columns[column + 1] - in.columns[column] != in.pointer[g + 1] - in.pointer[g]) {
      Rf_error("The output pointers do not match the matrix's columns in that order.");
    }
  }
  /* The bits a new row takes: rows are below n_rows. */
  in.bits = 0;
  while (in.bits < 63 && ((int64_t)1 << in.bits) < (int64_t)in.n_rows) in.bits++;
  return in;
}

/* The longest column of `in` reaching into positions low .. high - 1. */
static int longest_column(const struct renumbered_input *in, int64_t low, int64_t high) {
  int longest = 0;
  for (R_xlen_t g = 0; g < in->n_columns; g++) {
    if ((int64_t)in->pointer[g + 1] <= low || (int64_t)in->pointer[g] >= high) continue;
    int length = in->pointer[g + 1] - in->pointer[g];
    if (length > longest) longest = length;
  }
  return longest;
}

SEXP cellucid_renumbered_block(SEXP i, SEXP p, SEXP x, SEXP new_row, SEXP column_order, SEXP pointers,
                               SEXP first, SEXP last, SEXP index_size, SEXP data_dtype, SEXP part, SEXP padded) {
  struct renumbered_input in = renumbered_input_of(i, p, x, new_row, column_order, pointers, index_size, data_dtype);
  int64_t low = (int64_t)Rf_asReal(first);
  int64_t high = (int64_t)Rf_asReal(last);
  if (low < 0 || high < low || high > (int64_t)in.pointer[in.n_columns]) {
    Rf_error("The block lies outside the matrix's values.");
  }
  R_xlen_t n = raw_length((double)(high - low), 1, "A block");
  int which = Rf_asInteger(part);
  if (which < 0 || which > 2) Rf_error("A block's part is 0 (both), 1 (row indices) or 2 (values).");
  R_xlen_t length = which == 0 ? n : xlength_of(Rf_asReal(padded), "A chunk's length");
  if (length < n) Rf_error("A chunk of %.0f items cannot hold %.0f.", (double)length, (double)n);
  SEXP cells_out = PROTECT(which == 2 ? R_NilValue : Rf_allocVector(RAWSXP, raw_length((double)length, in.cell_size, "A block's row indices")));
  SEXP values_out = PROTECT(which == 1 ? R_NilValue : Rf_allocVector(RAWSXP, raw_length((double)length, in.value_size, "A block's values")));
  unsigned char *cell_bytes = which == 2 ? NULL : RAW(cells_out);
  unsigned char *value_bytes = which == 1 ? NULL : RAW(values_out);
  if (cell_bytes != NULL && length > n) memset(cell_bytes + n * in.cell_size, 0, (size_t)((length - n) * in.cell_size));
  if (value_bytes != NULL && length > n) memset(value_bytes + n * in.value_size, 0, (size_t)((length - n) * in.value_size));
  if (n > 0) {
    /* One buffer, as long as the longest column reaching into the block, used
       by each column in turn (R_alloc'd per column, every column's pairs would
       be held until the call returns: 16 bytes a value of the whole block). */
    int longest = longest_column(&in, low, high);
    struct renumbered *pairs = longest > 0 ? (struct renumbered *)R_alloc((size_t)longest, sizeof(struct renumbered)) : NULL;
    struct renumbered *spare = longest >= 64 ? (struct renumbered *)R_alloc((size_t)longest, sizeof(struct renumbered)) : NULL;
    if (renumbered_fill(&in, low, high, cell_bytes, value_bytes, pairs, spare) != 0) {
      Rf_error("A row index or its new row lies outside the %d rows renumbered.", (int)in.n_rows);
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


/*
 * A served store's expression matrices laid out gene-major once, in the
 * background (R/expression-input.R, .laid_out_arrays()): the gene-major
 * arrays of every matrix, made by the passes above over all its values,
 * into memory-mapped scratch files -- the operating system keeps their pages
 * or writes them to disk, never a second copy of the matrix in R's memory --
 * one matrix after another in the order the store planned them (X first),
 * each on every core (its cells, or its columns, in one segment a thread), so
 * R's thread answers the viewer meanwhile. Once laid out, every chunk is a
 * copy of the bytes it covers. Before, a chunk is made from the layout's
 * progress (cellucid_layout_block): the layouts stand while it is made, the
 * values each segment has laid out are copied, and only the rest is made,
 * on every core -- a gene asked halfway through costs half a pass over the
 * cells, shared among the threads, and none near the end. cellucid-python's
 * served store does the same (_served_store.py, GeneMajorLayout).
 *
 * The thread touches no R object: the matrices' i, p and x vectors are read
 * where R holds them (R never moves a vector), kept alive by the queue's
 * external pointer, and every argument is checked on R's thread first. A
 * scratch file is removed as soon as it is mapped (on Windows, when it is
 * closed), so nothing is left behind however R ends.
 */

void cellucid_scratch_init(struct scratch_map *m) {
  m->bytes = NULL;
  m->length = 0;
#ifdef _WIN32
  m->file = INVALID_HANDLE_VALUE;
  m->map = NULL;
#else
  m->fd = -1;
#endif
}

const char *cellucid_map_scratch(struct scratch_map *m, const char *path, size_t length) {
  m->bytes = NULL;
  m->length = length;
#ifdef _WIN32
  m->file = INVALID_HANDLE_VALUE;
  m->map = NULL;
  if (length == 0) return NULL;
  m->file = CreateFileA(path, GENERIC_READ | GENERIC_WRITE, 0, NULL, CREATE_NEW,
                        FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_DELETE_ON_CLOSE, NULL);
  if (m->file == INVALID_HANDLE_VALUE) return "could not create a scratch file";
  m->map = CreateFileMappingA(m->file, NULL, PAGE_READWRITE, (DWORD)((uint64_t)length >> 32), (DWORD)(length & 0xFFFFFFFFu), NULL);
  if (m->map == NULL) return "could not map a scratch file";
  m->bytes = (unsigned char *)MapViewOfFile(m->map, FILE_MAP_ALL_ACCESS, 0, 0, length);
  if (m->bytes == NULL) return "could not map a scratch file";
#else
  m->fd = -1;
  if (length == 0) return NULL;
  m->fd = open(path, O_RDWR | O_CREAT | O_EXCL, 0600);
  if (m->fd < 0) return "could not create a scratch file";
  unlink(path);
  if (ftruncate(m->fd, (off_t)length) != 0) return "could not size a scratch file (is the disk full?)";
  void *mapped = mmap(NULL, length, PROT_READ | PROT_WRITE, MAP_SHARED, m->fd, 0);
  if (mapped == MAP_FAILED) return "could not map a scratch file";
  m->bytes = (unsigned char *)mapped;
#endif
  return NULL;
}

void cellucid_unmap_scratch(struct scratch_map *m) {
#ifdef _WIN32
  if (m->bytes != NULL) UnmapViewOfFile(m->bytes);
  if (m->map != NULL) CloseHandle(m->map);
  if (m->file != INVALID_HANDLE_VALUE) CloseHandle(m->file);
  m->file = INVALID_HANDLE_VALUE;
  m->map = NULL;
#else
  if (m->bytes != NULL) munmap(m->bytes, m->length);
  if (m->fd >= 0) close(m->fd);
  m->fd = -1;
#endif
  m->bytes = NULL;
}

/* Monotonic seconds. On Windows the performance counter, which every
 * Windows has: GetTickCount64() is declared only for Vista's API level, which
 * Rtools 4.0's toolchains (R 4.1) do not target by default, so the 32-bit
 * build failed to link it. */
static double now_seconds(void) {
#ifdef _WIN32
  LARGE_INTEGER count, frequency;
  QueryPerformanceCounter(&count);
  QueryPerformanceFrequency(&frequency);
  return (double)count.QuadPart / (double)frequency.QuadPart;
#else
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + (double)t.tv_nsec / 1e9;
#endif
}

/* A short wait while another thread reaches a step (a pause, a counting's end). */
static void nap(void) {
#ifdef _WIN32
  SwitchToThread();
#else
  struct timespec t = {0, 20000};
  nanosleep(&t, NULL);
#endif
}

/* A layout thread's wait while it stands: longer, so the threads standing
 * leave the cores to the request (each wakes 5,000 times a second, not
 * 50,000). */
static void stand(void) {
#ifdef _WIN32
  Sleep(0);
#else
  struct timespec t = {0, 200000};
  nanosleep(&t, NULL);
#endif
}

enum { LAYOUT_WAITING = 0, LAYOUT_RUNNING = 1, LAYOUT_DONE = 2, LAYOUT_FAILED = 3 };

struct layout_queue;
struct layout_job;

/*
 * One thread's share of a layout: units first .. last - 1 (kind 0: cells in
 * output order; kind 1: output columns). Kind 0: `start` is where its values
 * of each gene begin, `next` its cursors (its counts while counting); the
 * values of gene g at start[g] .. next[g] - 1 are laid out, those of cells
 * first .. at - 1. Kind 1: the columns first .. at - 1 are laid out.
 */
struct layout_segment {
  struct layout_queue *queue;
  struct layout_job *job;
  int64_t first, last;
  int64_t *start, *next;
  int64_t at;      /* release: every unit before it laid out */
  int started;     /* kind 0: its placing began (release) */
  int paused;      /* the pause it stands at (release, with at and next) */
  int done;        /* its units are laid out (release) */
  int failed;
};

struct layout_job {
  int kind; /* 0: a genes x cells dgCMatrix (gene_major_cells), 1: renumbered (renumbered_fill) */
  struct gene_major_input gene_major;
  struct renumbered_input renumbered;
  struct gene_block block; /* kind 0: every value */
  int64_t n_values;
  int64_t units;     /* what progress counts: cells (kind 0) or columns (kind 1) */
  int state;         /* LAYOUT_* (release / acquire) */
  int placing;       /* kind 0: the counting is done, the segments' cursors set (release) */
  int n_segments;
  struct layout_segment *segments;
  int64_t *cursors;  /* the segments' start and next arrays (kind 0) */
  /* Kind 0 in an order: each column's output row (the order's inverse), made
   * while counting, let go on R's thread once the layout is done. */
  int *row_of;
  double started, finished;
  double placing_from; /* kind 0: when its placing began (written before `placing`) */
  /* What requests made before the layout ended (R's thread only): units they
   * passed over themselves (kind 0: cells, both passes of a count; kind 1:
   * values sorted), values they copied from the layout; the seconds they
   * held its placing, and took making its rest (all, and the last). */
  double request_units, copied;
  double held, spent, last;
  struct scratch_map parts[2];
};

struct layout_queue {
  struct layout_job *jobs;
  int n_jobs;
  int threads;
  int stop;
  /* A request's pause (R's thread sets it, 0 none): every running segment
   * stands at its next step until it is lifted, its progress published. */
  int hold;
  int generation;
  /* Units into each segment a layout stands at until lifted (-1 none): how
   * the tests catch a layout halfway. */
  int64_t park;
  int has_thread;
#ifdef _WIN32
  HANDLE thread;
#else
  pthread_t thread;
#endif
};

/*
 * A segment's step at unit `c`: its progress published; stands while a
 * request holds the layouts (its cursors then still, published with the
 * pause) or while parked. 1 when the layouts are stopped.
 */
static int segment_step(struct layout_segment *s, int64_t c) {
  struct layout_queue *queue = s->queue;
  __atomic_store_n(&s->at, c, __ATOMIC_RELEASE);
  for (;;) {
    if (__atomic_load_n(&queue->stop, __ATOMIC_RELAXED)) return 1;
    int hold = __atomic_load_n(&queue->hold, __ATOMIC_ACQUIRE);
    int64_t park = __atomic_load_n(&queue->park, __ATOMIC_ACQUIRE);
    int parked = park >= 0 && c - s->first >= park;
    if (hold == 0 && !parked) return 0;
    if (hold != 0 && __atomic_load_n(&s->paused, __ATOMIC_RELAXED) != hold) __atomic_store_n(&s->paused, hold, __ATOMIC_RELEASE);
    stand();
  }
}

static void count_segment(void *argument) {
  struct layout_segment *s = (struct layout_segment *)argument;
  struct layout_job *job = s->job;
  memset(s->next, 0, (size_t)job->gene_major.n_genes * sizeof(int64_t));
  if (job->row_of != NULL) {
    for (int64_t c = s->first; c < s->last; c++) job->row_of[job->gene_major.order[c]] = (int)c;
  }
  /* Counting does not pause: it is short, and a request waits for its end. */
  int64_t c = s->first;
  while (c < s->last && !__atomic_load_n(&s->queue->stop, __ATOMIC_RELAXED)) {
    int64_t until = c + 65536 < s->last ? c + 65536 : s->last;
    gene_major_cells(&job->gene_major, &job->block, c, until, s->next, NULL, NULL, NULL);
    c = until;
  }
}

static void place_segment(void *argument) {
  struct layout_segment *s = (struct layout_segment *)argument;
  struct layout_job *job = s->job;
  __atomic_store_n(&s->started, 1, __ATOMIC_RELEASE);
  if (job->kind == 0) {
    if (gene_major_cells(&job->gene_major, &job->block, s->first, s->last, s->next, job->parts[0].bytes, job->parts[1].bytes, s)) {
      s->failed = 1;
      return;
    }
  } else {
    /* Column by column: each a step, so a column is copied once laid out. */
    const struct renumbered_input *in = &job->renumbered;
    int longest = longest_column(in, in->pointer[s->first], in->pointer[s->last]);
    struct renumbered *pairs = (struct renumbered *)malloc((size_t)(longest > 0 ? longest : 1) * sizeof(struct renumbered));
    struct renumbered *spare = (struct renumbered *)malloc((size_t)(longest > 0 ? longest : 1) * sizeof(struct renumbered));
    if (pairs == NULL || spare == NULL) s->failed = 1;
    for (int64_t g = s->first; g < s->last && !s->failed; g++) {
      if (segment_step(s, g)) {
        s->failed = 1;
        break;
      }
      int64_t low = in->pointer[g], high = in->pointer[g + 1];
      if (high > low && renumbered_fill(in, low, high, job->parts[0].bytes + (size_t)low * (size_t)in->cell_size,
                                        job->parts[1].bytes + (size_t)low * (size_t)in->value_size, pairs, spare) != 0) {
        s->failed = 1;
      }
    }
    free(pairs);
    free(spare);
    if (s->failed) return;
  }
  __atomic_store_n(&s->at, s->last, __ATOMIC_RELEASE);
  __atomic_store_n(&s->done, 1, __ATOMIC_RELEASE);
}

/*
 * One matrix laid out by its segments, one a thread. Kind 0: each segment's
 * values of every gene counted, its cursors following the previous
 * segment's, then placed; kind 1: each segment's columns renumbered and
 * sorted. Either way each value lands where one pass would put it.
 */
static void run_layout(struct layout_queue *queue, struct layout_job *job) {
  job->started = now_seconds();
  __atomic_store_n(&job->state, LAYOUT_RUNNING, __ATOMIC_RELEASE);
  int failed = 0;
  if (job->kind == 0) {
    run_tasks(count_segment, job->segments, sizeof(struct layout_segment), job->n_segments);
    for (R_xlen_t g = 0; g < job->gene_major.n_genes && !failed; g++) {
      int64_t cursor = (int64_t)job->gene_major.pointer[g];
      for (int k = 0; k < job->n_segments; k++) {
        struct layout_segment *s = &job->segments[k];
        int64_t count = s->next[g];
        s->start[g] = s->next[g] = cursor;
        cursor += count;
      }
      if (cursor != (int64_t)job->gene_major.pointer[g + 1]) failed = 1;
    }
    job->placing_from = now_seconds();
    if (!failed) __atomic_store_n(&job->placing, 1, __ATOMIC_RELEASE);
  }
  if (!failed && !__atomic_load_n(&queue->stop, __ATOMIC_RELAXED)) {
    run_tasks(place_segment, job->segments, sizeof(struct layout_segment), job->n_segments);
  }
  for (int k = 0; k < job->n_segments; k++) failed = failed || job->segments[k].failed || !job->segments[k].done;
  job->finished = now_seconds();
  if (__atomic_load_n(&queue->stop, __ATOMIC_RELAXED)) failed = 1;
  __atomic_store_n(&job->state, failed ? LAYOUT_FAILED : LAYOUT_DONE, __ATOMIC_RELEASE);
}

#ifdef _WIN32
static DWORD WINAPI layout_worker(LPVOID argument) {
#else
static void *layout_worker(void *argument) {
#endif
  struct layout_queue *queue = (struct layout_queue *)argument;
  for (int j = 0; j < queue->n_jobs; j++) {
    if (__atomic_load_n(&queue->stop, __ATOMIC_RELAXED)) break;
    run_layout(queue, &queue->jobs[j]);
  }
#ifdef _WIN32
  return 0;
#else
  return NULL;
#endif
}

static void free_layout_queue(struct layout_queue *queue) {
  if (queue == NULL) return;
  __atomic_store_n(&queue->stop, 1, __ATOMIC_RELAXED);
  if (queue->has_thread) {
#ifdef _WIN32
    WaitForSingleObject(queue->thread, INFINITE);
    CloseHandle(queue->thread);
#else
    pthread_join(queue->thread, NULL);
#endif
    queue->has_thread = 0;
  }
  for (int j = 0; j < queue->n_jobs; j++) {
    cellucid_unmap_scratch(&queue->jobs[j].parts[0]);
    cellucid_unmap_scratch(&queue->jobs[j].parts[1]);
    free(queue->jobs[j].segments);
    free(queue->jobs[j].cursors);
    free(queue->jobs[j].row_of);
  }
  free(queue->jobs);
  free(queue);
}

static void layout_queue_finalizer(SEXP pointer) {
  struct layout_queue *queue = (struct layout_queue *)R_ExternalPtrAddr(pointer);
  R_ClearExternalPtr(pointer);
  free_layout_queue(queue);
}

static struct layout_queue *queue_of(SEXP pointer) {
  if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == NULL) Rf_error("The served store's layouts were stopped.");
  return (struct layout_queue *)R_ExternalPtrAddr(pointer);
}

/* A job's segments: `n` shares of its units (kind 1: of its values, whole columns each). */
static int layout_segments(struct layout_queue *queue, struct layout_job *job, int n) {
  if (job->kind == 0 && (int64_t)n > job->units) n = job->units > 0 ? (int)job->units : 1;
  if (job->kind == 1 && (int64_t)n > job->units) n = (int)job->units;
  job->n_segments = n;
  job->segments = (struct layout_segment *)calloc((size_t)n, sizeof(struct layout_segment));
  if (job->segments == NULL) return 0;
  if (job->kind == 0) {
    size_t genes = (size_t)job->gene_major.n_genes;
    job->cursors = (int64_t *)malloc(2 * (size_t)n * genes * sizeof(int64_t));
    if (job->cursors == NULL) return 0;
    if (job->gene_major.order != NULL) {
      job->row_of = (int *)malloc((size_t)job->gene_major.n_cells * sizeof(int));
      if (job->row_of == NULL) return 0;
    }
  }
  const int *pointer = job->kind == 1 ? job->renumbered.pointer : NULL;
  int64_t column = 0;
  for (int k = 0; k < n; k++) {
    struct layout_segment *s = &job->segments[k];
    s->queue = queue;
    s->job = job;
    if (job->kind == 0) {
      s->first = job->units * k / n;
      s->last = job->units * (k + 1) / n;
      s->start = job->cursors + (size_t)(2 * k) * (size_t)job->gene_major.n_genes;
      s->next = job->cursors + (size_t)(2 * k + 1) * (size_t)job->gene_major.n_genes;
    } else {
      /* Whole columns, about a share of the values each. */
      int64_t until = job->n_values * (k + 1) / n;
      s->first = column;
      while (column < job->units && (k == n - 1 || (int64_t)pointer[column + 1] <= until)) column++;
      s->last = column;
    }
    s->at = s->first;
  }
  return 1;
}

/*
 * Start laying out `jobs` (a list; each a list: kind 0 -- i, p, x, map,
 * origin, pointers (doubles), cells, index_size, dtype -- or kind 1 -- i, p,
 * x, new_row, column_order, pointers (integers), index_size, dtype), their
 * scratch files at `paths` (two a job: cell indices, values), each on
 * `threads` threads, one job after another. `park` (-1: none) stands every
 * segment at that many units in until cellucid_layouts_unpark() (tests).
 * Returns the queue, an external pointer that keeps `jobs` alive and, when
 * collected, stops the threads at their next step and lets the files go.
 */
SEXP cellucid_layouts_start(SEXP jobs, SEXP paths, SEXP threads, SEXP park) {
  if (TYPEOF(jobs) != VECSXP || TYPEOF(paths) != STRSXP || XLENGTH(paths) != 2 * XLENGTH(jobs)) {
    Rf_error("The layouts need their jobs and two scratch paths a job.");
  }
  int n_jobs = (int)XLENGTH(jobs);
  struct layout_queue *queue = (struct layout_queue *)calloc(1, sizeof(struct layout_queue));
  if (queue == NULL) Rf_error("Out of memory for the served store's layouts.");
  queue->jobs = (struct layout_job *)calloc((size_t)(n_jobs > 0 ? n_jobs : 1), sizeof(struct layout_job));
  if (queue->jobs == NULL) {
    free(queue);
    Rf_error("Out of memory for the served store's layouts.");
  }
  for (int j = 0; j < n_jobs; j++) {
    for (int part = 0; part < 2; part++) cellucid_scratch_init(&queue->jobs[j].parts[part]);
  }
  queue->n_jobs = n_jobs;
  queue->threads = threads_of(threads);
  double parked = Rf_asReal(park);
  queue->park = ISNAN(parked) || parked < 0 ? -1 : (int64_t)parked;
  /* Owned by the external pointer from here: an error below frees it. */
  SEXP pointer = PROTECT(R_MakeExternalPtr(queue, R_NilValue, jobs));
  R_RegisterCFinalizerEx(pointer, layout_queue_finalizer, TRUE);
  for (int j = 0; j < n_jobs; j++) {
    SEXP spec = VECTOR_ELT(jobs, j);
    struct layout_job *job = &queue->jobs[j];
    job->kind = Rf_asInteger(VECTOR_ELT(spec, 0));
    int cell_size, value_size;
    if (job->kind == 0) {
      job->gene_major = gene_major_input_of(VECTOR_ELT(spec, 1), VECTOR_ELT(spec, 2), VECTOR_ELT(spec, 3), VECTOR_ELT(spec, 4),
                                            VECTOR_ELT(spec, 5), VECTOR_ELT(spec, 6), VECTOR_ELT(spec, 7), VECTOR_ELT(spec, 8),
                                            VECTOR_ELT(spec, 9));
      job->n_values = (int64_t)job->gene_major.pointer[job->gene_major.n_genes];
      job->units = job->gene_major.n_cells;
      cell_size = job->gene_major.cell_size;
      value_size = job->gene_major.value_size;
      if (job->n_values > 0) job->block = gene_block_of(&job->gene_major, 0, job->n_values);
    } else if (job->kind == 1) {
      job->renumbered = renumbered_input_of(VECTOR_ELT(spec, 1), VECTOR_ELT(spec, 2), VECTOR_ELT(spec, 3), VECTOR_ELT(spec, 4),
                                            VECTOR_ELT(spec, 5), VECTOR_ELT(spec, 6), VECTOR_ELT(spec, 7), VECTOR_ELT(spec, 8));
      job->n_values = (int64_t)job->renumbered.pointer[job->renumbered.n_columns];
      job->units = job->renumbered.n_columns;
      cell_size = job->renumbered.cell_size;
      value_size = job->renumbered.value_size;
      for (R_xlen_t r = 0; r < job->renumbered.n_rows; r++) {
        if (job->renumbered.renumber[r] < 0 || job->renumbered.renumber[r] >= job->renumbered.n_rows) {
          Rf_error("A new row lies outside the %d rows.", (int)job->renumbered.n_rows);
        }
      }
    } else {
      Rf_error("A layout is of kind 0 or 1.");
    }
    if (!layout_segments(queue, job, queue->threads)) Rf_error("Out of memory for the served store's layouts.");
    size_t lengths[2] = {(size_t)job->n_values * (size_t)cell_size, (size_t)job->n_values * (size_t)value_size};
    for (int part = 0; part < 2; part++) {
      const char *failed = cellucid_map_scratch(&job->parts[part], Rf_translateChar(STRING_ELT(paths, 2 * j + part)), lengths[part]);
      if (failed != NULL) Rf_error("The served store's layout %s: %s.", failed, Rf_translateChar(STRING_ELT(paths, 2 * j + part)));
    }
  }
  if (n_jobs > 0) {
#ifdef _WIN32
    queue->thread = CreateThread(NULL, 0, layout_worker, queue, 0, NULL);
    if (queue->thread == NULL) Rf_error("Could not start the served store's layout thread.");
#else
    if (pthread_create(&queue->thread, NULL, layout_worker, queue) != 0) Rf_error("Could not start the served store's layout thread.");
#endif
    queue->has_thread = 1;
  }
  UNPROTECT(1);
  return pointer;
}

/* Stop the layouts now (a viewer stopped): the threads end at their next step, the files are let go. */
SEXP cellucid_layouts_stop(SEXP pointer) {
  if (TYPEOF(pointer) == EXTPTRSXP) layout_queue_finalizer(pointer);
  return R_NilValue;
}

/* Let parked layouts go on (cellucid_layouts_start(park =)). */
SEXP cellucid_layouts_unpark(SEXP pointer) {
  __atomic_store_n(&queue_of(pointer)->park, (int64_t)-1, __ATOMIC_RELEASE);
  return R_NilValue;
}

/* A job laid out: its order's inverse let go (R's thread; no request needs it again). */
static void release_done(struct layout_job *job) {
  if (job->row_of != NULL && __atomic_load_n(&job->state, __ATOMIC_ACQUIRE) == LAYOUT_DONE) {
    free(job->row_of);
    job->row_of = NULL;
  }
}

/* Units of a job laid out so far (relaxed: a progress report). */
static int64_t layout_units_done(struct layout_job *job) {
  int64_t done = 0;
  for (int k = 0; k < job->n_segments; k++) {
    struct layout_segment *s = &job->segments[k];
    done += __atomic_load_n(&s->at, __ATOMIC_RELAXED) - s->first;
  }
  return done;
}

/*
 * One job's state: c(state, units done, units, seconds running or taken,
 * units requests passed over themselves, values they copied from the layout).
 * state: 0 waiting, 1 running, 2 done, 3 failed (or stopped).
 */
SEXP cellucid_layout_state(SEXP pointer, SEXP index) {
  struct layout_queue *queue = queue_of(pointer);
  int j = Rf_asInteger(index);
  if (j < 0 || j >= queue->n_jobs) Rf_error("No layout %d.", j);
  struct layout_job *job = &queue->jobs[j];
  int state = __atomic_load_n(&job->state, __ATOMIC_ACQUIRE);
  release_done(job);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 6));
  REAL(out)[0] = (double)state;
  REAL(out)[1] = state == LAYOUT_DONE ? (double)job->units : (double)layout_units_done(job);
  REAL(out)[2] = (double)job->units;
  REAL(out)[3] = state == LAYOUT_WAITING ? 0.0 : state == LAYOUT_RUNNING ? now_seconds() - job->started : job->finished - job->started;
  REAL(out)[4] = job->request_units;
  REAL(out)[5] = job->copied;
  UNPROTECT(1);
  return out;
}

/* Copy positions from .. to - 1 of a job's layout into a block starting at position `low`. */
static void copy_laid_out(struct layout_job *job, int64_t from, int64_t to, int64_t low, unsigned char *cell_bytes, unsigned char *value_bytes) {
  if (to <= from) return;
  size_t cell_size = (size_t)(job->kind == 0 ? job->gene_major.cell_size : job->renumbered.cell_size);
  size_t value_size = (size_t)(job->kind == 0 ? job->gene_major.value_size : job->renumbered.value_size);
  memcpy(cell_bytes + (size_t)(from - low) * cell_size, job->parts[0].bytes + (size_t)from * cell_size, (size_t)(to - from) * cell_size);
  memcpy(value_bytes + (size_t)(from - low) * value_size, job->parts[1].bytes + (size_t)from * value_size, (size_t)(to - from) * value_size);
  job->copied += (double)(to - from);
}

/* A value of the block found by a scan: its gene (in the block), output row and value. */
struct scan_hit {
  int64_t row;
  double value;
  int gene;
};

/* One thread's share of a scan: matrix columns from .. to - 1. */
struct scan_piece {
  const struct layout_job *job;
  const struct gene_block *b;
  const int64_t *ats;
  R_xlen_t from, to;
  struct scan_hit *hits;
  int64_t n_hits, room;
  int failed;
};

/* The segment of a job holding output row (unit) `c`. */
static int segment_of(const struct layout_job *job, int64_t c) {
  int k = (int)((double)c * job->n_segments / (double)job->units);
  if (k >= job->n_segments) k = job->n_segments - 1;
  while (k > 0 && job->segments[k].first > c) k--;
  while (k < job->n_segments - 1 && job->segments[k].last <= c) k++;
  return k;
}

static void scan_task(void *argument) {
  struct scan_piece *piece = (struct scan_piece *)argument;
  const struct layout_job *job = piece->job;
  const struct gene_major_input *in = &job->gene_major;
  const struct gene_block *b = piece->b;
  const int *rows = in->rows, *columns = in->columns;
  const unsigned span = (unsigned)(b->r_high - b->r_low);
  /* Every value in turn, the matrix's own order: one comparison each, its
   * branch nearly always the same (a block's rows are a few of the genes'),
   * where a binary search a column mispredicts at every step. */
  int64_t k = columns[piece->from], end_all = columns[piece->to];
  R_xlen_t column = piece->from;
  const int r_low = b->r_low;
  for (; k < end_all; k++) {
    /* Eight values at a time with none in the block's rows passed over
     * together (the comparisons in one vector instruction). */
    while (k + 8 <= end_all) {
      unsigned any = 0;
      for (int j = 0; j < 8; j++) any |= (unsigned)(rows[k + j] - r_low) <= span;
      if (any) break;
      k += 8;
    }
    if (k >= end_all) break;
    if ((unsigned)(rows[k] - r_low) > span) continue;
    while (columns[column + 1] <= k) column++;
    {
      int gene = in->target[rows[k]];
      if (gene < b->g_low || gene > b->g_high) continue;
      int64_t c = job->row_of == NULL ? (int64_t)column : (int64_t)job->row_of[column];
      /* Laid out already (copied): not made again. */
      if (c < piece->ats[segment_of(job, c)]) continue;
      if (piece->n_hits == piece->room) {
        int64_t room = piece->room > 0 ? 2 * piece->room : 4096;
        struct scan_hit *grown = (struct scan_hit *)realloc(piece->hits, (size_t)room * sizeof(struct scan_hit));
        if (grown == NULL) {
          piece->failed = 1;
          return;
        }
        piece->hits = grown;
        piece->room = room;
      }
      struct scan_hit *hit = &piece->hits[piece->n_hits++];
      hit->row = c;
      hit->value = in->values[k];
      hit->gene = (int)(gene - b->g_low);
    }
  }
}

/* Placing the genes a scan found, each gene taken by the next free thread. */
struct scan_place {
  const struct layout_job *job;
  const struct gene_block *b;
  const int64_t *per_gene;
  struct renumbered *pairs, *spare;
  int64_t *cursors;
  unsigned char *cell_bytes, *value_bytes;
  int bits, n_genes;
  int next; /* the next gene to take (atomic) */
};

static void scan_place_task(void *argument) {
  struct scan_place *place = *(struct scan_place **)argument;
  const struct gene_major_input *in = &place->job->gene_major;
  const struct gene_block *b = place->b;
  for (;;) {
    int g = __atomic_fetch_add(&place->next, 1, __ATOMIC_RELAXED);
    if (g >= place->n_genes) return;
    struct renumbered *gene_pairs = place->pairs + place->per_gene[g];
    int64_t count = place->per_gene[g + 1] - place->per_gene[g];
    cellucid_sort_by_row(gene_pairs, place->spare + place->per_gene[g], (R_xlen_t)count, place->bits);
    for (int64_t h = 0; h < count; h++) {
      int segment = segment_of(place->job, gene_pairs[h].row);
      int64_t position = place->cursors[(size_t)segment * (size_t)place->n_genes + (size_t)g]++;
      if (position < b->low || position >= b->high) continue;
      R_xlen_t at = (R_xlen_t)(position - b->low);
      if (in->cell_size == 4) put32(place->cell_bytes + 4 * at, (uint32_t)gene_pairs[h].row);
      else put64(place->cell_bytes + 8 * at, (uint64_t)gene_pairs[h].row);
      if (in->value_size == 4) put_f4(place->value_bytes + 4 * at, gene_pairs[h].value);
      else put_f8(place->value_bytes + 8 * at, gene_pairs[h].value);
    }
  }
}

/*
 * The values of the block's genes in the cells not laid out yet (each
 * segment's cells ats[k] .. last - 1), placed at each segment's cursors of
 * those genes (`cursors`, a segment after another, the block's genes each):
 * one scan of the matrix's values in its own column order -- memory read in
 * sequence, on `threads` threads -- keeps the block's genes' values of those
 * cells; only they are ordered, by gene, then output row (each gene's few
 * sorted by a radix sort), and placed. A pass in output order instead visits
 * the columns at random (two cache misses a cell: 0.2-0.3 s for 10M cells on
 * ten threads). 0, or -1 out of memory.
 */
static int gene_major_scan(const struct layout_job *job, const struct gene_block *b, const int64_t *ats, int64_t *cursors,
                           unsigned char *cell_bytes, unsigned char *value_bytes, int threads) {
  const struct gene_major_input *in = &job->gene_major;
  R_xlen_t n_block_genes = b->g_high - b->g_low + 1;
  int n = threads < 1 ? 1 : threads;
  if ((R_xlen_t)n > in->n_cells) n = (int)in->n_cells;
  struct scan_piece *pieces = (struct scan_piece *)calloc((size_t)n, sizeof(struct scan_piece));
  if (pieces == NULL) return -1;
  /* Columns of about equal values each. */
  int64_t n_values = in->columns[in->n_cells];
  R_xlen_t column = 0;
  for (int k = 0; k < n; k++) {
    pieces[k].job = job;
    pieces[k].b = b;
    pieces[k].ats = ats;
    pieces[k].from = column;
    int64_t until = n_values * (k + 1) / n;
    while (column < in->n_cells && (k == n - 1 || (int64_t)in->columns[column + 1] <= until)) column++;
    pieces[k].to = column;
  }
  run_tasks(scan_task, pieces, sizeof(struct scan_piece), n);
  int failed = 0;
  int64_t total = 0;
  for (int k = 0; k < n; k++) {
    failed = failed || pieces[k].failed;
    total += pieces[k].n_hits;
  }
  int64_t *per_gene = (int64_t *)calloc((size_t)n_block_genes + 1, sizeof(int64_t));
  struct renumbered *pairs = (struct renumbered *)malloc((size_t)(total > 0 ? total : 1) * sizeof(struct renumbered));
  struct renumbered *spare = (struct renumbered *)malloc((size_t)(total > 0 ? total : 1) * sizeof(struct renumbered));
  if (failed || per_gene == NULL || pairs == NULL || spare == NULL) failed = 1;
  if (!failed) {
    /* By gene (a counting sort), then each gene by output row. */
    for (int k = 0; k < n; k++) {
      for (int64_t h = 0; h < pieces[k].n_hits; h++) per_gene[pieces[k].hits[h].gene + 1]++;
    }
    for (R_xlen_t g = 0; g < n_block_genes; g++) per_gene[g + 1] += per_gene[g];
    int64_t *fill = spare == NULL ? NULL : (int64_t *)malloc((size_t)n_block_genes * sizeof(int64_t));
    if (fill == NULL) {
      failed = 1;
    } else {
      memcpy(fill, per_gene, (size_t)n_block_genes * sizeof(int64_t));
      for (int k = 0; k < n; k++) {
        for (int64_t h = 0; h < pieces[k].n_hits; h++) {
          const struct scan_hit *hit = &pieces[k].hits[h];
          struct renumbered *pair = &pairs[fill[hit->gene]++];
          pair->row = hit->row;
          pair->value = hit->value;
        }
      }
      free(fill);
      int bits = 0;
      while (bits < 63 && ((int64_t)1 << bits) < (int64_t)in->n_cells) bits++;
      /* Each gene sorted and placed on its own, the genes shared among the threads. */
      struct scan_place place = {job, b, per_gene, pairs, spare, cursors, cell_bytes, value_bytes, bits, (int)n_block_genes, 0};
      int n_placers = n < (int)n_block_genes ? n : (int)n_block_genes;
      struct scan_place **placers = (struct scan_place **)malloc((size_t)n_placers * sizeof(struct scan_place *));
      if (placers == NULL) {
        failed = 1;
      } else {
        for (int k = 0; k < n_placers; k++) placers[k] = &place;
        run_tasks(scan_place_task, placers, sizeof(struct scan_place *), n_placers);
        free(placers);
      }
    }
  }
  for (int k = 0; k < n; k++) free(pieces[k].hits);
  free(pieces);
  free(per_gene);
  free(pairs);
  free(spare);
  return failed ? -1 : 0;
}

/*
 * Kind 0, laid out in part: the layouts held (every running segment stands
 * at its next step, its cursors still), each segment's values of the
 * block's genes laid out so far copied, and only those of its cells not laid
 * out yet made (gene_major_scan), on every core while the layouts' threads
 * stand, so the request has the machine (a few milliseconds: back-to-back
 * requests still leave the layout the time between them, and a burst waits
 * for it once that is sooner, layout_sooner). 0, -2 out of memory, -3 when
 * the layout failed.
 */
static int gene_major_shared(struct layout_queue *queue, struct layout_job *job, int64_t low, int64_t high,
                             unsigned char *cell_bytes, unsigned char *value_bytes) {
  const struct gene_major_input *in = &job->gene_major;
  struct gene_block b = gene_block_of(in, low, high);
  R_xlen_t n_block_genes = b.g_high - b.g_low + 1;
  int n = job->n_segments;
  int64_t *ats = (int64_t *)malloc((size_t)n * sizeof(int64_t));
  int64_t *cursors = (int64_t *)malloc((size_t)n * (size_t)n_block_genes * sizeof(int64_t));
  if (ats == NULL || cursors == NULL) {
    free(ats);
    free(cursors);
    return -2;
  }
  int failed = 0;
  int generation = ++queue->generation;
  if (generation == 0) generation = ++queue->generation;
  __atomic_store_n(&queue->hold, generation, __ATOMIC_RELEASE);
  /* The counting ends first (short, on every thread); then each segment
   * stands at a step, or has ended, or has not begun placing. */
  while (__atomic_load_n(&job->state, __ATOMIC_ACQUIRE) == LAYOUT_RUNNING && !__atomic_load_n(&job->placing, __ATOMIC_ACQUIRE)) nap();
  if (!__atomic_load_n(&job->placing, __ATOMIC_ACQUIRE) || __atomic_load_n(&job->state, __ATOMIC_ACQUIRE) == LAYOUT_FAILED) {
    __atomic_store_n(&queue->hold, 0, __ATOMIC_RELEASE);
    free(ats);
    free(cursors);
    return -3;
  }
  double holding = now_seconds();
  int n_pieces = 0;
  for (int k = 0; k < n; k++) {
    struct layout_segment *s = &job->segments[k];
    int64_t at = s->first;
    const int64_t *from = s->start;
    for (;;) {
      if (__atomic_load_n(&s->done, __ATOMIC_ACQUIRE)) {
        at = s->last;
        from = s->next;
        break;
      }
      if (!__atomic_load_n(&s->started, __ATOMIC_ACQUIRE)) break;
      if (__atomic_load_n(&s->paused, __ATOMIC_ACQUIRE) == generation) {
        at = __atomic_load_n(&s->at, __ATOMIC_RELAXED);
        from = s->next;
        break;
      }
      nap();
    }
    int64_t *next = cursors + (size_t)k * (size_t)n_block_genes;
    for (R_xlen_t g = 0; g < n_block_genes; g++) {
      next[g] = from[b.g_low + g];
      int64_t begin = s->start[b.g_low + g];
      copy_laid_out(job, begin > low ? begin : low, next[g] < high ? next[g] : high, low, cell_bytes, value_bytes);
    }
    ats[k] = at;
    if (at < s->last) {
      n_pieces++;
      job->request_units += (double)(s->last - at);
    }
  }
  if (n_pieces > 0 && gene_major_scan(job, &b, ats, cursors, cell_bytes, value_bytes, queue->threads) != 0) failed = -2;
  __atomic_store_n(&queue->hold, 0, __ATOMIC_RELEASE);
  job->held += now_seconds() - holding;
  free(ats);
  free(cursors);
  return failed;
}

/*
 * Kind 1, laid out in part: each column of the block a segment has laid out
 * copied, every other one renumbered and sorted here (its own values only).
 * 0, -1 when a row lies outside the rows renumbered, -2 out of memory.
 */
static int renumbered_shared(struct layout_job *job, int64_t low, int64_t high, unsigned char *cell_bytes, unsigned char *value_bytes) {
  const struct renumbered_input *in = &job->renumbered;
  int longest = longest_column(in, low, high);
  struct renumbered *pairs = (struct renumbered *)malloc((size_t)(longest > 0 ? longest : 1) * sizeof(struct renumbered));
  struct renumbered *spare = (struct renumbered *)malloc((size_t)(longest > 0 ? longest : 1) * sizeof(struct renumbered));
  int failed = pairs == NULL || spare == NULL ? -2 : 0;
  int k = 0;
  for (R_xlen_t g = 0; g < in->n_columns && !failed; g++) {
    int64_t from = in->pointer[g] > low ? in->pointer[g] : low;
    int64_t to = in->pointer[g + 1] < high ? in->pointer[g + 1] : high;
    if (to <= from) continue;
    while (k < job->n_segments - 1 && job->segments[k].last <= g) k++;
    if (g < __atomic_load_n(&job->segments[k].at, __ATOMIC_ACQUIRE)) {
      copy_laid_out(job, from, to, low, cell_bytes, value_bytes);
      continue;
    }
    if (renumbered_fill(in, from, to, cell_bytes + (size_t)(from - low) * (size_t)in->cell_size,
                        value_bytes + (size_t)(from - low) * (size_t)in->value_size, pairs, spare) != 0) {
      failed = -1;
    }
    job->request_units += (double)(in->pointer[g + 1] - in->pointer[g]);
  }
  free(pairs);
  free(spare);
  return failed;
}

/* The job being laid out now, if any (R's thread). */
static struct layout_job *running_job(struct layout_queue *queue) {
  for (int j = 0; j < queue->n_jobs; j++) {
    if (__atomic_load_n(&queue->jobs[j].state, __ATOMIC_ACQUIRE) == LAYOUT_RUNNING) return &queue->jobs[j];
  }
  return NULL;
}

/*
 * Whether a kind 0 layout running will end sooner than requests make its
 * rest again: at the pace it has placed while not held, its seconds left are
 * no more than the last request took, or than every request made from its
 * progress took together. The first keeps a gene asked near the end from
 * passing over the cells the layout is about to finish; the second bounds a
 * burst of genes (a hundred asked at once) to twice what waiting for the
 * layout costs, where each would otherwise pass over the cells left again.
 */
static int layout_sooner(struct layout_queue *queue, struct layout_job *job) {
  /* Parked (tests), it ends only when let go. */
  if (!__atomic_load_n(&job->placing, __ATOMIC_ACQUIRE) || __atomic_load_n(&queue->park, __ATOMIC_ACQUIRE) >= 0) return 0;
  int64_t done = layout_units_done(job);
  double working = now_seconds() - job->placing_from - job->held;
  if (done <= 0 || working <= 0) return 0;
  double left = working / (double)done * (double)(job->units - done);
  return left <= job->last || left <= job->spent;
}

/*
 * Positions first .. last - 1 of a job's gene-major arrays, as
 * list(cell index bytes, value bytes), whatever its layout's state: copied
 * once laid out; laid out in part, its progress copied and only the rest
 * made (gene_major_shared, renumbered_shared), unless the layout ends
 * sooner (layout_sooner); not begun, made directly on the layouts' threads'
 * count (the running job held meanwhile).
 */
SEXP cellucid_layout_block(SEXP pointer, SEXP index, SEXP first, SEXP last) {
  struct layout_queue *queue = queue_of(pointer);
  int j = Rf_asInteger(index);
  if (j < 0 || j >= queue->n_jobs) Rf_error("No layout %d.", j);
  struct layout_job *job = &queue->jobs[j];
  int64_t low = (int64_t)Rf_asReal(first);
  int64_t high = (int64_t)Rf_asReal(last);
  if (low < 0 || high < low || high > job->n_values) Rf_error("The block lies outside the layout.");
  int cell_size = job->kind == 0 ? job->gene_major.cell_size : job->renumbered.cell_size;
  int value_size = job->kind == 0 ? job->gene_major.value_size : job->renumbered.value_size;
  /* Allocated first: nothing below may leave the layouts held. */
  SEXP cells_out = PROTECT(Rf_allocVector(RAWSXP, raw_length((double)(high - low), cell_size, "A block's cell indices")));
  SEXP values_out = PROTECT(Rf_allocVector(RAWSXP, raw_length((double)(high - low), value_size, "A block's values")));
  int failed = 0;
  if (high > low) {
    int state = __atomic_load_n(&job->state, __ATOMIC_ACQUIRE);
    if (state == LAYOUT_FAILED) Rf_error("The layout failed.");
    if (state == LAYOUT_DONE) {
      release_done(job);
      copy_laid_out(job, low, high, low, RAW(cells_out), RAW(values_out));
    } else if (job->kind == 1) {
      failed = renumbered_shared(job, low, high, RAW(cells_out), RAW(values_out));
    } else if (state == LAYOUT_RUNNING && layout_sooner(queue, job)) {
      /* The layout ends sooner than its rest would be made again: waited for. */
      while ((state = __atomic_load_n(&job->state, __ATOMIC_ACQUIRE)) == LAYOUT_RUNNING) nap();
      if (state != LAYOUT_DONE) Rf_error("The layout failed.");
      copy_laid_out(job, low, high, low, RAW(cells_out), RAW(values_out));
    } else {
      double began = now_seconds();
      if (state == LAYOUT_RUNNING) {
        failed = gene_major_shared(queue, job, low, high, RAW(cells_out), RAW(values_out));
      } else {
        struct layout_job *running = running_job(queue);
        int generation = ++queue->generation;
        if (generation == 0) generation = ++queue->generation;
        __atomic_store_n(&queue->hold, generation, __ATOMIC_RELEASE);
        failed = gene_major_direct(&job->gene_major, low, high, queue->threads, RAW(cells_out), RAW(values_out), &job->request_units);
        __atomic_store_n(&queue->hold, 0, __ATOMIC_RELEASE);
        if (running != NULL && running->kind == 0 && __atomic_load_n(&running->placing, __ATOMIC_ACQUIRE)) {
          running->held += now_seconds() - began;
        }
      }
      double took = now_seconds() - began;
      job->spent += took;
      job->last = took;
    }
  }
  if (failed == -2) Rf_error("Out of memory for a gene-major block.");
  if (failed == -3) Rf_error("The layout failed.");
  if (failed != 0) Rf_error("The gene-major arrays do not match the matrix's values.");
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, cells_out);
  SET_VECTOR_ELT(out, 1, values_out);
  UNPROTECT(3);
  return out;
}

/*
 * Bytes offset .. offset + count - 1 of part `part` (1: cell indices, 2:
 * values) of a job laid out, then zeros up to `total` bytes: a chunk, or the
 * range of one a request asks for.
 */
SEXP cellucid_layout_slice(SEXP pointer, SEXP index, SEXP part, SEXP offset, SEXP count, SEXP total) {
  struct layout_queue *queue = queue_of(pointer);
  int j = Rf_asInteger(index);
  int which = Rf_asInteger(part);
  if (j < 0 || j >= queue->n_jobs || (which != 1 && which != 2)) Rf_error("No layout %d part %d.", j, which);
  struct layout_job *job = &queue->jobs[j];
  if (__atomic_load_n(&job->state, __ATOMIC_ACQUIRE) != LAYOUT_DONE) Rf_error("The layout is not done.");
  R_xlen_t from = xlength_of(Rf_asReal(offset), "A slice's offset");
  R_xlen_t n = xlength_of(Rf_asReal(count), "A slice's length");
  R_xlen_t size = xlength_of(Rf_asReal(total), "A slice's padded length");
  struct scratch_map *map = &job->parts[which - 1];
  if (n > size || (double)from + (double)n > (double)map->length) Rf_error("The slice lies outside the layout.");
  SEXP out = PROTECT(Rf_allocVector(RAWSXP, size));
  if (n > 0) memcpy(RAW(out), map->bytes + from, (size_t)n);
  if (size > n) memset(RAW(out) + n, 0, (size_t)(size - n));
  UNPROTECT(1);
  return out;
}
