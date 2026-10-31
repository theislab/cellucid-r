/*
 * The gene-major arrays of an expression matrix held outside R's memory
 * (BPCells, a DelayedArray backend such as HDF5Array; R/backed-matrix.R),
 * made from blocks of it read one at a time, never from the whole matrix:
 * nothing here holds more than one block of the matrix and the part of the
 * output being made.
 *
 * Each block arrives as compressed arrays (`index`, `pointers`, `values`),
 * its major lines its columns (cells, col0 ..) or its rows (genes) -- whole
 * cells, or whole genes, whichever the backend stores together -- and a
 * vector naming each of its rows (a source row, or an output gene; a block
 * may hold only the rows a pass needs). The passes:
 *
 *   1. counting (struct tally): each gene's stored values and their extent,
 *      once, before any order is known -- the store's metadata (pointers,
 *      dtype) and the refusal of a value that cannot be stored;
 *   2. prefix counts (struct prefix_tally): each gene's values in the cells
 *      before each prefix level, the cells in priority order;
 *   3. placing (struct target): the values of a run of genes at their
 *      genes' cursors, each value's cell its output row; then each gene put
 *      in ascending row (a radix sort, src/chunks.c). The genes a pass places
 *      are a share of the machine's memory (R/backed-matrix.R), so a matrix
 *      of any size is placed in as many passes as it takes.
 *
 * A target holds its genes' positions in memory (a written store's block,
 * cut into chunks) or in memory-mapped scratch files (a served store's
 * whole matrix, which the operating system keeps or writes to disk). Values
 * are stored as the store stores them: cell indices of `index_size` bytes,
 * values in "<f4" or "<f8", little-endian. Positions and counts are 64-bit
 * (doubles in R): a matrix may hold more than 2^31 - 1 values.
 */

#include "gene_major.h"
#include <stdlib.h>

#ifndef _WIN32
#include <pthread.h>
#endif

/* --- Arguments --------------------------------------------------------------- */

/* A block's arrays, checked: integer minor indices, its pointers and values
 * as doubles; `by_rows` 1 when its major lines are rows. `rows` names each
 * row of the block (its length is the block's row count). */
struct block {
  const int *index, *rows;
  const double *pointer, *values;
  R_xlen_t n_major, n_rows, col0;
  int by_rows;
};

static struct block block_of(SEXP index, SEXP pointers, SEXP values, SEXP by_rows, SEXP rows, SEXP col0) {
  if (TYPEOF(index) != INTSXP || TYPEOF(pointers) != REALSXP || TYPEOF(values) != REALSXP || TYPEOF(rows) != INTSXP ||
      XLENGTH(index) != XLENGTH(values) || XLENGTH(pointers) < 1) {
    Rf_error("A block is integer indices, and pointers and values as doubles.");
  }
  struct block b;
  b.index = INTEGER(index);
  b.rows = INTEGER(rows);
  b.pointer = REAL(pointers);
  b.values = REAL(values);
  b.n_major = XLENGTH(pointers) - 1;
  b.n_rows = XLENGTH(rows);
  b.by_rows = Rf_asLogical(by_rows) == TRUE;
  double c0 = Rf_asReal(col0);
  if (!(c0 >= 0) || c0 != floor(c0)) Rf_error("A block's first column is a whole number.");
  b.col0 = (R_xlen_t)c0;
  if (b.by_rows && b.n_major != b.n_rows) Rf_error("A block's rows do not match its lines.");
  double n = (double)XLENGTH(index), previous = 0;
  if (b.pointer[0] != 0 || b.pointer[b.n_major] != n) Rf_error("A block's pointers do not match its values.");
  for (R_xlen_t j = 1; j <= b.n_major; j++) {
    if (!(b.pointer[j] >= previous)) Rf_error("A block's pointers do not ascend.");
    previous = b.pointer[j];
  }
  return b;
}

/* The minor index at k, checked against `extent` positions. */
static inline R_xlen_t minor_at(const struct block *b, int64_t k, R_xlen_t extent) {
  int minor = b->index[k];
  if (minor < 0 || (R_xlen_t)minor >= extent) Rf_error("A block's index %d lies outside it.", minor);
  return (R_xlen_t)minor;
}

static void *external(SEXP pointer, const char *what) {
  if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == NULL) Rf_error("%s was released.", what);
  return R_ExternalPtrAddr(pointer);
}

/* --- 1. Counting --------------------------------------------------------------- */

struct tally {
  R_xlen_t n_rows;
  double *count, *low, *high; /* per row: stored non-zero values, smallest and largest finite */
  unsigned char *fraction;    /* per row: 1 once a finite value is not a whole number */
  /* values seen (zeros and non-finite ones included), then NaN, NA, +Inf, -Inf, zeros */
  double seen, nan, missing, positive, negative, zeros;
};

static void tally_free(struct tally *t) {
  if (t == NULL) return;
  free(t->count);
  free(t->low);
  free(t->high);
  free(t->fraction);
  free(t);
}

static void tally_finalizer(SEXP pointer) {
  struct tally *t = (struct tally *)R_ExternalPtrAddr(pointer);
  R_ClearExternalPtr(pointer);
  tally_free(t);
}

SEXP cellucid_backed_tally(SEXP n_rows) {
  double n = Rf_asReal(n_rows);
  if (!(n >= 0) || n != floor(n) || n > (double)R_XLEN_T_MAX) Rf_error("A tally's rows are a whole number.");
  struct tally *t = (struct tally *)calloc(1, sizeof(struct tally));
  if (t == NULL) Rf_error("Out of memory for a matrix's counts.");
  t->n_rows = (R_xlen_t)n;
  size_t rows = (size_t)(t->n_rows > 0 ? t->n_rows : 1);
  t->count = (double *)calloc(rows, sizeof(double));
  t->low = (double *)malloc(rows * sizeof(double));
  t->high = (double *)malloc(rows * sizeof(double));
  t->fraction = (unsigned char *)calloc(rows, 1);
  if (t->count == NULL || t->low == NULL || t->high == NULL || t->fraction == NULL) {
    tally_free(t);
    Rf_error("Out of memory for a matrix's counts.");
  }
  for (R_xlen_t r = 0; r < t->n_rows; r++) {
    t->low[r] = R_PosInf;
    t->high[r] = R_NegInf;
  }
  SEXP pointer = PROTECT(R_MakeExternalPtr(t, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(pointer, tally_finalizer, TRUE);
  UNPROTECT(1);
  return pointer;
}

/* Count a block whose rows are source rows `rows` (0-based). */
SEXP cellucid_backed_count(SEXP tally, SEXP index, SEXP pointers, SEXP values, SEXP by_rows, SEXP rows, SEXP col0) {
  struct tally *t = (struct tally *)external(tally, "A matrix's counts");
  struct block b = block_of(index, pointers, values, by_rows, rows, col0);
  for (R_xlen_t r = 0; r < b.n_rows; r++) {
    if (b.rows[r] < 0 || (R_xlen_t)b.rows[r] >= t->n_rows) Rf_error("A block's row %d lies outside the matrix's %.0f.", b.rows[r] + 1, (double)t->n_rows);
  }
  for (R_xlen_t j = 0; j < b.n_major; j++) {
    for (int64_t k = (int64_t)b.pointer[j]; k < (int64_t)b.pointer[j + 1]; k++) {
      R_xlen_t row = b.rows[b.by_rows ? j : minor_at(&b, k, b.n_rows)];
      double v = b.values[k];
      if (!isfinite(v)) {
        if (ISNA(v)) t->missing++;
        else if (isnan(v)) t->nan++;
        else if (v > 0) t->positive++;
        else t->negative++;
        continue;
      }
      if (v == 0) {
        t->zeros++;
        continue;
      }
      t->count[row] += 1;
      if (v < t->low[row]) t->low[row] = v;
      if (v > t->high[row]) t->high[row] = v;
      if (v != floor(v)) t->fraction[row] = 1;
    }
  }
  t->seen += (double)XLENGTH(index);
  return R_NilValue;
}

/* list(count, low, high, whole, c(seen, nan, missing, positive, negative, zeros)), per row. */
SEXP cellucid_backed_counted(SEXP tally) {
  struct tally *t = (struct tally *)external(tally, "A matrix's counts");
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP count = PROTECT(Rf_allocVector(REALSXP, t->n_rows));
  SEXP low = PROTECT(Rf_allocVector(REALSXP, t->n_rows));
  SEXP high = PROTECT(Rf_allocVector(REALSXP, t->n_rows));
  SEXP whole = PROTECT(Rf_allocVector(LGLSXP, t->n_rows));
  for (R_xlen_t r = 0; r < t->n_rows; r++) {
    REAL(count)[r] = t->count[r];
    REAL(low)[r] = t->low[r];
    REAL(high)[r] = t->high[r];
    LOGICAL(whole)[r] = t->fraction[r] ? FALSE : TRUE;
  }
  SEXP summary = PROTECT(Rf_allocVector(REALSXP, 6));
  double *s = REAL(summary);
  s[0] = t->seen;
  s[1] = t->nan;
  s[2] = t->missing;
  s[3] = t->positive;
  s[4] = t->negative;
  s[5] = t->zeros;
  SET_VECTOR_ELT(out, 0, count);
  SET_VECTOR_ELT(out, 1, low);
  SET_VECTOR_ELT(out, 2, high);
  SET_VECTOR_ELT(out, 3, whole);
  SET_VECTOR_ELT(out, 4, summary);
  UNPROTECT(6);
  return out;
}

/* --- What a placing or prefix pass reads: genes and cells mapped ------------------- */

/* A block's rows named by their output genes (`genes`, -1: left out), and
 * `new_row`, each source column's 0-based output row (-1: left out). */
static const int *new_rows_of(SEXP new_row, const struct block *b, R_xlen_t *n_columns) {
  if (TYPEOF(new_row) != INTSXP) Rf_error("The new rows are an integer vector.");
  *n_columns = XLENGTH(new_row);
  if (!b->by_rows && b->col0 + b->n_major > *n_columns) Rf_error("A block's columns lie outside the matrix.");
  return INTEGER(new_row);
}

/* The output row of the cell at position k of line j, or -1. */
static inline int row_of_cell(const struct block *b, const int *row_of, R_xlen_t n_columns, R_xlen_t j, int64_t k) {
  R_xlen_t column = b->col0 + (b->by_rows ? minor_at(b, k, n_columns - b->col0) : j);
  return row_of[column];
}

/* --- 2. Prefix counts ----------------------------------------------------------- */

struct prefix_tally {
  R_xlen_t n_genes, n_levels;
  double *levels;
  double *counts; /* n_genes x n_levels, column-major, by bucket until read */
};

static void prefix_free(struct prefix_tally *t) {
  if (t == NULL) return;
  free(t->levels);
  free(t->counts);
  free(t);
}

static void prefix_finalizer(SEXP pointer) {
  struct prefix_tally *t = (struct prefix_tally *)R_ExternalPtrAddr(pointer);
  R_ClearExternalPtr(pointer);
  prefix_free(t);
}

SEXP cellucid_backed_prefix_tally(SEXP n_genes, SEXP levels) {
  double genes = Rf_asReal(n_genes);
  if (TYPEOF(levels) != REALSXP || !(genes >= 0) || genes != floor(genes)) Rf_error("Prefix counts need a gene count and levels as doubles.");
  struct prefix_tally *t = (struct prefix_tally *)calloc(1, sizeof(struct prefix_tally));
  if (t == NULL) Rf_error("Out of memory for the prefix counts.");
  t->n_genes = (R_xlen_t)genes;
  t->n_levels = XLENGTH(levels);
  size_t cells = (size_t)(t->n_genes * t->n_levels);
  t->levels = (double *)malloc((size_t)(t->n_levels > 0 ? t->n_levels : 1) * sizeof(double));
  t->counts = (double *)calloc(cells > 0 ? cells : 1, sizeof(double));
  if (t->levels == NULL || t->counts == NULL) {
    prefix_free(t);
    Rf_error("Out of memory for the prefix counts.");
  }
  for (R_xlen_t j = 0; j < t->n_levels; j++) t->levels[j] = REAL(levels)[j];
  SEXP pointer = PROTECT(R_MakeExternalPtr(t, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(pointer, prefix_finalizer, TRUE);
  UNPROTECT(1);
  return pointer;
}

/* Each stored value counts toward every level above its cell's output row;
 * the block's rows are output genes `genes`. */
SEXP cellucid_backed_prefix_count(SEXP tally, SEXP index, SEXP pointers, SEXP values, SEXP by_rows, SEXP genes, SEXP col0,
                                  SEXP new_row) {
  struct prefix_tally *t = (struct prefix_tally *)external(tally, "The prefix counts");
  struct block b = block_of(index, pointers, values, by_rows, genes, col0);
  R_xlen_t n_columns;
  const int *row_of = new_rows_of(new_row, &b, &n_columns);
  for (R_xlen_t j = 0; j < b.n_major; j++) {
    for (int64_t k = (int64_t)b.pointer[j]; k < (int64_t)b.pointer[j + 1]; k++) {
      if (b.values[k] == 0) continue;
      int gene = b.rows[b.by_rows ? j : minor_at(&b, k, b.n_rows)];
      if (gene < 0 || gene >= t->n_genes) continue;
      int row = row_of_cell(&b, row_of, n_columns, j, k);
      if (row < 0) continue;
      /* The first level above the row. */
      R_xlen_t low = 0, high = t->n_levels;
      while (low < high) {
        R_xlen_t middle = low + (high - low) / 2;
        if (t->levels[middle] <= (double)row) low = middle + 1;
        else high = middle;
      }
      if (low < t->n_levels) t->counts[gene + low * t->n_genes] += 1;
    }
  }
  return R_NilValue;
}

/* The counts, n_genes x n_levels: values below each level, cumulative. */
SEXP cellucid_backed_prefix_counted(SEXP tally) {
  struct prefix_tally *t = (struct prefix_tally *)external(tally, "The prefix counts");
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)t->n_genes, (int)t->n_levels));
  double *counts = REAL(out);
  for (R_xlen_t g = 0; g < t->n_genes; g++) {
    double running = 0;
    for (R_xlen_t j = 0; j < t->n_levels; j++) {
      running += t->counts[g + j * t->n_genes];
      counts[g + j * t->n_genes] = running;
    }
  }
  UNPROTECT(1);
  return out;
}

/* --- 3. Placing ----------------------------------------------------------------- */

struct target {
  int64_t base, length;      /* the positions held: base .. base + length - 1 */
  R_xlen_t g_first, n_held;  /* the output genes held: g_first .. g_first + n_held - 1 */
  int64_t *start, *cursor;   /* per gene held: its first position and its next free one */
  int64_t n_cells;           /* output rows: every row is below it */
  int cell_size, value_size, bits;
  unsigned char *cells, *values;
  int mapped;
  struct scratch_map maps[2];
};

static void target_release(struct target *t) {
  if (t == NULL) return;
  if (t->mapped) {
    cellucid_unmap_scratch(&t->maps[0]);
    cellucid_unmap_scratch(&t->maps[1]);
  } else {
    free(t->cells);
    free(t->values);
  }
  free(t->start);
  free(t->cursor);
  free(t);
}

static void target_finalizer(SEXP pointer) {
  struct target *t = (struct target *)R_ExternalPtrAddr(pointer);
  R_ClearExternalPtr(pointer);
  target_release(t);
}

/*
 * A target for the output genes reaching into positions low .. high - 1 of
 * the gene-major arrays whose gene pointers are `pointers` (doubles, n_genes
 * + 1, the first 0): every value of those genes, whole (a gene is sorted
 * whole), in memory, or in the two scratch files `paths` (cell indices,
 * values) when given. `n_cells`: the output rows.
 */
SEXP cellucid_backed_target(SEXP pointers, SEXP low, SEXP high, SEXP n_cells, SEXP index_size, SEXP dtype, SEXP paths) {
  if (TYPEOF(pointers) != REALSXP || XLENGTH(pointers) < 2) Rf_error("A target needs the gene pointers as doubles.");
  if (TYPEOF(paths) != STRSXP || (XLENGTH(paths) != 0 && XLENGTH(paths) != 2)) Rf_error("A target's scratch paths are none or two.");
  const double *pointer = REAL(pointers);
  R_xlen_t n_genes = XLENGTH(pointers) - 1;
  double total = pointer[n_genes];
  double from = Rf_asReal(low), to = Rf_asReal(high);
  if (!(from >= 0) || !(to > from) || to > total || from != floor(from) || to != floor(to)) Rf_error("A target's positions lie outside the matrix's values.");
  if (TYPEOF(dtype) != STRSXP || XLENGTH(dtype) != 1) Rf_error("Expression values are stored as <f4 or <f8.");
  const char *name = CHAR(STRING_ELT(dtype, 0));
  int value_size = !strcmp(name, "<f4") ? 4 : !strcmp(name, "<f8") ? 8 : 0;
  if (value_size == 0) Rf_error("Expression values are stored as <f4 or <f8.");
  int cell_size = Rf_asInteger(index_size);
  if (cell_size != 4 && cell_size != 8) Rf_error("Cell indices are 4 or 8 bytes.");
  double cells = Rf_asReal(n_cells);
  if (!(cells >= 0) || cells != floor(cells)) Rf_error("A target's rows are a whole number.");
  /* The first gene with a value at `from`, and the last at `to - 1`. */
  R_xlen_t a = 0, b = n_genes;
  while (b - a > 1) {
    R_xlen_t m = a + (b - a) / 2;
    if (pointer[m] <= from) a = m; else b = m;
  }
  R_xlen_t g_first = a;
  a = g_first;
  b = n_genes;
  while (b - a > 1) {
    R_xlen_t m = a + (b - a) / 2;
    if (pointer[m] <= to - 1) a = m; else b = m;
  }
  R_xlen_t g_last = a;
  struct target *t = (struct target *)calloc(1, sizeof(struct target));
  if (t == NULL) Rf_error("Out of memory for a gene-major block.");
  cellucid_scratch_init(&t->maps[0]);
  cellucid_scratch_init(&t->maps[1]);
  SEXP pointer_out = PROTECT(R_MakeExternalPtr(t, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(pointer_out, target_finalizer, TRUE);
  t->g_first = g_first;
  t->n_held = g_last - g_first + 1;
  t->base = (int64_t)pointer[g_first];
  t->length = (int64_t)pointer[g_last + 1] - t->base;
  t->n_cells = (int64_t)cells;
  t->cell_size = cell_size;
  t->value_size = value_size;
  t->bits = 0;
  while (t->bits < 63 && ((int64_t)1 << t->bits) < t->n_cells) t->bits++;
  t->start = (int64_t *)malloc((size_t)(t->n_held + 1) * sizeof(int64_t));
  t->cursor = (int64_t *)malloc((size_t)t->n_held * sizeof(int64_t));
  if (t->start == NULL || t->cursor == NULL) Rf_error("Out of memory for a gene-major block.");
  for (R_xlen_t g = 0; g <= t->n_held; g++) t->start[g] = (int64_t)pointer[g_first + g];
  for (R_xlen_t g = 0; g < t->n_held; g++) t->cursor[g] = t->start[g];
  size_t lengths[2] = {(size_t)t->length * (size_t)cell_size, (size_t)t->length * (size_t)value_size};
  if (XLENGTH(paths) == 2) {
    t->mapped = 1;
    for (int part = 0; part < 2; part++) {
      const char *path = Rf_translateChar(STRING_ELT(paths, part));
      const char *failed = cellucid_map_scratch(&t->maps[part], path, lengths[part]);
      if (failed != NULL) Rf_error("Laying a matrix out gene-major %s: %s.", failed, path);
    }
    t->cells = t->maps[0].bytes;
    t->values = t->maps[1].bytes;
  } else {
    t->cells = (unsigned char *)malloc(lengths[0] > 0 ? lengths[0] : 1);
    t->values = (unsigned char *)malloc(lengths[1] > 0 ? lengths[1] : 1);
    if (t->cells == NULL || t->values == NULL) Rf_error("Out of memory for a gene-major block of %.0f values.", (double)t->length);
  }
  UNPROTECT(1);
  return pointer_out;
}

/* The output genes a target holds, c(first, last), 0-based. */
SEXP cellucid_backed_target_genes(SEXP target) {
  struct target *t = (struct target *)external(target, "A gene-major block");
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  REAL(out)[0] = (double)t->g_first;
  REAL(out)[1] = (double)(t->g_first + t->n_held - 1);
  UNPROTECT(1);
  return out;
}

static void gene_range(const struct target *t, SEXP first_gene, SEXP last_gene, R_xlen_t *first, R_xlen_t *last) {
  double a = Rf_asReal(first_gene), b = Rf_asReal(last_gene);
  if (!(a >= (double)t->g_first) || !(b >= a) || b > (double)(t->g_first + t->n_held - 1)) Rf_error("The genes lie outside the block.");
  *first = (R_xlen_t)a;
  *last = (R_xlen_t)b;
}

/* One thread's share of a placing: the output genes first .. last of a block. */
struct place_job {
  struct target *t;
  const struct block *b;
  const int *row_of;
  R_xlen_t first, last;
  int more_gene; /* 0, or the 1-based gene given more values than counted */
  int bad_row;   /* 1 when a cell's new row lies outside the rows */
};

static void place_some(struct place_job *job) {
  struct target *t = job->t;
  const struct block *b = job->b;
  for (R_xlen_t j = 0; j < b->n_major; j++) {
    int line_row = b->by_rows ? -1 : job->row_of[b->col0 + j];
    int line_gene = b->by_rows ? b->rows[j] : -1;
    if (!b->by_rows && line_row < 0) continue;
    if (b->by_rows && (line_gene < job->first || line_gene > job->last)) continue;
    for (int64_t k = (int64_t)b->pointer[j]; k < (int64_t)b->pointer[j + 1]; k++) {
      double v = b->values[k];
      if (v == 0) continue;
      int gene = b->by_rows ? line_gene : b->rows[b->index[k]];
      if (gene < job->first || gene > job->last) continue;
      int row = b->by_rows ? job->row_of[b->col0 + b->index[k]] : line_row;
      if (row < 0) continue;
      if ((int64_t)row >= t->n_cells) {
        job->bad_row = 1;
        return;
      }
      R_xlen_t h = (R_xlen_t)gene - t->g_first;
      int64_t position = t->cursor[h];
      if (position >= t->start[h + 1]) {
        job->more_gene = gene + 1;
        return;
      }
      t->cursor[h] = position + 1;
      size_t at = (size_t)(position - t->base);
      if (t->cell_size == 4) put32(t->cells + 4 * at, (uint32_t)row);
      else put64(t->cells + 8 * at, (uint64_t)row);
      if (t->value_size == 4) put_f4(t->values + 4 * at, v);
      else put_f8(t->values + 8 * at, v);
    }
  }
}

#ifdef _WIN32
static DWORD WINAPI place_worker(LPVOID argument) {
  place_some((struct place_job *)argument);
  return 0;
}
#else
static void *place_worker(void *argument) {
  place_some((struct place_job *)argument);
  return NULL;
}
#endif

/* Stored values below which a block is placed on one thread: a share of a
 * block's values is then cheaper than starting a thread (a scheduling
 * choice, not a limit). */
#define PLACE_THREADED_VALUES 65536

/*
 * Place the values of output genes first_gene .. last_gene that a block
 * holds (its rows named by their output genes, `genes`) at their genes'
 * cursors, each with its cell's output row. The genes are shared among
 * `threads` threads by their values (each thread reads the whole block and
 * writes only its genes, so the bytes never depend on their number: a
 * gene's values arrive in the block's order whichever thread places them).
 * A gene given more values than it was counted (the matrix changed since it
 * was counted) is refused.
 */
SEXP cellucid_backed_place(SEXP target, SEXP first_gene, SEXP last_gene, SEXP index, SEXP pointers, SEXP values,
                           SEXP by_rows, SEXP genes, SEXP col0, SEXP new_row, SEXP threads) {
  struct target *t = (struct target *)external(target, "A gene-major block");
  R_xlen_t first, last;
  gene_range(t, first_gene, last_gene, &first, &last);
  struct block b = block_of(index, pointers, values, by_rows, genes, col0);
  R_xlen_t n_columns;
  const int *row_of = new_rows_of(new_row, &b, &n_columns);
  /* Every index inside the block, checked here so the threads need not. */
  R_xlen_t extent = b.by_rows ? n_columns - b.col0 : b.n_rows;
  R_xlen_t n = XLENGTH(index);
  for (R_xlen_t k = 0; k < n; k++) {
    if (b.index[k] < 0 || (R_xlen_t)b.index[k] >= extent) Rf_error("A block's index %d lies outside it.", b.index[k]);
  }
  int n_threads = Rf_asInteger(threads);
  if (n_threads < 1 || n_threads == NA_INTEGER || n < PLACE_THREADED_VALUES) n_threads = 1;
  if ((R_xlen_t)n_threads > last - first + 1) n_threads = (int)(last - first + 1);
  struct place_job *jobs = (struct place_job *)R_alloc((size_t)n_threads, sizeof(struct place_job));
  /* Genes split by their positions: about as many values each (a share may be empty). */
  int64_t from = t->start[first - t->g_first], to = t->start[last - t->g_first + 1];
  R_xlen_t gene = first;
  for (int k = 0; k < n_threads; k++) {
    R_xlen_t end = last;
    if (k < n_threads - 1) {
      int64_t until = from + (to - from) * (k + 1) / n_threads;
      end = gene - 1;
      while (end < last && t->start[end + 1 - t->g_first] < until) end++;
    }
    jobs[k].t = t;
    jobs[k].b = &b;
    jobs[k].row_of = row_of;
    jobs[k].first = gene;
    jobs[k].last = end;
    jobs[k].more_gene = 0;
    jobs[k].bad_row = 0;
    gene = end + 1;
  }
#ifdef _WIN32
  HANDLE *handles = (HANDLE *)R_alloc((size_t)n_threads, sizeof(HANDLE));
  for (int k = 1; k < n_threads; k++) {
    handles[k] = CreateThread(NULL, 0, place_worker, &jobs[k], 0, NULL);
    /* A thread that could not start: its genes placed here. */
    if (handles[k] == NULL) place_some(&jobs[k]);
  }
  place_some(&jobs[0]);
  for (int k = 1; k < n_threads; k++) {
    if (handles[k] != NULL) {
      WaitForSingleObject(handles[k], INFINITE);
      CloseHandle(handles[k]);
    }
  }
#else
  pthread_t *handles = (pthread_t *)R_alloc((size_t)n_threads, sizeof(pthread_t));
  int *running = (int *)R_alloc((size_t)n_threads, sizeof(int));
  for (int k = 1; k < n_threads; k++) {
    running[k] = pthread_create(&handles[k], NULL, place_worker, &jobs[k]) == 0;
    /* A thread that could not start: its genes placed here. */
    if (!running[k]) place_some(&jobs[k]);
  }
  place_some(&jobs[0]);
  for (int k = 1; k < n_threads; k++) {
    if (running[k]) pthread_join(handles[k], NULL);
  }
#endif
  for (int k = 0; k < n_threads; k++) {
    if (jobs[k].bad_row) Rf_error("A new row lies outside the %.0f rows.", (double)t->n_cells);
    if (jobs[k].more_gene) Rf_error("Gene %d holds more values than were counted: the matrix changed while it was read.", jobs[k].more_gene);
  }
  return R_NilValue;
}

/* The genes a finish sorts, shared by its threads: each takes the next run
 * of genes; a gene holding one cell twice is recorded, never raised on a
 * worker thread (R's API is not thread-safe). */
struct finish_job {
  struct target *t;
  R_xlen_t first, last, next;
  const double *levels;
  R_xlen_t n_levels;
  double *counts; /* (last - first + 1) x n_levels, or NULL */
  R_xlen_t longest;
  int failed_gene;   /* 0, or the 1-based gene holding a cell twice */
  int out_of_memory;
#ifdef _WIN32
  CRITICAL_SECTION lock;
#else
  pthread_mutex_t lock;
#endif
};

/* Genes a worker takes at a time. */
#define FINISH_RUN 64

static R_xlen_t finish_take(struct finish_job *job) {
  R_xlen_t g;
#ifdef _WIN32
  EnterCriticalSection(&job->lock);
#else
  pthread_mutex_lock(&job->lock);
#endif
  g = job->next;
  job->next += FINISH_RUN;
#ifdef _WIN32
  LeaveCriticalSection(&job->lock);
#else
  pthread_mutex_unlock(&job->lock);
#endif
  return g;
}

static void finish_gene(struct finish_job *job, R_xlen_t g, struct renumbered *pairs, struct renumbered *spare) {
  struct target *t = job->t;
  R_xlen_t h = g - t->g_first;
  size_t at = (size_t)(t->start[h] - t->base);
  R_xlen_t n = (R_xlen_t)(t->start[h + 1] - t->start[h]);
  int cs = t->cell_size, vs = t->value_size;
  unsigned char *cells = t->cells + (size_t)cs * at;
  unsigned char *values = t->values + (size_t)vs * at;
  int ascending = 1;
  int64_t previous = -1;
  for (R_xlen_t k = 0; k < n; k++) {
    int64_t row = cs == 4 ? (int64_t)get32(cells + 4 * k) : (int64_t)get64(cells + 8 * k);
    if (row <= previous) {
      ascending = 0;
      break;
    }
    previous = row;
  }
  if (!ascending) {
    for (R_xlen_t k = 0; k < n; k++) {
      pairs[k].row = cs == 4 ? (int64_t)get32(cells + 4 * k) : (int64_t)get64(cells + 8 * k);
      pairs[k].value = vs == 4 ? get_f4(values + 4 * k) : get_f8(values + 8 * k);
    }
    cellucid_sort_by_row(pairs, spare, n, t->bits);
    for (R_xlen_t k = 0; k < n; k++) {
      if (k > 0 && pairs[k].row == pairs[k - 1].row) {
        __atomic_store_n(&job->failed_gene, (int)g + 1, __ATOMIC_RELAXED);
        return;
      }
      if (cs == 4) put32(cells + 4 * k, (uint32_t)pairs[k].row);
      else put64(cells + 8 * k, (uint64_t)pairs[k].row);
      if (vs == 4) put_f4(values + 4 * k, pairs[k].value);
      else put_f8(values + 8 * k, pairs[k].value);
    }
  }
  R_xlen_t n_genes = job->last - job->first + 1;
  for (R_xlen_t j = 0; j < job->n_levels; j++) {
    /* The first value at or past the level: rows ascend in a gene. */
    double level = job->levels[j];
    R_xlen_t low = 0, high = n;
    while (low < high) {
      R_xlen_t middle = low + (high - low) / 2;
      double row = cs == 4 ? (double)get32(cells + 4 * middle) : (double)get64(cells + 8 * middle);
      if (row < level) low = middle + 1;
      else high = middle;
    }
    job->counts[(g - job->first) + j * n_genes] = (double)low;
  }
}

static void finish_some(struct finish_job *job) {
  struct renumbered *pairs = NULL, *spare = NULL;
  if (job->longest > 0) {
    pairs = (struct renumbered *)malloc((size_t)job->longest * sizeof(struct renumbered));
    spare = (struct renumbered *)malloc((size_t)job->longest * sizeof(struct renumbered));
    if (pairs == NULL || spare == NULL) {
      free(pairs);
      free(spare);
      __atomic_store_n(&job->out_of_memory, 1, __ATOMIC_RELAXED);
      return;
    }
  }
  for (;;) {
    R_xlen_t from = finish_take(job);
    if (from > job->last || __atomic_load_n(&job->failed_gene, __ATOMIC_RELAXED)) break;
    R_xlen_t to = from + FINISH_RUN - 1 < job->last ? from + FINISH_RUN - 1 : job->last;
    for (R_xlen_t g = from; g <= to; g++) finish_gene(job, g, pairs, spare);
  }
  free(pairs);
  free(spare);
}

#ifdef _WIN32
static DWORD WINAPI finish_worker(LPVOID argument) {
  finish_some((struct finish_job *)argument);
  return 0;
}
#else
static void *finish_worker(void *argument) {
  finish_some((struct finish_job *)argument);
  return NULL;
}
#endif

/*
 * Output genes first_gene .. last_gene placed: each checked to hold every
 * value it was counted, then put in ascending row (one radix sort a gene,
 * none for a gene already ascending), the genes shared among `threads`
 * threads (each sorts its genes alone, so the bytes never depend on their
 * number). Returns the genes' prefix counts (rows below each of `levels`, a
 * genes x levels matrix), or NULL without levels.
 */
SEXP cellucid_backed_finish(SEXP target, SEXP first_gene, SEXP last_gene, SEXP levels, SEXP threads) {
  struct target *t = (struct target *)external(target, "A gene-major block");
  R_xlen_t first, last;
  gene_range(t, first_gene, last_gene, &first, &last);
  if (TYPEOF(levels) != REALSXP && !Rf_isNull(levels)) Rf_error("Prefix levels are doubles.");
  struct finish_job job;
  memset(&job, 0, sizeof(job));
  job.t = t;
  job.first = first;
  job.last = last;
  job.next = first;
  job.n_levels = Rf_isNull(levels) ? 0 : XLENGTH(levels);
  job.levels = job.n_levels > 0 ? REAL(levels) : NULL;
  for (R_xlen_t g = first; g <= last; g++) {
    R_xlen_t h = g - t->g_first;
    if (t->cursor[h] != t->start[h + 1]) Rf_error("Gene %d holds fewer values than were counted: the matrix changed while it was read.", (int)g + 1);
    R_xlen_t length = (R_xlen_t)(t->start[h + 1] - t->start[h]);
    if (length > job.longest) job.longest = length;
  }
  SEXP out = PROTECT(job.n_levels > 0 ? Rf_allocMatrix(REALSXP, (int)(last - first + 1), (int)job.n_levels) : R_NilValue);
  job.counts = job.n_levels > 0 ? REAL(out) : NULL;
  int n_threads = Rf_asInteger(threads);
  if (n_threads < 1 || n_threads == NA_INTEGER) n_threads = 1;
  R_xlen_t runs = (last - first) / FINISH_RUN + 1;
  if ((R_xlen_t)n_threads > runs) n_threads = (int)runs;
#ifdef _WIN32
  InitializeCriticalSection(&job.lock);
  HANDLE *handles = (HANDLE *)calloc((size_t)n_threads, sizeof(HANDLE));
  int started = 0;
  if (handles != NULL) {
    for (; started < n_threads - 1; started++) {
      handles[started] = CreateThread(NULL, 0, finish_worker, &job, 0, NULL);
      if (handles[started] == NULL) break;
    }
  }
  finish_some(&job);
  for (int k = 0; k < started; k++) {
    WaitForSingleObject(handles[k], INFINITE);
    CloseHandle(handles[k]);
  }
  free(handles);
  DeleteCriticalSection(&job.lock);
#else
  pthread_mutex_init(&job.lock, NULL);
  pthread_t *handles = (pthread_t *)calloc((size_t)n_threads, sizeof(pthread_t));
  int started = 0;
  if (handles != NULL) {
    for (; started < n_threads - 1; started++) {
      if (pthread_create(&handles[started], NULL, finish_worker, &job) != 0) break;
    }
  }
  /* R's thread works too; with fewer threads started, the rest of the genes are still taken. */
  finish_some(&job);
  for (int k = 0; k < started; k++) pthread_join(handles[k], NULL);
  free(handles);
  pthread_mutex_destroy(&job.lock);
#endif
  if (job.out_of_memory) Rf_error("Out of memory to sort a gene of %.0f values.", (double)job.longest);
  if (job.failed_gene) Rf_error("Gene %d holds two values for one cell.", job.failed_gene);
  UNPROTECT(1);
  return out;
}

/*
 * Positions first .. first + count - 1 (0-based, of the whole arrays) of part
 * `part` (1: cell indices, 2: values) of a target, then zeros up to `total`
 * items: a chunk, or the range of one a request asks for.
 */
SEXP cellucid_backed_slice(SEXP target, SEXP part, SEXP first, SEXP count, SEXP total) {
  struct target *t = (struct target *)external(target, "A gene-major block");
  int which = Rf_asInteger(part);
  if (which != 1 && which != 2) Rf_error("A gene-major block's part is 1 (cell indices) or 2 (values).");
  double from = Rf_asReal(first), n = Rf_asReal(count), size = Rf_asReal(total);
  if (!(from >= (double)t->base) || !(n >= 0) || !(size >= n) || from + n > (double)(t->base + t->length) ||
      from != floor(from) || n != floor(n) || size != floor(size)) {
    Rf_error("The slice lies outside the gene-major block.");
  }
  int item = which == 1 ? t->cell_size : t->value_size;
  double bytes = size * item;
  if (bytes > (double)R_XLEN_T_MAX) Rf_error("A slice of %.0f bytes is more than one R raw vector holds.", bytes);
  SEXP out = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)bytes));
  const unsigned char *source = (which == 1 ? t->cells : t->values) + (size_t)(from - (double)t->base) * (size_t)item;
  size_t copied = (size_t)n * (size_t)item;
  if (copied > 0) memcpy(RAW(out), source, copied);
  if ((size_t)bytes > copied) memset(RAW(out) + copied, 0, (size_t)bytes - copied);
  UNPROTECT(1);
  return out;
}

/* Let a target's memory or scratch files go now (it is not used again). */
SEXP cellucid_backed_release(SEXP target) {
  if (TYPEOF(target) == EXTPTRSXP) target_finalizer(target);
  return R_NilValue;
}
