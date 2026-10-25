/*
 * One pass over a numeric vector, at memory speed and with no copy: its
 * smallest and largest finite value, how many values are not finite (NA,
 * NaN, an infinity), whether one is zero, and whether every finite value is
 * a whole number (when asked: the check costs a floor() per value). cellucid_prepare() and the viewer check every expression
 * value before anything is written (R/matrix-input.R): base R's range(),
 * any(x == 0) and all(x == floor(x)) take three passes and allocate a
 * logical vector as long as the matrix for two of them.
 */

#include <R.h>
#include <Rinternals.h>

#include <math.h>

SEXP cellucid_extent(SEXP values, SEXP check_whole) {
  int want_whole = Rf_asLogical(check_whole) == TRUE;
  double low = R_PosInf, high = R_NegInf, not_finite = 0;
  int zero = 0, whole = 1;
  R_xlen_t n = XLENGTH(values);
  if (TYPEOF(values) == REALSXP) {
    const double *x = REAL(values);
    for (R_xlen_t i = 0; i < n; i++) {
      double v = x[i];
      if (!isfinite(v)) {
        not_finite++;
        continue;
      }
      if (v < low) low = v;
      if (v > high) high = v;
      if (v == 0) zero = 1;
      if (want_whole && whole && v != floor(v)) whole = 0;
    }
  } else if (TYPEOF(values) == INTSXP || TYPEOF(values) == LGLSXP) {
    const int *x = TYPEOF(values) == INTSXP ? INTEGER(values) : LOGICAL(values);
    for (R_xlen_t i = 0; i < n; i++) {
      int v = x[i];
      if (v == NA_INTEGER) {
        not_finite++;
        continue;
      }
      if (v < low) low = v;
      if (v > high) high = v;
      if (v == 0) zero = 1;
    }
  } else {
    Rf_error("The extent is computed over a numeric vector.");
  }
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 5));
  REAL(out)[0] = low;
  REAL(out)[1] = high;
  REAL(out)[2] = not_finite;
  REAL(out)[3] = zero;
  REAL(out)[4] = whole;
  UNPROTECT(1);
  return out;
}

/*
 * The compressed pointers of a column-compressed matrix's transpose, over the
 * rows it keeps: for row indices `i` (0-based) and `map` (one entry per row
 * of the matrix: the 0-based output row it becomes, or -1 when it is left
 * out), the cumulative count of stored values per output row, n_rows + 1 of
 * them starting at 0. A genes x cells dgCMatrix's transpose is the
 * gene-major csc_matrix the store holds, so this is its X/indptr, known
 * without transposing the matrix or copying the genes it keeps
 * (R/expression-input.R).
 */
SEXP cellucid_row_pointers(SEXP i, SEXP map, SEXP n_rows) {
  if (TYPEOF(i) != INTSXP || TYPEOF(map) != INTSXP) Rf_error("Row indices and the row map must be integer vectors.");
  /* Checked as a double before the cast: R_xlen_t is an int on 32-bit R. */
  double wanted = Rf_asReal(n_rows);
  if (!(wanted >= 0) || wanted >= (double)R_XLEN_T_MAX || wanted != floor(wanted)) {
    Rf_error("The row count must be a whole number below %.0f.", (double)R_XLEN_T_MAX);
  }
  R_xlen_t rows = (R_xlen_t)wanted;
  R_xlen_t matrix_rows = XLENGTH(map);
  const int *target = INTEGER(map);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, rows + 1));
  double *pointers = REAL(out);
  for (R_xlen_t r = 0; r <= rows; r++) pointers[r] = 0;
  const int *index = INTEGER(i);
  R_xlen_t n = XLENGTH(i);
  for (R_xlen_t k = 0; k < n; k++) {
    int row = index[k];
    if (row < 0 || row >= matrix_rows) Rf_error("Row index %d is outside the matrix.", row);
    int to = target[row];
    if (to < 0) continue;
    if (to >= rows) Rf_error("The row map names row %d of %d.", to, (int)rows);
    pointers[to + 1] += 1;
  }
  for (R_xlen_t r = 1; r <= rows; r++) pointers[r] += pointers[r - 1];
  UNPROTECT(1);
  return out;
}

/*
 * cellucid_extent over the values of a dgCMatrix (`x`, row indices `i`)
 * whose rows `map` keeps (map[row] >= 0): the extent of the genes stored,
 * without copying their values out first.
 */
SEXP cellucid_kept_extent(SEXP x, SEXP i, SEXP map, SEXP check_whole) {
  if (TYPEOF(x) != REALSXP || TYPEOF(i) != INTSXP || TYPEOF(map) != INTSXP || XLENGTH(x) != XLENGTH(i)) {
    Rf_error("The kept extent needs a dgCMatrix's x and i and an integer row map.");
  }
  int want_whole = Rf_asLogical(check_whole) == TRUE;
  double low = R_PosInf, high = R_NegInf, not_finite = 0;
  int zero = 0, whole = 1;
  const double *values = REAL(x);
  const int *rows = INTEGER(i);
  const int *target = INTEGER(map);
  R_xlen_t matrix_rows = XLENGTH(map);
  R_xlen_t n = XLENGTH(x);
  for (R_xlen_t k = 0; k < n; k++) {
    int row = rows[k];
    if (row < 0 || row >= matrix_rows) Rf_error("Row index %d is outside the matrix.", row);
    if (target[row] < 0) continue;
    double v = values[k];
    if (!isfinite(v)) {
      not_finite++;
      continue;
    }
    if (v < low) low = v;
    if (v > high) high = v;
    if (v == 0) zero = 1;
    if (want_whole && whole && v != floor(v)) whole = 0;
  }
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 5));
  REAL(out)[0] = low;
  REAL(out)[1] = high;
  REAL(out)[2] = not_finite;
  REAL(out)[3] = zero;
  REAL(out)[4] = whole;
  UNPROTECT(1);
  return out;
}
