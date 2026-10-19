/*
 * Latent-space outlier quantiles: the rule cellucid-python's prepare() and
 * cellucid serve, and the web app's data/latent-outliers.js, compute, bit for
 * bit. For each category with at least `min_points` cells, over the latent
 * values as stored (each rounded to float32, then widened exactly to double):
 *
 * 1. the centroid is the column sums, accumulated cell by cell in cell order,
 *    divided by the cell count;
 * 2. a cell's distance is the square root of its squared differences from the
 *    centroid, accumulated column by column in column order;
 * 3. its quantile is the number of the category's cells at that distance or
 *    closer, divided by the cell count, both as float32.
 *
 * Every other cell (a missing code, or a category under `min_points`) is NaN.
 *
 * Only IEEE double and float arithmetic is used, each operation rounded on
 * its own: no long double, and every product is stored before it is added so
 * that no compiler contracts the two into a fused multiply-add.
 */

#include <R.h>
#include <Rinternals.h>

#include <math.h>
#include <stdlib.h>

static int cellucid_compare_doubles(const void *left, const void *right) {
  const double a = *(const double *)left;
  const double b = *(const double *)right;
  return (a > b) - (a < b);
}

/* The number of values in `sorted` (ascending) at most `value`. */
static R_xlen_t cellucid_count_at_most(const double *sorted, R_xlen_t length,
                                       double value) {
  R_xlen_t low = 0;
  R_xlen_t high = length;
  while (low < high) {
    const R_xlen_t middle = low + (high - low) / 2;
    if (sorted[middle] <= value) {
      low = middle + 1;
    } else {
      high = middle;
    }
  }
  return low;
}

static double cellucid_as_float32(double value) { return (double)(float)value; }

SEXP cellucid_latent_outlier_quantiles(SEXP latent, SEXP codes,
                                       SEXP n_categories_sexp,
                                       SEXP min_points_sexp) {
  if (!Rf_isReal(latent) || !Rf_isMatrix(latent)) {
    Rf_error("The latent space must be a double matrix.");
  }
  if (!Rf_isInteger(codes)) {
    Rf_error("Category codes must be an integer vector.");
  }
  const R_xlen_t n_cells = (R_xlen_t)Rf_nrows(latent);
  const R_xlen_t width = (R_xlen_t)Rf_ncols(latent);
  if (XLENGTH(codes) != n_cells) {
    Rf_error("Category codes must be one integer per cell.");
  }
  const int n_categories = Rf_asInteger(n_categories_sexp);
  const int min_points = Rf_asInteger(min_points_sexp);
  if (n_categories == NA_INTEGER || n_categories < 0) {
    Rf_error("The category count must be a non-negative integer.");
  }
  if (min_points == NA_INTEGER || min_points < 1) {
    Rf_error("The fewest cells a category needs must be a positive integer.");
  }

  const double *values = REAL(latent);
  const int *code_of = INTEGER(codes);
  SEXP result = PROTECT(Rf_allocVector(REALSXP, n_cells));
  double *quantiles = REAL(result);
  for (R_xlen_t cell = 0; cell < n_cells; cell++) {
    quantiles[cell] = R_NaN;
  }

  /* The cells of each category, in cell order. */
  R_xlen_t *starts = (R_xlen_t *)R_alloc((size_t)n_categories + 1, sizeof(R_xlen_t));
  for (int code = 0; code <= n_categories; code++) {
    starts[code] = 0;
  }
  for (R_xlen_t cell = 0; cell < n_cells; cell++) {
    const int code = code_of[cell];
    if (code != NA_INTEGER && code >= 0 && code < n_categories) {
      starts[code + 1] += 1;
    }
  }
  for (int code = 0; code < n_categories; code++) {
    starts[code + 1] += starts[code];
  }
  const R_xlen_t n_members = starts[n_categories];
  R_xlen_t *members = (R_xlen_t *)R_alloc((size_t)(n_members > 0 ? n_members : 1), sizeof(R_xlen_t));
  R_xlen_t *next = (R_xlen_t *)R_alloc((size_t)n_categories + 1, sizeof(R_xlen_t));
  for (int code = 0; code < n_categories; code++) {
    next[code] = starts[code];
  }
  for (R_xlen_t cell = 0; cell < n_cells; cell++) {
    const int code = code_of[cell];
    if (code != NA_INTEGER && code >= 0 && code < n_categories) {
      members[next[code]++] = cell;
    }
  }

  for (int code = 0; code < n_categories; code++) {
    const R_xlen_t start = starts[code];
    const R_xlen_t n = starts[code + 1] - start;
    if (n < (R_xlen_t)min_points) {
      continue;
    }
    double *squared = (double *)R_alloc((size_t)n, sizeof(double));
    double *sorted = (double *)R_alloc((size_t)n, sizeof(double));
    for (R_xlen_t k = 0; k < n; k++) {
      squared[k] = 0.0;
    }
    for (R_xlen_t column = 0; column < width; column++) {
      const double *column_values = values + column * n_cells;
      double sum = 0.0;
      for (R_xlen_t k = 0; k < n; k++) {
        sum += cellucid_as_float32(column_values[members[start + k]]);
      }
      const double centroid = sum / (double)n;
      for (R_xlen_t k = 0; k < n; k++) {
        const double difference =
            cellucid_as_float32(column_values[members[start + k]]) - centroid;
        volatile double product = difference * difference;
        squared[k] += product;
      }
    }
    for (R_xlen_t k = 0; k < n; k++) {
      squared[k] = sqrt(squared[k]);
      sorted[k] = squared[k];
    }
    qsort(sorted, (size_t)n, sizeof(double), cellucid_compare_doubles);
    const double count = cellucid_as_float32((double)n);
    for (R_xlen_t k = 0; k < n; k++) {
      const double rank =
          cellucid_as_float32((double)cellucid_count_at_most(sorted, n, squared[k]));
      quantiles[members[start + k]] = cellucid_as_float32(rank / count);
    }
  }

  UNPROTECT(1);
  return result;
}
