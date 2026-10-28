/*
 * Priority order: the prefix counts of a gene-major X (docs: "Cellucid-ready
 * AnnData specification", priority order). For each gene and each prefix
 * level (ascending), how many of the gene's stored values lie in the cells
 * before that level: varm['cellucid_prefix_counts'], as cellucid-python's
 * prefix_counts() counts them. The cells are already in priority order.
 */

#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

/*
 * `p`, `i`: a dgCMatrix's column pointers and 0-based row indices.
 * `cells_by_column` TRUE: the matrix is genes x cells (a column per cell, in
 * priority order); `map` gives each row's output gene (0-based, -1 for a gene
 * not stored). FALSE: it is cells x genes (a column per gene, rows ascending).
 * `levels`: ascending cell counts. Returns an n_genes x n_levels double
 * matrix (counts reach 2^32 - 1, past R's integers).
 */
SEXP cellucid_prefix_counts(SEXP p, SEXP i, SEXP levels, SEXP n_genes, SEXP cells_by_column, SEXP map,
                            SEXP cells, SEXP new_row) {
  const int *pointers = INTEGER(p);
  const int *rows = INTEGER(i);
  R_xlen_t n_columns = XLENGTH(p) - 1;
  R_xlen_t n_levels = XLENGTH(levels);
  R_xlen_t genes = (R_xlen_t)asReal(n_genes);
  const double *level = REAL(levels);
  SEXP out = PROTECT(allocMatrix(REALSXP, (int)genes, (int)n_levels));
  double *counts = REAL(out);
  for (R_xlen_t k = 0; k < genes * n_levels; k++) counts[k] = 0;
  /* `cells`: the order a genes x cells matrix's columns are taken in (none:
   * their own); `new_row`: each cell's row in a cells x genes matrix
   * renumbered (none: its own rows, ascending in each gene). */
  const int *order = XLENGTH(cells) > 0 ? INTEGER(cells) : NULL;
  const int *renumber = XLENGTH(new_row) > 0 ? INTEGER(new_row) : NULL;
  if (asLogical(cells_by_column)) {
    const int *gene_of = INTEGER(map);
    double *running = (double *)R_alloc(genes > 0 ? genes : 1, sizeof(double));
    for (R_xlen_t g = 0; g < genes; g++) running[g] = 0;
    R_xlen_t next = 0;
    for (R_xlen_t cell = 0; cell <= n_columns && next < n_levels; cell++) {
      while (next < n_levels && level[next] == (double)cell) {
        for (R_xlen_t g = 0; g < genes; g++) counts[g + next * genes] = running[g];
        next++;
      }
      if (cell == n_columns) break;
      R_xlen_t column = order == NULL ? cell : order[cell];
      for (int e = pointers[column]; e < pointers[column + 1]; e++) {
        int g = gene_of[rows[e]];
        if (g >= 0) running[g] += 1;
      }
    }
  } else if (renumber != NULL) {
    /* Each value counts toward every level above its new row. */
    for (R_xlen_t g = 0; g < genes; g++) {
      for (int e = pointers[g]; e < pointers[g + 1]; e++) {
        double row = (double)renumber[rows[e]];
        R_xlen_t low = 0, high = n_levels; /* the first level above `row` */
        while (low < high) {
          R_xlen_t middle = low + (high - low) / 2;
          if (level[middle] <= row) low = middle + 1;
          else high = middle;
        }
        if (low < n_levels) counts[g + low * genes] += 1;
      }
      for (R_xlen_t j = 1; j < n_levels; j++) counts[g + j * genes] += counts[g + (j - 1) * genes];
    }
  } else {
    for (R_xlen_t g = 0; g < genes; g++) {
      int start = pointers[g];
      int end = pointers[g + 1];
      for (R_xlen_t j = 0; j < n_levels; j++) {
        /* The first entry at or past the level: rows ascend in a gene. */
        int low = start;
        int high = end;
        while (low < high) {
          int middle = low + (high - low) / 2;
          if ((double)rows[middle] < level[j]) low = middle + 1;
          else high = middle;
        }
        counts[g + j * genes] = (double)(low - start);
      }
    }
  }
  UNPROTECT(1);
  return out;
}
