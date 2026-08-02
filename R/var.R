# The gene axis: naming the genes to export and reading one gene's column
# out of whatever matrix the caller supplied.

.extract_gene_ids <- function(var, var_gene_id_column) {
  # The same four terms .validate_string() applies, plus the non-empty rule,
  # in the shape .resolve_created_at() uses for an optional argument: the
  # message has to name NULL, which the shared validator's does not. The class
  # term is here for the reason it is on every other string argument -- this
  # one is matched against names(var), used to subscript it, and published in
  # var_manifest.json, all of which assume a native string.
  if (
    !is.null(var_gene_id_column) &&
      (
        !is.character(var_gene_id_column) ||
          length(var_gene_id_column) != 1L ||
          is.na(var_gene_id_column) ||
          is.object(var_gene_id_column) ||
          !nzchar(var_gene_id_column)
      )
  ) {
    stop(
      "var_gene_id_column must be NULL or one non-empty string.",
      call. = FALSE
    )
  }

  if (is.null(var_gene_id_column)) {
    # A data frame always has row names, so rownames(var) is never NULL and
    # asking whether it is answers nothing. What a caller who never set them
    # has is the automatic sequence, and rownames() materializes that as the
    # strings "1" to "N" -- unique, drawable, and accepted by every check
    # downstream. Selecting them would publish an export whose genes are named
    # after their own row numbers, so a wet-lab user searching CD8A finds
    # nothing and there is no failure anywhere to say why. .row_names_info()
    # is what tells the two apart: it returns the row count negated while the
    # row names are still the automatic ones, and base R's own
    # as.matrix.data.frame() reads them the same way. The cellucid Python
    # package refuses the same input, because a default RangeIndex yields
    # integers and its gene identifiers must be strings.
    if (.row_names_info(var) < 0L) {
      stop(
        "var has only automatic row names, so rownames(var) would name the ",
        "genes '1' to '", nrow(var), "'. Set rownames(var) to the gene ",
        "identifiers, or pass var_gene_id_column.",
        call. = FALSE
      )
    }
    ids <- rownames(var)
    .validate_gene_ids(ids)
    return(ids)
  }
  if (!var_gene_id_column %in% names(var)) {
    stop(
      "var_gene_id_column '", var_gene_id_column, "' not found in var. Available columns: ",
      .format_value_list(names(var))
    )
  }
  ids <- var[[var_gene_id_column]]
  .validate_gene_ids(ids)
  ids
}

.validate_gene_ids <- function(gene_ids) {
  .require_unique_identifiers(gene_ids, what = "Gene")
}

.get_gene_column <- function(gene_expression, gene_idx, n_cells) {
  if (inherits(gene_expression, "dgCMatrix")) {
    p <- gene_expression@p
    i <- gene_expression@i
    x <- gene_expression@x

    start <- p[[gene_idx]] + 1L
    end <- p[[gene_idx + 1L]]
    values <- numeric(n_cells)
    if (end >= start) {
      rows <- i[start:end] + 1L
      values[rows] <- x[start:end]
    }
    return(values)
  }

  if (inherits(gene_expression, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Matrix package is required to export sparse gene_expression objects.")
    }
    dense <- as.matrix(gene_expression[, gene_idx, drop = FALSE])
    return(as.numeric(dense[, 1L]))
  }

  values <- gene_expression[, gene_idx]
  as.numeric(values)
}
