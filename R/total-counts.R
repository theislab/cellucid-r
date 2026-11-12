# Each cell's total count: obs['total_counts'] and uns['cellucid']['totalCounts']
# (the Cellucid-ready AnnData specification). cellucid-python computes and
# checks the same in src/cellucid/prepare_data/_totals.py.
#
# A file whose X holds counts records each cell's total count -- its library
# size, the t of library-size normalization (x * 10,000 / t) -- so the
# viewer normalizes a gene from one column instead of reading every gene to
# sum them. The record says what the totals cover: "all genes" for the
# caller's column (cellucid_prepare(total_counts = )), each cell's total over
# every gene measured, before any gene was selected (Seurat's nCount_RNA);
# "source panel" for the caller's column over every gene of the source's gene
# panel (total_counts_over = "source panel": a targeted or integration panel,
# such as the built-in atlases' 8,192 genes -- more genes than are stored, not
# the whole transcriptome), checked as "all genes" is; "stored genes" for the sum of each row of X as written, computed when no
# column is named.
#
# X holds counts when every value it stores is a whole number of at least 0.
# Summed in doubles: a sum of whole counts below 2^53 is exact, so both
# writers store the same bytes.

.total_counts_column <- "total_counts"
.total_counts_over_all <- "all genes"
.total_counts_over_panel <- "source panel"
# What a caller's named column may cover (total_counts_over =).
.total_counts_over_named <- c(.total_counts_over_all, .total_counts_over_panel)
.total_counts_over_stored <- "stored genes"
# How far below the stored genes' sum a named total may be and still be one.
.total_counts_tolerance <- 1e-6

# Whether an extent (.native_extent(whole = TRUE)'s five numbers) is of
# counts: no stored value, or every one whole and at least 0.
.extent_holds_counts <- function(extent) {
  !is.finite(extent[[1L]]) || (extent[[1L]] >= 0 && extent[[5L]] == 1)
}

# Each cell's sum over the genes kept (`index`: 0-based cell of each value,
# added into `totals`).
.add_by_cell <- function(totals, cells, values) {
  if (length(values) == 0L) return(totals)
  summed <- rowsum(values, cells, reorder = FALSE)
  at <- as.numeric(rownames(summed)) + 1
  totals[at] <- totals[at] + summed[, 1L]
  totals
}

# One block of a matrix held outside R's memory (.backed_blocks()) added to
# each cell's (column's) total.
.add_block_to_cell_totals <- function(totals, block) {
  if (block$by_rows) return(.add_by_cell(totals, block$index, block$values))
  # Each line a cell: its values' sum, from the running sum.
  lines <- length(block$pointers) - 1L
  if (lines <= 0L) return(totals)
  running <- c(0, cumsum(block$values))
  sums <- running[block$pointers[-1L] + 1] - running[block$pointers[-(lines + 1L)] + 1]
  at <- block$col0 + seq_len(lines)
  totals[at] <- totals[at] + sums
  totals
}

# Each row's (cell's) sum over the genes kept, in the matrix's cell order, as
# doubles, when every value of those genes is a whole number of at least 0;
# NULL otherwise (or without a matrix). `expression` is an
# .expression_input(), its genes chosen (.expression_select_genes()) and its
# cells aligned (.align_single_space()), not yet in priority order. A matrix
# held outside R's memory has its sums from its counting pass; with some of
# its genes kept, their rows are read once more, block by block, never whole.
.expression_row_totals <- function(expression) {
  if (is.null(expression)) return(NULL)
  matrix <- expression$matrix
  if (isTRUE(expression$backed)) {
    rows <- expression$genes %||% seq_len(nrow(matrix))
    if (!.extent_holds_counts(.backed_extent(expression$stats, rows))) return(NULL)
    # Summed over every gene by the counting pass (.backed_counts()). Only
    # when gene_identifiers keeps some of the genes are the kept genes' rows
    # read once more for their sums.
    totals <- if (is.null(expression$genes)) {
      expression$stats$cell_totals
    } else {
      blocks <- .backed_blocks(matrix, sort(expression$genes))
      summed <- numeric(ncol(matrix))
      for (k in seq_len(blocks$n)) summed <- .add_block_to_cell_totals(summed, blocks$read(k))
      summed
    }
    return(if (is.null(expression$columns)) totals else totals[expression$columns])
  }
  if (expression$sparse) {
    if (expression$feature_major && !is.null(expression$genes)) {
      if (!.extent_holds_counts(.native_kept_extent(matrix, .expression_gene_map(expression)$map, whole = TRUE))) {
        return(NULL)
      }
      kept <- numeric(nrow(matrix))
      kept[expression$genes] <- 1
      return(as.numeric(Matrix::crossprod(matrix, kept)))
    }
    if (!.extent_holds_counts(.native_extent(matrix@x, whole = TRUE))) return(NULL)
    return(as.numeric(if (expression$feature_major) Matrix::colSums(matrix) else Matrix::rowSums(matrix)))
  }
  if (!.extent_holds_counts(.native_extent(matrix, whole = TRUE))) return(NULL)
  as.numeric(if (expression$feature_major) colSums(matrix) else rowSums(matrix))
}

# total_counts = checked for its shape: one obs column's name (an AnnData),
# or for several modalities each one's column, named by modality.
.total_counts_argument <- function(value, modalities) {
  if (is.null(value)) return(NULL)
  if (is.null(modalities)) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value) || !is.null(names(value))) {
      .abort(
        "total_counts names the obs column holding each cell's total count over all genes ",
        "(computed before any gene was selected), for example total_counts = \"nCount_RNA\"."
      )
    }
    return(value)
  }
  value <- if (is.list(value)) unlist(value) else value
  if (!is.character(value) || length(value) == 0L || is.null(names(value)) || anyNA(value) ||
      any(!nzchar(value)) || any(!nzchar(names(value)))) {
    .abort(
      "With several modalities, total_counts names, for each modality, the obs column of ",
      "its cells' total count over all its genes, for example total_counts = c(RNA = \"nCount_RNA\")."
    )
  }
  unknown <- setdiff(names(value), modalities)
  if (length(unknown) > 0L) {
    .abort("total_counts names ", .format_value_list(unknown), ", which are not modalities written: ",
           .format_value_list(modalities), ".")
  }
  value
}

# total_counts_over = for one named column: what it covers ("all genes" when
# NULL); refused without a named column, since the totals the writer sums
# itself are over the stored genes.
.total_counts_over_argument <- function(value, named) {
  if (is.null(value)) return(.total_counts_over_all)
  if (!is.character(value) || length(value) != 1L || is.na(value) || !value %in% .total_counts_over_named) {
    .abort(
      "total_counts_over must be \"all genes\" (every gene measured) or \"source panel\" (every gene of ",
      "the source's gene panel, more than are stored), got ", .format_value_list(value), "."
    )
  }
  if (!named) {
    .abort(
      "total_counts_over says what the column total_counts names covers, and no total_counts was ",
      "given; name the column of each cell's total count with total_counts =."
    )
  }
  value
}

# total_counts_over = with several modalities: one value for every column
# total_counts names, or a named character vector by modality (a modality
# left out: "all genes").
.total_counts_over_by_modality <- function(value, named) {
  modalities <- names(named)
  if (is.null(value) || (is.character(value) && length(value) == 1L && is.null(names(value)))) {
    over <- .total_counts_over_argument(value, named = length(modalities) > 0L)
    return(stats::setNames(rep(over, length(modalities)), modalities))
  }
  value <- if (is.list(value)) unlist(value) else value
  if (!is.character(value) || is.null(names(value)) || any(!nzchar(names(value)))) {
    .abort(
      "With several modalities, total_counts_over is \"all genes\" or \"source panel\" for every column ",
      "total_counts names, or one of them by modality, for example total_counts_over = c(RNA = \"source panel\")."
    )
  }
  unnamed <- setdiff(names(value), modalities)
  if (length(unnamed) > 0L) {
    .abort("total_counts_over names ", .format_value_list(sort(unnamed)), ", for which total_counts names no column: ",
           .format_value_list(sort(modalities)), ".")
  }
  stats::setNames(vapply(modalities, function(modality) {
    .total_counts_over_argument(if (modality %in% names(value)) unname(value[[modality]]), named = TRUE)
  }, character(1L)), modalities)
}

# A named total needs a matrix of counts.
.require_counts_total <- function(sums, has_matrix) {
  if (!has_matrix) {
    .abort("total_counts names each cell's total count of gene_expression, and no gene_expression was supplied.")
  }
  if (is.null(sums)) {
    .abort(
      "total_counts names each cell's total count, and the expression matrix does not hold ",
      "counts (a value is not a whole number of at least 0): the totals normalize counts, so ",
      "write the counts as X (for a Seurat object, its counts layer)."
    )
  }
  invisible(sums)
}

# A named column holds each cell's total, at least the stored genes' sum.
.require_totals_column <- function(values, name, sums, cell_names) {
  label <- paste0("total_counts names obs column '", name, "', which")
  if (!is.numeric(values) || is.factor(values)) .abort(label, " is not a column of numbers.")
  values <- as.numeric(values)
  bad <- which(!is.finite(values) | values < 0)
  if (length(bad) > 0L) {
    .abort(label, " has ", format(length(bad), big.mark = ","), " ", .plural(length(bad), "cell"),
           " with no total of 0 or more (first '", cell_names[[bad[[1L]]]], "': ", values[[bad[[1L]]]], ").")
  }
  short <- which(values < sums * (1 - .total_counts_tolerance) - .total_counts_tolerance)
  if (length(short) > 0L) {
    cell <- short[[1L]]
    .abort(
      label, " is below the sum of the genes stored for ", format(length(short), big.mark = ","), " ",
      .plural(length(short), "cell"), " (first '", cell_names[[cell]], "': ", format(values[[cell]]),
      " against ", format(sums[[cell]]), "), so it is not each cell's total count over all genes. ",
      "Name the column of total counts (Seurat's nCount_RNA), not a log of it."
    )
  }
  invisible(NULL)
}

.total_counts_record <- function(field, over) {
  list(field = .uns_string(field), over = .uns_string(over))
}

# An AnnData's obs with its cells' totals, and the record of them (NULL:
# none). `sums` are .expression_row_totals() in obs's row order.
.with_total_counts <- function(obs, sums, named, cell_names, has_matrix, over = .total_counts_over_all) {
  if (!is.null(named)) {
    sums <- .require_counts_total(sums, has_matrix)
    if (!named %in% names(obs)) {
      .abort("total_counts names obs column '", named, "', which obs does not hold: ",
             .format_value_list(sort(names(obs))), ".")
    }
    .require_totals_column(obs[[named]], named, sums, cell_names)
    return(list(obs = obs, record = .total_counts_record(named, over)))
  }
  if (is.null(sums)) return(list(obs = obs, record = NULL))
  taken <- names(obs)[tolower(names(obs)) == .total_counts_column]
  if (length(taken) > 0L) {
    # The data's own column: recorded only when it is exactly the stored genes' sums.
    own <- obs[[taken[[1L]]]]
    same <- identical(taken[[1L]], .total_counts_column) && is.numeric(own) && !is.factor(own) &&
      !anyNA(own) && all(as.numeric(own) == sums)
    return(list(obs = obs, record = if (same) .total_counts_record(.total_counts_column, .total_counts_over_stored)))
  }
  obs[[.total_counts_column]] <- sums
  list(obs = obs, record = .total_counts_record(.total_counts_column, .total_counts_over_stored))
}

# Each modality of a MuData with its cells' totals (`space$totals`, written
# as mod/<name>/obs/total_counts, shown as <name>:total_counts) or named by a
# column of the cells' obs, and the record of each, by modality.
.modality_total_counts <- function(spaces, names, named, obs, obs_index, over = NULL) {
  records <- list()
  for (index in seq_along(spaces)) {
    space <- spaces[[index]]
    name <- names[[index]]
    sums <- .expression_row_totals(space$expression)
    column <- if (is.null(named)) NULL else unname(named[name])
    if (length(column) == 1L && !is.na(column)) {
      tryCatch({
        sums <- .require_counts_total(sums, !is.null(space$expression))
        if (!column %in% names(obs)) {
          .abort("total_counts names obs column '", column, "', which obs does not hold.")
        }
        rows <- space$positions %||% seq_along(obs_index)
        .require_totals_column(obs[[column]][rows], column, sums, obs_index[rows])
      }, error = function(error) .abort("Modality '", name, "': ", conditionMessage(error)))
      records[[name]] <- .total_counts_record(column, if (name %in% names(over)) over[[name]] else .total_counts_over_all)
    } else if (!is.null(sums)) {
      spaces[[index]]$totals <- sums
      records[[name]] <- .total_counts_record(paste0(name, ":", .total_counts_column), .total_counts_over_stored)
    }
  }
  list(spaces = spaces, records = records)
}
