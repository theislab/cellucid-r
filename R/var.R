# The gene axis: naming the genes to export, in one or more name sets, and
# reading one gene's column out of whatever matrix the caller supplied.

# Validate one name-set source: NULL for rownames(var), or one var column.
.validate_gene_name_source <- function(column, arg) {
  if (
    !is.character(column) ||
      length(column) != 1L ||
      is.na(column) ||
      is.object(column) ||
      !nzchar(column)
  ) {
    stop(arg, " must be NULL or one non-empty string.", call. = FALSE)
  }
  .validate_display_text(column, arg)
}

# Resolve one complete gene-name set over every var row. `column = NULL`
# names genes by rownames(var). A column's empty values (NA or "") take the
# row's name from rownames(var), and repeated names are made unique with the
# rule cellucid-python and the web app share (see make-unique.R). Uniqueness
# is resolved over the whole var rather than the exported subset, so a gene
# keeps its name however the export is narrowed.
.resolve_gene_name_set <- function(var, column, arg) {
  automatic <- .row_names_info(var) < 0L
  if (is.null(column)) {
    # A data frame always has row names, so rownames(var) is never NULL and
    # asking whether it is answers nothing. What a caller who never set them
    # has is the automatic sequence, and rownames() materializes that as the
    # strings "1" to "N" -- unique, drawable, and accepted by every check
    # downstream. Selecting them would publish an export whose genes are named
    # after their own row numbers, so a wet-lab user searching CD8A finds
    # nothing and there is no failure anywhere to say why. .row_names_info()
    # is what tells the two apart: it returns the row count negated while the
    # row names are still the automatic ones. The cellucid Python package
    # refuses the same input, because a default RangeIndex yields integers and
    # its gene names must be strings.
    if (automatic) {
      stop(
        "var has only automatic row names, so rownames(var) would name the ",
        "genes '1' to '", nrow(var), "'. Set rownames(var) to the gene ",
        "identifiers, or pass var_gene_id_column.",
        call. = FALSE
      )
    }
    return(.make_unique_names(rownames(var)))
  }
  column <- .validate_gene_name_source(column, arg)
  if (!column %in% names(var)) {
    stop(
      arg, " '", column, "' not found in var. Available columns: ",
      .format_value_list(names(var)),
      call. = FALSE
    )
  }
  values <- var[[column]]
  if (is.factor(values)) {
    values <- as.character(values)
  } else if (is.object(values) || !is.null(dim(values))) {
    # A classed column carries meaning in attributes that a name list drops,
    # and identity downstream is decided with base semantics that compare
    # those attributes, so the rule every identifier axis applies is applied
    # here: a gene name column is a native character (or factor) vector.
    stop("Gene keys must be a native character vector.", call. = FALSE)
  }
  # Automatic row names are no names at all, so an empty value has none to
  # take, exactly as a Python RangeIndex supplies integers rather than names.
  index <- if (automatic) rep(NA_character_, nrow(var)) else rownames(var)
  tryCatch(
    .make_unique_names(values, index = index),
    cellucid_data_mapping_error = function(condition) {
      stop(
        "var column '", column, "': ", conditionMessage(condition),
        call. = FALSE
      )
    }
  )
}

# Describe one gene-name set over the exported genes, in exported order. The
# default set's names are the manifest's field names, so its `names` is null;
# an alternative set lists one name per exported gene, with null wherever it
# equals the default name.
.gene_name_set_record <- function(source, names, exported_positions, default_names = NULL) {
  renamed <- list()
  for (position in exported_positions) {
    if (!identical(names$names[[position]], names$filled[[position]])) {
      renamed[[length(renamed) + 1L]] <- list(
        names$filled[[position]],
        names$names[[position]]
      )
    }
  }
  listed <- NULL
  if (!is.null(default_names)) {
    listed <- lapply(exported_positions, function(position) {
      if (identical(names$names[[position]], default_names$names[[position]])) {
        NULL
      } else {
        names$names[[position]]
      }
    })
  }
  list(
    source = source,
    filled_from_index = as.integer(sum(names$from_index[exported_positions])),
    renamed = renamed,
    names = listed
  )
}

.gene_name_sets <- function(
    var,
    default_source,
    default_names,
    alternative_columns,
    exported_positions
) {
  sets <- list(.gene_name_set_record(
    source = default_source,
    names = default_names,
    exported_positions = exported_positions
  ))
  if (is.null(alternative_columns)) {
    return(sets)
  }
  if (
    !is.character(alternative_columns) ||
      is.object(alternative_columns) ||
      !is.null(dim(alternative_columns))
  ) {
    stop(
      "var_gene_name_columns must be a character vector of var column names.",
      call. = FALSE
    )
  }
  seen <- character(0)
  for (column in alternative_columns) {
    if (!is.null(default_source) && identical(column, default_source)) {
      stop(
        "var_gene_name_columns lists '", column, "', which var_gene_id_column ",
        "already uses for the default gene names.",
        call. = FALSE
      )
    }
    names <- .resolve_gene_name_set(var, column, "var_gene_name_columns entry")
    if (column %in% seen) {
      stop(
        "var_gene_name_columns lists '", column, "' more than once.",
        call. = FALSE
      )
    }
    seen <- c(seen, column)
    .require_field_identities(
      names$names[exported_positions],
      what = paste0("var column '", column, "' gene name")
    )
    sets[[length(sets) + 1L]] <- .gene_name_set_record(
      source = column,
      names = names,
      exported_positions = exported_positions,
      default_names = default_names
    )
  }
  sets
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
