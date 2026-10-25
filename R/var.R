# The gene axis: which var column names the genes in the viewer, and the
# names it gives them.

# Validate the name of the var column that names the genes.
.validate_gene_name_source <- function(column, arg) {
  if (
    !is.character(column) ||
      length(column) != 1L ||
      is.na(column) ||
      is.object(column) ||
      !nzchar(column)
  ) {
    .abort(arg, " must be NULL or one non-empty string.")
  }
  .validate_display_text(column, arg)
}

# Resolve the gene names a var column gives, over every var row. A column's
# empty values (NA or "") take the row's name from rownames(var), and
# repeated names are made unique with the rule the web app applies when it
# reads the column (see make-unique.R). Uniqueness is resolved over the whole
# var rather than the written subset, so a gene keeps its name however
# gene_identifiers narrows the store.
.resolve_gene_name_set <- function(var, column, arg) {
  automatic <- .row_names_info(var) < 0L
  column <- .validate_gene_name_source(column, arg)
  if (!column %in% names(var)) {
    .abort(
      arg, " '", column, "' not found in var. Available columns: ",
      .format_value_list(names(var)), "."
    )
  }
  values <- var[[column]]
  # The rule the viewer applies to a gene-name column: text and categories as
  # they are, integers as their decimal digits; a floating-point or true/false
  # column cannot name genes.
  if (is.factor(values)) {
    values <- as.character(values)
  } else if (is.object(values) || !is.null(dim(values))) {
    # A classed column carries meaning in attributes that a name list drops,
    # and identity downstream is decided with base semantics that compare
    # those attributes, so the rule every identifier axis applies is applied
    # here: a gene name column is a native vector.
    .abort("Gene keys must be a native character vector.")
  } else if (is.integer(values)) {
    values <- as.character(values)
  } else if (is.double(values)) {
    .abort("var column '", column, "': floating-point values cannot name genes.")
  } else if (is.logical(values)) {
    # NA is R's untyped missing value: a column of nothing but NA names no
    # gene, and every row takes its row name.
    if (!all(is.na(values))) {
      .abort("var column '", column, "': true/false values cannot name genes.")
    }
    values <- as.character(values)
  }
  # Automatic row names are no names at all, so an empty value has none to
  # take, exactly as a Python RangeIndex supplies integers rather than names.
  index <- if (automatic) rep(NA_character_, nrow(var)) else rownames(var)
  tryCatch(
    .make_unique_names(values, index = index),
    cellucid_data_mapping_error = function(condition) {
      .abort(
        "var column '", column, "': ", conditionMessage(condition)
      )
    }
  )
}
