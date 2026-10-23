# AnnData's on-disk encodings (the anndata 0.8+ element specification), each
# written as Zarr groups and arrays (format 2 or 3): every element carries its
# `encoding-type` and `encoding-version` attributes, exactly as Python's
# anndata writes them, so anndata reads the store without a warning.

.encoding <- function(type, version, ...) {
  c(list(...), list(`encoding-type` = type, `encoding-version` = version))
}

.ad_dict <- function(store, path) {
  .zarr_group(store, path, .encoding("dict", "0.1.0"))
}

.ad_string_array <- function(store, path, strings) {
  .zarr_vector(store, path, strings, "|O", .encoding("string-array", "0.2.0"))
}

.ad_numeric_array <- function(store, path, values, dtype) {
  .zarr_vector(store, path, values, dtype, .encoding("array", "0.2.0"))
}

# The signed integer type of the category codes, exactly as pandas chooses
# it (pandas.core.dtypes.cast.coerce_indexer_dtype): int8 for fewer than 127
# categories, int16 for fewer than 32,767, int32 for fewer than 2^31 - 1,
# int64 beyond. pandas compares the category count with each type's maximum
# strictly, so 127 categories already take int16.
.ad_codes_dtype <- function(n_categories) {
  if (n_categories < 127) {
    "|i1"
  } else if (n_categories < 32767) {
    "<i2"
  } else if (n_categories < 2147483647) {
    "<i4"
  } else {
    "<i8"
  }
}

.ad_categorical <- function(store, path, codes, categories, ordered) {
  .zarr_group(store, path, .encoding("categorical", "0.2.0", ordered = ordered))
  .ad_numeric_array(store, .zarr_key(path, "codes"), codes, .ad_codes_dtype(length(categories)))
  .ad_string_array(store, .zarr_key(path, "categories"), categories)
}

.ad_nullable <- function(store, path, type, values, mask, dtype) {
  .zarr_group(store, path, .encoding(type, "0.1.0"))
  .ad_numeric_array(store, .zarr_key(path, "values"), values, dtype)
  .ad_numeric_array(store, .zarr_key(path, "mask"), mask, "|b1")
}

# How one data-frame column is stored. The rules are anndata's own:
#
# - a factor is a categorical with its levels, in their order (`ordered` for
#   an ordered factor); a missing value is code -1;
# - a character column is a categorical (levels in code-point order) when a
#   value repeats or is missing, and a string array otherwise, which is what
#   anndata's strings_to_categoricals() does on write;
# - a logical column is a boolean array, or a nullable boolean when it holds
#   NA;
# - an integer column is int32 (R's one integer type), or a nullable integer
#   with int32 values when it holds NA;
# - a double column is float64, NA and NaN both written as NaN, whether or not
#   its values are whole numbers.
#
# Each column keeps its own type, as pandas and anndata keep it: this is the
# rule cellucid's Python prepare() follows, so the same column written from R
# and from Python is the same element.
#
# Anything else -- a date, a time, a list -- is refused with the conversion
# that makes it storable, because writing its underlying numbers would
# store values the caller never saw.
.ad_column <- function(store, path, values, label) {
  if (is.factor(values)) {
    codes <- as.integer(values) - 1L
    codes[is.na(codes)] <- -1L
    .ad_categorical(store, path, codes, levels(values), is.ordered(values))
    return(invisible(NULL))
  }
  if (is.object(values) || !is.null(dim(values)) || is.list(values)) {
    stop(
      label, " is a ", class(values)[[1L]], " column, which has no AnnData ",
      "encoding Cellucid writes. Convert it first, for example with ",
      "as.character() or as.numeric().",
      call. = FALSE
    )
  }
  if (is.character(values)) {
    if (anyNA(values) || anyDuplicated(values)) {
      categories <- sort(unique(values[!is.na(values)]), method = "radix")
      codes <- match(values, categories) - 1L
      codes[is.na(codes)] <- -1L
      .ad_categorical(store, path, codes, categories, FALSE)
    } else {
      .ad_string_array(store, path, values)
    }
    return(invisible(NULL))
  }
  if (is.logical(values)) {
    if (anyNA(values)) {
      mask <- is.na(values)
      filled <- values
      filled[mask] <- FALSE
      .ad_nullable(store, path, "nullable-boolean", filled, mask, "|b1")
    } else {
      .ad_numeric_array(store, path, values, "|b1")
    }
    return(invisible(NULL))
  }
  if (is.integer(values)) {
    if (anyNA(values)) {
      mask <- is.na(values)
      filled <- values
      filled[mask] <- 0L
      .ad_nullable(store, path, "nullable-integer", filled, mask, "<i4")
    } else {
      .ad_numeric_array(store, path, values, "<i4")
    }
    return(invisible(NULL))
  }
  if (is.double(values)) {
    infinite <- which(is.infinite(values))
    if (length(infinite) > 0L) {
      stop(
        label, " holds ", .format_count(length(infinite)), " infinite ",
        .plural(length(infinite), "value", "values"), " (first at row ",
        infinite[[1L]], "). A colour scale has no position for an infinity; ",
        "replace it with NA, which is written as a missing value.",
        call. = FALSE
      )
    }
    .ad_numeric_array(store, path, values, "<f8")
    return(invisible(NULL))
  }
  stop(
    label, " is of type ", typeof(values), ", which has no AnnData encoding ",
    "Cellucid writes. Convert it first, for example with as.character() or ",
    "as.numeric().",
    call. = FALSE
  )
}

# A data frame: the index as a string array named `_index`, and one element
# per column, in `column-order`.
.ad_dataframe <- function(store, path, frame, index, axis) {
  columns <- names(frame)
  .zarr_group(store, path, .encoding(
    "dataframe", "0.2.0",
    `_index` = "_index",
    `column-order` = I(columns)
  ))
  .ad_string_array(store, .zarr_key(path, "_index"), index)
  for (column in columns) {
    .ad_column(
      store,
      .zarr_key(path, column),
      frame[[column]],
      paste0(axis, " column '", column, "'")
    )
  }
  invisible(NULL)
}

# A compressed sparse matrix. `pointers`, `indices` and `values` are the
# 0-based compressed arrays; indices and pointers are int32 unless the number
# of stored values reaches 2^31, when both are int64. `values` are stored in
# `dtype`. A `sliced` matrix -- the expression matrix, which the viewer reads
# one gene's range at a time -- chunks its data and indices in 2^16 values;
# any other (the neighbour graph, read whole) in 2^18.
.ad_sparse <- function(store, path, format, shape, pointers, indices, values,
                       dtype = "<f4", sliced = FALSE) {
  index_dtype <- if (length(values) >= 2147483648) "<i8" else "<i4"
  target <- if (sliced) .chunk_elements_sliced else .chunk_elements_whole
  .zarr_group(store, path, .encoding(format, "0.1.0", shape = I(as.numeric(shape))))
  # The expression values (sliced) are not shuffled (zarr-store.R).
  .zarr_vector(store, .zarr_key(path, "data"), values, dtype, target = target, shuffle = !sliced)
  .zarr_vector(store, .zarr_key(path, "indices"), indices, index_dtype, target = target)
  .zarr_vector(store, .zarr_key(path, "indptr"), pointers, index_dtype)
}

# A dense floating matrix with one cell per row, float32 unless `dtype` says
# otherwise. An embedding is read whole, so a chunk holds every column of a
# block of rows; dense expression is read one gene at a time, so a chunk
# holds a tall strip of few columns.
.ad_dense <- function(store, path, matrix, gene_major = FALSE, dtype = "<f4") {
  n_rows <- nrow(matrix)
  n_cols <- ncol(matrix)
  if (gene_major) {
    chunk_rows <- .chunk_length(n_rows, .chunk_elements_whole)
    chunk_cols <- .chunk_length(n_cols, floor(.chunk_elements_whole / chunk_rows))
  } else {
    chunk_cols <- max(1, n_cols)
    chunk_rows <- .chunk_length(n_rows, floor(.chunk_elements_whole / chunk_cols))
  }
  # Dense expression values are not shuffled (zarr-store.R).
  .zarr_matrix(store, path, matrix, dtype, chunk_rows, chunk_cols, .encoding("array", "0.2.0"),
    shuffle = !gene_major
  )
}

# `uns` values. A named list is a dict; every leaf is built with one of the
# constructors below, so a one-element array is never mistaken for a scalar.
# A leaf is told apart by its class, never by a member name: a dict keyed by
# obs column names (outlierQuantiles$fields) can hold any name.
.uns_leaf <- function(kind, value) structure(list(kind = kind, value = value), class = "cellucid_uns_leaf")
.uns_string <- function(value) .uns_leaf("string", value)
.uns_integer <- function(value) .uns_leaf("integer", value)
.uns_boolean <- function(value) .uns_leaf("boolean", value)
.uns_strings <- function(value) .uns_leaf("strings", value)
.uns_integers <- function(value) .uns_leaf("integers", value)
.uns_floats32 <- function(value) .uns_leaf("floats32", value)

.ad_uns_value <- function(store, path, value) {
  if (!inherits(value, "cellucid_uns_leaf")) {
    .ad_dict(store, path)
    for (name in names(value)) {
      .ad_uns_value(store, .zarr_key(path, name), value[[name]])
    }
    return(invisible(NULL))
  }
  kind <- value$kind
  switch(kind,
    string = .zarr_scalar(store, path, value$value, "|O", .encoding("string", "0.2.0")),
    integer = .zarr_scalar(store, path, value$value, "<i8", .encoding("numeric-scalar", "0.2.0")),
    boolean = .zarr_scalar(store, path, value$value, "|b1", .encoding("numeric-scalar", "0.2.0")),
    strings = .ad_string_array(store, path, value$value),
    integers = .ad_numeric_array(store, path, value$value, "<i8"),
    floats32 = .ad_numeric_array(store, path, value$value, "<f4")
  )
  invisible(NULL)
}
