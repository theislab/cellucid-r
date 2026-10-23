# Turning what an R caller passes for a matrix -- a base matrix, a data
# frame of numeric columns, or a Matrix::Matrix -- into the storage this
# package computes on.

# A data frame stands in for a numeric matrix when every one of its columns is
# a native numeric vector, because as.matrix() then produces exactly the matrix
# the caller would have built by hand. A classed column does not qualify: it
# carries meaning in its attributes that the matrix drops, so what the store
# would hold is not what the caller passed.
#
# The two terms catch different classes, and the second is the one that earns
# its place. `is.numeric()` is an internal generic, and base R answers FALSE
# for Date, POSIXct, difftime, and factor, so those never reach the second
# term at all. What reaches it is every classed numeric that has no
# is.numeric() method -- units::units, bit64::integer64, and any S3 class a
# caller wrote -- where the attribute is the only thing that says what the
# numbers mean. as.matrix() drops it silently: a data frame of one metre
# column and one kilometre column becomes a matrix of equal numbers, and the
# store holds a 45-degree vector where the caller's data points 0.06
# degrees off axis. There is no later check that can notice, because by then
# the values are ordinary finite doubles.
#
# This is the one rule for every input that accepts a data frame in place of
# a matrix -- gene_expression, the embeddings, and the vector fields -- so the
# same input cannot be refused on one of them and accepted on another.
.is_numeric_data_frame <- function(x) {
  is.data.frame(x) &&
    all(vapply(x, function(column) {
      is.numeric(column) && !is.object(column)
    }, logical(1)))
}

# The largest finite IEEE 754 binary32 value. Every matrix is stored as
# float32 (expression may be float64, but must fit the same range), and a
# larger magnitude would be stored as an infinity.
.float32_max_magnitude <- (2 - 2^-23) * 2^127

.require_finite_float32 <- function(values, name, positions = TRUE) {
  if (anyNA(values) || any(!is.finite(values))) {
    stop(
      name, " must contain only finite values. ",
      .describe_non_finite(values, name, positions = positions),
      call. = FALSE
    )
  }
  if (length(values) > 0L && max(abs(range(values))) > .float32_max_magnitude) {
    stop(
      name, " holds values beyond the float32 range (magnitude above ",
      "3.4e38), which would be stored as infinities.",
      call. = FALSE
    )
  }
  invisible(values)
}

# A numeric matrix with one row per cell: an embedding or a vector field. A
# plain numeric vector is the one column it declares.
.as_cell_matrix <- function(x, name) {
  if (is.numeric(x) && !is.object(x) && is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  } else if (.is_numeric_data_frame(x) || inherits(x, "Matrix")) {
    x <- as.matrix(x)
  } else if (!is.matrix(x) || !is.numeric(x) || is.object(x)) {
    stop(name, " must be a real numeric matrix or vector.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  .require_finite_float32(x, name)
  x
}

# The largest integer float32 holds exactly. Integer counts above it are
# stored as float64, so every count reads back as it was given.
.float32_exact_integer_limit <- 2^24

# The dtype the expression matrix is stored in: float32, unless every value
# is a whole number and one exceeds 2^24 in magnitude, which float32 cannot
# hold exactly -- then the whole matrix is float64. `values` are the stored
# values (a sparse matrix's non-zeros, a dense matrix's every value) of the
# genes that are stored, after gene_identifiers has chosen them, so a large
# count in a gene left out does not widen the genes kept. This is the rule of
# the Python writer (stored_expression in prepare_data/_inputs.py) for the
# types R has: an integer matrix is whole numbers by type, a double matrix
# when every value is, and a logical one never exceeds 1.
.expression_dtype <- function(values) {
  if (length(values) == 0L || max(abs(range(values))) <= .float32_exact_integer_limit) {
    return("<f4")
  }
  if (all(values == floor(values))) "<f8" else "<f4"
}

# The expression matrix, cells x genes, as either a column-compressed sparse
# matrix (dgCMatrix: its columns are genes, which is exactly AnnData's
# gene-major csc_matrix) or a dense base matrix of doubles. A stored zero of a
# sparse matrix is no value and is dropped.
.as_expression_matrix <- function(x, name) {
  if (inherits(x, "sparseMatrix")) {
    if (inherits(x, "nMatrix")) {
      stop(name, " must hold numbers, not a pattern of non-zeros.", call. = FALSE)
    }
    x <- methods::as(methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    .require_finite_float32(x@x, name, positions = FALSE)
    return(Matrix::drop0(x))
  }
  if (inherits(x, "Matrix") || .is_numeric_data_frame(x)) {
    x <- as.matrix(x)
  } else if (!is.matrix(x) || !(is.numeric(x) || is.logical(x)) || is.object(x)) {
    stop(name, " must be a real numeric matrix or a Matrix::Matrix.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  .require_finite_float32(x, name)
  x
}
