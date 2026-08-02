# Turning what an R caller passes for a matrix -- a base matrix, a data
# frame of numeric columns, or a Matrix::Matrix -- into the storage this
# package computes on.

# A data frame stands in for a numeric matrix when every one of its columns is
# a native numeric vector, because as.matrix() then produces exactly the matrix
# the caller would have built by hand. A classed column does not qualify: it
# carries meaning in its attributes that the matrix drops, so what the export
# would store is not what the caller passed.
#
# The two terms catch different classes, and the second is the one that earns
# its place. `is.numeric()` is an internal generic, and base R answers FALSE
# for Date, POSIXct, difftime, and factor, so those never reach the second
# term at all. What reaches it is every classed numeric that has no
# is.numeric() method -- units::units, bit64::integer64, and any S3 class a
# caller wrote -- where the attribute is the only thing that says what the
# numbers mean. as.matrix() drops it silently: a data frame of one metre
# column and one kilometre column becomes a matrix of equal numbers, and the
# export publishes a 45-degree vector where the caller's data points 0.06
# degrees off axis. There is no later check that can notice, because by then
# the values are ordinary finite doubles.
#
# This is the one rule for every axis that accepts a data frame in place of a
# matrix -- latent_space, gene_expression, the embeddings, and the vector
# fields -- so the same input cannot be refused on one of them and accepted on
# another.
.is_numeric_data_frame <- function(x) {
  is.data.frame(x) &&
    all(vapply(x, function(column) {
      is.numeric(column) && !is.object(column)
    }, logical(1)))
}

.as_dense_matrix <- function(x, name) {
  if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
  } else if (.is_numeric_data_frame(x)) {
    x <- as.matrix(x)
  } else if (
    !is.matrix(x) ||
      !is.numeric(x) ||
      is.object(x)
  ) {
    stop(name, " must be a finite real numeric matrix.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  x
}

.as_matrix_like <- function(x, name) {
  if (inherits(x, "Matrix")) {
    if (!inherits(x, "dMatrix")) {
      stop(name, " must contain real numeric values.", call. = FALSE)
    }
    values <- x@x
    if (anyNA(values) || any(!is.finite(values))) {
      stop(
        name,
        " must not contain infinite or missing sparse values.",
        call. = FALSE
      )
    }
    return(x)
  } else if (.is_numeric_data_frame(x)) {
    x <- as.matrix(x)
  } else if (
    !is.matrix(x) ||
      !is.numeric(x) ||
      is.object(x)
  ) {
    stop(name, " must be a real numeric matrix.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}
