# Transition drift: a vector field computed from a cell-cell transition
# matrix (CellRank's, for example), the counterpart of cellucid-python's
# compute_transition_drift(). The help page is hand-written in
# man/cellucid_transition_drift.Rd, as every page of this package is.
#
# Python's add_transition_drift_to_obsm() has no counterpart here, by design:
# R holds no AnnData in memory to add an obsm entry to. The drift is passed to
# cellucid_prepare(vector_fields = list(T_fwd = list(X_umap = drift))), which
# stores it as obsm[['T_fwd_umap']] and checks that key by the same
# vector-field key rule (.resolve_vector_fields()).
#
# The arithmetic is Python's, operation for operation: the product T %*% E in
# double precision, divided by the row sums of T when normalize_rows is TRUE,
# minus E. tests/testthat/fixtures/transition_drift_conformance.json, shared
# byte for byte with cellucid-python, pins the result.
cellucid_transition_drift <- function(transition_matrix, embedding, normalize_rows) {
  if (missing(transition_matrix)) {
    .abort("transition_matrix is required: a cells x cells transition matrix.")
  }
  if (missing(embedding)) {
    .abort("embedding is required: the cells x dimensions coordinates the drift is computed in.")
  }
  if (missing(normalize_rows)) {
    .abort(
      "normalize_rows is required: TRUE divides every row of transition_matrix ",
      "by its sum, FALSE uses it as given."
    )
  }
  normalize_rows <- .validate_flag(normalize_rows, "normalize_rows")
  coordinates <- .as_cell_matrix(embedding, "embedding")
  if (nrow(coordinates) < 1L) {
    .abort("embedding must contain at least one row.")
  }
  if (ncol(coordinates) < 1L) {
    .abort("embedding must have at least one column.")
  }
  .require_no_float32_underflow(coordinates, "embedding")
  transition <- .as_transition_matrix(transition_matrix, nrow(coordinates))

  product <- as.matrix(transition %*% coordinates)
  if (normalize_rows) {
    row_sums <- if (inherits(transition, "Matrix")) {
      Matrix::rowSums(transition)
    } else {
      rowSums(transition)
    }
    if (any(!is.finite(row_sums)) || any(row_sums == 0)) {
      .abort(
        "Every transition_matrix row must have a nonzero finite row sum ",
        "when normalize_rows = TRUE."
      )
    }
    product <- product / row_sums
  }
  drift <- product - coordinates
  .require_finite_float32(drift, "The computed transition drift")
  .require_no_float32_underflow(drift, "The computed transition drift")
  dimnames(drift) <- dimnames(coordinates)
  drift
}

# The transition matrix as the product needs it: a dgCMatrix for a Matrix
# sparse matrix, a base double matrix otherwise. Every value must be a finite,
# non-negative real number; logical and pattern matrices are refused, as
# Python refuses a boolean one.
.as_transition_matrix <- function(x, n_cells) {
  if (inherits(x, "sparseMatrix")) {
    if (!inherits(x, "dMatrix")) {
      .abort("transition_matrix must contain real numeric values.")
    }
    matrix <- methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix")
    values <- matrix@x
    positions <- FALSE
  } else {
    if (inherits(x, "Matrix")) {
      if (!inherits(x, "dMatrix")) {
        .abort("transition_matrix must contain real numeric values.")
      }
      x <- as.matrix(x)
    } else if (!is.matrix(x) || !is.numeric(x) || is.object(x)) {
      .abort(
        "transition_matrix must be a real numeric matrix or a Matrix::Matrix ",
        "matrix."
      )
    }
    matrix <- x
    storage.mode(matrix) <- "double"
    values <- matrix
    positions <- TRUE
  }
  if (nrow(matrix) != n_cells || ncol(matrix) != n_cells) {
    .abort(
      "transition_matrix must be ", n_cells, " x ", n_cells, " (one row and ",
      "one column per row of embedding), got ", nrow(matrix), " x ",
      ncol(matrix), "."
    )
  }
  if (anyNA(values) || any(!is.finite(values))) {
    .abort(
      "transition_matrix must contain only finite values. ",
      .describe_non_finite(values, "transition_matrix", positions = positions)
    )
  }
  if (any(values < 0)) {
    .abort("transition_matrix values must be non-negative.")
  }
  matrix
}

# Half the smallest float32 subnormal (2^-149): a nonzero magnitude at or
# below it rounds to zero in float32, so the stored value would no longer be
# the one computed.
.float32_underflow_magnitude <- 2^-150

.require_no_float32_underflow <- function(values, name) {
  underflow <- values != 0 & abs(values) <= .float32_underflow_magnitude
  if (any(underflow)) {
    .abort(
      name, " holds ", .format_count(sum(underflow)), " nonzero ",
      .plural(sum(underflow), "value", "values"), " below the smallest float32 ",
      "(magnitude at most 2^-150), which would be stored as zero."
    )
  }
  invisible(values)
}
