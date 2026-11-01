# cellucid_transition_drift(): the drift T %*% E - E, the counterpart of
# cellucid-python's compute_transition_drift().
#
# tests/testthat/fixtures/transition_drift_conformance.json is a byte-identical
# copy of cellucid-python's tests/fixtures/transition_drift_conformance.json,
# written by its scripts/generate_transition_drift_conformance.py with exact
# rational arithmetic, independently of both packages; both repositories pin
# the same SHA-256 and hold their implementation to it. An exact case (its
# product and row sums are exact in double precision in any summation order)
# must be reproduced bit for bit; a general case within the per-value bound
# the fixture states, twice the forward-error bound of the computation in any
# summation order, because R's %*% and rowSums() and numpy add in different
# orders.

.transition_drift_sha256 <- "0c3a15648c21a6441f2e530d25a2837298dacd5113942f3ef280b048a1ea81e8"

drift_case_matrix <- function(values, n_rows, n_cols) {
  matrix(vapply(values, as.numeric, numeric(1)), nrow = n_rows, ncol = n_cols, byrow = TRUE)
}

drift_case_transition <- function(case) {
  n <- case$nCells
  transition <- case$transition
  if (identical(transition$kind, "dense")) {
    return(drift_case_matrix(transition$values, n, n))
  }
  dense <- matrix(0, nrow = n, ncol = n)
  rows <- vapply(transition$rows, as.integer, integer(1)) + 1L
  columns <- vapply(transition$columns, as.integer, integer(1)) + 1L
  dense[cbind(rows, columns)] <- vapply(transition$values, as.numeric, numeric(1))
  dense
}

test_that("the shared conformance cases give the pinned drift, dense and sparse", {
  skip_if_not_installed("Matrix")
  expect_pinned_fixture("transition_drift_conformance.json", .transition_drift_sha256)
  document <- cellucid_read_fixture("transition_drift_conformance.json")
  expect_identical(document$version, 1L)
  cases <- document$cases
  expect_identical(sum(vapply(cases, function(case) isTRUE(case$exact), logical(1))), 5L)
  expect_length(cases, 9L)
  for (case in cases) {
    n <- case$nCells
    width <- case$width
    embedding <- drift_case_matrix(case$embedding, n, width)
    expected <- drift_case_matrix(case$drift, n, width)
    bound <- drift_case_matrix(case$bound, n, width)
    dense <- drift_case_transition(case)
    representations <- list(
      dense = dense,
      sparse = methods::as(Matrix::Matrix(dense, sparse = TRUE), "generalMatrix")
    )
    for (representation in names(representations)) {
      label <- paste0(case$name, " (", representation, ")")
      drift <- cellucid_transition_drift(
        representations[[representation]], embedding,
        normalize_rows = case$normalizeRows
      )
      expect_identical(typeof(drift), "double", label = label)
      expect_identical(dim(drift), c(n, width), label = label)
      if (isTRUE(case$exact)) {
        expect_identical(drift, expected, label = label)
      } else {
        expect_true(all(abs(drift - expected) <= bound), label = label)
      }
    }
  }
})

test_that("the drift is the transition-weighted mean position minus the cell's own", {
  embedding <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  transition <- matrix(c(0, 1, 1, 0, 1, 0, 0, 0, 2), nrow = 3, byrow = TRUE)
  expect_identical(
    cellucid_transition_drift(transition, embedding, normalize_rows = TRUE),
    matrix(c(0.5, 0.5, 0, 0, 0, 0), ncol = 2, byrow = TRUE)
  )
  # As given, the third row (weight 2) moves the cell to twice its position.
  expect_identical(
    cellucid_transition_drift(transition, embedding, normalize_rows = FALSE),
    matrix(c(1, 1, 0, 0, 0, 1), ncol = 2, byrow = TRUE)
  )
  # The identity, scaled or not, does not move any cell.
  expect_identical(
    cellucid_transition_drift(diag(3), embedding, normalize_rows = FALSE),
    matrix(0, nrow = 3, ncol = 2)
  )
  expect_identical(
    cellucid_transition_drift(3 * diag(3), embedding, normalize_rows = TRUE),
    matrix(0, nrow = 3, ncol = 2)
  )
})

test_that("every input form cellucid_prepare() takes for an embedding is accepted", {
  skip_if_not_installed("Matrix")
  embedding <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE,
                      dimnames = list(c("a", "b", "c"), c("UMAP_1", "UMAP_2")))
  transition <- matrix(c(0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L), nrow = 3, byrow = TRUE)
  expected <- cellucid_transition_drift(transition, embedding, normalize_rows = TRUE)
  expect_identical(dimnames(expected), dimnames(embedding))
  unnamed <- unname(expected)
  expect_identical(
    unname(cellucid_transition_drift(transition, as.data.frame(embedding), normalize_rows = TRUE)),
    unnamed
  )
  expect_identical(
    unname(cellucid_transition_drift(transition, Matrix::Matrix(embedding), normalize_rows = TRUE)),
    unnamed
  )
  for (form in list(
    Matrix::Matrix(transition * 1, sparse = TRUE),
    Matrix::Matrix(transition * 1, sparse = FALSE),
    methods::as(Matrix::Matrix(transition * 1, sparse = TRUE), "TsparseMatrix")
  )) {
    expect_identical(cellucid_transition_drift(form, embedding, normalize_rows = TRUE), expected)
  }
  expect_identical(
    cellucid_transition_drift(Matrix::Diagonal(3, 2), embedding, normalize_rows = TRUE),
    embedding - embedding
  )
  # A plain vector is the one column it declares, as for an embedding.
  expect_identical(
    cellucid_transition_drift(transition, c(0, 2, 4), normalize_rows = TRUE),
    matrix(c(3, 0, 0), ncol = 1L)
  )
  # Any width: every column is used.
  wide <- matrix(seq_len(30) / 4, nrow = 3)
  expect_identical(dim(cellucid_transition_drift(transition, wide, normalize_rows = TRUE)), c(3L, 10L))
})

test_that("the drift is a vector field cellucid_prepare() stores as <field>_<basis>", {
  embedding <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  transition <- matrix(c(0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1), nrow = 4, byrow = TRUE)
  drift <- cellucid_transition_drift(transition, embedding, normalize_rows = TRUE)
  store <- cellucid_prepare(
    priority_seed = NULL,
    out_path = file.path(tempfile("cellucid_drift_"), "drift.zarr"),
    format = "zarr",
    embeddings = list(X_umap = embedding),
    dataset_id = "drift",
    dataset_name = "Drift",
    vector_fields = list(T_fwd = list(X_umap = drift))
  )
  expect_identical(array_meta(store, "obsm/T_fwd_umap")$dtype, "<f4")
  expect_identical(store_array(store, "obsm/T_fwd_umap"), matrix(as_float32(drift), nrow = 4L))
})

test_that("invalid input is refused with the reason", {
  skip_if_not_installed("Matrix")
  embedding <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  expect_error(
    cellucid_transition_drift(diag(3), embedding),
    "normalize_rows is required: TRUE divides every row", fixed = TRUE
  )
  expect_error(cellucid_transition_drift(embedding = embedding, normalize_rows = TRUE), "transition_matrix is required")
  expect_error(cellucid_transition_drift(diag(3), normalize_rows = TRUE), "embedding is required")
  for (flag in list(NA, 1, "yes", c(TRUE, FALSE))) {
    expect_error(
      cellucid_transition_drift(diag(3), embedding, normalize_rows = flag),
      "normalize_rows must be exactly TRUE or FALSE.", fixed = TRUE
    )
  }
  expect_error(
    cellucid_transition_drift(diag(4), embedding, normalize_rows = FALSE),
    "transition_matrix must be 3 x 3 (one row and one column per row of embedding), got 4 x 4.",
    fixed = TRUE
  )
  expect_error(
    cellucid_transition_drift(diag(3)[, 1:2], embedding, normalize_rows = FALSE),
    "got 3 x 2", fixed = TRUE
  )
  expect_error(
    cellucid_transition_drift(diag(3) == 1, embedding, normalize_rows = FALSE),
    "transition_matrix must be a real numeric matrix or a Matrix::Matrix matrix.", fixed = TRUE
  )
  expect_error(
    cellucid_transition_drift(Matrix::Matrix(diag(3) == 1, sparse = TRUE), embedding, normalize_rows = FALSE),
    "transition_matrix must contain real numeric values.", fixed = TRUE
  )
  expect_error(
    cellucid_transition_drift(as.data.frame(diag(3)), embedding, normalize_rows = FALSE),
    "must be a real numeric matrix", fixed = TRUE
  )
  negative <- diag(3)
  negative[1, 2] <- -0.5
  expect_error(
    cellucid_transition_drift(negative, embedding, normalize_rows = FALSE),
    "transition_matrix values must be non-negative.", fixed = TRUE
  )
  expect_error(
    cellucid_transition_drift(Matrix::Matrix(negative, sparse = TRUE), embedding, normalize_rows = FALSE),
    "transition_matrix values must be non-negative.", fixed = TRUE
  )
  for (bad in c(NA, NaN, Inf)) {
    broken <- diag(3)
    broken[2, 3] <- bad
    expect_error(
      cellucid_transition_drift(broken, embedding, normalize_rows = FALSE),
      "transition_matrix must contain only finite values.", fixed = TRUE
    )
  }
  # A zero row is refused, never replaced, when rows are normalized; as given
  # it is a cell that moves to the origin.
  zero_row <- diag(3)
  zero_row[2, 2] <- 0
  expect_error(
    cellucid_transition_drift(zero_row, embedding, normalize_rows = TRUE),
    "Every transition_matrix row must have a nonzero finite row sum when normalize_rows = TRUE.",
    fixed = TRUE
  )
  expect_identical(
    cellucid_transition_drift(zero_row, embedding, normalize_rows = FALSE)[2, ],
    c(-1, 0)
  )
  expect_error(
    cellucid_transition_drift(diag(3), matrix(c(0, Inf, 1, 0, 0, 1), ncol = 2), normalize_rows = FALSE),
    "embedding must contain only finite values.", fixed = TRUE
  )
  expect_error(
    cellucid_transition_drift(diag(3), matrix(c(0, 1e39, 1, 0, 0, 1), ncol = 2), normalize_rows = FALSE),
    "beyond the float32 range", fixed = TRUE
  )
  expect_error(
    cellucid_transition_drift(diag(3), matrix(c(0, 2^-151, 1, 0, 0, 1), ncol = 2), normalize_rows = FALSE),
    "embedding holds 1 nonzero value below the smallest float32 (magnitude at most 2^-150), which would be stored as zero.",
    fixed = TRUE
  )
  expect_error(
    cellucid_transition_drift(diag(3), matrix(character(6), ncol = 2), normalize_rows = FALSE),
    "embedding must be a real numeric matrix or vector.", fixed = TRUE
  )
  expect_error(
    cellucid_transition_drift(matrix(0, 0, 0), matrix(numeric(0), ncol = 2), normalize_rows = FALSE),
    "embedding must contain at least one row.", fixed = TRUE
  )
  # A drift beyond the float32 range is refused, not stored as an infinity.
  expect_error(
    cellucid_transition_drift(matrix(c(1e10, 0, 0, 1), 2), matrix(c(3e30, 0), ncol = 1), normalize_rows = FALSE),
    "The computed transition drift holds values beyond the float32 range", fixed = TRUE
  )
  # And a drift that float32 would store as zero.
  expect_error(
    cellucid_transition_drift(matrix(c(1 + 2^-52, 0, 0, 1), 2), matrix(c(2^-126, 0), ncol = 1), normalize_rows = FALSE),
    "The computed transition drift holds 1 nonzero value below the smallest float32", fixed = TRUE
  )
})
