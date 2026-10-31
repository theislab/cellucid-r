# The dataset plan where an object's layout departs from the simple case: one
# assay holding its cells in another order or only some of them, a factor
# with an NA level, a dense genes x cells assay, and a matrix held outside
# R's memory read in many blocks.

# cellucid_prepare(data = ) of an object whose extraction gives `components`
# (as .dataset_from() would for a Seurat or SingleCellExperiment object).
prepare_components <- function(components) {
  local_mocked_bindings(.dataset_from = function(data) components, .package = "cellucid")
  out <- file.path(tempfile(), "plan.zarr")
  suppressMessages(cellucid_prepare(out, "zarr", data = "an object", priority_seed = NULL))
}

# X of a store as a dense cells x genes matrix, without dimnames (X stores
# none). unname(): Matrix before 1.5 gives as.matrix() of an unnamed sparse
# matrix dimnames list(NULL, NULL), later releases none.
store_x <- function(store) {
  shape <- jsonlite::fromJSON(file.path(store, "X", ".zattrs"))$shape
  p <- store_array(store, "X/indptr")
  i <- store_array(store, "X/indices")
  x <- store_array(store, "X/data")
  unname(as.matrix(Matrix::sparseMatrix(i = i, p = p, x = x, dims = shape, index1 = FALSE)))
}

one_assay_components <- function(assay_cells) {
  cells <- paste0("c", 1:6)
  genes <- c("g1", "g2")
  # genes x cells, its columns in the assay's own cell order: the value of
  # gene g in cell cK is K (+ 100 for g2).
  values <- rbind(as.numeric(sub("c", "", assay_cells)), as.numeric(sub("c", "", assay_cells)) + 100)
  dimnames(values) <- list(genes, assay_cells)
  list(
    dataset_id = "plan", dataset_name = "Plan",
    embeddings = list(X_umap = matrix(seq_len(12) / 10, 6, 2, dimnames = list(cells, NULL))),
    obs = data.frame(row.names = cells, group = rep(c("a", "b"), 3)),
    modalities = list(RNA = list(
      name = "RNA",
      gene_expression = cellucid:::.feature_major(Matrix::Matrix(values, sparse = TRUE)),
      var = data.frame(row.names = genes),
      cells = if (identical(assay_cells, cells)) NULL else assay_cells
    )),
    notes = character(),
    kind = "Seurat"
  )
}

test_that("one assay holding the cells in another order is written in the cells' order", {
  skip_if_not_installed("Matrix")
  store <- prepare_components(one_assay_components(paste0("c", 6:1)))
  expect_identical(store_array(store, "obs/_index"), paste0("c", 1:6))
  x <- store_x(store)
  expect_identical(x[, 1], as.double(1:6))
  expect_identical(x[, 2], as.double(101:106))
})

test_that("one assay holding only some of the cells is refused with what to do", {
  skip_if_not_installed("Matrix")
  expect_error(
    prepare_components(one_assay_components(paste0("c", c(2, 4, 6, 1)))),
    "Modality 'RNA' holds 4 of the dataset's 6 cells.*subset\\(data, cells = Cells\\(data\\[\\[\"RNA\"\\]\\]\\)\\).*modalities = c\\("
  )
})

test_that("a factor's NA level is a missing value, not a category", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(file.path(tempfile(), "f.zarr"))
  n <- nrow(args$obs)
  args$obs$tagged <- addNA(factor(rep(c("x", NA, "y"), length.out = n), levels = c("x", "y")))
  args$latent_key <- "X_pca"
  args$centroid_min_points <- 1L
  store <- suppressMessages(do.call(cellucid_prepare, args))
  expect_identical(store_array(store, "obs/tagged/categories"), c("x", "y"))
  codes <- store_array(store, "obs/tagged/codes")
  expect_identical(as.integer(codes), as.integer(rep(c(0, -1, 1), length.out = n)))
  # The outlier quantiles give the missing cells none.
  quantiles <- store_array(store, "uns/cellucid/outlierQuantiles/fields/tagged")
  expect_true(all(is.nan(quantiles[codes == -1])))
  expect_false(any(is.nan(quantiles[codes != -1])))
})

test_that("a dense genes x cells assay is read through its transpose, never transposed whole", {
  matrix <- matrix(as.double(1:12), nrow = 3, dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
  view <- cellucid:::.transposed(matrix)
  expect_identical(cellucid:::.value_dim(view), c(4L, 3L))
  expect_identical(cellucid:::.matrix_block(view, 2:4, c(1, 3)), t(matrix)[2:4, c(1, 3), drop = FALSE])
  # The planned X of a feature-major dense matrix is that view.
  expression <- cellucid:::.expression_input(cellucid:::.feature_major(matrix), "gene_expression")
  captured <- NULL
  local_mocked_bindings(
    .ad_dense = function(store, path, matrix, gene_major = FALSE, dtype = "<f4") captured <<- matrix,
    .package = "cellucid"
  )
  cellucid:::.expression_put(NULL, "X", expression, "<f4")
  expect_s3_class(captured, "cellucid_transposed")
})

test_that("a matrix held outside memory is read in blocks of whole cells that cover it once", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("Matrix")
  set.seed(5)
  expected <- Matrix::rsparsematrix(40L, 23L, 0.2)
  expected[, 7] <- 0
  expected <- Matrix::drop0(expected)
  dimnames(expected) <- list(paste0("g", 1:40), paste0("c", 1:23))
  delayed <- DelayedArray::DelayedArray(expected)
  # Blocks of 80 values: two columns at a time, twelve blocks, never joined.
  local_mocked_bindings(.backed_block_values = 80, .package = "cellucid")
  blocks <- cellucid:::.backed_blocks(delayed)
  expect_identical(blocks$n, 12L)
  for (k in seq_len(blocks$n)) {
    read <- blocks$read(k)
    expect_false(read$by_rows)
    expect_identical(read$rows, 1:40)
    columns <- read$col0 + seq_len(length(read$pointers) - 1L)
    block <- methods::new("dgCMatrix", i = read$index, p = as.integer(read$pointers), x = read$values, Dim = c(40L, length(columns)))
    expect_equal(unname(as.matrix(block)), unname(as.matrix(expected[, columns, drop = FALSE])))
  }
  # Only the rows a pass needs.
  some <- cellucid:::.backed_blocks(delayed, c(3L, 9L, 40L))$read(1L)
  expect_identical(some$rows, c(3L, 9L, 40L))
  block <- methods::new("dgCMatrix", i = some$index, p = as.integer(some$pointers), x = some$values, Dim = c(3L, length(some$pointers) - 1L))
  expect_equal(unname(as.matrix(block)), unname(as.matrix(expected[c(3, 9, 40), some$col0 + seq_len(ncol(block)), drop = FALSE])))
})

test_that("genes kept from a genes x cells assay are read in place, in their order, never copied", {
  skip_if_not_installed("Matrix")
  set.seed(11)
  genes <- paste0("g", 1:30)
  cells <- paste0("c", 1:40)
  values <- Matrix::rsparsematrix(30L, 40L, 0.3, rand.x = function(n) round(stats::runif(n, 1, 50)))
  dimnames(values) <- list(genes, cells)
  # The row subset gene_identifiers used to make copied every value: the
  # kept genes are positions the gene-major arrays read.
  expression <- cellucid:::.expression_input(cellucid:::.feature_major(values), "gene_expression")
  kept <- cellucid:::.expression_select_genes(expression, c(3L, 7L, 8L))
  expect_identical(kept$matrix, expression$matrix)
  expect_identical(cellucid:::.expression_dims(kept), c(40L, 3L))
  expect_identical(cellucid:::.expression_gene_names(kept), c("g3", "g7", "g8"))
  for (chosen in list(c("g2", "g5", "g9", "g30"), c("g30", "g4", "g17", "g1", "g12"))) {
    components <- list(
      dataset_id = "kept", dataset_name = "Kept",
      embeddings = list(X_umap = matrix(seq_len(80) / 10, 40, 2, dimnames = list(cells, NULL))),
      obs = data.frame(row.names = cells),
      modalities = list(RNA = list(
        name = "RNA",
        gene_expression = cellucid:::.feature_major(values),
        var = data.frame(row.names = genes),
        cells = NULL
      )),
      gene_identifiers = chosen,
      notes = character(),
      kind = "Seurat"
    )
    store <- prepare_components(components)
    expect_identical(store_array(store, "var/_index"), chosen)
    expect_identical(store_x(store), unname(t(as.matrix(values[chosen, , drop = FALSE]))))
  }
})
