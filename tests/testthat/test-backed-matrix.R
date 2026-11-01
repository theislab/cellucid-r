# Expression matrices held outside R's memory (BPCells, DelayedArray
# backends; R/backed-matrix.R) are read block by block -- counted, then
# placed gene by gene -- and give the bytes a dgCMatrix of the same values
# gives, written or served, in priority order or not, whatever the blocks and
# the passes.

store_path <- function(name = "data.zarr") file.path(tempfile("cellucid_backed_"), name)

backed_values <- function(n_genes = 50L, n_cells = 400L, seed = 3) {
  set.seed(seed)
  values <- Matrix::rsparsematrix(n_genes, n_cells, 0.2, rand.x = function(n) round(stats::runif(n, 1, 30)))
  values[, 9] <- 0
  values <- Matrix::drop0(values)
  dimnames(values) <- list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells)))
  values
}

backed_sce <- function(matrix, umap) {
  object <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = matrix))
  SingleCellExperiment::reducedDim(object, "UMAP") <- umap
  object
}

# A BPCells function, looked up when called: BPCells is not on CRAN, so the
# package declares it nowhere, and these tests skip without it.
bpcells <- function(name) getExportedValue("BPCells", name)

# A BPCells matrix directory of `values`, in its storage order ("col": cells
# together, as a 10x or Seurat matrix is; "row": genes together).
bpcells_of <- function(values, order = "col") {
  written <- bpcells("write_matrix_dir")(
    bpcells("convert_matrix_type")(methods::as(values, "IterableMatrix"), "uint32_t"), tempfile("cellucid_bpcells_")
  )
  if (identical(order, "row")) bpcells("transpose_storage_order")(written, tempfile("cellucid_bpcells_")) else written
}

# Every file of two stores (a folder, or one .zarr.zip) identical, byte for byte.
expect_same_store <- function(actual, expected) {
  if (!dir.exists(expected)) {
    expect_identical(readBin(actual, "raw", file.size(actual)), readBin(expected, "raw", file.size(expected)))
    return(invisible(NULL))
  }
  files <- sort(list.files(expected, recursive = TRUE, all.files = TRUE))
  expect_identical(sort(list.files(actual, recursive = TRUE, all.files = TRUE)), files)
  differ <- files[!vapply(files, function(file) {
    identical(
      readBin(file.path(actual, file), "raw", file.size(file.path(actual, file))),
      readBin(file.path(expected, file), "raw", file.size(file.path(expected, file)))
    )
  }, logical(1))]
  expect_identical(differ, character(0))
}

backed_kinds <- function(values) {
  kinds <- list()
  if (requireNamespace("BPCells", quietly = TRUE)) {
    kinds$`BPCells, cells together` <- bpcells_of(values, "col")
    kinds$`BPCells, genes together` <- bpcells_of(values, "row")
  }
  if (requireNamespace("DelayedArray", quietly = TRUE)) {
    kinds$`DelayedArray, sparse` <- DelayedArray::DelayedArray(values)
    kinds$`DelayedArray, dense` <- DelayedArray::DelayedArray(as.matrix(values))
  }
  kinds
}

test_that("a backed assay is written byte for byte as its dgCMatrix, in priority order or not", {
  skip_without_sce()
  values <- backed_values()
  umap <- matrix(stats::rnorm(2L * ncol(values)), ncol = 2L, dimnames = list(colnames(values), NULL))
  kinds <- backed_kinds(values)
  if (length(kinds) == 0L) skip("Neither BPCells nor DelayedArray is installed.")
  prepare <- function(matrix, format = "zarr", ...) {
    path <- store_path(if (identical(format, "zarr")) "data.zarr" else "data.zarr.zip")
    suppressMessages(cellucid_prepare(path, format, data = backed_sce(matrix, umap), dataset_id = "d", dataset_name = "D", ...))
  }
  for (arguments in list(list(), list(priority_seed = NULL), list(gene_identifiers = c("g5", "g2", "g40", "g33")))) {
    expected <- do.call(prepare, c(list(values), arguments))
    for (kind in names(kinds)) {
      expect_same_store(do.call(prepare, c(list(kinds[[kind]]), arguments)), expected)
    }
    # Blocks of a few cells or genes, and passes of a few hundred values: X
    # placed in many passes, each block of positions let go after its chunks.
    local_mocked_bindings(.backed_block_values = 120, .backed_group_values = function(bytes_per_value) 300)
    for (kind in names(kinds)) {
      expect_same_store(do.call(prepare, c(list(kinds[[kind]]), arguments)), expected)
    }
  }
  expect_same_store(prepare(kinds[[1L]], "zarr.zip"), prepare(values, "zarr.zip"))
})

test_that("a backed assay is served byte for byte as its dgCMatrix, laid out in passes of a share of memory", {
  skip_without_sce()
  values <- backed_values(n_genes = 40L, n_cells = 600L, seed = 4)
  umap <- matrix(stats::rnorm(2L * ncol(values)), ncol = 2L, dimnames = list(colnames(values), NULL))
  kinds <- backed_kinds(values)
  if (length(kinds) == 0L) skip("Neither BPCells nor DelayedArray is installed.")
  served <- function(matrix) {
    components <- cellucid:::.dataset_components(backed_sce(matrix, umap), list(dataset_id = "d", dataset_name = "D"))
    components$priority_seed <- cellucid:::.served_priority_seed
    store <- suppressMessages(memory_store_of(components))
    list(folder = dump_store(store), store = store)
  }
  expected <- served(values)$folder
  for (kind in names(kinds)) {
    expect_same_store(served(kinds[[kind]])$folder, expected)
  }
  local_mocked_bindings(.backed_block_values = 200, .backed_group_values = function(bytes_per_value) 500)
  for (kind in names(kinds)) {
    run <- served(kinds[[kind]])
    expect_same_store(run$folder, expected)
    # One laid-out matrix, in scratch files, let go when the viewer stops.
    expect_length(run$store$backed_layouts, 1L)
    cellucid:::.backed_release_layouts(run$store)
    expect_null(run$store$backed_layouts)
  }
})

test_that("a backed matrix is never read whole: one pass to count, then one or two a share to place, each over its genes", {
  skip_without_sce()
  skip_if_not_installed("BPCells")
  set.seed(5)
  values <- Matrix::rsparsematrix(100L, 1000L, 0.4, rand.x = function(n) round(stats::runif(n, 1, 30)))
  dimnames(values) <- list(paste0("g", 1:100), paste0("c", 1:1000))
  umap <- matrix(stats::rnorm(2L * ncol(values)), ncol = 2L, dimnames = list(colnames(values), NULL))
  matrix <- bpcells_of(values)
  read <- new.env()
  blocks <- cellucid:::.backed_blocks
  local_mocked_bindings(
    .backed_block_values = 2000,
    .backed_blocks = function(matrix, rows = NULL) {
      plan <- blocks(matrix, rows)
      list(n = plan$n, read = function(k) {
        block <- plan$read(k)
        read$blocks[[length(read$blocks) + 1L]] <- c(length(block$rows), length(block$pointers) - 1)
        block
      })
    }
  )
  passes <- function(...) {
    read$blocks <- list()
    suppressMessages(cellucid_prepare(store_path(), "zarr", data = backed_sce(matrix, umap), dataset_id = "d", dataset_name = "D", ...))
    # Whole cells, at most 2000 values of the block's rows at a time.
    expect_true(all(vapply(read$blocks, function(dims) dims[[1L]] * dims[[2L]] <= 2000, logical(1))))
    sum(vapply(read$blocks, function(dims) dims[[2L]], numeric(1))) / ncol(values)
  }
  # Count, then place; the prefix counts are counted as X is placed.
  expect_identical(passes(priority_seed = NULL), 2)
  expect_identical(passes(), 2)
  # Passes of one chunk of values (2^14): X placed once a block of positions
  # for its values, then again for its cells, each pass over its genes' rows.
  local_mocked_bindings(.backed_group_values = function(bytes_per_value) 16384)
  n_blocks <- ceiling(length(values@x) / 16384)
  expect_gt(n_blocks, 1)
  expect_identical(passes(), 1 + 2 * n_blocks)
  expect_true(all(vapply(read$blocks, function(dims) dims[[1L]], numeric(1))[-seq_len(50)] < 100))
})

test_that("a backed assay keeps integer counts above 2^24 exact, decided on the genes stored", {
  skip_without_sce()
  skip_if_not_installed("BPCells")
  values <- backed_values(n_genes = 20L, n_cells = 100L, seed = 6)
  values[3, 5] <- 2^25 + 1
  umap <- matrix(stats::rnorm(2L * ncol(values)), ncol = 2L, dimnames = list(colnames(values), NULL))
  matrix <- bpcells_of(values)
  dtype_of <- function(store) store_json(store, "X/data/.zarray")$dtype
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = backed_sce(matrix, umap), dataset_id = "d", dataset_name = "D", priority_seed = NULL))
  expect_identical(dtype_of(store), "<f8")
  expect_true((2^25 + 1) %in% store_array(store, "X/data"))
  kept <- suppressMessages(cellucid_prepare(
    store_path(), "zarr", data = backed_sce(matrix, umap), dataset_id = "d", dataset_name = "D",
    priority_seed = NULL, gene_identifiers = c("g1", "g2", "g4")
  ))
  expect_identical(dtype_of(kept), "<f4")
})

test_that("preparing a DelayedArray assay warns nothing; only the OpenMP pause warning is muffled", {
  skip_without_sce()
  skip_if_not_installed("DelayedArray")
  values <- backed_values()
  umap <- matrix(stats::rnorm(2L * ncol(values)), ncol = 2L, dimnames = list(colnames(values), NULL))
  # Blocks of a few cells: many conversions, each of which warned once where
  # the OpenMP runtime cannot release its threads (macOS's libomp).
  local_mocked_bindings(.backed_block_values = 400)
  expect_no_warning(suppressMessages(cellucid_prepare(
    store_path(), "zarr", data = backed_sce(DelayedArray::DelayedArray(values), umap), dataset_id = "d", dataset_name = "D"
  )))
  expect_no_warning(.without_omp_pause_warning(warning("omp_pause_resource() failed to relinquish resources on device 0")))
  expect_warning(.without_omp_pause_warning(warning("another warning")), "another warning")
})

test_that("a backed assay holding a value that cannot be stored is refused when it is counted", {
  skip_without_sce()
  skip_if_not_installed("DelayedArray")
  values <- backed_values(n_genes = 10L, n_cells = 40L, seed = 7)
  umap <- matrix(stats::rnorm(2L * ncol(values)), ncol = 2L, dimnames = list(colnames(values), NULL))
  broken <- values
  broken[2, 3] <- NA
  broken[4, 5] <- Inf
  broken[6, 7] <- NaN
  # Every NA and NaN is counted, on a sparse backend and on a dense one
  # (DelayedArray's own coercion to a dgCMatrix dropped both before
  # Bioconductor 3.19: R/backed-matrix.R, .delayed_block_dgcmatrix()), and
  # in blocks of a few cells too.
  whole <- cellucid:::.backed_block_values
  for (seed in list(broken, as.matrix(broken))) {
    for (block_values in c(whole, 30)) {
      local_mocked_bindings(.backed_block_values = block_values)
      expect_error(
        suppressMessages(cellucid_prepare(store_path(), "zarr", data = backed_sce(DelayedArray::DelayedArray(seed), umap), dataset_id = "d", dataset_name = "D")),
        "must contain only finite values.*1 NaN, 1 NA, 1 \\+Inf"
      )
    }
  }
  broken <- values
  broken[2, 3] <- 1e39
  expect_error(
    suppressMessages(cellucid_prepare(store_path(), "zarr", data = backed_sce(DelayedArray::DelayedArray(broken), umap), dataset_id = "d", dataset_name = "D")),
    "beyond the float32 range"
  )
})

test_that("a cells x genes backed matrix passed as gene_expression is read through its transpose", {
  skip_if_not_installed("BPCells")
  values <- backed_values(n_genes = 25L, n_cells = 200L, seed = 8)
  umap <- matrix(stats::rnorm(2L * ncol(values)), ncol = 2L)
  prepare <- function(matrix) {
    suppressMessages(cellucid_prepare(
      store_path(), "zarr", embeddings = list(X_umap = umap), gene_expression = matrix,
      obs = data.frame(row.names = colnames(values)), var = data.frame(row.names = rownames(values)),
      dataset_id = "d", dataset_name = "D"
    ))
  }
  expect_same_store(prepare(Matrix::t(bpcells_of(values))), prepare(Matrix::t(values)))
})

test_that("a Seurat v5 assay whose layer is a BPCells matrix is written as its dgCMatrix", {
  skip_without_seurat()
  skip_if_not_installed("BPCells")
  skip_if_not(seurat_v5())
  object <- pbmc_seurat()
  data <- methods::as(SeuratObject::LayerData(object, assay = "RNA", layer = "data"), "dgCMatrix")
  expected <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "p", dataset_name = "P"))
  backed <- suppressWarnings(suppressMessages(bpcells("write_matrix_dir")(methods::as(data, "IterableMatrix"), tempfile("cellucid_bpcells_"))))
  SeuratObject::LayerData(object, assay = "RNA", layer = "data") <- backed
  expect_true(inherits(SeuratObject::LayerData(object, assay = "RNA", layer = "data"), "IterableMatrix"))
  expect_message(
    actual <- cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "p", dataset_name = "P"),
    "block by block, never whole"
  )
  expect_same_store(actual, expected)
})

test_that("a gene given more or fewer values than it was counted is refused, not written", {
  skip_if_not_installed("BPCells")
  values <- backed_values(n_genes = 12L, n_cells = 60L, seed = 9)
  expression <- suppressMessages(cellucid:::.expression_input(cellucid:::.feature_major(bpcells_of(values)), "gene_expression"))
  # The matrix changes after it was counted: a value more in gene 2.
  changed <- values
  changed[2, which(values[2, ] == 0)[[1L]]] <- 7
  expression$matrix <- bpcells_of(changed)
  store <- cellucid:::.store_open(tempfile(), zip = FALSE, zarr_format = 2L)
  arrays <- cellucid:::.backed_arrays(store, "X", expression, c(60L, 12L), "<f4")
  expect_error(arrays$values$chunk(1, 10, 16384), "holds more values than were counted")
})
