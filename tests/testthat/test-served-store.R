# The store cellucid_view() and cellucid_serve() serve from memory holds the
# arrays of the store cellucid_prepare(format = "zarr") writes, key for key
# and value for value: both are made by the one dataset plan (dataset-plan.R),
# the served one making each chunk only when it is asked for and leaving it
# uncompressed (it is read over loopback). Every array is read from both and
# compared -- dtype, shape and values -- by this package's test reader and by
# zarr-python (expect_same_arrays(), helper-objects.R).

store_path <- function(name = "data.zarr") file.path(tempfile("cellucid_served_"), name)

test_that("the served store of explicit inputs holds the prepared store's arrays", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  args$latent_key <- "X_pca"
  folder <- suppressMessages(do.call(cellucid_prepare, args))
  given <- args[setdiff(names(args), c("out_path", "format"))]
  store <- suppressMessages(memory_store_of(cellucid:::.dataset_components(NULL, given)))
  expect_same_arrays(folder, store)
})

test_that("a Seurat object's served store holds its prepared store's arrays", {
  skip_without_seurat()
  object <- pbmc_seurat()
  folder <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "pbmc", dataset_name = "PBMC"))
  store <- suppressMessages(memory_store_of(
    cellucid:::.dataset_components(object, list(dataset_id = "pbmc", dataset_name = "PBMC"))
  ))
  expect_same_arrays(folder, store)
})

test_that("a multimodal object's served MuData store holds its prepared store's arrays", {
  skip_without_seurat()
  object <- pbmc_multimodal()
  folder <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "pbmc", dataset_name = "PBMC"))
  store <- suppressMessages(memory_store_of(
    cellucid:::.dataset_components(object, list(dataset_id = "pbmc", dataset_name = "PBMC"))
  ))
  expect_same_arrays(folder, store)
  expect_true("mod/ADT/X/data/0" %in% store$keys)
})

test_that("a SingleCellExperiment's served store holds its prepared store's arrays", {
  skip_without_sce()
  object <- small_sce()
  folder <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "sce", dataset_name = "SCE"))
  store <- suppressMessages(memory_store_of(
    cellucid:::.dataset_components(object, list(dataset_id = "sce", dataset_name = "SCE"))
  ))
  expect_same_arrays(folder, store)
})

test_that("an object's genes x cells matrix is stored as the same cells x genes matrix", {
  # The object path never transposes the matrix (the indptr counted from the
  # row indices, the indices and data made block by block); the caller's
  # path takes cells x genes as it is. Both must write the same arrays.
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  from_object <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A"))
  logcounts <- SummarizedExperiment::assay(object, "logcounts")
  explicit <- suppressMessages(cellucid_prepare(
    store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A",
    gene_expression = Matrix::t(logcounts),
    var = as.data.frame(SummarizedExperiment::rowData(object))
  ))
  for (key in c("X/data", "X/indices", "X/indptr", "var/_index", "var/symbol")) {
    expect_identical(store_array(from_object, key), store_array(explicit, key), info = key)
  }
})

test_that("the served store holds what the object held when it was planned", {
  # R copies on change: an object changed after the viewer started is a new
  # object, and the served values -- even of the chunks first made after the
  # change -- stay those of the object as planned.
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  folder <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A"))
  store <- suppressMessages(memory_store_of(
    cellucid:::.dataset_components(object, list(dataset_id = "a", dataset_name = "A"))
  ))
  changed <- SummarizedExperiment::assay(object, "logcounts")
  changed@x[] <- 0
  SummarizedExperiment::assay(object, "logcounts") <- changed
  object$cluster <- NULL
  expect_same_arrays(folder, store)
})

test_that("a served store declares no compressor and no shuffle, and keeps its chunks' values as they are", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  store <- suppressMessages(memory_store_of(
    cellucid:::.dataset_components(object, list(dataset_id = "a", dataset_name = "A"))
  ))
  served <- dump_store(store)
  meta <- jsonlite::read_json(file.path(served, "obsm", "X_UMAP", ".zarray"))
  expect_null(meta$compressor)
  expect_null(meta$filters)
  umap <- SingleCellExperiment::reducedDim(object, "UMAP")
  chunk <- cellucid:::.store_get(store, "obsm/X_UMAP/0.0")
  expect_identical(length(chunk), as.integer(4L * prod(unlist(meta$chunks))))
  # Row-major float32, the cells' coordinates in order.
  values <- readBin(chunk, "double", n = length(umap), size = 4L, endian = "little")
  expect_equal(values, as.vector(t(umap)), tolerance = 1e-6)
})

# A genes x cells matrix of `n_cells` cells, each holding `per_cell` values
# at genes spread over `n_genes`, every value distinct.
spread_matrix <- function(n_genes, n_cells, per_cell) {
  rows <- unlist(lapply(seq_len(n_cells), function(cell) {
    sort(((cell * 7L + seq_len(per_cell) * 13L) %% n_genes))
  }))
  Matrix::sparseMatrix(
    i = rows, p = seq(0, by = per_cell, length.out = n_cells + 1L), x = seq_along(rows) / 7,
    dims = c(n_genes, n_cells), index1 = FALSE
  )
}

test_that("every gene-major block is the matching slice of the transpose, at any boundary", {
  skip_if_not_installed("Matrix")
  matrix <- spread_matrix(n_genes = 37L, n_cells = 400L, per_cell = 9L)
  # Every gene kept in its order, then some genes kept in another order (the
  # rows of gene_identifiers, read in place): each block is the slice of the
  # transpose of exactly those rows.
  for (kept in list(seq_len(nrow(matrix)), c(30L, 2L, 17L, 3L, 36L, 9L))) {
    map <- rep.int(-1L, nrow(matrix))
    map[kept] <- seq_along(kept) - 1L
    origin <- kept - 1L
    transposed <- Matrix::t(matrix[kept, , drop = FALSE])
    pointers <- cellucid:::.native_row_pointers(matrix@i, map, length(kept))
    expect_identical(pointers, as.numeric(transposed@p))
    n <- length(transposed@x)
    # Boundaries inside genes, on gene starts, the first and last value, empty.
    for (range in list(c(0, n), c(0, 1), c(n - 1, n), c(5, 5), c(17, min(n, 950)), c(pointers[[4]], pointers[[5]]), c(min(101, n), n))) {
      block <- cellucid:::.native_gene_major_block(matrix, map, origin, pointers, range[[1]], range[[2]], 4L, "<f8")
      positions <- seq_len(range[[2]] - range[[1]]) + range[[1]]
      expect_identical(readBin(block[[1]], "integer", n = length(positions), size = 4L, endian = "little"),
                       transposed@i[positions], info = paste(c(length(kept), range), collapse = "-"))
      expect_identical(readBin(block[[2]], "double", n = length(positions), size = 8L, endian = "little"),
                       transposed@x[positions], info = paste(c(length(kept), range), collapse = "-"))
    }
    block <- cellucid:::.native_gene_major_block(matrix, map, origin, pointers, 30, 90, 8L, "<f4")
    expect_identical(decode_chunk(block[[1]], "<i8", 60L), as.numeric(transposed@i[31:90]))
  }
})

test_that("a served store made of many small gene-major blocks, evicted from its cache, holds the prepared arrays", {
  # Over 2^16 values, so X's data and indices span several chunks; a cache
  # of a few kilobytes makes every block one chunk, and evicts it.
  skip_if_not_installed("Matrix")
  matrix <- spread_matrix(n_genes = 211L, n_cells = 9000L, per_cell = 17L)
  colnames(matrix) <- paste0("c", seq_len(ncol(matrix)))
  rownames(matrix) <- paste0("g", seq_len(nrow(matrix)))
  components <- list(
    embeddings = list(X_umap = matrix(seq_len(2L * ncol(matrix)) / 3, ncol = 2L, dimnames = list(colnames(matrix), NULL))),
    gene_expression = cellucid:::.feature_major(matrix),
    dataset_id = "blocks", dataset_name = "Blocks"
  )
  plan_components <- cellucid:::.dataset_components(NULL, components)
  folder <- file.path(tempfile("cellucid_blocks_"), "blocks.zarr")
  plan <- suppressMessages(cellucid:::.plan_dataset(plan_components))
  disk <- cellucid:::.store_open(folder, zip = FALSE, zarr_format = 2L)
  plan$write(disk)
  cellucid:::.store_close(disk)
  for (cache_bytes in c(0, 600000, 268435456)) {
    store <- suppressMessages(memory_store_of(plan_components, cache_bytes = cache_bytes))
    expect_same_arrays(folder, store)
    expect_lte(store$cache$used, cache_bytes)
  }
  transposed <- Matrix::t(matrix)
  expect_identical(store_array(folder, "X/indices"), transposed@i)
  expect_equal(store_array(folder, "X/data"), transposed@x, tolerance = 1e-7)
})
