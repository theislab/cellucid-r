# Small single-cell objects of every class the object front door reads, made
# with the installed SeuratObject (v4 or v5: the suite runs against both) and
# SingleCellExperiment. Each helper skips, naming what is missing.

skip_without_seurat <- function() {
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("Matrix")
}

skip_without_sce <- function() {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("Matrix")
}

seurat_v5 <- function() utils::packageVersion("SeuratObject") >= "5.0.0"

# SeuratObject's pbmc_small, brought to the installed version's layout (a v5
# Assay5 with SeuratObject 5, a v4 Assay with SeuratObject 4).
pbmc_seurat <- function() {
  env <- new.env()
  utils::data("pbmc_small", package = "SeuratObject", envir = env)
  object <- suppressWarnings(suppressMessages(SeuratObject::UpdateSeuratObject(env$pbmc_small)))
  if (seurat_v5() && !methods::is(object[["RNA"]], "Assay5")) {
    suppressWarnings(object[["RNA"]] <- methods::as(object[["RNA"]], "Assay5"))
  }
  object
}

# pbmc_small with a second assay (ADT-like: 5 features over every cell).
pbmc_multimodal <- function() {
  object <- pbmc_seurat()
  counts <- Matrix::Matrix(
    matrix(rpois(5L * ncol(object), 3), nrow = 5L, dimnames = list(paste0("prot", 1:5), colnames(object))),
    sparse = TRUE
  )
  create <- if (seurat_v5()) SeuratObject::CreateAssay5Object else SeuratObject::CreateAssayObject
  object[["ADT"]] <- create(counts = counts)
  object
}

small_sce <- function(alt = TRUE) {
  set.seed(2)
  n_cells <- 12L
  counts <- Matrix::rsparsematrix(30L, n_cells, 0.3, rand.x = function(n) as.double(rpois(n, 4) + 1))
  dimnames(counts) <- list(paste0("gene", 1:30), paste0("cell", 1:n_cells))
  logcounts <- counts
  logcounts@x <- log1p(logcounts@x)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts, logcounts = logcounts),
    colData = S4Vectors::DataFrame(
      cluster = factor(rep(c("a", "b", "c"), 4L)),
      depth = seq_len(n_cells) * 10,
      label = paste0("L", seq_len(n_cells)),
      row.names = colnames(counts)
    ),
    rowData = S4Vectors::DataFrame(symbol = paste0("S", 1:30), row.names = rownames(counts)),
    reducedDims = list(
      UMAP = matrix(rnorm(n_cells * 2L), ncol = 2L),
      PCA = matrix(rnorm(n_cells * 4L), ncol = 4L)
    )
  )
  SingleCellExperiment::sizeFactors(sce) <- seq(0.5, 1.6, length.out = n_cells)
  if (alt) {
    adt <- Matrix::Matrix(matrix(rpois(4L * n_cells, 5), 4L, dimnames = list(paste0("ab", 1:4), colnames(counts))), sparse = TRUE)
    SingleCellExperiment::altExp(sce, "ADT") <- SummarizedExperiment::SummarizedExperiment(list(counts = adt))
  }
  graph <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4), j = c(2, 1, 4, 3), x = c(0.5, 0.5, 1, 1), dims = c(n_cells, n_cells)
  )
  SingleCellExperiment::colPair(sce, "snn") <- graph
  sce
}

# A memory store written key by key to a folder (each chunk made as the
# server makes it), so it is read as any store on disk.
dump_store <- function(store, folder = tempfile("cellucid_dumped_")) {
  for (key in store$keys) {
    target <- file.path(folder, key)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    writeBin(cellucid:::.store_get(store, key), target)
  }
  folder
}

# The metadata of every node of a Zarr format 2 store (its .zmetadata
# members), with the codecs a served store leaves out removed: each array's
# compressor, and its filters (the shuffle filter; text keeps vlen-utf8).
metadata_without_codecs <- function(folder) {
  members <- jsonlite::read_json(file.path(folder, ".zmetadata"), simplifyVector = FALSE)$metadata
  lapply(members, function(member) {
    if (!is.null(member$zarr_format) && !is.null(member$dtype)) {
      member$compressor <- NULL
      member$filters <- NULL
    }
    member
  })
}

# The store cellucid_view() serves (a memory store) holds the arrays of the
# store cellucid_prepare() wrote: the same keys, every metadata document the
# same but for the codecs, which the served one leaves out (it is read over
# loopback, uncompressed), and every array the same dtype, shape and values
# -- read by this package's own test reader, and by zarr-python when one is
# found (compare_stores.py). `served` is the memory store, or a folder of the
# bytes its server answered with.
expect_same_arrays <- function(folder, served) {
  if (is.environment(served)) served <- dump_store(served)
  files <- list.files(folder, recursive = TRUE, all.files = TRUE)
  testthat::expect_setequal(list.files(served, recursive = TRUE, all.files = TRUE), files)
  testthat::expect_identical(metadata_without_codecs(served), metadata_without_codecs(folder))
  served_members <- jsonlite::read_json(file.path(served, ".zmetadata"), simplifyVector = FALSE)$metadata
  arrays <- sub("/\\.zarray$", "", grep("(^|/)\\.zarray$", names(served_members), value = TRUE))
  for (key in arrays) {
    meta <- served_members[[paste0(key, "/.zarray")]]
    testthat::expect_null(meta$compressor, info = key)
    text <- identical(meta$dtype, "|O")
    testthat::expect_identical(meta$filters, if (text) list(list(id = "vlen-utf8")), info = key)
    testthat::expect_identical(store_array(served, key), store_array(folder, key), info = key)
  }
  python <- zarr_python()
  if (!is.null(python)) {
    output <- suppressWarnings(system2(
      python, c(shQuote(testthat::test_path("python", "compare_stores.py")), shQuote(folder), shQuote(served)),
      stdout = TRUE, stderr = TRUE
    ))
    testthat::expect_identical(attr(output, "status"), NULL, info = paste(output, collapse = "\n"))
    testthat::expect_identical(output, paste("ok", length(arrays)))
  }
  invisible(served)
}

# A Python with zarr-python 3 (CELLUCID_ANNDATA_PYTHON, CELLUCID_PYTHON or
# python3), or NULL; never on CRAN.
zarr_python <- function() {
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) return(NULL)
  candidates <- c(Sys.getenv("CELLUCID_ANNDATA_PYTHON"), Sys.getenv("CELLUCID_PYTHON"), Sys.which("python3"))
  for (candidate in candidates[nzchar(candidates)]) {
    status <- suppressWarnings(system2(
      candidate, c("-c", shQuote("import zarr, sys; sys.exit(int(zarr.__version__.split('.')[0]) < 3)")),
      stdout = FALSE, stderr = FALSE
    ))
    if (identical(status, 0L)) return(candidate)
  }
  NULL
}

memory_store_of <- function(components, cache_bytes = 268435456) {
  plan <- cellucid:::.plan_dataset(components)
  store <- cellucid:::.store_open_memory(2L, cache_bytes)
  plan$write(store)
  cellucid:::.store_close(store)
  store
}
