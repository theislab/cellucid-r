# obs['total_counts'] and uns['cellucid']['totalCounts']: each cell's total
# count (the Cellucid-ready AnnData specification, totalCounts; R/total-counts.R).
# cellucid-python's tests/test_prepare_total_counts_contract.py pins the same,
# and test-cross-writer.R checks both writers store the same bytes.

total_counts_matrix <- function() {
  set.seed(4)
  x <- Matrix::rsparsematrix(20L, 8L, 0.5, rand.x = function(n) as.double(rpois(n, 3) + 1))
  x[4L, ] <- 0
  dimnames(x) <- list(paste0("c", 1:20), paste0("g", 1:8))
  Matrix::drop0(x)
}

write_totals <- function(gene_expression, ..., obs = NULL, name = "a.zarr") {
  store <- file.path(tempfile("cellucid_totals_"), name)
  suppressMessages(cellucid_prepare(
    store, "zarr", embeddings = list(X_umap = matrix(seq_len(40) / 7, ncol = 2L)),
    gene_expression = gene_expression, obs = obs, dataset_id = "t", dataset_name = "T",
    priority_seed = NULL, ...
  ))
  store
}

totals_record <- function(store, key = "uns/cellucid/totalCounts") {
  if (!dir.exists(file.path(store, key))) return(NULL)
  list(field = store_array(store, file.path(key, "field")), over = store_array(store, file.path(key, "over")))
}

test_that("counts are summed over the stored genes, in every matrix form", {
  skip_if_not_installed("DelayedArray")
  x <- total_counts_matrix()
  sums <- as.numeric(Matrix::rowSums(x))
  forms <- list(
    sparse = x,
    dense = as.matrix(x),
    genes_by_cells = cellucid:::.feature_major(Matrix::t(x)),
    delayed = DelayedArray::DelayedArray(x)
  )
  for (form in names(forms)) {
    store <- write_totals(forms[[form]])
    expect_identical(totals_record(store), list(field = "total_counts", over = "stored genes"), info = form)
    expect_identical(store_json(store, "obs/total_counts/.zarray")$dtype, "<f8", info = form)
    expect_identical(store_array(store, "obs/total_counts"), sums, info = form)
  }
  expect_identical(store_array(write_totals(x), "obs/total_counts")[[4L]], 0)
})

test_that("the sums are of the genes stored after gene_identifiers", {
  skip_if_not_installed("DelayedArray")
  x <- total_counts_matrix()
  for (matrix in list(x, cellucid:::.feature_major(Matrix::t(x)), DelayedArray::DelayedArray(x))) {
    store <- write_totals(matrix, gene_identifiers = c("g1", "g5", "g8"))
    expect_identical(store_array(store, "obs/total_counts"), as.numeric(Matrix::rowSums(x[, c(1L, 5L, 8L)])))
  }
})

test_that("a matrix that holds no counts records no totals", {
  logged <- total_counts_matrix()
  logged@x <- log1p(logged@x)
  store <- write_totals(logged)
  expect_null(totals_record(store))
  expect_false(dir.exists(file.path(store, "obs/total_counts")))
})

test_that("the caller's column is recorded over all genes", {
  x <- total_counts_matrix()
  obs <- data.frame(nCount_RNA = as.numeric(Matrix::rowSums(x)) + 7, row.names = rownames(x))
  store <- write_totals(x, obs = obs, total_counts = "nCount_RNA")
  expect_identical(totals_record(store), list(field = "nCount_RNA", over = "all genes"))
  expect_false(dir.exists(file.path(store, "obs/total_counts")))
})

test_that("a column that is not the totals is refused", {
  x <- total_counts_matrix()
  sums <- as.numeric(Matrix::rowSums(x))
  expect_error(write_totals(x, obs = data.frame(tc = log1p(sums), row.names = rownames(x)), total_counts = "tc"),
               "below the sum of the genes stored for")
  expect_error(write_totals(x, obs = data.frame(tc = replace(sums, 5L, NA), row.names = rownames(x)), total_counts = "tc"),
               "with no total of 0 or more")
  expect_error(write_totals(x, obs = data.frame(tc = rep("a", 20L), row.names = rownames(x)), total_counts = "tc"),
               "not a column of numbers")
  expect_error(write_totals(x, total_counts = "tc"), "which obs does not hold")
  logged <- x
  logged@x <- log1p(logged@x)
  expect_error(write_totals(logged, obs = data.frame(tc = sums + 1, row.names = rownames(x)), total_counts = "tc"),
               "does not hold counts")
  expect_error(write_totals(x, total_counts = c("a", "b")), "names the obs column")
})

test_that("the data's own total_counts column is kept", {
  x <- total_counts_matrix()
  sums <- as.numeric(Matrix::rowSums(x))
  same <- write_totals(x, obs = data.frame(total_counts = sums, row.names = rownames(x)))
  expect_identical(totals_record(same), list(field = "total_counts", over = "stored genes"))
  other <- write_totals(x, obs = data.frame(total_counts = sums + 3, row.names = rownames(x)))
  expect_null(totals_record(other))
  expect_identical(store_array(other, "obs/total_counts"), sums + 3)
})

test_that("a multimodal object records each counts modality's totals", {
  skip_without_sce()
  store <- file.path(tempfile("cellucid_totals_"), "m.zarr")
  sce <- small_sce()
  suppressMessages(cellucid_prepare(store, "zarr", data = sce, dataset_id = "m", dataset_name = "M"))
  expect_identical(totals_record(store, "uns/cellucid/totalCounts/ADT"), list(field = "ADT:total_counts", over = "stored genes"))
  expect_false(dir.exists(file.path(store, "uns/cellucid/totalCounts/main")))
  adt <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, "ADT"))
  rows <- store_array(store, "mod/ADT/obs/cellucid_original_position") + 1
  expect_identical(store_array(store, "mod/ADT/obs/total_counts"), as.numeric(Matrix::colSums(adt))[rows])
  named <- file.path(tempfile("cellucid_totals_"), "n.zarr")
  SummarizedExperiment::colData(sce)$adt_total <- as.numeric(Matrix::colSums(adt)) + 1
  suppressMessages(cellucid_prepare(named, "zarr", data = sce, dataset_id = "m", dataset_name = "M",
                                    total_counts = c(ADT = "adt_total")))
  expect_identical(totals_record(named, "uns/cellucid/totalCounts/ADT"), list(field = "adt_total", over = "all genes"))
  expect_error(suppressMessages(cellucid_prepare(file.path(tempfile("cellucid_totals_"), "x.zarr"), "zarr", data = sce,
                                                 dataset_id = "m", dataset_name = "M", total_counts = c(main = "adt_total"))),
               "Modality 'main'.*does not hold counts")
  expect_error(suppressMessages(cellucid_prepare(file.path(tempfile("cellucid_totals_"), "y.zarr"), "zarr", data = sce,
                                                 dataset_id = "m", dataset_name = "M", total_counts = c(atac = "x"))),
               "not modalities written")
  expect_error(suppressMessages(cellucid_prepare(file.path(tempfile("cellucid_totals_"), "z.zarr"), "zarr", data = sce,
                                                 dataset_id = "m", dataset_name = "M", total_counts = "adt_total")),
               "With several modalities")
})

test_that("a column over the source's gene panel is recorded as a source panel (F154)", {
  # The built-in atlases' totals cover their 8,192-gene integration panels:
  # more genes than are stored, not the whole transcriptome.
  x <- total_counts_matrix()
  sums <- as.numeric(Matrix::rowSums(x))
  obs <- data.frame(tc = sums + 5, row.names = rownames(x))
  store <- write_totals(x, obs = obs, total_counts = "tc", total_counts_over = "source panel")
  expect_identical(totals_record(store), list(field = "tc", over = "source panel"))
  # Checked as "all genes" is.
  expect_error(write_totals(x, obs = data.frame(tc = sums / 2, row.names = rownames(x)), total_counts = "tc",
                            total_counts_over = "source panel"), "below the sum of the genes stored")
  expect_error(write_totals(x, total_counts_over = "source panel"), "no total_counts was given")
  expect_error(write_totals(x, obs = obs, total_counts = "tc", total_counts_over = "stored genes"),
               "total_counts_over must be \"all genes\".*or \"source panel\"")
})

test_that("an object's record of a source panel's totals is kept, and a column named anew says its own", {
  skip_without_sce()
  sce <- small_sce(alt = FALSE)
  counts <- SummarizedExperiment::assay(sce, "counts")
  SummarizedExperiment::colData(sce)$tc <- as.numeric(Matrix::colSums(counts)) + 3
  S4Vectors::metadata(sce)$cellucid <- list(id = "p", name = "P", totalCounts = list(field = "tc", over = "source panel"))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts),
    colData = SummarizedExperiment::colData(sce), reducedDims = SingleCellExperiment::reducedDims(sce),
    metadata = S4Vectors::metadata(sce))
  kept <- file.path(tempfile("cellucid_totals_"), "k.zarr")
  suppressMessages(cellucid_prepare(kept, "zarr", data = sce))
  expect_identical(totals_record(kept), list(field = "tc", over = "source panel"))
  renamed <- file.path(tempfile("cellucid_totals_"), "r.zarr")
  suppressMessages(cellucid_prepare(renamed, "zarr", data = sce, total_counts = "tc"))
  expect_identical(totals_record(renamed), list(field = "tc", over = "all genes"))
})

test_that("several modalities say what each named column covers", {
  skip_without_sce()
  sce <- small_sce()
  adt <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, "ADT"))
  SummarizedExperiment::colData(sce)$adt_total <- as.numeric(Matrix::colSums(adt)) + 1
  store <- file.path(tempfile("cellucid_totals_"), "m.zarr")
  suppressMessages(cellucid_prepare(store, "zarr", data = sce, dataset_id = "m", dataset_name = "M",
                                    total_counts = c(ADT = "adt_total"), total_counts_over = c(ADT = "source panel")))
  expect_identical(totals_record(store, "uns/cellucid/totalCounts/ADT"), list(field = "adt_total", over = "source panel"))
  expect_error(suppressMessages(cellucid_prepare(file.path(tempfile("cellucid_totals_"), "x.zarr"), "zarr", data = sce,
                                                 dataset_id = "m", dataset_name = "M", total_counts = c(ADT = "adt_total"),
                                                 total_counts_over = c(main = "source panel"))),
               "total_counts_over names c\\(\"main\"\\)")
})
