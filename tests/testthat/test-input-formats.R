# The input format registry Cellucid reads by, shared with the web app and
# cellucid-python.
#
# inst/input_formats.json (installed as the package's input_formats.json,
# which cellucid_view() and cellucid_serve() classify paths with) is a
# byte-identical copy of the
# web app's tests/fixtures/input-formats.json (written from
# assets/js/data/input-formats.js) and of cellucid-python's packaged
# input_formats.json; all three repositories pin the same SHA-256. It names
# every format the viewer opens by its id, the name endings of each layout it
# arrives in, the container signatures and the folder markers. The id is what
# a dataset that states none is named after ("<format>_<name>") and what a
# server reports as the store's kind.

.input_formats_sha256 <- "ed44816fc5299aea7b87a8f7e805a704779517e04a2dd5b9b6340ced511ba7e7"

registry_path <- function() system.file("input_formats.json", package = "cellucid", mustWork = TRUE)

test_that("the three repositories hold the same pinned input format registry", {
  expect_identical(cellucid_file_sha256(registry_path()), .input_formats_sha256)
  registry <- jsonlite::fromJSON(registry_path(), simplifyVector = FALSE)
  expect_identical(registry$schema, 1L)
  ids <- vapply(registry$formats, function(entry) entry$id, character(1))
  expect_true(all(c("h5ad", "zarr", "loom", "tenx-h5", "tenx-outs", "h5mu", "mudata-zarr") %in% ids))
})

test_that("cellucid_prepare() writes the formats and layouts the registry names", {
  registry <- jsonlite::fromJSON(registry_path(), simplifyVector = FALSE)
  zarr <- Filter(function(entry) identical(entry$id, "zarr"), registry$formats)[[1L]]
  # format = "zarr" writes a .zarr folder, format = "zarr.zip" a .zarr.zip file.
  expect_true(".zarr" %in% unlist(zarr$names$folder))
  expect_true(".zarr.zip" %in% unlist(zarr$names$zip))
  markers <- unlist(registry$folderMarkers$zarr)
  expect_true(all(c(".zgroup", "zarr.json") %in% markers))
})

test_that("a path is classified as the registry and the web app classify it", {
  classify <- cellucid:::.classify_path
  folder <- tempfile("cellucid_paths_")
  dir.create(folder)
  hdf5 <- file.path(folder, "data.h5ad")
  writeBin(c(as.raw(c(0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a)), raw(100)), hdf5)
  expect_identical(classify(hdf5)[c("kind", "layout")], list(kind = "h5ad", layout = "file"))
  # Content decides the container, the name the format within it.
  loom <- file.path(folder, "data.loom")
  file.copy(hdf5, loom)
  expect_identical(classify(loom)$kind, "loom")
  # An HDF5 file named .zip is an HDF5 format, not a zipped store.
  zipped <- file.path(folder, "odd.zip")
  file.copy(hdf5, zipped)
  expect_identical(classify(zipped)$layout, "file")
  rds <- file.path(folder, "object.rds")
  saveRDS(list(1), rds)
  expect_identical(classify(rds)[c("kind", "layout")], list(kind = "rds", layout = "file"))
  text <- file.path(folder, "notes.txt")
  writeLines("hello", text)
  expect_error(classify(text), "A file Cellucid opens is one of .h5ad")
  expect_error(classify(file.path(folder, "missing.h5ad")), "Nothing exists at")
  expect_error(classify(folder), "is not a Zarr store")
  store <- suppressMessages(cellucid_prepare(
    priority_seed = NULL,
    file.path(folder, "s.zarr"), "zarr",
    embeddings = list(X_umap = diag(2)), dataset_id = "s", dataset_name = "S"
  ))
  expect_identical(classify(store)[c("kind", "layout", "zarr_format")], list(kind = "zarr", layout = "folder", zarr_format = 2L))
  zip <- suppressMessages(cellucid_prepare(
    priority_seed = NULL,
    file.path(folder, "s.zarr.zip"), "zarr.zip",
    embeddings = list(X_umap = diag(2)), dataset_id = "s", dataset_name = "S"
  ))
  expect_identical(classify(zip)[c("kind", "layout")], list(kind = "zarr", layout = "zip"))
})

test_that("an h5Seurat file is refused in the registry's words before it is served", {
  # The web app and cellucid serve (Python) refuse it with the same words:
  # the registry's refusal, which names the conversion route.
  registry <- jsonlite::fromJSON(registry_path(), simplifyVector = FALSE)
  h5seurat <- Filter(function(entry) identical(entry$id, "h5seurat"), registry$formats)[[1L]]
  expect_match(h5seurat$refusal, "SeuratDisk::LoadH5Seurat", fixed = TRUE)
  expect_match(h5seurat$refusal, "saveRDS", fixed = TRUE)
  folder <- tempfile("cellucid_h5seurat_")
  dir.create(folder)
  path <- file.path(folder, "pbmc.h5Seurat")
  writeBin(c(as.raw(c(0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a)), raw(100)), path)
  expect_error(cellucid:::.classify_path(path), paste0("pbmc.h5Seurat: ", h5seurat$refusal), fixed = TRUE)
})
