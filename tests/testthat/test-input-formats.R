# The input format registry Cellucid reads by, shared with the web app and
# cellucid-python.
#
# inst/input_formats.json (installed as the package's input_formats.json,
# which cellucid_view() and cellucid_serve() classify paths with) is a
# byte-identical copy of the
# web app's tests/fixtures/input-formats.json (its home, written from
# assets/js/data/input-formats.js) and of cellucid-python's packaged
# input_formats.json; each package's test compares its copy with the web
# repository's file. It names
# every format the viewer opens by its id, the name endings of each layout it
# arrives in, the container signatures and the folder markers. The id is what
# a dataset that states none is named after ("<format>_<name>_<structure>")
# and what a server reports as the store's kind.

registry_path <- function() system.file("input_formats.json", package = "cellucid", mustWork = TRUE)

test_that("the input format registry is the web app's, byte for byte", {
  expect_shared_file("inst/input_formats.json", "tests/fixtures/input-formats.json")
})

test_that("the packaged registry names the formats the viewer opens", {
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
  # One sentence in the three implementations (the web app's openableNamesText,
  # cellucid-python's openable_files_text()): the ending that stands alone,
  # h5Seurat (refused) left out.
  expected <- "an .h5ad file, a .loom file, an .h5 file, an .h5mu file, an .rds file, an .rda/.rdata file or a .zarr.zip archive"
  expect_identical(cellucid:::.openable_files_text(), expected)
  expect_error(classify(text), paste0(" is not a file Cellucid opens (", expected, ")."), fixed = TRUE)
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

test_that("a gzip-compressed copy is refused with how to serve it, as cellucid serve refuses it", {
  # GEO ships many .h5ad files gzipped: the web app reads such a copy whole,
  # but a served file is read in parts, which a gzip stream has none of.
  folder <- tempfile("cellucid_gz_")
  dir.create(folder)
  copy <- file.path(folder, "GSE1_cite.h5ad.gz")
  connection <- gzfile(copy, "wb")
  writeBin(c(as.raw(c(0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a)), raw(64)), connection)
  close(connection)
  expect_error(
    cellucid:::.classify_path(copy),
    "GSE1_cite.h5ad.gz is a gzip-compressed copy of AnnData (.h5ad): it cannot be read in parts. Decompress it (gunzip GSE1_cite.h5ad.gz)",
    fixed = TRUE
  )
  # A format read as one stream names its compressed endings itself.
  rds <- file.path(folder, "object.rds.gz")
  saveRDS(list(1), rds)
  expect_identical(cellucid:::.classify_path(rds)$kind, "rds")
})

test_that("Cellucid's own files are said for what they are, in the registry's words", {
  # The web app and cellucid serve (Python) say the same sentence (cellucidFiles' said).
  registry <- jsonlite::fromJSON(registry_path(), simplifyVector = FALSE)
  expect_identical(vapply(registry$cellucidFiles, function(entry) entry$id, ""), c("library", "session"))
  folder <- tempfile("cellucid_own_")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  for (name in c("annotation.cellucid-session", "LAB.Cellucid-Library")) {
    entry <- cellucid:::.cellucid_file_of_name(name)
    expect_identical(entry$id, if (endsWith(name, "session")) "session" else "library")
    path <- file.path(folder, name)
    writeBin(charToRaw("CELLUCID_SESSION/4\n"), path)
    expect_error(cellucid:::.classify_path(path), paste0(name, ": ", entry$said), fixed = TRUE)
    expect_error(cellucid_serve(path), entry$said, fixed = TRUE)
  }
  expect_null(cellucid:::.cellucid_file_of_name("cells.h5ad"))
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

test_that("a vlen-utf8 chunk's counts and lengths are checked against its bytes before anything is allocated", {
  # Unchecked, four bytes stating 2^31 - 1 strings would allocate a character
  # vector that long before a string is read.
  le <- function(value) writeBin(as.integer(value), raw(), size = 4L, endian = "little")
  good <- c(le(2L), le(2L), charToRaw("ab"), le(0L))
  expect_identical(cellucid:::.decode_vlen_utf8(good), c("ab", ""))
  expect_error(cellucid:::.decode_vlen_utf8(le(.Machine$integer.max)), "more strings than its bytes hold")
  expect_error(cellucid:::.decode_vlen_utf8(c(le(-1L), le(0L))), "more strings than its bytes hold")
  expect_error(cellucid:::.decode_vlen_utf8(raw(2)), "shorter than its count")
  expect_error(cellucid:::.decode_vlen_utf8(c(le(1L), le(1000L), charToRaw("ab"))), "runs past the end")
  expect_error(cellucid:::.decode_vlen_utf8(c(le(1L), le(-5L))), "runs past the end")
  expect_error(cellucid:::.decode_vlen_utf8(c(le(2L), le(3L), charToRaw("abc"), raw(2))), "ends inside a string.s length")
})
