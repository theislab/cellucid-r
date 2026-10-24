# The input format registry Cellucid reads by, shared with the web app and
# cellucid-python.
#
# tests/testthat/fixtures/input_formats.json is a byte-identical copy of the
# web app's tests/fixtures/input-formats.json (written from
# assets/js/data/input-formats.js) and of cellucid-python's packaged
# input_formats.json; all three repositories pin the same SHA-256. It names
# every format the viewer opens by its id, the name endings of each layout it
# arrives in, the container signatures and the folder markers. The id is what
# a dataset that states none is named after ("<format>_<name>") and what a
# server reports as the store's kind.

.input_formats_sha256 <- "926cbd0c4690e2b9e2616e088141ba06c5a54c35d41f94adbaf6bbf89b31208b"

test_that("the three repositories hold the same pinned input format registry", {
  expect_pinned_fixture("input_formats.json", .input_formats_sha256)
  registry <- cellucid_read_fixture("input_formats.json")
  expect_identical(registry$schema, 1L)
  ids <- vapply(registry$formats, function(entry) entry$id, character(1))
  expect_identical(ids, c("h5ad", "zarr"))
})

test_that("cellucid_prepare() writes the formats and layouts the registry names", {
  registry <- cellucid_read_fixture("input_formats.json")
  zarr <- Filter(function(entry) identical(entry$id, "zarr"), registry$formats)[[1L]]
  # format = "zarr" writes a .zarr folder, format = "zarr.zip" a .zarr.zip file.
  expect_true(".zarr" %in% unlist(zarr$names$folder))
  expect_true(".zarr.zip" %in% unlist(zarr$names$zip))
  markers <- unlist(registry$folderMarkers$zarr)
  expect_true(all(c(".zgroup", "zarr.json") %in% markers))
})
