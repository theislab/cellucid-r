# The gene sample, as R reads it: a large dataset stored cell by cell opens
# every cell on the map in the viewer and reads a gene for a sample of the
# open cells (the first in their order). A session's `genes` record and the
# `color_by` event's `cells_read` are checked as the web app
# (decodeGeneSampleRecord, notifyColorBy) and cellucid-python
# (tests/test_gene_sample_contract.py) check them.

genes_record <- list(cells = 40, order = list(method = "blocks", seed = 0, blockSize = 19), digest = "0123456789abcdef")

test_that("a session's gene sample is read as recorded", {
  read <- cellucid:::.session_genes(genes_record, 100)
  expect_identical(read$cells, 40)
  expect_identical(read$order$blockSize, 19)
  expect_identical(read$digest, "0123456789abcdef")
  expect_null(cellucid:::.session_genes(utils::modifyList(genes_record, list(order = NULL), keep.null = TRUE), 100)$order)
})

test_that("a session's gene sample out of its contract is refused", {
  refused <- function(genes, message) {
    expect_error(cellucid:::.session_genes(genes, 100), message, fixed = TRUE)
  }
  refused(utils::modifyList(genes_record, list(cells = 100)), "genes.cells must be a count of open cells below")
  refused(utils::modifyList(genes_record, list(cells = 0)), "genes.cells must be a count of open cells below")
  refused(c(genes_record, list(extra = 1)), "genes must have exactly cells, order and digest")
  refused(utils::modifyList(genes_record, list(order = list(method = "uniform", seed = 0, blockSize = 19))), "genes.order must be")
  refused(utils::modifyList(genes_record, list(order = list(method = "blocks", seed = -1, blockSize = 19))), "genes.order.seed must be a uint32")
  refused(utils::modifyList(genes_record, list(order = list(method = "blocks", seed = 0, blockSize = 0))), "genes.order.blockSize must be a positive count")
  refused(utils::modifyList(genes_record, list(digest = "ABC")), "genes.digest must be 16 lowercase hex digits")
})

test_that("color_by says the cells a gene is read for, and only for a gene", {
  event <- cellucid:::.parse_viewer_event(
    "{\"type\":\"color_by\",\"field\":\"ENSG00000168484\",\"source\":\"var\",\"symbol\":\"SFTPC\",\"cells_read\":40000}"
  )
  expect_identical(as.numeric(event$payload$cells_read), 40000)
  refused <- function(text) expect_error(cellucid:::.parse_viewer_event(text), class = "cellucid_protocol_error")
  refused("{\"type\":\"color_by\",\"field\":\"leiden\",\"source\":\"obs\",\"symbol\":null,\"cells_read\":3}")
  refused("{\"type\":\"color_by\",\"field\":\"g\",\"source\":\"var\",\"symbol\":null,\"cells_read\":-1}")
  refused("{\"type\":\"color_by\",\"field\":\"g\",\"source\":\"var\",\"symbol\":null}")
})
