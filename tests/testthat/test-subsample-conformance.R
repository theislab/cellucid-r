# cellucid_subsample_cells() follows the subsampling specification: every case
# of the shared conformance fixture keeps exactly the cells it lists.
#
# tests/testthat/fixtures/subsample_conformance.json is a byte-identical copy
# of the web app's tests/fixtures/subsample-conformance.json, computed in
# plain Python by its scripts/generate-subsample-conformance.py independently
# of every implementation; the web app and cellucid-python pin the same
# SHA-256. Positions in the fixture are 0-based; the function returns R's
# 1-based positions.

.subsample_sha256 <- "a805a302a88e77e0c7e07e998811ce79825b5a34890db79a9462e801ab84a60d"

.subsample_case <- function(case) {
  if (identical(case$method, "uniform")) {
    return(cellucid_subsample_cells(case$nObs, case$k, seed = case$seed))
  }
  if (identical(case$method, "blocks")) {
    return(cellucid_subsample_cells(case$nObs, case$k, method = "blocks", seed = case$seed, block_size = case$blockSize))
  }
  if (identical(case$method, "stratified")) {
    codes <- vapply(case$codes, function(code) if (is.null(code)) -1L else as.integer(code), integer(1))
    groups <- factor(ifelse(codes < 0L, NA, paste0("c", codes)), levels = paste0("c", seq_len(case$nCategories) - 1L))
    return(cellucid_subsample_cells(case$nObs, case$k, method = "stratified", seed = case$seed, groupby = groups,
                                    allocation = case$allocation, min_per_group = case$minPerCategory))
  }
  values <- vapply(case$coordinates, function(value) if (is.null(value)) NaN else as.numeric(value), numeric(1))  # "Infinity", "-0" too
  coordinates <- matrix(values, ncol = case$dimensions, byrow = TRUE)
  cellucid_subsample_cells(case$nObs, case$k, method = "spatial", seed = case$seed, basis = coordinates,
                           columns = seq_len(case$dimensions), grid_size = case$gridSize)
}

test_that("the fixture is the pinned one", {
  expect_pinned_fixture("subsample_conformance.json", .subsample_sha256)
  expect_identical(cellucid_read_fixture("subsample_conformance.json")$version, 1L)
})

test_that("ranks and sizes are the specification's", {
  document <- cellucid_read_fixture("subsample_conformance.json")
  for (item in document$rank) {
    positions <- vapply(item$positions, as.numeric, numeric(1))
    expect_identical(cellucid:::.subsample_rank(item$seed, positions), vapply(item$ranks, as.numeric, numeric(1)))
  }
  for (item in document$refusals) {
    expect_error(.subsample_case(item), "cannot hold", label = item$name)
  }
  for (item in document$fraction) {
    expect_identical(cellucid:::.subsample_count_from_fraction(item$fraction, item$nObs), as.numeric(item$k))
  }
})

test_that("every case keeps the cells the fixture lists, with its digest", {
  document <- cellucid_read_fixture("subsample_conformance.json")
  for (case in document$cases) {
    kept <- .subsample_case(case) - 1
    if (is.null(case$kept) && is.null(case$keptCount)) {
      expect_identical(kept, as.numeric(seq_len(case$nObs) - 1L), label = case$name)
      expect_identical(cellucid:::.subsample_digest(case$nObs, NULL), case$digest, label = case$name)
      next
    }
    if (!is.null(case$kept)) {
      expect_identical(kept, vapply(case$kept, as.numeric, numeric(1)), label = case$name)
    } else {
      expect_identical(length(kept), as.integer(case$keptCount), label = case$name)
      expect_identical(head(kept, 8), vapply(case$keptHead, as.numeric, numeric(1)), label = case$name)
      expect_identical(tail(kept, 8), vapply(case$keptTail, as.numeric, numeric(1)), label = case$name)
    }
    expect_identical(cellucid:::.subsample_digest(case$nObs, kept), case$digest, label = case$name)
  }
})

test_that("the session encoding of a cell set is the fixture's", {
  document <- cellucid_read_fixture("subsample_conformance.json")
  for (item in document$encoding) {
    if (is.null(item$kept)) {
      expect_identical(cellucid:::.subsample_digest(item$nObs, NULL), item$digest)
      next
    }
    kept <- vapply(item$kept, as.numeric, numeric(1))
    bytes <- cellucid:::.subsample_encode_cell_set(kept)
    expect_identical(paste(format(as.hexmode(as.integer(bytes)), width = 2L), collapse = ""), item$bytesHex)
    expect_identical(cellucid:::.subsample_decode_cell_set(bytes, item$nObs), kept)
    expect_identical(cellucid:::.subsample_digest(item$nObs, kept), item$digest)
  }
})

test_that("uniform and blocks are nested, and requests outside the specification are refused", {
  previous <- NULL
  for (k in c(3, 40, 900, 4999)) {
    kept <- cellucid_subsample_cells(5000, k, seed = 11)
    if (!is.null(previous)) expect_true(all(previous %in% kept))
    previous <- kept
  }
  expect_error(cellucid_subsample_cells(10, 0))
  expect_error(cellucid_subsample_cells(10, 3, seed = -1))
  expect_error(cellucid_subsample_cells(10, 3, seed = 2^32))
  expect_error(cellucid_subsample_cells(10, fraction = 0))
  expect_error(cellucid_subsample_cells(10, 3, fraction = 0.5))
  expect_error(cellucid_subsample_cells(10, 3, method = "blocks"))
  expect_error(cellucid_subsample_cells(4, 2, method = "stratified", groupby = factor(c("a", "b", "c", "a")), min_per_group = 1))
  expect_identical(cellucid_subsample_cells(5, 9), 1:5)
})

test_that("a session's subsample cells are read as the web app writes them", {
  kept <- cellucid_subsample_cells(10000, 1234, method = "blocks", seed = 9, block_size = 64) - 1
  record <- list(
    spec = list(method = "blocks", nObs = 10000L, k = 1234L, seed = 9L, blockSize = 64L),
    kept = gsub("\n", "", jsonlite::base64_enc(cellucid:::.subsample_encode_cell_set(kept)), fixed = TRUE),
    digest = cellucid:::.subsample_digest(10000, kept)
  )
  decoded <- cellucid:::.session_cells(record, 1234L)
  expect_identical(decoded$kept, kept)
  record$digest <- strrep("0", 16)
  expect_error(cellucid:::.session_cells(record, 1234L), "digest")
})

test_that("categories of plain values are ordered by code point, never by the locale", {
  values <- rep(c("a", "B"), each = 5)
  by_code_point <- factor(values, levels = c("B", "a"))
  expect_identical(
    cellucid_subsample_cells(10, 3, method = "stratified", groupby = values, allocation = "balanced"),
    cellucid_subsample_cells(10, 3, method = "stratified", groupby = by_code_point, allocation = "balanced")
  )
})
