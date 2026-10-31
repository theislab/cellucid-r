# cellucid_subsample_cells() follows the subsampling specification: every case
# of the shared conformance fixture keeps exactly the cells it lists.
#
# tests/testthat/fixtures/subsample_conformance.json is a byte-identical copy
# of the web app's tests/fixtures/subsample-conformance.json, computed in
# plain Python by its scripts/generate-subsample-conformance.py independently
# of every implementation; the web app and cellucid-python pin the same
# SHA-256. Positions in the fixture are 0-based; the function returns R's
# 1-based positions.

.subsample_sha256 <- "ee4b8462f66e6e441c2ab9f0948ca26d04451aa48017c38e577e3f13c8f1f48b"

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

test_that("a fraction's shortest decimal is the one Python's repr gives, ties to the even digit, on every platform", {
  # Found in exact integer arithmetic (src/subsample.c), not by the C
  # library's printf, which breaks a tie as each likes: Windows' wrote 2^-25
  # as "...313", and the fixture's 2^-25 of 50,331,648 cells then came out 2.
  shortest <- function(x) cellucid:::.native_shortest_decimal(x)
  expect_identical(shortest(2^-25), "2.9802322387695312e-08")
  expect_identical(shortest(2^-24), "5.960464477539063e-08")
  expect_identical(shortest(2^-53), "1.1102230246251565e-16")
  expect_identical(shortest(2^-60), "8.673617379884035e-19")
  expect_identical(shortest(2^-1074), "5e-324")
  expect_identical(shortest(2^-1022), "2.2250738585072014e-308")
  expect_identical(shortest(2^-1022 - 2^-1074), "2.225073858507201e-308")
  expect_identical(shortest(0.1), "1e-01")
  expect_identical(shortest(1 / 3), "3.333333333333333e-01")
  expect_identical(shortest(0x1.9769d41p-4), "9.946616017259657e-02") # R reads "0.0994...57" one ulp low on some builds
  expect_identical(shortest(1 - 2^-53), "9.999999999999999e-01")
  expect_identical(shortest(1), "1e+00")
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

# An object in a Cellucid-ready file's stored order, by hand as a reader gives
# it: rows (columns) by ascending rank of `seed`, each one's original position
# recorded, and the record beside them.
.stored_cells <- function(frame, seed) {
  n <- nrow(frame)
  positions <- order(cellucid:::.subsample_rank(seed, seq_len(n) - 1), method = "radix") - 1
  stored <- frame[positions + 1, , drop = FALSE]
  stored$cellucid_original_position <- positions
  attr(stored, "cellucid") <- list(priority = list(seed = seed, specVersion = 1L, nObs = n, levels = cellucid:::.prefix_levels(n)))
  stored
}

test_that("an object in a file's stored order is refused until restored, then keeps the fixture's cells", {
  # The viewer opens the fixture's cells (original positions) of a file in
  # priority order; an object holding that file's rows as stored would take
  # other cells by position, as pull, subset and a session would.
  document <- cellucid_read_fixture("subsample_conformance.json")
  for (case in document$cases) {
    if (is.null(case$kept) || identical(case$method, "spatial") || case$nObs < 2) next
    n <- case$nObs
    frame <- data.frame(row.names = paste0("cell", seq_len(n) - 1L))
    request <- list(n_cells = case$k, method = case$method, seed = case$seed)
    if (identical(case$method, "stratified")) {
      codes <- vapply(case$codes, function(code) if (is.null(code)) -1L else as.integer(code), integer(1))
      frame$group <- factor(ifelse(codes < 0L, NA, paste0("c", codes)), levels = paste0("c", seq_len(case$nCategories) - 1L))
      request <- c(request, list(groupby = "group", allocation = case$allocation, min_per_group = case$minPerCategory))
    }
    if (identical(case$method, "blocks")) request <- c(request, list(block_size = case$blockSize))
    stored <- .stored_cells(frame, 7)
    expect_error(do.call(cellucid_subsample_cells, c(list(stored), request)), "cellucid_restore_original_order", label = case$name)
    restored <- cellucid_restore_original_order(stored)
    kept <- do.call(cellucid_subsample_cells, c(list(restored), request))
    expected <- vapply(case$kept, as.numeric, numeric(1))
    expect_identical(as.numeric(kept) - 1, expected, label = case$name)
    expect_identical(rownames(restored)[kept], paste0("cell", expected), label = case$name)
  }
  skip_without_sce()
  case <- Filter(function(case) identical(case$name, "spatial a constant axis"), document$cases)[[1L]]
  coordinates <- matrix(vapply(case$coordinates, as.numeric, numeric(1)), ncol = case$dimensions, byrow = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 1, ncol = case$nObs, dimnames = list("g", paste0("cell", seq_len(case$nObs) - 1L)))),
    reducedDims = list(SPACE = coordinates)
  )
  positions <- order(cellucid:::.subsample_rank(7, seq_len(case$nObs) - 1), method = "radix") - 1
  stored <- sce[, positions + 1]
  SummarizedExperiment::colData(stored)$cellucid_original_position <- positions
  S4Vectors::metadata(stored)$cellucid <- list(priority = list(seed = 7, specVersion = 1L, nObs = case$nObs, levels = cellucid:::.prefix_levels(case$nObs)))
  request <- list(n_cells = case$k, method = "spatial", seed = case$seed, basis = "SPACE", columns = seq_len(case$dimensions), grid_size = case$gridSize)
  expect_error(do.call(cellucid_subsample_cells, c(list(stored), request)), "cellucid_restore_original_order")
  kept <- do.call(cellucid_subsample_cells, c(list(cellucid_restore_original_order(stored)), request))
  expect_identical(as.numeric(kept) - 1, vapply(case$kept, as.numeric, numeric(1)))
})

test_that("a record that no longer describes the cells warns and takes each row as its own position", {
  stored <- .stored_cells(data.frame(row.names = paste0("c", 1:200)), 0)
  subset <- stored[1:150, , drop = FALSE]
  attr(subset, "cellucid") <- attr(stored, "cellucid")
  expect_warning(kept <- cellucid_subsample_cells(subset, 15, seed = 3), "written for 200 cells")
  expect_identical(kept, cellucid_subsample_cells(150, 15, seed = 3))
})
