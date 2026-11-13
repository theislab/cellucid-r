# cellucid_subsample_cells() follows the subsampling specification: every case
# of the shared conformance fixture keeps exactly the cells it lists.
#
# tests/testthat/fixtures/subsample_conformance.json is a byte-identical copy
# of the web app's tests/fixtures/subsample-conformance.json, computed in
# plain Python by its scripts/generate-subsample-conformance.py independently
# of every implementation; the test compares the copy with the web
# repository's file, as cellucid-python's does. Positions in the fixture are
# 0-based; the function returns R's 1-based positions.

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

test_that("the fixture is the web app's, byte for byte", {
  expect_shared_file("tests/testthat/fixtures/subsample_conformance.json", "tests/fixtures/subsample-conformance.json")
})

test_that("the fixture is of the version this suite reads", {
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

test_that("a number's text is the one JavaScript's String() and Python's repr give, ties to the even digit, on every platform", {
  # Found in exact integer arithmetic (src/number_text.c), not by the C
  # library's printf, which breaks a tie as each likes: Windows' writes 2^-25
  # as "...313", and the fixture's 2^-25 of 50,331,648 cells would come out 2.
  # The texts are node's (JSON.stringify); the same digits a share's
  # subsample size is taken from, and a session value is hashed by.
  text <- function(x) cellucid:::.native_number_text(x)
  cases <- list(
    list(2^-25, "2.9802322387695312e-8"), list(2^-24, "5.960464477539063e-8"),
    list(2^-53, "1.1102230246251565e-16"), list(2^-60, "8.673617379884035e-19"), list(2^-1074, "5e-324"),
    list(2^-1022, "2.2250738585072014e-308"), list(2^-1022 - 2^-1074, "2.225073858507201e-308"),
    list(0.1, "0.1"), list(1 / 3, "0.3333333333333333"),
    list(0x1.9769d41p-4, "0.09946616017259657"), # R reads "0.0994...57" one ulp low on some builds
    list(1 - 2^-53, "0.9999999999999999"), list(1, "1"), list(1e21, "1e+21"), list(1e-7, "1e-7"),
    list(1.5e-7, "1.5e-7"), list(2^66, "73786976294838210000"), list(2^70, "1.1805916207174113e+21"),
    list(2^66 + 2^14, "73786976294838220000"), list(-0, "0"),
    list(.Machine$double.xmax, "1.7976931348623157e+308"), list(0.000001, "0.000001"), list(-3.875, "-3.875"),
    list(4294967295, "4294967295"), list(2^53 + 2, "9007199254740994")
  )
  for (case in cases) expect_identical(text(case[[1L]]), case[[2L]], info = case[[2L]])
  expect_identical(text(c(5L, -2L, 0L, .Machine$integer.max, -.Machine$integer.max)),
                   c("5", "-2", "0", "2147483647", "-2147483647"))
  expect_error(text(NA_real_), "finite")
  expect_error(text(Inf), "finite")
})

test_that("a number's text at the extremes of a double is node's, each double named by its IEEE 754 bits", {
  # Named by bits, since R's own reading of a literal is not correctly rounded
  # on every build. The texts are node's, generated from these bits; among
  # them the doubles a decimal one digit shorter reaches by rounding up to the
  # next power of ten (1e23, 1e-7, 1e-323: the double lies just below it).
  double_of <- function(bits) {
    bytes <- as.raw(strtoi(substring(bits, seq(1L, 15L, 2L), seq(2L, 16L, 2L)), 16L))
    readBin(bytes, "double", 1L, size = 8L, endian = "big")
  }
  text <- function(x) cellucid:::.native_number_text(x)
  cases <- c(
    "7fefffffffffffff" = "1.7976931348623157e+308",
    "ffefffffffffffff" = "-1.7976931348623157e+308",
    "0000000000000001" = "5e-324",
    "8000000000000001" = "-5e-324",
    "0010000000000000" = "2.2250738585072014e-308",
    "000fffffffffffff" = "2.225073858507201e-308",
    "444b1ae4d6e2ef50" = "1e+21",
    "3e7ad7f29abcaf48" = "1e-7",
    "3eb0c6f7a0b5ed8d" = "0.000001",
    "444b1ae4d6e2ef4f" = "999999999999999900000",
    "444b1ae4d6e2ef51" = "1.0000000000000001e+21",
    "44b52d02c7e14af6" = "1e+23",
    "4480f0cf064dd592" = "1e+22",
    "0000000000000002" = "1e-323",
    "0000000000000004" = "2e-323",
    "0000000000000014" = "1e-322",
    "0000000000000013" = "9.4e-323",
    "7e37e43c8800759c" = "1e+300",
    "7fe1ccf385ebc8a0" = "1e+308",
    "01a56e1fc2f8f359" = "1e-300",
    "000730d67819e8d2" = "1e-308",
    "000012688b70e62b" = "1e-310",
    "44b52d02c7e14af7" = "1.0000000000000001e+23",
    "3ea0c6f7a0b5ed8d" = "5e-7",
    "3e8421f5f40d8376" = "1.5e-7",
    "3fb999999999999a" = "0.1",
    "3fc999999999999a" = "0.2",
    "3fd3333333333334" = "0.30000000000000004",
    "3fd5555555555555" = "0.3333333333333333",
    "3fe5555555555555" = "0.6666666666666666",
    "40fe240c9fbe76c9" = "123456.789",
    "4011666666666666" = "4.35",
    "4004000000000000" = "2.5",
    "3eb4b3fd5942cd96" = "0.000001234",
    "7feffffffffffffe" = "1.7976931348623155e+308",
    "3e60000000000000" = "2.9802322387695312e-8",
    "0008000000000000" = "1.1125369292536007e-308",
    "7fe0000000000000" = "8.98846567431158e+307",
    "4450000000000000" = "1.1805916207174113e+21",
    "4340000000000000" = "9007199254740992",
    "4340000000000001" = "9007199254740994",
    "4350000000000000" = "18014398509481984",
    "43b0000000000000" = "1152921504606847000",
    "43e0000000000000" = "9223372036854776000",
    "43f0000000000000" = "18446744073709552000",
    "2b2bff2ee48e0530" = "1e-100",
    "358dee7a4ad4b81f" = "1e-50",
    "3bc79ca10c924223" = "1e-20",
    "3ddb7cdfd9d7bdbb" = "1e-10",
    "3ee4f8b588e368f1" = "0.00001",
    "3f50624dd2f1a9fc" = "0.001",
    "408f400000000000" = "1000",
    "430c6bf526340000" = "1000000000000000",
    "4341c37937e08000" = "10000000000000000",
    "4376345785d8a000" = "100000000000000000",
    "43abc16d674ec800" = "1000000000000000000",
    "43e158e460913d00" = "10000000000000000000",
    "4415af1d78b58c40" = "100000000000000000000",
    "45208b2a2c280291" = "1e+25",
    "46293e5939a08cea" = "1e+30",
    "483d6329f1c35ca5" = "1e+40",
    "54b249ad2594c37d" = "1e+100",
    "6974e718d7d7625a" = "1e+200"
  )
  for (bits in names(cases)) {
    expect_identical(text(double_of(bits)), cases[[bits]], info = bits)
  }
  # The extremes as R names them.
  expect_identical(text(.Machine$double.xmax), "1.7976931348623157e+308")
  expect_identical(text(-.Machine$double.xmax), "-1.7976931348623157e+308")
  expect_identical(text(.Machine$double.xmin), "2.2250738585072014e-308")
  expect_identical(text(.Machine$double.eps), "2.220446049250313e-16")
  # A whole number written by hand: every length up to 2^53's, and its ends.
  powers <- cumprod(rep(10, 15L))
  for (digits in seq_along(powers)) {
    expect_identical(text(powers[[digits]] - 1), strrep("9", digits), info = digits)
    expect_identical(text(powers[[digits]]), paste0("1", strrep("0", digits)), info = digits)
  }
  expect_identical(text(2^53 - 1), "9007199254740991")
  expect_identical(text(2^53), "9007199254740992")
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
