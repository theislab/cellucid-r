# Latent-space outlier quantiles are the bytes cellucid-python's prepare()
# stores and the web app computes, bit for bit.
#
# tests/testthat/fixtures/latent_outlier_conformance.json is a byte-identical
# copy of the web app's tests/fixtures/latent-outlier-conformance.json, written
# by its scripts/generate-latent-outlier-conformance.py with numpy,
# independently of every implementation; both repositories pin the same
# SHA-256.

.latent_outlier_sha256 <- "a812db4dd25c9b5b44549c001a0e3927043042ca12ad7d8567c73f583ecf498d"

.float32_hex <- function(values) {
  bytes <- writeBin(values, raw(), size = 4L, endian = "little")
  # The fixture writes numpy's NaN; R's NaN is the platform's 0/0.
  nan <- which(is.nan(values))
  if (length(nan) > 0L) {
    bytes[as.vector(outer(1:4, (nan - 1L) * 4L, "+"))] <- as.raw(c(0x00, 0x00, 0xc0, 0x7f))
  }
  paste(format(as.hexmode(as.integer(bytes)), width = 2L), collapse = "")
}

.quantiles <- function(latent, codes, n_categories, min_points) {
  cellucid:::.native_latent_outlier_quantiles(latent, codes, as.integer(n_categories), as.integer(min_points))
}

test_that("the shared conformance cases give the pinned bytes", {
  expect_pinned_fixture("latent_outlier_conformance.json", .latent_outlier_sha256)
  document <- cellucid_read_fixture("latent_outlier_conformance.json")
  expect_identical(document$version, 2L)
  expect_length(document$cases, 5L)
  for (case in document$cases) {
    # null: a NaN row's value (a cell the latent space does not place).
    values <- vapply(case$latent, function(value) if (is.null(value)) NaN else as.numeric(value), numeric(1))
    latent <- matrix(values, ncol = case$width, byrow = TRUE)
    codes <- vapply(case$codes, as.integer, integer(1))
    codes[codes < 0L] <- NA_integer_
    quantiles <- .quantiles(latent, codes, case$nCategories, case$minPoints)
    expect_identical(.float32_hex(quantiles), case$quantilesHex, label = case$name)
  }
})

test_that("ties share the larger rank; a small category and a missing code are NaN", {
  latent <- matrix(c(0, 1, 1, 4, 7, 8, 3, 0, 0, 0, 0, 7, 8, 3), ncol = 2L)
  quantiles <- .quantiles(latent, c(0L, 0L, 0L, 0L, 1L, 1L, NA_integer_), 2L, 3L)
  expect_identical(quantiles[1:4], c(0.75, 0.5, 0.5, 1))
  expect_true(all(is.nan(quantiles[5:7])))
})

test_that("the quantiles use the latent values as float32 stores them", {
  # 1 + 2^-30 is 1 in float32: the first two cells are at the same distance
  # (2/3) from the centre (5/3), so both have 2 of the 3 cells at it or closer.
  # In double, their distances would differ and so would their quantiles.
  latent <- matrix(c(1 + 2^-30, 1, 3), ncol = 1L)
  expect_identical(.quantiles(latent, c(0L, 0L, 0L), 1L, 1L), as_float32(c(2, 2, 3) / 3))
})

test_that("an obs column is partitioned as the viewer shows it", {
  partition <- cellucid:::.outlier_partition
  expect_identical(
    partition(factor(c("b", NA, "a"), levels = c("a", "b", "c"))),
    list(codes = c(1L, NA, 0L), n_categories = 3L)
  )
  expect_identical(partition(c(TRUE, NA, FALSE)), list(codes = c(1L, NA, 0L), n_categories = 2L))
  expect_identical(partition(c("x", NA, "y", "x")), list(codes = c(0L, NA, 1L, 0L), n_categories = 2L))
  expect_null(partition(c(1, 2, 3)))
  expect_null(partition(1:3))
})
