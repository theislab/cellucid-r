# Latent-space outlier quantiles are the bytes cellucid-python and the web app
# compute, bit for bit. The fixture is written by cellucid-python's
# scripts/generate_latent_outlier_conformance.py and kept byte for byte here
# and in the web app under the same pin.

.latent_outlier_sha256 <- "46ea30291fc2d1fc51fc2021837eb8914062f8080b6ad51995cacaa1875b9d80"

.float32_hex <- function(values) {
  bytes <- writeBin(values, raw(), size = 4L, endian = "little")
  paste(format(as.hexmode(as.integer(bytes)), width = 2L), collapse = "")
}

test_that("the shared conformance cases give the pinned bytes", {
  expect_pinned_fixture("latent_outlier_conformance.json", .latent_outlier_sha256)
  document <- cellucid_read_fixture("latent_outlier_conformance.json")
  expect_length(document$cases, 4L)
  for (case in document$cases) {
    values <- vapply(case$latent, as.numeric, numeric(1))
    latent <- matrix(values, ncol = case$width, byrow = TRUE)
    codes <- vapply(case$codes, as.integer, integer(1))
    quantiles <- cellucid:::.compute_latent_space_quantiles(
      latent,
      codes,
      categories = seq_len(case$nCategories),
      min_points = case$minPoints
    )
    expect_identical(.float32_hex(quantiles), case$quantilesHex, label = case$name)
  }
})

test_that("ties share the larger rank; a small category and a missing code are NaN", {
  latent <- matrix(c(0, 1, 1, 4, 7, 8, 3, 0, 0, 0, 0, 7, 8, 3), ncol = 2L)
  quantiles <- cellucid:::.compute_latent_space_quantiles(
    latent,
    c(0L, 0L, 0L, 0L, 1L, 1L, NA_integer_),
    categories = c("a", "b"),
    min_points = 3L
  )
  expect_identical(quantiles[1:4], c(0.75, 0.5, 0.5, 1))
  expect_true(all(is.nan(quantiles[5:7])))
})
