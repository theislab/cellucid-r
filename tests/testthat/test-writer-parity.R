# Every float32 payload is rounded exactly once, at the write, and the
# metadata that names it is spelled the way the other writer spells it.
#
# cellucid-python and cellucid-r publish one format from one input, so a value
# that differs between them is a defect in whichever one rounded twice. These
# tests pin this writer's half of that contract: the caller's own numbers are
# scaled at their own precision, the product is rounded once, and the centroid
# published as a JSON double is measured from coordinates that were never
# rounded.

.parity_float32 <- function(values) {
  readBin(
    writeBin(as.double(values), raw(), size = 4L, endian = "little"),
    what = "numeric",
    n = length(values),
    size = 4L,
    endian = "little"
  )
}

.parity_read_float32 <- function(path, n) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = "numeric", size = 4L, n = n, endian = "little")
}

# Coordinates and vectors with every float64 mantissa bit occupied, so
# rounding the input before scaling changes the published value for a large
# share of the components instead of for none.
.parity_fixture <- function(n_cells = 96L) {
  set.seed(20260730L)
  umap2 <- matrix(
    stats::rnorm(n_cells * 2L, mean = 3.7182818284590452, sd = 4.1231056256176605),
    ncol = 2L
  )
  vectors <- matrix(
    stats::rnorm(n_cells * 2L, sd = 0.3183098861837907),
    ncol = 2L
  )
  list(
    n_cells = n_cells,
    umap2 = umap2,
    vectors = vectors,
    obs = data.frame(
      cluster = factor(rep(c("A", "B", "C", "D"), length.out = n_cells))
    ),
    latent = matrix(
      stats::rnorm(n_cells * 3L),
      ncol = 3L
    )
  )
}

.parity_scale_factor <- function(coords) {
  axis_mins <- apply(coords, 2, min)
  axis_maxs <- apply(coords, 2, max)
  2 / max(axis_maxs - axis_mins)
}

.parity_export <- function(fixture, out) {
  unlink(out, recursive = TRUE, force = TRUE)
  cellucid_prepare(
    dataset_id = "writer-parity",
    dataset_name = "Writer parity",
    latent_space = fixture$latent,
    obs = fixture$obs,
    X_umap_2d = fixture$umap2,
    vector_fields = list(velocity_umap_2d = fixture$vectors),
    out_dir = out,
    centroid_min_points = 4L,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )
  out
}

test_that("a vector payload is the scaled value rounded to float32 once", {
  fixture <- .parity_fixture()
  out <- .parity_export(fixture, file.path(tempdir(), "cellucid_parity_vectors"))

  scale_factor <- .parity_scale_factor(fixture$umap2)
  raw_values <- as.vector(t(fixture$vectors))
  rounded_once <- .parity_float32(raw_values * scale_factor)
  rounded_twice <- .parity_float32(.parity_float32(raw_values) * scale_factor)

  published <- .parity_read_float32(
    file.path(out, "vectors", "0_2d.bin"),
    length(raw_values)
  )

  expect_identical(published, rounded_once)
  # The fixture is only evidence if the two rounding models disagree on it.
  expect_gt(sum(rounded_once != rounded_twice), 0L)
  expect_false(identical(published, rounded_twice))
})

test_that("a points payload is the normalized value rounded to float32 once", {
  fixture <- .parity_fixture()
  out <- .parity_export(fixture, file.path(tempdir(), "cellucid_parity_points"))

  axis_mins <- apply(fixture$umap2, 2, min)
  axis_maxs <- apply(fixture$umap2, 2, max)
  center <- (axis_mins + axis_maxs) / 2
  scale_factor <- 2 / max(axis_maxs - axis_mins)
  normalized <- sweep(fixture$umap2, 2, center, FUN = "-") * scale_factor

  published <- .parity_read_float32(
    file.path(out, "points_2d.bin"),
    length(normalized)
  )
  expect_identical(published, .parity_float32(as.vector(t(normalized))))
})

test_that("a centroid is measured from coordinates that were never rounded", {
  fixture <- .parity_fixture()
  out <- .parity_export(fixture, file.path(tempdir(), "cellucid_parity_centroids"))

  axis_mins <- apply(fixture$umap2, 2, min)
  axis_maxs <- apply(fixture$umap2, 2, max)
  center <- (axis_mins + axis_maxs) / 2
  scale_factor <- 2 / max(axis_maxs - axis_mins)
  coords <- sweep(fixture$umap2, 2, center, FUN = "-") * scale_factor

  manifest <- jsonlite::fromJSON(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  centroids <- manifest[["_categoricalFields"]][[1L]][[6L]][["2"]]
  expect_length(centroids, 4L)

  codes <- as.integer(fixture$obs$cluster) - 1L
  for (code in seq_along(centroids) - 1L) {
    points <- coords[codes == code, , drop = FALSE]
    exact <- colMeans(points)
    distances <- sqrt(rowSums(sweep(points, 2, exact, FUN = "-")^2))
    threshold <- as.numeric(
      stats::quantile(distances, probs = 0.95, names = FALSE, type = 7)
    )
    inliers <- points[distances <= threshold, , drop = FALSE]
    expected <- if (nrow(inliers) >= 4L) colMeans(inliers) else exact

    published <- as.numeric(unlist(centroids[[code + 1L]][["position"]]))
    expect_identical(published, as.numeric(expected))

    # The same measurement taken from float32-rounded coordinates is a
    # different number, so the assertion above is not vacuous.
    rounded <- matrix(
      .parity_float32(as.vector(points)),
      nrow = nrow(points)
    )
    rounded_exact <- colMeans(rounded)
    rounded_distances <- sqrt(
      rowSums(sweep(rounded, 2, rounded_exact, FUN = "-")^2)
    )
    rounded_threshold <- as.numeric(
      stats::quantile(rounded_distances, probs = 0.95, names = FALSE, type = 7)
    )
    rounded_inliers <- rounded[rounded_distances <= rounded_threshold, , drop = FALSE]
    expect_false(identical(published, as.numeric(colMeans(rounded_inliers))))
  }
})

test_that("a vector field entry uses the key order both writers publish", {
  fixture <- .parity_fixture()
  out <- .parity_export(fixture, file.path(tempdir(), "cellucid_parity_identity"))

  identity <- jsonlite::fromJSON(
    file.path(out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  entry <- identity$vector_fields$fields$velocity_umap
  # cellucid-python emits exactly this order, and the output format
  # specification prints exactly this order.
  expect_identical(
    names(entry),
    c("label", "available_dimensions", "default_dimension", "files", "basis")
  )
})

test_that("an identifier defect is reported in the words the other writer uses", {
  # `what` is a singular noun in every caller, so the position frame composes
  # as a sentence: "Gene identifier at position 1 ...", not "Gene identifiers
  # identifier at position 1 ...". cellucid-python composes the identical
  # sentence from the identical nouns.
  expect_error(
    cellucid:::.require_field_identities(c("ok", "Liver "), "Gene"),
    paste0(
      "^Gene identifier at position 1 is displayed verbatim, so it must not ",
      "carry characters that have no glyph"
    )
  )
  expect_error(
    cellucid:::.require_field_identities(c("a", "a"), "Observation field"),
    "^Observation field key 'a' is duplicated\\.$"
  )
  expect_error(
    cellucid:::.require_field_identities(c("a", "a"), "Vector field"),
    "^Vector field key 'a' is duplicated\\.$"
  )
  # A set of values is shown as a list, with the boundary that lets a reader
  # tell where it ends. cellucid-python prints the same list.
  expect_error(
    cellucid:::.require_dense_payload_indices(list(0L, 1L, 1L), axis = "Gene"),
    "^Gene payload indices must be exactly 0\\.\\.2, each used once; got \\[0, 1, 1\\]\\.$"
  )
  expect_identical(cellucid:::.format_value_list(c("a", "b")), "['a', 'b']")
  expect_identical(cellucid:::.format_value_list(c(0L, 2L)), "[0, 2]")
})
