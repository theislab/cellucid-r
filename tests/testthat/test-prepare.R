test_that("cellucid_prepare writes expected core files", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(
    cluster = factor(c("A", "A", "B")),
    score = c(0.1, 0.2, 0.3)
  )
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  expr <- matrix(c(0, 1, 2, 3, 4, 5), nrow = 3, ncol = 2)
  var <- data.frame(symbol = c("G1", "G2"))
  rownames(var) <- var$symbol

  out <- file.path(tempdir(), "cellucid_test_export")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    latent_space = latent,
    obs = obs,
    var = var,
    gene_expression = expr,
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE
  )

  expect_true(file.exists(file.path(out, "points_2d.bin")))
  expect_true(file.exists(file.path(out, "obs_manifest.json")))
  expect_true(file.exists(file.path(out, "var_manifest.json")))
  expect_true(file.exists(file.path(out, "dataset_identity.json")))

  obs_dir <- file.path(out, "obs")
  expect_true(file.exists(file.path(obs_dir, "cluster.codes.u8")))
  expect_true(file.exists(file.path(obs_dir, "cluster.outliers.f32")))
  expect_true(file.exists(file.path(obs_dir, "score.values.f32")))

  var_dir <- file.path(out, "var")
  expect_true(file.exists(file.path(var_dir, "G1.values.f32")))
  expect_true(file.exists(file.path(var_dir, "G2.values.f32")))
})

test_that("points_2d.bin is float32 row-major and normalized", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))
  umap2 <- matrix(c(0, 0, 2, 0, 0, 1), ncol = 2, byrow = TRUE)

  out <- file.path(tempdir(), "cellucid_test_points")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    latent_space = latent,
    obs = obs,
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE
  )

  # Expected normalization (matches cellucid-python prepare_data.py)
  axis_mins <- apply(umap2, 2, min)
  axis_maxs <- apply(umap2, 2, max)
  axis_ranges <- axis_maxs - axis_mins
  max_range <- max(axis_ranges)
  center <- (axis_mins + axis_maxs) / 2
  scale_factor <- 2 / max_range
  expected <- sweep(umap2, 2, center, "-") * scale_factor
  expected_flat <- as.vector(t(expected))

  path <- file.path(out, "points_2d.bin")
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  got <- readBin(con, what = "numeric", size = 4, n = length(expected_flat), endian = "little")

  expect_equal(got, expected_flat, tolerance = 1e-6)
})

test_that("continuous quantization uses reserved missing marker", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(score = c(0, NA, Inf))
  umap1 <- matrix(c(0, 1, 2), ncol = 1)

  out <- file.path(tempdir(), "cellucid_test_quant")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = out,
    obs_continuous_quantization = 8,
    centroid_min_points = 1,
    force = TRUE
  )

  vals_path <- file.path(out, "obs", "score.values.u8")
  expect_true(file.exists(vals_path))
  con <- file(vals_path, open = "rb")
  on.exit(close(con), add = TRUE)
  got <- readBin(con, what = "raw", n = 3)
  expect_equal(got, as.raw(c(0, 255, 255)))
})

test_that("connectivity export writes uint16 edge pairs", {
  skip_if_not_installed("Matrix")

  latent <- matrix(c(0, 0, 1, 1, 2, 2, 3, 3), ncol = 2, byrow = TRUE)
  obs <- data.frame(a = 1:4)
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)

  conn <- matrix(0, nrow = 4, ncol = 4)
  conn[1, 2] <- 1
  conn[1, 3] <- 1
  conn[4, 3] <- 1

  out <- file.path(tempdir(), "cellucid_test_connectivity")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    latent_space = latent,
    obs = obs,
    connectivities = conn,
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE
  )

  src_path <- file.path(out, "connectivity", "edges.src.bin")
  dst_path <- file.path(out, "connectivity", "edges.dst.bin")
  expect_true(file.exists(src_path))
  expect_true(file.exists(dst_path))

  read_u16 <- function(path) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    readBin(con, what = "integer", size = 2, endian = "little", n = 1000)
  }

  # Expected unique undirected edges: (0,1), (0,2), (2,3)
  expect_equal(read_u16(src_path)[1:3], c(0, 0, 2))
  expect_equal(read_u16(dst_path)[1:3], c(1, 2, 3))
})

test_that("vector fields are exported and scaled with embedding normalization", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))

  # Range is 4 on both axes -> scale_factor = 2 / 4 = 0.5
  umap2 <- matrix(c(0, 0,
                    4, 0,
                    0, 4),
                  ncol = 2, byrow = TRUE)

  vector_fields <- list(
    velocity_umap_2d = matrix(c(2, 0,
                                2, 0,
                                2, 0),
                              ncol = 2, byrow = TRUE)
  )

  out <- file.path(tempdir(), "cellucid_test_vectors")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    latent_space = latent,
    obs = obs,
    X_umap_2d = umap2,
    vector_fields = vector_fields,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE
  )

  vec_path <- file.path(out, "vectors", "velocity_umap_2d.bin")
  expect_true(file.exists(vec_path))

  con <- file(vec_path, open = "rb")
  on.exit(close(con), add = TRUE)
  got <- readBin(con, what = "numeric", size = 4, endian = "little", n = 3 * 2)

  # Input vectors are (2, 0); scale_factor=0.5 -> expected (1, 0)
  expect_equal(got, rep(c(1, 0), 3), tolerance = 1e-6)

  ident <- jsonlite::read_json(file.path(out, "dataset_identity.json"), simplifyVector = TRUE)
  expect_equal(ident$vector_fields$default_field, "velocity_umap")
  expect_true("velocity_umap" %in% names(ident$vector_fields$fields))
  expect_true("2d" %in% names(ident$vector_fields$fields$velocity_umap$files))
})

test_that("obs manifest preserves min/max precision in JSON", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(score = c(1 / 3, 2 / 3, NA_real_))
  umap1 <- matrix(c(0, 1, 2), ncol = 1)

  out <- file.path(tempdir(), "cellucid_test_json_precision")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = out,
    obs_continuous_quantization = 8,
    centroid_min_points = 1,
    force = TRUE
  )

  manifest <- jsonlite::read_json(file.path(out, "obs_manifest.json"), simplifyVector = FALSE)
  fields <- manifest[["_continuousFields"]]
  expect_equal(length(fields), 1L)

  entry <- fields[[1]]
  expect_equal(entry[[1]], "score")
  expect_equal(entry[[2]], 1 / 3, tolerance = 1e-15)
  expect_equal(entry[[3]], 2 / 3, tolerance = 1e-15)
})

test_that("cellucid_prepare errors on non-finite embedding values", {
  latent <- matrix(c(0, 0), ncol = 2)
  obs <- data.frame(cluster = factor("A"))
  umap2 <- matrix(c(0, NA), ncol = 2)

  expect_error(
    cellucid_prepare(
      latent_space = latent,
      obs = obs,
      X_umap_2d = umap2,
      out_dir = file.path(tempdir(), "cellucid_test_bad_embedding"),
      centroid_min_points = 1,
      force = TRUE
    ),
    "non-finite values"
  )
})

test_that("cellucid_prepare errors without embeddings", {
  latent <- matrix(c(0, 0), ncol = 2)
  obs <- data.frame(a = 1)
  expect_error(
    cellucid_prepare(latent_space = latent, obs = obs, X_umap_1d = NULL, X_umap_2d = NULL, X_umap_3d = NULL),
    "At least one embedding"
  )
})

test_that("X_umap_4d is rejected", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(a = 1:2)
  umap2 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  umap4 <- matrix(rep(0, 8), ncol = 4)
  expect_error(
    cellucid_prepare(latent_space = latent, obs = obs, X_umap_2d = umap2, X_umap_4d = umap4),
    "4D visualization is not yet implemented"
  )
})

test_that("duplicate gene identifiers are rejected", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)

  expr <- matrix(0, nrow = 3, ncol = 2)
  var <- data.frame(symbol = c("G1", "G1"))

  out <- file.path(tempdir(), "cellucid_test_dup_gene_ids")
  unlink(out, recursive = TRUE, force = TRUE)

  expect_error(
    cellucid_prepare(
      latent_space = latent,
      obs = obs,
      var = var,
      gene_expression = expr,
      var_gene_id_column = "symbol",
      X_umap_2d = umap2,
      out_dir = out,
      centroid_min_points = 1,
      force = TRUE
    ),
    "Gene identifiers must be unique"
  )
})

test_that("gene ids that collide after sanitization are rejected", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)

  expr <- matrix(0, nrow = 3, ncol = 2)
  var <- data.frame(symbol = c("A/B", "A B"))
  rownames(var) <- var$symbol

  out <- file.path(tempdir(), "cellucid_test_gene_id_collision")
  unlink(out, recursive = TRUE, force = TRUE)

  expect_error(
    cellucid_prepare(
      latent_space = latent,
      obs = obs,
      var = var,
      gene_expression = expr,
      X_umap_2d = umap2,
      out_dir = out,
      centroid_min_points = 1,
      force = TRUE
    ),
    "collide after filename sanitization"
  )
})

test_that("gene_identifiers subsets var export and dataset identity counts exported genes", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)

  expr <- matrix(c(0, 1, 2, 3, 4, 5), nrow = 3, ncol = 2)
  var <- data.frame(symbol = c("G1", "G2"))
  rownames(var) <- var$symbol

  out <- file.path(tempdir(), "cellucid_test_gene_subset")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    latent_space = latent,
    obs = obs,
    var = var,
    gene_expression = expr,
    gene_identifiers = c("G2"),
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE
  )

  var_dir <- file.path(out, "var")
  expect_false(file.exists(file.path(var_dir, "G1.values.f32")))
  expect_true(file.exists(file.path(var_dir, "G2.values.f32")))

  ident <- jsonlite::read_json(file.path(out, "dataset_identity.json"), simplifyVector = TRUE)
  expect_equal(ident$stats$n_genes, 1L)
})

test_that("obs keys that collide after sanitization are rejected", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(
    "a/b" = factor(c("A", "A", "B")),
    "a b" = factor(c("A", "A", "B")),
    check.names = FALSE
  )
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)

  out <- file.path(tempdir(), "cellucid_test_obs_key_collision")
  unlink(out, recursive = TRUE, force = TRUE)

  expect_error(
    cellucid_prepare(
      latent_space = latent,
      obs = obs,
      X_umap_2d = umap2,
      out_dir = out,
      centroid_min_points = 1,
      force = TRUE
    ),
    "obs_keys contains names that collide after filename sanitization"
  )
})
