test_that("only the canonical preparation API is exported", {
  expect_true("cellucid_prepare" %in% getNamespaceExports("cellucid"))
  expect_false("prepare" %in% getNamespaceExports("cellucid"))
  expect_false(
    exists("prepare", envir = asNamespace("cellucid"), inherits = FALSE)
  )
})

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
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    var = var,
    gene_expression = expr,
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint8"
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

test_that("identity obs fields follow emitted compact manifest order", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(
    group = factor(c("A", "A", "B")),
    score = c(0.1, 0.2, 0.3),
    selected = c(TRUE, FALSE, TRUE),
    quality = c(1, 2, 3)
  )
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  out <- file.path(tempdir(), "cellucid_test_identity_obs_order")
  unlink(out, recursive = TRUE, force = TRUE)

  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  identity <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  manifest_order <- c(
    lapply(
      manifest[["_continuousFields"]],
      function(field) list(key = field[[1]], kind = "continuous")
    ),
    lapply(
      manifest[["_categoricalFields"]],
      function(field) {
        list(
          key = field[[1]],
          kind = "category",
          n_categories = length(field[[2]])
        )
      }
    )
  )

  expect_identical(identity$obs_fields, manifest_order)
  expect_identical(
    vapply(identity$obs_fields, `[[`, character(1), "key"),
    c("score", "quality", "group", "selected")
  )
  expect_identical(identity$stats$n_obs_fields, 4L)
  expect_identical(identity$stats$n_continuous_fields, 2L)
  expect_identical(identity$stats$n_categorical_fields, 2L)
})

test_that("force FALSE refuses to mix with an existing generation", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  original_obs <- data.frame(
    group = factor(c("A", "A", "B")),
    score = c(0.1, 0.2, 0.3)
  )
  changed_obs <- data.frame(
    group = factor(c("A", "A", "B")),
    different_score = c(0.1, 0.2, 0.3)
  )
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  out <- file.path(tempdir(), "cellucid_test_existing_obs_contract")
  unlink(out, recursive = TRUE, force = TRUE)

  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = original_obs,
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  expect_error(
    cellucid_prepare(
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = changed_obs,
      X_umap_2d = umap2,
      out_dir = out,
      centroid_min_points = 1,
      force = FALSE,
      obs_categorical_dtype = "uint16"
    ),
    "out_dir already exists.*force = TRUE"
  )
})

test_that("force TRUE publishes one complete replacement generation", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap1 <- matrix(c(0, 1), ncol = 1)
  out <- tempfile("cellucid_r_atomic_replacement_")

  cellucid_prepare(
    dataset_id = "first-generation",
    dataset_name = "First generation",
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )
  writeLines("stale", file.path(out, "stale-generation-file"))

  cellucid_prepare(
    dataset_id = "second-generation",
    dataset_name = "Second generation",
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  identity <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  expect_identical(identity$id, "second-generation")
  expect_false(file.exists(file.path(out, "stale-generation-file")))
  stage_prefix <- paste0(".", basename(out), ".stage-")
  previous_prefix <- paste0(".", basename(out), ".previous-")
  siblings <- list.files(dirname(out), all.files = TRUE)
  expect_false(any(startsWith(siblings, stage_prefix)))
  expect_false(any(startsWith(siblings, previous_prefix)))
})

test_that("manifest field keys must already be portable filenames", {
  expect_error(
    cellucid:::.safe_filename_component(".."),
    "portable identifier"
  )
  expect_error(
    cellucid:::.safe_filename_component(1),
    "exactly one string"
  )
  expect_error(
    cellucid:::.safe_filename_component("gene."),
    "portable identifier"
  )
  for (value in c("CON", "con.txt", "AUX", "LPT1", "com9.csv")) {
    expect_error(
      cellucid:::.safe_filename_component(value),
      "reserved on Windows"
    )
  }
  expect_error(
    cellucid:::.assert_unique_filename_components(
      c("Gene", "gene"),
      "Gene identifiers"
    ),
    "case-insensitive filesystems"
  )
})

test_that("categorical storage uses the exact portable code capacities", {
  compact <- factor(sprintf("category_%03d", seq_len(255L)))
  compact_summary <- cellucid:::.summarize_obs_field(
    compact,
    key = "compact",
    obs_continuous_quantization = NULL,
    obs_categorical_dtype = "uint8"
  )
  expect_identical(compact_summary$codes_dtype, "uint8")

  wide <- factor(sprintf("category_%03d", seq_len(256L)))
  expect_error(
    cellucid:::.summarize_obs_field(
      wide,
      key = "wide",
      obs_continuous_quantization = NULL,
      obs_categorical_dtype = "uint8"
    ),
    "uint8 supports at most 255"
  )
  wide_summary <- cellucid:::.summarize_obs_field(
    wide,
    key = "wide",
    obs_continuous_quantization = NULL,
    obs_categorical_dtype = "uint16"
  )
  expect_identical(wide_summary$codes_dtype, "uint16")

  oversized <- factor(sprintf("category_%05d", seq_len(65536L)))
  expect_error(
    cellucid:::.summarize_obs_field(
      oversized,
      key = "oversized",
      obs_continuous_quantization = NULL,
      obs_categorical_dtype = "uint16"
    ),
    "uint16 supports at most 65,535"
  )
})

test_that("logical categories retain boolean JSON identity", {
  payload <- cellucid:::.category_values_and_codes(c(FALSE, TRUE, NA, FALSE))
  expect_identical(payload$categories, c(FALSE, TRUE))
  expect_identical(payload$codes, c(0L, 1L, -1L, 0L))

  centroids <- cellucid:::.compute_centroids_for_field(
    coords = matrix(c(0, 2), ncol = 1),
    codes = c(0L, 0L),
    categories = FALSE,
    outlier_quantile = 0.95,
    min_points = 1
  )
  expect_identical(centroids[[1]]$category, FALSE)
})

test_that("compression and quantization settings require exact scalar values", {
  expect_null(cellucid:::.normalize_compression(NULL))
  expect_identical(cellucid:::.normalize_compression(6), 6)
  expect_identical(
    cellucid:::.validate_quantization_bits(8, "bits"),
    8
  )

  for (value in list(0, 6.5, "6", TRUE, NA_real_, c(5, 6))) {
    expect_error(
      cellucid:::.normalize_compression(value),
      "NULL or one integer from 1 to 9"
    )
  }
  for (value in list(0, 8.5, "8", TRUE, NA_real_, c(8, 16))) {
    expect_error(
      cellucid:::.validate_quantization_bits(value, "bits"),
      "exactly 8 or 16"
    )
  }
})

test_that("public writer options reject partial matches and non-boolean force values", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap2 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)

  for (value in list("a", "u", "uint", "", NA_character_, c("auto", "uint8"))) {
    out <- tempfile("cellucid_r_exact_category_option_")
    expect_error(
      cellucid_prepare(
        dataset_id = "test-dataset",
        dataset_name = "Test dataset",
        latent_space = latent,
        obs = obs,
        X_umap_2d = umap2,
        out_dir = out,
        obs_categorical_dtype = value,
        force = TRUE
      ),
      "obs_categorical_dtype must be exactly one of"
    )
    expect_false(dir.exists(out))
  }

  for (value in list("TRUE", 1, 0, NA, c(TRUE, FALSE))) {
    out <- tempfile("cellucid_r_exact_force_")
    expect_error(
      cellucid_prepare(
        dataset_id = "test-dataset",
        dataset_name = "Test dataset",
        latent_space = latent,
        obs = obs,
        X_umap_2d = umap2,
        out_dir = out,
        force = value,
        obs_categorical_dtype = "uint16"
      ),
      "force must be exactly TRUE or FALSE"
    )
    expect_false(dir.exists(out))
  }
})

test_that("dataset and source identity require exact caller-owned strings", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap2 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)

  invalid_identity <- list(
    list(dataset_id = NULL, dataset_name = "Dataset"),
    list(dataset_id = "dataset", dataset_name = NULL),
    list(dataset_id = 1, dataset_name = "Dataset"),
    list(dataset_id = "dataset", dataset_name = 1),
    list(dataset_id = "", dataset_name = "Dataset"),
    list(dataset_id = "data/set", dataset_name = "Dataset"),
    list(dataset_id = "dataset", dataset_name = ""),
    list(dataset_id = "dataset", dataset_name = "   "),
    list(dataset_id = "dataset", dataset_name = " Dataset"),
    list(dataset_id = "dataset", dataset_name = "Dataset ")
  )
  for (identity in invalid_identity) {
    out <- tempfile("cellucid_r_invalid_identity_")
    expect_error(
      cellucid_prepare(
        dataset_id = identity$dataset_id,
        dataset_name = identity$dataset_name,
        latent_space = latent,
        obs = obs,
        X_umap_2d = umap2,
        out_dir = out,
        force = TRUE,
        obs_categorical_dtype = "uint16"
      ),
      "dataset_id|dataset_name"
    )
    expect_false(dir.exists(out))
  }

  partial_sources <- list(
    list(source_name = "Source"),
    list(source_url = "https://example.test/data"),
    list(source_citation = "Citation"),
    list(
      source_name = "Source",
      source_url = "https://example.test/data"
    )
  )
  for (source in partial_sources) {
    out <- tempfile("cellucid_r_partial_source_")
    expect_error(
      do.call(
        cellucid_prepare,
        c(
          list(
            dataset_id = "test-dataset",
            dataset_name = "Test dataset",
            latent_space = latent,
            obs = obs,
            X_umap_2d = umap2,
            out_dir = out,
            obs_categorical_dtype = "uint16",
            force = TRUE
          ),
          source
        )
      ),
      "source_name, source_url, and source_citation"
    )
    expect_false(dir.exists(out))
  }
})

test_that("the public writer exposes only the current fixed artifact layout", {
  removed <- c(
    "X_umap_4d",
    "obs_manifest_filename",
    "obs_binary_dirname",
    "var_manifest_filename",
    "var_binary_dirname",
    "connectivity_manifest_filename",
    "connectivity_binary_dirname"
  )
  expect_false(any(removed %in% names(formals(cellucid_prepare))))
})

test_that("uint8 serialization rejects wrapping and truncation", {
  invalid <- list(
    -1,
    256,
    1.5,
    NaN,
    Inf,
    "1",
    TRUE,
    c(0, NA_real_)
  )
  for (value in invalid) {
    path <- tempfile("cellucid_r_invalid_uint8_")
    expect_error(
      cellucid:::.write_uint8(path, value),
      "uint8 values"
    )
    expect_false(file.exists(path))
  }
})

test_that("unsupported observation label types are rejected without stringification", {
  expect_error(
    cellucid:::.category_values_and_codes(as.Date("2026-01-01")),
    "factor, logical, or character"
  )
  expect_error(
    cellucid:::.category_values_and_codes(I(list("A", "B"))),
    "factor, logical, or character"
  )
})

test_that("embedding normalization rejects undefined and non-finite domains", {
  expect_error(
    cellucid:::.normalize_embedding(matrix(c(2, 2, 2), ncol = 1)),
    "positive finite coordinate range"
  )
  expect_error(
    cellucid:::.normalize_embedding(matrix(c(0, Inf), ncol = 1)),
    "finite"
  )
})

test_that("matrix and observation inputs reject coercive representations", {
  expect_error(
    cellucid:::.as_dense_matrix(
      matrix(c("1", "2"), ncol = 1),
      "latent_space"
    ),
    "finite real numeric matrix"
  )
  expect_error(
    cellucid:::.as_matrix_like(
      matrix(c("1", "2"), ncol = 1),
      "gene_expression"
    ),
    "real numeric matrix"
  )
  expect_error(
    cellucid:::.collect_embeddings(
      list(`1` = c("0", "1"), `2` = NULL, `3` = NULL)
    ),
    "finite real numeric matrix"
  )
  expect_error(
    cellucid:::.observation_kind(
      as.Date("2026-01-01"),
      "date"
    ),
    "native numeric, logical, character, or factor"
  )
})

test_that("centroid_min_points requires one positive integer", {
  for (value in list(0, -1, 1.5, "1", TRUE, NA_real_, c(1, 2))) {
    expect_error(
      cellucid:::.validate_positive_integer(
        value,
        "centroid_min_points"
      ),
      "one positive integer"
    )
  }
  expect_identical(
    cellucid:::.validate_positive_integer(
      1,
      "centroid_min_points"
    ),
    1L
  )
})

test_that("singleton scientific arrays remain JSON arrays", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "A")))
  umap1 <- matrix(c(0, 1), ncol = 1)
  out <- tempfile("cellucid_r_singleton_arrays_")

  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  identity <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )

  expect_identical(identity$embeddings$available_dimensions, list(1L))
  expect_identical(
    manifest[["_categoricalFields"]][[1]][[2]],
    list("A")
  )
  expect_identical(
    manifest[["_categoricalFields"]][[1]][[5]][["1"]][[1]]$position,
    list(0L)
  )
})

test_that("invalid centroid and quantization domains are rejected", {
  expect_error(
    cellucid:::.compute_centroids_for_field(
      coords = matrix(c(0, 1), ncol = 1),
      codes = c(0L, 0L),
      categories = "A",
      outlier_quantile = 0.5,
      min_points = 1
    ),
    "strictly between 0.5 and 1"
  )
  expect_error(
    cellucid:::.quantize_continuous(
      numeric(0),
      bits = 8,
      field_name = "all-invalid"
    ),
    "contains no finite values"
  )
  expect_error(
    cellucid:::.quantize_continuous(
      c(3, 3),
      bits = 8,
      field_name = "constant"
    ),
    "all values are constant"
  )
  expect_error(
    cellucid:::.quantize_continuous(
      c(0, NA_real_, 1),
      bits = 8,
      field_name = "missing"
    ),
    "only finite values"
  )
  expect_error(
    cellucid:::.quantize_continuous(
      c(0, Inf),
      bits = 8,
      field_name = "infinite"
    ),
    "only finite values"
  )
})

test_that("only generated outlier NaN uses the exact reserved marker", {
  q <- cellucid:::.quantize_nullable_outlier_quantiles(
    c(2, NaN, 4),
    bits = 8,
    field_name = "cluster_outliers"
  )

  expect_identical(q$quantized, c(0L, 255L, 254L))
  expect_identical(q$min_val, 2)
  expect_identical(q$max_val, 4)
  expect_identical(q$scale, 127)

  expect_error(
    cellucid:::.quantize_nullable_outlier_quantiles(
      c(0, NA_real_, 1),
      bits = 8,
      field_name = "cluster_outliers"
    ),
    "only generated NaN"
  )
  expect_error(
    cellucid:::.quantize_nullable_outlier_quantiles(
      c(0, Inf, 1),
      bits = 8,
      field_name = "cluster_outliers"
    ),
    "finite values or generated NaN"
  )
})

test_that("manifest identifiers require exact character values", {
  expect_error(
    cellucid:::.validate_character_vector(1:2, "obs_keys"),
    "character vector"
  )
  expect_error(
    cellucid:::.validate_gene_ids(c(1, 2)),
    "character vector"
  )

  var <- data.frame(numeric_id = c(101, 102))
  expect_error(
    cellucid:::.extract_gene_ids(var, "numeric_id"),
    "character vector"
  )
  expect_error(
    cellucid:::.extract_gene_ids(var, 1),
    "NULL or one non-empty string"
  )
})

test_that("points_2d.bin is float32 row-major and normalized", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))
  umap2 <- matrix(c(0, 0, 2, 0, 0, 1), ncol = 2, byrow = TRUE)

  out <- file.path(tempdir(), "cellucid_test_points")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
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

test_that("continuous observations reject missing values for every codec", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(score = c(0, 1, NA_real_))
  umap1 <- matrix(c(0, 1, 2), ncol = 1)

  for (quantization in list(NULL, 8L)) {
    out <- tempfile("cellucid_r_nonfinite_obs_")

    expect_error(
      cellucid_prepare(
        dataset_id = "test-dataset",
        dataset_name = "Test dataset",
        latent_space = latent,
        obs = obs,
        X_umap_1d = umap1,
        out_dir = out,
        obs_continuous_quantization = quantization,
        centroid_min_points = 1,
        force = TRUE,
        obs_categorical_dtype = "uint16"
      ),
      "score.*only finite values"
    )
    expect_false(dir.exists(out))
  }
})

test_that("categorical outlier quantiles preserve generated missing markers", {
  latent <- matrix(
    c(0, 0, 1, 0, 4, 0, 8, 0),
    ncol = 2,
    byrow = TRUE
  )
  obs <- data.frame(
    cluster = factor(c("large", "large", "large", "small"))
  )
  umap2 <- matrix(
    c(0, 0, 1, 0, 0, 1, 1, 1),
    ncol = 2,
    byrow = TRUE
  )
  out <- tempfile("cellucid_r_outlier_missing_")

  cellucid_prepare(
    dataset_id = "outlier-missing-marker",
    dataset_name = "Outlier missing marker",
    latent_space = latent,
    obs = obs,
    X_umap_2d = umap2,
    out_dir = out,
    obs_continuous_quantization = 8,
    centroid_min_points = 2,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  con <- file(
    file.path(out, "obs", "cluster.outliers.u8"),
    open = "rb"
  )
  on.exit(close(con), add = TRUE)
  got <- readBin(con, what = "raw", n = 4)
  expect_identical(as.integer(got), c(127L, 0L, 254L, 255L))
})

test_that("connectivity export writes uint16 edge pairs", {

  latent <- matrix(c(0, 0, 1, 1, 2, 2, 3, 3), ncol = 2, byrow = TRUE)
  obs <- data.frame(a = 1:4)
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)

  conn <- matrix(0, nrow = 4, ncol = 4)
  conn[1, 2] <- 1
  conn[2, 1] <- 1
  conn[1, 3] <- 1
  conn[3, 1] <- 1
  conn[4, 3] <- 1
  conn[3, 4] <- 1

  out <- file.path(tempdir(), "cellucid_test_connectivity")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    connectivities = conn,
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  src_path <- file.path(out, "connectivity", "edges.src.bin")
  dst_path <- file.path(out, "connectivity", "edges.dst.bin")
  weights_path <- file.path(
    out,
    "connectivity",
    "edges.weights.f64.bin"
  )
  expect_true(file.exists(src_path))
  expect_true(file.exists(dst_path))
  expect_true(file.exists(weights_path))

  read_u16 <- function(path) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    readBin(con, what = "integer", size = 2, endian = "little", n = 1000)
  }

  # Expected unique undirected edges: (0,1), (0,2), (2,3)
  expect_equal(read_u16(src_path)[1:3], c(0, 0, 2))
  expect_equal(read_u16(dst_path)[1:3], c(1, 2, 3))
  weights_connection <- file(weights_path, open = "rb")
  on.exit(close(weights_connection), add = TRUE)
  expect_identical(
    readBin(
      weights_connection,
      what = "double",
      size = 8L,
      endian = "little",
      n = 3L
    ),
    rep(1, 3)
  )
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
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    X_umap_2d = umap2,
    vector_fields = vector_fields,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
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

test_that("vector fields use one exact finite UMAP naming and ownership contract", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap2 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  vector2 <- matrix(c(1, 0, 1, 0), ncol = 2, byrow = TRUE)

  invalid_vectors <- list(
    list(
      value = list(velocity = vector2),
      message = "must use '<field>_umap'"
    ),
    list(
      value = list(velocity_umap = NULL),
      message = "cannot be NULL"
    ),
    list(
      value = list(
        velocity_umap = vector2,
        velocity_umap_2d = vector2
      ),
      message = "same field and dimension"
    ),
    list(
      value = list(
        velocity_umap_3d = matrix(
          c(1, 0, 0, 1, 0, 0),
          ncol = 3,
          byrow = TRUE
        )
      ),
      message = "matching 3D embedding"
    ),
    list(
      value = list(
        velocity_umap_2d = matrix(
          c(1, 0, Inf, 0),
          ncol = 2,
          byrow = TRUE
        )
      ),
      message = "finite"
    )
  )

  for (case in invalid_vectors) {
    out <- tempfile("cellucid_r_invalid_vector_")
    expect_error(
      cellucid_prepare(
        dataset_id = "test-dataset",
        dataset_name = "Test dataset",
        latent_space = latent,
        obs = obs,
        X_umap_2d = umap2,
        vector_fields = case$value,
        out_dir = out,
        centroid_min_points = 1,
        force = TRUE,
        obs_categorical_dtype = "uint16"
      ),
      case$message
    )
    expect_false(dir.exists(out))
  }
})

test_that("multiple vector fields require one explicit exact default", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap2 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  vector2 <- matrix(c(1, 0, 1, 0), ncol = 2, byrow = TRUE)
  fields <- list(
    velocity_umap = vector2,
    displacement_umap = vector2
  )

  base_args <- list(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    X_umap_2d = umap2,
    vector_fields = fields,
    centroid_min_points = 1,
    obs_categorical_dtype = "uint16",
    force = TRUE
  )
  out <- tempfile("cellucid_r_vector_default_missing_")
  expect_error(
    do.call(cellucid_prepare, c(base_args, list(out_dir = out))),
    "vector_field_default is required"
  )
  expect_false(dir.exists(out))

  out <- tempfile("cellucid_r_vector_default_unknown_")
  expect_error(
    do.call(
      cellucid_prepare,
      c(
        base_args,
        list(
          out_dir = out,
          vector_field_default = "unknown_umap"
        )
      )
    ),
    "vector_field_default.*not present"
  )
  expect_false(dir.exists(out))

  out <- tempfile("cellucid_r_vector_default_exact_")
  do.call(
    cellucid_prepare,
    c(
      base_args,
      list(
        out_dir = out,
        vector_field_default = "displacement_umap"
      )
    )
  )
  identity <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  expect_identical(
    identity$vector_fields$default_field,
    "displacement_umap"
  )
  expect_identical(
    identity$vector_fields$fields$velocity_umap$label,
    "velocity_umap"
  )
})

test_that("obs manifest preserves min/max precision in JSON", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(score = c(1 / 3, 2 / 3, 1 / 2))
  umap1 <- matrix(c(0, 1, 2), ncol = 1)

  out <- file.path(tempdir(), "cellucid_test_json_precision")
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = out,
    obs_continuous_quantization = 8,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  manifest <- jsonlite::read_json(file.path(out, "obs_manifest.json"), simplifyVector = FALSE)
  fields <- manifest[["_continuousFields"]]
  expect_equal(length(fields), 1L)

  entry <- fields[[1]]
  expect_equal(entry[[1]], "score")
  expect_equal(entry[[2]], 1 / 3, tolerance = 1e-15)
  expect_equal(entry[[3]], 2 / 3, tolerance = 1e-15)
})

test_that("gene expression rejects missing values for every codec", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  expression <- matrix(c(0, NA_real_, 2), ncol = 1)
  var <- data.frame(row.names = "GeneA")

  for (quantization in list(NULL, 8L)) {
    out <- tempfile("cellucid_r_nonfinite_gene_")

    expect_error(
      cellucid_prepare(
        dataset_id = "test-dataset",
        dataset_name = "Test dataset",
        latent_space = latent,
        obs = obs,
        var = var,
        gene_expression = expression,
        X_umap_2d = umap2,
        out_dir = out,
        var_quantization = quantization,
        centroid_min_points = 1,
        force = TRUE,
        obs_categorical_dtype = "uint16"
      ),
      "gene_expression.*only finite values"
    )
    expect_false(dir.exists(out))
  }
})

test_that("cellucid_prepare errors on non-finite embedding values", {
  latent <- matrix(c(0, 0), ncol = 2)
  obs <- data.frame(cluster = factor("A"))
  umap2 <- matrix(c(0, NA), ncol = 2)

  expect_error(
    cellucid_prepare(
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = obs,
      X_umap_2d = umap2,
      out_dir = file.path(tempdir(), "cellucid_test_bad_embedding"),
      centroid_min_points = 1,
      force = TRUE,
      obs_categorical_dtype = "uint16"
    ),
    "non-finite values"
  )
})

test_that("cellucid_prepare errors without embeddings", {
  latent <- matrix(c(0, 0), ncol = 2)
  obs <- data.frame(a = 1)
  expect_error(
    cellucid_prepare(
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = obs,
      X_umap_1d = NULL,
      X_umap_2d = NULL,
      X_umap_3d = NULL,
      out_dir = tempfile("cellucid_r_missing_embedding_"),
      force = TRUE,
      obs_categorical_dtype = "uint16"
    ),
    "At least one embedding"
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
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = obs,
      var = var,
      gene_expression = expr,
      var_gene_id_column = "symbol",
      X_umap_2d = umap2,
      out_dir = out,
      centroid_min_points = 1,
      force = TRUE,
      obs_categorical_dtype = "uint16"
    ),
    "Gene identifiers must be unique"
  )
})

test_that("gene ids must already be portable filename components", {
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
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = obs,
      var = var,
      gene_expression = expr,
      X_umap_2d = umap2,
      out_dir = out,
      centroid_min_points = 1,
      force = TRUE,
      obs_categorical_dtype = "uint16"
    ),
    "portable identifier"
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
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    var = var,
    gene_expression = expr,
    gene_identifiers = c("G2"),
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  var_dir <- file.path(out, "var")
  expect_false(file.exists(file.path(var_dir, "G1.values.f32")))
  expect_true(file.exists(file.path(var_dir, "G2.values.f32")))

  ident <- jsonlite::read_json(file.path(out, "dataset_identity.json"), simplifyVector = TRUE)
  expect_equal(ident$stats$n_genes, 1L)
})

test_that("gene_identifiers rejects missing and non-character selections", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  expr <- matrix(c(0, 1, 2, 3, 4, 5), nrow = 3, ncol = 2)
  var <- data.frame(symbol = c("G1", "G2"))
  rownames(var) <- var$symbol

  expect_error(
    cellucid_prepare(
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = obs,
      var = var,
      gene_expression = expr,
      gene_identifiers = c("G2", "missing"),
      X_umap_2d = umap2,
      out_dir = file.path(tempdir(), "cellucid_test_missing_gene_subset"),
      centroid_min_points = 1,
      force = TRUE,
      obs_categorical_dtype = "uint16"
    ),
    "identifiers not found in var: missing"
  )

  expect_error(
    cellucid_prepare(
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = obs,
      var = var,
      gene_expression = expr,
      gene_identifiers = 2,
      X_umap_2d = umap2,
      out_dir = file.path(tempdir(), "cellucid_test_numeric_gene_subset"),
      centroid_min_points = 1,
      force = TRUE,
      obs_categorical_dtype = "uint16"
    ),
    "gene_identifiers must be a character vector"
  )

  expect_error(
    cellucid_prepare(
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = obs,
      var = var,
      gene_expression = expr,
      gene_identifiers = c("G2", "G2"),
      X_umap_2d = umap2,
      out_dir = file.path(tempdir(), "cellucid_test_duplicate_gene_subset"),
      centroid_min_points = 1,
      force = TRUE,
      obs_categorical_dtype = "uint16"
    ),
    "must not contain duplicate identifiers"
  )
})

test_that("obs keys must already be portable filename components", {
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
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = obs,
      X_umap_2d = umap2,
      out_dir = out,
      centroid_min_points = 1,
      force = TRUE,
      obs_categorical_dtype = "uint16"
    ),
    "portable identifier"
  )

  expect_error(
    cellucid_prepare(
      dataset_id = "test-dataset",
      dataset_name = "Test dataset",
      latent_space = latent,
      obs = obs,
      obs_keys = 1,
      X_umap_2d = umap2,
      out_dir = file.path(tempdir(), "cellucid_test_numeric_obs_key"),
      centroid_min_points = 1,
      force = TRUE,
      obs_categorical_dtype = "uint16"
    ),
    "obs_keys must be a character vector"
  )
})

test_that("float32 writer rejects overflow and underflow before publication", {
  float32_max <- (2 - 2^-23) * 2^127
  float32_min_subnormal <- 2^-149
  path <- tempfile("cellucid_r_float32_boundary_")

  cellucid:::.write_float32_vector(
    path,
    c(
      -float32_max,
      -float32_min_subnormal,
      0,
      float32_min_subnormal,
      float32_max
    )
  )
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  encoded <- readBin(
    connection,
    what = "numeric",
    size = 4L,
    endian = "little",
    n = 5L
  )
  expect_true(all(is.finite(encoded)))
  expect_true(all(encoded[c(1L, 2L)] < 0))
  expect_identical(encoded[[3L]], 0)
  expect_true(all(encoded[c(4L, 5L)] > 0))

  expect_error(
    cellucid:::.write_float32_vector(tempfile(), 2^128),
    "finite float32"
  )
  expect_error(
    cellucid:::.write_float32_vector(tempfile(), 2^-150),
    "finite float32"
  )
  expect_error(
    cellucid:::.write_float32_vector(tempfile(), NaN),
    "generated NaN"
  )

  nullable_path <- tempfile("cellucid_r_float32_nullable_")
  cellucid:::.write_float32_vector(
    nullable_path,
    c(0, NaN),
    allow_generated_nan = TRUE
  )
  nullable_connection <- file(nullable_path, open = "rb")
  on.exit(close(nullable_connection), add = TRUE)
  nullable <- readBin(
    nullable_connection,
    what = "numeric",
    size = 4L,
    endian = "little",
    n = 2L
  )
  expect_identical(nullable[[1L]], 0)
  expect_true(is.nan(nullable[[2L]]))
})

test_that("all public float32 artifacts reject unrepresentable finite values", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  umap2 <- latent
  categorical_obs <- data.frame(group = factor(c("A", "B")))
  base_args <- list(
    dataset_id = "float32-contract",
    dataset_name = "Float32 contract",
    latent_space = latent,
    X_umap_2d = umap2,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )
  cases <- list(
    continuous_observation = list(
      obs = data.frame(score = c(0, 1e300))
    ),
    gene_expression = list(
      obs = categorical_obs,
      var = data.frame(
        symbol = "GENE",
        row.names = "GENE"
      ),
      gene_expression = matrix(c(0, 1e300), ncol = 1)
    ),
    vector_field = list(
      obs = categorical_obs,
      vector_fields = list(
        velocity_umap_2d = matrix(
          c(0, 0, 1e300, 0),
          ncol = 2,
          byrow = TRUE
        )
      )
    )
  )

  for (case_name in names(cases)) {
    out <- tempfile(paste0("cellucid_r_", case_name, "_"))
    args <- c(
      base_args,
      cases[[case_name]],
      list(out_dir = out)
    )
    expect_error(
      do.call(cellucid_prepare, args),
      "finite float32",
      info = case_name
    )
    expect_false(dir.exists(out), info = case_name)
  }
})
