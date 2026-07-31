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

  # Every payload path is an index, so the directory listing is the same on
  # every dataset and carries none of this dataset's vocabulary.
  obs_dir <- file.path(out, "obs")
  expect_identical(
    sort(list.files(obs_dir)),
    c("0.codes.u8", "0.outliers.f32", "1.values.f32")
  )

  var_dir <- file.path(out, "var")
  expect_identical(
    sort(list.files(var_dir)),
    c("0.values.f32", "1.values.f32")
  )

  obs_manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(
    obs_manifest[["_obsSchemas"]]$continuous$pathPattern,
    "obs/{index}.values.f32"
  )
  expect_identical(
    obs_manifest[["_obsSchemas"]]$categorical$codesPathPattern,
    "obs/{index}.codes.{ext}"
  )
  expect_identical(
    obs_manifest[["_obsSchemas"]]$categorical$outlierPathPattern,
    "obs/{index}.outliers.f32"
  )
  expect_identical(obs_manifest[["_continuousFields"]][[1L]], list(1L, "score"))
  expect_identical(obs_manifest[["_categoricalFields"]][[1L]][[1L]], 0L)
  expect_identical(obs_manifest[["_categoricalFields"]][[1L]][[2L]], "cluster")

  var_manifest <- jsonlite::read_json(
    file.path(out, "var_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(
    var_manifest[["_varSchema"]]$pathPattern,
    "var/{index}.values.f32"
  )
  expect_identical(
    var_manifest$fields,
    list(list(0L, "G1"), list(1L, "G2"))
  )
})

# Identifiers that no filesystem could carry, and that no longer have to. The
# cellucid Python package holds the identical set in its own contract test.
.hostile_obs_keys <- c("% mito", "cell type", "CON", "Field", "field")
.hostile_gene_names <- c(
  "HLA-DRB1/2", "NUL.txt", "trailing.", "細胞", "Gene", "gene"
)

.hostile_export <- function(out, ...) {
  n_cells <- 6L
  coordinates <- matrix(
    c(0, 0, 1, 0, 0, 1, 1, 1, 0.5, 0.25, 0.25, 0.75),
    ncol = 2,
    byrow = TRUE
  )
  # Distinct within-category distances, so the generated outlier quantiles vary
  # and the quantized variant of this export is encodable.
  latent <- matrix(
    c(0, 0, 3, 0, 0.5, 0, 7, 0, 1.5, 0, 9, 0),
    ncol = 2,
    byrow = TRUE
  )
  connectivities <- matrix(0, nrow = n_cells, ncol = n_cells)
  for (edge in list(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 6))) {
    connectivities[edge[[1]], edge[[2]]] <- 0.5
    connectivities[edge[[2]], edge[[1]]] <- 0.5
  }

  obs <- data.frame(
    `% mito` = seq(0, 1, length.out = n_cells),
    `cell type` = factor(c("a", "a", "a", "b", "b", "b")),
    CON = factor(c("x", "y", "y", "y", "x", "x")),
    Field = seq(1, 2, length.out = n_cells),
    field = seq(3, 4, length.out = n_cells),
    check.names = FALSE
  )
  expect_identical(names(obs), .hostile_obs_keys)

  var <- data.frame(row.names = .hostile_gene_names)
  expression <- matrix(
    as.numeric(seq_len(n_cells * length(.hostile_gene_names))),
    nrow = n_cells,
    ncol = length(.hostile_gene_names)
  )

  arguments <- list(
    dataset_id = "payload-index",
    dataset_name = "Payload index",
    created_at = "2026-01-01T00:00:00Z",
    latent_space = latent,
    obs = obs,
    var = var,
    gene_expression = expression,
    connectivities = connectivities,
    X_umap_2d = coordinates,
    vector_fields = list(
      Velocity_umap_2d = coordinates * 0.25,
      velocity_umap_2d = coordinates * 0.5
    ),
    vector_field_default = "Velocity_umap",
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )
  do.call(cellucid_prepare, utils::modifyList(arguments, list(...)))
  out
}

test_that("identifiers that can never be filenames export without complaint", {
  out <- .hostile_export(tempfile("cellucid_r_hostile_identifiers_"))

  obs_manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  var_manifest <- jsonlite::read_json(
    file.path(out, "var_manifest.json"),
    simplifyVector = FALSE
  )
  exported_obs_keys <- c(
    vapply(obs_manifest[["_continuousFields"]], `[[`, character(1), 2L),
    vapply(obs_manifest[["_categoricalFields"]], `[[`, character(1), 2L)
  )
  expect_setequal(exported_obs_keys, .hostile_obs_keys)
  expect_identical(
    vapply(var_manifest$fields, `[[`, character(1), 2L),
    .hostile_gene_names
  )
})

test_that("no payload filename carries any identifier", {
  out <- .hostile_export(tempfile("cellucid_r_neutral_paths_"))

  indexed_payload <- "^(0|[1-9][0-9]*)\\.(values|codes|outliers)\\.(f32|u8|u16)$"
  for (directory in c("obs", "var")) {
    names <- sort(list.files(file.path(out, directory)))
    expect_gt(length(names), 0L)
    expect_true(all(grepl(indexed_payload, names)), info = directory)
  }
  expect_true(
    all(grepl(
      "^(0|[1-9][0-9]*)_[123]d\\.bin$",
      list.files(file.path(out, "vectors"))
    ))
  )
  expect_identical(
    sort(list.files(file.path(out, "connectivity"))),
    c("edges.dst.bin", "edges.src.bin", "edges.weights.f64.bin")
  )

  # Nothing anywhere in the tree spells any identifier, including the sanitized
  # forms an escaping scheme would have produced.
  written <- list.files(out, recursive = TRUE, all.files = TRUE)
  joined <- paste(written, collapse = "\n")
  identifiers <- c(
    .hostile_obs_keys,
    .hostile_gene_names,
    "Velocity_umap",
    "velocity_umap"
  )
  for (identifier in identifiers) {
    expect_false(grepl(identifier, joined, fixed = TRUE), info = identifier)
    expect_false(
      grepl(
        gsub("/", "_", gsub(" ", "_", identifier, fixed = TRUE), fixed = TRUE),
        joined,
        fixed = TRUE
      ),
      info = identifier
    )
  }
  expect_setequal(
    written,
    c(
      "points_2d.bin",
      "obs_manifest.json",
      "var_manifest.json",
      "connectivity_manifest.json",
      "dataset_identity.json",
      file.path("obs", c(
        "0.values.f32",
        "1.codes.u16", "1.outliers.f32",
        "2.codes.u16", "2.outliers.f32",
        "3.values.f32",
        "4.values.f32"
      )),
      file.path("var", sprintf("%d.values.f32", 0:5)),
      file.path(
        "connectivity",
        c("edges.src.bin", "edges.dst.bin", "edges.weights.f64.bin")
      ),
      file.path("vectors", c("0_2d.bin", "1_2d.bin"))
    )
  )
})

test_that("obs indices are one dense space shared by both field arrays", {
  out <- .hostile_export(tempfile("cellucid_r_shared_obs_index_"))
  manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  continuous <- manifest[["_continuousFields"]]
  categorical <- manifest[["_categoricalFields"]]

  # The space follows the obs column order, so it interleaves the two arrays.
  expect_identical(
    lapply(continuous, function(field) field[1:2]),
    list(list(0L, "% mito"), list(3L, "Field"), list(4L, "field"))
  )
  expect_identical(
    lapply(categorical, function(field) field[1:2]),
    list(list(1L, "cell type"), list(2L, "CON"))
  )
  indices <- c(
    vapply(continuous, `[[`, integer(1), 1L),
    vapply(categorical, `[[`, integer(1), 1L)
  )
  expect_identical(sort(indices), 0:4)
})

test_that("var and vector indices are dense and position based", {
  out <- .hostile_export(tempfile("cellucid_r_dense_var_vector_index_"))
  var_manifest <- jsonlite::read_json(
    file.path(out, "var_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(
    vapply(var_manifest$fields, `[[`, integer(1), 1L),
    seq_along(.hostile_gene_names) - 1L
  )

  identity <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  # Both writers emit vector fields in code-point order of their ids, so the
  # same input receives the same payload index in either language.
  expect_identical(
    lapply(identity$vector_fields$fields, function(field) field$files),
    list(
      Velocity_umap = list(`2d` = "vectors/0_2d.bin"),
      velocity_umap = list(`2d` = "vectors/1_2d.bin")
    )
  )
})

test_that("manifest path patterns substitute the index", {
  out <- .hostile_export(tempfile("cellucid_r_index_path_patterns_"))
  obs_manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  var_manifest <- jsonlite::read_json(
    file.path(out, "var_manifest.json"),
    simplifyVector = FALSE
  )

  expect_identical(
    obs_manifest[["_obsSchemas"]]$continuous$pathPattern,
    "obs/{index}.values.f32"
  )
  expect_identical(
    obs_manifest[["_obsSchemas"]]$categorical,
    list(
      codesPathPattern = "obs/{index}.codes.{ext}",
      outlierPathPattern = "obs/{index}.outliers.f32",
      outlierExt = "f32",
      outlierDtype = "float32",
      outlierQuantized = FALSE
    )
  )
  expect_identical(
    var_manifest[["_varSchema"]]$pathPattern,
    "var/{index}.values.f32"
  )
})

test_that("quantized entries keep the index and gain their bounds", {
  out <- .hostile_export(
    tempfile("cellucid_r_quantized_index_entries_"),
    var_quantization = 8,
    obs_continuous_quantization = 8,
    compression = 6
  )

  var_manifest <- jsonlite::read_json(
    file.path(out, "var_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(
    var_manifest[["_varSchema"]]$pathPattern,
    "var/{index}.values.u8.gz"
  )
  for (position in seq_along(var_manifest$fields)) {
    field <- var_manifest$fields[[position]]
    expect_length(field, 4L)
    expect_identical(field[[1L]], position - 1L)
    expect_identical(field[[2L]], .hostile_gene_names[[position]])
    expect_lt(field[[3L]], field[[4L]])
  }

  obs_manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(
    obs_manifest[["_obsSchemas"]]$continuous$pathPattern,
    "obs/{index}.values.u8.gz"
  )
  expect_identical(
    obs_manifest[["_obsSchemas"]]$categorical$outlierPathPattern,
    "obs/{index}.outliers.u8.gz"
  )
  for (field in obs_manifest[["_continuousFields"]]) {
    expect_length(field, 4L)
  }
  for (field in obs_manifest[["_categoricalFields"]]) {
    expect_length(field, 8L)
  }
  expect_true(all(endsWith(list.files(file.path(out, "obs")), ".gz")))
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
      function(field) list(key = field[[2]], kind = "continuous")
    ),
    lapply(
      manifest[["_categoricalFields"]],
      function(field) {
        list(
          key = field[[2]],
          kind = "category",
          n_categories = length(field[[3]])
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

  # obs/ is written by both arrays, so their declared indices are one shared
  # space: exactly the obs_keys positions, 0..3, each used once.
  expect_identical(
    vapply(manifest[["_continuousFields"]], `[[`, integer(1), 1L),
    c(1L, 3L)
  )
  expect_identical(
    vapply(manifest[["_categoricalFields"]], `[[`, integer(1), 1L),
    c(0L, 2L)
  )
  expect_identical(
    sort(list.files(file.path(out, "obs"))),
    c(
      "0.codes.u16", "0.outliers.f32",
      "1.values.f32",
      "2.codes.u16", "2.outliers.f32",
      "3.values.f32"
    )
  )
})

test_that("the writer refuses a payload index space it cannot write once", {
  expect_error(
    cellucid:::.require_dense_payload_indices(
      list(0L, 1L, 1L),
      axis = "Observation"
    ),
    paste0(
      "^Observation payload indices must be exactly 0\\.\\.2, each used once; ",
      "got \\[0, 1, 1\\]\\.$"
    )
  )
  expect_error(
    cellucid:::.require_dense_payload_indices(list(0L, 2L), axis = "Gene"),
    "^Gene payload indices must be exactly 0\\.\\.1, each used once; got \\[0, 2\\]\\.$"
  )
  expect_error(
    cellucid:::.require_dense_payload_indices(list(0L, 1), axis = "Gene"),
    "^Gene payload index at position 1 must be a native integer\\.$"
  )
  expect_identical(
    cellucid:::.require_dense_payload_indices(list(1L, 0L), axis = "Gene"),
    c(1L, 0L)
  )
  expect_identical(
    cellucid:::.require_dense_payload_indices(list(), axis = "Gene"),
    integer(0)
  )
})

test_that("an emitted manifest whose entries collide on one index is rejected", {
  colliding_obs <- list(
    `_format` = "compact_v1",
    n_points = 2L,
    centroid_outlier_quantile = NULL,
    latent_key = "latent_space",
    compression = NULL,
    `_obsSchemas` = list(),
    `_continuousFields` = list(list(0L, "score")),
    `_categoricalFields` = list(
      list(0L, "cluster", I("A"), "uint8", 255L, list())
    )
  )
  expect_error(
    cellucid:::.identity_obs_fields_from_compact_manifest(colliding_obs),
    "^Observation payload indices must be exactly 0\\.\\.1, each used once"
  )

  expect_error(
    cellucid:::.gene_names_from_compact_manifest(
      list(fields = list(list(0L, "TNMD"), list(0L, "CIITA")))
    ),
    "^Gene payload indices must be exactly 0\\.\\.1, each used once"
  )
  expect_error(
    cellucid:::.gene_names_from_compact_manifest(
      list(fields = list(list(0L, "TNMD"), list(1L, "TNMD")))
    ),
    "^Gene key 'TNMD' is duplicated\\.$"
  )
  expect_error(
    cellucid:::.gene_names_from_compact_manifest(
      list(fields = list(list("0", "TNMD")))
    ),
    "^Gene payload index at position 0 must be a native integer\\.$"
  )
  expect_identical(
    cellucid:::.gene_names_from_compact_manifest(
      list(fields = list(list(0L, "TNMD"), list(1L, "CIITA")))
    ),
    c("TNMD", "CIITA")
  )
})

test_that("a payload pattern must resolve to one complete path", {
  expect_identical(
    cellucid:::.expand_payload_pattern(
      "obs/{index}.codes.{ext}.gz",
      index = 3L,
      label = "obs categorical codesPathPattern",
      ext = "u16"
    ),
    "obs/3.codes.u16.gz"
  )
  expect_error(
    cellucid:::.expand_payload_pattern(
      "obs/{index}.codes.{ext}",
      index = 3L,
      label = "obs categorical codesPathPattern"
    ),
    paste0(
      "^obs categorical codesPathPattern retains an unsubstituted ",
      "placeholder: 'obs/3\\.codes\\.\\{ext\\}'\\.$"
    )
  )
  expect_error(
    cellucid:::.expand_payload_pattern(NULL, index = 0L, label = "var pathPattern"),
    "^var pathPattern must be a non-empty path pattern\\.$"
  )
})

test_that("an axis directory must hold exactly the payloads its manifest declares", {
  root <- tempfile("cellucid_r_declared_payloads_")
  dir.create(file.path(root, "var"), recursive = TRUE)
  writeBin(raw(1L), file.path(root, "var", "0.values.f32"))

  expect_identical(
    cellucid:::.require_declared_payloads_on_disk(
      root,
      directory_name = "var",
      declared = "var/0.values.f32",
      axis = "Gene"
    ),
    "var/0.values.f32"
  )
  expect_error(
    cellucid:::.require_declared_payloads_on_disk(
      root,
      directory_name = "var",
      declared = c("var/0.values.f32", "var/1.values.f32"),
      axis = "Gene"
    ),
    paste0(
      "^Gene manifest does not describe the payloads that were written\\. ",
      "Declared but absent: var/1\\.values\\.f32\\. Written but undeclared: \\.$"
    )
  )
  writeBin(raw(1L), file.path(root, "var", "stray.bin"))
  expect_error(
    cellucid:::.require_declared_payloads_on_disk(
      root,
      directory_name = "var",
      declared = "var/0.values.f32",
      axis = "Gene"
    ),
    "Written but undeclared: var/stray\\.bin\\.$"
  )
})

test_that("payload positions are one-based and become zero-based indices", {
  expect_identical(cellucid:::.payload_index(1L), 0L)
  expect_identical(cellucid:::.payload_index(4), 3L)
  for (invalid in list(0L, -1L, 1.5, NA_integer_, c(1L, 2L), "1", Inf)) {
    expect_error(
      cellucid:::.payload_index(invalid),
      "^Payload position must be one positive integer\\.$"
    )
  }
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
  transaction_prefix <- paste0(".", basename(out), ".cellucid-")
  siblings <- list.files(dirname(out), all.files = TRUE)
  expect_false(any(startsWith(siblings, transaction_prefix)))
})

test_that("a live independent exporter excludes replacement until process exit", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap1 <- matrix(c(0, 1), ncol = 1)
  out <- tempfile("cellucid_r_live_export_lock_")

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
  sentinel <- file.path(out, "first-generation-sentinel")
  writeLines("owned by the first generation", sentinel)
  lock_path <- file.path(
    dirname(out),
    paste0(".", basename(out), ".cellucid.lock")
  )

  cluster <- parallel::makeCluster(1L)
  on.exit({
    if (!is.null(cluster)) {
      parallel::stopCluster(cluster)
    }
  }, add = TRUE)
  held <- parallel::clusterCall(
    cluster,
    function(path) {
      library(cellucid)
      lock <- cellucid:::.native_export_lock_acquire(path)
      stopifnot(!is.null(lock))
      assign(".cellucid_test_export_lock", lock, envir = .GlobalEnv)
      TRUE
    },
    lock_path
  )
  expect_identical(held, list(TRUE))

  expect_error(
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
    ),
    "generation.*already active"
  )
  expect_true(file.exists(sentinel))

  parallel::stopCluster(cluster)
  cluster <- NULL
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
  expect_false(file.exists(sentinel))
  expect_true(file.exists(lock_path))
  expect_identical(unname(file.info(lock_path)$size), 0)
  siblings <- list.files(dirname(out), all.files = TRUE)
  expect_false(any(startsWith(
    siblings,
    paste0(".", basename(out), ".cellucid-")
  )))
})

test_that("rejected native lock attempts do not leak descriptors or handles", {
  parent <- tempfile("cellucid_r_native_lock_resources_")
  dir.create(parent)
  lock_path <- file.path(parent, ".target.cellucid.lock")
  cluster <- parallel::makeCluster(1L)
  on.exit({
    if (!is.null(cluster)) {
      parallel::stopCluster(cluster)
    }
  }, add = TRUE)
  held <- parallel::clusterCall(
    cluster,
    function(path) {
      library(cellucid)
      lock <- cellucid:::.native_export_lock_acquire(path)
      stopifnot(!is.null(lock))
      assign(".cellucid_test_export_lock", lock, envir = .GlobalEnv)
      TRUE
    },
    lock_path
  )
  expect_identical(held, list(TRUE))
  expect_null(cellucid:::.native_export_lock_acquire(lock_path))

  resource_count <- function() {
    if (.Platform$OS.type == "windows") {
      return(as.numeric(cellucid:::.native_process_handle_count()))
    }
    descriptor_root <- if (dir.exists("/proc/self/fd")) {
      "/proc/self/fd"
    } else if (dir.exists("/dev/fd")) {
      "/dev/fd"
    } else {
      return(NA_real_)
    }
    as.numeric(length(list.files(descriptor_root)))
  }

  before <- resource_count()
  skip_if(is.na(before), "native process resource count is unavailable")
  for (attempt in seq_len(1000L)) {
    stopifnot(is.null(cellucid:::.native_export_lock_acquire(lock_path)))
  }
  gc()
  after <- resource_count()
  expect_identical(after, before)
})

test_that("garbage collection finalizes an abandoned native lock handle", {
  parent <- tempfile("cellucid_r_native_lock_gc_")
  dir.create(parent)
  lock_path <- file.path(parent, ".target.cellucid.lock")
  lock <- cellucid:::.native_export_lock_acquire(lock_path)
  expect_false(is.null(lock))

  cluster <- parallel::makeCluster(1L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  worker_can_lock <- function() {
    parallel::clusterCall(
      cluster,
      function(path) {
        library(cellucid)
        worker_lock <- cellucid:::.native_export_lock_acquire(path)
        if (is.null(worker_lock)) {
          return(FALSE)
        }
        stopifnot(identical(
          cellucid:::.native_export_lock_release(worker_lock),
          c(0L, 0L, 1L)
        ))
        TRUE
      },
      lock_path
    )[[1]]
  }

  expect_false(worker_can_lock())
  rm(lock)
  invisible(gc())
  expect_true(worker_can_lock())
})

test_that("escaped native handles remain inert after namespace unload and gc", {
  script_path <- tempfile("cellucid_r_unload_lock_", fileext = ".R")
  library_expression <- paste(deparse(.libPaths()), collapse = " ")
  writeLines(
    c(
      paste0(".libPaths(", library_expression, ")"),
      "library(cellucid)",
      "lock_path <- tempfile('cellucid-r-unload-lock-')",
      "handle <- cellucid:::.native_export_lock_acquire(lock_path)",
      "stopifnot(!is.null(handle))",
      "unloadNamespace('cellucid')",
      "rm(handle)",
      "invisible(gc())",
      "library(cellucid)",
      "replacement <- cellucid:::.native_export_lock_acquire(lock_path)",
      "stopifnot(!is.null(replacement))",
      paste0(
        "stopifnot(identical(",
        "cellucid:::.native_export_lock_release(replacement), ",
        "c(0L, 0L, 1L)))"
      ),
      "unloadNamespace('cellucid')",
      "invisible(gc())",
      "cat('unload-gc-reload-ok\\n')"
    ),
    script_path
  )
  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  output <- suppressWarnings(system2(
    rscript,
    c("--vanilla", shQuote(script_path)),
    stdout = TRUE,
    stderr = TRUE,
    timeout = 30
  ))
  status <- attr(output, "status", exact = TRUE)
  if (is.null(status)) {
    status <- 0L
  }

  expect_identical(as.integer(status), 0L, info = paste(output, collapse = "\n"))
  expect_true(
    any(grepl("unload-gc-reload-ok", output, fixed = TRUE)),
    info = paste(output, collapse = "\n")
  )
})

test_that("abrupt native lock-owner death releases the persistent sidecar", {
  parent <- tempfile("cellucid_r_native_lock_death_")
  dir.create(parent)
  lock_path <- file.path(parent, ".target.cellucid.lock")
  cluster <- parallel::makeCluster(1L)
  cluster_is_live <- TRUE
  on.exit({
    if (cluster_is_live) {
      suppressWarnings(try(parallel::stopCluster(cluster), silent = TRUE))
    }
  }, add = TRUE)
  owner_pid <- parallel::clusterCall(
    cluster,
    function(path) {
      library(cellucid)
      lock <- cellucid:::.native_export_lock_acquire(path)
      stopifnot(!is.null(lock))
      assign(".cellucid_test_export_lock", lock, envir = .GlobalEnv)
      Sys.getpid()
    },
    lock_path
  )[[1]]

  # Windows only implements process termination for SIGTERM. Both branches
  # bypass cluster cleanup so this still exercises abrupt owner death.
  owner_signal <- if (.Platform$OS.type == "windows") {
    tools::SIGTERM
  } else {
    tools::SIGKILL
  }
  expect_true(tools::pskill(owner_pid, owner_signal))
  suppressWarnings(try(parallel::stopCluster(cluster), silent = TRUE))
  cluster_is_live <- FALSE

  recovered_lock <- NULL
  for (attempt in seq_len(500L)) {
    recovered_lock <- cellucid:::.native_export_lock_acquire(lock_path)
    if (!is.null(recovered_lock)) {
      break
    }
    Sys.sleep(0.01)
  }
  expect_false(is.null(recovered_lock))
  expect_identical(
    cellucid:::.native_export_lock_release(recovered_lock),
    c(0L, 0L, 1L)
  )
  expect_true(file.exists(lock_path))
})

test_that("same-process export locks reject recursion without dropping ownership", {
  parent <- tempfile("cellucid_r_process_lock_parent_")
  dir.create(parent)
  target <- file.path(parent, "target")
  other_target <- file.path(parent, "other-target")
  lock_path <- file.path(parent, ".target.cellucid.lock")

  first_lock <- cellucid:::.acquire_export_generation_lock(target)
  first_lock_is_owned <- TRUE
  on.exit({
    if (first_lock_is_owned) {
      cellucid:::.release_export_generation_lock(first_lock)
    }
  }, add = TRUE)

  cluster <- parallel::makeCluster(1L)
  on.exit({
    if (!is.null(cluster)) {
      parallel::stopCluster(cluster)
    }
  }, add = TRUE)
  worker_can_lock <- function(cluster, path) {
    parallel::clusterCall(
      cluster,
      function(lock_path) {
        library(cellucid)
        lock <- cellucid:::.native_export_lock_acquire(lock_path)
        if (is.null(lock)) {
          return(FALSE)
        }
        identical(
          cellucid:::.native_export_lock_release(lock),
          c(0L, 0L, 1L)
        )
      },
      path
    )[[1]]
  }

  expect_false(worker_can_lock(cluster, lock_path))
  expect_error(
    cellucid:::.acquire_export_generation_lock(target),
    "generation.*already active"
  )
  expect_false(worker_can_lock(cluster, lock_path))

  other_lock <- cellucid:::.acquire_export_generation_lock(other_target)
  expect_true(cellucid:::.release_export_generation_lock(other_lock))
  expect_true(cellucid:::.release_export_generation_lock(first_lock))
  first_lock_is_owned <- FALSE
  expect_true(worker_can_lock(cluster, lock_path))
})

test_that("case aliases and Unicode targets retain exact lock ownership", {
  parent <- tempfile("cellucid_r_case_lock_parent_")
  dir.create(parent)
  canonical_target <- file.path(parent, "CaseSensitiveTarget")
  alias_target <- file.path(parent, "casesensitivetarget")
  canonical_lock_path <- file.path(
    parent,
    ".CaseSensitiveTarget.cellucid.lock"
  )
  alias_lock_path <- file.path(
    parent,
    ".casesensitivetarget.cellucid.lock"
  )
  lock <- cellucid:::.acquire_export_generation_lock(canonical_target)
  lock_is_owned <- TRUE
  on.exit({
    if (lock_is_owned) {
      cellucid:::.release_export_generation_lock(lock)
    }
  }, add = TRUE)

  if (file.exists(alias_lock_path)) {
    expect_error(
      cellucid:::.acquire_export_generation_lock(alias_target),
      "generation.*already active"
    )
    cluster <- parallel::makeCluster(1L)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    alias_contended <- parallel::clusterCall(
      cluster,
      function(path) {
        library(cellucid)
        is.null(cellucid:::.native_export_lock_acquire(path))
      },
      alias_lock_path
    )[[1]]
    expect_true(alias_contended)
  }

  expect_true(cellucid:::.release_export_generation_lock(lock))
  lock_is_owned <- FALSE
  expect_true(file.exists(canonical_lock_path))

  unicode_target <- file.path(parent, "células-细胞")
  unicode_lock <- cellucid:::.acquire_export_generation_lock(unicode_target)
  expect_true(cellucid:::.release_export_generation_lock(unicode_lock))
  unicode_lock_path <- file.path(parent, ".células-细胞.cellucid.lock")
  expect_true(file.exists(unicode_lock_path))
  expect_identical(unname(file.info(unicode_lock_path)$size), 0)
})

test_that("native lock handles reject foreign pointers and real contention", {
  parent <- tempfile("cellucid_r_native_lock_alias_")
  dir.create(parent)
  target <- file.path(parent, "target")
  lock_path <- file.path(parent, ".target.cellucid.lock")
  lock <- cellucid:::.acquire_export_generation_lock(target)
  lock_is_owned <- TRUE
  on.exit({
    if (lock_is_owned) {
      cellucid:::.release_export_generation_lock(lock)
    }
  }, add = TRUE)

  expect_error(
    cellucid:::.native_export_lock_release(new("externalptr")),
    "not owned by Cellucid"
  )

  cluster <- parallel::makeCluster(1L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  alias_contended <- parallel::clusterCall(
    cluster,
    function(path) {
      library(cellucid)
      is.null(cellucid:::.native_export_lock_acquire(path))
    },
    lock_path
  )[[1]]
  expect_true(alias_contended)

  expect_true(cellucid:::.release_export_generation_lock(lock))
  lock_is_owned <- FALSE
})

test_that("hard-linked lock aliases are rejected without mutating ownership", {
  skip_on_os("windows")
  parent <- tempfile("cellucid_r_native_hardlink_alias_")
  dir.create(parent)
  target <- file.path(parent, "target")
  lock_path <- file.path(parent, ".target.cellucid.lock")
  alias_lock_path <- file.path(parent, ".target-alias.cellucid.lock")
  lock <- cellucid:::.acquire_export_generation_lock(target)
  lock_is_owned <- TRUE
  on.exit({
    if (lock_is_owned) {
      cellucid:::.release_export_generation_lock(lock)
    }
    unlink(alias_lock_path)
  }, add = TRUE)

  expect_true(file.link(lock_path, alias_lock_path))
  expect_error(
    cellucid:::.native_export_lock_acquire(alias_lock_path),
    "non-linked regular file"
  )
  expect_true(file.exists(lock_path))
  expect_true(file.exists(alias_lock_path))
  expect_identical(unlink(alias_lock_path), 0L)
  expect_false(file.exists(alias_lock_path))
  expect_true(cellucid:::.release_export_generation_lock(lock))
  lock_is_owned <- FALSE
})

test_that("forked children discard inherited process-local lock claims", {
  skip_on_os("windows")
  parent <- tempfile("cellucid_r_fork_lock_parent_")
  dir.create(parent)
  target <- file.path(parent, "target")
  first_lock <- cellucid:::.acquire_export_generation_lock(target)
  on.exit(
    cellucid:::.release_export_generation_lock(first_lock),
    add = TRUE
  )

  child <- parallel::mcparallel({
    result <- tryCatch(
      cellucid:::.acquire_export_generation_lock(target),
      error = conditionMessage
    )
    list(
      result = result,
      process_id = Sys.getpid(),
      registry_pid = cellucid:::.export_generation_lock_registry$pid,
      claims = ls(
        envir = cellucid:::.export_generation_lock_registry$paths,
        all.names = TRUE
      )
    )
  })
  child_result <- parallel::mccollect(child)[[1]]

  expect_match(child_result$result, "generation.*already active")
  expect_identical(child_result$registry_pid, child_result$process_id)
  expect_identical(child_result$claims, character())
})

test_that("a forked child obtains a real lock after the parent releases", {
  skip_on_os("windows")
  parent <- tempfile("cellucid_r_fork_transfer_parent_")
  dir.create(parent)
  target <- file.path(parent, "target")
  lock_path <- file.path(parent, ".target.cellucid.lock")
  start_path <- file.path(parent, "start-child")
  ready_path <- file.path(parent, "child-ready")
  release_path <- file.path(parent, "release-child")
  wait_for_file <- function(path) {
    for (attempt in seq_len(500L)) {
      if (file.exists(path)) {
        return(TRUE)
      }
      Sys.sleep(0.01)
    }
    FALSE
  }

  parent_lock <- cellucid:::.acquire_export_generation_lock(target)
  parent_lock_is_owned <- TRUE
  on.exit({
    if (parent_lock_is_owned) {
      cellucid:::.release_export_generation_lock(parent_lock)
    }
  }, add = TRUE)

  child <- parallel::mcparallel({
    if (!wait_for_file(start_path)) {
      stop("Parent never released the child start gate.")
    }
    child_lock <- cellucid:::.acquire_export_generation_lock(target)
    on.exit(
      cellucid:::.release_export_generation_lock(child_lock),
      add = TRUE
    )
    file.create(ready_path)
    if (!wait_for_file(release_path)) {
      stop("Parent never released the child lock gate.")
    }
    "real-lock-held"
  })
  child_is_active <- TRUE
  on.exit({
    if (child_is_active) {
      file.create(start_path)
      file.create(release_path)
      parallel::mccollect(child, wait = FALSE)
    }
  }, add = TRUE)

  expect_true(cellucid:::.release_export_generation_lock(parent_lock))
  parent_lock_is_owned <- FALSE
  file.create(start_path)
  expect_true(wait_for_file(ready_path))

  cluster <- parallel::makeCluster(1L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  worker_can_lock <- parallel::clusterCall(
    cluster,
    function(path) {
      library(cellucid)
      lock <- cellucid:::.native_export_lock_acquire(path)
      if (is.null(lock)) {
        return(FALSE)
      }
      identical(
        cellucid:::.native_export_lock_release(lock),
        c(0L, 0L, 1L)
      )
    },
    lock_path
  )[[1]]
  expect_false(worker_can_lock)

  file.create(release_path)
  child_result <- parallel::mccollect(child)[[1]]
  child_is_active <- FALSE
  expect_identical(child_result, "real-lock-held")
})

test_that("stage cleanup failure cannot strand export lock ownership", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap1 <- matrix(c(0, 1), ncol = 1)
  out <- tempfile("cellucid_r_cleanup_failure_")
  cleanup_was_attempted <- FALSE

  testthat::with_mocked_bindings(
    captured_error <- tryCatch(
      cellucid_prepare(
        dataset_id = "cleanup-failure",
        dataset_name = "Cleanup failure",
        latent_space = latent,
        obs = obs,
        X_umap_1d = umap1,
        out_dir = out,
        centroid_min_points = 1,
        force = TRUE,
        obs_categorical_dtype = "uint16"
      ),
      error = identity
    ),
    .normalize_embedding = function(...) {
      stop("synthetic generation failure")
    },
    .remove_export_tree = function(...) {
      cleanup_was_attempted <<- TRUE
      stop("synthetic stage cleanup failure")
    },
    .package = "cellucid"
  )

  expect_true(cleanup_was_attempted)
  expect_s3_class(captured_error, "error")
  expect_match(
    conditionMessage(captured_error),
    "synthetic stage cleanup failure"
  )
  retry_lock <- cellucid:::.acquire_export_generation_lock(out)
  expect_true(cellucid:::.release_export_generation_lock(retry_lock))
  stage_prefix <- paste0(".", basename(out), ".cellucid-stage-")
  stages <- list.files(
    dirname(out),
    pattern = paste0("^", stage_prefix),
    full.names = TRUE,
    all.files = TRUE
  )
  expect_length(stages, 1L)
  cellucid:::.recover_export_transaction(out)
  expect_false(file.exists(stages))
  controls <- cellucid:::.export_transaction_control_paths(out)
  expect_false(file.exists(controls$journal))
  expect_false(file.exists(controls$journal_temp))
})

test_that("legacy lock contents survive successful publication", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap1 <- matrix(c(0, 1), ncol = 1)
  out <- tempfile("cellucid_r_legacy_lock_")
  lock_path <- file.path(
    dirname(out),
    paste0(".", basename(out), ".cellucid.lock")
  )
  writeLines("legacy-owner-pid=424242", lock_path)

  cellucid_prepare(
    dataset_id = "legacy-lock-generation",
    dataset_name = "Legacy lock generation",
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  expect_identical(readLines(lock_path), "legacy-owner-pid=424242")
})

test_that("unsafe export lock paths are rejected without mutation", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap1 <- matrix(c(0, 1), ncol = 1)
  prepare_target <- function(out) {
    cellucid_prepare(
      dataset_id = "unsafe-lock-generation",
      dataset_name = "Unsafe lock generation",
      latent_space = latent,
      obs = obs,
      X_umap_1d = umap1,
      out_dir = out,
      centroid_min_points = 1,
      force = TRUE,
      obs_categorical_dtype = "uint16"
    )
  }

  directory_target <- tempfile("cellucid_r_directory_lock_")
  directory_lock <- file.path(
    dirname(directory_target),
    paste0(".", basename(directory_target), ".cellucid.lock")
  )
  dir.create(directory_lock)
  expect_error(prepare_target(directory_target), "regular file")
  expect_false(file.exists(directory_target))
  expect_true(dir.exists(directory_lock))

  symlink_target <- tempfile("cellucid_r_symlink_lock_")
  symlink_lock <- file.path(
    dirname(symlink_target),
    paste0(".", basename(symlink_target), ".cellucid.lock")
  )
  victim <- tempfile("cellucid_r_lock_victim_")
  writeLines("must remain unchanged", victim)
  skip_if_not(
    file.symlink(victim, symlink_lock),
    "symbolic links unavailable"
  )
  expect_identical(cellucid:::.export_path_info(symlink_lock)$kind, 3L)
  if (.Platform$OS.type != "windows") {
    expect_true(nzchar(Sys.readlink(symlink_lock)))
  }
  expect_error(prepare_target(symlink_target), "symbolic link")
  expect_false(file.exists(symlink_target))
  expect_identical(cellucid:::.export_path_info(symlink_lock)$kind, 3L)
  expect_identical(readLines(victim), "must remain unchanged")
  cleanup_status <- NULL
  expect_warning(
    cleanup_status <- .remove_test_reparse(symlink_lock),
    regexp = NA
  )
  expect_identical(cleanup_status, 0L)
  expect_identical(cellucid:::.export_path_info(symlink_lock)$kind, 0L)
  expect_identical(readLines(victim), "must remain unchanged")
})

test_that("dataset_id is the one identifier that still names a directory", {
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
    cellucid:::.validate_dataset_id("bad id!"),
    "^dataset_id 'bad id!' is not a portable identifier\\. "
  )
  expect_identical(cellucid:::.validate_dataset_id("neutral-paths"), "neutral-paths")
})

test_that("identifier axes are held to distinctness and to being drawable", {
  # The nouns cellucid_prepare() actually passes. They are singular because
  # they are read as the subject of "<what> identifier at position N ...",
  # which is the sentence cellucid-python composes from the same nouns.
  callers <- c(
    "Gene",
    "Observation field",
    "Vector field"
  )

  for (what in callers) {
    expect_error(
      cellucid:::.require_field_identities(c("Gene", "Gene"), what),
      paste0("^", what, " key 'Gene' is duplicated\\.$")
    )
    expect_error(
      cellucid:::.require_field_identities(c("ok", ""), what),
      paste0("^", what, " keys must contain only non-missing, non-empty strings\\.$")
    )
    # A gene name and a category label are drawn in the same legend, so an
    # identifier obeys the rule every category label obeys.
    expect_error(
      cellucid:::.require_field_identities(c("ok", "Liver "), what),
      paste0(
        "^", what, " identifier at position 1 is displayed verbatim, so it ",
        "must not carry characters that have no glyph"
      )
    )
    expect_error(
      cellucid:::.require_field_identities(paste0("a", "\u200b", "b"), what),
      paste0("^", what, " identifier at position 0 is displayed verbatim")
    )
    # Payload paths are indices, so an identifier's spelling can no longer
    # collide with another identifier's file on any filesystem.
    expect_identical(
      cellucid:::.require_field_identities(
        c("Gene", "gene", "bad id!", "CON", "trailing.", "HLA-DRB1/2"),
        what
      ),
      c("Gene", "gene", "bad id!", "CON", "trailing.", "HLA-DRB1/2")
    )
  }
})

test_that("identifiers no filesystem could hold are exported verbatim", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(
    CON = factor(c("A", "A", "B")),
    `sample score (%)` = c(0.1, 0.2, 0.3),
    check.names = FALSE
  )
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  expr <- matrix(c(0, 1, 2, 3, 4, 5), nrow = 3, ncol = 2)
  var <- data.frame(symbol = c("Gene", "gene"))
  rownames(var) <- c("row_1", "row_2")

  out <- tempfile("cellucid_r_unportable_identifiers_")
  cellucid_prepare(
    dataset_id = "unportable-identifiers",
    dataset_name = "Unportable identifiers",
    latent_space = latent,
    obs = obs,
    var = var,
    var_gene_id_column = "symbol",
    gene_expression = expr,
    X_umap_2d = umap2,
    vector_fields = list(
      CON.velocity_umap_2d = matrix(c(1, 0, 1, 0, 1, 0), ncol = 2, byrow = TRUE)
    ),
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  obs_manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(obs_manifest[["_categoricalFields"]][[1L]][[2L]], "CON")
  expect_identical(
    obs_manifest[["_continuousFields"]][[1L]],
    list(1L, "sample score (%)")
  )

  var_manifest <- jsonlite::read_json(
    file.path(out, "var_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(
    var_manifest$fields,
    list(list(0L, "Gene"), list(1L, "gene"))
  )

  identity <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  expect_identical(identity$vector_fields$default_field, "CON.velocity_umap")
  expect_identical(
    identity$vector_fields$fields[["CON.velocity_umap"]]$files[["2d"]],
    "vectors/0_2d.bin"
  )
  expect_identical(list.files(file.path(out, "vectors")), "0_2d.bin")
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

test_that("created_at is exact, reproducible, and validated before mutation", {
  expect_null(formals(cellucid_prepare)$created_at)
  invalid <- list(
    123,
    TRUE,
    NA_character_,
    factor("2026-07-26T12:34:56Z"),
    c("2026-07-26T12:34:56Z", "2026-07-26T12:34:57Z"),
    "",
    "2026-7-26T12:34:56Z",
    "2026-07-26 12:34:56Z",
    "2026-07-26T12:34:56+00:00",
    "2026-07-26T12:34:56.000Z",
    "0000-01-01T00:00:00Z",
    "2026-02-29T12:34:56Z",
    "2026-13-01T12:34:56Z",
    "2026-01-00T12:34:56Z",
    "2026-07-26T24:00:00Z",
    "2026-07-26T12:60:00Z",
    "2026-07-26T12:34:60Z"
  )
  for (value in invalid) {
    parent <- tempfile("cellucid_r_invalid_created_at_")
    expect_error(
      cellucid_prepare(
        dataset_id = "invalid-created-at",
        dataset_name = "Invalid created at",
        out_dir = file.path(parent, "generation"),
        obs_categorical_dtype = "uint16",
        created_at = value
      ),
      "created_at"
    )
    expect_false(dir.exists(parent))
  }

  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap1 <- matrix(c(0, 1), ncol = 1)
  fixed <- "2024-02-29T23:59:59Z"
  out <- tempfile("cellucid_r_exact_created_at_")
  cellucid_prepare(
    dataset_id = "exact-created-at",
    dataset_name = "Exact created at",
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16",
    created_at = fixed
  )
  identity <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  expect_identical(identity$created_at, fixed)

  default_out <- tempfile("cellucid_r_default_created_at_")
  cellucid_prepare(
    dataset_id = "default-created-at",
    dataset_name = "Default created at",
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = default_out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )
  default_identity <- jsonlite::read_json(
    file.path(default_out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  expect_match(
    default_identity$created_at,
    "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$"
  )
  expect_identical(
    cellucid:::.resolve_created_at(default_identity$created_at),
    default_identity$created_at
  )
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
  invalid_identity <- c(
    invalid_identity,
    lapply(c(seq_len(31L), 127L), function(code) {
      list(
        dataset_id = "dataset",
        dataset_name = paste0("Data", intToUtf8(code), "set")
      )
    })
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

  invalid_sources <- list(
    list(
      arguments = list(source_url = "https://example.test/data"),
      message = "source_name"
    ),
    list(
      arguments = list(source_citation = "Citation"),
      message = "source_name"
    ),
    list(
      arguments = list(
        source_url = "https://example.test/data",
        source_citation = "Citation"
      ),
      message = "source_name"
    ),
    list(
      arguments = list(source_name = 1),
      message = "source_name"
    ),
    list(
      arguments = list(source_name = ""),
      message = "source_name"
    ),
    list(
      arguments = list(source_name = "Source", source_url = 1),
      message = "source_url"
    ),
    list(
      arguments = list(source_name = "Source", source_url = ""),
      message = "source_url"
    ),
    list(
      arguments = list(source_name = "Source", source_citation = 1),
      message = "source_citation"
    ),
    list(
      arguments = list(source_name = "Source", source_citation = ""),
      message = "source_citation"
    )
  )
  for (case in invalid_sources) {
    parent <- tempfile("cellucid_r_invalid_source_")
    out <- file.path(parent, "generation")
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
          case$arguments
        )
      ),
      case$message
    )
    expect_false(dir.exists(parent))
  }

  valid_sources <- list(
    list(
      arguments = list(source_name = "Source"),
      expected = list(name = "Source")
    ),
    list(
      arguments = list(
        source_name = "Source",
        source_url = "https://example.test/data"
      ),
      expected = list(
        name = "Source",
        url = "https://example.test/data"
      )
    ),
    list(
      arguments = list(
        source_name = "Source",
        source_citation = "Citation"
      ),
      expected = list(
        name = "Source",
        citation = "Citation"
      )
    ),
    list(
      arguments = list(
        source_name = "Source",
        source_url = "https://example.test/data",
        source_citation = "Citation"
      ),
      expected = list(
        name = "Source",
        url = "https://example.test/data",
        citation = "Citation"
      )
    )
  )
  for (case in valid_sources) {
    out <- tempfile("cellucid_r_valid_source_")
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
        case$arguments
      )
    )
    identity <- jsonlite::read_json(
      file.path(out, "dataset_identity.json"),
      simplifyVector = FALSE
    )
    expect_identical(identity$source, case$expected)
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
    manifest[["_categoricalFields"]][[1]][[3]],
    list("A")
  )
  expect_identical(
    manifest[["_categoricalFields"]][[1]][[6]][["1"]][[1]]$position,
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

test_that("quantization owns the viewer-visible float32 value domain", {
  roundtrip_float32 <- function(values) {
    readBin(
      writeBin(
        as.double(values),
        raw(),
        size = 4L,
        endian = "little"
      ),
      what = double(),
      n = length(values),
      size = 4L,
      endian = "little"
    )
  }

  # Native-double variation finer than float32 resolution is one float32 value,
  # which is one constant, and constants have an encoding.
  collapsed <- cellucid:::.quantize_continuous(
    c(1, 1 + 1e-8, 1 + 2e-8),
    bits = 8L,
    field_name = "collapsed"
  )
  expect_identical(collapsed$quantized, c(0L, 0L, 0L))
  expect_identical(collapsed$min_val, 1)
  expect_identical(collapsed$max_val, 1)
  expect_identical(collapsed$scale, 0)

  expect_error(
    cellucid:::.quantize_continuous(
      c(0, 2^128),
      bits = 8L,
      field_name = "overflow"
    ),
    "remain finite.*float32 domain"
  )

  visible_values <- roundtrip_float32(
    c(1, 1 + 1e-7, 1 + 2e-7)
  )
  q <- cellucid:::.quantize_continuous(
    visible_values,
    bits = 8L,
    field_name = "visible"
  )
  expect_identical(q$quantized, c(0L, 127L, 254L))
  expect_identical(q$min_val, visible_values[[1L]])
  expect_identical(q$max_val, visible_values[[3L]])
  decoded <- q$min_val +
    q$quantized / 254 * (q$max_val - q$min_val)
  expect_identical(roundtrip_float32(decoded), visible_values)

  underflow_source <- c(2^-150, 2^-149, 2^-148)
  underflow_visible <- roundtrip_float32(underflow_source)
  expect_identical(underflow_visible, c(0, 2^-149, 2^-148))
  underflow_q <- cellucid:::.quantize_continuous(
    underflow_source,
    bits = 8L,
    field_name = "subnormal"
  )
  expect_identical(underflow_q$quantized, c(0L, 127L, 254L))
  expect_identical(underflow_q$min_val, 0)
  expect_identical(underflow_q$max_val, 2^-148)
  underflow_decoded <- underflow_q$min_val +
    underflow_q$quantized / 254 *
      (underflow_q$max_val - underflow_q$min_val)
  expect_identical(
    roundtrip_float32(underflow_decoded),
    underflow_visible
  )

  observed_chunk_lengths <- integer()
  original_roundtrip <- cellucid:::.roundtrip_finite_float32_chunk
  chunked <- testthat::with_mocked_bindings(
    cellucid:::.quantize_continuous(
      visible_values,
      bits = 8L,
      field_name = "chunked",
      .chunk_size = 2L
    ),
    .roundtrip_finite_float32_chunk = function(values, field_name) {
      observed_chunk_lengths <<- c(
        observed_chunk_lengths,
        length(values)
      )
      original_roundtrip(values, field_name)
    },
    .package = "cellucid"
  )
  expect_identical(chunked, q)
  expect_identical(observed_chunk_lengths, c(2L, 1L, 2L, 1L))
})

# A gene expressed at one level in every published cell -- very often zero, once
# an atlas is subset to one lineage -- is ordinary scientific data, and so is an
# obs column a subset flattened. compact_v1 publishes such a field as equal
# bounds with every code 0, so dequantization returns the exact constant
# instead of an approximation and nothing divides by (maxValue - minValue).
# The Python exporter is asserted against the same contract in
# cellucid-python/tests/test_prepare_data_constant_field_contract.py.
.constant_field_export <- function(out, ...) {
  n_cells <- 6L
  cellucid_prepare(
    out_dir = out,
    dataset_id = "constant-fields",
    dataset_name = "Constant fields",
    latent_space = matrix(
      c(0, 0, 3, 0, 0.5, 0, 7, 0, 1.5, 0, 9, 0),
      ncol = 2,
      byrow = TRUE
    ),
    X_umap_2d = matrix(
      c(0, 0, 1, 0, 0, 1, 1, 1, 0.5, 0.25, 0.25, 0.75),
      ncol = 2,
      byrow = TRUE
    ),
    obs = data.frame(
      group = factor(c("a", "a", "a", "b", "b", "b")),
      flat = rep(1.5, n_cells),
      zeroed = rep(0, n_cells),
      varies = seq(0, 1, length.out = n_cells)
    ),
    var = data.frame(row.names = c("ZERO", "CONST", "NORM")),
    gene_expression = cbind(
      rep(0, n_cells),
      rep(2.5, n_cells),
      seq(0.5, 1.5, length.out = n_cells)
    ),
    centroid_min_points = 1L,
    obs_categorical_dtype = "uint16",
    force = TRUE,
    ...
  )
  out
}

.read_payload_codes <- function(out, path_pattern, payload_index, size, n) {
  relative <- sub("{index}", payload_index, path_pattern, fixed = TRUE)
  connection <- gzfile(file.path(out, relative), "rb")
  on.exit(close(connection))
  readBin(
    connection,
    what = "integer",
    size = size,
    n = n,
    signed = FALSE,
    endian = "little"
  )
}

.read_payload_float32 <- function(out, path_pattern, payload_index, n) {
  relative <- sub("{index}", payload_index, path_pattern, fixed = TRUE)
  connection <- gzfile(file.path(out, relative), "rb")
  on.exit(close(connection))
  readBin(
    connection,
    what = "double",
    size = 4L,
    n = n,
    endian = "little"
  )
}

# Decode exactly as dequantize() in the viewer's data-loaders.js does: the
# expression is evaluated in double precision and stored into a Float32Array.
.viewer_dequantize <- function(codes, minimum, maximum, bits) {
  max_quant <- if (bits == 8L) 254 else 65534
  decoded <- minimum + codes * ((maximum - minimum) / max_quant)
  readBin(
    writeBin(decoded, raw(), size = 4L, endian = "little"),
    what = double(),
    n = length(decoded),
    size = 4L,
    endian = "little"
  )
}

.manifest_entry <- function(entries, name) {
  for (entry in entries) {
    if (identical(entry[[2L]], name)) {
      return(entry)
    }
  }
  stop("no manifest entry named '", name, "'")
}

# JSON carries one number type, so a whole-valued bound parses back as an R
# integer. The published value is the number, not the R storage mode it lands
# in, and the viewer reads both as one Number.
.manifest_bounds <- function(entry) {
  c(as.double(entry[[3L]]), as.double(entry[[4L]]))
}

for (bits in c(8L, 16L)) {
  local({
    bits <- bits
    test_that(
      paste0(
        "constant genes and obs fields round-trip exactly at ",
        bits,
        " bits"
      ),
      {
        n_cells <- 6L
        out <- .constant_field_export(
          tempfile(paste0("cellucid_r_constant_", bits, "bit_")),
          var_quantization = bits,
          obs_continuous_quantization = bits,
          compression = 6
        )
        code_size <- if (bits == 8L) 1L else 2L

        var_manifest <- jsonlite::read_json(
          file.path(out, "var_manifest.json"),
          simplifyVector = FALSE
        )
        var_pattern <- var_manifest[["_varSchema"]]$pathPattern
        obs_manifest <- jsonlite::read_json(
          file.path(out, "obs_manifest.json"),
          simplifyVector = FALSE
        )
        obs_pattern <- obs_manifest[["_obsSchemas"]]$continuous$pathPattern

        constants <- list(
          list(entries = var_manifest$fields, pattern = var_pattern,
               name = "ZERO", value = 0),
          list(entries = var_manifest$fields, pattern = var_pattern,
               name = "CONST", value = 2.5),
          list(entries = obs_manifest[["_continuousFields"]],
               pattern = obs_pattern, name = "flat", value = 1.5),
          list(entries = obs_manifest[["_continuousFields"]],
               pattern = obs_pattern, name = "zeroed", value = 0)
        )
        for (constant in constants) {
          entry <- .manifest_entry(constant$entries, constant$name)
          expect_length(entry, 4L)
          bounds <- .manifest_bounds(entry)
          # The manifest declares the constant case: equal bounds, both the
          # constant itself.
          expect_identical(bounds[[1L]], constant$value, info = constant$name)
          expect_identical(bounds[[2L]], constant$value, info = constant$name)
          codes <- .read_payload_codes(
            out,
            constant$pattern,
            entry[[1L]],
            size = code_size,
            n = n_cells
          )
          # Every code is 0: the writer never derived a scale for this field.
          expect_identical(codes, integer(n_cells), info = constant$name)
          # The general dequantization arithmetic stays exact and finite: it
          # divides by maxQuant, never by (maxValue - minValue).
          decoded <- .viewer_dequantize(
            codes,
            bounds[[1L]],
            bounds[[2L]],
            bits
          )
          expect_identical(
            decoded,
            rep(constant$value, n_cells),
            info = constant$name
          )
          expect_true(all(is.finite(decoded)), info = constant$name)
        }

        # A field the writer took the general path for still reaches both
        # terminal codes and keeps distinct bounds.
        max_quant <- if (bits == 8L) 254L else 65534L
        for (varying in list(
          list(entries = var_manifest$fields, pattern = var_pattern,
               name = "NORM"),
          list(entries = obs_manifest[["_continuousFields"]],
               pattern = obs_pattern, name = "varies")
        )) {
          entry <- .manifest_entry(varying$entries, varying$name)
          bounds <- .manifest_bounds(entry)
          expect_lt(bounds[[1L]], bounds[[2L]])
          codes <- .read_payload_codes(
            out,
            varying$pattern,
            entry[[1L]],
            size = code_size,
            n = n_cells
          )
          expect_identical(min(codes), 0L, info = varying$name)
          expect_identical(max(codes), max_quant, info = varying$name)
        }
      }
    )
  })
}

test_that("constant genes and obs fields round-trip exactly unquantized", {
  n_cells <- 6L
  out <- .constant_field_export(
    tempfile("cellucid_r_constant_f32_"),
    var_quantization = NULL,
    obs_continuous_quantization = NULL
  )

  var_manifest <- jsonlite::read_json(
    file.path(out, "var_manifest.json"),
    simplifyVector = FALSE
  )
  obs_manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  cases <- list(
    list(entries = var_manifest$fields,
         pattern = var_manifest[["_varSchema"]]$pathPattern,
         name = "ZERO", value = 0),
    list(entries = var_manifest$fields,
         pattern = var_manifest[["_varSchema"]]$pathPattern,
         name = "CONST", value = 2.5),
    list(entries = obs_manifest[["_continuousFields"]],
         pattern = obs_manifest[["_obsSchemas"]]$continuous$pathPattern,
         name = "flat", value = 1.5)
  )
  for (case in cases) {
    entry <- .manifest_entry(case$entries, case$name)
    # A full-precision entry carries no bounds at all.
    expect_length(entry, 2L)
    expect_identical(
      .read_payload_float32(out, case$pattern, entry[[1L]], n_cells),
      rep(case$value, n_cells),
      info = case$name
    )
  }
})

test_that("float32-collapsed variation is published as a constant field", {
  n_cells <- 6L
  values <- 1 + (seq_len(n_cells) - 1) * 1e-8
  out <- tempfile("cellucid_r_float32_collapsed_")
  cellucid_prepare(
    out_dir = out,
    dataset_id = "float32-quantization",
    dataset_name = "Float32 quantization",
    latent_space = matrix(
      c(0, 0, 3, 0, 0.5, 0, 7, 0, 1.5, 0, 9, 0),
      ncol = 2,
      byrow = TRUE
    ),
    X_umap_2d = matrix(
      c(0, 0, 1, 0, 0, 1, 1, 1, 0.5, 0.25, 0.25, 0.75),
      ncol = 2,
      byrow = TRUE
    ),
    obs = data.frame(
      group = factor(c("a", "a", "a", "b", "b", "b")),
      score = values
    ),
    centroid_min_points = 1L,
    obs_categorical_dtype = "uint16",
    obs_continuous_quantization = 8L,
    force = TRUE
  )

  obs_manifest <- jsonlite::read_json(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  entry <- .manifest_entry(obs_manifest[["_continuousFields"]], "score")
  bounds <- .manifest_bounds(entry)
  expect_identical(bounds[[1L]], 1)
  expect_identical(bounds[[2L]], 1)
  codes <- .read_payload_codes(
    out,
    obs_manifest[["_obsSchemas"]]$continuous$pathPattern,
    entry[[1L]],
    size = 1L,
    n = n_cells
  )
  expect_identical(codes, integer(n_cells))
  expect_identical(
    .viewer_dequantize(codes, bounds[[1L]], bounds[[2L]], 8L),
    rep(1, n_cells)
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
    file.path(out, "obs", "0.outliers.u8"),
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

  vec_path <- file.path(out, "vectors", "0_2d.bin")
  expect_true(file.exists(vec_path))

  con <- file(vec_path, open = "rb")
  on.exit(close(con), add = TRUE)
  got <- readBin(con, what = "numeric", size = 4, endian = "little", n = 3 * 2)

  # Input vectors are (2, 0); scale_factor=0.5 -> expected (1, 0)
  expect_equal(got, rep(c(1, 0), 3), tolerance = 1e-6)

  ident <- jsonlite::read_json(file.path(out, "dataset_identity.json"), simplifyVector = TRUE)
  expect_equal(ident$vector_fields$default_field, "velocity_umap")
  expect_true("velocity_umap" %in% names(ident$vector_fields$fields))
  expect_equal(
    ident$vector_fields$fields$velocity_umap$files$`2d`,
    "vectors/0_2d.bin"
  )
})

test_that("a _1d vector field key accepts the plain vector it declares", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))

  # Range is 4 -> scale_factor = 2 / 4 = 0.5
  umap1 <- matrix(c(0, 2, 4), ncol = 1)

  out <- tempfile("cellucid_r_vector_1d_")
  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    vector_fields = list(pseudotime_umap_1d = c(2, 2, 2)),
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  vec_path <- file.path(out, "vectors", "0_1d.bin")
  expect_true(file.exists(vec_path))

  con <- file(vec_path, open = "rb")
  on.exit(close(con), add = TRUE)
  got <- readBin(con, what = "numeric", size = 4, endian = "little", n = 3)
  expect_equal(got, rep(1, 3), tolerance = 1e-6)

  ident <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = TRUE
  )
  expect_equal(ident$vector_fields$default_field, "pseudotime_umap")
  expect_equal(
    ident$vector_fields$fields$pseudotime_umap$available_dimensions,
    1L
  )
})

test_that("vector fields use one exact finite UMAP naming and ownership contract", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap2 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  vector2 <- matrix(c(1, 0, 1, 0), ncol = 2, byrow = TRUE)

  invalid_vectors <- list(
    list(
      value = list(velocity = vector2),
      message = "key 'velocity' must exactly match '<field>_umap_<1\\|2\\|3>d'"
    ),
    list(
      value = list(velocity_umap = vector2),
      message = paste0(
        "key 'velocity_umap' must exactly match ",
        "'<field>_umap_<1\\|2\\|3>d'"
      )
    ),
    list(
      value = list(velocity_umap_4d = vector2),
      message = "key 'velocity_umap_4d' must exactly match"
    ),
    list(
      value = list(velocity_umap_2d = NULL),
      message = "cannot be NULL"
    ),
    list(
      value = list(
        velocity_umap_2d = vector2,
        velocity_umap_2d = vector2
      ),
      message = "names must be unique"
    ),
    list(
      value = list(velocity_umap_2d = c(1, 1)),
      message = "declares 2D but contains 1 components"
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

test_that("the vector field key grammar constrains only the shape", {
  # cellucid-python's _VECTOR_KEY_PATTERN is '^(?P<field>.+_umap)_([123])d$'.
  # The part before _umap is otherwise unconstrained, because it names a field
  # in dataset_identity.json and never a file, so R must not stay stricter.
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap2 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  vector2 <- matrix(c(1, 0, 1, 0), ncol = 2, byrow = TRUE)
  fields <- list(vector2, vector2, vector2, vector2)
  names(fields) <- c(
    "RNA velocity/latent_umap_2d",
    ".leading_dot_umap_2d",
    "-leading_dash_umap_2d",
    "細胞_umap_2d"
  )

  out <- tempfile("cellucid_r_loose_vector_keys_")
  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    X_umap_2d = umap2,
    vector_fields = fields,
    vector_field_default = ".leading_dot_umap",
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  identity <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = FALSE
  )
  # Code-point order: '-' (0x2D) < '.' (0x2E) < 'R' (0x52) < non-ASCII.
  expect_identical(
    names(identity$vector_fields$fields),
    c(
      "-leading_dash_umap",
      ".leading_dot_umap",
      "RNA velocity/latent_umap",
      "細胞_umap"
    )
  )
  expect_identical(
    vapply(
      identity$vector_fields$fields,
      function(field) field$files[["2d"]],
      character(1),
      USE.NAMES = FALSE
    ),
    sprintf("vectors/%d_2d.bin", 0:3)
  )
  expect_identical(
    sort(list.files(file.path(out, "vectors"))),
    sprintf("%d_2d.bin", 0:3)
  )
})

test_that("multiple vector fields require one explicit exact default", {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap2 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  vector2 <- matrix(c(1, 0, 1, 0), ncol = 2, byrow = TRUE)
  fields <- list(
    velocity_umap_2d = vector2,
    displacement_umap_2d = vector2
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

test_that("obs manifest preserves viewer float32 min/max precision in JSON", {
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
  expect_identical(entry[[1]], 0L)
  expect_equal(entry[[2]], "score")
  viewer_bounds <- readBin(
    writeBin(
      c(1 / 3, 2 / 3),
      raw(),
      size = 4L,
      endian = "little"
    ),
    what = double(),
    n = 2L,
    size = 4L,
    endian = "little"
  )
  expect_identical(entry[[3]], viewer_bounds[[1L]])
  expect_identical(entry[[4]], viewer_bounds[[2L]])
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
    "Gene key 'G1' is duplicated\\."
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

  # The subset is renumbered from zero: the index is a position in the export,
  # never a position in var.
  var_dir <- file.path(out, "var")
  expect_identical(list.files(var_dir), "0.values.f32")
  manifest <- jsonlite::read_json(
    file.path(out, "var_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(manifest$fields, list(list(0L, "G2")))

  ident <- jsonlite::read_json(file.path(out, "dataset_identity.json"), simplifyVector = TRUE)
  expect_equal(ident$stats$n_genes, 1L)
})

test_that("gene_identifiers order is the exported payload index order", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(cluster = factor(c("A", "A", "B")))
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  expr <- matrix(as.numeric(seq_len(9L)), nrow = 3, ncol = 3)
  var <- data.frame(symbol = c("G1", "G2", "G3"))
  rownames(var) <- var$symbol

  out <- tempfile("cellucid_r_gene_subset_order_")
  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = latent,
    obs = obs,
    var = var,
    gene_expression = expr,
    gene_identifiers = c("G3", "G1"),
    X_umap_2d = umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  manifest <- jsonlite::read_json(
    file.path(out, "var_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(manifest$fields, list(list(0L, "G3"), list(1L, "G1")))
  expect_identical(
    sort(list.files(file.path(out, "var"))),
    c("0.values.f32", "1.values.f32")
  )

  # var/0 must hold G3's column, not var's first column.
  con <- file(file.path(out, "var", "0.values.f32"), open = "rb")
  on.exit(close(con), add = TRUE)
  expect_equal(
    readBin(con, what = "numeric", size = 4, endian = "little", n = 3L),
    expr[, 3L],
    tolerance = 1e-6
  )
})

.gene_scope_fixture <- function(gene_ids, id_column = NULL) {
  var <- data.frame(symbol = gene_ids, stringsAsFactors = FALSE)
  if (is.null(id_column)) {
    rownames(var) <- gene_ids
  } else {
    names(var) <- id_column
    rownames(var) <- paste0("row_", seq_along(gene_ids))
  }
  list(
    latent = matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE),
    obs = data.frame(cluster = factor(c("A", "A", "B"))),
    umap2 = matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE),
    expr = matrix(
      as.numeric(seq_len(3L * length(gene_ids))),
      nrow = 3L,
      ncol = length(gene_ids)
    ),
    var = var
  )
}

.prepare_gene_scope <- function(fixture, out, gene_identifiers, id_column = NULL) {
  unlink(out, recursive = TRUE, force = TRUE)
  cellucid_prepare(
    dataset_id = "test-dataset",
    dataset_name = "Test dataset",
    latent_space = fixture$latent,
    obs = fixture$obs,
    var = fixture$var,
    gene_expression = fixture$expr,
    var_gene_id_column = id_column,
    gene_identifiers = gene_identifiers,
    X_umap_2d = fixture$umap2,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )
}

test_that("a gene identifier no path could hold is exported and recorded", {
  # Payload paths are indices, so nothing about a gene's spelling can name a
  # file any more. Every one of these was rejected while paths carried names.
  for (gene_id in c("HLA-DRB1/2", "CON", "trailing.", "GENE_B")) {
    fixture <- .gene_scope_fixture(c("gene_a", gene_id, "gene_b"))
    out <- file.path(tempdir(), "cellucid_test_unportable_gene_id")

    .prepare_gene_scope(fixture, out, gene_identifiers = c("gene_a", gene_id))

    manifest <- jsonlite::read_json(
      file.path(out, "var_manifest.json"),
      simplifyVector = FALSE
    )
    expect_identical(
      manifest$fields,
      list(list(0L, "gene_a"), list(1L, gene_id)),
      info = gene_id
    )
    expect_equal(
      sort(list.files(file.path(out, "var"))),
      c("0.values.f32", "1.values.f32"),
      info = gene_id
    )
    unlink(out, recursive = TRUE, force = TRUE)
  }
})

test_that("a duplicate var gene id is rejected even when it is not exported", {
  # Every var row is addressable through gene_identifiers=, so a repeated
  # identifier makes that lookup ambiguous whichever genes are selected.
  fixture <- .gene_scope_fixture(
    c("gene_a", "gene_c", "gene_c"),
    id_column = "gene_id"
  )
  out <- file.path(tempdir(), "cellucid_test_unexported_duplicate_gene_id")

  expect_error(
    .prepare_gene_scope(
      fixture,
      out,
      gene_identifiers = c("gene_a"),
      id_column = "gene_id"
    ),
    "Gene key 'gene_c' is duplicated\\."
  )
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

test_that("obs keys no path could hold are exported and recorded", {
  latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  obs <- data.frame(
    "a/b" = factor(c("A", "A", "B")),
    "a b" = factor(c("A", "A", "B")),
    check.names = FALSE
  )
  umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)

  out <- file.path(tempdir(), "cellucid_test_obs_key_collision")
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
  expect_identical(
    vapply(manifest[["_categoricalFields"]], `[[`, character(1), 2L),
    c("a/b", "a b")
  )
  expect_identical(
    vapply(manifest[["_categoricalFields"]], `[[`, integer(1), 1L),
    c(0L, 1L)
  )
  expect_identical(
    sort(list.files(file.path(out, "obs"))),
    c("0.codes.u16", "0.outliers.f32", "1.codes.u16", "1.outliers.f32")
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

# Every user-supplied string the viewer prints must read as the value it stores.
# A category label, a dataset name, a description, and a source line are drawn
# verbatim in the legend, the field selector, and every exported figure, so a
# character with no glyph changes the value without changing the picture. The
# writer rejects those strings instead of trimming them: trimming rewrites an
# annotation nobody asked to change, and merges "Liver " into a separate "Liver"
# category, moving cells between them without saying so.

.display_text_cases <- list(
  list(
    labels = c("Liver ", "Gut ", "Liver "),
    expected = "ends with U+0020 SPACE"
  ),
  list(
    labels = c(" Liver", "Gut", "Gut"),
    expected = "starts with U+0020 SPACE"
  ),
  list(
    labels = c("   ", "Gut", "Gut"),
    expected = "starts with U+0020 SPACE"
  ),
  list(
    labels = c("", "Gut", "Gut"),
    expected = "'' is empty"
  ),
  list(
    labels = c(paste0("Liver", "\u00a0"), "Gut", "Gut"),
    expected = "ends with U+00A0 NO-BREAK SPACE"
  ),
  list(
    labels = c(paste0("Liver", "\u3000"), "Gut", "Gut"),
    expected = "ends with U+3000 IDEOGRAPHIC SPACE"
  ),
  list(
    labels = c("Li\tver", "Gut", "Gut"),
    expected = "contains U+0009 (control character)"
  ),
  list(
    labels = c("Li\nver", "Gut", "Gut"),
    expected = "contains U+000A (control character)"
  ),
  list(
    labels = c(paste0("Liver", "\u0085"), "Gut", "Gut"),
    expected = "contains U+0085 (control character)"
  ),
  list(
    labels = c(paste0("Liver", "\u200b"), "Gut", "Gut"),
    expected = "contains U+200B ZERO WIDTH SPACE"
  ),
  list(
    labels = c(paste0("\ufeff", "Liver"), "Gut", "Gut"),
    expected = "contains U+FEFF ZERO WIDTH NO-BREAK SPACE"
  ),
  list(
    labels = c(paste0("Liver", "\u2060"), "Gut", "Gut"),
    expected = "contains U+2060 WORD JOINER"
  )
)

.display_text_prepare <- function(out, obs, ...) {
  coordinates <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  args <- list(
    dataset_id = "display-text",
    dataset_name = "Display text",
    latent_space = coordinates,
    obs = obs,
    X_umap_2d = coordinates,
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )
  overrides <- list(...)
  for (name in names(overrides)) {
    args[[name]] <- overrides[[name]]
  }
  do.call(cellucid_prepare, args)
}

test_that("category labels the viewer cannot show verbatim are rejected", {
  for (case in .display_text_cases) {
    out <- tempfile("cellucid_r_display_label_")
    expect_error(
      .display_text_prepare(out, data.frame(organ = case$labels)),
      "Categorical field 'organ'",
      fixed = TRUE,
      info = case$expected
    )
    expect_error(
      .display_text_prepare(out, data.frame(organ = case$labels)),
      case$expected,
      fixed = TRUE
    )
    expect_false(dir.exists(out), info = case$expected)
  }
})

test_that("the published suo organ labels are rejected in one message", {
  published <- c(
    "Bone_marrow ",
    "Gut ",
    "Kidney ",
    "Liver ",
    "Mesenteric_lymph_node",
    "Skin ",
    "Spleen ",
    "Thymus ",
    "Yolk_sac "
  )
  obs <- data.frame(organ = factor(published[c(1L, 4L, 5L)], levels = published))
  out <- tempfile("cellucid_r_display_suo_")

  failure <- tryCatch(
    .display_text_prepare(out, obs),
    error = conditionMessage
  )

  expect_true(grepl("8 labels", failure, fixed = TRUE))
  for (label in published[c(1L, 2L, 3L, 4L)]) {
    expect_true(grepl(paste0("'", label, "'"), failure, fixed = TRUE))
  }
  expect_false(grepl("'Mesenteric_lymph_node'", failure, fixed = TRUE))
  expect_false(dir.exists(out))
})

test_that("a padded label colliding with its clean twin is rejected", {
  out <- tempfile("cellucid_r_display_collision_")
  failure <- tryCatch(
    .display_text_prepare(
      out,
      data.frame(organ = c("Liver", "Liver ", "Liver"))
    ),
    error = conditionMessage
  )

  expect_true(grepl("'Liver '", failure, fixed = TRUE))
  expect_true(grepl("ends with U+0020 SPACE", failure, fixed = TRUE))
  expect_true(grepl("can merge two categories into one", failure, fixed = TRUE))
  expect_false(dir.exists(out))
})

test_that("labels a whitespace-collapsing renderer draws alike are rejected", {
  out <- tempfile("cellucid_r_display_render_")
  failure <- tryCatch(
    .display_text_prepare(
      out,
      data.frame(organ = c("T cell", "T  cell", "T cell"))
    ),
    error = conditionMessage
  )

  expect_true(grepl("'T  cell'", failure, fixed = TRUE))
  expect_true(grepl("'T cell'", failure, fixed = TRUE))
  expect_true(grepl("drawn identically", failure, fixed = TRUE))
  expect_false(dir.exists(out))
})

test_that("the rejection names the field, the reason, and the repair", {
  out <- tempfile("cellucid_r_display_repair_")
  failure <- tryCatch(
    .display_text_prepare(out, data.frame(organ = c("Liver ", "Gut", "Gut"))),
    error = conditionMessage
  )

  expect_true(grepl(
    "the legend, the field selector, and exported figures",
    failure,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Cellucid does not clean them for you",
    failure,
    fixed = TRUE
  ))
  expect_true(grepl(
    "obs[['organ']] <- factor(trimws(as.character(obs[['organ']])",
    failure,
    fixed = TRUE
  ))
})

test_that("labels that read as they are stored survive unchanged", {
  out <- tempfile("cellucid_r_display_clean_")
  obs <- data.frame(
    organ = c("T cell", "Ganglion", "Islet"),
    case = c("Liver", "liver", "LIVER"),
    flag = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  .display_text_prepare(out, obs)

  manifest <- jsonlite::fromJSON(
    file.path(out, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  categories <- lapply(
    manifest[["_categoricalFields"]],
    function(field) unlist(field[[3L]])
  )
  names(categories) <- vapply(
    manifest[["_categoricalFields"]],
    function(field) field[[2L]],
    character(1)
  )
  expect_identical(
    categories$organ,
    c("Ganglion", "Islet", "T cell")
  )
  # Case-only differences are visible on screen, so they stay two categories.
  expect_identical(categories$case, c("LIVER", "Liver", "liver"))
  expect_identical(categories$flag, c(FALSE, TRUE))
})

test_that("identity text the viewer cannot show verbatim is rejected", {
  arguments <- c(
    "dataset_name",
    "dataset_description",
    "source_name",
    "source_url",
    "source_citation"
  )
  values <- list(
    list(value = " padded ", expected = "starts with U+0020 SPACE"),
    list(
      value = paste0("padded", "\u00a0"),
      expected = "ends with U+00A0 NO-BREAK SPACE"
    ),
    list(
      value = paste0("hid", "\u200b", "den"),
      expected = "contains U+200B ZERO WIDTH SPACE"
    ),
    list(
      value = paste0("cont", "\u0001", "rol"),
      expected = "contains U+0001 (control character)"
    )
  )

  for (argument in arguments) {
    for (case in values) {
      out <- tempfile("cellucid_r_display_identity_")
      extra <- list()
      if (argument %in% c("source_url", "source_citation")) {
        extra$source_name <- "Exact source"
      }
      extra[[argument]] <- case$value
      failure <- tryCatch(
        do.call(
          .display_text_prepare,
          c(
            list(out = out, obs = data.frame(group = factor(c("A", "A", "B")))),
            extra
          )
        ),
        error = conditionMessage
      )
      expect_true(
        startsWith(failure, paste0(argument, " is displayed verbatim")),
        info = paste(argument, case$expected)
      )
      expect_true(
        grepl(case$expected, failure, fixed = TRUE),
        info = paste(argument, case$expected)
      )
      expect_false(dir.exists(out), info = paste(argument, case$expected))
    }
  }
})
