# A manifest is validated in memory and then serialized, and until the bytes
# are read back nothing proves the two describe the same thing. jsonlite is not
# a total function on the payloads this package builds: `auto_unbox = TRUE`
# turns a one-element vector into a scalar, a non-finite double into a JSON
# *string*, and an empty unnamed list into `[]` where the format specifies an
# object. Each of those publishes cleanly and fails only when the viewer opens
# the export, so the writer re-reads every manifest it writes and requires it to
# parse back to exactly the payload that was validated.

.readback_fixture <- function(n_cells = 40L) {
  list(
    latent = matrix(seq_len(n_cells * 3L) / n_cells, ncol = 3L),
    umap2 = matrix(
      c(seq_len(n_cells) / n_cells, rev(seq_len(n_cells)) / n_cells),
      ncol = 2L
    )
  )
}

.readback_prepare <- function(out, obs, ...) {
  fixture <- .readback_fixture()
  args <- list(
    dataset_id = "readback",
    dataset_name = "Read back",
    latent_space = fixture$latent,
    obs = obs,
    X_umap_2d = fixture$umap2,
    out_dir = out,
    centroid_min_points = 2L,
    force = TRUE,
    obs_categorical_dtype = "uint8",
    created_at = "2020-01-01T00:00:00Z"
  )
  # A case that replaces an embedding passes it here, so `...` overrides the
  # default rather than colliding with it.
  overrides <- list(...)
  for (name in names(overrides)) {
    args[name] <- list(overrides[[name]])
  }
  do.call(cellucid_prepare, args)
}

.readback_out_dir <- function(name) {
  out <- file.path(tempdir(), paste0("cellucid_readback_", name))
  unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

test_that("an export with no observation fields declares _obsSchemas as an object", {
  # The web reader calls requireRecord() on _obsSchemas and refuses an array,
  # so an export with no obs field at all -- an obs with no columns, or an
  # explicit empty obs_keys -- must still spell the empty schema map `{}`.
  # jsonlite renders an empty *unnamed* list as `[]`, which published a dataset
  # the viewer then rejected with "expected an object".
  for (variant in c("no-columns", "empty-obs-keys")) {
    out <- .readback_out_dir(variant)
    if (identical(variant, "no-columns")) {
      suppressMessages(capture.output(
        .readback_prepare(out, obs = data.frame(row.names = seq_len(40L)))
      ))
    } else {
      suppressMessages(capture.output(
        .readback_prepare(
          out,
          obs = data.frame(grp = factor(rep(c("A", "B"), each = 20L))),
          obs_keys = character(0)
        )
      ))
    }

    text <- paste(
      readLines(file.path(out, "obs_manifest.json"), warn = FALSE),
      collapse = ""
    )
    expect_true(grepl('"_obsSchemas":{}', text, fixed = TRUE))
    expect_false(grepl('"_obsSchemas":[]', text, fixed = TRUE))

    parsed <- jsonlite::fromJSON(text, simplifyVector = FALSE)
    # An empty JSON object parses to a list carrying character(0) names; an
    # empty JSON array parses to one carrying NULL names. That distinction is
    # the whole defect, so it is asserted rather than the length alone.
    expect_identical(names(parsed[["_obsSchemas"]]), character(0))
    expect_identical(length(parsed[["_obsSchemas"]]), 0L)
    expect_identical(parsed[["_continuousFields"]], list())
    expect_identical(names(parsed[["_continuousFields"]]), NULL)
    expect_identical(parsed[["_categoricalFields"]], list())
    expect_identical(names(parsed[["_categoricalFields"]]), NULL)
  }
})

test_that("a non-finite number in a manifest is refused rather than written as a string", {
  # jsonlite writes NaN, Inf and NA as JSON *strings*, so a bound that went
  # non-finite would publish a well-formed manifest in which minValue is the
  # text "NaN". Nothing downstream can tell that from a number, so the writer
  # refuses it at the read-back.
  path <- file.path(tempdir(), "cellucid_readback_probe.json")
  on.exit(unlink(path, force = TRUE), add = TRUE)

  for (value in list(NaN, Inf, -Inf, NA_real_)) {
    expect_error(
      cellucid:::.write_json(path, list(minValue = value)),
      "does not read back",
      fixed = TRUE
    )
  }
  expect_silent(cellucid:::.write_json(path, list(minValue = 0.5)))
})

test_that("the read-back refuses a file that says something else than the payload", {
  # The guard is proved against the bytes, not against the encoder: given the
  # payload the writer validated, a file carrying any other shape is refused.
  # The two cases are the losses jsonlite's mapping can produce -- an object
  # collapsed to an array, and a one-element array collapsed to a scalar.
  path <- file.path(tempdir(), "cellucid_readback_tampered.json")
  on.exit(unlink(path, force = TRUE), add = TRUE)

  payload <- list(
    schemas = cellucid:::.json_object(),
    categories = I(c("only")),
    n_points = 3L
  )
  expect_silent(cellucid:::.write_json(path, payload))
  expect_error(cellucid:::.require_manifest_reads_back(path, payload), NA)

  tampered <- list(
    object_became_array = '{"schemas":[],"categories":["only"],"n_points":3}',
    array_became_scalar = '{"schemas":{},"categories":"only","n_points":3}',
    key_dropped = '{"schemas":{},"categories":["only"]}',
    key_reordered = '{"categories":["only"],"schemas":{},"n_points":3}',
    value_changed = '{"schemas":{},"categories":["only"],"n_points":4}',
    not_json = '{"schemas":{},'
  )
  for (name in names(tampered)) {
    writeLines(tampered[[name]], con = path, useBytes = TRUE)
    expect_error(
      cellucid:::.require_manifest_reads_back(path, payload),
      "does not read back",
      fixed = TRUE,
      label = name
    )
  }

  # A file that is not valid UTF-8 is refused before it is parsed.
  writeBin(as.raw(c(0x7b, 0x22, 0xff, 0x22, 0x3a, 0x31, 0x7d)), path)
  expect_error(
    cellucid:::.require_manifest_reads_back(path, payload),
    "not valid UTF-8",
    fixed = TRUE
  )
})

test_that("the read-back accepts every manifest the writer actually emits", {
  # The guard must not be so strict that ordinary exports fail: one case per
  # optional part of the format.
  fixture <- .readback_fixture()
  obs <- data.frame(
    grp = factor(rep(c("A", "B"), each = 20L)),
    val = seq_len(40L) / 40
  )
  expr <- matrix(seq_len(120L) / 40, nrow = 40L)
  var <- data.frame(symbol = c("g1", "g2", "g3"))
  rownames(var) <- var$symbol
  connectivities <- matrix(0, 40L, 40L)
  connectivities[1L, 2L] <- connectivities[2L, 1L] <- 0.5
  connectivities[3L, 4L] <- connectivities[4L, 3L] <- 1

  cases <- list(
    plain = list(),
    quantized = list(obs_continuous_quantization = 8L, var_quantization = 8L),
    genes = list(gene_expression = expr, var = var),
    connectivity = list(connectivities = connectivities),
    compressed = list(compression = 6L),
    no_centroids = list(centroid_outlier_quantile = NULL),
    single_category = list(),
    one_dimension = list(
      X_umap_2d = NULL,
      X_umap_1d = matrix(seq_len(40L) / 40, ncol = 1L)
    ),
    vectors = list(
      vector_fields = list(
        velocity_umap_2d = matrix(seq_len(80L) / 40, ncol = 2L)
      ),
      vector_field_default = "velocity_umap"
    )
  )

  for (name in names(cases)) {
    out <- .readback_out_dir(name)
    case_obs <- if (identical(name, "single_category")) {
      data.frame(grp = factor(rep("solo", 40L)), val = seq_len(40L) / 40)
    } else {
      obs
    }
    args <- c(list(out = out, obs = case_obs), cases[[name]])
    expect_error(
      suppressMessages(capture.output(do.call(.readback_prepare, args))),
      NA
    )
    for (manifest in c(
      "obs_manifest.json",
      "var_manifest.json",
      "connectivity_manifest.json",
      "dataset_identity.json"
    )) {
      manifest_path <- file.path(out, manifest)
      if (!file.exists(manifest_path)) {
        next
      }
      expect_error(
        jsonlite::fromJSON(manifest_path, simplifyVector = FALSE),
        NA
      )
    }
  }
})

test_that("every manifest the writer emits is valid UTF-8 whatever the labels are", {
  # Category labels are the one part of a manifest a caller fully controls, and
  # writeLines(useBytes = TRUE) copies bytes rather than re-encoding them.
  latin1 <- c("caf\xe9", "na\xefve")
  Encoding(latin1) <- "latin1"
  labels <- list(
    ascii = c("A", "B"),
    utf8 = c("caf\u00e9", "\u4e2d\u6587"),
    latin1 = latin1
  )
  for (name in names(labels)) {
    out <- .readback_out_dir(paste0("utf8_", name))
    values <- labels[[name]]
    obs <- data.frame(
      grp = factor(values[rep_len(seq_along(values), 40L)], levels = values),
      val = seq_len(40L) / 40
    )
    suppressMessages(capture.output(.readback_prepare(out, obs = obs)))
    path <- file.path(out, "obs_manifest.json")
    bytes <- readBin(path, "raw", file.size(path))
    text <- rawToChar(bytes)
    Encoding(text) <- "UTF-8"
    expect_false(is.na(iconv(text, "UTF-8", "UTF-8")))
    expect_error(jsonlite::fromJSON(text, simplifyVector = FALSE), NA)
  }
})
