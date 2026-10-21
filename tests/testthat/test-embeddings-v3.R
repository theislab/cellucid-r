# Export format v3: any number of embeddings per dimension, each an item with
# an id, label, dimension, file and source; vector fields bound to embedding
# ids; and gene-name sets in var_manifest.json. cellucid-python writes the same
# shapes, and test-writer-parity-v3.R pins the two writers byte for byte.

.v3_base <- function(out, ...) {
  args <- list(
    dataset_id = "format-v3",
    dataset_name = "Format v3",
    latent_space = matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE),
    obs = data.frame(cluster = factor(c("A", "A", "B"))),
    embeddings = list(X_umap = matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)),
    out_dir = out,
    centroid_min_points = 1,
    force = TRUE,
    obs_categorical_dtype = "uint16",
    created_at = "2026-01-01T00:00:00Z"
  )
  overrides <- list(...)
  for (name in names(overrides)) {
    args[name] <- list(overrides[[name]])
  }
  do.call(cellucid_prepare, args)
  out
}

.read_identity <- function(out) {
  jsonlite::read_json(file.path(out, "dataset_identity.json"), simplifyVector = FALSE)
}

test_that("cellucid_embedding validates what it describes", {
  wide <- matrix(as.numeric(seq_len(15L)), nrow = 3L)
  spec <- cellucid_embedding(wide, columns = c(4, 2), label = "PCA", source_key = "X_pca")
  expect_s3_class(spec, "cellucid_embedding")
  expect_identical(spec$columns, c(4L, 2L))

  expect_error(cellucid_embedding(wide), "values has 5 columns, and an embedding plots 1, 2, or 3")
  expect_error(cellucid_embedding(wide, columns = c(1, 2, 3, 4)), "must name 1, 2, or 3 columns")
  expect_error(cellucid_embedding(wide, columns = c(1, 1)), "each axis needs a different column")
  expect_error(cellucid_embedding(wide, columns = c(0, 1)), "names a column that does not exist")
  expect_error(cellucid_embedding(wide, columns = c(1, 6)), "the array has columns 1 to 5")
  expect_error(cellucid_embedding(wide, columns = c(1.5, 2)), "whole column numbers")
  expect_error(cellucid_embedding(wide, columns = c(1, NA)), "numeric vector of 1-based")
  expect_error(cellucid_embedding(c(1, 2, 3), columns = c(1, 2)), "needs at least 2 columns")
  expect_error(cellucid_embedding(wide, columns = 1:2, label = " PCA"), "^label is displayed verbatim")
  expect_error(cellucid_embedding(wide, columns = 1:2, source_key = ""), "source_key must be NULL")
  expect_error(cellucid_embedding(wide, columns = 1:2, source_key = 1), "source_key must be exactly one string")
})

test_that("every embedding is one item, ordered by dimension and id", {
  coords3 <- matrix(c(0, 0, 0, 2, 1, 4, 1, 3, 2), ncol = 3, byrow = TRUE)
  wide <- cbind(coords3, c(5, 1, 3), c(9, 9, 8))
  out <- .v3_base(
    tempfile("cellucid_r_v3_items_"),
    embeddings = list(
      X_umap_3d = coords3,
      X_umap = coords3[, 1:2],
      X_pca = cellucid_embedding(wide, columns = c(5, 4), label = "PCA", source_key = "X_pca"),
      X_line = coords3[, 3]
    ),
    embedding_defaults = c("2d" = "X_umap")
  )
  identity <- .read_identity(out)
  expect_identical(identity$version, 3L)
  embeddings <- identity$embeddings
  expect_identical(names(embeddings), c("default_dimension", "defaults", "items"))
  expect_identical(embeddings$default_dimension, 3L)
  expect_identical(
    embeddings$defaults,
    list(`1d` = "X_line", `2d` = "X_umap", `3d` = "X_umap_3d")
  )
  expect_identical(names(embeddings$items), c("X_line", "X_pca", "X_umap", "X_umap_3d"))
  expect_identical(
    embeddings$items$X_pca,
    list(
      label = "PCA",
      dimension = 2L,
      file = "embeddings/1_2d.bin",
      source = list(key = "X_pca", columns = list(4L, 3L))
    )
  )
  expect_identical(embeddings$items$X_umap$label, "X_umap")
  expect_null(embeddings$items$X_umap$source$key)
  expect_identical(
    sort(list.files(file.path(out, "embeddings"))),
    c("0_1d.bin", "1_2d.bin", "2_2d.bin", "3_3d.bin")
  )

  # Every categorical field carries one centroid list per embedding.
  manifest <- jsonlite::read_json(file.path(out, "obs_manifest.json"), simplifyVector = FALSE)
  expect_identical(manifest[["_format"]], "compact_v2")
  expect_identical(
    names(manifest[["_categoricalFields"]][[1L]][[6L]]),
    c("X_line", "X_pca", "X_umap", "X_umap_3d")
  )
  expect_length(manifest[["_categoricalFields"]][[1L]][[6L]]$X_pca[[1L]]$position, 2L)

  # The selected columns are normalized on their own: column 5 spans 9..8 and
  # column 4 spans 1..5, so the wider one sets the scale 2 / 4.
  connection <- file(file.path(out, "embeddings", "1_2d.bin"), open = "rb")
  on.exit(close(connection), add = TRUE)
  expect_identical(
    readBin(connection, what = "numeric", size = 4L, endian = "little", n = 6L),
    c(0.25, 1, 0.25, -1, -0.25, 0)
  )

  no_centroids <- .v3_base(
    tempfile("cellucid_r_v3_no_centroids_"),
    embeddings = list(X_umap = coords3[, 1:2], X_line = coords3[, 3]),
    centroid_outlier_quantile = NULL
  )
  manifest <- jsonlite::read_json(
    file.path(no_centroids, "obs_manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(
    manifest[["_categoricalFields"]][[1L]][[6L]],
    list(X_line = list(), X_umap = list())
  )
})

test_that("each dimension with several embeddings names its default", {
  umap <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  several <- list(X_umap = umap, X_tsne = umap[, 2:1])
  out <- tempfile("cellucid_r_v3_defaults_")
  expect_error(
    .v3_base(out, embeddings = several),
    paste0(
      "^2 2D embeddings were supplied \\(c\\(\"X_tsne\", \"X_umap\"\\)\\), so ",
      "embedding_defaults must name the one shown first"
    )
  )
  expect_false(dir.exists(out))
  expect_error(
    .v3_base(out, embeddings = several, embedding_defaults = c("3d" = "X_umap")),
    "names a 3d default, and no 3D embedding was supplied"
  )
  expect_error(
    .v3_base(out, embeddings = several, embedding_defaults = c("2D" = "X_umap")),
    "key '2D' must be \"1d\", \"2d\", or \"3d\""
  )
  expect_error(
    .v3_base(out, embeddings = several, embedding_defaults = c("2d" = "X_pca")),
    "is not a 2D embedding"
  )
  expect_error(
    .v3_base(out, embeddings = several, embedding_defaults = "X_umap"),
    "embedding_defaults must map"
  )
  expect_false(dir.exists(out))
  .v3_base(out, embeddings = several, embedding_defaults = list(`2d` = "X_tsne"))
  expect_identical(.read_identity(out)$embeddings$defaults, list(`2d` = "X_tsne"))

  expect_error(
    .v3_base(tempfile(), embeddings = list(X_umap = umap, X_short = umap[1:2, ])),
    "All embeddings must have the same number of cells: 'X_umap' has 3 and 'X_short' has 2"
  )
  expect_error(
    .v3_base(tempfile(), embeddings = list(umap)),
    "^embeddings must be a named list: every embedding needs an id\\.$"
  )
  expect_error(
    .v3_base(tempfile(), embeddings = list(X_umap = umap, umap)),
    "Embedding keys must contain only non-missing, non-empty strings"
  )
  expect_error(
    .v3_base(tempfile(), embeddings = list(X_wide = cbind(umap, umap))),
    "has 4 columns, and an embedding plots 1, 2, or 3"
  )
})

test_that("a vector field binds to several embeddings, in item order", {
  coords3 <- matrix(c(0, 0, 0, 2, 1, 4, 1, 3, 2), ncol = 3, byrow = TRUE)
  out <- .v3_base(
    tempfile("cellucid_r_v3_vector_bindings_"),
    embeddings = list(X_umap_3d = coords3, X_umap = coords3[, 1:2]),
    vector_fields = list(
      velocity = list(X_umap_3d = coords3 / 8, X_umap = coords3[, 1:2] / 8)
    )
  )
  identity <- .read_identity(out)
  expect_identical(
    identity$vector_fields,
    list(
      default_field = "velocity",
      fields = list(velocity = list(
        label = "velocity",
        files = list(X_umap = "vectors/0_0.bin", X_umap_3d = "vectors/0_1.bin")
      ))
    )
  )
  expect_identical(sort(list.files(file.path(out, "vectors"))), c("0_0.bin", "0_1.bin"))
  expect_error(
    .v3_base(tempfile(), vector_field_default = "velocity"),
    "^vector_field_default was provided but no vector fields exist\\.$"
  )
})

test_that("alternative gene names are metadata-only name sets", {
  expr <- matrix(as.numeric(seq_len(9L)), nrow = 3L)
  var <- data.frame(
    symbol = c("CD8A", NA, "CD8A"),
    alias = factor(c("ENSG1", "ENSG2", "ENSG3")),
    row.names = c("ENSG1", "ENSG2", "ENSG3")
  )
  out <- tempfile("cellucid_r_v3_gene_name_sets_")
  expect_message(
    .v3_base(
      out,
      var = var,
      gene_expression = expr,
      var_gene_name_columns = c("symbol", "alias")
    ),
    paste0(
      "^Gene names from var column 'symbol': 1 gene has no 'symbol' value and ",
      "uses its var index name\\. 1 duplicate name was made unique, e\\.g\\. ",
      "'CD8A' \u2192 'CD8A-1'\\.\\n$"
    )
  )
  manifest <- jsonlite::read_json(file.path(out, "var_manifest.json"), simplifyVector = FALSE)
  expect_identical(
    names(manifest),
    c("_format", "n_points", "compression", "quantization", "_varSchema", "fields", "gene_name_sets")
  )
  expect_identical(manifest[["_format"]], "compact_v2")
  expect_identical(
    manifest$fields,
    list(list(0L, "ENSG1"), list(1L, "ENSG2"), list(2L, "ENSG3"))
  )
  expect_identical(
    manifest$gene_name_sets,
    list(
      list(source = NULL, filled_from_index = 0L, renamed = list(), names = NULL),
      list(
        source = "symbol",
        filled_from_index = 1L,
        renamed = list(list("CD8A", "CD8A-1")),
        names = list("CD8A", NULL, "CD8A-1")
      ),
      list(
        source = "alias",
        filled_from_index = 0L,
        renamed = list(),
        names = list(NULL, NULL, NULL)
      )
    )
  )
  expect_identical(sort(list.files(file.path(out, "var"))), sprintf("%d.values.f32", 0:2))

  refusals <- list(
    list(columns = c("symbol", "symbol"), message = "lists 'symbol' more than once"),
    list(columns = "missing", message = "var_gene_name_columns entry 'missing' not found in var"),
    list(columns = 1, message = "must be a character vector of var column names"),
    list(columns = "symbol ", message = "is displayed verbatim")
  )
  for (refusal in refusals) {
    failed <- tempfile("cellucid_r_v3_gene_name_refusal_")
    expect_error(
      .v3_base(
        failed,
        var = var,
        gene_expression = expr,
        var_gene_name_columns = refusal$columns
      ),
      refusal$message,
      info = paste(refusal$columns, collapse = ",")
    )
    expect_false(dir.exists(failed))
  }
  expect_error(
    .v3_base(
      tempfile(),
      var = var,
      gene_expression = expr,
      var_gene_id_column = "alias",
      var_gene_name_columns = "alias"
    ),
    "which var_gene_id_column already uses for the default gene names"
  )
  expect_error(
    .v3_base(tempfile(), var = var, var_gene_name_columns = "symbol"),
    "and no gene_expression was supplied"
  )
})

test_that("an empty gene name needs a row name to take", {
  expr <- matrix(as.numeric(seq_len(6L)), nrow = 3L)
  automatic <- data.frame(symbol = c("CD8A", ""))
  expect_error(
    .v3_base(
      tempfile(),
      var = automatic,
      gene_expression = expr,
      var_gene_id_column = "symbol"
    ),
    paste0(
      "^var column 'symbol': Name at position 1 is empty, and its index name ",
      "is not a non-empty string: NA\\.$"
    )
  )
  out <- tempfile("cellucid_r_v3_filled_default_")
  expect_message(
    .v3_base(
      out,
      var = data.frame(symbol = c("CD8A", ""), row.names = c("ENSG1", "ENSG2")),
      gene_expression = expr,
      var_gene_id_column = "symbol"
    ),
    "1 gene has no 'symbol' value and uses its var index name\\."
  )
  manifest <- jsonlite::read_json(file.path(out, "var_manifest.json"), simplifyVector = FALSE)
  expect_identical(manifest$fields, list(list(0L, "CD8A"), list(1L, "ENSG2")))
})

test_that("a gene-name set manifest is validated against its field names", {
  default_names <- c("A", "B")
  set0 <- list(source = NULL, filled_from_index = 0L, renamed = list(), names = NULL)
  alternative <- function(names, source = "symbol") {
    list(source = source, filled_from_index = 0L, renamed = list(), names = names)
  }
  validate <- function(sets) {
    cellucid:::.require_gene_name_sets(sets, default_names)
  }
  expect_silent(validate(list(set0, alternative(list("x", NULL)))))
  expect_error(validate(list()), "must be a non-empty list")
  expect_error(validate(list(alternative(list("x", NULL)))), "names must be null")
  expect_error(validate(list(set0, alternative(list("x")))), "one entry per gene")
  expect_error(validate(list(set0, alternative(list("A", NULL)))), "null wherever a name equals")
  expect_error(validate(list(set0, alternative(list("B", NULL)))), "gene key 'B' is duplicated")
  expect_error(
    validate(list(set0, alternative(list("x", NULL)), alternative(list("y", NULL)))),
    "must read each source once"
  )
  expect_error(
    validate(list(set0, alternative(list("x", NULL), source = NULL))),
    "must name a var column"
  )
})
