# The help page for cellucid_prepare() is hand-written in
# man/cellucid_prepare.Rd and is the only source for it. Nothing in this
# package is generated: NAMESPACE is hand-written too, because
# useDynLib(cellucid, .registration = TRUE, .fixes = "C_") has no generator
# input anywhere under R/ and no generator would emit it. Edit the .Rd file
# directly when this signature changes; tests/testthat/test-current-contract.R
# rebuilds the documented \usage block into a function and compares it to
# these formals, so the two cannot drift apart.
cellucid_prepare <- function(
    latent_space = NULL,
    obs = NULL,
    var = NULL,
    gene_expression = NULL,
    var_gene_id_column = NULL,
    gene_identifiers = NULL,
    connectivities = NULL,
    out_dir = file.path(getwd(), "exports"),
    obs_keys = NULL,
    centroid_outlier_quantile = 0.95,
    centroid_min_points = 10,
    force = FALSE,
    var_quantization = NULL,
    obs_continuous_quantization = NULL,
    obs_categorical_dtype,
    compression = NULL,
    dataset_name = NULL,
    dataset_description = "",
    dataset_id = NULL,
    source_name = NULL,
    source_url = NULL,
    source_citation = NULL,
    X_umap_1d = NULL,
    X_umap_2d = NULL,
    X_umap_3d = NULL,
    vector_fields = NULL,
    vector_field_default = NULL,
    created_at = NULL
) {
  manifest_format_version <- "compact_v1"
  obs_binary_dirname <- "obs"
  var_binary_dirname <- "var"
  connectivity_binary_dirname <- "connectivity"

  force <- .validate_force(force)
  compression <- .normalize_compression(compression)
  obs_categorical_dtype <- .validate_obs_categorical_dtype(
    obs_categorical_dtype
  )
  dataset_id <- .validate_dataset_id(dataset_id)
  dataset_name <- .validate_dataset_name(dataset_name)
  dataset_description <- .validate_display_text(
    dataset_description,
    "dataset_description",
    allow_empty = TRUE
  )
  created_at <- .resolve_created_at(created_at)
  source_info <- .validate_source_identity(
    source_name,
    source_url,
    source_citation
  )

  if (!is.null(var_quantization)) {
    var_quantization <- .validate_quantization_bits(
      var_quantization,
      arg = "var_quantization"
    )
  }
  if (!is.null(obs_continuous_quantization)) {
    obs_continuous_quantization <- .validate_quantization_bits(
      obs_continuous_quantization,
      arg = "obs_continuous_quantization"
    )
  }
  centroid_outlier_quantile <- .validate_centroid_outlier_quantile(
    centroid_outlier_quantile
  )
  centroid_min_points <- .validate_positive_integer(
    centroid_min_points,
    "centroid_min_points"
  )

  final_out_dir <- .validate_output_path(out_dir)
  .dir_create(dirname(final_out_dir))
  out_dir <- NULL
  transaction_cleanup_armed <- FALSE
  export_lock <- .acquire_export_generation_lock(final_out_dir)
  on.exit(
    .release_export_generation_lock(export_lock),
    add = TRUE
  )
  on.exit({
    if (transaction_cleanup_armed) {
      recovery_error <- tryCatch(
        {
          .recover_export_transaction(final_out_dir)
          NULL
        },
        error = identity
      )
      if (inherits(recovery_error, "error")) {
        stop(
          "Failed to recover rejected export transaction for ",
          final_out_dir,
          ": ",
          conditionMessage(recovery_error),
          call. = FALSE
        )
      }
    }
  }, add = TRUE, after = FALSE)
  .recover_export_transaction(final_out_dir)
  transaction_cleanup_armed <- TRUE
  transaction <- .begin_export_transaction(final_out_dir, force)
  out_dir <- transaction$stage
  .dir_create(out_dir)

  obs_binary_dir <- file.path(out_dir, obs_binary_dirname)
  .dir_create(obs_binary_dir)

  embeddings_in <- list(
    `1` = X_umap_1d,
    `2` = X_umap_2d,
    `3` = X_umap_3d
  )
  embeddings <- .collect_embeddings(embeddings_in)
  n_cells <- .validate_embeddings(embeddings)
  connectivity_edges <- NULL
  if (!is.null(connectivities)) {
    connectivity_edges <- .validate_connectivity_edges(
      connectivities = connectivities,
      n_cells = n_cells
    )
  }

  normalization_info <- list()
  for (dim in names(embeddings)) {
    if (any(!is.finite(embeddings[[dim]]))) {
      stop(
        "X_umap_", dim, "d contains non-finite values (NA/NaN/Inf). ",
        "Please filter or impute before export."
      )
    }
    normalized <- .normalize_embedding(embeddings[[dim]])
    embeddings[[dim]] <- normalized$coords
    normalization_info[[dim]] <- normalized$info
  }

  available_dimensions <- sort(as.integer(names(embeddings)))
  default_dimension <- if (3L %in% available_dimensions) {
    3L
  } else if (2L %in% available_dimensions) {
    2L
  } else {
    1L
  }

  # Save embeddings
  for (dim in available_dimensions) {
    file_base <- file.path(out_dir, sprintf("points_%dd.bin", dim))
    .write_float32_matrix_row_major(file_base, embeddings[[as.character(dim)]], compression = compression)
  }

  if (is.null(latent_space)) {
    stop("latent_space is required for outlier quantile calculation.")
  }
  latent <- .as_dense_matrix(latent_space, name = "latent_space")
  if (nrow(latent) != n_cells) {
    stop(
      "Latent space has ", nrow(latent), " cells, but embeddings have ", n_cells, " cells."
    )
  }

  if (!is.data.frame(obs)) {
    stop("obs must be a data.frame.", call. = FALSE)
  }
  if (nrow(obs) != n_cells) {
    stop("obs has ", nrow(obs), " rows, but embeddings have ", n_cells, " cells.")
  }

  if (is.null(obs_keys)) {
    obs_keys <- names(obs)
  }
  .validate_character_vector(obs_keys, what = "obs_keys")
  missing <- setdiff(obs_keys, names(obs))
  if (length(missing) > 0) {
    stop(
      "obs_keys contain columns not in obs: ",
      .format_value_list(missing),
      ". Available columns: ",
      .format_value_list(names(obs))
    )
  }
  .require_field_identities(obs_keys, what = "Observation field")

  # Collect lightweight metadata for obs fields (used in dataset identity)
  obs_field_summaries <- lapply(
    obs_keys,
    function(key) .summarize_obs_field(
      obs[[key]],
      key = key,
      obs_continuous_quantization = obs_continuous_quantization,
      obs_categorical_dtype = obs_categorical_dtype
    )
  )

  # Vector fields (optional)
  vector_fields_identity <- NULL
  vector_payload_indices <- list()
  if (!is.null(vector_fields)) {
    vector_result <- .export_vector_fields(
      vector_fields = vector_fields,
      embeddings = embeddings,
      normalization_info = normalization_info,
      out_dir = out_dir,
      compression = compression,
      vector_field_default = vector_field_default
    )
    vector_fields_identity <- vector_result$identity
    vector_payload_indices <- vector_result$payload_indices
  } else if (!is.null(vector_field_default)) {
    stop(
      "vector_field_default requires vector_fields.",
      call. = FALSE
    )
  }

  # ===========================================================================
  # OBS EXPORT
  # ===========================================================================
  obs_manifest_path <- file.path(out_dir, "obs_manifest.json")
  obs_export <- .export_obs(
    obs = obs,
    obs_keys = obs_keys,
    latent = latent,
    embeddings = embeddings,
    obs_binary_dir = obs_binary_dir,
    obs_binary_dirname = obs_binary_dirname,
    centroid_outlier_quantile = centroid_outlier_quantile,
    centroid_min_points = centroid_min_points,
    obs_continuous_quantization = obs_continuous_quantization,
    obs_categorical_dtype = obs_categorical_dtype,
    compression = compression
  )

  obs_manifest_payload <- list(
    `_format` = manifest_format_version,
    n_points = as.integer(n_cells),
    centroid_outlier_quantile = if (!is.null(centroid_outlier_quantile)) {
      centroid_outlier_quantile
    } else {
      NULL
    },
    latent_key = "latent_space",
    compression = if (!is.null(compression)) compression else NULL,
    `_obsSchemas` = obs_export$schemas,
    `_continuousFields` = obs_export$continuous_fields,
    `_categoricalFields` = obs_export$categorical_fields
  )
  .write_json(obs_manifest_path, obs_manifest_payload, pretty = FALSE)

  identity_obs_fields <- .identity_obs_fields_from_compact_manifest(
    obs_manifest_payload
  )
  centroid_quantile_matches <- if (is.null(centroid_outlier_quantile)) {
    is.null(obs_manifest_payload$centroid_outlier_quantile)
  } else {
    is.numeric(obs_manifest_payload$centroid_outlier_quantile) &&
      length(obs_manifest_payload$centroid_outlier_quantile) == 1L &&
      !is.na(obs_manifest_payload$centroid_outlier_quantile) &&
      obs_manifest_payload$centroid_outlier_quantile ==
        centroid_outlier_quantile
  }
  compression_matches <- if (is.null(compression)) {
    is.null(obs_manifest_payload$compression)
  } else {
    is.numeric(obs_manifest_payload$compression) &&
      length(obs_manifest_payload$compression) == 1L &&
      !is.na(obs_manifest_payload$compression) &&
      obs_manifest_payload$compression == compression
  }
  if (
    !is.numeric(obs_manifest_payload$n_points) ||
      length(obs_manifest_payload$n_points) != 1L ||
      is.na(obs_manifest_payload$n_points) ||
      obs_manifest_payload$n_points != n_cells ||
      !centroid_quantile_matches ||
      !identical(obs_manifest_payload$latent_key, "latent_space") ||
      !compression_matches
  ) {
    stop(
      "Observation manifest ",
      obs_manifest_path,
      " does not match the staged export settings.",
      call. = FALSE
    )
  }
  ordered_obs_field_summaries <- c(
    Filter(
      function(field_info) identical(field_info$kind, "continuous"),
      obs_field_summaries
    ),
    Filter(
      function(field_info) identical(field_info$kind, "category"),
      obs_field_summaries
    )
  )
  expected_identity_obs_fields <- lapply(
    ordered_obs_field_summaries,
    function(field_info) {
      entry <- list(key = field_info$key, kind = field_info$kind)
      if (identical(field_info$kind, "category")) {
        entry$n_categories <- as.integer(field_info$category_count)
      }
      entry
    }
  )
  if (!identical(identity_obs_fields, expected_identity_obs_fields)) {
    stop(
      "Observation manifest ",
      obs_manifest_path,
      " fields do not match the staged observation fields.",
      call. = FALSE
    )
  }
  .require_declared_payloads_on_disk(
    out_dir,
    directory_name = obs_binary_dirname,
    declared = .declared_obs_payload_paths(obs_manifest_payload),
    axis = "Observation"
  )

  # ===========================================================================
  # VAR (GENE EXPRESSION) EXPORT
  # ===========================================================================
  genes_to_export <- character(0)
  if (!is.null(gene_expression)) {
    if (!is.data.frame(var)) {
      stop(
        "var must be a data.frame when gene_expression is supplied.",
        call. = FALSE
      )
    }

    gene_expression <- .as_matrix_like(gene_expression, name = "gene_expression")
    if (nrow(gene_expression) != n_cells) {
      stop("gene_expression has ", nrow(gene_expression), " cells, but embeddings have ", n_cells, " cells.")
    }
    if (nrow(var) != ncol(gene_expression)) {
      stop("var has ", nrow(var), " rows, but gene_expression has ", ncol(gene_expression), " genes.")
    }

    all_gene_ids <- .extract_gene_ids(var, var_gene_id_column)
    gene_id_to_idx <- stats::setNames(seq_along(all_gene_ids), all_gene_ids)

    if (is.null(gene_identifiers)) {
      genes_to_export <- all_gene_ids
    } else {
      .validate_character_vector(
        gene_identifiers,
        what = "gene_identifiers"
      )
      if (anyDuplicated(gene_identifiers)) {
        stop(
          "gene_identifiers must not contain duplicate identifiers.",
          call. = FALSE
        )
      }
      missing_genes <- setdiff(gene_identifiers, names(gene_id_to_idx))
      if (length(missing_genes) > 0) {
        stop(
          "gene_identifiers contains identifiers not found in var: ",
          paste(utils::head(missing_genes, 5), collapse = ", "),
          if (length(missing_genes) > 5) "..." else "",
          call. = FALSE
        )
      }
      genes_to_export <- gene_identifiers
    }
    # Uniqueness spans the whole var, because every row is addressable through
    # gene_identifiers=. Being drawable is a property of a name the viewer
    # shows, so it is checked on the exported genes: a var row left out reaches
    # no manifest, exactly as an obs column left out of obs_keys= does.
    .require_field_identities(genes_to_export, what = "Gene")

    var_manifest_path <- file.path(out_dir, "var_manifest.json")
    var_binary_dir <- file.path(out_dir, var_binary_dirname)
    .dir_create(var_binary_dir)

    var_manifest_fields <- vector("list", length(genes_to_export))

    for (idx in seq_along(genes_to_export)) {
      gene_id <- genes_to_export[[idx]]
      gene_idx <- unname(gene_id_to_idx[[gene_id]])
      payload_index <- .payload_index(idx)

      values <- .get_gene_column(gene_expression, gene_idx, n_cells = n_cells)

      if (!is.null(var_quantization)) {
        q <- .quantize_continuous(
          values,
          bits = var_quantization,
          field_name = gene_id
        )
        ext <- if (var_quantization == 8L) "u8" else "u16"
        value_path <- file.path(
          var_binary_dir,
          sprintf("%d.values.%s", payload_index, ext)
        )
        if (var_quantization == 8L) {
          .write_uint8(value_path, q$quantized, compression = compression)
        } else {
          .write_uint16(value_path, q$quantized, compression = compression)
        }
        var_manifest_fields[[idx]] <- list(
          payload_index,
          gene_id,
          q$min_val,
          q$max_val
        )
      } else {
        value_path <- file.path(
          var_binary_dir,
          sprintf("%d.values.f32", payload_index)
        )
        .write_float32_vector(value_path, values, compression = compression)
        var_manifest_fields[[idx]] <- list(payload_index, gene_id)
      }
    }

    gz_suffix <- if (!is.null(compression)) ".gz" else ""
    var_schema <- if (!is.null(var_quantization)) {
      ext <- if (var_quantization == 8L) "u8" else "u16"
      dtype_str <- if (var_quantization == 8L) "uint8" else "uint16"
      list(
        kind = "continuous",
        pathPattern = sprintf("%s/{index}.values.%s%s", var_binary_dirname, ext, gz_suffix),
        ext = ext,
        dtype = dtype_str,
        quantized = TRUE,
        quantizationBits = var_quantization
      )
    } else {
      list(
        kind = "continuous",
        pathPattern = sprintf("%s/{index}.values.f32%s", var_binary_dirname, gz_suffix),
        ext = "f32",
        dtype = "float32",
        quantized = FALSE
      )
    }

    var_manifest_payload <- list(
      `_format` = manifest_format_version,
      n_points = as.integer(n_cells),
      var_gene_id_column = var_gene_id_column,
      compression = if (!is.null(compression)) compression else NULL,
      quantization = if (!is.null(var_quantization)) var_quantization else NULL,
      `_varSchema` = var_schema,
      fields = var_manifest_fields
    )
    .write_json(var_manifest_path, var_manifest_payload, pretty = FALSE)

    exported_gene_names <- .gene_names_from_compact_manifest(
      var_manifest_payload
    )
    if (!identical(exported_gene_names, unname(genes_to_export))) {
      stop(
        "Gene manifest ",
        var_manifest_path,
        " names do not match the staged genes.",
        call. = FALSE
      )
    }
    .require_declared_payloads_on_disk(
      out_dir,
      directory_name = var_binary_dirname,
      declared = .declared_var_payload_paths(var_manifest_payload),
      axis = "Gene"
    )
  }

  # ===========================================================================
  # CONNECTIVITY EXPORT
  # ===========================================================================
  connectivity_meta <- list(
    n_edges = NULL,
    max_neighbors = NULL,
    index_dtype = NULL
  )

  if (!is.null(connectivities)) {
    if (is.null(connectivity_edges)) {
      stop("Validated connectivity edge pairs are unavailable.", call. = FALSE)
    }
    connectivity_manifest_path <- file.path(
      out_dir,
      "connectivity_manifest.json"
    )
    connectivity_meta$n_edges <- connectivity_edges$n_edges
    connectivity_meta$max_neighbors <- connectivity_edges$max_neighbors
    connectivity_meta$index_dtype <- connectivity_edges$index_dtype
    connectivity_binary_dir <- file.path(out_dir, connectivity_binary_dirname)
    .dir_create(connectivity_binary_dir)

    .write_connectivity_edges(
      connectivity_edges = connectivity_edges,
      out_dir = connectivity_binary_dir,
      compression = compression
    )

    sources_fname <- "edges.src.bin"
    dests_fname <- "edges.dst.bin"
    weights_fname <- "edges.weights.f64.bin"
    manifest_sources <- paste0(connectivity_binary_dirname, "/", sources_fname)
    manifest_dests <- paste0(connectivity_binary_dirname, "/", dests_fname)
    manifest_weights <- paste0(connectivity_binary_dirname, "/", weights_fname)
    if (!is.null(compression)) {
      manifest_sources <- paste0(manifest_sources, ".gz")
      manifest_dests <- paste0(manifest_dests, ".gz")
      manifest_weights <- paste0(manifest_weights, ".gz")
    }

    connectivity_manifest_payload <- list(
      format = "edge_pairs",
      n_cells = as.integer(n_cells),
      n_edges = connectivity_edges$n_edges,
      max_neighbors = connectivity_edges$max_neighbors,
      index_bytes = connectivity_edges$index_bytes,
      index_dtype = connectivity_edges$index_dtype,
      sourcesPath = manifest_sources,
      destinationsPath = manifest_dests,
      weightsPath = manifest_weights,
      weight_dtype = "float64",
      weight_bytes = 8L,
      compression = if (!is.null(compression)) compression else NULL
    )
    .write_json(connectivity_manifest_path, connectivity_manifest_payload, pretty = FALSE)

    declared_connectivity_paths <- character(0)
    for (path_key in c("sourcesPath", "destinationsPath", "weightsPath")) {
      declared_path <- connectivity_manifest_payload[[path_key]]
      if (
        !is.character(declared_path) ||
          length(declared_path) != 1L ||
          is.na(declared_path) ||
          !nzchar(declared_path)
      ) {
        stop(
          "connectivity_manifest.json ",
          path_key,
          " must be one payload path.",
          call. = FALSE
        )
      }
      declared_connectivity_paths <- c(
        declared_connectivity_paths,
        declared_path
      )
    }
    .require_declared_payloads_on_disk(
      out_dir,
      directory_name = connectivity_binary_dirname,
      declared = unique(declared_connectivity_paths),
      axis = "Connectivity"
    )
  }

  if (!is.null(vector_fields_identity)) {
    declared_vector_paths <- unlist(
      lapply(
        vector_fields_identity$fields,
        function(field_metadata) unlist(field_metadata$files, use.names = FALSE)
      ),
      use.names = FALSE
    )
    .require_dense_payload_indices(
      vector_payload_indices,
      axis = "Vector field"
    )
    .require_declared_payloads_on_disk(
      out_dir,
      directory_name = "vectors",
      declared = unique(declared_vector_paths),
      axis = "Vector field"
    )
  }

  # ===========================================================================
  # DATASET IDENTITY
  # ===========================================================================
  identity_path <- file.path(out_dir, "dataset_identity.json")
  cellucid_version <- as.character(utils::packageVersion("cellucid"))

  n_genes <- if (!is.null(gene_expression)) {
    as.integer(length(genes_to_export))
  } else {
    0L
  }

  n_obs_fields <- length(identity_obs_fields)
  n_categorical_fields <- sum(
    vapply(
      identity_obs_fields,
      function(field_info) identical(field_info$kind, "category"),
      logical(1)
    )
  )
  n_continuous_fields <- sum(
    vapply(
      identity_obs_fields,
      function(field_info) identical(field_info$kind, "continuous"),
      logical(1)
    )
  )

  export_settings <- list(
    compression = if (!is.null(compression)) compression else NULL,
    var_quantization = if (!is.null(var_quantization)) var_quantization else NULL,
    obs_continuous_quantization = if (!is.null(obs_continuous_quantization)) {
      obs_continuous_quantization
    } else {
      NULL
    },
    obs_categorical_dtype = obs_categorical_dtype
  )

  gz_suffix <- if (!is.null(compression)) ".gz" else ""
  embeddings_meta <- list(
    available_dimensions = I(available_dimensions),
    default_dimension = default_dimension,
    files = list()
  )
  for (dim in available_dimensions) {
    embeddings_meta$files[[sprintf("%dd", dim)]] <- sprintf("points_%dd.bin%s", dim, gz_suffix)
  }

  identity_payload <- list(
    version = 2,
    id = dataset_id,
    name = dataset_name,
    description = dataset_description,
    created_at = created_at,
    cellucid_data_version = cellucid_version,
    stats = list(
      n_cells = as.integer(n_cells),
      n_genes = as.integer(n_genes),
      n_obs_fields = as.integer(n_obs_fields),
      n_categorical_fields = as.integer(n_categorical_fields),
      n_continuous_fields = as.integer(n_continuous_fields),
      has_connectivity = !is.null(connectivity_meta$n_edges),
      n_edges = connectivity_meta$n_edges
    ),
    embeddings = embeddings_meta,
    obs_fields = identity_obs_fields,
    export_settings = export_settings
  )

  if (!is.null(source_info)) {
    identity_payload$source <- source_info
  }
  if (!is.null(vector_fields_identity)) {
    identity_payload$vector_fields <- vector_fields_identity
  }

  .write_json(identity_path, identity_payload, pretty = 2)
  .publish_export_generation(
    out_dir,
    final_out_dir,
    transaction_id = transaction$transaction_id,
    had_target = transaction$had_target
  )
  transaction_cleanup_armed <- FALSE

  invisible(NULL)
}

.validate_force <- function(force) {
  if (
    !is.logical(force) ||
      length(force) != 1L ||
      is.na(force) ||
      is.object(force)
  ) {
    stop("force must be exactly TRUE or FALSE.", call. = FALSE)
  }
  force
}

.validate_string <- function(value, arg) {
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value)
  ) {
    stop(arg, " must be exactly one string.", call. = FALSE)
  }
  value
}

.validate_required_string <- function(value, arg) {
  value <- .validate_string(value, arg)
  if (!nzchar(value) || !identical(value, trimws(value))) {
    stop(
      arg,
      " must be one non-empty string without leading or trailing whitespace.",
      call. = FALSE
    )
  }
  value
}

# Text the viewer prints verbatim -- a category label, a dataset name, a
# description, a source line -- is drawn with the browser's default white-space
# handling, in the legend, the field selector, and every exported figure.
# Characters that carry no glyph therefore change the value without changing
# the picture: "Liver " and "Liver" are two distinct categories that look
# identical on screen and in an exported SVG.
#
# All three classes below are rejected rather than repaired. Trimming would
# rewrite a scientific label the caller never asked to change, and would merge
# two distinct categories into one, silently moving cells between them. The
# Python exporter enforces the identical rule in _display_text_defect().
#
# R character vectors cannot hold U+0000, so the control class starts at U+0001.
.control_character_pattern <- "[\u0001-\u001f\u007f-\u009f]"
# Zero-width characters with no meaning of their own: ZERO WIDTH SPACE, WORD
# JOINER, and the byte-order mark a spreadsheet leaves at the front of a UTF-8
# CSV. U+200C and U+200D are deliberately absent: they join Indic, Persian, and
# emoji sequences, so banning them would reject real text.
.invisible_character_pattern <- "[\u200b\u2060\ufeff]"
# Every Unicode whitespace character that is not already a control character.
.edge_whitespace_class <- "[\u0020\u00a0\u1680\u2000-\u200a\u2028\u2029\u202f\u205f\u3000]"
.edge_whitespace_pattern <- paste0(
  "^",
  .edge_whitespace_class,
  "|",
  .edge_whitespace_class,
  "$"
)
.whitespace_run_pattern <- paste0(.edge_whitespace_class, "+")
# Every character the two non-control rejected classes can hold, so an
# invisible defect is named on screen exactly as the Python exporter names it.
.display_character_names <- c(
  "0020" = "U+0020 SPACE",
  "00A0" = "U+00A0 NO-BREAK SPACE",
  "1680" = "U+1680 OGHAM SPACE MARK",
  "2000" = "U+2000 EN QUAD",
  "2001" = "U+2001 EM QUAD",
  "2002" = "U+2002 EN SPACE",
  "2003" = "U+2003 EM SPACE",
  "2004" = "U+2004 THREE-PER-EM SPACE",
  "2005" = "U+2005 FOUR-PER-EM SPACE",
  "2006" = "U+2006 SIX-PER-EM SPACE",
  "2007" = "U+2007 FIGURE SPACE",
  "2008" = "U+2008 PUNCTUATION SPACE",
  "2009" = "U+2009 THIN SPACE",
  "200A" = "U+200A HAIR SPACE",
  "2028" = "U+2028 LINE SEPARATOR",
  "2029" = "U+2029 PARAGRAPH SEPARATOR",
  "202F" = "U+202F NARROW NO-BREAK SPACE",
  "205F" = "U+205F MEDIUM MATHEMATICAL SPACE",
  "3000" = "U+3000 IDEOGRAPHIC SPACE",
  "200B" = "U+200B ZERO WIDTH SPACE",
  "2060" = "U+2060 WORD JOINER",
  "FEFF" = "U+FEFF ZERO WIDTH NO-BREAK SPACE"
)

.describe_character <- function(character) {
  key <- sprintf("%04X", utf8ToInt(character))
  named <- unname(.display_character_names[key])
  if (!is.na(named)) {
    return(named)
  }
  # The C0 and C1 ranges carry no Unicode name at all, so say what they are.
  paste0("U+", key, " (control character)")
}

.escape_display_codepoint <- function(codepoint) {
  if (codepoint == 92L) {
    return("\\\\")
  }
  if (codepoint == 39L) {
    return("\\'")
  }
  if (codepoint == 9L) {
    return("\\t")
  }
  if (codepoint == 10L) {
    return("\\n")
  }
  if (codepoint == 13L) {
    return("\\r")
  }
  character <- intToUtf8(codepoint)
  if (
    codepoint == 32L ||
      !(
        grepl(.control_character_pattern, character, perl = TRUE) ||
          grepl(.invisible_character_pattern, character, perl = TRUE) ||
          grepl(.edge_whitespace_class, character, perl = TRUE)
      )
  ) {
    return(character)
  }
  if (codepoint <= 255L) {
    return(sprintf("\\x%02x", codepoint))
  }
  if (codepoint <= 65535L) {
    return(sprintf("\\u%04x", codepoint))
  }
  sprintf("\\U%08x", codepoint)
}

.escape_display_text <- function(text) {
  codepoints <- utf8ToInt(text)
  if (length(codepoints) == 0L) {
    return("''")
  }
  paste0(
    "'",
    paste(
      vapply(codepoints, .escape_display_codepoint, character(1)),
      collapse = ""
    ),
    "'"
  )
}

.display_text_defect <- function(text) {
  control <- regexpr(.control_character_pattern, text, perl = TRUE)
  if (control[[1L]] > 0L) {
    position <- control[[1L]]
    return(paste0(
      "contains ",
      .describe_character(substring(text, position, position))
    ))
  }
  hidden <- regexpr(.invisible_character_pattern, text, perl = TRUE)
  if (hidden[[1L]] > 0L) {
    position <- hidden[[1L]]
    return(paste0(
      "contains ",
      .describe_character(substring(text, position, position))
    ))
  }
  edge <- regexpr(.edge_whitespace_pattern, text, perl = TRUE)
  if (edge[[1L]] > 0L) {
    position <- edge[[1L]]
    return(paste0(
      if (position == 1L) "starts with" else "ends with",
      " ",
      .describe_character(substring(text, position, position))
    ))
  }
  NULL
}

.validate_display_text <- function(value, arg, allow_empty = FALSE) {
  value <- .validate_string(value, arg)
  if (!nzchar(value)) {
    if (allow_empty) {
      return(value)
    }
    stop(arg, " must be a non-empty string.", call. = FALSE)
  }
  defect <- .display_text_defect(value)
  if (!is.null(defect)) {
    stop(
      arg,
      " is displayed verbatim, so it must not carry characters that have no ",
      "glyph: ",
      .escape_display_text(value),
      " ",
      defect,
      ". Cellucid does not remove them for you, because that would change ",
      "text you did not ask to change. Pass the exact text you want shown.",
      call. = FALSE
    )
  }
  value
}

.validate_category_labels <- function(categories, key) {
  # Logical and numeric category labels carry no display text to inspect.
  if (!is.character(categories) || length(categories) == 0L) {
    return(invisible(categories))
  }
  suspect <- !nzchar(categories) |
    grepl(.control_character_pattern, categories, perl = TRUE) |
    grepl(.invisible_character_pattern, categories, perl = TRUE) |
    grepl(.edge_whitespace_pattern, categories, perl = TRUE)
  if (any(suspect)) {
    defects <- vapply(
      categories[suspect],
      function(label) {
        if (!nzchar(label)) {
          return("'' is empty")
        }
        paste0(.escape_display_text(label), " ", .display_text_defect(label))
      },
      character(1),
      USE.NAMES = FALSE
    )
    shown <- utils::head(defects, 10L)
    listed <- paste0("  - ", shown, collapse = "\n")
    remainder <- length(defects) - length(shown)
    if (remainder > 0L) {
      listed <- paste0(listed, "\n  - ... and ", remainder, " more")
    }
    counted <- if (length(defects) == 1L) {
      "1 label"
    } else {
      paste0(length(defects), " labels")
    }
    stop(
      "Categorical field '",
      key,
      "' has ",
      counted,
      " the viewer cannot show as written:\n",
      listed,
      "\nLabels are drawn verbatim in the legend, the field selector, and ",
      "exported figures, so these characters are invisible on screen while ",
      "still making the label a different value from the one it looks like. ",
      "Cellucid does not clean them for you: trimming a label rewrites your ",
      "annotation and can merge two categories into one, moving cells between ",
      "them. Clean the column and export again, for example:\n",
      "    obs[['",
      key,
      "']] <- factor(trimws(as.character(obs[['",
      key,
      "']]), whitespace = ",
      '"[\\\\h\\\\v]"',
      "))",
      call. = FALSE
    )
  }

  collapsed <- gsub(.whitespace_run_pattern, " ", categories, perl = TRUE)
  collision <- which(duplicated(collapsed))
  if (length(collision) > 0L) {
    later <- collision[[1L]]
    first <- match(collapsed[[later]], collapsed)
    stop(
      "Categorical field '",
      key,
      "' labels ",
      .escape_display_text(categories[[first]]),
      " and ",
      .escape_display_text(categories[[later]]),
      " are stored as different categories but are drawn identically: a run ",
      "of whitespace collapses to a single space in the legend, the field ",
      "selector, and exported figures. Rename one of them so the two ",
      "categories can be told apart, then export again.",
      call. = FALSE
    )
  }
  invisible(categories)
}

.validate_dataset_name <- function(value) {
  .validate_display_text(value, "dataset_name")
}


.validate_dataset_id <- function(value) {
  .safe_filename_component(value, what = "dataset_id")
}

.resolve_created_at <- function(value) {
  if (is.null(value)) {
    return(format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC",
      usetz = FALSE
    ))
  }
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value)
  ) {
    stop(
      "created_at must be NULL or exactly one string in ",
      "'YYYY-MM-DDTHH:MM:SSZ' UTC format.",
      call. = FALSE
    )
  }
  if (!grepl(
    "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$",
    value,
    perl = TRUE,
    useBytes = TRUE
  )) {
    stop(
      "created_at must use exact 'YYYY-MM-DDTHH:MM:SSZ' UTC format.",
      call. = FALSE
    )
  }

  year <- as.integer(substr(value, 1L, 4L))
  month <- as.integer(substr(value, 6L, 7L))
  day <- as.integer(substr(value, 9L, 10L))
  hour <- as.integer(substr(value, 12L, 13L))
  minute <- as.integer(substr(value, 15L, 16L))
  second <- as.integer(substr(value, 18L, 19L))
  leap_year <- year %% 400L == 0L ||
    (year %% 4L == 0L && year %% 100L != 0L)
  days_in_month <- c(
    31L,
    if (leap_year) 29L else 28L,
    31L,
    30L,
    31L,
    30L,
    31L,
    31L,
    30L,
    31L,
    30L,
    31L
  )
  if (
    year < 1L ||
      month < 1L ||
      month > 12L ||
      day < 1L ||
      day > days_in_month[[month]] ||
      hour < 0L ||
      hour > 23L ||
      minute < 0L ||
      minute > 59L ||
      second < 0L ||
      second > 59L
  ) {
    stop(
      "created_at must be a valid UTC calendar timestamp in ",
      "'YYYY-MM-DDTHH:MM:SSZ' format.",
      call. = FALSE
    )
  }
  value
}

.validate_source_identity <- function(name, url, citation) {
  if (!is.null(name)) {
    name <- .validate_display_text(name, "source_name")
  }
  if (!is.null(url)) {
    url <- .validate_display_text(url, "source_url")
  }
  if (!is.null(citation)) {
    citation <- .validate_display_text(
      citation,
      "source_citation"
    )
  }
  if (is.null(name) && is.null(url) && is.null(citation)) {
    return(NULL)
  }
  if (is.null(name)) {
    stop(
      "source_name is required whenever source_url or source_citation ",
      "is supplied.",
      call. = FALSE
    )
  }
  source <- list(name = name)
  if (!is.null(url)) {
    source$url <- url
  }
  if (!is.null(citation)) {
    source$citation <- citation
  }
  source
}

.validate_obs_categorical_dtype <- function(value) {
  supported <- c("uint8", "uint16")
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value) ||
      !value %in% supported
  ) {
    stop(
      "obs_categorical_dtype must be exactly one of \"uint8\" or \"uint16\".",
      call. = FALSE
    )
  }
  value
}

.normalize_compression <- function(compression) {
  if (is.null(compression)) {
    return(NULL)
  }
  if (
    !is.numeric(compression) ||
      length(compression) != 1L ||
      is.na(compression) ||
      !is.finite(compression) ||
      compression != floor(compression) ||
      compression < 1 ||
      compression > 9
  ) {
    stop(
      "compression must be NULL or one integer from 1 to 9.",
      call. = FALSE
    )
  }
  compression
}

.validate_quantization_bits <- function(bits, arg) {
  if (
    !is.numeric(bits) ||
      length(bits) != 1L ||
      is.na(bits) ||
      !is.finite(bits) ||
      !(bits %in% c(8, 16))
  ) {
    stop(arg, " must be exactly 8 or 16.", call. = FALSE)
  }
  bits
}

.validate_centroid_outlier_quantile <- function(value) {
  if (is.null(value)) {
    return(NULL)
  }
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value <= 0.5 ||
      value >= 1
  ) {
    stop(
      "centroid_outlier_quantile must be NULL or one finite number ",
      "strictly between 0.5 and 1.",
      call. = FALSE
    )
  }
  value
}

.validate_positive_integer <- function(value, arg) {
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value != floor(value) ||
      value < 1 ||
      value > .Machine$integer.max
  ) {
    stop(arg, " must be one positive integer.", call. = FALSE)
  }
  as.integer(value)
}

.validate_character_vector <- function(values, what) {
  if (!is.character(values) || !is.null(dim(values))) {
    stop(what, " must be a character vector.", call. = FALSE)
  }
  if (anyNA(values) || any(!nzchar(values))) {
    stop(
      what,
      " must contain only non-missing, non-empty strings.",
      call. = FALSE
    )
  }
  invisible(values)
}

.validate_output_path <- function(path) {
  path <- .validate_required_string(path, "out_dir")
  expanded <- path.expand(path)
  resolved <- normalizePath(
    expanded,
    winslash = "/",
    mustWork = FALSE
  )
  protected <- unique(c(
    "/",
    normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    normalizePath(path.expand("~"), winslash = "/", mustWork = TRUE)
  ))
  if (resolved %in% protected) {
    stop(
      "out_dir must name a dedicated dataset output directory.",
      call. = FALSE
    )
  }
  leaf <- basename(expanded)
  if (!nzchar(leaf) || leaf %in% c(".", "..")) {
    stop(
      "out_dir must name a dedicated dataset output directory.",
      call. = FALSE
    )
  }
  parent <- normalizePath(
    dirname(expanded),
    winslash = "/",
    mustWork = FALSE
  )
  final_path <- file.path(parent, leaf)
  if (final_path %in% protected) {
    stop(
      "out_dir must name a dedicated dataset output directory.",
      call. = FALSE
    )
  }
  final_path
}

.export_generation_lock_registry <- new.env(parent = emptyenv())
.export_generation_lock_registry$pid <- Sys.getpid()
.export_generation_lock_registry$paths <- new.env(parent = emptyenv())

.native_export_lock_acquire <- function(lock_path) {
  .Call(C_cellucid_export_lock_acquire, lock_path)
}

.native_export_lock_release <- function(lock_handle) {
  .Call(C_cellucid_export_lock_release, lock_handle)
}

.native_process_handle_count <- function() {
  .Call(C_cellucid_process_handle_count)
}

.native_export_path_info <- function(path) {
  .Call(C_cellucid_export_path_info, path)
}

.native_export_transaction_id <- function() {
  .Call(C_cellucid_export_transaction_id)
}

.native_export_write_journal <- function(path, contents) {
  .Call(C_cellucid_export_write_journal, path, contents)
}

.native_export_sync_directory <- function(path) {
  .Call(C_cellucid_export_sync_directory, path)
}

.native_export_lock_cleanup_status <- function(lock_handle) {
  status <- .native_export_lock_release(lock_handle)
  if (
    !is.integer(status) ||
      length(status) != 3L ||
      anyNA(status)
  ) {
    stop("Native export lock cleanup returned an invalid status.", call. = FALSE)
  }
  status
}

.assert_native_export_lock_cleanup <- function(status, lock_path) {
  if (status[[1L]] != 0L || status[[2L]] != 0L) {
    stop(
      "Could not completely release the export generation lock ",
      lock_path,
      " (unlock error ",
      status[[1L]],
      ", close error ",
      status[[2L]],
      ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.refresh_export_generation_lock_registry <- function() {
  process_id <- Sys.getpid()
  if (!identical(.export_generation_lock_registry$pid, process_id)) {
    inherited_registry <- .export_generation_lock_registry$paths
    inherited_keys <- ls(envir = inherited_registry, all.names = TRUE)
    .export_generation_lock_registry$pid <- process_id
    .export_generation_lock_registry$paths <- new.env(parent = emptyenv())
    cleanup_failures <- character()
    for (lock_key in inherited_keys) {
      record <- get(lock_key, envir = inherited_registry, inherits = FALSE)
      if (!is.null(record$handle)) {
        status <- tryCatch(
          .native_export_lock_cleanup_status(record$handle),
          error = identity
        )
        if (inherits(status, "error")) {
          cleanup_failures <- c(
            cleanup_failures,
            paste0(record$path, ": ", conditionMessage(status))
          )
        } else if (status[[1L]] != 0L || status[[2L]] != 0L) {
          cleanup_failures <- c(
            cleanup_failures,
            paste0(
              record$path,
              ": unlock error ",
              status[[1L]],
              ", close error ",
              status[[2L]]
            )
          )
        }
      }
    }
    if (length(cleanup_failures) != 0L) {
      stop(
        "Could not discard every fork-inherited export lock handle: ",
        paste(cleanup_failures, collapse = "; "),
        call. = FALSE
      )
    }
  }
  process_id
}

.export_generation_lock_path <- function(final_path) {
  parent <- normalizePath(
    dirname(final_path),
    winslash = "/",
    mustWork = TRUE
  )
  file.path(
    parent,
    paste0(".", basename(final_path), ".cellucid.lock")
  )
}

.validate_export_generation_lock_path <- function(lock_path) {
  path_info <- .export_path_info(lock_path)
  if (path_info$kind == 3L) {
    stop(
      "Export lock path must not be a symbolic link or reparse point: ",
      lock_path,
      call. = FALSE
    )
  }
  if (.Platform$OS.type != "windows") {
    link_target <- suppressWarnings(Sys.readlink(lock_path))
    if (
      length(link_target) == 1L &&
        !is.na(link_target) &&
        nzchar(link_target)
    ) {
      stop(
        "Export lock path must not be a symbolic link: ",
        lock_path,
        call. = FALSE
      )
    }
  }
  if (
    path_info$kind != 0L &&
      (path_info$kind != 1L || !utils::file_test("-f", lock_path))
  ) {
    stop(
      "Export lock path must identify a regular file: ",
      lock_path,
      call. = FALSE
    )
  }
  invisible(lock_path)
}

.canonical_export_generation_lock_key <- function(lock_path) {
  key <- normalizePath(lock_path, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows") {
    key <- tolower(key)
  }
  key
}

.acquire_export_generation_lock <- function(final_path) {
  lock_path <- .export_generation_lock_path(final_path)
  .validate_export_generation_lock_path(lock_path)
  lock_key <- .canonical_export_generation_lock_key(lock_path)
  process_id <- .refresh_export_generation_lock_registry()
  registry <- .export_generation_lock_registry$paths

  if (exists(lock_key, envir = registry, inherits = FALSE)) {
    stop(
      "An export generation is already active for ",
      final_path,
      ".",
      call. = FALSE
    )
  }

  assign(
    lock_key,
    list(handle = NULL, path = lock_path, pid = process_id),
    envir = registry
  )
  reservation_is_owned <- TRUE
  on.exit({
    if (
      reservation_is_owned &&
        exists(lock_key, envir = registry, inherits = FALSE)
    ) {
      rm(list = lock_key, envir = registry)
    }
  }, add = TRUE)

  lock <- .native_export_lock_acquire(lock_path)
  if (is.null(lock)) {
    stop(
      "An export generation is already active for ",
      final_path,
      ".",
      call. = FALSE
    )
  }

  export_lock <- list(
    handle = lock,
    key = lock_key,
    path = lock_path,
    pid = process_id
  )
  assign(lock_key, export_lock, envir = registry)
  reservation_is_owned <- FALSE
  export_lock
}

.release_export_generation_lock <- function(export_lock) {
  process_id <- .refresh_export_generation_lock_registry()
  if (!identical(export_lock$pid, process_id)) {
    return(invisible(FALSE))
  }

  registry <- .export_generation_lock_registry$paths
  if (!exists(export_lock$key, envir = registry, inherits = FALSE)) {
    stop(
      "Export generation lock ownership is not active: ",
      export_lock$path,
      call. = FALSE
    )
  }
  active_lock <- get(export_lock$key, envir = registry, inherits = FALSE)
  if (!identical(active_lock$handle, export_lock$handle)) {
    stop(
      "Export generation lock ownership does not match: ",
      export_lock$path,
      call. = FALSE
    )
  }

  status <- .native_export_lock_cleanup_status(export_lock$handle)
  if (status[[3L]] == 1L) {
    rm(list = export_lock$key, envir = registry)
  }
  .assert_native_export_lock_cleanup(status, export_lock$path)
  invisible(TRUE)
}

.export_transaction_format <- "cellucid-export-transaction"
.export_transaction_version <- 1L
.export_transaction_id_pattern <- "^[0-9a-f]{32}$"

.validate_export_transaction_id <- function(transaction_id) {
  if (
    !is.character(transaction_id) ||
      length(transaction_id) != 1L ||
      is.na(transaction_id) ||
      is.object(transaction_id) ||
      !grepl(
        .export_transaction_id_pattern,
        transaction_id,
        perl = TRUE,
        useBytes = TRUE
      )
  ) {
    stop("Export transaction identity is not canonical.", call. = FALSE)
  }
  transaction_id
}

.export_transaction_paths <- function(final_path, transaction_id) {
  transaction_id <- .validate_export_transaction_id(transaction_id)
  parent <- dirname(final_path)
  stem <- paste0(".", basename(final_path))
  list(
    journal = file.path(
      parent,
      paste0(stem, ".cellucid-transaction.json")
    ),
    journal_temp = file.path(
      parent,
      paste0(stem, ".cellucid-transaction.json.tmp")
    ),
    stage = file.path(
      parent,
      paste0(stem, ".cellucid-stage-", transaction_id)
    ),
    backup = file.path(
      parent,
      paste0(stem, ".cellucid-backup-", transaction_id)
    )
  )
}

.export_transaction_control_paths <- function(final_path) {
  paths <- .export_transaction_paths(final_path, strrep("0", 32L))
  paths[c("journal", "journal_temp")]
}

.export_path_info <- function(path) {
  info <- .native_export_path_info(path)
  if (
    !is.double(info) ||
      length(info) != 2L ||
      anyNA(info) ||
      info[[1L]] < 0 ||
      info[[1L]] > 4 ||
      info[[1L]] != floor(info[[1L]]) ||
      info[[2L]] < 0 ||
      info[[2L]] != floor(info[[2L]])
  ) {
    stop("Native export path inspection returned invalid state.", call. = FALSE)
  }
  list(
    kind = as.integer(info[[1L]]),
    links = info[[2L]]
  )
}

.require_export_directory_or_absent <- function(path, label) {
  info <- .export_path_info(path)
  if (info$kind == 0L) {
    return(FALSE)
  }
  if (info$kind != 2L) {
    stop(
      label,
      " must be an ordinary non-symbolic directory or absent: ",
      path,
      call. = FALSE
    )
  }
  TRUE
}

.require_export_regular_file <- function(path, label) {
  info <- .export_path_info(path)
  if (info$kind == 0L) {
    stop(label, " is missing: ", path, call. = FALSE)
  }
  if (info$kind != 1L || info$links != 1) {
    stop(
      label,
      " must be one non-linked, non-symbolic regular file: ",
      path,
      call. = FALSE
    )
  }
  invisible(info)
}

.fsync_export_directory <- function(path) {
  status <- .native_export_sync_directory(path)
  if (
    !is.logical(status) ||
      length(status) != 1L ||
      is.na(status)
  ) {
    stop(
      "Native export directory synchronization returned invalid state.",
      call. = FALSE
    )
  }
  invisible(status)
}

.rename_export_path <- function(source, destination) {
  if (!isTRUE(file.rename(source, destination))) {
    stop(
      "Could not rename export transaction path ",
      source,
      " to ",
      destination,
      ".",
      call. = FALSE
    )
  }
  .fsync_export_directory(dirname(destination))
  invisible(destination)
}

.remove_export_tree <- function(path) {
  if (!.require_export_directory_or_absent(
    path,
    "Export transaction directory"
  )) {
    return(invisible(path))
  }
  status <- unlink(path, recursive = TRUE, force = TRUE)
  if (status != 0L || .export_path_info(path)$kind != 0L) {
    stop(
      "Could not remove export transaction directory: ",
      path,
      call. = FALSE
    )
  }
  .fsync_export_directory(dirname(path))
  invisible(path)
}

.remove_export_control_file <- function(path) {
  .require_export_regular_file(
    path,
    "Export transaction control file"
  )
  status <- unlink(path, recursive = FALSE, force = TRUE)
  if (status != 0L || .export_path_info(path)$kind != 0L) {
    stop(
      "Could not remove export transaction control file: ",
      path,
      call. = FALSE
    )
  }
  .fsync_export_directory(dirname(path))
  invisible(path)
}

.serialize_export_transaction <- function(transaction_id, had_target) {
  transaction_id <- .validate_export_transaction_id(transaction_id)
  if (
    !is.logical(had_target) ||
      length(had_target) != 1L ||
      is.na(had_target) ||
      is.object(had_target)
  ) {
    stop(
      "Export transaction had_target must be exactly boolean.",
      call. = FALSE
    )
  }
  paste0(
    "{\"format\":\"",
    .export_transaction_format,
    "\",\"version\":",
    .export_transaction_version,
    ",\"transaction_id\":\"",
    transaction_id,
    "\",\"had_target\":",
    if (had_target) "true" else "false",
    "}\n"
  )
}

.read_export_transaction <- function(journal_path) {
  .require_export_regular_file(
    journal_path,
    "Export transaction journal"
  )
  journal_bytes <- readBin(
    journal_path,
    what = "raw",
    n = 513L
  )
  if (length(journal_bytes) > 512L) {
    stop(
      "Export transaction journal is unexpectedly large: ",
      journal_path,
      call. = FALSE
    )
  }
  if (
    length(journal_bytes) == 0L ||
      any(as.integer(journal_bytes) > 127L)
  ) {
    stop(
      "Export transaction journal is malformed: ",
      journal_path,
      call. = FALSE
    )
  }
  journal <- tryCatch(
    rawToChar(journal_bytes),
    error = function(error) NULL
  )
  if (is.null(journal)) {
    stop(
      "Export transaction journal is malformed: ",
      journal_path,
      call. = FALSE
    )
  }

  prefix <- paste0(
    "{\"format\":\"",
    .export_transaction_format,
    "\",\"version\":",
    .export_transaction_version,
    ",\"transaction_id\":\""
  )
  for (had_target in c(TRUE, FALSE)) {
    suffix <- paste0(
      "\",\"had_target\":",
      if (had_target) "true" else "false",
      "}\n"
    )
    expected_length <- nchar(prefix, type = "bytes") +
      32L +
      nchar(suffix, type = "bytes")
    if (
      nchar(journal, type = "bytes") == expected_length &&
        startsWith(journal, prefix) &&
        endsWith(journal, suffix)
    ) {
      transaction_id <- substr(
        journal,
        nchar(prefix, type = "chars") + 1L,
        nchar(prefix, type = "chars") + 32L
      )
      if (
        grepl(
          .export_transaction_id_pattern,
          transaction_id,
          perl = TRUE,
          useBytes = TRUE
        ) &&
          identical(
            journal_bytes,
            charToRaw(.serialize_export_transaction(
              transaction_id,
              had_target
            ))
          )
      ) {
        return(list(
          transaction_id = transaction_id,
          had_target = had_target
        ))
      }
    }
  }
  stop(
    "Export transaction journal is not canonical: ",
    journal_path,
    call. = FALSE
  )
}

.require_active_export_transaction <- function(
    journal_path,
    transaction_id,
    had_target
) {
  transaction <- .read_export_transaction(journal_path)
  if (
    !identical(transaction$transaction_id, transaction_id) ||
      !identical(transaction$had_target, had_target)
  ) {
    stop(
      "Export transaction journal does not describe the active transaction.",
      call. = FALSE
    )
  }
  invisible(transaction)
}

.discard_export_transaction_temp <- function(journal_temp) {
  if (.export_path_info(journal_temp)$kind == 0L) {
    return(invisible(journal_temp))
  }
  .remove_export_control_file(journal_temp)
}

.recover_export_transaction <- function(final_path) {
  controls <- .export_transaction_control_paths(final_path)
  .discard_export_transaction_temp(controls$journal_temp)
  if (.export_path_info(controls$journal)$kind == 0L) {
    return(invisible(final_path))
  }

  transaction <- .read_export_transaction(controls$journal)
  paths <- .export_transaction_paths(
    final_path,
    transaction$transaction_id
  )
  target_exists <- .require_export_directory_or_absent(
    final_path,
    "Export target"
  )
  stage_exists <- .require_export_directory_or_absent(
    paths$stage,
    "Staged export generation"
  )
  backup_exists <- .require_export_directory_or_absent(
    paths$backup,
    "Prior export generation"
  )
  state <- c(target_exists, stage_exists, backup_exists)

  if (transaction$had_target) {
    if (identical(state, c(TRUE, FALSE, FALSE))) {
      NULL
    } else if (identical(state, c(TRUE, TRUE, FALSE))) {
      .remove_export_tree(paths$stage)
    } else if (identical(state, c(FALSE, TRUE, TRUE))) {
      .rename_export_path(paths$backup, final_path)
      .remove_export_tree(paths$stage)
    } else if (identical(state, c(TRUE, FALSE, TRUE))) {
      .remove_export_tree(paths$backup)
    } else {
      stop(
        "Export transaction cannot be recovered without guessing whether ",
        "to commit or roll back: target/stage/backup state is ",
        paste(state, collapse = "/"),
        ".",
        call. = FALSE
      )
    }
  } else {
    if (identical(state, c(FALSE, FALSE, FALSE))) {
      NULL
    } else if (identical(state, c(FALSE, TRUE, FALSE))) {
      .remove_export_tree(paths$stage)
    } else if (identical(state, c(TRUE, FALSE, FALSE))) {
      NULL
    } else {
      stop(
        "Initial export transaction cannot be recovered without guessing ",
        "whether to commit or roll back: target/stage/backup state is ",
        paste(state, collapse = "/"),
        ".",
        call. = FALSE
      )
    }
  }

  .require_active_export_transaction(
    paths$journal,
    transaction$transaction_id,
    transaction$had_target
  )
  .remove_export_control_file(paths$journal)
  invisible(final_path)
}

.new_export_transaction_id <- function(final_path) {
  for (attempt in seq_len(128L)) {
    transaction_id <- .native_export_transaction_id()
    transaction_id <- .validate_export_transaction_id(transaction_id)
    paths <- .export_transaction_paths(final_path, transaction_id)
    if (
      .export_path_info(paths$stage)$kind == 0L &&
        .export_path_info(paths$backup)$kind == 0L
    ) {
      return(transaction_id)
    }
  }
  stop(
    "Could not allocate a unique export transaction identity.",
    call. = FALSE
  )
}

.write_export_transaction <- function(
    final_path,
    transaction_id,
    had_target
) {
  paths <- .export_transaction_paths(final_path, transaction_id)
  if (
    .export_path_info(paths$journal)$kind != 0L ||
      .export_path_info(paths$journal_temp)$kind != 0L
  ) {
    stop(
      "Export transaction control path is already occupied for ",
      final_path,
      ".",
      call. = FALSE
    )
  }
  journal_bytes <- charToRaw(.serialize_export_transaction(
    transaction_id,
    had_target
  ))
  status <- .native_export_write_journal(
    paths$journal_temp,
    journal_bytes
  )
  if (!identical(status, TRUE)) {
    stop(
      "Native export transaction journal write returned invalid state.",
      call. = FALSE
    )
  }
  .rename_export_path(paths$journal_temp, paths$journal)
  invisible(paths$journal)
}

.begin_export_transaction <- function(final_path, force) {
  had_target <- .require_export_directory_or_absent(
    final_path,
    "Export target"
  )
  if (had_target && !force) {
    stop(
      "out_dir already exists; set force = TRUE to replace the complete ",
      "generation.",
      call. = FALSE
    )
  }
  transaction_id <- .new_export_transaction_id(final_path)
  .write_export_transaction(
    final_path,
    transaction_id,
    had_target
  )
  paths <- .export_transaction_paths(final_path, transaction_id)
  if (
    !dir.create(
      paths$stage,
      recursive = FALSE,
      showWarnings = FALSE
    )
  ) {
    stop("Could not create the staged output directory.", call. = FALSE)
  }
  .fsync_export_directory(dirname(final_path))
  list(
    transaction_id = transaction_id,
    had_target = had_target,
    stage = paths$stage,
    backup = paths$backup
  )
}

.publish_export_generation <- function(
    stage,
    final_path,
    transaction_id,
    had_target
) {
  paths <- .export_transaction_paths(final_path, transaction_id)
  if (!identical(stage, paths$stage)) {
    stop(
      "Staged export path does not belong to the active transaction.",
      call. = FALSE
    )
  }
  .require_active_export_transaction(
    paths$journal,
    transaction_id,
    had_target
  )
  if (.export_path_info(paths$journal_temp)$kind != 0L) {
    stop(
      "Export transaction journal temporary path reappeared during publication.",
      call. = FALSE
    )
  }
  if (!.require_export_directory_or_absent(
    stage,
    "Staged export generation"
  )) {
    stop("Staged export directory is missing: ", stage, call. = FALSE)
  }

  if (had_target) {
    if (!.require_export_directory_or_absent(
      final_path,
      "Export target"
    )) {
      stop(
        "Prior export generation is missing: ",
        final_path,
        call. = FALSE
      )
    }
    if (.export_path_info(paths$backup)$kind != 0L) {
      stop(
        "Export backup path is already occupied: ",
        paths$backup,
        call. = FALSE
      )
    }
    .rename_export_path(final_path, paths$backup)
  } else if (.export_path_info(final_path)$kind != 0L) {
    stop(
      "Initial export target appeared during publication: ",
      final_path,
      call. = FALSE
    )
  }

  .rename_export_path(stage, final_path)
  if (had_target) {
    .remove_export_tree(paths$backup)
  }
  .require_active_export_transaction(
    paths$journal,
    transaction_id,
    had_target
  )
  .remove_export_control_file(paths$journal)
  invisible(final_path)
}

.dir_create <- function(path) {
  if (!dir.exists(path)) {
    created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!created && !dir.exists(path)) {
      stop("Could not create directory: ", path, call. = FALSE)
    }
  }
  invisible(path)
}

# dataset_id is the one identifier that still names a path: the export
# directory a producer publishes and every data path addresses it by. Payload
# names on the obs, var, and vectors axes are integer indices, so no field
# identifier reaches this rule.
.safe_filename_component <- function(name, what = "Filename component") {
  name <- .validate_required_string(name, what)
  if (
    !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", name) ||
      grepl("\\.$", name) ||
      nchar(name, type = "bytes") > 180L
  ) {
    stop(
      what, " '", name, "' is not a portable identifier. Use ",
      "1-180 ASCII letters, numbers, '.', '_', or '-', beginning with a ",
      "letter or number and not ending with '.'.",
      call. = FALSE
    )
  }
  device_base <- toupper(strsplit(name, ".", fixed = TRUE)[[1]][[1]])
  if (
    device_base %in% c("CON", "PRN", "AUX", "NUL") ||
      grepl("^(COM|LPT)[1-9]$", device_base)
  ) {
    stop(
      what, " '", name, "' is reserved on Windows.",
      call. = FALSE
    )
  }
  name
}

# A payload filename is an integer index, so an identifier is never a path and
# carries no filename rule at all: no ASCII restriction, no case-collision
# rule, no reserved Windows name. What survives is what the identity is
# actually for. It names one field in the manifest, so it must be a non-empty
# string, distinct within its axis so a lookup resolves one field, and text the
# viewer can draw exactly as it is stored -- the same rule every category label
# obeys, because a gene name and a category label are shown in the same legend.
# The cellucid Python package enforces the identical rule.
#
# `what` names one field in the singular -- "Gene", "Observation field",
# "Vector field" -- because it is read as the subject of a sentence about a
# single identifier: "Gene identifier at position 1 ...". The checks that speak
# about the whole axis add their own plural, so no caller has to guess which
# number a message will need. cellucid-python passes the same singular nouns to
# _require_field_identities().
.require_field_identities <- function(values, what) {
  .validate_character_vector(values, what = paste0(what, " keys"))
  for (position in seq_along(values)) {
    .validate_display_text(
      values[[position]],
      paste0(what, " identifier at position ", position - 1L),
      allow_empty = FALSE
    )
  }
  .require_unique_identifiers(values, what = what)
}

# Identity, not path. A repeated key names no single field, so a lookup by key
# resolves to more than one row whatever the payloads are called. The reported
# key is the first repeat in supplied order, which is the one
# cellucid-python's _require_unique_identifiers() reports.
.require_unique_identifiers <- function(keys, what) {
  .validate_character_vector(keys, what = paste0(what, " keys"))
  repeated <- which(duplicated(keys))
  if (length(repeated) > 0L) {
    stop(
      what,
      " key '",
      keys[[repeated[[1L]]]],
      "' is duplicated.",
      call. = FALSE
    )
  }
  invisible(keys)
}

# A message that shows the reader a set of values shows it as a list, not as a
# sentence that happens to continue after a comma: without a boundary,
# "columns not in obs: a, b. Available columns: x" cannot be read back as two
# lists. cellucid-python interpolates a Python list, which prints as [0, 1] for
# integers and ['a', 'b'] for strings, and the two writers report one defect to
# one reader, so they print it the same way.
.format_value_list <- function(values) {
  if (is.character(values)) {
    items <- sprintf("'%s'", values)
  } else {
    items <- format(values, trim = TRUE)
  }
  paste0("[", paste(items, collapse = ", "), "]")
}

# A payload filename is the position of its field on its axis, never the
# field's own name, so no exported path depends on a dataset's vocabulary. The
# manifest entry carries the same integer as its first element, and the viewer
# substitutes it into the schema path pattern to reach the bytes.
.payload_index <- function(position) {
  if (
    !is.numeric(position) ||
      length(position) != 1L ||
      is.na(position) ||
      !is.finite(position) ||
      position != floor(position) ||
      position < 1 ||
      position > .Machine$integer.max
  ) {
    stop("Payload position must be one positive integer.", call. = FALSE)
  }
  as.integer(position) - 1L
}

# Within one axis directory the indices must be exactly 0 to N-1, each used
# once. The index *is* the filename, so two fields holding one index write into
# one file: the second overwrites the first and the viewer then draws one
# field's values under the other field's name. Nothing downstream can detect
# that, so it is asserted here in the writer, against the manifest that was
# just built, rather than only in a test.
.require_dense_payload_indices <- function(indices, axis) {
  for (position in seq_along(indices)) {
    index <- indices[[position]]
    if (
      !is.integer(index) ||
        length(index) != 1L ||
        is.na(index)
    ) {
      stop(
        axis,
        " payload index at position ",
        position - 1L,
        " must be a native integer.",
        call. = FALSE
      )
    }
  }
  resolved <- as.integer(unlist(indices, use.names = FALSE))
  if (!identical(sort(resolved), seq_along(resolved) - 1L)) {
    stop(
      axis,
      " payload indices must be exactly 0..",
      length(resolved) - 1L,
      ", each used once; got ",
      .format_value_list(sort(resolved)),
      ".",
      call. = FALSE
    )
  }
  invisible(resolved)
}

.expand_payload_pattern <- function(pattern, index, label, ext = NULL) {
  if (
    !is.character(pattern) ||
      length(pattern) != 1L ||
      is.na(pattern) ||
      !nzchar(pattern)
  ) {
    stop(label, " must be a non-empty path pattern.", call. = FALSE)
  }
  expanded <- gsub("{index}", as.character(index), pattern, fixed = TRUE)
  if (!is.null(ext)) {
    expanded <- gsub("{ext}", ext, expanded, fixed = TRUE)
  }
  if (grepl("[{}]", expanded)) {
    stop(
      label,
      " retains an unsubstituted placeholder: '",
      expanded,
      "'.",
      call. = FALSE
    )
  }
  expanded
}

# The manifest is the only index the viewer has: a payload it does not declare
# is invisible, and a payload it declares but that was never written fails the
# dataset at read time, in the browser, long after the export succeeded. Both
# are caught here by re-expanding the emitted path patterns and comparing them
# against the directory that was actually written.
.require_declared_payloads_on_disk <- function(
    out_dir,
    directory_name,
    declared,
    axis
) {
  directory <- file.path(out_dir, directory_name)
  on_disk <- character(0)
  if (dir.exists(directory)) {
    for (entry in sort(list.files(directory, all.files = TRUE, no.. = TRUE))) {
      if (!utils::file_test("-f", file.path(directory, entry))) {
        stop(
          axis,
          " payload directory holds a non-file entry: ",
          file.path(directory, entry),
          call. = FALSE
        )
      }
      on_disk <- c(on_disk, paste0(directory_name, "/", entry))
    }
  }
  missing_payloads <- sort(setdiff(declared, on_disk))
  undeclared <- sort(setdiff(on_disk, declared))
  if (length(missing_payloads) > 0L || length(undeclared) > 0L) {
    stop(
      axis,
      " manifest does not describe the payloads that were written. ",
      "Declared but absent: ",
      paste(missing_payloads, collapse = ", "),
      ". Written but undeclared: ",
      paste(undeclared, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(declared)
}

.codes_extension_by_dtype <- c(uint8 = "u8", uint16 = "u16")

.declared_obs_payload_paths <- function(manifest) {
  schemas <- manifest[["_obsSchemas"]]
  declared <- character(0)
  for (field in manifest[["_continuousFields"]]) {
    declared <- c(declared, .expand_payload_pattern(
      schemas$continuous$pathPattern,
      index = field[[1L]],
      label = "obs continuous pathPattern"
    ))
  }
  for (field in manifest[["_categoricalFields"]]) {
    extension <- unname(.codes_extension_by_dtype[field[[4L]]])
    if (is.na(extension)) {
      stop(
        "Categorical obs field '",
        field[[2L]],
        "' declares an unknown codes dtype '",
        field[[4L]],
        "'.",
        call. = FALSE
      )
    }
    declared <- c(
      declared,
      .expand_payload_pattern(
        schemas$categorical$codesPathPattern,
        index = field[[1L]],
        label = "obs categorical codesPathPattern",
        ext = extension
      ),
      .expand_payload_pattern(
        schemas$categorical$outlierPathPattern,
        index = field[[1L]],
        label = "obs categorical outlierPathPattern"
      )
    )
  }
  unique(declared)
}

.declared_var_payload_paths <- function(manifest) {
  schema <- manifest[["_varSchema"]]
  declared <- vapply(
    manifest$fields,
    function(field) {
      .expand_payload_pattern(
        schema$pathPattern,
        index = field[[1L]],
        label = "var pathPattern"
      )
    },
    character(1)
  )
  unique(declared)
}

.gene_names_from_compact_manifest <- function(manifest) {
  fields <- manifest$fields
  if (!is.list(fields)) {
    stop("compact_v1 var manifest fields must be a list.", call. = FALSE)
  }
  gene_names <- character(0)
  payload_indices <- list()
  for (field in fields) {
    if (
      !is.list(field) ||
        !(length(field) %in% c(2L, 4L)) ||
        !is.character(field[[2L]]) ||
        length(field[[2L]]) != 1L ||
        is.na(field[[2L]]) ||
        !nzchar(field[[2L]])
    ) {
      stop(
        "compact_v1 var fields must be exact [index, name] or ",
        "[index, name, minValue, maxValue] tuples.",
        call. = FALSE
      )
    }
    payload_indices[[length(payload_indices) + 1L]] <- field[[1L]]
    gene_names <- c(gene_names, field[[2L]])
  }
  .require_unique_identifiers(gene_names, what = "Gene")
  .require_dense_payload_indices(payload_indices, axis = "Gene")
  gene_names
}

.as_dense_matrix <- function(x, name) {
  if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
  } else if (
    is.data.frame(x) &&
      all(vapply(x, function(column) {
        is.numeric(column) && !is.object(column)
      }, logical(1)))
  ) {
    x <- as.matrix(x)
  } else if (
    !is.matrix(x) ||
      !is.numeric(x) ||
      is.object(x)
  ) {
    stop(name, " must be a finite real numeric matrix.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  x
}

.as_matrix_like <- function(x, name) {
  if (inherits(x, "Matrix")) {
    if (!inherits(x, "dMatrix")) {
      stop(name, " must contain real numeric values.", call. = FALSE)
    }
    values <- x@x
    if (anyNA(values) || any(!is.finite(values))) {
      stop(
        name,
        " must not contain infinite or missing sparse values.",
        call. = FALSE
      )
    }
    return(x)
  } else if (
    is.data.frame(x) &&
      all(vapply(x, function(column) {
        is.numeric(column) && !is.object(column)
      }, logical(1)))
  ) {
    x <- as.matrix(x)
  } else if (
    !is.matrix(x) ||
      !is.numeric(x) ||
      is.object(x)
  ) {
    stop(name, " must be a real numeric matrix.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.collect_embeddings <- function(embeddings_in) {
  embeddings <- list()
  for (dim_chr in names(embeddings_in)) {
    x <- embeddings_in[[dim_chr]]
    if (is.null(x)) {
      next
    }
    dim_int <- as.integer(dim_chr)
    if (
      is.numeric(x) &&
        !is.object(x) &&
        is.null(dim(x)) &&
        dim_int == 1L
    ) {
      x <- matrix(x, ncol = 1L)
    } else if (
      is.data.frame(x) &&
        all(vapply(x, function(column) {
          is.numeric(column) && !is.object(column)
        }, logical(1)))
    ) {
      x <- as.matrix(x)
    } else if (inherits(x, "Matrix")) {
      x <- as.matrix(x)
    } else if (
      !is.matrix(x) ||
        !is.numeric(x) ||
        is.object(x)
    ) {
      stop(
        "X_umap_", dim_chr,
        "d must be a finite real numeric matrix.",
        call. = FALSE
      )
    } else {
      x <- unclass(x)
    }
    storage.mode(x) <- "double"
    embeddings[[dim_chr]] <- x
  }
  if (length(embeddings) == 0) {
    stop("At least one embedding must be provided: X_umap_1d, X_umap_2d, or X_umap_3d.")
  }
  embeddings
}

.validate_embeddings <- function(embeddings) {
  n_cells <- NULL
  for (dim_chr in names(embeddings)) {
    dim_int <- as.integer(dim_chr)
    arr <- embeddings[[dim_chr]]
    if (length(dim(arr)) != 2L) {
      stop("X_umap_", dim_int, "d must be a 2D array.")
    }
    if (ncol(arr) != dim_int) {
      stop("X_umap_", dim_int, "d must have exactly ", dim_int, " columns, got ", ncol(arr), ".")
    }
    if (is.null(n_cells)) {
      n_cells <- nrow(arr)
    } else if (nrow(arr) != n_cells) {
      stop("All embeddings must have the same number of cells (rows).")
    }
  }
  if (
    is.null(n_cells) ||
      !is.numeric(n_cells) ||
      length(n_cells) != 1L ||
      n_cells < 1L ||
      n_cells > .Machine$integer.max
  ) {
    stop(
      "Embeddings must contain between 1 and 2,147,483,647 cells.",
      call. = FALSE
    )
  }
  as.integer(n_cells)
}

.normalize_embedding <- function(arr) {
  if (
    !is.matrix(arr) ||
      !is.numeric(arr) ||
      anyNA(arr) ||
      any(!is.finite(arr))
  ) {
    stop(
      "Embedding coordinates must be a finite numeric matrix.",
      call. = FALSE
    )
  }
  axis_mins <- apply(arr, 2, min)
  axis_maxs <- apply(arr, 2, max)
  axis_ranges <- axis_maxs - axis_mins
  max_range <- max(axis_ranges)
  if (!is.finite(max_range) || max_range <= 0) {
    stop(
      "Embedding coordinates must span a positive finite coordinate range.",
      call. = FALSE
    )
  }
  center <- (axis_mins + axis_maxs) / 2
  scale_factor <- 2 / max_range

  centered <- sweep(arr, 2, center, FUN = "-")
  coords <- centered * scale_factor

  list(
    coords = coords,
    info = list(
      original_range = as.numeric(max_range),
      center = as.numeric(center),
      scale_factor = as.numeric(scale_factor)
    )
  )
}

.write_float32_matrix_row_major <- function(path, mat, compression = NULL) {
  .write_float32_vector(path, as.vector(t(mat)), compression = compression)
}

.write_float32_vector <- function(
    path,
    values,
    compression = NULL,
    allow_generated_nan = FALSE
) {
  if (
    !is.logical(allow_generated_nan) ||
      length(allow_generated_nan) != 1L ||
      is.na(allow_generated_nan)
  ) {
    stop("allow_generated_nan must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.numeric(values) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop(
      "float32 values must be a native numeric vector.",
      call. = FALSE
    )
  }
  has_forbidden_missing <- any(is.na(values) & !is.nan(values))
  if (
    has_forbidden_missing ||
      any(is.infinite(values)) ||
      (!allow_generated_nan && any(is.nan(values)))
  ) {
    stop(
      "float32 values must contain only finite float32 values; generated ",
      "NaN is accepted only when allow_generated_nan is TRUE.",
      call. = FALSE
    )
  }

  float32_max <- (2 - 2^-23) * 2^127
  float32_min_subnormal <- 2^-149
  outside_finite_float32 <- (
    abs(values) > float32_max |
      (values != 0 & abs(values) < float32_min_subnormal)
  )
  if (any(outside_finite_float32, na.rm = TRUE)) {
    stop(
      "Every nonzero finite float32 value must be representable between ",
      "2^-149 and (2 - 2^-23) * 2^127 in absolute magnitude.",
      call. = FALSE
    )
  }

  con <- .binary_connection(path, compression = compression)
  on.exit(close(con), add = TRUE)
  writeBin(as.double(values), con, size = 4L, endian = "little")
  invisible(.final_binary_path(path, compression))
}

.write_float64_vector <- function(path, values, compression = NULL) {
  if (
    !is.double(values) ||
      is.object(values) ||
      !is.null(dim(values)) ||
      anyNA(values) ||
      any(!is.finite(values))
  ) {
    stop(
      "float64 values must be a finite native double vector.",
      call. = FALSE
    )
  }
  con <- .binary_connection(path, compression = compression)
  on.exit(close(con), add = TRUE)
  writeBin(values, con, size = 8L, endian = "little")
  invisible(.final_binary_path(path, compression))
}

.write_uint8 <- function(path, values, compression = NULL) {
  if (
    !(is.integer(values) || is.double(values)) ||
      is.object(values) ||
      !is.null(dim(values)) ||
      anyNA(values) ||
      any(!is.finite(values)) ||
      any(values != floor(values)) ||
      any(values < 0 | values > 255)
  ) {
    stop(
      "uint8 values must be a numeric vector of finite integers in the ",
      "inclusive range 0 to 255.",
      call. = FALSE
    )
  }
  con <- .binary_connection(path, compression = compression)
  on.exit(close(con), add = TRUE)
  raw_vals <- as.raw(as.integer(values))
  writeBin(raw_vals, con, endian = "little")
  invisible(.final_binary_path(path, compression))
}

.write_uint16 <- function(path, values, compression = NULL) {
  con <- .binary_connection(path, compression = compression)
  on.exit(close(con), add = TRUE)
  raw_vals <- .pack_uint16(values)
  writeBin(raw_vals, con, endian = "little")
  invisible(.final_binary_path(path, compression))
}

.write_uint32 <- function(path, values, compression = NULL) {
  con <- .binary_connection(path, compression = compression)
  on.exit(close(con), add = TRUE)
  raw_vals <- .pack_uint32(values)
  writeBin(raw_vals, con, endian = "little")
  invisible(.final_binary_path(path, compression))
}

.final_binary_path <- function(path, compression) {
  if (!is.null(compression)) paste0(path, ".gz") else path
}

.binary_connection <- function(path, compression = NULL) {
  final <- .final_binary_path(path, compression)
  if (!is.null(compression)) {
    return(gzfile(final, open = "wb", compression = compression))
  }
  file(final, open = "wb")
}

.identity_obs_fields_from_compact_manifest <- function(manifest) {
  expected_names <- c(
    "_format",
    "n_points",
    "centroid_outlier_quantile",
    "latent_key",
    "compression",
    "_obsSchemas",
    "_continuousFields",
    "_categoricalFields"
  )
  if (
    !is.list(manifest) ||
      is.null(names(manifest)) ||
      length(names(manifest)) != length(expected_names) ||
      !setequal(names(manifest), expected_names)
  ) {
    stop(
      "obs manifest must contain exactly the current compact_v1 fields.",
      call. = FALSE
    )
  }
  if (!identical(manifest[["_format"]], "compact_v1")) {
    stop("obs manifest must use the current compact_v1 format.", call. = FALSE)
  }

  continuous_fields <- manifest[["_continuousFields"]]
  categorical_fields <- manifest[["_categoricalFields"]]
  if (!is.list(continuous_fields)) {
    stop(
      "compact_v1 obs manifest _continuousFields must be a list.",
      call. = FALSE
    )
  }
  if (!is.list(categorical_fields)) {
    stop(
      "compact_v1 obs manifest _categoricalFields must be a list.",
      call. = FALSE
    )
  }

  identity_fields <- list()
  manifest_keys <- character(0)
  # Both arrays write into obs/, so their payload indices share one space.
  payload_indices <- list()
  for (field in continuous_fields) {
    if (
      !is.list(field) ||
        !(length(field) %in% c(2L, 4L)) ||
        !is.character(field[[2]]) ||
        length(field[[2]]) != 1L ||
        is.na(field[[2]]) ||
        !nzchar(field[[2]])
    ) {
      stop(
        "compact_v1 continuous observation fields must be exact ",
        "[index, key] or [index, key, minValue, maxValue] tuples.",
        call. = FALSE
      )
    }
    payload_indices[[length(payload_indices) + 1L]] <- field[[1]]
    key <- field[[2]]
    manifest_keys <- c(manifest_keys, key)
    identity_fields[[length(identity_fields) + 1L]] <- list(
      key = key,
      kind = "continuous"
    )
  }

  for (field in categorical_fields) {
    if (
      !is.list(field) ||
        !(length(field) %in% c(6L, 8L)) ||
        !is.character(field[[2]]) ||
        length(field[[2]]) != 1L ||
        is.na(field[[2]]) ||
        !nzchar(field[[2]]) ||
        !(is.atomic(field[[3]]) || is.list(field[[3]])) ||
        !is.null(dim(field[[3]]))
    ) {
      stop(
        "compact_v1 categorical observation fields must be exact ",
        "six- or eight-member tuples with a category array.",
        call. = FALSE
      )
    }
    payload_indices[[length(payload_indices) + 1L]] <- field[[1]]
    key <- field[[2]]
    manifest_keys <- c(manifest_keys, key)
    identity_fields[[length(identity_fields) + 1L]] <- list(
      key = key,
      kind = "category",
      n_categories = as.integer(length(field[[3]]))
    )
  }

  .require_unique_identifiers(
    manifest_keys,
    what = "Observation field"
  )
  .require_dense_payload_indices(payload_indices, axis = "Observation")
  identity_fields
}

.write_json <- function(path, payload, pretty = FALSE) {
  json <- jsonlite::toJSON(
    payload,
    auto_unbox = TRUE,
    null = "null",
    digits = 17
  )

  if (isTRUE(pretty)) {
    json <- jsonlite::prettify(json)
  } else if (is.numeric(pretty) && length(pretty) == 1L) {
    json <- jsonlite::prettify(json, indent = as.integer(pretty))
  }

  writeLines(json, con = path, useBytes = TRUE)
  invisible(path)
}

.float32_quantization_chunk_size <- 1048576L

.roundtrip_finite_float32_chunk <- function(values, field_name) {
  encoded <- writeBin(
    as.double(values),
    raw(),
    size = 4L,
    endian = "little"
  )
  rounded <- readBin(
    encoded,
    what = double(),
    n = length(values),
    size = 4L,
    endian = "little"
  )
  if (
    length(rounded) != length(values) ||
      anyNA(rounded) ||
      any(!is.finite(rounded))
  ) {
    stop(
      "Cannot quantize field '",
      field_name,
      "': values must remain finite in the viewer's float32 domain.",
      call. = FALSE
    )
  }
  rounded
}

# Name compact_v1's constant-field case from one payload's bounds.
#
# A gene detected in no published cell, or an obs column a subset flattened, is
# ordinary scientific data, so the format encodes it rather than rejecting it.
# The case is declared by equal bounds -- the entry stays exactly
# [index, key, minValue, maxValue], so neither its shape nor its length moves --
# and every code is written as 0. Writer and reader then take the same named
# branch instead of the general arithmetic: the writer never divides by
# (maxValue - minValue), and the reader fills minValue straight through. The
# value returns bit-exact, not within a quantization step.
#
# Sole derivation point on this writer, and the exact peer of
# _is_constant_continuous_range() in the Python exporter.
.is_constant_continuous_range <- function(min_val, max_val) {
  min_val == max_val
}

.quantize_continuous <- function(
    values,
    bits,
    field_name,
    .chunk_size = .float32_quantization_chunk_size
) {
  bits <- .validate_quantization_bits(bits, arg = "bits")
  field_name <- .validate_required_string(field_name, "field_name")
  .chunk_size <- .validate_positive_integer(
    .chunk_size,
    ".chunk_size"
  )
  if (
    !is.numeric(values) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop(
      "Cannot quantize field '",
      field_name,
      "': values must be a native numeric vector.",
      call. = FALSE
    )
  }
  if (length(values) == 0L) {
    stop(
      "Cannot quantize field '",
      field_name,
      "': it contains no finite values.",
      call. = FALSE
    )
  }

  min_val <- Inf
  max_val <- -Inf
  starts <- seq.int(1L, length(values), by = .chunk_size)
  for (start in starts) {
    end <- min(
      length(values),
      start + as.double(.chunk_size) - 1
    )
    source_chunk <- values[start:end]
    if (anyNA(source_chunk) || any(!is.finite(source_chunk))) {
      stop(
        "Cannot quantize field '",
        field_name,
        "': values must contain only finite values.",
        call. = FALSE
      )
    }
    rounded <- .roundtrip_finite_float32_chunk(
      source_chunk,
      field_name
    )
    min_val <- min(min_val, min(rounded))
    max_val <- max(max_val, max(rounded))
  }

  if (bits == 8L) {
    max_quant <- 254L
  } else {
    max_quant <- 65534L
  }

  if (.is_constant_continuous_range(min_val, max_val)) {
    # The constant case, taken before any scale exists. The general path
    # divides by (max_val - min_val), which is exactly zero here, so a constant
    # field must never reach it. This also covers native-double variation that
    # collapses to one value in the viewer's float32 domain: one float32 value
    # is one constant. Mirrors _quantize_continuous() in the Python exporter.
    return(list(
      quantized = integer(length(values)),
      min_val = min_val,
      max_val = max_val,
      scale = 0
    ))
  }

  scale <- max_quant / (max_val - min_val)
  quantized <- integer(length(values))
  for (start in starts) {
    end <- min(
      length(values),
      start + as.double(.chunk_size) - 1
    )
    rounded <- .roundtrip_finite_float32_chunk(
      values[start:end],
      field_name
    )
    normalized <- (rounded - min_val) * scale
    normalized[normalized < 0] <- 0
    normalized[normalized > max_quant] <- max_quant
    quantized[start:end] <- as.integer(
      floor(normalized + 1e-8)
    )
  }

  list(
    quantized = quantized,
    min_val = min_val,
    max_val = max_val,
    scale = scale
  )
}

.quantize_nullable_outlier_quantiles <- function(values, bits, field_name) {
  bits <- .validate_quantization_bits(bits, arg = "bits")
  field_name <- .validate_required_string(field_name, "field_name")
  if (
    !is.numeric(values) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop(
      "Cannot quantize outlier field '",
      field_name,
      "': values must be a native numeric vector.",
      call. = FALSE
    )
  }
  if (any(is.na(values) & !is.nan(values))) {
    stop(
      "Cannot quantize outlier field '",
      field_name,
      "': only generated NaN values may represent missing quantiles.",
      call. = FALSE
    )
  }
  if (any(is.infinite(values))) {
    stop(
      "Cannot quantize outlier field '",
      field_name,
      "': values must contain only finite values or generated NaN.",
      call. = FALSE
    )
  }

  missing_mask <- is.nan(values)
  finite_quantized <- .quantize_continuous(
    values[!missing_mask],
    bits = bits,
    field_name = field_name
  )
  missing_value <- if (bits == 8L) 255L else 65535L
  quantized <- rep(missing_value, length(values))
  quantized[!missing_mask] <- finite_quantized$quantized

  list(
    quantized = quantized,
    min_val = finite_quantized$min_val,
    max_val = finite_quantized$max_val,
    scale = finite_quantized$scale
  )
}

.category_values_and_codes <- function(values, key) {
  if (is.logical(values) && !is.object(values)) {
    present <- !is.na(values)
    categories <- c(FALSE, TRUE)
    categories <- categories[vapply(
      categories,
      function(category) any(present & values == category),
      logical(1)
    )]
    codes <- match(values, categories) - 1L
    codes[!present] <- -1L
    return(list(categories = categories, codes = codes))
  }

  if (is.factor(values)) {
    categorical <- values
  } else if (is.character(values) && !is.object(values)) {
    categorical <- factor(values)
  } else {
    stop(
      "Categorical observations must be a factor, logical, or character ",
      "vector.",
      call. = FALSE
    )
  }
  codes <- as.integer(categorical) - 1L
  codes[is.na(codes)] <- -1L
  categories <- as.character(levels(categorical))
  .validate_category_labels(categories, key)
  list(
    categories = categories,
    codes = codes
  )
}

.observation_kind <- function(values, key) {
  if (is.factor(values)) {
    return("category")
  }
  if (
    (is.logical(values) || is.character(values)) &&
      !is.object(values)
  ) {
    return("category")
  }
  if (is.numeric(values) && !is.object(values)) {
    return("continuous")
  }
  stop(
    "Observation field '", key,
    "' must be a native numeric, logical, character, or factor vector.",
    call. = FALSE
  )
}

.summarize_obs_field <- function(s, key, obs_continuous_quantization, obs_categorical_dtype) {
  kind <- .observation_kind(s, key)

  if (identical(kind, "continuous")) {
    dtype_str <- if (is.null(obs_continuous_quantization)) {
      "float32"
    } else if (obs_continuous_quantization == 8L) {
      "uint8"
    } else {
      "uint16"
    }
    return(
      list(
        key = key,
        kind = "continuous",
        quantized = !is.null(obs_continuous_quantization),
        quantization_bits = if (!is.null(obs_continuous_quantization)) {
          obs_continuous_quantization
        } else {
          NULL
        },
        dtype = dtype_str
      )
    )
  }

  categories <- .category_values_and_codes(s, key)$categories
  n_categories <- length(categories)
  dtype_str <- if (identical(obs_categorical_dtype, "uint8")) {
    if (n_categories > 255L) {
      stop(
        "Field '", key, "' has ", n_categories,
        " categories, but uint8 supports at most 255."
      )
    }
    "uint8"
  } else {
    if (n_categories > 65535L) {
      stop(
        "Field '", key, "' has ", n_categories,
        " categories, but uint16 supports at most 65,535."
      )
    }
    "uint16"
  }

  list(
    key = key,
    kind = "category",
    category_count = as.integer(n_categories),
    codes_dtype = dtype_str,
    outlier_quantized = !is.null(obs_continuous_quantization),
    outlier_quantization_bits = if (!is.null(obs_continuous_quantization)) {
      obs_continuous_quantization
    } else {
      NULL
    }
  )
}

.export_obs <- function(
    obs,
    obs_keys,
    latent,
    embeddings,
    obs_binary_dir,
    obs_binary_dirname,
    centroid_outlier_quantile,
    centroid_min_points,
    obs_continuous_quantization,
    obs_categorical_dtype,
    compression
) {
  continuous_fields <- list()
  categorical_fields <- list()
  continuous_dtype_info <- list()
  categorical_dtype_info <- list()

  for (idx in seq_along(obs_keys)) {
    key <- obs_keys[[idx]]
    s <- obs[[key]]
    # Continuous and categorical payloads share one obs/ directory, so they
    # share one index space: the position of the key in obs_keys.
    payload_index <- .payload_index(idx)

    kind <- .observation_kind(s, key)

    if (identical(kind, "continuous")) {
      values <- s
      storage.mode(values) <- "double"
      if (length(values) != nrow(obs)) {
        stop("Length mismatch for obs[['", key, "']]: ", length(values), " vs ", nrow(obs))
      }
      if (anyNA(values) || any(!is.finite(values))) {
        stop(
          "Observation field '", key,
          "' must contain only finite values.",
          call. = FALSE
        )
      }

      if (!is.null(obs_continuous_quantization)) {
        q <- .quantize_continuous(values, bits = obs_continuous_quantization, field_name = key)
        ext <- if (obs_continuous_quantization == 8L) "u8" else "u16"
        dtype_str <- if (obs_continuous_quantization == 8L) "uint8" else "uint16"
        value_path <- file.path(
          obs_binary_dir,
          sprintf("%d.values.%s", payload_index, ext)
        )
        if (obs_continuous_quantization == 8L) {
          .write_uint8(value_path, q$quantized, compression = compression)
        } else {
          .write_uint16(value_path, q$quantized, compression = compression)
        }
        continuous_fields[[length(continuous_fields) + 1L]] <- list(
          payload_index,
          key,
          q$min_val,
          q$max_val
        )
        if (length(continuous_dtype_info) == 0L) {
          continuous_dtype_info <- list(
            ext = ext,
            dtype = dtype_str,
            quantized = TRUE,
            quantizationBits = obs_continuous_quantization
          )
        }
      } else {
        value_path <- file.path(
          obs_binary_dir,
          sprintf("%d.values.f32", payload_index)
        )
        .write_float32_vector(value_path, values, compression = compression)
        continuous_fields[[length(continuous_fields) + 1L]] <- list(
          payload_index,
          key
        )
        if (length(continuous_dtype_info) == 0L) {
          continuous_dtype_info <- list(ext = "f32", dtype = "float32", quantized = FALSE)
        }
      }

      next
    }

    categorical <- .category_values_and_codes(s, key)
    categories <- categorical$categories
    codes <- categorical$codes

    n_categories <- length(categories)
    dtype_choice <- if (identical(obs_categorical_dtype, "uint8")) {
      if (n_categories > 255L) {
        stop(
          "Field '", key, "' has ", n_categories,
          " categories, but uint8 supports at most 255."
        )
      }
      list(dtype = "uint8", missing_value = 255L)
    } else {
      if (n_categories > 65535L) {
        stop(
          "Field '", key, "' has ", n_categories,
          " categories, but uint16 supports at most 65,535."
        )
      }
      list(dtype = "uint16", missing_value = 65535L)
    }

    missing_value <- dtype_choice$missing_value
    dtype_str <- dtype_choice$dtype

    codes_typed <- rep(missing_value, length(codes))
    valid_mask <- codes >= 0L
    codes_typed[valid_mask] <- codes[valid_mask]

    if (identical(dtype_str, "uint8")) {
      codes_fname <- sprintf("%d.codes.u8", payload_index)
      codes_path <- file.path(obs_binary_dir, codes_fname)
      .write_uint8(codes_path, codes_typed, compression = compression)
    } else {
      codes_fname <- sprintf("%d.codes.u16", payload_index)
      codes_path <- file.path(obs_binary_dir, codes_fname)
      .write_uint16(codes_path, codes_typed, compression = compression)
    }

    if (is.null(centroid_outlier_quantile)) {
      centroids_by_dim <- stats::setNames(
        lapply(names(embeddings), function(...) list()),
        names(embeddings)
      )
    } else {
      centroids_by_dim <- .compute_centroids_for_all_dimensions(
        embeddings = embeddings,
        codes = codes,
        categories = categories,
        outlier_quantile = centroid_outlier_quantile,
        min_points = centroid_min_points
      )
    }

    outlier_quantiles <- .compute_latent_space_quantiles(
      latent = latent,
      codes = codes,
      categories = categories,
      min_points = centroid_min_points
    )

    if (!is.null(obs_continuous_quantization)) {
      oq <- .quantize_nullable_outlier_quantiles(
        outlier_quantiles,
        bits = obs_continuous_quantization,
        field_name = paste0(key, "_outliers")
      )
      oq_ext <- if (obs_continuous_quantization == 8L) "u8" else "u16"
      oq_dtype_str <- if (obs_continuous_quantization == 8L) "uint8" else "uint16"

      outlier_fname <- sprintf("%d.outliers.%s", payload_index, oq_ext)
      outlier_path <- file.path(obs_binary_dir, outlier_fname)
      if (obs_continuous_quantization == 8L) {
        .write_uint8(outlier_path, oq$quantized, compression = compression)
      } else {
        .write_uint16(outlier_path, oq$quantized, compression = compression)
      }

      categorical_fields[[length(categorical_fields) + 1L]] <- list(
        payload_index,
        key,
        I(categories),
        dtype_str,
        as.integer(missing_value),
        centroids_by_dim,
        oq$min_val,
        oq$max_val
      )

      if (length(categorical_dtype_info) == 0L) {
        categorical_dtype_info <- list(
          codesExt = if (identical(dtype_str, "uint8")) "u8" else "u16",
          outlierExt = oq_ext,
          outlierDtype = oq_dtype_str,
          outlierQuantized = TRUE
        )
      }
    } else {
      outlier_fname <- sprintf("%d.outliers.f32", payload_index)
      outlier_path <- file.path(obs_binary_dir, outlier_fname)
      .write_float32_vector(
        outlier_path,
        outlier_quantiles,
        compression = compression,
        allow_generated_nan = TRUE
      )

      categorical_fields[[length(categorical_fields) + 1L]] <- list(
        payload_index,
        key,
        I(categories),
        dtype_str,
        as.integer(missing_value),
        centroids_by_dim
      )

      if (length(categorical_dtype_info) == 0L) {
        categorical_dtype_info <- list(
          codesExt = if (identical(dtype_str, "uint8")) "u8" else "u16",
          outlierExt = "f32",
          outlierDtype = "float32",
          outlierQuantized = FALSE
        )
      }
    }
  }

  gz_suffix <- if (!is.null(compression)) ".gz" else ""
  schemas <- list()
  if (length(continuous_dtype_info) > 0L) {
    cont <- list(
      pathPattern = sprintf(
        "%s/{index}.values.%s%s",
        obs_binary_dirname,
        continuous_dtype_info$ext,
        gz_suffix
      ),
      ext = continuous_dtype_info$ext,
      dtype = continuous_dtype_info$dtype,
      quantized = isTRUE(continuous_dtype_info$quantized)
    )
    if (isTRUE(continuous_dtype_info$quantized)) {
      cont$quantizationBits <- as.integer(continuous_dtype_info$quantizationBits)
    }
    schemas$continuous <- cont
  }

  if (length(categorical_dtype_info) > 0L) {
    schemas$categorical <- list(
      codesPathPattern = sprintf("%s/{index}.codes.{ext}%s", obs_binary_dirname, gz_suffix),
      outlierPathPattern = sprintf(
        "%s/{index}.outliers.%s%s",
        obs_binary_dirname,
        categorical_dtype_info$outlierExt,
        gz_suffix
      ),
      outlierExt = categorical_dtype_info$outlierExt,
      outlierDtype = categorical_dtype_info$outlierDtype,
      outlierQuantized = isTRUE(categorical_dtype_info$outlierQuantized)
    )
  }

  list(
    schemas = schemas,
    continuous_fields = continuous_fields,
    categorical_fields = categorical_fields
  )
}

.compute_centroids_for_field <- function(coords, codes, categories, outlier_quantile = 0.95, min_points = 10) {
  if (nrow(coords) != length(codes)) {
    stop("coords and codes must have the same length.")
  }

  outlier_quantile <- .validate_centroid_outlier_quantile(
    outlier_quantile
  )
  if (is.null(outlier_quantile)) {
    stop(
      "outlier_quantile must be a finite number when computing centroids.",
      call. = FALSE
    )
  }

  centroids <- list()
  for (code in seq_along(categories) - 1L) {
    label <- categories[[code + 1L]]
    idx <- which(codes == code)
    n <- length(idx)
    if (n < min_points) {
      next
    }

    pts <- coords[idx, , drop = FALSE]
    center <- colMeans(pts)

    if (n > min_points) {
      diffs <- sweep(pts, 2, center, FUN = "-")
      dists <- sqrt(rowSums(diffs * diffs))
      thr <- as.numeric(stats::quantile(dists, probs = outlier_quantile, names = FALSE, type = 7))
      inlier_mask <- dists <= thr
      n_in <- sum(inlier_mask)
      if (n_in >= min_points) {
        pts_in <- pts[inlier_mask, , drop = FALSE]
        center <- colMeans(pts_in)
        used_count <- n_in
      } else {
        used_count <- n
      }
    } else {
      used_count <- n
    }

    centroids[[length(centroids) + 1L]] <- list(
      category = label,
      position = I(as.numeric(center)),
      n_points = as.integer(used_count)
    )
  }

  centroids
}

.compute_centroids_for_all_dimensions <- function(embeddings, codes, categories, outlier_quantile = 0.95, min_points = 10) {
  result <- list()
  for (dim in names(embeddings)) {
    coords <- embeddings[[dim]]
    result[[dim]] <- .compute_centroids_for_field(
      coords = coords,
      codes = codes,
      categories = categories,
      outlier_quantile = outlier_quantile,
      min_points = min_points
    )
  }
  result
}

.compute_latent_space_quantiles <- function(latent, codes, categories, min_points = 10) {
  n_cells <- nrow(latent)
  quantiles <- rep(NaN, n_cells)

  for (code in seq_along(categories) - 1L) {
    idx <- which(codes == code)
    n <- length(idx)
    if (n < min_points) {
      next
    }
    pts <- latent[idx, , drop = FALSE]
    centroid <- colMeans(pts)
    diffs <- sweep(pts, 2, centroid, FUN = "-")
    dists <- sqrt(rowSums(diffs * diffs))
    ranks <- rank(dists, ties.method = "max")
    quantiles[idx] <- as.numeric(ranks) / n
  }

  quantiles
}

.extract_gene_ids <- function(var, var_gene_id_column) {
  if (
    !is.null(var_gene_id_column) &&
      (
        !is.character(var_gene_id_column) ||
          length(var_gene_id_column) != 1L ||
          is.na(var_gene_id_column) ||
          !nzchar(var_gene_id_column)
      )
  ) {
    stop(
      "var_gene_id_column must be NULL or one non-empty string.",
      call. = FALSE
    )
  }

  if (is.null(var_gene_id_column)) {
    ids <- rownames(var)
    if (is.null(ids)) {
      stop(
        "var has no row names; provide an explicit character gene ID column.",
        call. = FALSE
      )
    }
    .validate_gene_ids(ids)
    return(ids)
  }
  if (!var_gene_id_column %in% names(var)) {
    stop(
      "var_gene_id_column '", var_gene_id_column, "' not found in var. Available columns: ",
      paste(names(var), collapse = ", ")
    )
  }
  ids <- var[[var_gene_id_column]]
  .validate_gene_ids(ids)
  ids
}

.validate_gene_ids <- function(gene_ids) {
  .require_unique_identifiers(gene_ids, what = "Gene")
}

.get_gene_column <- function(gene_expression, gene_idx, n_cells) {
  if (inherits(gene_expression, "dgCMatrix")) {
    p <- gene_expression@p
    i <- gene_expression@i
    x <- gene_expression@x

    start <- p[[gene_idx]] + 1L
    end <- p[[gene_idx + 1L]]
    values <- numeric(n_cells)
    if (end >= start) {
      rows <- i[start:end] + 1L
      values[rows] <- x[start:end]
    }
    return(values)
  }

  if (inherits(gene_expression, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Matrix package is required to export sparse gene_expression objects.")
    }
    dense <- as.matrix(gene_expression[, gene_idx, drop = FALSE])
    return(as.numeric(dense[, 1L]))
  }

  values <- gene_expression[, gene_idx]
  as.numeric(values)
}

.pack_uint16 <- function(values) {
  if (
    !(is.integer(values) || is.double(values)) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop("uint16 values must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(values) || any(!is.finite(values))) {
    stop("uint16 values must all be finite.", call. = FALSE)
  }
  if (any(values != floor(values))) {
    stop("uint16 values must all be integers.", call. = FALSE)
  }
  if (any(values < 0 | values > 65535)) {
    stop(
      "uint16 values must all be in the inclusive range 0 to 65,535.",
      call. = FALSE
    )
  }
  integer_values <- as.integer(values)
  lo <- as.raw(bitwAnd(integer_values, 0xFFL))
  hi <- as.raw(bitwAnd(bitwShiftR(integer_values, 8L), 0xFFL))
  as.raw(rbind(lo, hi))
}

.pack_uint32 <- function(values) {
  if (
    !(is.integer(values) || is.double(values)) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop("uint32 values must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(values) || any(!is.finite(values))) {
    stop("uint32 values must all be finite.", call. = FALSE)
  }
  if (any(values != floor(values))) {
    stop("uint32 values must all be integers.", call. = FALSE)
  }
  if (any(values < 0 | values > 4294967295)) {
    stop(
      "uint32 values must all be in the inclusive range 0 to 4,294,967,295.",
      call. = FALSE
    )
  }
  b0 <- as.raw(floor(values) %% 256)
  b1 <- as.raw(floor(values / 256) %% 256)
  b2 <- as.raw(floor(values / 65536) %% 256)
  b3 <- as.raw(floor(values / 16777216) %% 256)
  as.raw(rbind(b0, b1, b2, b3))
}

.connectivity_triplet_values <- function(triplet) {
  if ("x" %in% names(triplet)) {
    return(triplet$x)
  }
  rep.int(1, nrow(triplet))
}

.connectivity_triplet <- function(matrix) {
  triplet_matrix <- if (inherits(matrix, "TsparseMatrix")) {
    matrix
  } else {
    methods::as(matrix, "TsparseMatrix")
  }
  rows <- triplet_matrix@i
  columns <- triplet_matrix@j
  if (length(rows) > 1L) {
    coordinate_order <- order(rows, columns, method = "radix")
    ordered_rows <- rows[coordinate_order]
    ordered_columns <- columns[coordinate_order]
    duplicate_offsets <- which(
      ordered_rows[-1L] == ordered_rows[-length(ordered_rows)] &
        ordered_columns[-1L] == ordered_columns[-length(ordered_columns)]
    )
    if (length(duplicate_offsets) > 0L) {
      duplicate_offset <- duplicate_offsets[[1L]]
      stop(
        "Connectivity contains duplicate sparse coordinates at (",
        ordered_rows[[duplicate_offset]], ", ",
        ordered_columns[[duplicate_offset]], ").",
        call. = FALSE
      )
    }
  }
  Matrix::summary(triplet_matrix)
}

.validate_connectivity_edges <- function(connectivities, n_cells) {
  info <- .connectivity_index_dtype(n_cells)
  expected_shape <- c(n_cells, n_cells)

  if (inherits(connectivities, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Matrix package is required to validate sparse connectivity matrices.")
    }
    if (!(
      inherits(connectivities, "dMatrix") ||
        inherits(connectivities, "lMatrix") ||
        inherits(connectivities, "nMatrix")
    )) {
      stop(
        "Connectivity values must be real numeric or logical values.",
        call. = FALSE
      )
    }
    if (!identical(as.numeric(dim(connectivities)), as.numeric(expected_shape))) {
      stop(
        "Connectivity matrix shape must exactly match the cell axis: expected (",
        n_cells, ", ", n_cells, "), got (",
        nrow(connectivities), ", ", ncol(connectivities), ").",
        call. = FALSE
      )
    }

    triplet <- .connectivity_triplet(connectivities)
    values <- .connectivity_triplet_values(triplet)
    if (anyNA(values) || any(!is.finite(values))) {
      stop("Connectivity values must all be finite.", call. = FALSE)
    }
    if (any(values < 0)) {
      stop(
        "Connectivity weights must all be non-negative.",
        call. = FALSE
      )
    }
    if (any(values == 0)) {
      stop(
        "Sparse connectivity storage must omit zero-weight coordinates.",
        call. = FALSE
      )
    }
    diagonal_index <- seq_len(n_cells)
    diagonal_values <- connectivities[cbind(diagonal_index, diagonal_index)]
    if (any(diagonal_values != 0)) {
      stop(
        "Connectivity diagonal values must all be exactly zero.",
        call. = FALSE
      )
    }
    if (length(Matrix::which(connectivities != Matrix::t(connectivities))) != 0L) {
      stop("Connectivity matrix must be exactly symmetric.", call. = FALSE)
    }

    upper_triplet <- .connectivity_triplet(
      Matrix::triu(connectivities, k = 1L)
    )
    weights <- as.double(.connectivity_triplet_values(upper_triplet))
    sources <- as.numeric(upper_triplet$i - 1L)
    destinations <- as.numeric(upper_triplet$j - 1L)
  } else {
    if (
      !is.matrix(connectivities) ||
        !(is.numeric(connectivities) || is.logical(connectivities)) ||
        is.object(connectivities) ||
        is.complex(connectivities)
    ) {
      stop(
        "connectivities must be a real numeric or logical matrix or a ",
        "Matrix::Matrix sparse matrix.",
        call. = FALSE
      )
    }
    if (!identical(as.numeric(dim(connectivities)), as.numeric(expected_shape))) {
      stop(
        "Connectivity matrix shape must exactly match the cell axis: expected (",
        n_cells, ", ", n_cells, "), got (",
        nrow(connectivities), ", ", ncol(connectivities), ").",
        call. = FALSE
      )
    }
    if (anyNA(connectivities) || any(!is.finite(connectivities))) {
      stop("Connectivity values must all be finite.", call. = FALSE)
    }
    if (any(connectivities < 0)) {
      stop(
        "Connectivity weights must all be non-negative.",
        call. = FALSE
      )
    }
    if (any(diag(connectivities) != 0)) {
      stop(
        "Connectivity diagonal values must all be exactly zero.",
        call. = FALSE
      )
    }
    if (!all(connectivities == t(connectivities))) {
      stop("Connectivity matrix must be exactly symmetric.", call. = FALSE)
    }

    coordinates <- which(
      upper.tri(connectivities) & connectivities > 0,
      arr.ind = TRUE
    )
    sources <- as.numeric(coordinates[, 1L] - 1L)
    destinations <- as.numeric(coordinates[, 2L] - 1L)
    weights <- as.double(connectivities[coordinates])
  }

  if (length(sources) > 0L) {
    order_index <- order(sources, destinations)
    sources <- sources[order_index]
    destinations <- destinations[order_index]
    weights <- weights[order_index]
  }
  degrees <- tabulate(
    c(sources + 1, destinations + 1),
    nbins = n_cells
  )
  max_neighbors <- if (length(degrees) > 0L) max(degrees) else 0L

  list(
    sources = sources,
    destinations = destinations,
    weights = weights,
    n_edges = length(sources),
    max_neighbors = max_neighbors,
    index_dtype = info$index_dtype,
    index_bytes = info$index_bytes
  )
}

.write_connectivity_edges <- function(connectivity_edges, out_dir, compression) {
  sources_fname <- "edges.src.bin"
  dests_fname <- "edges.dst.bin"
  weights_fname <- "edges.weights.f64.bin"
  sources_path <- file.path(out_dir, sources_fname)
  dests_path <- file.path(out_dir, dests_fname)
  weights_path <- file.path(out_dir, weights_fname)

  if (identical(connectivity_edges$index_dtype, "uint16")) {
    .write_uint16(
      sources_path,
      connectivity_edges$sources,
      compression = compression
    )
    .write_uint16(
      dests_path,
      connectivity_edges$destinations,
      compression = compression
    )
  } else {
    .write_uint32(
      sources_path,
      connectivity_edges$sources,
      compression = compression
    )
    .write_uint32(
      dests_path,
      connectivity_edges$destinations,
      compression = compression
    )
  }
  .write_float64_vector(
    weights_path,
    connectivity_edges$weights,
    compression = compression
  )
  invisible(NULL)
}

.connectivity_index_dtype <- function(n_cells) {
  if (
    !is.numeric(n_cells) ||
      length(n_cells) != 1L ||
      is.na(n_cells) ||
      !is.finite(n_cells) ||
      n_cells <= 0 ||
      n_cells != floor(n_cells)
  ) {
    stop("n_cells must be one positive integer.", call. = FALSE)
  }
  if (n_cells <= 65536) {
    list(index_dtype = "uint16", index_bytes = 2L)
  } else if (n_cells <= 4294967296) {
    list(index_dtype = "uint32", index_bytes = 4L)
  } else {
    stop(
      "Connectivity cannot exceed 4,294,967,296 cells because the current ",
      "browser contract supports at most the complete uint32 index domain.",
      call. = FALSE
    )
  }
}

.export_vector_fields <- function(
    vector_fields,
    embeddings,
    normalization_info,
    out_dir,
    compression,
    vector_field_default
) {
  if (
    !is.list(vector_fields) ||
      is.object(vector_fields) ||
      is.null(names(vector_fields)) ||
      length(vector_fields) == 0L
  ) {
    stop("vector_fields must be one non-empty named list.", call. = FALSE)
  }
  input_names <- names(vector_fields)
  if (anyNA(input_names) || any(!nzchar(input_names))) {
    stop(
      "vector_fields must have native non-empty string names.",
      call. = FALSE
    )
  }
  if (anyDuplicated(input_names)) {
    stop("vector_fields names must be unique.", call. = FALSE)
  }

  grouped <- list()
  for (input_name in input_names) {
    value <- vector_fields[[input_name]]
    if (is.null(value)) {
      stop(
        "Vector field '", input_name, "' cannot be NULL.",
        call. = FALSE
      )
    }

    # The key names no file any more, so the grammar constrains only the shape
    # the caller declares. cellucid-python's _VECTOR_KEY_PATTERN is exactly
    # this, and a key it accepts must not be refused here.
    explicit_match <- regexec(
      "^(.+_umap)_([123])d$",
      input_name,
      perl = TRUE
    )
    explicit_parts <- regmatches(input_name, explicit_match)[[1]]
    if (length(explicit_parts) == 0L) {
      stop(
        "Vector field key '", input_name,
        "' must exactly match '<field>_umap_<1|2|3>d'.",
        call. = FALSE
      )
    }
    field_id <- explicit_parts[[2]]
    declared_dimension <- as.integer(explicit_parts[[3]])

    inferred <- .infer_vector_shape(value, input_name)
    dimension <- inferred$dim
    if (declared_dimension != dimension) {
      stop(
        "Vector field '", input_name, "' declares ", declared_dimension,
        "D but contains ", dimension, " components.",
        call. = FALSE
      )
    }
    dimension_key <- as.character(dimension)
    if (!dimension_key %in% names(embeddings)) {
      stop(
        "Vector field '", input_name, "' requires a matching ",
        dimension, "D embedding.",
        call. = FALSE
      )
    }
    if (nrow(inferred$dense) != nrow(embeddings[[dimension_key]])) {
      stop(
        "Vector field '", input_name, "' has ", nrow(inferred$dense),
        " rows; expected ", nrow(embeddings[[dimension_key]]), ".",
        call. = FALSE
      )
    }

    if (is.null(grouped[[field_id]])) {
      grouped[[field_id]] <- list()
    }
    if (!is.null(grouped[[field_id]][[dimension_key]])) {
      stop(
        "vector_fields declares the same field and dimension more than once: ",
        field_id, " ", dimension, "D.",
        call. = FALSE
      )
    }

    normalization <- normalization_info[[dimension_key]]
    scale_factor <- normalization$scale_factor
    if (
      !is.numeric(scale_factor) ||
        length(scale_factor) != 1L ||
        is.na(scale_factor) ||
        !is.finite(scale_factor) ||
        scale_factor <= 0
    ) {
      stop(
        "Vector field '", input_name,
        "' has no exact positive finite embedding scale.",
        call. = FALSE
      )
    }
    grouped[[field_id]][[dimension_key]] <- (
      inferred$dense * scale_factor
    )
  }

  # cellucid-python emits its validated fields in sorted key order, and the
  # payload index is that order, so the two writers must sort identically.
  # method = "radix" orders by code point exactly as Python's sorted() does,
  # instead of by the session's collation locale.
  field_ids <- sort(names(grouped), method = "radix")
  .require_field_identities(field_ids, what = "Vector field")

  if (is.null(vector_field_default)) {
    if (length(field_ids) != 1L) {
      stop(
        "vector_field_default is required when more than one vector field ",
        "is supplied.",
        call. = FALSE
      )
    }
    default_field <- field_ids[[1]]
  } else {
    default_field <- .validate_required_string(
      vector_field_default,
      "vector_field_default"
    )
    if (!default_field %in% field_ids) {
      stop(
        "vector_field_default '", default_field,
        "' is not present in vector_fields.",
        call. = FALSE
      )
    }
  }

  vectors_dir <- file.path(out_dir, "vectors")
  .dir_create(vectors_dir)
  fields_meta <- list()
  gz_suffix <- if (!is.null(compression)) ".gz" else ""

  payload_indices <- list()
  for (idx in seq_along(field_ids)) {
    field_id <- field_ids[[idx]]
    payload_index <- .payload_index(idx)
    payload_indices[[length(payload_indices) + 1L]] <- payload_index
    dimensions <- sort(as.integer(names(grouped[[field_id]])))
    files <- list()
    for (dimension in dimensions) {
      filename <- sprintf("%d_%dd.bin", payload_index, dimension)
      .write_float32_matrix_row_major(
        file.path(vectors_dir, filename),
        grouped[[field_id]][[as.character(dimension)]],
        compression = compression
      )
      files[[sprintf("%dd", dimension)]] <- sprintf(
        "vectors/%s%s",
        filename,
        gz_suffix
      )
    }
    fields_meta[[field_id]] <- list(
      label = field_id,
      available_dimensions = I(dimensions),
      default_dimension = max(dimensions),
      files = files,
      basis = "umap"
    )
  }
  list(
    identity = list(
      default_field = default_field,
      fields = fields_meta
    ),
    payload_indices = payload_indices
  )
}

.infer_vector_shape <- function(arr, name) {
  if (inherits(arr, "Matrix")) {
    dense <- as.matrix(arr)
  } else if (
    is.numeric(arr) &&
      !is.object(arr) &&
      is.null(dim(arr))
  ) {
    dense <- matrix(arr, ncol = 1L)
  } else if (
    is.matrix(arr) &&
      is.numeric(arr) &&
      !is.object(arr)
  ) {
    dense <- arr
  } else if (
    is.data.frame(arr) &&
      all(vapply(arr, is.numeric, logical(1)))
  ) {
    dense <- as.matrix(arr)
  } else {
    stop(
      "Vector field '", name,
      "' must be a finite real vector or matrix.",
      call. = FALSE
    )
  }
  storage.mode(dense) <- "double"
  if (anyNA(dense) || any(!is.finite(dense))) {
    stop(
      "Vector field '", name, "' must contain only finite values.",
      call. = FALSE
    )
  }
  if (!ncol(dense) %in% c(1L, 2L, 3L)) {
    stop(
      "Vector field '", name, "' must have exactly 1, 2, or 3 components.",
      call. = FALSE
    )
  }
  list(dim = as.integer(ncol(dense)), dense = dense)
}
