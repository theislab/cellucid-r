# The help pages for cellucid_prepare() and cellucid_embedding() are
# hand-written in man/cellucid_prepare.Rd and man/cellucid_embedding.Rd and
# are the only source for them. Nothing in this
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
    var_gene_name_columns = NULL,
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
    dataset_name,
    dataset_description = NULL,
    dataset_id,
    source_name = NULL,
    source_url = NULL,
    source_citation = NULL,
    embeddings,
    embedding_defaults = NULL,
    vector_fields = NULL,
    vector_field_default = NULL,
    created_at = NULL
) {
  manifest_format_version <- .manifest_format_version
  identity_version <- 3L
  embeddings_dirname <- "embeddings"
  obs_binary_dirname <- "obs"
  var_binary_dirname <- "var"
  connectivity_binary_dirname <- "connectivity"

  # An argument that was never supplied and an argument supplied with the wrong
  # value are two different mistakes, and only the message tells the caller
  # which one they made. Every required argument is therefore reported here,
  # before any value is inspected, so a missing one is never described as a
  # value of the wrong type. Each `missing()` sits directly in an `if`
  # condition in this body, which is the one place R guarantees it reads the
  # argument's own supplied-ness.
  absent_required <- character()
  if (missing(obs_categorical_dtype)) {
    absent_required <- c(absent_required, "obs_categorical_dtype")
  }
  if (missing(dataset_name)) {
    absent_required <- c(absent_required, "dataset_name")
  }
  if (missing(dataset_id)) {
    absent_required <- c(absent_required, "dataset_id")
  }
  if (missing(embeddings)) {
    absent_required <- c(absent_required, "embeddings")
  }
  .require_supplied_arguments(absent_required)

  force <- .validate_force(force)
  compression <- .normalize_compression(compression)
  obs_categorical_dtype <- .validate_obs_categorical_dtype(
    obs_categorical_dtype
  )
  dataset_id <- .validate_dataset_id(dataset_id)
  dataset_name <- .validate_dataset_name(dataset_name)
  dataset_description <- .validate_dataset_description(dataset_description)
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

  # Every embedding and each dimension's default are validated, and each
  # item is normalized on its own, before any embedding payload is written.
  resolved_embeddings <- .resolve_embeddings(embeddings, embedding_defaults)
  embedding_items <- resolved_embeddings$items
  n_cells <- resolved_embeddings$n_cells
  default_dimension <- max(vapply(
    embedding_items,
    function(item) item$dimension,
    integer(1)
  ))
  connectivity_edges <- NULL
  if (!is.null(connectivities)) {
    connectivity_edges <- .validate_connectivity_edges(
      connectivities = connectivities,
      n_cells = n_cells
    )
  }

  # Save embeddings: embeddings/{index}_{dimension}d.bin, rounded to float32
  # exactly once, at the write.
  embeddings_dir <- file.path(out_dir, embeddings_dirname)
  .dir_create(embeddings_dir)
  for (item in embedding_items) {
    .write_float32_matrix_row_major(
      file.path(embeddings_dir, sprintf("%d_%dd.bin", item$index, item$dimension)),
      item$coords,
      compression = compression
    )
  }

  if (is.null(latent_space)) {
    stop("latent_space is required for outlier quantile calculation.")
  }
  latent <- .as_dense_matrix(latent_space, name = "latent_space")
  # The latent space is never written, only measured: the outlier quantiles
  # come from distances computed in it. Those distances are published as
  # float32, so a latent coordinate that cannot survive the conversion cannot
  # produce a quantile that means anything, and cellucid-python refuses it for
  # the same reason.
  .require_finite_float32_source(latent, "latent_space")
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

  # Vector fields (optional), each bound to exported embedding ids
  vector_fields_identity <- NULL
  vector_payload_indices <- list()
  resolved_vector_fields <- .resolve_vector_fields(
    vector_fields = vector_fields,
    items = embedding_items,
    vector_field_default = vector_field_default
  )
  if (length(resolved_vector_fields$fields) > 0L) {
    vector_result <- .write_vector_fields(
      resolved_vector_fields,
      items = embedding_items,
      out_dir = out_dir,
      compression = compression
    )
    vector_fields_identity <- vector_result$identity
    vector_payload_indices <- vector_result$payload_indices
  }

  # ===========================================================================
  # OBS EXPORT
  # ===========================================================================
  obs_manifest_path <- file.path(out_dir, "obs_manifest.json")
  obs_export <- .export_obs(
    obs = obs,
    obs_keys = obs_keys,
    latent = latent,
    embedding_items = embedding_items,
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

    if (!is.null(var_gene_id_column)) {
      var_gene_id_column <- .validate_gene_name_source(
        var_gene_id_column,
        "var_gene_id_column"
      )
    }
    default_names <- .resolve_gene_name_set(
      var,
      var_gene_id_column,
      "var_gene_id_column"
    )
    all_gene_ids <- default_names$names
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
        shown <- utils::head(missing_genes, 5L)
        remainder <- length(missing_genes) - length(shown)
        stop(
          "gene_identifiers contains identifiers not found in var: ",
          .format_value_list(shown),
          if (remainder > 0L) paste0(" and ", remainder, " more") else "",
          ". Identifiers are matched against the default gene names after ",
          "repeated names were made unique.",
          call. = FALSE
        )
      }
      genes_to_export <- gene_identifiers
    }
    # Being drawable is a property of a name the viewer shows, so it is
    # checked on the exported genes: a var row left out reaches no manifest,
    # exactly as an obs column left out of obs_keys= does.
    .require_field_identities(genes_to_export, what = "Gene")
    exported_positions <- unname(gene_id_to_idx[genes_to_export])
    gene_name_sets <- .gene_name_sets(
      var,
      default_source = var_gene_id_column,
      default_names = default_names,
      alternative_columns = var_gene_name_columns,
      exported_positions = exported_positions
    )

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
        ext <- .quantization_ext(var_quantization)
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

    gz_suffix <- .gz_suffix(compression)
    var_schema <- if (!is.null(var_quantization)) {
      ext <- .quantization_ext(var_quantization)
      dtype_str <- .quantization_dtype(var_quantization)
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
      compression = if (!is.null(compression)) compression else NULL,
      quantization = if (!is.null(var_quantization)) var_quantization else NULL,
      `_varSchema` = var_schema,
      fields = var_manifest_fields,
      gene_name_sets = gene_name_sets
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
    for (name_set in gene_name_sets) {
      note <- .unique_names_note(
        name_set$filled_from_index,
        name_set$renamed,
        name_set$source
      )
      if (!is.null(note)) {
        source_label <- if (is.null(name_set$source)) {
          "var index"
        } else {
          paste0("var column '", name_set$source, "'")
        }
        message("Gene names from ", source_label, ": ", note)
      }
    }
  } else if (!is.null(var_gene_name_columns)) {
    stop(
      "var_gene_name_columns names alternative gene names, and no ",
      "gene_expression was supplied, so the export has no genes to name.",
      call. = FALSE
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

  # One item per exported embedding, in payload index order, each naming the
  # exact file the format derives for it.
  gz_suffix <- .gz_suffix(compression)
  items_meta <- list()
  for (item in embedding_items) {
    items_meta[[item$id]] <- list(
      label = item$label,
      dimension = as.integer(item$dimension),
      file = sprintf(
        "%s/%d_%dd.bin%s",
        embeddings_dirname,
        item$index,
        item$dimension,
        gz_suffix
      ),
      source = list(
        key = item$source_key,
        columns = I(as.integer(item$columns))
      )
    )
  }
  embeddings_meta <- list(
    default_dimension = as.integer(default_dimension),
    defaults = resolved_embeddings$defaults,
    items = items_meta
  )

  identity_payload <- list(
    version = identity_version,
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

  # Every embedding item is declared by path in dataset_identity.json, so
  # the embeddings directory is reconciled against that declaration -- the
  # declared name and the written name are two expressions of one compression
  # setting, and a generation whose coordinates disagree would publish and
  # then fail in the browser with nothing to draw.
  .require_dense_payload_indices(
    lapply(embedding_items, function(item) item$index),
    axis = "Embedding"
  )
  .require_declared_payloads_on_disk(
    out_dir,
    directory_name = embeddings_dirname,
    declared = .declared_embedding_payload_paths(identity_payload),
    axis = "Embedding"
  )

  # The root is the last surface where a declaration and a file can disagree,
  # so it is reconciled whole: the manifests and the payload directories this
  # call created are named, and anything else that reached the root is
  # refused rather than published.
  declared_root_files <- c("dataset_identity.json", "obs_manifest.json")
  declared_root_directories <- c(obs_binary_dirname, embeddings_dirname)
  if (!is.null(gene_expression)) {
    declared_root_files <- c(declared_root_files, "var_manifest.json")
    declared_root_directories <- c(
      declared_root_directories,
      var_binary_dirname
    )
  }
  if (!is.null(connectivities)) {
    declared_root_files <- c(
      declared_root_files,
      "connectivity_manifest.json"
    )
    declared_root_directories <- c(
      declared_root_directories,
      connectivity_binary_dirname
    )
  }
  if (!is.null(vector_fields_identity)) {
    declared_root_directories <- c(declared_root_directories, "vectors")
  }
  .require_declared_payloads_on_disk(
    out_dir,
    directory_name = NULL,
    declared = unique(declared_root_files),
    axis = "Export",
    declared_directories = unique(declared_root_directories)
  )

  .publish_export_generation(
    out_dir,
    final_out_dir,
    transaction_id = transaction$transaction_id,
    had_target = transaction$had_target
  )
  transaction_cleanup_armed <- FALSE

  invisible(NULL)
}
