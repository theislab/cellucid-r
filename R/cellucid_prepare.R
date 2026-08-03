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
    dataset_name,
    dataset_description = NULL,
    dataset_id,
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
    # Before normalization, because normalization hides it. Centring and
    # scaling pull every coordinate into roughly [-1, 1], so a source
    # coordinate of 1e300 or 1e-320 reaches the writer already inside the
    # float32 range and would be published as a number the caller never
    # supplied. cellucid-python gates the same values at the same point.
    .require_finite_float32_source(
      embeddings[[dim]],
      paste0("X_umap_", dim, "d")
    )
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
        shown <- utils::head(missing_genes, 5L)
        remainder <- length(missing_genes) - length(shown)
        stop(
          "gene_identifiers contains identifiers not found in var: ",
          .format_value_list(shown),
          if (remainder > 0L) paste0(" and ", remainder, " more") else "",
          ".",
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

  gz_suffix <- .gz_suffix(compression)
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

  # Each axis directory has now been reconciled against the manifest that
  # declares it, but the point payloads sit at the export root and are declared
  # in dataset_identity.json, so nothing above compared the coordinate files
  # the viewer is told to fetch against the files that were written. The
  # declared name and the written name are two independent spellings of the
  # compression setting -- the sprintf() above and .final_binary_path() at the
  # write -- and a generation whose points disagree publishes successfully and
  # then fails in the browser with no coordinates. The root is the last surface
  # where a declaration and a file can disagree, so it is reconciled whole:
  # the manifests and the axis directories this call created are named here
  # too, and anything else that reached the root is refused rather than
  # published.
  declared_root_files <- c("dataset_identity.json", "obs_manifest.json")
  declared_root_directories <- obs_binary_dirname
  for (dimension_key in names(identity_payload$embeddings$files)) {
    declared_point_path <- identity_payload$embeddings$files[[dimension_key]]
    if (
      !is.character(declared_point_path) ||
        length(declared_point_path) != 1L ||
        is.na(declared_point_path) ||
        !nzchar(declared_point_path)
    ) {
      stop(
        "dataset_identity.json embeddings.files ",
        dimension_key,
        " must be one payload path.",
        call. = FALSE
      )
    }
    declared_root_files <- c(declared_root_files, declared_point_path)
  }
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
