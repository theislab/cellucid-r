#' Export single-cell data to the Cellucid viewer format
#'
#' `cellucid_prepare()` mirrors the `prepare()` function from the `cellucid`
#' Python package and writes a directory of binary (`.bin`, `.f32`, `.u8`, `.u16`)
#' and JSON manifest files consumed by the Cellucid WebGL viewer.
#'
#' Gene identifiers and `obs` keys are sanitized for filenames; they must be
#' unique after sanitization, otherwise export fails to prevent silent
#' overwrites.
#'
#' @param latent_space Numeric matrix-like of shape `(n_cells, n_dims)`.
#'   Required for per-cell outlier quantile calculations for categorical fields.
#' @param obs A `data.frame` of cell metadata with `n_cells` rows.
#' @param var Optional `data.frame` of feature metadata with `n_genes` rows.
#'   Required if `gene_expression` is provided.
#' @param gene_expression Optional numeric matrix-like of shape
#'   `(n_cells, n_genes)`.
#' @param var_gene_id_column Column name in `var` containing gene identifiers, or
#'   `"index"` (default) to use `rownames(var)`.
#' @param gene_identifiers Optional character vector specifying which genes to
#'   export. If `NULL`, all genes are exported.
#' @param connectivities Optional adjacency/connectivity matrix of shape
#'   `(n_cells, n_cells)`. Sparse matrices from the `Matrix` package are
#'   supported.
#' @param out_dir Output directory path. Defaults to `file.path(getwd(),
#'   "exports")`.
#' @param obs_keys Optional character vector of column names in `obs` to export.
#'   If `NULL`, all columns are exported.
#' @param centroid_outlier_quantile Quantile of distances to keep as inliers
#'   when computing centroids for categorical fields. Set to `NULL` to disable
#'   centroid computation.
#' @param centroid_min_points Minimum number of points in a category required to
#'   compute centroids and outlier quantiles.
#' @param obs_manifest_filename Filename for the obs manifest JSON.
#' @param obs_binary_dirname Directory name (under `out_dir`) for obs binaries.
#' @param var_manifest_filename Filename for the var manifest JSON.
#' @param var_binary_dirname Directory name (under `out_dir`) for var binaries.
#' @param connectivity_manifest_filename Filename for the connectivity manifest
#'   JSON.
#' @param connectivity_binary_dirname Directory name (under `out_dir`) for
#'   connectivity binaries.
#' @param force If `TRUE`, overwrite existing files; if `FALSE`, skip writing
#'   files whose manifest already exists (and skip embedding/vector files that
#'   already exist).
#' @param var_quantization Integer bits for gene expression quantization (`8` or
#'   `16`), or `NULL` for full float32 export.
#' @param obs_continuous_quantization Integer bits for continuous obs field
#'   quantization (`8` or `16`), or `NULL` for full float32 export.
#' @param obs_categorical_dtype One of `"auto"`, `"uint8"`, or `"uint16"` for
#'   categorical obs codes.
#' @param compression Optional gzip compression level (`1`-`9`). Use `NULL` or
#'   `0` to disable compression.
#' @param dataset_name Optional dataset name for `dataset_identity.json`.
#' @param dataset_description Optional dataset description.
#' @param dataset_id Optional dataset identifier. If `NULL`, a filesystem-safe
#'   identifier is derived from `dataset_name` or the `out_dir` name.
#' @param source_name Optional data source name for metadata.
#' @param source_url Optional data source URL for metadata.
#' @param source_citation Optional citation text for metadata.
#' @param X_umap_1d Optional 1D embedding matrix of shape `(n_cells, 1)`.
#' @param X_umap_2d Optional 2D embedding matrix of shape `(n_cells, 2)`.
#' @param X_umap_3d Optional 3D embedding matrix of shape `(n_cells, 3)`.
#' @param X_umap_4d Reserved for future development. Must be `NULL`.
#' @param vector_fields Optional named list of per-cell displacement vectors.
#'   Each element must be a numeric vector (1D) or matrix with 1/2/3 columns and
#'   `n_cells` rows.
#'
#' @return Invisibly returns `NULL`. Called for its side effects (file export).
#' @export
cellucid_prepare <- function(
    latent_space = NULL,
    obs = NULL,
    var = NULL,
    gene_expression = NULL,
    var_gene_id_column = "index",
    gene_identifiers = NULL,
    connectivities = NULL,
    out_dir = file.path(getwd(), "exports"),
    obs_keys = NULL,
    centroid_outlier_quantile = 0.95,
    centroid_min_points = 10,
    obs_manifest_filename = "obs_manifest.json",
    obs_binary_dirname = "obs",
    var_manifest_filename = "var_manifest.json",
    var_binary_dirname = "var",
    connectivity_manifest_filename = "connectivity_manifest.json",
    connectivity_binary_dirname = "connectivity",
    force = FALSE,
    var_quantization = NULL,
    obs_continuous_quantization = NULL,
    obs_categorical_dtype = "auto",
    compression = NULL,
    dataset_name = NULL,
    dataset_description = NULL,
    dataset_id = NULL,
    source_name = NULL,
    source_url = NULL,
    source_citation = NULL,
    X_umap_1d = NULL,
    X_umap_2d = NULL,
    X_umap_3d = NULL,
    X_umap_4d = NULL,
    vector_fields = NULL
) {
  manifest_format_version <- "compact_v1"

  if (!is.null(X_umap_4d)) {
    stop(
      "4D visualization is not yet implemented. ",
      "The X_umap_4d parameter is reserved for future development. ",
      "Please use X_umap_1d, X_umap_2d, or X_umap_3d."
    )
  }

  compression <- .normalize_compression(compression)
  obs_categorical_dtype <- match.arg(obs_categorical_dtype, c("auto", "uint8", "uint16"))

  if (!is.null(var_quantization)) {
    .validate_quantization_bits(var_quantization, arg = "var_quantization")
  }
  if (!is.null(obs_continuous_quantization)) {
    .validate_quantization_bits(obs_continuous_quantization, arg = "obs_continuous_quantization")
  }

  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
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
    check_path <- if (!is.null(compression)) paste0(file_base, ".gz") else file_base
    if (.file_exists_skip(check_path, description = basename(check_path), force = force)) {
      next
    }
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

  if (is.null(obs)) {
    stop("obs data.frame is required.")
  }
  obs <- as.data.frame(obs, stringsAsFactors = FALSE)
  if (nrow(obs) != n_cells) {
    stop("obs has ", nrow(obs), " rows, but embeddings have ", n_cells, " cells.")
  }

  if (is.null(obs_keys)) {
    obs_keys <- names(obs)
  } else {
    obs_keys <- as.character(obs_keys)
    missing <- setdiff(obs_keys, names(obs))
    if (length(missing) > 0) {
      stop(
        "obs_keys contain columns not in obs: ",
        paste(missing, collapse = ", "),
        ". Available columns: ",
        paste(names(obs), collapse = ", ")
      )
    }
  }

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
  if (!is.null(vector_fields)) {
    vector_result <- .export_vector_fields(
      vector_fields = vector_fields,
      embeddings = embeddings,
      normalization_info = normalization_info,
      out_dir = out_dir,
      compression = compression,
      force = force
    )
    vector_fields_identity <- vector_result$identity
  }

  # ===========================================================================
  # OBS EXPORT
  # ===========================================================================
  obs_manifest_path <- file.path(out_dir, obs_manifest_filename)
  if (!.file_exists_skip(obs_manifest_path, description = "obs manifest", force = force)) {
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
        as.numeric(centroid_outlier_quantile)
      } else {
        NULL
      },
      latent_key = "latent_space",
      compression = if (!is.null(compression)) as.integer(compression) else NULL,
      `_obsSchemas` = obs_export$schemas,
      `_continuousFields` = obs_export$continuous_fields,
      `_categoricalFields` = obs_export$categorical_fields
    )
    .write_json(obs_manifest_path, obs_manifest_payload, pretty = FALSE)
  }

  # ===========================================================================
  # VAR (GENE EXPRESSION) EXPORT
  # ===========================================================================
  genes_to_export <- character(0)
  if (!is.null(gene_expression)) {
    if (is.null(var)) {
      stop("var data.frame must be provided when gene_expression is given.")
    }
    var <- as.data.frame(var, stringsAsFactors = FALSE)

    gene_expression <- .as_matrix_like(gene_expression, name = "gene_expression")
    if (nrow(gene_expression) != n_cells) {
      stop("gene_expression has ", nrow(gene_expression), " cells, but embeddings have ", n_cells, " cells.")
    }
    if (nrow(var) != ncol(gene_expression)) {
      stop("var has ", nrow(var), " rows, but gene_expression has ", ncol(gene_expression), " genes.")
    }

    all_gene_ids <- .extract_gene_ids(var, var_gene_id_column)
    .validate_gene_ids(all_gene_ids)
    gene_id_to_idx <- stats::setNames(seq_along(all_gene_ids), all_gene_ids)

    if (is.null(gene_identifiers)) {
      genes_to_export <- all_gene_ids
    } else {
      genes_to_export <- as.character(gene_identifiers)
      missing_genes <- setdiff(genes_to_export, names(gene_id_to_idx))
      if (length(missing_genes) > 0) {
        warning(
          length(missing_genes),
          " gene identifiers not found in var: ",
          paste(utils::head(missing_genes, 5), collapse = ", "),
          if (length(missing_genes) > 5) "..."
        )
      }
      genes_to_export <- intersect(genes_to_export, names(gene_id_to_idx))
    }

    var_manifest_path <- file.path(out_dir, var_manifest_filename)
    if (!.file_exists_skip(var_manifest_path, description = "var manifest", force = force)) {
      var_binary_dir <- file.path(out_dir, var_binary_dirname)
      .dir_create(var_binary_dir)

      var_manifest_fields <- vector("list", length(genes_to_export))
      safe_gene_ids <- .assert_unique_filename_components(genes_to_export, what = "Gene identifiers")

      for (idx in seq_along(genes_to_export)) {
        gene_id <- genes_to_export[[idx]]
        gene_idx <- unname(gene_id_to_idx[[gene_id]])
        safe_gene_id <- safe_gene_ids[[idx]]

        values <- .get_gene_column(gene_expression, gene_idx, n_cells = n_cells)

        if (!is.null(var_quantization)) {
          q <- .quantize_continuous(values, bits = var_quantization)
          ext <- if (var_quantization == 8L) "u8" else "u16"
          value_path <- file.path(var_binary_dir, sprintf("%s.values.%s", safe_gene_id, ext))
          if (var_quantization == 8L) {
            .write_uint8(value_path, q$quantized, compression = compression)
          } else {
            .write_uint16(value_path, q$quantized, compression = compression)
          }
          var_manifest_fields[[idx]] <- list(gene_id, q$min_val, q$max_val)
        } else {
          value_path <- file.path(var_binary_dir, sprintf("%s.values.f32", safe_gene_id))
          .write_float32_vector(value_path, values, compression = compression)
          var_manifest_fields[[idx]] <- list(gene_id)
        }
      }

      gz_suffix <- if (!is.null(compression)) ".gz" else ""
      var_schema <- if (!is.null(var_quantization)) {
        ext <- if (var_quantization == 8L) "u8" else "u16"
        dtype_str <- if (var_quantization == 8L) "uint8" else "uint16"
        list(
          kind = "continuous",
          pathPattern = sprintf("%s/{key}.values.%s%s", var_binary_dirname, ext, gz_suffix),
          ext = ext,
          dtype = dtype_str,
          quantized = TRUE,
          quantizationBits = as.integer(var_quantization)
        )
      } else {
        list(
          kind = "continuous",
          pathPattern = sprintf("%s/{key}.values.f32%s", var_binary_dirname, gz_suffix),
          ext = "f32",
          dtype = "float32",
          quantized = FALSE
        )
      }

      var_manifest_payload <- list(
        `_format` = manifest_format_version,
        n_points = as.integer(n_cells),
        var_gene_id_column = var_gene_id_column,
        compression = if (!is.null(compression)) as.integer(compression) else NULL,
        quantization = if (!is.null(var_quantization)) as.integer(var_quantization) else NULL,
        `_varSchema` = var_schema,
        fields = var_manifest_fields
      )
      .write_json(var_manifest_path, var_manifest_payload, pretty = FALSE)
    }
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
    connectivity_manifest_path <- file.path(out_dir, connectivity_manifest_filename)
    if (!.file_exists_skip(connectivity_manifest_path, description = "connectivity manifest", force = force)) {
      connectivity_binary_dir <- file.path(out_dir, connectivity_binary_dirname)
      .dir_create(connectivity_binary_dir)

      conn_export <- .export_connectivity(
        connectivities = connectivities,
        n_cells = n_cells,
        out_dir = connectivity_binary_dir,
        connectivity_binary_dirname = connectivity_binary_dirname,
        compression = compression
      )

      connectivity_meta$n_edges <- conn_export$n_edges
      connectivity_meta$max_neighbors <- conn_export$max_neighbors
      connectivity_meta$index_dtype <- conn_export$index_dtype

      manifest_sources <- paste0(connectivity_binary_dirname, "/", conn_export$sources_fname)
      manifest_dests <- paste0(connectivity_binary_dirname, "/", conn_export$dests_fname)
      if (!is.null(compression)) {
        manifest_sources <- paste0(manifest_sources, ".gz")
        manifest_dests <- paste0(manifest_dests, ".gz")
      }

      connectivity_manifest_payload <- list(
        format = "edge_pairs",
        n_cells = as.integer(n_cells),
        n_edges = as.integer(conn_export$n_edges),
        max_neighbors = as.integer(conn_export$max_neighbors),
        index_bytes = as.integer(conn_export$index_bytes),
        index_dtype = conn_export$index_dtype,
        sourcesPath = manifest_sources,
        destinationsPath = manifest_dests,
        compression = if (!is.null(compression)) as.integer(compression) else NULL
      )
      .write_json(connectivity_manifest_path, connectivity_manifest_payload, pretty = FALSE)
    } else {
      conn_info <- .connectivity_index_dtype(n_cells)
      connectivity_meta$index_dtype <- conn_info$index_dtype
      connectivity_meta$n_edges <- .approximate_edge_count(connectivities)
    }
  }

  # ===========================================================================
  # DATASET IDENTITY
  # ===========================================================================
  identity_path <- file.path(out_dir, "dataset_identity.json")
  cellucid_version <- tryCatch(
    as.character(utils::packageVersion("cellucid")),
    error = function(e) "unknown"
  )

  if (is.null(dataset_id)) {
    if (!is.null(dataset_name)) {
      dataset_id <- .safe_filename_component(dataset_name)
    } else {
      dataset_id <- .safe_filename_component(basename(out_dir))
    }
  }
  if (is.null(dataset_name)) {
    dataset_name <- basename(out_dir)
  }

  n_genes <- if (!is.null(gene_expression)) {
    as.integer(length(genes_to_export))
  } else {
    0L
  }

  n_obs_fields <- length(obs_field_summaries)
  n_categorical_fields <- sum(vapply(obs_field_summaries, function(x) x$kind == "category", logical(1)))
  n_continuous_fields <- sum(vapply(obs_field_summaries, function(x) x$kind == "continuous", logical(1)))

  identity_obs_fields <- lapply(
    obs_field_summaries,
    function(field_info) {
      entry <- list(key = field_info$key, kind = field_info$kind)
      if (identical(field_info$kind, "category") && !is.null(field_info$category_count)) {
        entry$n_categories <- as.integer(field_info$category_count)
      }
      entry
    }
  )

  source_info <- NULL
  if (!is.null(source_name) || !is.null(source_url) || !is.null(source_citation)) {
    source_info <- list()
    if (!is.null(source_name)) source_info$name <- as.character(source_name)
    if (!is.null(source_url)) source_info$url <- as.character(source_url)
    if (!is.null(source_citation)) source_info$citation <- as.character(source_citation)
  }

  export_settings <- list(
    compression = if (!is.null(compression)) as.integer(compression) else NULL,
    var_quantization = if (!is.null(var_quantization)) as.integer(var_quantization) else NULL,
    obs_continuous_quantization = if (!is.null(obs_continuous_quantization)) as.integer(obs_continuous_quantization) else NULL,
    obs_categorical_dtype = obs_categorical_dtype
  )

  gz_suffix <- if (!is.null(compression)) ".gz" else ""
  embeddings_meta <- list(
    available_dimensions = available_dimensions,
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
    description = if (!is.null(dataset_description)) as.character(dataset_description) else "",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
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

  invisible(NULL)
}

#' @rdname cellucid_prepare
#' @export
prepare <- cellucid_prepare

.normalize_compression <- function(compression) {
  if (is.null(compression)) {
    return(NULL)
  }
  compression <- as.integer(compression)
  if (is.na(compression) || compression <= 0L) {
    return(NULL)
  }
  if (compression < 1L || compression > 9L) {
    stop("compression must be an integer between 1 and 9, or NULL/0 to disable.")
  }
  compression
}

.validate_quantization_bits <- function(bits, arg) {
  bits <- as.integer(bits)
  if (is.na(bits) || !(bits %in% c(8L, 16L))) {
    stop(arg, " must be 8, 16, or NULL.")
  }
  bits
}

.dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

.file_exists_skip <- function(path, description, force = FALSE) {
  if (file.exists(path) && !isTRUE(force)) {
    return(TRUE)
  }
  FALSE
}

.safe_filename_component <- function(name) {
  safe <- gsub("[^A-Za-z0-9._-]+", "_", as.character(name))
  safe <- gsub("^[._]+|[._]+$", "", safe)
  if (nchar(safe) == 0) {
    safe <- "field"
  }
  safe
}

.assert_unique_filename_components <- function(keys, what) {
  keys <- as.character(keys)
  safe <- vapply(keys, .safe_filename_component, character(1))

  dup_safe <- unique(safe[duplicated(safe)])
  if (length(dup_safe) == 0L) {
    return(invisible(safe))
  }

  collision_lines <- character(0)
  for (safe_key in dup_safe) {
    originals <- unique(keys[safe == safe_key])
    preview <- paste(sprintf("'%s'", utils::head(originals, 5)), collapse = ", ")
    if (length(originals) > 5L) {
      preview <- paste0(preview, ", ...")
    }
    collision_lines <- c(collision_lines, sprintf("  - '%s' <- %s", safe_key, preview))
  }

  stop(
    what,
    " contains names that collide after filename sanitization.\n",
    "Please rename to avoid collisions. Collisions:\n",
    paste(collision_lines, collapse = "\n"),
    call. = FALSE
  )
}

.as_dense_matrix <- function(x, name) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(dim(x)) || length(dim(x)) != 2L) {
    stop(name, " must be a 2D matrix-like object.")
  }
  x
}

.as_matrix_like <- function(x, name) {
  if (inherits(x, "Matrix")) {
    return(x)
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  x <- as.matrix(x)
  if (is.null(dim(x)) || length(dim(x)) != 2L) {
    stop(name, " must be a 2D matrix-like object.")
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
    if (is.vector(x) && dim_int == 1L) {
      x <- matrix(as.numeric(x), ncol = 1L)
    } else if (is.data.frame(x)) {
      x <- as.matrix(x)
    } else {
      x <- as.matrix(x)
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
  as.integer(n_cells)
}

.normalize_embedding <- function(arr) {
  axis_mins <- apply(arr, 2, min)
  axis_maxs <- apply(arr, 2, max)
  axis_ranges <- axis_maxs - axis_mins
  max_range <- max(axis_ranges)
  if (max_range < 1e-8) {
    max_range <- 1
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

.write_float32_vector <- function(path, values, compression = NULL) {
  con <- .binary_connection(path, compression = compression)
  on.exit(close(con), add = TRUE)
  writeBin(as.numeric(values), con, size = 4L, endian = "little")
  invisible(.final_binary_path(path, compression))
}

.write_uint8 <- function(path, values, compression = NULL) {
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

.write_uint64 <- function(path, values, compression = NULL) {
  con <- .binary_connection(path, compression = compression)
  on.exit(close(con), add = TRUE)
  raw_vals <- .pack_uint64(values)
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

.quantize_continuous <- function(values, bits, field_name = "unknown") {
  bits <- .validate_quantization_bits(bits, arg = "bits")
  values <- as.numeric(values)

  invalid_mask <- is.na(values) | is.infinite(values)
  valid_mask <- !invalid_mask

  n_total <- length(values)
  n_nan <- sum(is.na(values))
  n_inf <- sum(is.infinite(values))
  n_valid <- sum(valid_mask)

  stats <- list(
    n_total = n_total,
    n_valid = n_valid,
    n_nan = n_nan,
    n_inf = n_inf
  )

  if (n_valid == 0L) {
    min_val <- 0
    max_val <- 1
    stats$warning <- "all_invalid"
  } else {
    valid_values <- values[valid_mask]
    min_val <- min(valid_values)
    max_val <- max(valid_values)
    stats$data_min <- min_val
    stats$data_max <- max_val
  }

  if (identical(max_val, min_val)) {
    max_val <- min_val + 1
    stats$constant_value <- TRUE
  }

  if (bits == 8L) {
    max_quant <- 254L
    nan_value <- 255L
  } else {
    max_quant <- 65534L
    nan_value <- 65535L
  }

  scale <- max_quant / (max_val - min_val)

  quantized <- integer(n_total)
  if (n_valid > 0L) {
    normalized <- (values[valid_mask] - min_val) * scale
    clipped <- pmin(pmax(normalized, 0), max_quant)
    quantized[valid_mask] <- as.integer(floor(clipped + 1e-8))
  }
  quantized[invalid_mask] <- nan_value

  list(
    quantized = quantized,
    min_val = as.numeric(min_val),
    max_val = as.numeric(max_val),
    scale = as.numeric(scale),
    stats = stats,
    field_name = field_name
  )
}

.select_category_dtype <- function(n_categories) {
  if (n_categories <= 254L) {
    list(dtype = "uint8", missing_value = 255L)
  } else {
    list(dtype = "uint16", missing_value = 65535L)
  }
}

.summarize_obs_field <- function(s, key, obs_continuous_quantization, obs_categorical_dtype) {
  kind <- if (is.factor(s)) {
    "category"
  } else if (is.logical(s)) {
    "category"
  } else if (is.numeric(s)) {
    "continuous"
  } else {
    "category"
  }

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
        key = as.character(key),
        kind = "continuous",
        quantized = !is.null(obs_continuous_quantization),
        quantization_bits = if (!is.null(obs_continuous_quantization)) as.integer(obs_continuous_quantization) else NULL,
        dtype = dtype_str
      )
    )
  }

  cat <- as.factor(s)
  categories <- levels(cat)
  n_categories <- length(categories)
  dtype_str <- if (identical(obs_categorical_dtype, "auto")) {
    if (n_categories <= 254L) "uint8" else "uint16"
  } else if (identical(obs_categorical_dtype, "uint8")) {
    "uint8"
  } else {
    "uint16"
  }

  list(
    key = as.character(key),
    kind = "category",
    category_count = as.integer(n_categories),
    codes_dtype = dtype_str,
    outlier_quantized = !is.null(obs_continuous_quantization),
    outlier_quantization_bits = if (!is.null(obs_continuous_quantization)) as.integer(obs_continuous_quantization) else NULL
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
  safe_keys <- .assert_unique_filename_components(obs_keys, what = "obs_keys")
  continuous_fields <- list()
  categorical_fields <- list()
  continuous_dtype_info <- list()
  categorical_dtype_info <- list()

  for (idx in seq_along(obs_keys)) {
    key <- obs_keys[[idx]]
    s <- obs[[key]]
    safe_key <- safe_keys[[idx]]

    kind <- if (is.factor(s)) {
      "category"
    } else if (is.logical(s)) {
      "category"
    } else if (is.numeric(s)) {
      "continuous"
    } else {
      "category"
    }

    if (identical(kind, "continuous")) {
      values <- suppressWarnings(as.numeric(s))
      if (length(values) != nrow(obs)) {
        stop("Length mismatch for obs[['", key, "']]: ", length(values), " vs ", nrow(obs))
      }

      if (!is.null(obs_continuous_quantization)) {
        q <- .quantize_continuous(values, bits = obs_continuous_quantization, field_name = key)
        ext <- if (obs_continuous_quantization == 8L) "u8" else "u16"
        dtype_str <- if (obs_continuous_quantization == 8L) "uint8" else "uint16"
        value_path <- file.path(obs_binary_dir, sprintf("%s.values.%s", safe_key, ext))
        if (obs_continuous_quantization == 8L) {
          .write_uint8(value_path, q$quantized, compression = compression)
        } else {
          .write_uint16(value_path, q$quantized, compression = compression)
        }
        continuous_fields[[length(continuous_fields) + 1L]] <- list(key, q$min_val, q$max_val)
        if (length(continuous_dtype_info) == 0L) {
          continuous_dtype_info <- list(
            ext = ext,
            dtype = dtype_str,
            quantized = TRUE,
            quantizationBits = as.integer(obs_continuous_quantization)
          )
        }
      } else {
        value_path <- file.path(obs_binary_dir, sprintf("%s.values.f32", safe_key))
        .write_float32_vector(value_path, values, compression = compression)
        continuous_fields[[length(continuous_fields) + 1L]] <- list(key)
        if (length(continuous_dtype_info) == 0L) {
          continuous_dtype_info <- list(ext = "f32", dtype = "float32", quantized = FALSE)
        }
      }

      next
    }

    cat <- as.factor(s)
    categories <- as.character(levels(cat))
    codes <- as.integer(cat) - 1L
    codes[is.na(codes)] <- -1L

    n_categories <- length(categories)
    dtype_choice <- if (identical(obs_categorical_dtype, "auto")) {
      .select_category_dtype(n_categories)
    } else if (identical(obs_categorical_dtype, "uint8")) {
      if (n_categories > 254L) {
        stop(
          "Field '", key, "' has ", n_categories, " categories, but uint8 can only hold 254. ",
          "Use 'auto' or 'uint16'."
        )
      }
      list(dtype = "uint8", missing_value = 255L)
    } else {
      list(dtype = "uint16", missing_value = 65535L)
    }

    missing_value <- dtype_choice$missing_value
    dtype_str <- dtype_choice$dtype

    codes_typed <- rep(missing_value, length(codes))
    valid_mask <- codes >= 0L
    codes_typed[valid_mask] <- codes[valid_mask]

    if (identical(dtype_str, "uint8")) {
      codes_fname <- sprintf("%s.codes.u8", safe_key)
      codes_path <- file.path(obs_binary_dir, codes_fname)
      .write_uint8(codes_path, codes_typed, compression = compression)
    } else {
      codes_fname <- sprintf("%s.codes.u16", safe_key)
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
      oq <- .quantize_continuous(outlier_quantiles, bits = obs_continuous_quantization, field_name = paste0(key, "_outliers"))
      oq_ext <- if (obs_continuous_quantization == 8L) "u8" else "u16"
      oq_dtype_str <- if (obs_continuous_quantization == 8L) "uint8" else "uint16"

      outlier_fname <- sprintf("%s.outliers.%s", safe_key, oq_ext)
      outlier_path <- file.path(obs_binary_dir, outlier_fname)
      if (obs_continuous_quantization == 8L) {
        .write_uint8(outlier_path, oq$quantized, compression = compression)
      } else {
        .write_uint16(outlier_path, oq$quantized, compression = compression)
      }

      categorical_fields[[length(categorical_fields) + 1L]] <- list(
        key,
        categories,
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
      outlier_fname <- sprintf("%s.outliers.f32", safe_key)
      outlier_path <- file.path(obs_binary_dir, outlier_fname)
      .write_float32_vector(outlier_path, outlier_quantiles, compression = compression)

      categorical_fields[[length(categorical_fields) + 1L]] <- list(
        key,
        categories,
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
        "%s/{key}.values.%s%s",
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
      codesPathPattern = sprintf("%s/{key}.codes.{ext}%s", obs_binary_dirname, gz_suffix),
      outlierPathPattern = sprintf(
        "%s/{key}.outliers.%s%s",
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

  if (!(outlier_quantile > 0.5 && outlier_quantile < 1)) {
    outlier_quantile <- 0.95
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
      category = as.character(label),
      position = as.numeric(center),
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
  if (is.null(var_gene_id_column) || identical(var_gene_id_column, "index")) {
    ids <- rownames(var)
    if (is.null(ids)) {
      ids <- as.character(seq_len(nrow(var)) - 1L)
    }
    return(as.character(ids))
  }
  if (!var_gene_id_column %in% names(var)) {
    stop(
      "var_gene_id_column '", var_gene_id_column, "' not found in var. Available columns: ",
      paste(names(var), collapse = ", ")
    )
  }
  as.character(var[[var_gene_id_column]])
}

.validate_gene_ids <- function(gene_ids) {
  gene_ids <- as.character(gene_ids)

  if (anyNA(gene_ids)) {
    stop("Gene identifiers contain NA. Provide a non-missing ID for every gene.", call. = FALSE)
  }
  if (any(gene_ids == "")) {
    stop("Gene identifiers contain empty strings. Provide a non-empty ID for every gene.", call. = FALSE)
  }

  dup <- unique(gene_ids[duplicated(gene_ids)])
  if (length(dup) > 0L) {
    preview <- paste(sprintf("'%s'", utils::head(dup, 5)), collapse = ", ")
    if (length(dup) > 5L) {
      preview <- paste0(preview, ", ...")
    }
    stop("Gene identifiers must be unique. Duplicates: ", preview, call. = FALSE)
  }

  invisible(gene_ids)
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
  values <- as.integer(values)
  if (any(is.na(values))) {
    stop("uint16 values contain NA.")
  }
  values <- values %% 65536L
  lo <- as.raw(bitwAnd(values, 0xFFL))
  hi <- as.raw(bitwAnd(bitwShiftR(values, 8L), 0xFFL))
  as.raw(rbind(lo, hi))
}

.pack_uint32 <- function(values) {
  values <- as.numeric(values)
  if (anyNA(values)) {
    stop("uint32 values contain NA.")
  }
  values <- values %% 4294967296
  b0 <- as.raw(floor(values) %% 256)
  b1 <- as.raw(floor(values / 256) %% 256)
  b2 <- as.raw(floor(values / 65536) %% 256)
  b3 <- as.raw(floor(values / 16777216) %% 256)
  as.raw(rbind(b0, b1, b2, b3))
}

.pack_uint64 <- function(values) {
  values <- as.numeric(values)
  if (anyNA(values)) {
    stop("uint64 values contain NA.")
  }
  two32 <- 4294967296
  lo <- values %% two32
  hi <- floor(values / two32)

  lo_raw <- .pack_uint32(lo)
  hi_raw <- .pack_uint32(hi)

  raw <- as.raw(rbind(matrix(lo_raw, nrow = 4), matrix(hi_raw, nrow = 4)))
  raw
}

.export_connectivity <- function(connectivities, n_cells, out_dir, connectivity_binary_dirname, compression) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required to export connectivity matrices.")
  }

  if (inherits(connectivities, "Matrix")) {
    conn <- connectivities
  } else if (is.matrix(connectivities)) {
    conn <- Matrix::Matrix(connectivities, sparse = TRUE)
  } else {
    stop("connectivities must be a matrix or a Matrix::Matrix sparse matrix.")
  }

  if (nrow(conn) != n_cells || ncol(conn) != n_cells) {
    stop(
      "Connectivity matrix shape (", nrow(conn), ", ", ncol(conn), ") does not match number of cells ", n_cells, "."
    )
  }

  info <- .connectivity_index_dtype(n_cells)
  index_dtype <- info$index_dtype
  index_bytes <- info$index_bytes

  conn_sym <- conn + Matrix::t(conn)
  if (length(conn_sym@x) > 0L) {
    conn_sym@x[] <- 1
  }

  trip <- Matrix::summary(conn_sym)
  max_neighbors <- max(tabulate(trip$i, nbins = n_cells))

  src <- trip$i - 1L
  dst <- trip$j - 1L
  keep <- src < dst
  src <- src[keep]
  dst <- dst[keep]
  n_edges <- length(src)

  ord <- order(src, dst)
  src <- src[ord]
  dst <- dst[ord]

  sources_fname <- "edges.src.bin"
  dests_fname <- "edges.dst.bin"
  sources_path <- file.path(out_dir, sources_fname)
  dests_path <- file.path(out_dir, dests_fname)

  if (identical(index_dtype, "uint16")) {
    .write_uint16(sources_path, src, compression = compression)
    .write_uint16(dests_path, dst, compression = compression)
  } else if (identical(index_dtype, "uint32")) {
    .write_uint32(sources_path, src, compression = compression)
    .write_uint32(dests_path, dst, compression = compression)
  } else {
    .write_uint64(sources_path, src, compression = compression)
    .write_uint64(dests_path, dst, compression = compression)
  }

  list(
    n_edges = as.integer(n_edges),
    max_neighbors = as.integer(max_neighbors),
    index_dtype = index_dtype,
    index_bytes = as.integer(index_bytes),
    sources_fname = sources_fname,
    dests_fname = dests_fname
  )
}

.connectivity_index_dtype <- function(n_cells) {
  n_cells <- as.numeric(n_cells)
  if (n_cells <= 65535) {
    list(index_dtype = "uint16", index_bytes = 2L)
  } else if (n_cells <= 4294967295) {
    list(index_dtype = "uint32", index_bytes = 4L)
  } else {
    list(index_dtype = "uint64", index_bytes = 8L)
  }
}

.approximate_edge_count <- function(connectivities) {
  if (inherits(connectivities, "Matrix")) {
    return(as.integer(length(connectivities@x) %/% 2L))
  }
  if (is.matrix(connectivities)) {
    return(as.integer(sum(connectivities != 0) %/% 2L))
  }
  NULL
}

.export_vector_fields <- function(vector_fields, embeddings, normalization_info, out_dir, compression, force) {
  if (!is.list(vector_fields) || is.null(names(vector_fields))) {
    stop("vector_fields must be a named list of arrays.")
  }
  nms <- names(vector_fields)
  if (anyNA(nms) || any(nms == "")) {
    stop("vector_fields must be a named list with non-empty names.")
  }
  if (anyDuplicated(nms)) {
    dup <- unique(nms[duplicated(nms)])
    preview <- paste(sprintf("'%s'", utils::head(dup, 5)), collapse = ", ")
    if (length(dup) > 5L) {
      preview <- paste0(preview, ", ...")
    }
    stop("vector_fields names must be unique. Duplicates: ", preview, call. = FALSE)
  }

  vectors_dir <- file.path(out_dir, "vectors")
  .dir_create(vectors_dir)

  suffix_re <- "^(.+)_([123])d$"
  grouped <- list()
  explicit_dims <- list()

  # Pass 1: explicit keys (<field>_<dim>d)
  for (nm in names(vector_fields)) {
    arr <- vector_fields[[nm]]
    if (is.null(arr)) next

    m <- regexec(suffix_re, nm)
    parts <- regmatches(nm, m)[[1]]
    if (length(parts) == 0) next

    field_id <- parts[[2]]
    dim <- as.integer(parts[[3]])

    explicit_dims[[field_id]] <- union(explicit_dims[[field_id]], dim)
    if (is.null(grouped[[field_id]])) grouped[[field_id]] <- list()
    grouped[[field_id]][[as.character(dim)]] <- arr
  }

  # Pass 2: implicit keys (<field>_umap), infer dim from array shape
  for (nm in names(vector_fields)) {
    arr <- vector_fields[[nm]]
    if (is.null(arr)) next

    m <- regexec(suffix_re, nm)
    parts <- regmatches(nm, m)[[1]]
    if (length(parts) > 0) next

    inferred <- .infer_vector_shape(arr, nm)
    inferred_dim <- inferred$dim
    field_id <- nm

    if (inferred_dim %in% (explicit_dims[[field_id]] %||% integer(0))) {
      next
    }
    if (is.null(grouped[[field_id]])) grouped[[field_id]] <- list()
    grouped[[field_id]][[as.character(inferred_dim)]] <- arr
  }

  fields_meta <- list()
  gz_suffix <- if (!is.null(compression)) ".gz" else ""

  for (field_id in names(grouped)) {
    safe_id <- .safe_filename_component(field_id)
    if (!identical(safe_id, field_id)) {
      stop(
        "Vector field id '", field_id, "' contains unsupported characters. Use '", safe_id, "' instead."
      )
    }

    by_dim <- grouped[[field_id]]
    files <- list()
    dims <- integer(0)

    for (dim in sort(as.integer(names(by_dim)))) {
      if (!as.character(dim) %in% names(embeddings)) {
        next
      }

      inferred <- .infer_vector_shape(by_dim[[as.character(dim)]], paste0(field_id, "_", dim, "d"))
      if (inferred$dim != dim) {
        stop("Vector field '", field_id, "' declared as ", dim, "D but has shape mismatch.")
      }
      vec <- inferred$dense
      if (nrow(vec) != nrow(embeddings[[as.character(dim)]])) {
        stop(
          "Vector field '", field_id, "' ", dim, "D has ", nrow(vec), " rows, expected ",
          nrow(embeddings[[as.character(dim)]])
        )
      }

      scale_factor <- normalization_info[[as.character(dim)]]$scale_factor %||% 1
      if (!identical(scale_factor, 1)) {
        vec <- vec * as.numeric(scale_factor)
      }

      filename <- sprintf("%s_%dd.bin", field_id, dim)
      path <- file.path(vectors_dir, filename)
      check_path <- if (!is.null(compression)) paste0(path, ".gz") else path
      if (.file_exists_skip(check_path, description = basename(check_path), force = force)) {
        next
      }

      .write_float32_matrix_row_major(path, vec, compression = compression)

      files[[sprintf("%dd", dim)]] <- sprintf("vectors/%s%s", filename, gz_suffix)
      dims <- c(dims, dim)
    }

    if (length(dims) == 0L) next

    entry <- list(
      label = .vector_field_label(field_id),
      available_dimensions = dims,
      default_dimension = max(dims),
      files = files
    )
    if (endsWith(field_id, "_umap")) {
      entry$basis <- "umap"
    }

    fields_meta[[field_id]] <- entry
  }

  if (length(fields_meta) == 0L) {
    return(list(identity = NULL))
  }

  default_field <- if ("velocity_umap" %in% names(fields_meta)) {
    "velocity_umap"
  } else {
    sort(names(fields_meta))[[1]]
  }

  list(
    identity = list(
      default_field = default_field,
      fields = fields_meta
    )
  )
}

.infer_vector_shape <- function(arr, name) {
  if (is.vector(arr)) {
    dense <- matrix(as.numeric(arr), ncol = 1L)
  } else if (is.data.frame(arr)) {
    dense <- as.matrix(arr)
  } else {
    dense <- as.matrix(arr)
  }
  storage.mode(dense) <- "double"

  if (length(dim(dense)) != 2L) {
    stop("Vector field '", name, "' must be 1D or 2D array.")
  }
  if (ncol(dense) %in% c(1L, 2L, 3L)) {
    return(list(dim = as.integer(ncol(dense)), dense = dense))
  }
  stop("Vector field '", name, "' must have 1/2/3 components, got ", ncol(dense), ".")
}

.vector_field_label <- function(field_id) {
  base <- if (endsWith(field_id, "_umap")) substr(field_id, 1L, nchar(field_id) - 5L) else field_id
  base <- trimws(gsub("_", " ", base))
  titled <- if (nchar(base) > 0) paste0(toupper(substr(base, 1L, 1L)), substr(base, 2L, nchar(base))) else field_id
  if (endsWith(field_id, "_umap")) {
    paste0(titled, " (UMAP)")
  } else {
    titled
  }
}

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}
