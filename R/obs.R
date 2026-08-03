# The observation axis: deciding what kind of field a column is, turning a
# categorical column into categories and codes, and writing every obs
# payload with the schemas that describe them.

# How wide a categorical field's codes are stored, and the code value that
# stands for a missing one. The width the caller chose fixes both, and it also
# fixes how many categories the field may carry: the top code of the domain is
# reserved for "missing", so a uint8 field addresses 255 categories and a
# uint16 field 65,535. Every place that needs a width asks here, so the limit
# and the reserved code cannot come apart.
.categorical_code_dtype <- function(n_categories, key, obs_categorical_dtype) {
  if (identical(obs_categorical_dtype, "uint8")) {
    if (n_categories > 255L) {
      stop(
        "Field '", key, "' has ", n_categories,
        " categories, but uint8 supports at most 255."
      )
    }
    return(list(dtype = "uint8", missing_value = 255L))
  }
  if (n_categories > 65535L) {
    stop(
      "Field '", key, "' has ", n_categories,
      " categories, but uint16 supports at most 65,535."
    )
  }
  list(dtype = "uint16", missing_value = 65535L)
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
    } else {
      .quantization_dtype(obs_continuous_quantization)
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
  dtype_str <- .categorical_code_dtype(
    n_categories,
    key,
    obs_categorical_dtype
  )$dtype

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
          "' must contain only finite values. ",
          .describe_non_finite(
            values,
            paste0("Observation field '", key, "'")
          ),
          call. = FALSE
        )
      }

      if (!is.null(obs_continuous_quantization)) {
        q <- .quantize_continuous(values, bits = obs_continuous_quantization, field_name = key)
        ext <- .quantization_ext(obs_continuous_quantization)
        dtype_str <- .quantization_dtype(obs_continuous_quantization)
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
    dtype_choice <- .categorical_code_dtype(
      n_categories,
      key,
      obs_categorical_dtype
    )

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
      oq_ext <- .quantization_ext(obs_continuous_quantization)
      oq_dtype_str <- .quantization_dtype(obs_continuous_quantization)

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
          codesExt = .codes_extension_by_dtype[[dtype_str]],
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
          codesExt = .codes_extension_by_dtype[[dtype_str]],
          outlierExt = "f32",
          outlierDtype = "float32",
          outlierQuantized = FALSE
        )
      }
    }
  }

  gz_suffix <- .gz_suffix(compression)
  # The schema map is a JSON object even when it is empty, which it is for an
  # export carrying no observation field at all. An empty unnamed list would
  # serialize as `[]`, and the viewer refuses a non-object here.
  schemas <- .json_object()
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
