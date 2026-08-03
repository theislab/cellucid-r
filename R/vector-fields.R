# Per-cell vectors drawn over an embedding. A key declares the field and the
# dimension it belongs to, the components are scaled by that embedding's own
# scale so the arrows keep their length relative to the points, and the
# fields are ordered by code point so both writers assign the same indices.

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
  # The whole identity rule is applied here, not the parts of it this caller
  # can still fail. `grouped` is built by name assignment, so its names are
  # already distinct and the uniqueness clause inside
  # .require_field_identities() cannot fire from this call site; the
  # display-text clause can, because the key grammar constrains only the shape
  # and lets any characters through in front of `_umap`. Splitting the rule
  # here to call only the reachable half would put a second spelling of "what
  # a field identity is" in the package, which is the drift this file is
  # supposed to be free of. cellucid-python calls its whole
  # _require_field_identities() on this axis for the same reason, over a dict
  # whose keys are equally unique by construction.
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
  gz_suffix <- .gz_suffix(compression)

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
  } else if (.is_numeric_data_frame(arr)) {
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
      "Vector field '", name, "' must contain only finite values. ",
      .describe_non_finite(
        dense,
        paste0("Vector field '", name, "'")
      ),
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
