# Per-cell vectors drawn over an embedding. A field is bound to embedding ids:
# list(velocity = list(X_umap = v2, X_umap_3d = v3)) binds the field
# `velocity` to two exported embeddings, and each array has one row per cell
# and one column per dimension of the embedding it is drawn over. Components
# are scaled by that embedding's own normalization factor, so the arrows keep
# their length relative to the points. Fields are ordered by code point and
# each field's embeddings follow the item order, so both writers assign the
# same indices: vectors/{fieldIndex}_{embeddingIndex}.bin.

# Validate and scale every binding before anything is written. Returns the
# scaled arrays as list(fields = list(<field> = list(<embedding id> = matrix)),
# default_field = <id or NULL>).
.resolve_vector_fields <- function(vector_fields, items, vector_field_default) {
  if (is.null(vector_fields)) {
    if (!is.null(vector_field_default)) {
      stop(
        "vector_field_default was provided but no vector fields exist.",
        call. = FALSE
      )
    }
    return(list(fields = list(), default_field = NULL))
  }
  if (
    !is.list(vector_fields) ||
      is.object(vector_fields) ||
      is.data.frame(vector_fields)
  ) {
    stop(
      "vector_fields must map each field id to a named list of embedding ids ",
      "and arrays, for example list(velocity = list(X_umap = velocity_2d)).",
      call. = FALSE
    )
  }
  if (!is.null(vector_field_default)) {
    vector_field_default <- .validate_required_string(
      vector_field_default,
      "vector_field_default"
    )
  }
  if (length(vector_fields) == 0L) {
    if (!is.null(vector_field_default)) {
      stop(
        "vector_field_default was provided but no vector fields exist.",
        call. = FALSE
      )
    }
    return(list(fields = list(), default_field = NULL))
  }
  field_names <- names(vector_fields)
  if (is.null(field_names) || anyNA(field_names) || any(!nzchar(field_names))) {
    stop(
      "vector_fields must have native non-empty string names.",
      call. = FALSE
    )
  }
  if (anyDuplicated(field_names)) {
    stop("vector_fields names must be unique.", call. = FALSE)
  }

  item_ids <- vapply(items, function(item) item$id, character(1))
  dimension_by_id <- stats::setNames(
    vapply(items, function(item) item$dimension, integer(1)),
    item_ids
  )
  scale_by_id <- stats::setNames(
    vapply(items, function(item) item$scale_factor, numeric(1)),
    item_ids
  )

  fields <- list()
  for (field_id in field_names) {
    bindings <- vector_fields[[field_id]]
    if (
      !is.list(bindings) ||
        is.object(bindings) ||
        is.data.frame(bindings)
    ) {
      stop(
        "Vector field '", field_id, "' must map embedding ids to arrays, for ",
        "example list(", field_id, " = list(X_umap = vectors)).",
        call. = FALSE
      )
    }
    if (length(bindings) == 0L) {
      stop("Vector field '", field_id, "' binds no embedding.", call. = FALSE)
    }
    embedding_names <- names(bindings)
    if (
      is.null(embedding_names) ||
        anyNA(embedding_names) ||
        any(!nzchar(embedding_names))
    ) {
      stop(
        "Vector field '", field_id, "' must name the embedding of every array.",
        call. = FALSE
      )
    }
    if (anyDuplicated(embedding_names)) {
      stop(
        "Vector field '", field_id, "' binds one embedding more than once.",
        call. = FALSE
      )
    }
    bound <- list()
    for (embedding_id in embedding_names) {
      if (!embedding_id %in% item_ids) {
        stop(
          "Vector field '", field_id, "' is bound to embedding '", embedding_id,
          "', which is not one of the exported embeddings: ",
          .format_value_list(item_ids), ".",
          call. = FALSE
        )
      }
      label <- paste0(
        "Vector field '", field_id, "' on embedding '", embedding_id, "'"
      )
      inferred <- .infer_vector_shape(bindings[[embedding_id]], label)
      dimension <- dimension_by_id[[embedding_id]]
      if (inferred$dim != dimension) {
        stop(
          label, " must have exactly ", dimension, " columns, got ",
          inferred$dim, ".",
          call. = FALSE
        )
      }
      if (nrow(inferred$dense) != nrow(items[[1L]]$coords)) {
        stop(
          label, " must have exactly ", nrow(items[[1L]]$coords),
          " rows, got ", nrow(inferred$dense), ".",
          call. = FALSE
        )
      }
      .require_finite_float32_source(inferred$dense, label)
      scale_factor <- scale_by_id[[embedding_id]]
      if (
        !is.numeric(scale_factor) ||
          length(scale_factor) != 1L ||
          !is.finite(scale_factor) ||
          scale_factor <= 0
      ) {
        stop(label, " requires a finite positive scale factor.", call. = FALSE)
      }
      bound[[embedding_id]] <- inferred$dense * scale_factor
    }
    # Each field's embeddings follow the item order, which is the order the
    # identity lists them in and the order cellucid-python writes them.
    fields[[field_id]] <- bound[item_ids[item_ids %in% names(bound)]]
  }

  # method = "radix" orders by code point exactly as Python's sorted() does,
  # instead of by the session's collation locale.
  field_ids <- sort(names(fields), method = "radix")
  fields <- fields[field_ids]

  if (is.null(vector_field_default)) {
    if (length(field_ids) != 1L) {
      stop(
        "vector_field_default is required when more than one vector field ",
        "is declared.",
        call. = FALSE
      )
    }
    default_field <- field_ids[[1L]]
  } else {
    if (!vector_field_default %in% field_ids) {
      stop(
        "vector_field_default '", vector_field_default,
        "' does not match an available field: ",
        .format_value_list(field_ids), ".",
        call. = FALSE
      )
    }
    default_field <- vector_field_default
  }
  # The whole identity rule is applied to the field ids, as cellucid-python
  # applies its _require_field_identities() on this axis: the display-text
  # clause is the one this caller can fail.
  .require_field_identities(field_ids, what = "Vector field")
  list(fields = fields, default_field = default_field)
}

# Write the validated fields and return their identity entry and indices.
.write_vector_fields <- function(resolved, items, out_dir, compression) {
  index_by_id <- stats::setNames(
    vapply(items, function(item) item$index, integer(1)),
    vapply(items, function(item) item$id, character(1))
  )
  vectors_dir <- file.path(out_dir, "vectors")
  .dir_create(vectors_dir)
  gz_suffix <- .gz_suffix(compression)
  fields_meta <- list()
  payload_indices <- list()
  field_ids <- names(resolved$fields)
  for (position in seq_along(field_ids)) {
    field_id <- field_ids[[position]]
    payload_index <- .payload_index(position)
    payload_indices[[length(payload_indices) + 1L]] <- payload_index
    files <- list()
    for (embedding_id in names(resolved$fields[[field_id]])) {
      filename <- sprintf(
        "%d_%d.bin",
        payload_index,
        index_by_id[[embedding_id]]
      )
      .write_float32_matrix_row_major(
        file.path(vectors_dir, filename),
        resolved$fields[[field_id]][[embedding_id]],
        compression = compression
      )
      files[[embedding_id]] <- paste0("vectors/", filename, gz_suffix)
    }
    # Key order is the one the output format specification prints and the one
    # cellucid-python emits.
    fields_meta[[field_id]] <- list(
      label = field_id,
      files = files
    )
  }
  list(
    identity = list(
      default_field = resolved$default_field,
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
      name,
      " must be a finite real vector or matrix.",
      call. = FALSE
    )
  }
  storage.mode(dense) <- "double"
  if (anyNA(dense) || any(!is.finite(dense))) {
    stop(
      name, " must contain only finite values. ",
      .describe_non_finite(dense, name),
      call. = FALSE
    )
  }
  list(dim = as.integer(ncol(dense)), dense = unname(dense))
}
