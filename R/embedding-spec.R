# The embeddings an export carries. Any number of embeddings per dimension can
# be exported; each one is an item with an id, a display label, a dimension,
# and the source columns it was read from. Items are ordered by dimension and
# then by id in code-point order, and an item's position in that order is its
# payload index: its coordinates are written to embeddings/{index}_{dim}d.bin.
# cellucid-python's prepare_data/_embeddings.py resolves the same items in the
# same order, so both writers name every file identically.

# One embedding to export, with the columns it is read from. Columns are
# 1-based here, as everywhere in R, and are published 0-based.
cellucid_embedding <- function(
    values,
    columns = NULL,
    label = NULL,
    source_key = NULL
) {
  dense <- .embedding_values_matrix(values, "values")
  width <- ncol(dense)
  if (is.null(columns)) {
    if (!width %in% c(1L, 2L, 3L)) {
      stop(
        "values has ", width, " columns, and an embedding plots 1, 2, or 3 of ",
        "them. Pass columns = c(...) naming the 1-based columns to plot, for ",
        "example columns = c(1, 2).",
        call. = FALSE
      )
    }
  } else {
    columns <- .validate_embedding_columns(columns, width, "columns")
  }
  if (!is.null(label)) {
    label <- .validate_display_text(label, "label")
  }
  if (!is.null(source_key)) {
    source_key <- .validate_string(source_key, "source_key")
    if (!nzchar(source_key)) {
      stop("source_key must be NULL or one non-empty string.", call. = FALSE)
    }
  }
  structure(
    list(
      values = dense,
      columns = columns,
      label = label,
      source_key = source_key
    ),
    class = "cellucid_embedding"
  )
}

# Coordinates as a double matrix: a base numeric matrix, a numeric data frame,
# a Matrix, or a plain numeric vector, which is the one column it declares.
.embedding_values_matrix <- function(x, label) {
  if (is.numeric(x) && !is.object(x) && is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  } else if (.is_numeric_data_frame(x)) {
    x <- as.matrix(x)
  } else if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
  } else if (!is.matrix(x) || !is.numeric(x) || is.object(x)) {
    stop(
      label, " must be a finite real numeric matrix or vector.",
      call. = FALSE
    )
  } else {
    x <- unclass(x)
  }
  storage.mode(x) <- "double"
  x
}

# The 1-based columns one embedding reads. The checks run in the order
# cellucid-python's validate_embedding_columns() runs them, so one input is
# refused for one reason by both writers.
.validate_embedding_columns <- function(columns, width, label) {
  if (
    !is.numeric(columns) ||
      is.object(columns) ||
      !is.null(dim(columns)) ||
      anyNA(columns)
  ) {
    stop(label, " must be a numeric vector of 1-based column numbers.", call. = FALSE)
  }
  dimension <- length(columns)
  if (!dimension %in% c(1L, 2L, 3L)) {
    stop(
      label, " must name 1, 2, or 3 columns, one per axis; got ", dimension, ".",
      call. = FALSE
    )
  }
  if (width < dimension) {
    stop(
      "A ", dimension, "D embedding needs at least ", dimension,
      " columns; this array has ", width, ".",
      call. = FALSE
    )
  }
  if (any(columns != floor(columns))) {
    stop(label, " must hold whole column numbers.", call. = FALSE)
  }
  if (any(columns < 1) || any(columns > width)) {
    stop(
      label, " names a column that does not exist; the array has columns 1 to ",
      width, ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(columns)) {
    stop(label, ": each axis needs a different column.", call. = FALSE)
  }
  as.integer(columns)
}

.resolve_one_embedding <- function(embedding_id, value) {
  label <- paste0("Embedding '", embedding_id, "'")
  if (inherits(value, "cellucid_embedding")) {
    spec <- value
  } else {
    if (is.null(value)) {
      stop(label, " has no values.", call. = FALSE)
    }
    spec <- list(
      values = .embedding_values_matrix(value, label),
      columns = NULL,
      label = NULL,
      source_key = NULL
    )
    width <- ncol(spec$values)
    if (!width %in% c(1L, 2L, 3L)) {
      stop(
        label, " has ", width, " columns, and an embedding plots 1, 2, or 3 ",
        "of them. Pass cellucid_embedding(values, columns = c(...)) naming the ",
        "1-based columns to plot, for example columns = c(1, 2).",
        call. = FALSE
      )
    }
  }
  columns <- if (is.null(spec$columns)) {
    seq_len(ncol(spec$values))
  } else {
    spec$columns
  }
  selected <- spec$values[, columns, drop = FALSE]
  if (anyNA(selected) || any(!is.finite(selected))) {
    stop(
      label, " must contain only finite values. ",
      .describe_non_finite(selected, label),
      call. = FALSE
    )
  }
  # Before normalization, because normalization hides it. Centring and
  # scaling pull every coordinate into roughly [-1, 1], so a source
  # coordinate of 1e300 or 1e-320 reaches the writer already inside the
  # float32 range and would be published as a number the caller never
  # supplied. cellucid-python gates the same values at the same point.
  .require_finite_float32_source(selected, label)
  list(
    coords = unname(selected),
    columns = as.integer(columns) - 1L,
    label = if (is.null(spec$label)) embedding_id else spec$label,
    source_key = spec$source_key
  )
}

# Validate every embedding and the per-dimension defaults, then normalize each
# item on its own. Everything is checked before anything is written.
.resolve_embeddings <- function(embeddings, embedding_defaults) {
  if (
    !is.list(embeddings) ||
      is.object(embeddings) ||
      is.data.frame(embeddings)
  ) {
    stop(
      "embeddings must be a named list mapping each embedding id to a matrix ",
      "or a cellucid_embedding(), for example list(X_umap = umap).",
      call. = FALSE
    )
  }
  if (length(embeddings) == 0L) {
    stop("embeddings must contain at least one embedding.", call. = FALSE)
  }
  ids <- names(embeddings)
  if (is.null(ids)) {
    stop("embeddings must be a named list: every embedding needs an id.", call. = FALSE)
  }
  .require_field_identities(ids, what = "Embedding")

  resolved <- list()
  n_cells <- NULL
  first_id <- NULL
  for (embedding_id in ids) {
    item <- .resolve_one_embedding(embedding_id, embeddings[[embedding_id]])
    rows <- nrow(item$coords)
    if (is.null(n_cells)) {
      n_cells <- rows
      first_id <- embedding_id
    } else if (rows != n_cells) {
      stop(
        "All embeddings must have the same number of cells: '", first_id,
        "' has ", n_cells, " and '", embedding_id, "' has ", rows, ".",
        call. = FALSE
      )
    }
    resolved[[embedding_id]] <- item
  }
  if (is.null(n_cells) || n_cells < 1L || n_cells > .Machine$integer.max) {
    stop(
      "Embeddings must contain between 1 and 2,147,483,647 cells.",
      call. = FALSE
    )
  }

  dimensions <- vapply(resolved, function(item) length(item$columns), integer(1))
  # Code-point order within a dimension: method = "radix" sorts by bytes, which
  # for UTF-8 is code-point order, exactly as Python's sorted() does.
  order_ids <- character(0)
  ids_by_dimension <- list()
  for (dimension in sort(unique(dimensions))) {
    in_dimension <- sort(ids[dimensions[ids] == dimension], method = "radix")
    ids_by_dimension[[as.character(dimension)]] <- in_dimension
    order_ids <- c(order_ids, in_dimension)
  }
  defaults <- .resolve_embedding_defaults(embedding_defaults, ids_by_dimension)

  items <- list()
  for (index in seq_along(order_ids)) {
    embedding_id <- order_ids[[index]]
    item <- resolved[[embedding_id]]
    normalized <- .normalize_embedding(item$coords)
    items[[index]] <- list(
      id = embedding_id,
      index = index - 1L,
      label = item$label,
      dimension = length(item$columns),
      columns = item$columns,
      source_key = item$source_key,
      coords = normalized$coords,
      scale_factor = normalized$info$scale_factor,
      original_range = normalized$info$original_range
    )
  }
  list(items = items, defaults = defaults, n_cells = as.integer(n_cells))
}

.resolve_embedding_defaults <- function(embedding_defaults, ids_by_dimension) {
  requested <- list()
  if (!is.null(embedding_defaults)) {
    if (
      !(is.character(embedding_defaults) || is.list(embedding_defaults)) ||
        is.object(embedding_defaults) ||
        is.null(names(embedding_defaults)) ||
        anyNA(names(embedding_defaults))
    ) {
      stop(
        "embedding_defaults must map \"1d\", \"2d\", or \"3d\" to an embedding ",
        "id, for example c(\"2d\" = \"X_umap\").",
        call. = FALSE
      )
    }
    for (slot in names(embedding_defaults)) {
      if (!slot %in% c("1d", "2d", "3d")) {
        stop(
          "embedding_defaults key '", slot, "' must be \"1d\", \"2d\", or \"3d\".",
          call. = FALSE
        )
      }
      dimension <- substr(slot, 1L, 1L)
      if (is.null(ids_by_dimension[[dimension]])) {
        stop(
          "embedding_defaults names a ", slot, " default, and no ", dimension,
          "D embedding was supplied.",
          call. = FALSE
        )
      }
      embedding_id <- embedding_defaults[[slot]]
      if (
        !is.character(embedding_id) ||
          length(embedding_id) != 1L ||
          is.na(embedding_id) ||
          !embedding_id %in% ids_by_dimension[[dimension]]
      ) {
        stop(
          "embedding_defaults[[\"", slot, "\"]] is not a ", dimension,
          "D embedding. The ", dimension, "D embeddings are: ",
          .format_value_list(ids_by_dimension[[dimension]]), ".",
          call. = FALSE
        )
      }
      requested[[slot]] <- embedding_id
    }
  }
  defaults <- list()
  for (dimension in names(ids_by_dimension)) {
    slot <- paste0(dimension, "d")
    candidates <- ids_by_dimension[[dimension]]
    if (!is.null(requested[[slot]])) {
      defaults[[slot]] <- requested[[slot]]
    } else if (length(candidates) == 1L) {
      defaults[[slot]] <- candidates[[1L]]
    } else {
      stop(
        length(candidates), " ", dimension, "D embeddings were supplied (",
        .format_value_list(candidates), "), so embedding_defaults must name the ",
        "one shown first, for example embedding_defaults = c(\"", slot,
        "\" = \"", candidates[[1L]], "\").",
        call. = FALSE
      )
    }
  }
  defaults
}
