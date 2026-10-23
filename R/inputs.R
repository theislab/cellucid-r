# Checking what cellucid_prepare() was given and turning it into the elements
# of the AnnData store: the obsm arrays, the vector fields, the neighbour
# graph, and the names of the cells. Everything here runs before anything is
# written.

# --- Embeddings -------------------------------------------------------------

# Each embedding becomes obsm[[id]], whole: every column, at its own
# coordinates. The viewer picks the columns it plots and normalizes them.
.resolve_embeddings <- function(embeddings) {
  if (!is.list(embeddings) || is.object(embeddings) || is.data.frame(embeddings)) {
    stop(
      "embeddings must be a named list mapping each obsm key to a matrix, ",
      "for example list(X_umap = umap).",
      call. = FALSE
    )
  }
  if (length(embeddings) == 0L) {
    stop("embeddings must contain at least one embedding.", call. = FALSE)
  }
  ids <- names(embeddings)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids))) {
    stop("embeddings must be a named list: every embedding needs an obsm key.", call. = FALSE)
  }
  for (id in ids) {
    .validate_element_key(id, "Embedding")
  }
  .require_distinct_element_keys(ids, "Embedding")
  matrices <- list()
  n_cells <- NULL
  for (id in ids) {
    matrix <- .as_cell_matrix(embeddings[[id]], paste0("Embedding '", id, "'"))
    if (ncol(matrix) < 1L) {
      stop("Embedding '", id, "' has no columns.", call. = FALSE)
    }
    if (is.null(n_cells)) {
      n_cells <- nrow(matrix)
      first <- id
    } else if (nrow(matrix) != n_cells) {
      stop(
        "All embeddings must have one row per cell: '", first, "' has ",
        n_cells, " rows and '", id, "' has ", nrow(matrix), ".",
        call. = FALSE
      )
    }
    matrices[[id]] <- matrix
  }
  if (n_cells < 1L) {
    stop("Embeddings must contain at least one cell.", call. = FALSE)
  }
  list(matrices = matrices, n_cells = as.integer(n_cells))
}

# embedding_defaults names the embedding each dimension shows first: an obsm
# key, plotted by its first columns, or list(key = , columns = ) naming the
# 1-based columns to plot. The store records the columns 0-based, as AnnData
# indexes them.
.resolve_embedding_defaults <- function(embedding_defaults, matrices) {
  if (is.null(embedding_defaults)) {
    return(list())
  }
  if (
    !(is.character(embedding_defaults) || is.list(embedding_defaults)) ||
      is.object(embedding_defaults) ||
      is.null(names(embedding_defaults)) ||
      anyNA(names(embedding_defaults))
  ) {
    stop(
      "embedding_defaults must map \"1d\", \"2d\", or \"3d\" to an obsm key, ",
      "for example c(\"2d\" = \"X_umap\").",
      call. = FALSE
    )
  }
  slots <- names(embedding_defaults)
  if (anyDuplicated(slots)) {
    stop("embedding_defaults names one dimension more than once.", call. = FALSE)
  }
  resolved <- list()
  for (slot in slots) {
    if (!slot %in% c("1d", "2d", "3d")) {
      stop("embedding_defaults key '", slot, "' must be \"1d\", \"2d\", or \"3d\".", call. = FALSE)
    }
    dimension <- as.integer(substr(slot, 1L, 1L))
    spec <- embedding_defaults[[slot]]
    label <- paste0("embedding_defaults[[\"", slot, "\"]]")
    columns <- NULL
    if (is.list(spec)) {
      if (!setequal(names(spec), c("key", "columns")) || length(spec) != 2L) {
        stop(label, " must be an obsm key or list(key = , columns = ).", call. = FALSE)
      }
      key <- spec$key
      columns <- spec$columns
    } else {
      key <- spec
    }
    if (!is.character(key) || length(key) != 1L || is.na(key) || !key %in% names(matrices)) {
      stop(
        label, " must name one of the embeddings: ",
        .format_value_list(names(matrices)), ".",
        call. = FALSE
      )
    }
    width <- ncol(matrices[[key]])
    if (is.null(columns)) {
      if (width < dimension) {
        stop(
          label, " names '", key, "', which has ", width, " ",
          .plural(width, "column", "columns"), "; a ", slot, " view needs ",
          dimension, ".",
          call. = FALSE
        )
      }
      resolved[[paste0("embedding", slot)]] <- .uns_string(key)
      next
    }
    if (
      !is.numeric(columns) || is.object(columns) || anyNA(columns) ||
        length(columns) != dimension || any(columns != floor(columns)) ||
        any(columns < 1) || any(columns > width) || anyDuplicated(columns)
    ) {
      stop(
        label, "$columns must name ", dimension, " different 1-based ",
        .plural(dimension, "column", "columns"), " of '", key,
        "', which has columns 1 to ", width, ".",
        call. = FALSE
      )
    }
    resolved[[paste0("embedding", slot)]] <- list(
      key = .uns_string(key),
      columns = .uns_integers(as.integer(columns) - 1L)
    )
  }
  resolved
}

# --- Vector fields ----------------------------------------------------------

# A vector field is stored under scVelo's name for it: bound to the embedding
# 'X_<basis>', the field '<field>' becomes obsm[['<field>_<basis>']]
# (list(velocity = list(X_umap = v)) writes 'velocity_umap', on X_branches
# 'velocity_branches'), which the viewer draws as the field '<field>' over that
# embedding by the vector-field key rule (.parse_vector_field_key). Each array
# has one column per column of its embedding. Every key must read back as its
# own field and embedding, so a field named like an embedding ('X_...'), or a
# key another embedding's longer basis claims, is refused.
.resolve_vector_fields <- function(vector_fields, matrices, n_cells) {
  if (is.null(vector_fields)) {
    return(list(arrays = list(), fields = character()))
  }
  if (!is.list(vector_fields) || is.object(vector_fields) || is.data.frame(vector_fields)) {
    stop(
      "vector_fields must map each field to a named list of embeddings and ",
      "arrays, for example list(velocity = list(X_umap = velocity_umap)).",
      call. = FALSE
    )
  }
  fields <- names(vector_fields)
  if (length(vector_fields) > 0L && (is.null(fields) || anyNA(fields) || any(!nzchar(fields)))) {
    stop("vector_fields must name every field.", call. = FALSE)
  }
  if (anyDuplicated(fields)) {
    stop("vector_fields names one field more than once.", call. = FALSE)
  }
  arrays <- list()
  named <- list()
  for (field in fields) {
    bindings <- vector_fields[[field]]
    if (!is.list(bindings) || is.object(bindings) || is.data.frame(bindings) ||
        length(bindings) == 0L || is.null(names(bindings))) {
      stop(
        "Vector field '", field, "' must map embeddings to arrays, for ",
        "example list(", field, " = list(X_umap = vectors)).",
        call. = FALSE
      )
    }
    for (embedding in names(bindings)) {
      label <- paste0("Vector field '", field, "' on embedding '", embedding, "'")
      if (!embedding %in% names(matrices)) {
        stop(
          label, ": '", embedding, "' is not one of the embeddings: ",
          .format_value_list(names(matrices)), ".",
          call. = FALSE
        )
      }
      if (!startsWith(embedding, .embedding_key_prefix) || embedding == .embedding_key_prefix) {
        stop(
          label, ": the viewer binds a vector field to the embedding 'X_<basis>' ",
          "it is named for ('velocity_umap' to 'X_umap'), so the embedding needs ",
          "an 'X_<basis>' key. Rename '", embedding, "', for example to 'X_",
          embedding, "'.",
          call. = FALSE
        )
      }
      values <- .as_cell_matrix(bindings[[embedding]], label)
      width <- ncol(matrices[[embedding]])
      if (ncol(values) != width) {
        stop(
          label, " has ", ncol(values), " columns, and '", embedding, "' has ",
          width, ": a vector field has one component per embedding column.",
          call. = FALSE
        )
      }
      if (nrow(values) != n_cells) {
        stop(label, " must have one row per cell (", n_cells, "), got ", nrow(values), ".", call. = FALSE)
      }
      basis <- substring(embedding, nchar(.embedding_key_prefix) + 1L)
      key <- paste0(field, "_", basis)
      if (!is.null(arrays[[key]])) {
        stop(
          label, " and vector field '", named[[key]]$field, "' on '",
          named[[key]]$embedding, "' would both be stored as obsm[['", key,
          "']]. Rename one field.",
          call. = FALSE
        )
      }
      .validate_element_key(key, "Vector field obsm key")
      arrays[[key]] <- values
      named[[key]] <- list(field = field, embedding = embedding)
    }
  }
  collisions <- intersect(tolower(names(arrays)), tolower(names(matrices)))
  if (length(collisions) > 0L) {
    stop(
      "A vector field would be stored as obsm[['", collisions[[1L]], "']], ",
      "which is also the key of an embedding. Rename the field or the embedding.",
      call. = FALSE
    )
  }
  obsm_keys <- c(names(matrices), names(arrays))
  for (key in names(named)) {
    field <- named[[key]]$field
    embedding <- named[[key]]$embedding
    parsed <- .parse_vector_field_key(key, obsm_keys)
    if (is.null(parsed) || !identical(parsed$field, field) ||
        !identical(paste0(.embedding_key_prefix, parsed$basis), embedding)) {
      reads_as <- if (is.null(parsed)) {
        "an embedding: a key that starts with 'X_' is never a vector field"
      } else {
        paste0("the field '", parsed$field, "' on '", .embedding_key_prefix, parsed$basis, "'")
      }
      stop(
        "Vector field '", field, "' on embedding '", embedding, "' would be ",
        "stored as obsm[['", key, "']], which the viewer reads as ", reads_as,
        ". Rename the field.",
        call. = FALSE
      )
    }
  }
  list(arrays = arrays, fields = fields, bound = named)
}

# The web app's best-guess ranking of obsm keys (bestGuessTier in
# data/data-mapping.js): the first needle a lower-cased key contains decides
# its tier, a lower tier is a better guess, and a key matching none ranks 7.
.best_guess_tiers <- list(
  list(0L, "umap"),
  list(1L, c("tsne", "t-sne", "t_sne")),
  list(2L, c("draw_graph", "force")),
  list(3L, "phate"),
  list(4L, c("pacmap", "trimap")),
  list(5L, "spatial"),
  list(6L, "diffmap"),
  list(8L, c("pca", "scvi", "latent", "harmony", "x_emb"))
)

# Keys in best-guess order: by tier, then by code point.
.best_guess_order <- function(keys) {
  tier <- vapply(keys, function(key) {
    name <- tolower(key)
    for (entry in .best_guess_tiers) {
      if (any(vapply(entry[[2L]], function(needle) grepl(needle, name, fixed = TRUE), logical(1L)))) {
        return(entry[[1L]])
      }
    }
    7L
  }, integer(1L))
  # Radix order compares strings byte by byte: code-point order for UTF-8.
  keys[order(tier, keys, method = "radix")]
}

# The embeddings the web app shows a store with when it opens
# (resolveEmbeddingSlots and defaultEnabledEmbeddings, from metadata): the
# slots the defaults name over the default rule (X_umap_1d/2d/3d, else a
# 1-3 column X_umap); when nothing names a slot and the rule fills none, the
# best guess (per dimension the first key of that exact width, else the first
# wider key in 2D); and every other 1-3 column embedding. The latent space is
# never guessed or shown beside.
.shown_embedding_keys <- function(matrices, named_slots, latent_key) {
  widths <- vapply(matrices, ncol, integer(1L))
  candidates <- widths[names(widths) != if (is.null(latent_key)) "" else latent_key]
  dimensional <- intersect(c("X_umap_1d", "X_umap_2d", "X_umap_3d"), names(widths))
  slots <- if (length(dimensional) > 0L) {
    dimensional
  } else if ("X_umap" %in% names(widths) && widths[["X_umap"]] %in% 1:3) {
    "X_umap"
  } else {
    character()
  }
  if (length(named_slots) > 0L) {
    slots <- c(slots, unlist(named_slots, use.names = FALSE))
  } else if (length(slots) == 0L) {
    for (dimension in 3:1) {
      exact <- .best_guess_order(names(candidates)[candidates == dimension])
      if (length(exact) > 0L) slots <- c(slots, exact[[1L]])
    }
    if (length(slots) == 0L) {
      wide <- .best_guess_order(names(candidates)[candidates > 3L])
      if (length(wide) > 0L) slots <- wide[[1L]]
    }
  }
  unique(c(slots, names(candidates)[candidates <= 3L]))
}

# --- The neighbour graph ----------------------------------------------------

# The graph is stored as obsp['connectivities'], a CSR matrix. It must be a
# square cells x cells matrix, exactly symmetric, with a zero diagonal and
# finite non-negative weights; a stored zero is not an edge and is dropped.
# For a symmetric matrix the column-compressed arrays of a dgCMatrix are the
# row-compressed arrays of the same matrix, so they are written as they are.
.resolve_connectivities <- function(connectivities, n_cells) {
  if (inherits(connectivities, "Matrix")) {
    graph <- methods::as(
      methods::as(methods::as(connectivities, "dMatrix"), "generalMatrix"),
      "CsparseMatrix"
    )
  } else if (
    is.matrix(connectivities) &&
      (is.numeric(connectivities) || is.logical(connectivities)) &&
      !is.object(connectivities)
  ) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("The Matrix package is required to store connectivities.", call. = FALSE)
    }
    storage.mode(connectivities) <- "double"
    graph <- methods::as(Matrix::Matrix(connectivities, sparse = TRUE), "generalMatrix")
    graph <- methods::as(graph, "CsparseMatrix")
  } else {
    stop(
      "connectivities must be a real numeric or logical matrix or a ",
      "Matrix::Matrix sparse matrix.",
      call. = FALSE
    )
  }
  if (nrow(graph) != n_cells || ncol(graph) != n_cells) {
    stop(
      "connectivities must be ", n_cells, " x ", n_cells, " (one row and one ",
      "column per cell), got ", nrow(graph), " x ", ncol(graph), ".",
      call. = FALSE
    )
  }
  weights <- graph@x
  if (anyNA(weights) || any(!is.finite(weights))) {
    stop("Connectivity weights must all be finite.", call. = FALSE)
  }
  if (any(weights < 0)) {
    stop("Connectivity weights must all be non-negative.", call. = FALSE)
  }
  graph <- Matrix::drop0(graph)
  if (any(Matrix::diag(graph) != 0)) {
    stop("Connectivity diagonal values must all be exactly zero.", call. = FALSE)
  }
  if (!isTRUE(Matrix::isSymmetric(graph, tol = 0))) {
    stop("Connectivity matrix must be exactly symmetric.", call. = FALSE)
  }
  .require_finite_float32(graph@x, "connectivities", positions = FALSE)
  graph
}

# --- Cell and gene names ----------------------------------------------------

.supplied_row_names <- function(x) {
  if (is.data.frame(x)) {
    if (.row_names_info(x) < 0L) NULL else rownames(x)
  } else {
    rownames(x)
  }
}

# The obs index: the cell names every input agrees on. Row names may come
# from obs, from gene_expression and from any embedding; where several are
# given they must be the same names in the same order, because a different
# order means the rows describe different cells. With no names at all the
# cells are named "0" to "n - 1", as anndata names them.
.resolve_cell_names <- function(sources, n_cells) {
  named <- Filter(Negate(is.null), sources)
  if (length(named) == 0L) {
    return(as.character(seq_len(n_cells) - 1L))
  }
  reference <- named[[1L]]
  for (source in names(named)[-1L]) {
    names_here <- named[[source]]
    if (!identical(as.character(names_here), as.character(reference))) {
      differs <- which(as.character(names_here) != as.character(reference))[[1L]]
      stop(
        "The row names of ", source, " and ", names(named)[[1L]], " disagree ",
        "at row ", differs, " ('", names_here[[differs]], "' and '",
        reference[[differs]], "'). Every input must list the same cells in ",
        "the same order.",
        call. = FALSE
      )
    }
  }
  if (anyNA(reference)) {
    stop("Cell names must not be missing (", names(named)[[1L]], ").", call. = FALSE)
  }
  unique_names <- .make_unique_names(as.character(reference))
  if (length(unique_names$renamed) > 0L) {
    message("Cell names: ", .unique_names_note(0, unique_names$renamed, NULL))
  }
  unique_names$names
}
