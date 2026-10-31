# Checking what cellucid_prepare() was given and turning it into the elements
# of the AnnData store: the obsm arrays, the vector fields, the neighbour
# graph, and the names of the cells. Everything here runs before anything is
# written.

# --- Embeddings -------------------------------------------------------------

# Each embedding becomes obsm[[id]], whole: every column, at its own
# coordinates. The viewer picks the columns it plots and normalizes them.
.resolve_embeddings <- function(embeddings) {
  if (!is.list(embeddings) || is.object(embeddings) || is.data.frame(embeddings)) {
    .abort(
      "embeddings must be a named list mapping each obsm key to a matrix, ",
      "for example list(X_umap = umap)."
    )
  }
  if (length(embeddings) == 0L) {
    .abort("embeddings must contain at least one embedding.")
  }
  ids <- names(embeddings)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids))) {
    .abort("embeddings must be a named list: every embedding needs an obsm key.")
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
      .abort("Embedding '", id, "' has no columns.")
    }
    if (is.null(n_cells)) {
      n_cells <- nrow(matrix)
      first <- id
    } else if (nrow(matrix) != n_cells) {
      .abort(
        "All embeddings must have one row per cell: '", first, "' has ",
        n_cells, " ", .plural(n_cells, "row"), " and '", id, "' has ", nrow(matrix), "."
      )
    }
    matrices[[id]] <- matrix
  }
  if (n_cells < 1L) {
    .abort("Embeddings must contain at least one cell.")
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
    .abort(
      "embedding_defaults must map \"1d\", \"2d\", or \"3d\" to an obsm key, ",
      "for example c(\"2d\" = \"X_umap\")."
    )
  }
  slots <- names(embedding_defaults)
  if (anyDuplicated(slots)) {
    .abort("embedding_defaults names one dimension more than once.")
  }
  resolved <- list()
  for (slot in slots) {
    if (!slot %in% c("1d", "2d", "3d")) {
      .abort("embedding_defaults key '", slot, "' must be \"1d\", \"2d\", or \"3d\".")
    }
    dimension <- as.integer(substr(slot, 1L, 1L))
    spec <- embedding_defaults[[slot]]
    label <- paste0("embedding_defaults[[\"", slot, "\"]]")
    columns <- NULL
    if (is.list(spec)) {
      if (!setequal(names(spec), c("key", "columns")) || length(spec) != 2L) {
        .abort(label, " must be an obsm key or list(key = , columns = ).")
      }
      key <- spec$key
      columns <- spec$columns
    } else {
      key <- spec
    }
    if (!is.character(key) || length(key) != 1L || is.na(key) || !key %in% names(matrices)) {
      .abort(
        label, " must name one of the embeddings: ",
        .format_value_list(names(matrices)), "."
      )
    }
    width <- ncol(matrices[[key]])
    if (is.null(columns)) {
      if (width < dimension) {
        .abort(
          label, " names '", key, "', which has ", width, " ",
          .plural(width, "column", "columns"), "; a ", slot, " view needs ",
          dimension, "."
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
      .abort(
        label, "$columns must name ", dimension, " different 1-based ",
        .plural(dimension, "column", "columns"), " of '", key,
        "', which has columns 1 to ", width, "."
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
    .abort(
      "vector_fields must map each field to a named list of embeddings and ",
      "arrays, for example list(velocity = list(X_umap = velocity_umap))."
    )
  }
  fields <- names(vector_fields)
  if (length(vector_fields) > 0L && (is.null(fields) || anyNA(fields) || any(!nzchar(fields)))) {
    .abort("vector_fields must name every field.")
  }
  if (anyDuplicated(fields)) {
    .abort("vector_fields names one field more than once.")
  }
  arrays <- list()
  named <- list()
  for (field in fields) {
    bindings <- vector_fields[[field]]
    if (!is.list(bindings) || is.object(bindings) || is.data.frame(bindings) ||
        length(bindings) == 0L || is.null(names(bindings))) {
      .abort(
        "Vector field '", field, "' must map embeddings to arrays, for ",
        "example list(", field, " = list(X_umap = vectors))."
      )
    }
    for (embedding in names(bindings)) {
      label <- paste0("Vector field '", field, "' on embedding '", embedding, "'")
      if (!embedding %in% names(matrices)) {
        .abort(
          label, ": '", embedding, "' is not one of the embeddings: ",
          .format_value_list(names(matrices)), "."
        )
      }
      if (!startsWith(embedding, .embedding_key_prefix) || embedding == .embedding_key_prefix) {
        .abort(
          label, ": the viewer binds a vector field to the embedding 'X_<basis>' ",
          "it is named for ('velocity_umap' to 'X_umap'), so the embedding needs ",
          "an 'X_<basis>' key. Rename '", embedding, "', for example to 'X_",
          embedding, "'."
        )
      }
      values <- .as_cell_matrix(bindings[[embedding]], label)
      width <- ncol(matrices[[embedding]])
      if (ncol(values) != width) {
        .abort(
          label, " has ", ncol(values), " ", .plural(ncol(values), "column"), ", and '", embedding, "' has ",
          width, ": a vector field has one component per embedding column."
        )
      }
      if (nrow(values) != n_cells) {
        .abort(label, " must have one row per cell (", n_cells, "), got ", nrow(values), ".")
      }
      basis <- substring(embedding, nchar(.embedding_key_prefix) + 1L)
      key <- paste0(field, "_", basis)
      if (!is.null(arrays[[key]])) {
        .abort(
          label, " and vector field '", named[[key]]$field, "' on '",
          named[[key]]$embedding, "' would both be stored as obsm[['", key,
          "']]. Rename one field."
        )
      }
      .validate_element_key(key, "Vector field obsm key")
      arrays[[key]] <- values
      named[[key]] <- list(field = field, embedding = embedding)
    }
  }
  collisions <- intersect(tolower(names(arrays)), tolower(names(matrices)))
  if (length(collisions) > 0L) {
    .abort(
      "A vector field would be stored as obsm[['", collisions[[1L]], "']], ",
      "which is also the key of an embedding. Rename the field or the embedding."
    )
  }
  obsm_keys <- c(names(matrices), names(arrays))
  for (key in names(named)) {
    field <- named[[key]]$field
    embedding <- named[[key]]$embedding
    parsed <- .parse_vector_field_key(key, obsm_keys)
    if (is.null(parsed) || !identical(parsed$field, field) ||
        !identical(parsed$embedding, embedding)) {
      reads_as <- if (is.null(parsed)) {
        "an embedding: a key that starts with 'X_' is never a vector field"
      } else {
        paste0("the field '", parsed$field, "' on '", parsed$embedding, "'")
      }
      .abort(
        "Vector field '", field, "' on embedding '", embedding, "' would be ",
        "stored as obsm[['", key, "']], which the viewer reads as ", reads_as,
        ". Rename the field."
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

# Keys in best-guess order: by tier (of a key's own name, .obsm_namespace),
# then namespace (the dataset's own keys, then the modality it opens with,
# then the other modalities), then code point.
.best_guess_order <- function(keys, modality = NULL) {
  tier <- vapply(keys, function(key) {
    name <- tolower(.obsm_namespace(key)$local)
    for (entry in .best_guess_tiers) {
      if (any(vapply(entry[[2L]], function(needle) grepl(needle, name, fixed = TRUE), logical(1L)))) {
        return(entry[[1L]])
      }
    }
    7L
  }, integer(1L))
  namespace <- vapply(keys, function(key) {
    prefix <- .obsm_namespace(key)$prefix
    if (identical(prefix, "")) {
      0L
    } else if (!is.null(modality) && identical(prefix, paste0(modality, ":"))) {
      1L
    } else {
      2L
    }
  }, integer(1L))
  # Radix order compares strings byte by byte: code-point order for UTF-8.
  keys[order(tier, namespace, enc2utf8(keys), method = "radix")]
}

# The conventional names of one UMAP, in the order the default rule takes
# them (UMAP_KEYS in the web app's data-mapping.js): Scanpy's X_umap (also a
# Seurat object's umap reduction), a SingleCellExperiment's UMAP reducedDim
# (scater's runUMAP()), then a Seurat object's weighted-nearest-neighbour and
# Azimuth reference UMAPs. Reductions and reducedDims are obsm['X_<name>']
# with the name's case kept.
.umap_keys <- c("X_umap", "X_UMAP", "X_wnn.umap", "X_ref.umap")

# The embeddings the web app shows a store with when it opens
# (resolveEmbeddingSlots and defaultEnabledEmbeddings, from metadata): the
# slots the defaults name over the default rule (X_umap_1d/2d/3d, else the
# first 1-3 column one of .umap_keys); when nothing names a slot and the rule fills none, the
# best guess (per dimension the first key of that exact width, else the first
# wider key in 2D); and every other 1-3 column embedding. The latent space is
# never guessed or shown beside.
.shown_embedding_keys <- function(matrices, named_slots, latent_key, modality = NULL) {
  widths <- vapply(matrices, ncol, integer(1L))
  candidates <- widths[names(widths) != if (is.null(latent_key)) "" else latent_key]
  slots <- character()
  # The dataset's own names, then (several modalities) the same names in the
  # namespace of the modality it opens with.
  for (prefix in if (is.null(modality)) "" else c("", paste0(modality, ":"))) {
    dimensional <- intersect(paste0(prefix, c("X_umap_1d", "X_umap_2d", "X_umap_3d")), names(widths))
    slots <- if (length(dimensional) > 0L) {
      dimensional
    } else {
      single <- paste0(prefix, .umap_keys)
      single <- single[single %in% names(widths)]
      single <- single[widths[single] %in% 1:3]
      single[seq_len(min(1L, length(single)))]
    }
    if (length(slots) > 0L) break
  }
  if (length(named_slots) > 0L) {
    slots <- c(slots, unlist(named_slots, use.names = FALSE))
  } else if (length(slots) == 0L) {
    for (dimension in 3:1) {
      exact <- .best_guess_order(names(candidates)[candidates == dimension], modality)
      if (length(exact) > 0L) slots <- c(slots, exact[[1L]])
    }
    if (length(slots) == 0L) {
      wide <- .best_guess_order(names(candidates)[candidates > 3L], modality)
      if (length(wide) > 0L) slots <- wide[[1L]]
    }
  }
  unique(c(slots, names(candidates)[candidates <= 3L]))
}

# --- The neighbour graph ----------------------------------------------------

# The neighbour graph a dataset opens with, of its graph names (the web app's
# defaultGraph and cellucid-python's default_graph, pinned by the shared
# conformance fixture): Scanpy's 'connectivities'; else Seurat's
# shared-nearest-neighbour graph of the modality (assay) it opens with
# ('<modality>_snn'); else the first other shared-nearest-neighbour graph
# ('<assay>_snn', or 'SNN', in any case) in code-point order; else NULL.
# Seurat's '<assay>_nn' (directed) and Scanpy's 'distances' are never the
# default. .graph_candidates lists them all in that order.
.graph_candidates <- function(keys, modality = NULL) {
  keys <- as.character(keys)
  shared <- keys[tolower(keys) == "snn" | (endsWith(tolower(keys), "_snn") & nchar(keys) > 4L)]
  shared <- shared[order(enc2utf8(shared), method = "radix")]
  own <- if (is.null(modality)) character() else intersect(paste0(modality, "_snn"), keys)
  unique(c(intersect("connectivities", keys), own, shared))
}

.default_graph <- function(keys, modality = NULL) {
  candidates <- .graph_candidates(keys, modality)
  if (length(candidates) == 0L) NULL else candidates[[1L]]
}

# The graph is stored as obsp['connectivities'], a CSR matrix. It must be a
# square cells x cells matrix, exactly symmetric, with finite non-negative
# weights; a stored zero is not an edge and is dropped, and so is the
# diagonal: a self-loop is no edge anyone can draw (Seurat's FindNeighbors
# counts each cell among its own neighbours, so its '<assay>_snn' has 1 on
# the diagonal), as the web app and cellucid-python drop it. A directed graph
# (Seurat's '<assay>_nn', Scanpy's 'distances') is refused.
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
      .abort("The Matrix package is required to store connectivities.")
    }
    storage.mode(connectivities) <- "double"
    graph <- methods::as(Matrix::Matrix(connectivities, sparse = TRUE), "generalMatrix")
    graph <- methods::as(graph, "CsparseMatrix")
  } else {
    .abort(
      "connectivities must be a real numeric or logical matrix or a ",
      "Matrix::Matrix sparse matrix."
    )
  }
  if (nrow(graph) != n_cells || ncol(graph) != n_cells) {
    .abort(
      "connectivities must be ", n_cells, " x ", n_cells, " (one row and one ",
      "column per cell), got ", nrow(graph), " x ", ncol(graph), "."
    )
  }
  weights <- graph@x
  if (anyNA(weights) || any(!is.finite(weights))) {
    .abort("Connectivity weights must all be finite.")
  }
  if (any(weights < 0)) {
    .abort("Connectivity weights must all be non-negative.")
  }
  Matrix::diag(graph) <- 0
  graph <- Matrix::drop0(graph)
  if (!isTRUE(Matrix::isSymmetric(graph, tol = 0))) {
    .abort(
      "Connectivity matrix must be exactly symmetric: a directed neighbour graph ",
      "(Seurat's '<assay>_nn', Scanpy's 'distances') is not one the viewer draws; ",
      "pass a symmetric one (obsp[['connectivities']], or Seurat's '<assay>_snn')."
    )
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
      .abort(
        "The row names of ", source, " and ", names(named)[[1L]], " disagree ",
        "at row ", differs, " ('", names_here[[differs]], "' and '",
        reference[[differs]], "'). Every input must list the same cells in ",
        "the same order."
      )
    }
  }
  if (anyNA(reference)) {
    .abort("Cell names must not be missing (", names(named)[[1L]], ").")
  }
  unique_names <- .make_unique_names(as.character(reference))
  if (length(unique_names$renamed) > 0L) {
    .inform("Cell names: ", .unique_names_note(0, unique_names$renamed, NULL))
  }
  unique_names$names
}
