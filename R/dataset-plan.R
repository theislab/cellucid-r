# The dataset plan: every component checked, and the function that puts the
# whole store into a store backend. cellucid_prepare() writes it to a folder
# or a ZIP file; cellucid_view() and cellucid_serve() hold it in memory and
# serve it (zarr-store.R). One plan for both is what makes the served store
# the same bytes as the prepared one.
#
# A dataset of one feature space is a Cellucid-ready AnnData. A dataset of
# several (Seurat assays, a SingleCellExperiment's alternative experiments) is
# a Cellucid-ready MuData: the cells' obs, obsm, obsp and uns['cellucid'] at
# the root, and each modality an AnnData under mod/<name> (mudata.R).

# The components the plan is made from: what `data` holds (dataset-from.R),
# with every argument the caller gave in its place. `given` is a named list
# of the caller's arguments, NULL where not given. The arguments that
# describe genes -- gene_expression, var and gene_identifiers --
# replace an object's only feature space, and are refused for an object of
# several (each has its own genes); var_gene_id_column and modality then name
# a column of, and the modality, the viewer opens with.
.dataset_components <- function(data, given) {
  if (is.null(data)) {
    if (!is.null(given$modalities)) {
      .abort("modalities chooses among the modalities of data, and no data was supplied.")
    }
    components <- given[setdiff(names(given), c("gene_expression", "var", "modalities"))]
    components$modalities <- list(list(
      name = NULL,
      gene_expression = given$gene_expression,
      var = given$var,
      cells = NULL
    ))
    components$notes <- character()
    return(components)
  }
  components <- .dataset_from(data)
  for (name in setdiff(names(given), c("gene_expression", "var", "modalities"))) {
    if (!is.null(given[[name]])) components[[name]] <- given[[name]]
  }
  if (length(components$embeddings) == 0L) {
    compute <- if (identical(components$kind, "Seurat")) {
      "data <- Seurat::RunUMAP(data, dims = 1:30)"
    } else {
      "data <- scater::runUMAP(data)"
    }
    .abort(
      "This ", components$kind, " holds no embedding (no ",
      if (identical(components$kind, "Seurat")) "reductions" else "reducedDims",
      ") for Cellucid to place its cells at. Compute one first, for example ",
      compute, ", or pass embeddings = list(X_umap = coordinates)."
    )
  }
  modalities <- components$modalities
  if (!is.null(given$modalities)) {
    .validate_character_vector(given$modalities, "modalities")
    if (length(given$modalities) == 0L || anyDuplicated(given$modalities)) {
      .abort("modalities must name at least one modality of data, each once.")
    }
    unknown <- setdiff(given$modalities, names(modalities))
    if (length(unknown) > 0L) {
      .abort(
        "modalities names ", .format_value_list(unknown), ", which data does not hold. ",
        "Its modalities: ", .format_value_list(names(modalities)), "."
      )
    }
    modalities <- modalities[given$modalities]
  }
  if (length(modalities) == 0L) {
    .abort("data holds no expression matrix Cellucid reads.")
  }
  gene_arguments <- c("gene_expression", "var", "gene_identifiers")
  given_genes <- gene_arguments[!vapply(given[gene_arguments], is.null, logical(1L))]
  if (length(modalities) > 1L && length(given_genes) > 0L) {
    # Each modality has its own genes: change the object's assay first, or
    # write one modality (cellucid-python refuses them with a MuData too).
    .abort(
      .format_value_list(given_genes), " describe the genes of one feature space, and data ",
      "holds several (", .format_value_list(names(modalities)), "). Change the assay in ",
      "the object first, or write one modality with modalities = \"", names(modalities)[[1L]], "\"."
    )
  }
  if (length(modalities) == 1L) {
    for (name in c("gene_expression", "var")) {
      if (!is.null(given[[name]])) modalities[[1L]][[name]] <- given[[name]]
    }
    if (!is.null(given$gene_expression)) {
      # A new matrix has other genes than the object's var describes.
      if (is.null(given$var)) modalities[[1L]]["var"] <- list(NULL)
      modalities[[1L]]$cells <- NULL
    }
  }
  if (length(modalities) == 1L) {
    # One feature space is an AnnData, whatever the object called it; the
    # name is kept to say which one a message is about.
    modalities[[1L]]$source <- modalities[[1L]]$name
    modalities[[1L]]$name <- NULL
  }
  components$modalities <- unname(modalities)
  components
}

# The modality the viewer opens with: the one named, else the shared rule.
.launch_modality <- function(names, modality) {
  if (length(names) == 1L) {
    return(1L)
  }
  chosen <- if (is.null(modality)) .default_modality(names) else modality
  if (!is.character(chosen) || length(chosen) != 1L || !chosen %in% names) {
    .abort(
      "modality must name one of the modalities: ",
      .format_value_list(names), "."
    )
  }
  match(chosen, names)
}

# The names a gene-expression modality goes by, lower-cased: MuData's rna,
# Cell Ranger's Gene Expression, Seurat's RNA and SCT, GEX.
.gene_expression_modalities <- c("rna", "gex", "gene expression", "gene_expression", "sct")

# The modality a dataset opens with, of its modalities in stored order: the
# first whose name is a gene-expression name (any letter case), else the
# first. The web app (defaultModality in data/data-mapping.js) and
# cellucid-python (default_modality) apply the same rule, which the shared
# data-mapping conformance fixture pins (defaultModality). A Seurat object's
# active assay is listed first, a SingleCellExperiment's main experiment.
.default_modality <- function(names) {
  if (length(names) == 0L) {
    return(NULL)
  }
  gene <- which(tolower(names) %in% .gene_expression_modalities)
  names[[if (length(gene) > 0L) gene[[1L]] else 1L]]
}

# --- The plan ----------------------------------------------------------------

.plan_dataset <- function(components) {
  if (length(components$notes) > 0L) {
    .inform(paste(components$notes, collapse = "\n"))
  }
  .require_identity(components)
  dataset_id <- .validate_dataset_id(components$dataset_id)
  dataset_name <- .validate_dataset_name(components$dataset_name)
  dataset_description <- .validate_optional_text(
    components$dataset_description, "dataset_description", allow_empty = TRUE
  )
  citation <- .validate_linked_text(
    components$citation_text, components$citation_url, "citation_text", "citation_url", "text"
  )
  source <- .validate_linked_text(
    components$source_name, components$source_url, "source_name", "source_url", "name"
  )

  # --- obsm: the embeddings and the vector fields -------------------------
  resolved <- .resolve_embeddings(components$embeddings)
  matrices <- resolved$matrices
  n_cells <- resolved$n_cells
  vectors <- .resolve_vector_fields(components$vector_fields, matrices, n_cells)

  # --- obs ----------------------------------------------------------------
  obs <- components$obs
  if (is.null(obs)) {
    obs <- data.frame(row = seq_len(n_cells))[, 0L, drop = FALSE]
  } else if (!is.data.frame(obs)) {
    .abort("obs must be a data.frame with one row per cell.")
  }
  if (nrow(obs) != n_cells) {
    .abort(
      "obs has ", nrow(obs), " ", .plural(nrow(obs), "row", "rows"),
      ", but the embeddings have ", n_cells, " ", .plural(n_cells, "cell", "cells"), "."
    )
  }
  .validate_frame_columns(obs, "obs column")

  # --- The feature spaces -------------------------------------------------
  modalities <- components$modalities
  multimodal <- length(modalities) > 1L
  names <- vapply(modalities, function(modality) {
    if (is.null(modality$name)) "" else modality$name
  }, character(1L))
  if (multimodal) {
    for (name in names) .validate_element_key(name, "Modality")
    .require_distinct_element_keys(names, "Modality")
  }
  launch <- if (multimodal) .launch_modality(names, components$modality) else 1L
  if (!multimodal && !is.null(components$modality)) {
    .abort(
      "modality names a modality, and the data has one feature space, ",
      "which is written as an AnnData with no modalities."
    )
  }
  spaces <- vector("list", length(modalities))
  for (index in seq_along(modalities)) {
    is_launch <- index == launch
    spaces[[index]] <- .resolve_feature_space(
      modalities[[index]],
      label = if (multimodal) paste0("modality '", names[[index]], "'") else NULL,
      gene_identifiers = if (is_launch) components$gene_identifiers else NULL,
      var_gene_id_column = if (is_launch) components$var_gene_id_column else NULL
    )
  }
  if (!is.null(spaces[[launch]]$unused) && length(modalities) > 0L) {
    supplied <- spaces[[launch]]$unused
    .abort(
      supplied, " describes the genes of gene_expression, and no gene_expression was supplied."
    )
  }

  if (multimodal) {
    taken <- intersect(tolower(names), tolower(c(names(matrices), names(vectors$arrays))))
    if (length(taken) > 0L) {
      .abort(
        "The obsm key '", taken[[1L]], "' is also a modality's name, under which MuData ",
        "keeps the modality's mask of measured cells. Rename the embedding."
      )
    }
  }

  # --- Cell names -----------------------------------------------------------
  whole <- Filter(function(space) is.null(space$cells), spaces)
  cell_sources <- c(
    list(obs = .supplied_row_names(obs)),
    if (length(whole) > 0L && !is.null(whole[[1L]]$expression)) {
      list(gene_expression = .expression_cell_names(whole[[1L]]$expression))
    },
    `names<-`(
      lapply(matrices, rownames),
      paste0("embedding '", names(matrices), "'")
    )
  )
  obs_index <- .resolve_cell_names(cell_sources, n_cells)
  for (index in seq_along(spaces)) {
    spaces[[index]] <- .place_feature_space(spaces[[index]], obs_index, n_cells, names[[index]])
  }
  if (!multimodal) {
    # An AnnData's X has one row per cell, in the cells' order: MuData's
    # obsmap places a modality's rows, an AnnData has none.
    spaces[[1L]] <- .align_single_space(spaces[[1L]], n_cells)
  }

  # --- obsp -----------------------------------------------------------------
  graph <- NULL
  if (!is.null(components$connectivities)) {
    graph <- .resolve_connectivities(components$connectivities, n_cells)
  }

  # --- Priority order (R/priority.R) ---------------------------------------
  # Every per-cell element in priority order before anything is derived from
  # it (the outlier quantiles), so the file holds what its rows give.
  priority <- .cell_priority(n_cells, components$priority_seed)
  # The cells by original position, which the viewer's positions name: the
  # cells' order as given (an object in a verified stored order was put back
  # in its original order first, .in_original_order()).
  cell_names <- obs_index
  if (!is.null(priority)) {
    rows <- priority$order
    matrices <- lapply(matrices, function(values) values[rows, , drop = FALSE])
    vectors$arrays <- lapply(vectors$arrays, function(values) values[rows, , drop = FALSE])
    obs <- .obs_in_priority(obs, priority)
    obs_index <- obs_index[rows]
    if (!is.null(graph)) graph <- .graph_in_priority(graph, priority)
    for (index in seq_along(spaces)) spaces[[index]] <- .space_in_priority(spaces[[index]], priority)
  }

  # --- uns['cellucid'] ------------------------------------------------------
  defaults <- .resolve_defaults(
    matrices = matrices,
    obs = obs,
    vectors = vectors,
    embedding_defaults = components$embedding_defaults,
    latent_key = components$latent_key,
    obs_keys = components$obs_keys,
    var_gene_id_column = components$var_gene_id_column,
    vector_field_default = components$vector_field_default,
    normalize_embeddings = components$normalize_embeddings,
    centroid_min_points = components$centroid_min_points,
    modality = if (multimodal && !is.null(components$modality)) names[[launch]]
  )
  cellucid <- list(schema = .uns_integer(1L), id = .uns_string(dataset_id), name = .uns_string(dataset_name))
  if (!is.null(dataset_description)) cellucid$description <- .uns_string(dataset_description)
  if (!is.null(citation)) cellucid$citation <- citation
  if (!is.null(source)) cellucid$source <- source
  if (length(defaults) > 0L) cellucid$defaults <- defaults
  if (!is.null(defaults[["latentKey"]])) {
    latent <- defaults[["latentKey"]]$value
    min_points <- defaults[["centroidMinPoints"]]$value
    if (is.null(min_points)) min_points <- .default_outlier_min_points
    cellucid$outlierQuantiles <- .stored_outlier_quantiles(obs, latent, matrices[[latent]], min_points)
  }
  if (!is.null(priority)) cellucid$priority <- .priority_record(priority)
  uns <- list(title = .uns_string(dataset_name), cellucid = cellucid)

  shared <- list(
    obs = obs, obs_index = obs_index, matrices = matrices, vectors = vectors,
    graph = graph, uns = uns
  )
  write <- if (multimodal) {
    function(store) .write_mudata(store, shared, spaces, names)
  } else {
    function(store) .write_anndata(store, shared, spaces[[1L]])
  }
  launch_space <- spaces[[launch]]
  list(
    write = write,
    kind = if (multimodal) "mudata-zarr" else "zarr",
    dataset_id = dataset_id,
    dataset_name = dataset_name,
    n_cells = n_cells,
    obs_index = obs_index,
    cell_names = cell_names,
    n_genes = if (is.null(launch_space$expression)) 0L else length(launch_space$index),
    modalities = if (multimodal) names else NULL,
    launch_modality = if (multimodal) names[[launch]] else NULL
  )
}

# dataset_id and dataset_name, reported together when either is missing.
.require_identity <- function(components) {
  absent <- c(
    if (is.null(components$dataset_id)) "dataset_id",
    if (is.null(components$dataset_name)) "dataset_name"
  )
  if (length(absent) == 0L) {
    return(invisible(NULL))
  }
  .abort(
    paste(absent, collapse = " and "), if (length(absent) > 1L) " are" else " is",
    " required: ", .required_argument_guidance[[absent[[1L]]]], ". ",
    "An object can declare its own in misc$cellucid (Seurat) or ",
    "metadata(x)$cellucid (SingleCellExperiment): list(id = , name = )."
  )
}

# One feature space: its expression matrix (or none), the genes that name its
# columns and its var, checked against one another.
.resolve_feature_space <- function(modality, label, gene_identifiers, var_gene_id_column) {
  what <- function(name) if (is.null(label)) name else paste0(label, " ", name)
  if (is.null(modality$gene_expression)) {
    supplied <- c(
      var = !is.null(modality$var),
      gene_identifiers = !is.null(gene_identifiers),
      var_gene_id_column = !is.null(var_gene_id_column)
    )
    return(list(
      expression = NULL,
      var = data.frame(row = integer(0))[, 0L, drop = FALSE],
      index = character(0),
      cells = NULL,
      unused = if (any(supplied)) names(supplied)[supplied][[1L]]
    ))
  }
  expression <- .expression_input(modality$gene_expression, what("gene_expression"))
  genes <- .resolve_genes(expression, modality$var, var_gene_id_column, gene_identifiers)
  expression <- genes$expression
  # The dtype is decided on the genes stored, after gene_identifiers chose
  # them, as the Python writer decides it.
  list(
    expression = expression,
    dtype = .expression_value_dtype(expression),
    var = genes$var,
    index = genes$index,
    cells = modality$cells,
    label = label,
    source = modality$source
  )
}

# Where a modality's cells are among all cells: NULL when it holds every cell
# in the same order, else its rows' 1-based positions (MuData's obsmap).
.place_feature_space <- function(space, obs_index, n_cells, name) {
  if (is.null(space$expression)) {
    return(space)
  }
  n_rows <- .expression_dims(space$expression)[[1L]]
  if (is.null(space$cells)) {
    if (n_rows != n_cells) {
      .abort(
        if (is.null(space$label)) "gene_expression" else space$label, " has ", n_rows, " ",
        .plural(n_rows, "row", "rows"), ", but the embeddings have ", n_cells, " ",
        .plural(n_cells, "cell", "cells"), ". gene_expression is cells x genes; ",
        "transpose a genes x cells matrix with t()."
      )
    }
    space$obs_index <- obs_index
    space$positions <- NULL
    return(space)
  }
  positions <- match(space$cells, obs_index)
  if (anyNA(positions) || anyDuplicated(positions)) {
    .abort(
      "Modality '", name, "' names cells that are not the dataset's cells, ",
      "or one cell twice."
    )
  }
  if (length(positions) == n_cells && identical(positions, seq_len(n_cells))) {
    space$positions <- NULL
  } else {
    space$positions <- positions
  }
  space$obs_index <- obs_index[positions]
  space
}

# The one feature space of an AnnData, its rows in the dataset's cell order:
# an assay that holds every cell in another order is reordered (a copy of
# its matrix); one that holds only some of the cells is refused, since an
# AnnData's X has a row for every cell.
.align_single_space <- function(space, n_cells) {
  positions <- space$positions
  if (is.null(space$expression) || is.null(positions)) {
    return(space)
  }
  what <- if (!is.null(space$source)) paste0("Modality '", space$source, "'") else "gene_expression"
  if (length(positions) != n_cells) {
    .abort(
      what, " holds ", length(positions), " of the dataset's ", n_cells, " ", .plural(n_cells, "cell"), ", and one ",
      "feature space is written as an AnnData, whose X has a row for every cell. Keep only ",
      "those cells first (for a Seurat object: subset(data, cells = Cells(data[[\"",
      if (is.null(space$source)) "<assay>" else space$source, "\"]]))), or write it ",
      "beside an assay of every cell as a MuData (modalities = c(\"<assay>\", \"",
      if (is.null(space$source)) "<assay>" else space$source, "\"))."
    )
  }
  space$expression <- .expression_select_cells(space$expression, order(positions))
  space$positions <- NULL
  space$obs_index <- space$obs_index[order(positions)]
  space
}

# --- Writing -------------------------------------------------------------------

# The elements of one AnnData group at `prefix` ("" for the root): X, obs,
# var, obsm, obsp, varm, varp, layers (empty) and uns, in the order anndata writes
# them. `obs`, `obsm` and the graph are the shared cells' (the root of an
# AnnData), or none (a MuData modality, whose cells' annotations are the
# root's).
.write_anndata_group <- function(store, prefix, space, obs, obs_index, matrices, vectors, graph, uns) {
  key <- function(...) .zarr_key(prefix, ...)
  if (!is.null(space$expression)) {
    .expression_put(store, key("X"), space$expression, space$dtype)
  }
  .ad_dataframe(store, key("obs"), obs, obs_index, "obs")
  .ad_dataframe(store, key("var"), space$var, space$index, "var")
  .ad_dict(store, key("obsm"))
  for (name in names(matrices)) {
    .ad_dense(store, key("obsm", name), matrices[[name]])
  }
  for (name in names(vectors$arrays)) {
    .ad_dense(store, key("obsm", name), vectors$arrays[[name]])
  }
  .ad_dict(store, key("obsp"))
  if (!is.null(graph)) .ad_graph(store, key("obsp/connectivities"), graph)
  .ad_dict(store, key("varm"))
  # A sparse X in priority order: each gene's values below each prefix
  # level, one column per level (R/priority.R).
  if (!is.null(space$prefix_levels) && !is.null(space$expression)) {
    counts <- .expression_prefix_counts(space$expression, space$prefix_levels)
    if (!is.null(counts)) {
      .zarr_matrix(
        store, key("varm", .prefix_counts_key), counts, "<u4",
        .chunk_length(nrow(counts), .chunk_elements_whole), 1L, .encoding("array", "0.2.0")
      )
    }
  }
  .ad_dict(store, key("varp"))
  .ad_dict(store, key("layers"))
  .ad_uns_value(store, key("uns"), uns)
  invisible(NULL)
}

# The neighbour graph, obsp['connectivities'], as a csr_matrix: a symmetric
# dgCMatrix's column-compressed arrays are its row-compressed ones, written as
# they are, or a .graph_in_priority() one, both axes in priority order: its
# pointers counted once, when first needed (one integer a cell), and its
# indices and values made one chunk at a time from the graph as it is --
# row r is column order[r], its columns renumbered and sorted (src/chunks.c)
# -- each chunk made when it is written or read and then let go, as
# cellucid-python serves a graph. Nothing as large as the graph is held
# beside it.
.ad_graph <- function(store, path, graph) {
  if (!inherits(graph, "cellucid_graph_in_order")) {
    return(.ad_sparse(store, path, "csr_matrix", dim(graph), pointers = graph@p, indices = graph@i, values = graph@x))
  }
  matrix <- graph$matrix
  shape <- dim(matrix)
  order <- graph$order
  pointers <- .deferred(function() c(0L, cumsum(diff(matrix@p)[order])), length = shape[[1L]] + 1)
  index_size <- .zarr_item_size(.sparse_index_dtype("csr_matrix", shape))
  new_row <- graph$new_row - 1L
  column_order <- order - 1L
  part <- function(which) {
    structure(list(
      length = length(matrix@x),
      chunk = function(first, count, chunk_length, ...) {
        .native_renumbered_block(matrix, new_row, first - 1, first - 1 + count, index_size, "<f4",
                                 column_order = column_order, pointers = .realize(pointers),
                                 part = which, padded = chunk_length)
      }
    ), class = "cellucid_chunked")
  }
  .ad_sparse(store, path, "csr_matrix", shape, pointers = pointers, indices = part(1L), values = part(2L))
}

.write_anndata <- function(store, shared, space) {
  if (!is.null(space$expression) && .expression_dims(space$expression)[[1L]] != length(shared$obs_index)) {
    n_rows <- .expression_dims(space$expression)[[1L]]
    .abort("Internal error: X has ", n_rows, " ", .plural(n_rows, "row"), " for ",
           length(shared$obs_index), " ", .plural(length(shared$obs_index), "cell"), ".")
  }
  .zarr_group(store, "", .encoding("anndata", "0.1.0"))
  .write_anndata_group(
    store, "", space,
    obs = shared$obs, obs_index = shared$obs_index,
    matrices = shared$matrices, vectors = shared$vectors,
    graph = shared$graph, uns = shared$uns
  )
}

# The genes: the expression matrix narrowed to gene_identifiers, and the var
# data frame and index that describe its columns. `expression` is an
# .expression_input().
.resolve_genes <- function(expression, var, var_gene_id_column, gene_identifiers) {
  n_genes <- .expression_dims(expression)[[2L]]
  if (is.null(var)) {
    var <- data.frame(row = seq_len(n_genes))[, 0L, drop = FALSE]
  } else if (!is.data.frame(var)) {
    .abort("var must be a data.frame with one row per gene.")
  }
  if (nrow(var) != n_genes) {
    .abort(
      "var has ", nrow(var), " ", .plural(nrow(var), "row", "rows"),
      ", but gene_expression has ", n_genes, " ", .plural(n_genes, "gene", "genes"),
      " (columns)."
    )
  }
  .validate_frame_columns(var, "var column")
  if (!is.null(var_gene_id_column)) {
    var_gene_id_column <- .validate_gene_name_source(var_gene_id_column, "var_gene_id_column")
  }

  row_names <- .supplied_row_names(var)
  column_names <- .expression_gene_names(expression)
  if (!is.null(row_names) && !is.null(column_names) && !identical(row_names, column_names)) {
    differs <- which(row_names != column_names)[[1L]]
    .abort(
      "rownames(var) and colnames(gene_expression) disagree at gene ", differs,
      " ('", row_names[[differs]], "' and '", column_names[[differs]], "'). ",
      "var must describe the columns of gene_expression in their order."
    )
  }
  index <- if (!is.null(row_names)) row_names else column_names
  if (is.null(index)) {
    if (is.null(var_gene_id_column)) {
      .abort(
        "The genes have no names: set rownames(var) or ",
        "colnames(gene_expression) to the gene identifiers, or pass ",
        "var_gene_id_column."
      )
    }
    index <- .resolve_gene_name_set(var, var_gene_id_column, "var_gene_id_column")$names
  }
  if (anyNA(index) || any(!nzchar(index))) {
    .abort("Gene names must not be missing or empty.")
  }
  unique_index <- .make_unique_names(as.character(index))
  if (length(unique_index$renamed) > 0L) {
    .inform("Gene names: ", .unique_names_note(0, unique_index$renamed, NULL))
  }
  index <- unique_index$names
  rownames(var) <- index

  if (!is.null(gene_identifiers)) {
    default_names <- if (is.null(var_gene_id_column)) {
      index
    } else {
      .resolve_gene_name_set(var, var_gene_id_column, "var_gene_id_column")$names
    }
    .validate_character_vector(gene_identifiers, "gene_identifiers")
    if (anyDuplicated(gene_identifiers)) {
      .abort("gene_identifiers must not contain duplicate identifiers.")
    }
    positions <- match(gene_identifiers, default_names)
    missing_genes <- gene_identifiers[is.na(positions)]
    if (length(missing_genes) > 0L) {
      shown <- utils::head(missing_genes, 5L)
      remainder <- length(missing_genes) - length(shown)
      .abort(
        "gene_identifiers contains identifiers not found in var: ",
        .format_value_list(shown),
        if (remainder > 0L) paste0(" and ", remainder, " more") else "",
        ". Identifiers are matched against the gene names the viewer shows ",
        "(var_gene_id_column, or the var row names), after repeated names ",
        "were made unique."
      )
    }
    expression <- .expression_select_genes(expression, positions)
    var <- var[positions, , drop = FALSE]
    index <- index[positions]
  }
  if (!is.null(var_gene_id_column)) {
    # The column is stored as it is; the viewer fills its empty values from
    # the index and makes repeated names unique when it reads it, and this
    # says ahead of time what that will change.
    shown <- .resolve_gene_name_set(var, var_gene_id_column, "var_gene_id_column")
    note <- .unique_names_note(shown$filled_from_index, shown$renamed, var_gene_id_column)
    if (!is.null(note)) {
      .inform("Gene names from var column '", var_gene_id_column, "': ", note)
    }
  }
  rownames(var) <- NULL
  list(expression = expression, var = var, index = index)
}

# The launch defaults of uns['cellucid']['defaults']: only what the caller
# set, so every other choice stays the viewer's own. `modality` is the
# modality a multimodal dataset opens with, when modality named one.
.resolve_defaults <- function(
    matrices,
    obs,
    vectors,
    embedding_defaults,
    latent_key,
    obs_keys,
    var_gene_id_column,
    vector_field_default,
    normalize_embeddings,
    centroid_min_points,
    modality = NULL
) {
  defaults <- .resolve_embedding_defaults(embedding_defaults, matrices)
  if (!is.null(modality)) {
    defaults$modality <- .uns_string(modality)
  }
  if (!is.null(var_gene_id_column)) {
    defaults$geneIdColumn <- .uns_string(var_gene_id_column)
  }
  if (!is.null(latent_key)) {
    latent_key <- .validate_string(latent_key, "latent_key")
    if (!latent_key %in% names(matrices)) {
      .abort(
        "latent_key must name one of the embeddings: ",
        .format_value_list(names(matrices)), "."
      )
    }
    defaults$latentKey <- .uns_string(latent_key)
  }
  if (!is.null(vector_field_default)) {
    vector_field_default <- .validate_string(vector_field_default, "vector_field_default")
    if (length(vectors$fields) == 0L) {
      .abort("vector_field_default names a field, and no vector_fields were supplied.")
    }
    if (!vector_field_default %in% vectors$fields) {
      .abort(
        "vector_field_default must name one of vector_fields: ",
        .format_value_list(vectors$fields), "."
      )
    }
    # The overlay opens on a field only where the viewer shows one of its
    # embeddings as the file opens (the web app's defaultVectorBindings).
    named_slots <- lapply(
      defaults[startsWith(names(defaults), "embedding")],
      function(ref) if (is.null(ref$key)) ref$value else ref$key$value
    )
    shown <- .shown_embedding_keys(matrices, named_slots, latent_key)
    drawn <- unique(vapply(
      Filter(function(binding) binding$embedding %in% shown, vectors$bound),
      function(binding) binding$field,
      character(1L)
    ))
    if (!vector_field_default %in% drawn) {
      .abort(
        "vector_field_default names '", vector_field_default, "', and the viewer ",
        "shows none of its embeddings when the file opens, so it is not drawn. ",
        "Show one of them (embedding_defaults), or choose a field drawn on a ",
        "shown embedding: ", .format_value_list(sort(drawn)), "."
      )
    }
    defaults$vectorFieldDefault <- .uns_string(vector_field_default)
  }
  if (!is.null(obs_keys)) {
    .validate_character_vector(obs_keys, "obs_keys")
    if (length(obs_keys) == 0L) {
      .abort("obs_keys must name at least one obs column.")
    }
    if (anyDuplicated(obs_keys)) {
      .abort("obs_keys names '", obs_keys[duplicated(obs_keys)][[1L]], "' twice.")
    }
    absent <- setdiff(obs_keys, names(obs))
    if (length(absent) > 0L) {
      .abort(
        "obs_keys contain columns not in obs: ", .format_value_list(absent),
        ". Available columns: ", .format_value_list(names(obs)), "."
      )
    }
    defaults$obsKeys <- .uns_strings(obs_keys)
  }
  if (!is.null(normalize_embeddings)) {
    defaults$normalizeEmbeddings <- .uns_boolean(
      .validate_flag(normalize_embeddings, "normalize_embeddings")
    )
  }
  if (!is.null(centroid_min_points)) {
    defaults$centroidMinPoints <- .uns_integer(
      .validate_positive_integer(centroid_min_points, "centroid_min_points")
    )
  }
  order <- c(
    "embedding1d", "embedding2d", "embedding3d", "modality", "geneIdColumn", "latentKey",
    "vectorFieldDefault", "obsKeys", "normalizeEmbeddings", "centroidMinPoints"
  )
  defaults[intersect(order, names(defaults))]
}
