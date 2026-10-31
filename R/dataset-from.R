# What an R single-cell object holds, as the components of the dataset plan
# (dataset-plan.R): cellucid_prepare(data = ), cellucid_view() and
# cellucid_serve() read an object through this one extraction. The mapping is
# the web app's for the same object saved with saveRDS() (the transcoders in
# data/anndata/transcode/seurat.js and sce.js), so the obs columns, the obsm
# keys and the modality names agree whichever way the object reaches the
# viewer:
#
# | Seurat                              | SingleCellExperiment / SE          | becomes                       |
# | ----------------------------------- | ---------------------------------- | ----------------------------- |
# | meta.data, active.ident             | colData, sizeFactors               | obs (active_ident, size_factor) |
# | reductions$<r>@cell.embeddings      | reducedDims                        | obsm['X_<name>']              |
# | each assay, the active one first    | the main experiment, then altExps  | one modality each             |
# | assay data (else counts)            | logcounts (else counts, else first) | X                            |
# | meta.features / Assay5 meta.data    | rowData                            | var                           |
# | var.features                        |                                    | var['highly_variable']        |
# | graphs (<active>_snn first)         | colPairs (connectivities first)    | obsp['connectivities']        |
# | misc$cellucid                       | metadata(x)$cellucid               | dataset identity and defaults |
#
# Seurat v5 split layers (counts.1, data.ctrl) are joined first, on a copy;
# the caller's object is never changed. A Cellucid-ready AnnData holds one
# matrix per modality (its layers group is empty, the shared contract), so
# counts beside data, and further assays, are not read. Matrices kept outside
# R's memory (BPCells, HDF5Array and other DelayedArray backends) stay there:
# they are read block by block as they are written or served
# (backed-matrix.R), never whole. scale.data is not read: it is dense, holds
# only the variable genes, and is derived from data.
#
# SeuratObject, SingleCellExperiment and SummarizedExperiment are suggested,
# not imported: each is loaded only when an object of its class arrives.

.dataset_from <- function(data) {
  classes <- class(data)
  package <- attr(classes, "package")
  if (inherits(data, "seurat") && !inherits(data, "Seurat")) {
    .abort(
      "data is a Seurat v2 object (class 'seurat'), which Cellucid does not read. ",
      "Update it first: data <- SeuratObject::UpdateSeuratObject(data)."
    )
  }
  if (!is.null(package) && !identical(package, "methods") && !identical(package, ".GlobalEnv")) {
    .require_object_package(package, classes[[1L]])
  }
  if (methods::is(data, "Seurat")) {
    return(.dataset_from_seurat(data))
  }
  if (methods::is(data, "SummarizedExperiment")) {
    return(.dataset_from_experiment(data))
  }
  .abort(
    "data must be a Seurat, SingleCellExperiment or SummarizedExperiment object; ",
    "got ", .describe_class(data), ". Pass matrices and data frames as the ",
    "other arguments instead (embeddings, obs, gene_expression, ...)."
  )
}

.describe_class <- function(x) {
  package <- attr(class(x), "package")
  paste0("a ", class(x)[[1L]], if (!is.null(package)) paste0(" (", package, ")"), " object")
}

.require_object_package <- function(package, class_name) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install <- if (package %in% c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "DelayedArray", "HDF5Array")) {
      paste0("BiocManager::install(\"", package, "\")")
    } else if (identical(package, "BPCells")) {
      "install.packages(\"BPCells\", repos = c(\"https://bnprks.r-universe.dev\", getOption(\"repos\")))"
    } else {
      paste0("install.packages(\"", package, "\")")
    }
    .abort(
      "data is a ", class_name, " object from the ", package, " package, which is ",
      "not installed, so it cannot be read. Install it with ", install, "."
    )
  }
  invisible(TRUE)
}

# --- Shared rules -------------------------------------------------------------

# A data frame of the columns Cellucid stores (factor, character, logical,
# integer, double; an Rle expanded), and a note naming every other column.
.storable_columns <- function(frame, what) {
  keep <- list()
  skipped <- character()
  for (name in names(frame)) {
    column <- frame[[name]]
    if (methods::is(column, "Rle")) column <- as.vector(column)
    storable <- is.null(dim(column)) && (
      is.factor(column) ||
        (!is.object(column) && (is.character(column) || is.logical(column) ||
          is.integer(column) || is.double(column)))
    )
    if (!storable || !nzchar(name)) {
      skipped <- c(skipped, name)
      next
    }
    keep[[name]] <- column
  }
  # The columns kept as they are. as.data.frame() would pass them to
  # data.frame() as named arguments, and an argument's name is a symbol, which
  # an R session that is not UTF-8 (Windows R before 4.2) makes native: a
  # column named in Greek became "<U+03B1>", and was then refused as a path.
  out <- if (length(keep) == 0L) {
    data.frame(row = seq_len(nrow(frame)))[, 0L, drop = FALSE]
  } else {
    structure(keep, class = "data.frame", row.names = .set_row_names(nrow(frame)))
  }
  note <- if (length(skipped) > 0L) {
    paste0(
      what, ": not stored, no AnnData encoding Cellucid writes: ",
      .format_value_list(skipped), " (convert them to character or numbers to keep them)."
    )
  }
  list(frame = out, note = note)
}

# A matrix as the plan reads it: a base matrix, a Matrix sparse matrix, or one
# held outside R's memory (BPCells, a DelayedArray backend such as
# HDF5Array), which stays where it is and is read block by block when it is
# written or served (backed-matrix.R). Anything else is refused.
.readable_matrix <- function(matrix, where) {
  if (is.matrix(matrix) || inherits(matrix, "Matrix") || !is.null(.backed_package(matrix))) {
    return(matrix)
  }
  .abort(
    where, " is ", .describe_class(matrix), ", which Cellucid does not read. ",
    "Convert it to a dgCMatrix first: as(x, \"dgCMatrix\")."
  )
}

# The obsm embeddings: each matrix's rows placed in the dataset's cell order
# by name; one that covers other cells is left out and named.
.place_embeddings <- function(embeddings, cells, what) {
  out <- list()
  notes <- character()
  for (name in names(embeddings)) {
    matrix <- embeddings[[name]]
    if (is.data.frame(matrix)) matrix <- as.matrix(matrix)
    matrix <- .readable_matrix(matrix, paste0(what, " '", name, "'"))
    # An embedding is a few columns a cell, read whole.
    if (!is.matrix(matrix)) matrix <- as.matrix(matrix)
    if (length(dim(matrix)) != 2L || ncol(matrix) == 0L || !is.numeric(matrix)) {
      notes <- c(notes, paste0(what, " '", name, "' holds no numeric coordinates and was not read."))
      next
    }
    rows <- if (is.null(rownames(matrix))) {
      if (nrow(matrix) == length(cells)) seq_along(cells) else NULL
    } else {
      match(cells, rownames(matrix))
    }
    if (is.null(rows) || anyNA(rows)) {
      notes <- c(notes, paste0(
        what, " '", name, "' covers ", nrow(matrix), " of the ", length(cells), " ",
        .plural(length(cells), "cell"), " and was not read: an embedding has a position for every cell."
      ))
      next
    }
    matrix <- matrix[rows, , drop = FALSE]
    rownames(matrix) <- cells
    key <- paste0(.embedding_key_prefix, name)
    problem <- tryCatch({
      .validate_element_key(key, "Embedding")
      NULL
    }, error = conditionMessage)
    if (!is.null(problem)) {
      notes <- c(notes, problem)
      next
    }
    out[[key]] <- matrix
  }
  list(embeddings = out, notes = notes)
}

# The neighbour graph: the first of `graphs` (in the order given) that is a
# cells x cells matrix over the dataset's cells, exactly symmetric once its
# diagonal is removed. Each graph passed over is named.
.choose_graph <- function(graphs, cells, what) {
  notes <- character()
  for (name in names(graphs)) {
    graph <- graphs[[name]]
    if (!inherits(graph, "Matrix") && !is.matrix(graph)) next
    if (nrow(graph) != length(cells) || ncol(graph) != length(cells)) {
      notes <- c(notes, paste0(what, " '", name, "' is not ", length(cells), " x ", length(cells), " and was not read."))
      next
    }
    if (!is.null(rownames(graph)) && !identical(as.character(rownames(graph)), as.character(cells))) {
      order <- match(cells, rownames(graph))
      if (anyNA(order)) {
        notes <- c(notes, paste0(what, " '", name, "' names other cells and was not read."))
        next
      }
      graph <- graph[order, order, drop = FALSE]
    }
    graph <- methods::as(methods::as(methods::as(graph, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    Matrix::diag(graph) <- 0
    graph <- Matrix::drop0(graph)
    dimnames(graph) <- NULL
    if (!isTRUE(Matrix::isSymmetric(graph, tol = 0))) {
      notes <- c(notes, paste0(
        what, " '", name, "' is not symmetric and was not read: the viewer draws ",
        "an undirected neighbour graph."
      ))
      next
    }
    return(list(graph = graph, name = name, notes = notes))
  }
  list(graph = NULL, name = NULL, notes = notes)
}

# The identity and launch defaults an object declares in its own record
# (Seurat misc$cellucid, SingleCellExperiment metadata(x)$cellucid), shaped
# as uns['cellucid'] is: id, name, description, citation (text, url), source
# (name, url) and defaults (embedding1d/2d/3d, modality, geneIdColumn,
# latentKey, vectorFieldDefault, obsKeys, normalizeEmbeddings,
# centroidMinPoints). The web app reads the same record when the object is
# opened from an .rds file.
.declared_components <- function(record, where) {
  if (is.null(record)) {
    return(list(components = list(), notes = character()))
  }
  if (!is.list(record)) {
    .abort(where, " must be a list such as list(id = \"my-study\", name = \"My study\").")
  }
  scalar <- function(value) if (is.null(value)) NULL else unlist(value, use.names = FALSE)
  notes <- character()
  # citation and source are list(text = , url = ) and list(name = , url = );
  # any other value is named and not read (cellucid-python reads only those).
  member <- function(value, name) if (is.list(value)) value[[name]] else NULL
  for (key in c("citation", "source")) {
    if (!is.null(record[[key]]) && !is.list(record[[key]])) {
      notes <- c(notes, paste0(
        where, "$", key, ": not read, it must be ",
        if (key == "citation") "list(text = , url = )" else "list(name = , url = )", "."
      ))
    }
  }
  out <- list(
    dataset_id = scalar(record$id),
    dataset_name = scalar(record$name),
    dataset_description = scalar(record$description),
    citation_text = scalar(member(record$citation, "text")),
    citation_url = scalar(member(record$citation, "url")),
    source_name = scalar(member(record$source, "name")),
    source_url = scalar(member(record$source, "url"))
  )
  defaults <- record$defaults
  if (!is.null(defaults)) {
    slots <- list()
    for (slot in c("1d", "2d", "3d")) {
      spec <- defaults[[paste0("embedding", slot)]]
      if (is.null(spec)) next
      slots[[slot]] <- if (is.list(spec)) {
        list(key = scalar(spec$key), columns = as.integer(unlist(spec$columns)) + 1L)
      } else {
        scalar(spec)
      }
    }
    if (length(slots) > 0L) out$embedding_defaults <- slots
    out$latent_key <- scalar(defaults$latentKey)
    out$obs_keys <- if (!is.null(defaults$obsKeys)) as.character(unlist(defaults$obsKeys))
    out$var_gene_id_column <- scalar(defaults$geneIdColumn)
    out$vector_field_default <- scalar(defaults$vectorFieldDefault)
    out$modality <- scalar(defaults$modality)
    out$normalize_embeddings <- scalar(defaults$normalizeEmbeddings)
    out$centroid_min_points <- scalar(defaults$centroidMinPoints)
    unknown <- setdiff(names(defaults), c(
      paste0("embedding", c("1d", "2d", "3d")), "latentKey", "obsKeys", "geneIdColumn",
      "vectorFieldDefault", "modality", "normalizeEmbeddings", "centroidMinPoints"
    ))
    if (length(unknown) > 0L) {
      notes <- c(notes, paste0(where, "$defaults: not written, no such launch default: ", .format_value_list(unknown), "."))
    }
  }
  list(components = Filter(Negate(is.null), out), notes = notes)
}

# --- Seurat --------------------------------------------------------------------

.dataset_from_seurat <- function(object) {
  version <- tryCatch(as.character(methods::slot(object, "version")), error = function(error) NA_character_)
  if (is.na(version) || utils::compareVersion(version, "3.0.0") < 0) {
    .abort(
      "This Seurat object was made by Seurat ", if (is.na(version)) "before 3.0" else version,
      ", which Cellucid does not read. Update it first: ",
      "data <- SeuratObject::UpdateSeuratObject(data)."
    )
  }
  assays <- methods::slot(object, "assays")
  v5 <- vapply(assays, function(assay) methods::.hasSlot(assay, "layers"), logical(1L))
  if (any(v5) && utils::packageVersion("SeuratObject") < "5.0.0") {
    .abort(
      "This Seurat object has Seurat v5 assays (", .format_value_list(names(assays)[v5]),
      "), which SeuratObject ", as.character(utils::packageVersion("SeuratObject")),
      " cannot read. Install SeuratObject 5 or later: install.packages(\"SeuratObject\")."
    )
  }
  meta <- methods::slot(object, "meta.data")
  cells <- rownames(meta)
  if (is.null(cells) || length(cells) == 0L) cells <- colnames(object)
  notes <- character()

  stored <- .storable_columns(meta, "Seurat meta.data")
  obs <- stored$frame
  rownames(obs) <- cells
  notes <- c(notes, stored$note)
  ident <- methods::slot(object, "active.ident")
  if (length(ident) == length(cells) && !"active_ident" %in% names(obs)) {
    if (!is.null(names(ident))) ident <- ident[cells]
    obs$active_ident <- unname(ident)
  }

  embeddings <- list()
  for (name in names(methods::slot(object, "reductions"))) {
    embeddings[[name]] <- methods::slot(methods::slot(object, "reductions")[[name]], "cell.embeddings")
  }
  placed <- .place_embeddings(embeddings, cells, "Seurat reduction")
  notes <- c(notes, placed$notes)

  active <- as.character(methods::slot(object, "active.assay"))
  order <- c(intersect(active, names(assays)), setdiff(names(assays), active))
  modalities <- list()
  for (name in order) {
    read <- .seurat_assay(assays[[name]], name, cells)
    notes <- c(notes, read$notes)
    if (!is.null(read$modality)) modalities[[name]] <- read$modality
  }

  graphs <- methods::slot(object, "graphs")
  # The graphs the viewer's default rule would open with, in its order
  # (.graph_candidates: connectivities, <active>_snn, any other *_snn); the
  # first that is a cells x cells symmetric matrix is written.
  candidates <- .graph_candidates(names(graphs), active)
  chosen <- .choose_graph(graphs[candidates], cells, "Seurat graph")
  notes <- c(notes, chosen$notes)
  passed_over <- setdiff(names(graphs), candidates)
  if (length(passed_over) > 0L) {
    notes <- c(notes, paste0(
      "Seurat graph ", paste0("'", passed_over, "'", collapse = ", "),
      " not read: the viewer draws one undirected neighbour graph, a shared-nearest-neighbour one."
    ))
  }

  misc <- tryCatch(methods::slot(object, "misc"), error = function(error) list())
  declared <- .declared_components(misc$cellucid, "misc$cellucid")
  notes <- c(notes, declared$notes)
  components <- declared$components
  project <- tryCatch(as.character(methods::slot(object, "project.name")), error = function(error) "")
  if (is.null(components$dataset_name) && length(project) == 1L && nzchar(project) &&
      !identical(project, "SeuratProject")) {
    components$dataset_name <- project
  }
  components$embeddings <- placed$embeddings
  components$obs <- obs
  components$connectivities <- chosen$graph
  components$modalities <- modalities
  components$notes <- if (length(notes) > 0L) paste0("From the Seurat object: ", notes)
  components$kind <- "Seurat"
  components
}

# One Seurat assay as a modality: its expression (data, else counts), counts
# as a layer beside data, its var and its cells.
.seurat_assay <- function(assay, name, cells) {
  where <- paste0("assays$", name)
  notes <- character()
  if (methods::.hasSlot(assay, "layers")) {
    layer_names <- names(methods::slot(assay, "layers"))
    split <- grepl("^(counts|data)\\.", layer_names)
    if (any(split)) {
      # JoinLayers() returns a new assay; the caller's object is untouched.
      # It exists from SeuratObject 5, which a v5 assay implies (checked in
      # .dataset_from_seurat()); it is looked up when called, so R CMD check
      # with SeuratObject 4 installed does not report it missing.
      assay <- getExportedValue("SeuratObject", "JoinLayers")(assay)
      layer_names <- names(methods::slot(assay, "layers"))
    }
    # The layer as stored, named by the assay's own record of its cells and
    # features (LayerData() would subset it, which copies it).
    layer <- function(kind) {
      if (!kind %in% layer_names) return(NULL)
      matrix <- .readable_matrix(methods::slot(assay, "layers")[[kind]], paste0(where, " layer ", kind))
      dimnames(matrix) <- list(
        SeuratObject::Features(assay, layer = kind),
        SeuratObject::Cells(assay, layer = kind)
      )
      matrix
    }
    data <- layer("data")
    # counts is read only when there is no data: the store holds one matrix
    # per modality.
    counts <- if (is.null(data)) layer("counts")
    other <- setdiff(layer_names, c("data", "counts"))
    if (length(other) > 0L) {
      notes <- c(notes, paste0(where, " layers not read: ", .format_value_list(other), "."))
    }
    features <- rownames(if (!is.null(data)) data else counts)
    meta <- methods::slot(assay, "meta.data")
    var_frame <- if (nrow(meta) == length(SeuratObject::Features(assay))) {
      rownames(meta) <- SeuratObject::Features(assay)
      meta[features, , drop = FALSE]
    } else {
      data.frame(row = seq_along(features))[, 0L, drop = FALSE]
    }
    variable <- tryCatch(SeuratObject::VariableFeatures(assay), error = function(error) {
      notes <<- c(notes, paste0(where, " variable features were not read (", conditionMessage(error), "); ",
                                "var has no highly_variable column."))
      character()
    })
  } else {
    present <- function(matrix) !is.null(matrix) && length(dim(matrix)) == 2L && all(dim(matrix) > 0L)
    data <- methods::slot(assay, "data")
    counts <- if (methods::.hasSlot(assay, "counts")) methods::slot(assay, "counts") else NULL
    if (!present(data)) data <- NULL
    if (!present(counts)) counts <- NULL
    features <- rownames(if (!is.null(data)) data else counts)
    meta <- methods::slot(assay, "meta.features")
    var_frame <- if (nrow(meta) == length(features)) {
      meta[features, , drop = FALSE]
    } else {
      data.frame(row = seq_along(features))[, 0L, drop = FALSE]
    }
    variable <- methods::slot(assay, "var.features")
  }
  primary <- if (!is.null(data)) data else counts
  if (is.null(primary)) {
    return(list(modality = NULL, notes = paste0(where, " has neither data nor counts and was not read.")))
  }
  stored <- .storable_columns(var_frame, paste0(where, " feature metadata"))
  var <- stored$frame
  notes <- c(notes, stored$note)
  rownames(var) <- features
  if (!"highly_variable" %in% names(var) && length(variable) > 0L) {
    var$highly_variable <- features %in% variable
  }
  assay_cells <- colnames(primary)
  list(
    modality = list(
      name = name,
      gene_expression = .feature_major(primary),
      var = var,
      cells = if (identical(as.character(assay_cells), as.character(cells))) NULL else assay_cells
    ),
    notes = notes
  )
}

# --- SingleCellExperiment and SummarizedExperiment -------------------------------

.dataset_from_experiment <- function(object) {
  is_sce <- methods::is(object, "SingleCellExperiment")
  if (is_sce) {
    .require_object_package("SingleCellExperiment", "SingleCellExperiment")
    version <- tryCatch(object@int_metadata$version, error = function(error) NULL)
    if (is.null(version) || !methods::.hasSlot(object, "int_colData") ||
        utils::compareVersion(as.character(version), "1.7.1") < 0) {
      .abort(
        "This SingleCellExperiment was made by SingleCellExperiment ",
        if (is.null(version)) "before 1.8" else as.character(version),
        ", which Cellucid does not read. Update it first: data <- updateObject(data)."
      )
    }
  }
  .require_object_package("SummarizedExperiment", "SummarizedExperiment")
  kind <- if (is_sce) "SingleCellExperiment" else "SummarizedExperiment"
  cells <- colnames(object)
  if (is.null(cells)) cells <- as.character(seq_len(ncol(object)) - 1L)
  notes <- character()
  if (anyDuplicated(cells)) {
    # Made unique as anndata's obs_names_make_unique() does (cellucid-python's
    # resolve_obs): every part of the object follows its columns, so they are
    # renamed in place, and the assays and embeddings keep their cells.
    unique_names <- .make_unique_names(as.character(cells))
    .inform("Cell names: ", .unique_names_note(0, unique_names$renamed, NULL))
    cells <- unique_names$names
    colnames(object) <- cells
  }

  col_data <- SummarizedExperiment::colData(object)
  frame <- if (ncol(col_data) == 0L) {
    data.frame(row = seq_len(ncol(object)))[, 0L, drop = FALSE]
  } else {
    structure(
      lapply(stats::setNames(nm = names(col_data)), function(name) col_data[[name]]),
      class = "data.frame", row.names = .set_row_names(ncol(object))
    )
  }
  stored <- .storable_columns(frame, paste0(kind, " colData"))
  obs <- stored$frame
  notes <- c(notes, stored$note)
  rownames(obs) <- cells

  embeddings <- list()
  if (is_sce) {
    factors <- SingleCellExperiment::sizeFactors(object)
    if (!is.null(factors) && !any(c("size_factor", "sizeFactor") %in% names(obs))) {
      obs$size_factor <- as.double(factors)
    }
    embeddings <- as.list(SingleCellExperiment::reducedDims(object))
  }
  placed <- .place_embeddings(embeddings, cells, "reducedDim")
  notes <- c(notes, placed$notes)

  main <- if (is_sce) {
    name <- tryCatch(SingleCellExperiment::mainExpName(object), error = function(error) NULL)
    if (is.null(name) || is.na(name) || !nzchar(name)) "main" else name
  } else {
    "main"
  }
  modalities <- list()
  read <- .experiment_modality(object, main, cells, kind)
  notes <- c(notes, read$notes)
  if (!is.null(read$modality)) modalities[[main]] <- read$modality
  if (is_sce) {
    for (name in SingleCellExperiment::altExpNames(object)) {
      read <- .experiment_modality(SingleCellExperiment::altExp(object, name), name, cells, paste0("altExp '", name, "'"))
      notes <- c(notes, read$notes)
      if (!is.null(read$modality)) modalities[[name]] <- read$modality
    }
  }

  graph <- list(graph = NULL, notes = character())
  if (is_sce) {
    pair_names <- SingleCellExperiment::colPairNames(object)
    # The colPairs the viewer's default rule would open with, in its order (.graph_candidates).
    preferred <- .graph_candidates(pair_names)
    pairs <- lapply(stats::setNames(nm = preferred), function(name) {
      SingleCellExperiment::colPair(object, name, asSparse = TRUE)
    })
    graph <- .choose_graph(pairs, cells, "colPair")
    notes <- c(notes, graph$notes)
  }

  declared <- .declared_components(S4Vectors::metadata(object)$cellucid, "metadata(x)$cellucid")
  notes <- c(notes, declared$notes)
  components <- declared$components
  components$embeddings <- placed$embeddings
  components$obs <- obs
  components$connectivities <- graph$graph
  components$modalities <- modalities
  components$notes <- if (length(notes) > 0L) paste0("From the ", kind, ": ", notes)
  components$kind <- kind
  components
}

# One experiment (the main one or an alternative one) as a modality: logcounts
# (else counts, else the first assay) as X; the store holds one matrix per
# modality, so the other assays are not read.
.experiment_modality <- function(experiment, name, cells, what) {
  assay_names <- SummarizedExperiment::assayNames(experiment)
  n_assays <- length(SummarizedExperiment::assays(experiment))
  if (n_assays == 0L) {
    return(list(modality = NULL, notes = paste0(what, " has no assay and was not read.")))
  }
  if (is.null(assay_names) || any(!nzchar(assay_names))) {
    return(list(modality = NULL, notes = paste0(what, " has unnamed assays and was not read; name them with assayNames().")))
  }
  primary <- if ("logcounts" %in% assay_names) {
    "logcounts"
  } else if ("counts" %in% assay_names) {
    "counts"
  } else {
    assay_names[[1L]]
  }
  read <- function(assay_name) {
    matrix <- SummarizedExperiment::assay(experiment, assay_name, withDimnames = TRUE)
    .readable_matrix(matrix, paste0(what, " assay '", assay_name, "'"))
  }
  row_data <- SummarizedExperiment::rowData(experiment)
  frame <- if (ncol(row_data) == 0L) {
    data.frame(row = seq_len(nrow(experiment)))[, 0L, drop = FALSE]
  } else {
    structure(
      lapply(stats::setNames(nm = names(row_data)), function(column) row_data[[column]]),
      class = "data.frame", row.names = .set_row_names(nrow(experiment))
    )
  }
  stored <- .storable_columns(frame, paste0(what, " rowData"))
  var <- stored$frame
  features <- rownames(experiment)
  if (!is.null(features)) rownames(var) <- features
  experiment_cells <- colnames(experiment)
  list(
    modality = list(
      name = name,
      gene_expression = .feature_major(read(primary)),
      var = var,
      cells = if (is.null(experiment_cells) || identical(as.character(experiment_cells), as.character(cells))) {
        NULL
      } else {
        experiment_cells
      }
    ),
    notes = stored$note
  )
}
