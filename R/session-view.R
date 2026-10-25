# An R object as the viewer showed it, for applying a session to it: the
# same extraction cellucid_view() serves (dataset-from.R, dataset-plan.R),
# so a session is checked against exactly what it was saved on -- the cells
# and their obs columns, the obsm embeddings (X_<reduction>), the neighbour
# graph, and each modality's genes and matrix. cellucid-python's
# anndata_session (AnnData, and MuDataSessionView for several modalities)
# checks a session against the same picture.

.session_view <- function(data, extra) {
  components <- .dataset_components(data, extra)
  obs <- components$obs
  cells <- rownames(obs)
  obsm <- components$embeddings
  for (field in names(components$vector_fields)) {
    for (basis in names(components$vector_fields[[field]])) {
      obsm[[paste0(field, "_", sub("^X_", "", basis))]] <- components$vector_fields[[field]][[basis]]
    }
  }
  modalities <- components$modalities
  names(modalities) <- vapply(modalities, function(item) item$name %||% "", character(1L))
  list(
    n_cells = length(cells),
    cells = cells,
    obs = obs,
    obsm = obsm,
    obsp = if (is.null(components$connectivities)) character() else "connectivities",
    modalities = modalities,
    modality_names = if (length(modalities) > 1L) names(modalities) else NULL,
    gene_identifiers = components$gene_identifiers
  )
}

# The modality a mapping reads: the named one, or the one feature space.
.session_modality <- function(view, mapping) {
  if (is.null(view$modality_names)) view$modalities[[1L]] else view$modalities[[mapping[["modality"]]]]
}

.modality_matrix <- function(modality) {
  matrix <- modality$gene_expression
  if (inherits(matrix, "cellucid_feature_major")) matrix$matrix else matrix
}

# The genes a mapping shows, named by its gene-name set and made unique by
# the shared rule, as the viewer named them; narrowed to gene_identifiers
# when cellucid_view() was given them.
.session_gene_names <- function(view, mapping) {
  modality <- .session_modality(view, mapping)
  var <- modality$var
  if (is.null(var)) var <- data.frame(row.names = rownames(.modality_matrix(modality)))
  gene_names <- mapping[["geneNames"]]
  names <- if (identical(gene_names, "index")) {
    .make_unique_names(rownames(var))$names
  } else {
    .resolve_gene_name_set(var, substring(gene_names, 8L), "geneNames")$names
  }
  if (!is.null(view$gene_identifiers)) names <- names[match(view$gene_identifiers, names)]
  names
}

# Every key a mapping names exists on the view, and every embedding's
# columns fit it (cellucid-python's require_mapping_matches_anndata).
.require_mapping_matches_view <- function(view, mapping) {
  problems <- character()
  refs <- c(
    lapply(seq_along(mapping$embeddings$enabled), function(position) {
      list(label = paste0("embeddings.enabled[", position - 1L, "]"), ref = mapping$embeddings$enabled[[position]])
    }),
    lapply(mapping$vectorFields, function(binding) {
      list(label = paste0("vectorFields embedding of '", binding$key, "'"), ref = binding$embedding)
    })
  )
  for (item in refs) {
    key <- item$ref$key
    array <- view$obsm[[key]]
    if (is.null(array)) {
      problems <- c(problems, paste0(item$label, " uses obsm['", key, "'], which this dataset does not have"))
      next
    }
    width <- if (is.null(dim(array))) 1L else ncol(array)
    columns <- item$ref$columns
    if (width < length(columns)) {
      problems <- c(problems, paste0(
        item$label, " uses obsm['", key, "']: A ", length(columns), "D embedding needs at least ",
        length(columns), " ", .plural(length(columns), "column"), "; this array has ", width, "."
      ))
    } else if (any(columns >= width)) {
      problems <- c(problems, paste0(
        item$label, " uses obsm['", key, "']: Embedding column ", columns[columns >= width][[1L]],
        " does not exist; the array has columns 0 to ", width - 1L, "."
      ))
    }
  }
  for (binding in mapping$vectorFields) {
    if (is.null(view$obsm[[binding$key]])) {
      problems <- c(problems, paste0("vector field obsm['", binding$key, "'] is not in this dataset"))
    }
  }
  if (!is.null(mapping$knn) && !mapping$knn %in% view$obsp) {
    problems <- c(problems, paste0("knn uses obsp['", mapping$knn, "'], which this dataset does not have"))
  }
  if (!is.null(mapping$latent) && is.null(view$obsm[[mapping$latent]])) {
    problems <- c(problems, paste0("latent uses obsm['", mapping$latent, "'], which this dataset does not have"))
  }
  # An object's view holds one matrix per modality, as X: no raw, no layers.
  missing_matrix <- function(matrix) {
    if (identical(matrix, "X")) return(NULL)
    if (identical(matrix, "raw")) return("raw")
    paste0("layers['", substring(matrix, 7L), "']")
  }
  if (!is.null(mapping$expression) && !is.null(missing_matrix(mapping$expression))) {
    problems <- c(problems, paste0("expression uses ", missing_matrix(mapping$expression), ", which this dataset does not have"))
  }
  var <- .session_modality(view, mapping)$var
  columns <- if (is.null(var)) character() else names(var)
  if (!identical(mapping$geneNames, "index") && !substring(mapping$geneNames, 8L) %in% columns) {
    problems <- c(problems, paste0("geneNames uses var column '", substring(mapping$geneNames, 8L), "', which this dataset does not have"))
  }
  for (matrix in names(mapping$log1p)) {
    absent <- missing_matrix(matrix)
    if (!is.null(absent)) problems <- c(problems, paste0("log1p names ", absent, ", which this dataset does not have"))
  }
  hv <- mapping$highlyVariable
  if (!is.null(hv)) {
    if (!hv %in% columns) {
      problems <- c(problems, paste0("highlyVariable uses var column '", hv, "', which this dataset does not have"))
    } else if (!is.logical(var[[hv]])) {
      problems <- c(problems, paste0("highlyVariable uses var column '", hv, "', which is not a true/false column"))
    }
  }
  if (length(problems) > 0L) {
    .abort("The session's data mapping does not match this dataset: ", paste(problems, collapse = "; "), ".")
  }
  invisible(NULL)
}

# How many genes the mapping exposes (0 without an expression matrix).
.session_var_count <- function(view, mapping) {
  if (is.null(mapping$expression)) return(0L)
  length(.session_gene_names(view, mapping))
}

# One gene's values over every cell of the view: from the modality's matrix
# (genes x cells), NA where the modality did not measure a cell.
.session_gene_values <- function(view, mapping, gene) {
  modality <- .session_modality(view, mapping)
  matrix <- .modality_matrix(modality)
  row <- if (!is.null(view$gene_identifiers)) match(view$gene_identifiers[[gene]], rownames(matrix)) else gene
  values <- as.numeric(matrix[row, , drop = TRUE])
  if (is.null(modality$cells)) return(values)
  placed <- rep(NA_real_, view$n_cells)
  placed[match(modality$cells, view$cells)] <- values
  placed
}
