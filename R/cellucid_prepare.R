# The help page for cellucid_prepare() is hand-written in
# man/cellucid_prepare.Rd and is the only source for it. Nothing in this
# package is generated: NAMESPACE is hand-written too, because
# useDynLib(cellucid, .registration = TRUE, .fixes = "C_") has no generator
# input anywhere under R/. Edit the .Rd file directly when this signature
# changes; tests/testthat/test-current-contract.R rebuilds the documented
# \usage block into a function and compares it to these formals, so the two
# cannot drift apart.
cellucid_prepare <- function(
    out_path,
    format,
    embeddings,
    dataset_id,
    dataset_name,
    obs = NULL,
    var = NULL,
    gene_expression = NULL,
    connectivities = NULL,
    vector_fields = NULL,
    gene_identifiers = NULL,
    dataset_description = NULL,
    citation_text = NULL,
    citation_url = NULL,
    source_name = NULL,
    source_url = NULL,
    embedding_defaults = NULL,
    latent_key = NULL,
    obs_keys = NULL,
    var_gene_id_column = NULL,
    vector_field_default = NULL,
    normalize_embeddings = NULL,
    centroid_min_points = NULL,
    zarr_format = NULL,
    force = FALSE
) {
  # An argument that was never supplied and an argument supplied with the wrong
  # value are two different mistakes, and only the message tells the caller
  # which one they made. Every required argument is therefore reported here,
  # before any value is inspected. Each `missing()` sits directly in an `if`
  # condition in this body, which is the one place R guarantees it reads the
  # argument's own supplied-ness.
  absent_required <- character()
  if (missing(out_path)) {
    absent_required <- c(absent_required, "out_path")
  }
  if (missing(format)) {
    absent_required <- c(absent_required, "format")
  }
  if (missing(embeddings)) {
    absent_required <- c(absent_required, "embeddings")
  }
  if (missing(dataset_id)) {
    absent_required <- c(absent_required, "dataset_id")
  }
  if (missing(dataset_name)) {
    absent_required <- c(absent_required, "dataset_name")
  }
  .require_supplied_arguments(absent_required)

  # --- Scalars: everything the caller can get wrong in one value ------------
  format <- .validate_format(format)
  zarr_format <- .validate_zarr_format(zarr_format)
  force <- .validate_flag(force, "force")
  target <- .resolve_output_path(out_path, format = format, force = force)
  target$zarr_format <- zarr_format
  dataset_id <- .validate_dataset_id(dataset_id)
  dataset_name <- .validate_dataset_name(dataset_name)
  dataset_description <- .validate_optional_text(
    dataset_description, "dataset_description", allow_empty = TRUE
  )
  citation <- .validate_linked_text(
    citation_text, citation_url, "citation_text", "citation_url", "text"
  )
  source <- .validate_linked_text(
    source_name, source_url, "source_name", "source_url", "name"
  )

  # --- obsm: the embeddings and the vector fields ---------------------------
  resolved <- .resolve_embeddings(embeddings)
  matrices <- resolved$matrices
  n_cells <- resolved$n_cells
  vectors <- .resolve_vector_fields(vector_fields, matrices, n_cells)

  # --- obs ------------------------------------------------------------------
  if (is.null(obs)) {
    obs <- data.frame(row = seq_len(n_cells))[, 0L, drop = FALSE]
  } else if (!is.data.frame(obs)) {
    stop("obs must be a data.frame with one row per cell.", call. = FALSE)
  }
  if (nrow(obs) != n_cells) {
    stop(
      "obs has ", nrow(obs), " ", .plural(nrow(obs), "row", "rows"),
      ", but the embeddings have ", n_cells, " ", .plural(n_cells, "cell", "cells"), ".",
      call. = FALSE
    )
  }
  .validate_frame_columns(obs, "obs column")

  # --- X and var ------------------------------------------------------------
  if (is.null(gene_expression)) {
    supplied <- c(
      var = !is.null(var),
      gene_identifiers = !is.null(gene_identifiers),
      var_gene_id_column = !is.null(var_gene_id_column)
    )
    if (any(supplied)) {
      stop(
        names(supplied)[supplied][[1L]], " describes the genes of ",
        "gene_expression, and no gene_expression was supplied.",
        call. = FALSE
      )
    }
    expression <- NULL
    var <- data.frame(row = integer(0))[, 0L, drop = FALSE]
    var_index <- character(0)
  } else {
    expression <- .as_expression_matrix(gene_expression, "gene_expression")
    if (nrow(expression) != n_cells) {
      stop(
        "gene_expression has ", nrow(expression), " ",
        .plural(nrow(expression), "row", "rows"), ", but the embeddings have ",
        n_cells, " ", .plural(n_cells, "cell", "cells"), ". gene_expression is ",
        "cells x genes; transpose a genes x cells matrix with t().",
        call. = FALSE
      )
    }
    genes <- .resolve_genes(expression, var, var_gene_id_column, gene_identifiers)
    expression <- genes$expression
    var <- genes$var
    var_index <- genes$index
    # The dtype is decided on the genes stored, after gene_identifiers chose
    # them, as the Python writer decides it.
    expression_dtype <- .expression_dtype(
      if (inherits(expression, "sparseMatrix")) expression@x else expression
    )
  }

  cell_sources <- c(
    list(obs = .supplied_row_names(obs)),
    if (!is.null(expression)) list(gene_expression = rownames(expression)),
    `names<-`(
      lapply(matrices, rownames),
      paste0("embedding '", names(matrices), "'")
    )
  )
  obs_index <- .resolve_cell_names(cell_sources, n_cells)

  # --- obsp -----------------------------------------------------------------
  graph <- NULL
  if (!is.null(connectivities)) {
    graph <- .resolve_connectivities(connectivities, n_cells)
  }

  # --- uns['cellucid'] ------------------------------------------------------
  defaults <- .resolve_defaults(
    matrices = matrices,
    obs = obs,
    vectors = vectors,
    embedding_defaults = embedding_defaults,
    latent_key = latent_key,
    obs_keys = obs_keys,
    var_gene_id_column = var_gene_id_column,
    vector_field_default = vector_field_default,
    normalize_embeddings = normalize_embeddings,
    centroid_min_points = centroid_min_points
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

  # --- Write the store, then move it into place ------------------------------
  .write_output(target, write = function(store) {
    .zarr_group(store, "", .encoding("anndata", "0.1.0"))
    if (!is.null(expression)) {
      if (inherits(expression, "sparseMatrix")) {
        .ad_sparse(
          store, "X", "csc_matrix", dim(expression),
          pointers = expression@p, indices = expression@i, values = expression@x,
          dtype = expression_dtype, sliced = TRUE
        )
      } else {
        .ad_dense(store, "X", expression, gene_major = TRUE, dtype = expression_dtype)
      }
    }
    .ad_dataframe(store, "obs", obs, obs_index, "obs")
    .ad_dataframe(store, "var", var, var_index, "var")
    .ad_dict(store, "obsm")
    for (key in names(matrices)) {
      .ad_dense(store, .zarr_key("obsm", key), matrices[[key]])
    }
    for (key in names(vectors$arrays)) {
      .ad_dense(store, .zarr_key("obsm", key), vectors$arrays[[key]])
    }
    .ad_dict(store, "obsp")
    if (!is.null(graph)) {
      .ad_sparse(
        store, "obsp/connectivities", "csr_matrix", dim(graph),
        pointers = graph@p, indices = graph@i, values = graph@x
      )
    }
    .ad_dict(store, "varm")
    .ad_dict(store, "varp")
    .ad_dict(store, "layers")
    .ad_uns_value(store, "uns", list(title = .uns_string(dataset_name), cellucid = cellucid))
  })
  invisible(target$path)
}

# The genes: the expression matrix narrowed to gene_identifiers, and the var
# data frame and index that describe its columns.
.resolve_genes <- function(expression, var, var_gene_id_column, gene_identifiers) {
  n_genes <- ncol(expression)
  if (is.null(var)) {
    var <- data.frame(row = seq_len(n_genes))[, 0L, drop = FALSE]
  } else if (!is.data.frame(var)) {
    stop("var must be a data.frame with one row per gene.", call. = FALSE)
  }
  if (nrow(var) != n_genes) {
    stop(
      "var has ", nrow(var), " ", .plural(nrow(var), "row", "rows"),
      ", but gene_expression has ", n_genes, " ", .plural(n_genes, "gene", "genes"),
      " (columns).",
      call. = FALSE
    )
  }
  .validate_frame_columns(var, "var column")
  if (!is.null(var_gene_id_column)) {
    var_gene_id_column <- .validate_gene_name_source(var_gene_id_column, "var_gene_id_column")
  }

  row_names <- .supplied_row_names(var)
  column_names <- colnames(expression)
  if (!is.null(row_names) && !is.null(column_names) && !identical(row_names, column_names)) {
    differs <- which(row_names != column_names)[[1L]]
    stop(
      "rownames(var) and colnames(gene_expression) disagree at gene ", differs,
      " ('", row_names[[differs]], "' and '", column_names[[differs]], "'). ",
      "var must describe the columns of gene_expression in their order.",
      call. = FALSE
    )
  }
  index <- if (!is.null(row_names)) row_names else column_names
  if (is.null(index)) {
    if (is.null(var_gene_id_column)) {
      stop(
        "The genes have no names: set rownames(var) or ",
        "colnames(gene_expression) to the gene identifiers, or pass ",
        "var_gene_id_column.",
        call. = FALSE
      )
    }
    index <- .resolve_gene_name_set(var, var_gene_id_column, "var_gene_id_column")$names
  }
  if (anyNA(index) || any(!nzchar(index))) {
    stop("Gene names must not be missing or empty.", call. = FALSE)
  }
  unique_index <- .make_unique_names(as.character(index))
  if (length(unique_index$renamed) > 0L) {
    message("Gene names: ", .unique_names_note(0, unique_index$renamed, NULL))
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
      stop("gene_identifiers must not contain duplicate identifiers.", call. = FALSE)
    }
    positions <- match(gene_identifiers, default_names)
    missing_genes <- gene_identifiers[is.na(positions)]
    if (length(missing_genes) > 0L) {
      shown <- utils::head(missing_genes, 5L)
      remainder <- length(missing_genes) - length(shown)
      stop(
        "gene_identifiers contains identifiers not found in var: ",
        .format_value_list(shown),
        if (remainder > 0L) paste0(" and ", remainder, " more") else "",
        ". Identifiers are matched against the gene names the viewer shows ",
        "(var_gene_id_column, or the var row names), after repeated names ",
        "were made unique.",
        call. = FALSE
      )
    }
    expression <- expression[, positions, drop = FALSE]
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
      message("Gene names from var column '", var_gene_id_column, "': ", note)
    }
  }
  rownames(var) <- NULL
  list(expression = expression, var = var, index = index)
}

# The launch defaults of uns['cellucid']['defaults']: only what the caller
# set, so every other choice stays the viewer's own.
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
    centroid_min_points
) {
  defaults <- .resolve_embedding_defaults(embedding_defaults, matrices)
  if (!is.null(var_gene_id_column)) {
    defaults$geneIdColumn <- .uns_string(var_gene_id_column)
  }
  if (!is.null(latent_key)) {
    latent_key <- .validate_string(latent_key, "latent_key")
    if (!latent_key %in% names(matrices)) {
      stop(
        "latent_key must name one of the embeddings: ",
        .format_value_list(names(matrices)), ".",
        call. = FALSE
      )
    }
    defaults$latentKey <- .uns_string(latent_key)
  }
  if (!is.null(vector_field_default)) {
    vector_field_default <- .validate_string(vector_field_default, "vector_field_default")
    if (length(vectors$fields) == 0L) {
      stop("vector_field_default names a field, and no vector_fields were supplied.", call. = FALSE)
    }
    if (!vector_field_default %in% vectors$fields) {
      stop(
        "vector_field_default must name one of vector_fields: ",
        .format_value_list(vectors$fields), ".",
        call. = FALSE
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
      stop(
        "vector_field_default names '", vector_field_default, "', and the viewer ",
        "shows none of its embeddings when the file opens, so it is not drawn. ",
        "Show one of them (embedding_defaults), or choose a field drawn on a ",
        "shown embedding: ", .format_value_list(sort(drawn)), ".",
        call. = FALSE
      )
    }
    defaults$vectorFieldDefault <- .uns_string(vector_field_default)
  }
  if (!is.null(obs_keys)) {
    .validate_character_vector(obs_keys, "obs_keys")
    if (length(obs_keys) == 0L) {
      stop("obs_keys must name at least one obs column.", call. = FALSE)
    }
    if (anyDuplicated(obs_keys)) {
      stop("obs_keys names '", obs_keys[duplicated(obs_keys)][[1L]], "' twice.", call. = FALSE)
    }
    absent <- setdiff(obs_keys, names(obs))
    if (length(absent) > 0L) {
      stop(
        "obs_keys contain columns not in obs: ", .format_value_list(absent),
        ". Available columns: ", .format_value_list(names(obs)), ".",
        call. = FALSE
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
    "embedding1d", "embedding2d", "embedding3d", "geneIdColumn", "latentKey",
    "vectorFieldDefault", "obsKeys", "normalizeEmbeddings", "centroidMinPoints"
  )
  defaults[intersect(order, names(defaults))]
}
