# A data mapping: what a session's viewer plotted and read (the embeddings
# shown, the neighbour graph, the gene names, the expression matrix, the
# latent space, the vector fields, and for a dataset of several modalities
# which one). A session file records the mapping it was saved with, and
# cellucid_apply_session() checks it against the object before applying
# anything. The rules are the web app's (data/data-mapping.js) and
# cellucid-python's (data_mapping.validate_data_mapping), case for case: the
# shared conformance fixture tests/testthat/fixtures/data_mapping_conformance.json
# pins the reason code of every refusal (sections dataMapping and
# mappingModality).
#
# A mapping arrives as jsonlite::parse_json(simplifyVector = FALSE) parses
# JSON: an object is a named list, an array an unnamed list, null NULL.

.mapping_keys <- c(
  "kind", "embeddings", "knn", "geneNames", "expression", "log1p", "highlyVariable",
  "latent", "vectorFields"
)
.dimension_slots <- c("1d", "2d", "3d")

# A JSON object: a named list (jsonlite's `{}` is a list with empty names).
.is_json_object <- function(value) {
  is.list(value) && !is.null(names(value)) && !is.data.frame(value)
}

.is_json_array <- function(value) {
  is.list(value) && is.null(names(value))
}

# A JSON integer: jsonlite reads one that fits as an integer and a larger one
# as a whole double. A boolean is not one, and neither is a fraction.
.is_json_integer <- function(value) {
  length(value) == 1L && !is.logical(value) && (
    is.integer(value) && !is.na(value) ||
      is.double(value) && is.finite(value) && value == floor(value) && abs(value) > .Machine$integer.max
  )
}

.is_json_string <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value)
}

.require_mapping_keys <- function(value, expected, label) {
  if (!.is_json_object(value)) {
    stop(.data_mapping_error("not_object", paste0(label, " must be an object.")))
  }
  missing <- .sort_code_points(setdiff(expected, names(value)))
  unknown <- .sort_code_points(setdiff(names(value), expected))
  if (length(missing) > 0L || length(unknown) > 0L) {
    details <- c(
      if (length(missing) > 0L) paste0("missing ", paste(missing, collapse = ", ")),
      if (length(unknown) > 0L) paste0("unknown ", paste(unknown, collapse = ", "))
    )
    stop(.data_mapping_error("keys", paste0(label, " fields: ", paste(details, collapse = "; "), ".")))
  }
  value
}

.require_mapping_string <- function(value, label) {
  if (!.is_json_string(value) || !nzchar(value)) {
    stop(.data_mapping_error("not_string", paste0(label, " must be a non-empty string.")))
  }
  value
}

.optional_mapping_string <- function(value, label) {
  if (is.null(value)) NULL else .require_mapping_string(value, label)
}

# Strings in code-point order, as Python and JavaScript compare them (never
# the locale's collation): -1, 0 or 1.
.compare_code_points <- function(a, b) {
  x <- .code_points(a)
  y <- .code_points(b)
  shared <- min(length(x), length(y))
  if (shared > 0L) {
    differ <- which(x[seq_len(shared)] != y[seq_len(shared)])
    if (length(differ) > 0L) {
      first <- differ[[1L]]
      return(if (x[[first]] < y[[first]]) -1L else 1L)
    }
  }
  sign(length(x) - length(y))
}

# One embedding reference's sort key compared with another's: key in
# code-point order, then columns element by element (a shorter prefix first).
.compare_embedding_refs <- function(a, b) {
  by_key <- .compare_code_points(a[["key"]], b[["key"]])
  if (by_key != 0L) return(by_key)
  x <- unlist(a[["columns"]])
  y <- unlist(b[["columns"]])
  shared <- min(length(x), length(y))
  if (shared > 0L) {
    differ <- which(x[seq_len(shared)] != y[seq_len(shared)])
    if (length(differ) > 0L) {
      first <- differ[[1L]]
      return(if (x[[first]] < y[[first]]) -1L else 1L)
    }
  }
  sign(length(x) - length(y))
}

# Validate one {key, columns} reference; returns it with integer columns.
.validate_embedding_ref <- function(value, dimension, label) {
  ref <- .require_mapping_keys(value, c("key", "columns"), label)
  key <- .require_mapping_string(ref[["key"]], paste0(label, ".key"))
  columns <- ref[["columns"]]
  if (!.is_json_array(columns)) {
    stop(.data_mapping_error("columns_not_list", paste0(label, ".columns must be a list.")))
  }
  if (!is.null(dimension) && length(columns) != dimension) {
    stop(.data_mapping_error("column_count", paste0(label, ".columns must list exactly ", dimension, " ", .plural(dimension, "column"), ".")))
  }
  if (length(columns) < 1L || length(columns) > 3L) {
    stop(.data_mapping_error("column_count", paste0(label, ".columns must list 1 to 3 columns.")))
  }
  for (column in columns) {
    if (!.is_json_integer(column) || column < 0) {
      stop(.data_mapping_error("column_not_integer", paste0(label, ".columns must hold non-negative integers.")))
    }
  }
  columns <- vapply(columns, as.integer, integer(1L))
  if (anyDuplicated(columns)) {
    stop(.data_mapping_error("column_duplicate", paste0(label, ": each axis needs a different column.")))
  }
  list(key = key, columns = columns)
}

.is_expression_matrix <- function(value) {
  .is_json_string(value) && (value %in% c("X", "raw") || (startsWith(value, "layer:") && nchar(value) > 6L))
}

.validate_embeddings_group <- function(value, label) {
  group <- .require_mapping_keys(value, c("enabled", "defaults"), label)
  if (!.is_json_array(group[["enabled"]])) {
    stop(.data_mapping_error("not_list", paste0(label, ".enabled must be a list.")))
  }
  if (length(group[["enabled"]]) == 0L) {
    stop(.data_mapping_error("no_embedding", paste0(label, ".enabled must list at least one embedding.")))
  }
  enabled <- lapply(seq_along(group[["enabled"]]), function(position) {
    .validate_embedding_ref(group[["enabled"]][[position]], NULL, paste0(label, ".enabled[", position - 1L, "]"))
  })
  for (index in seq_along(enabled)[-1L]) {
    if (.compare_embedding_refs(enabled[[index - 1L]], enabled[[index]]) > 0L) {
      stop(.data_mapping_error("embedding_order", paste0(label, ".enabled must be sorted by embedding key, then columns.")))
    }
  }
  for (index in seq_along(enabled)[-1L]) {
    if (.compare_embedding_refs(enabled[[index - 1L]], enabled[[index]]) == 0L) {
      stop(.data_mapping_error("embedding_duplicate", paste0(label, ".enabled repeats an embedding.")))
    }
  }
  slots <- .require_mapping_keys(group[["defaults"]], .dimension_slots, paste0(label, ".defaults"))
  defaults <- stats::setNames(vector("list", 3L), .dimension_slots)
  for (slot in .dimension_slots) {
    ref <- slots[[slot]]
    if (is.null(ref)) next
    valid <- .validate_embedding_ref(ref, as.integer(substr(slot, 1L, 1L)), paste0(label, ".defaults.", slot))
    known <- any(vapply(enabled, function(item) .compare_embedding_refs(item, valid) == 0L, logical(1L)))
    if (!known) {
      stop(.data_mapping_error("default_not_enabled", paste0(label, ".defaults.", slot, " must be one of the enabled embeddings.")))
    }
    defaults[slot] <- list(valid)
  }
  if (all(vapply(defaults, is.null, logical(1L)))) {
    stop(.data_mapping_error("no_embedding", paste0(label, ".defaults must name at least one default.")))
  }
  list(enabled = enabled, defaults = defaults)
}

.validate_log1p <- function(value, label) {
  if (!.is_json_object(value)) {
    stop(.data_mapping_error("log1p", paste0(label, " must be an object of expression matrices.")))
  }
  out <- list()
  matrices <- names(value)
  for (matrix in .sort_code_points(matrices)) {
    if (!.is_expression_matrix(matrix)) {
      stop(.data_mapping_error("log1p", paste0(label, ' keys must be expression matrices: "X", "raw" or "layer:<name>".')))
    }
    entry <- .require_mapping_keys(value[[matrix]], "base", paste0(label, ".", matrix))
    base <- entry[["base"]]
    if (!is.null(base) && !(is.numeric(base) && length(base) == 1L && is.finite(base) && base > 0 && base != 1)) {
      stop(.data_mapping_error("log1p_base", paste0(label, ".", matrix, ".base must be null (base e) or a positive number other than 1.")))
    }
    out[[matrix]] <- list(base = base)
  }
  if (length(out) == 0L) structure(list(), names = character()) else out
}

# Validate one data mapping and return its canonical form; a refusal is a
# cellucid_data_mapping_error carrying the reason code the conformance
# fixture pins.
.validate_data_mapping <- function(value, label = "dataMapping") {
  if (!.is_json_object(value)) {
    stop(.data_mapping_error("not_object", paste0(label, " must be an object.")))
  }
  if (!identical(value[["kind"]], "anndata")) {
    stop(.data_mapping_error("kind", paste0(label, '.kind must be "anndata".')))
  }
  has_modality <- "modality" %in% names(value)
  mapping <- .require_mapping_keys(value, c(.mapping_keys, if (has_modality) "modality"), label)
  embeddings <- .validate_embeddings_group(mapping[["embeddings"]], paste0(label, ".embeddings"))

  if (!.is_json_array(mapping[["vectorFields"]])) {
    stop(.data_mapping_error("not_list", paste0(label, ".vectorFields must be a list.")))
  }
  bindings <- lapply(seq_along(mapping[["vectorFields"]]), function(position) {
    binding_label <- paste0(label, ".vectorFields[", position - 1L, "]")
    binding <- .require_mapping_keys(mapping[["vectorFields"]][[position]], c("key", "embedding"), binding_label)
    list(
      key = .require_mapping_string(binding[["key"]], paste0(binding_label, ".key")),
      embedding = .validate_embedding_ref(binding[["embedding"]], NULL, paste0(binding_label, ".embedding"))
    )
  })
  compare_bindings <- function(a, b) {
    by_embedding <- .compare_embedding_refs(a[["embedding"]], b[["embedding"]])
    if (by_embedding != 0L) by_embedding else .compare_code_points(a[["key"]], b[["key"]])
  }
  for (index in seq_along(bindings)[-1L]) {
    if (compare_bindings(bindings[[index - 1L]], bindings[[index]]) > 0L) {
      stop(.data_mapping_error("vector_order", paste0(
        label, ".vectorFields must be sorted by embedding key, embedding columns, then key."
      )))
    }
  }
  for (index in seq_along(bindings)[-1L]) {
    if (compare_bindings(bindings[[index - 1L]], bindings[[index]]) == 0L) {
      stop(.data_mapping_error("vector_duplicate", paste0(label, ".vectorFields repeats a binding.")))
    }
  }
  for (position in seq_along(bindings)) {
    bound <- any(vapply(embeddings$enabled, function(ref) {
      .compare_embedding_refs(ref, bindings[[position]]$embedding) == 0L
    }, logical(1L)))
    if (!bound) {
      stop(.data_mapping_error("vector_not_enabled", paste0(
        label, ".vectorFields[", position - 1L, "] is bound to an embedding that is not enabled."
      )))
    }
  }

  knn <- .optional_mapping_string(mapping[["knn"]], paste0(label, ".knn"))
  modality <- NULL
  if (has_modality) {
    modality <- mapping[["modality"]]
    if (!.is_json_string(modality) || !nzchar(modality)) {
      stop(.data_mapping_error("modality", paste0(
        label, ".modality must name one of the dataset's modalities; a dataset with one feature ",
        "space has no modality key."
      )))
    }
  }
  gene_names <- mapping[["geneNames"]]
  if (!(identical(gene_names, "index") ||
        (.is_json_string(gene_names) && startsWith(gene_names, "column:") && nchar(gene_names) > 7L))) {
    stop(.data_mapping_error("gene_names", paste0(label, '.geneNames must be "index" or "column:<var column>".')))
  }
  expression <- mapping[["expression"]]
  if (!is.null(expression) && !.is_expression_matrix(expression)) {
    stop(.data_mapping_error("expression", paste0(label, '.expression must be "X", "raw", "layer:<name>" or null.')))
  }
  log1p <- .validate_log1p(mapping[["log1p"]], paste0(label, ".log1p"))
  highly_variable <- .optional_mapping_string(mapping[["highlyVariable"]], paste0(label, ".highlyVariable"))
  latent <- .optional_mapping_string(mapping[["latent"]], paste0(label, ".latent"))

  canonical <- list(kind = "anndata", embeddings = embeddings, knn = knn)
  if (has_modality) canonical$modality <- modality
  canonical <- c(canonical, list(
    geneNames = gene_names, expression = expression, log1p = log1p,
    highlyVariable = highly_variable, latent = latent, vectorFields = bindings
  ))
  canonical
}

# A mapping's modality checked against the dataset's modalities (NULL for a
# dataset of one feature space): codes modality_unexpected, modality_missing,
# modality_unknown.
.require_mapping_modality <- function(mapping, modalities) {
  named <- "modality" %in% names(mapping)
  if (is.null(modalities)) {
    if (named) {
      stop(.data_mapping_error("modality_unexpected", paste0(
        "This mapping reads the modality '", mapping[["modality"]], "', and this dataset has one feature space."
      )))
    }
    return(invisible(NULL))
  }
  if (!named) {
    stop(.data_mapping_error("modality_missing", paste0(
      "This dataset has several modalities (", paste(modalities, collapse = ", "),
      "); its mapping must name the one it reads."
    )))
  }
  if (!mapping[["modality"]] %in% modalities) {
    stop(.data_mapping_error("modality_unknown", paste0(
      "This mapping reads the modality '", mapping[["modality"]], "', which this dataset does not have (",
      paste(modalities, collapse = ", "), ")."
    )))
  }
  invisible(NULL)
}
