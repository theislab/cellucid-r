# The manifests, and the accounting that keeps them true. A manifest is the
# only index the viewer has, so what it declares is re-expanded and compared
# against the payloads that were actually written, on every axis, before the
# generation is published. The manifest's own bytes get the same treatment:
# every one is read back and required to parse to exactly the payload that was
# validated, because validating a payload proves nothing about the file a
# reader will actually open.

# An empty JSON object. jsonlite renders an empty *unnamed* list as `[]`, so a
# map that happens to be empty -- the schema map of an export carrying no
# observation field -- silently changes kind on the way to disk, and the
# viewer refuses it with "expected an object". Naming the empty case here
# keeps the payload, rather than jsonlite's default, the place that decides
# whether a node is an object.
.json_object <- function() {
  structure(list(), names = character(0))
}

# What one payload node becomes in JSON. jsonlite's mapping is not injective:
# `auto_unbox = TRUE` unboxes a one-element vector, an unnamed list and a
# named list are two different kinds, and a non-finite double is written as a
# string. Naming the kind is what lets the read-back compare a parsed manifest
# against the payload it came from without re-implementing the encoder.
.json_kind <- function(x) {
  if (is.null(x)) {
    return("null")
  }
  if (inherits(x, "AsIs")) {
    return("array")
  }
  if (is.list(x)) {
    return(if (is.null(names(x))) "array" else "object")
  }
  if (!is.atomic(x)) {
    return("unsupported")
  }
  # A named atomic vector still renders as an array: jsonlite drops the names.
  if (length(x) != 1L) {
    return("array")
  }
  if (is.character(x)) {
    return("string")
  }
  if (is.logical(x)) {
    return("boolean")
  }
  if (is.numeric(x)) {
    return("number")
  }
  "unsupported"
}

.json_elements <- function(x) {
  if (inherits(x, "AsIs")) {
    x <- unclass(x)
  }
  if (is.list(x)) {
    return(unname(x))
  }
  lapply(seq_along(x), function(index) x[[index]])
}

# Compare one payload node against the same node parsed back out of the file.
# Returns NULL when they agree, or the path of the first place they do not.
# Numbers are compared by value rather than by storage type, because jsonlite
# parses a whole number back as an integer whatever it was written from; a
# number that became a *string* is caught by the kind, which is the case that
# matters.
.first_json_difference <- function(expected, actual, path = "$") {
  expected_kind <- .json_kind(expected)
  actual_kind <- .json_kind(actual)
  if (identical(expected_kind, "unsupported")) {
    return(paste0(path, " is not a JSON value"))
  }
  if (!identical(expected_kind, actual_kind)) {
    return(paste0(
      path, " was written as ", actual_kind, " but validated as ", expected_kind
    ))
  }
  if (identical(expected_kind, "null")) {
    return(NULL)
  }
  if (identical(expected_kind, "object")) {
    if (!identical(names(expected), names(actual))) {
      return(paste0(
        path, " keys are [", paste(names(actual), collapse = ", "),
        "] but were validated as [", paste(names(expected), collapse = ", "), "]"
      ))
    }
    for (index in seq_along(expected)) {
      difference <- .first_json_difference(
        expected[[index]],
        actual[[index]],
        paste0(path, ".", names(expected)[[index]])
      )
      if (!is.null(difference)) {
        return(difference)
      }
    }
    return(NULL)
  }
  if (identical(expected_kind, "array")) {
    expected_elements <- .json_elements(expected)
    actual_elements <- .json_elements(actual)
    if (length(expected_elements) != length(actual_elements)) {
      return(paste0(
        path, " holds ", length(actual_elements), " elements but was validated ",
        "with ", length(expected_elements)
      ))
    }
    for (index in seq_along(expected_elements)) {
      difference <- .first_json_difference(
        expected_elements[[index]],
        actual_elements[[index]],
        paste0(path, "[", index - 1L, "]")
      )
      if (!is.null(difference)) {
        return(difference)
      }
    }
    return(NULL)
  }
  same <- if (identical(expected_kind, "number")) {
    isTRUE(as.numeric(expected) == as.numeric(actual))
  } else {
    identical(as.vector(expected), as.vector(actual))
  }
  if (!same) {
    return(paste0(
      path, " is ", format(actual), " but was validated as ", format(expected)
    ))
  }
  NULL
}

# Read one manifest back out of the staging directory and require it to say
# what the payload said. Everything upstream validates the payload; this is the
# only place that validates the file, and a generation whose manifest lost a
# shape on the way to disk is rejected here instead of publishing and failing
# in the browser.
.require_manifest_reads_back <- function(path, payload) {
  bytes <- readBin(path, "raw", file.size(path))
  text <- rawToChar(bytes)
  Encoding(text) <- "UTF-8"
  if (is.na(iconv(text, "UTF-8", "UTF-8"))) {
    stop(
      basename(path),
      " does not read back as the manifest that was validated: ",
      "the file is not valid UTF-8.",
      call. = FALSE
    )
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(text, simplifyVector = FALSE),
    error = function(condition) {
      stop(
        basename(path),
        " does not read back as the manifest that was validated: ",
        conditionMessage(condition),
        call. = FALSE
      )
    }
  )
  difference <- .first_json_difference(payload, parsed)
  if (!is.null(difference)) {
    stop(
      basename(path),
      " does not read back as the manifest that was validated: ",
      difference,
      ".",
      call. = FALSE
    )
  }
  invisible(path)
}

.write_json <- function(path, payload, pretty = FALSE) {
  json <- jsonlite::toJSON(
    payload,
    auto_unbox = TRUE,
    null = "null",
    digits = 17
  )

  if (isTRUE(pretty)) {
    json <- jsonlite::prettify(json)
  } else if (is.numeric(pretty) && length(pretty) == 1L) {
    json <- jsonlite::prettify(json, indent = as.integer(pretty))
  }

  writeLines(json, con = path, useBytes = TRUE)
  .require_manifest_reads_back(path, payload)
  invisible(path)
}

# A payload filename is the position of its field on its axis, never the
# field's own name, so no exported path depends on a dataset's vocabulary. The
# manifest entry carries the same integer as its first element, and the viewer
# substitutes it into the schema path pattern to reach the bytes.
.payload_index <- function(position) {
  if (
    !is.numeric(position) ||
      length(position) != 1L ||
      is.na(position) ||
      !is.finite(position) ||
      position != floor(position) ||
      position < 1 ||
      position > .Machine$integer.max
  ) {
    stop("Payload position must be one positive integer.", call. = FALSE)
  }
  as.integer(position) - 1L
}

# Within one axis directory the indices must be exactly 0 to N-1, each used
# once. The index *is* the filename, so two fields holding one index write into
# one file: the second overwrites the first and the viewer then draws one
# field's values under the other field's name. Nothing downstream can detect
# that, so it is asserted here in the writer, against the manifest that was
# just built, rather than only in a test.
.require_dense_payload_indices <- function(indices, axis) {
  for (position in seq_along(indices)) {
    index <- indices[[position]]
    if (
      !is.integer(index) ||
        length(index) != 1L ||
        is.na(index)
    ) {
      stop(
        axis,
        " payload index at position ",
        position - 1L,
        " must be a native integer.",
        call. = FALSE
      )
    }
  }
  resolved <- as.integer(unlist(indices, use.names = FALSE))
  if (!identical(sort(resolved), seq_along(resolved) - 1L)) {
    stop(
      axis,
      " payload indices must be exactly 0..",
      length(resolved) - 1L,
      ", each used once; got ",
      .format_value_list(sort(resolved)),
      ".",
      call. = FALSE
    )
  }
  invisible(resolved)
}

.expand_payload_pattern <- function(pattern, index, label, ext = NULL) {
  if (
    !is.character(pattern) ||
      length(pattern) != 1L ||
      is.na(pattern) ||
      !nzchar(pattern)
  ) {
    stop(label, " must be a non-empty path pattern.", call. = FALSE)
  }
  expanded <- gsub("{index}", as.character(index), pattern, fixed = TRUE)
  if (!is.null(ext)) {
    expanded <- gsub("{ext}", ext, expanded, fixed = TRUE)
  }
  if (grepl("[{}]", expanded)) {
    stop(
      label,
      " retains an unsubstituted placeholder: '",
      expanded,
      "'.",
      call. = FALSE
    )
  }
  expanded
}

# The manifest is the only index the viewer has: a payload it does not declare
# is invisible, and a payload it declares but that was never written fails the
# dataset at read time, in the browser, long after the export succeeded. Both
# are caught here by re-expanding the emitted path patterns and comparing them
# against the directory that was actually written.
#
# `directory_name = NULL` reconciles the export root itself rather than one
# axis directory. The root is where the point payloads live, and they are the
# only artifact this format declares by path from there -- in
# dataset_identity.json -- so nothing else in the package can catch a
# coordinate file whose name on disk is not the name the viewer is told to
# fetch. The root also holds the manifests and the axis directories, so
# `declared_directories` names the subdirectories that belong to the
# generation: a directory that is not one of them is refused by kind, which
# also keeps a declared payload from being satisfied by a directory of the
# same name.
.require_declared_payloads_on_disk <- function(
    out_dir,
    directory_name,
    declared,
    axis,
    declared_directories = character(0)
) {
  directory <- if (is.null(directory_name)) {
    out_dir
  } else {
    file.path(out_dir, directory_name)
  }
  prefix <- if (is.null(directory_name)) "" else paste0(directory_name, "/")
  on_disk <- character(0)
  on_disk_directories <- character(0)
  if (dir.exists(directory)) {
    for (entry in sort(list.files(directory, all.files = TRUE, no.. = TRUE))) {
      entry_path <- file.path(directory, entry)
      relative <- paste0(prefix, entry)
      if (utils::file_test("-d", entry_path)) {
        if (!relative %in% declared_directories) {
          stop(
            axis,
            " payload directory holds a non-file entry: ",
            entry_path,
            call. = FALSE
          )
        }
        on_disk_directories <- c(on_disk_directories, relative)
        next
      }
      if (!utils::file_test("-f", entry_path)) {
        stop(
          axis,
          " payload directory holds a non-file entry: ",
          entry_path,
          call. = FALSE
        )
      }
      on_disk <- c(on_disk, relative)
    }
  }
  missing_payloads <- sort(union(
    setdiff(declared, on_disk),
    setdiff(declared_directories, on_disk_directories)
  ))
  undeclared <- sort(setdiff(on_disk, declared))
  if (length(missing_payloads) > 0L || length(undeclared) > 0L) {
    stop(
      axis,
      " manifest does not describe the payloads that were written. ",
      "Declared but absent: ",
      .format_value_list(missing_payloads),
      ". Written but undeclared: ",
      .format_value_list(undeclared),
      ".",
      call. = FALSE
    )
  }
  invisible(declared)
}

.declared_obs_payload_paths <- function(manifest) {
  schemas <- manifest[["_obsSchemas"]]
  declared <- character(0)
  for (field in manifest[["_continuousFields"]]) {
    declared <- c(declared, .expand_payload_pattern(
      schemas$continuous$pathPattern,
      index = field[[1L]],
      label = "obs continuous pathPattern"
    ))
  }
  for (field in manifest[["_categoricalFields"]]) {
    extension <- unname(.codes_extension_by_dtype[field[[4L]]])
    if (is.na(extension)) {
      stop(
        "Categorical obs field '",
        field[[2L]],
        "' declares an unknown codes dtype '",
        field[[4L]],
        "'.",
        call. = FALSE
      )
    }
    declared <- c(
      declared,
      .expand_payload_pattern(
        schemas$categorical$codesPathPattern,
        index = field[[1L]],
        label = "obs categorical codesPathPattern",
        ext = extension
      ),
      .expand_payload_pattern(
        schemas$categorical$outlierPathPattern,
        index = field[[1L]],
        label = "obs categorical outlierPathPattern"
      )
    )
  }
  unique(declared)
}

.declared_var_payload_paths <- function(manifest) {
  schema <- manifest[["_varSchema"]]
  declared <- vapply(
    manifest$fields,
    function(field) {
      .expand_payload_pattern(
        schema$pathPattern,
        index = field[[1L]],
        label = "var pathPattern"
      )
    },
    character(1)
  )
  unique(declared)
}

.gene_names_from_compact_manifest <- function(manifest) {
  fields <- manifest$fields
  if (!is.list(fields)) {
    stop("compact_v1 var manifest fields must be a list.", call. = FALSE)
  }
  gene_names <- character(0)
  payload_indices <- list()
  for (field in fields) {
    if (
      !is.list(field) ||
        !(length(field) %in% c(2L, 4L)) ||
        !is.character(field[[2L]]) ||
        length(field[[2L]]) != 1L ||
        is.na(field[[2L]]) ||
        !nzchar(field[[2L]])
    ) {
      stop(
        "compact_v1 var fields must be exact [index, name] or ",
        "[index, name, minValue, maxValue] tuples.",
        call. = FALSE
      )
    }
    payload_indices[[length(payload_indices) + 1L]] <- field[[1L]]
    gene_names <- c(gene_names, field[[2L]])
  }
  .require_unique_identifiers(gene_names, what = "Gene")
  .require_dense_payload_indices(payload_indices, axis = "Gene")
  gene_names
}

.identity_obs_fields_from_compact_manifest <- function(manifest) {
  expected_names <- c(
    "_format",
    "n_points",
    "centroid_outlier_quantile",
    "latent_key",
    "compression",
    "_obsSchemas",
    "_continuousFields",
    "_categoricalFields"
  )
  if (
    !is.list(manifest) ||
      is.null(names(manifest)) ||
      length(names(manifest)) != length(expected_names) ||
      !setequal(names(manifest), expected_names)
  ) {
    stop(
      "obs manifest must contain exactly the current compact_v1 fields.",
      call. = FALSE
    )
  }
  if (!identical(manifest[["_format"]], "compact_v1")) {
    stop("obs manifest must use the current compact_v1 format.", call. = FALSE)
  }
  # The schema map is an object in every generation, including one with no
  # observation field at all, where it is empty. An empty unnamed list is a
  # JSON array once serialized, and the viewer refuses a non-object here, so
  # the kind is asserted rather than only the contents.
  if (!identical(.json_kind(manifest[["_obsSchemas"]]), "object")) {
    stop(
      "obs manifest _obsSchemas must be a JSON object.",
      call. = FALSE
    )
  }

  continuous_fields <- manifest[["_continuousFields"]]
  categorical_fields <- manifest[["_categoricalFields"]]
  if (!is.list(continuous_fields)) {
    stop(
      "compact_v1 obs manifest _continuousFields must be a list.",
      call. = FALSE
    )
  }
  if (!is.list(categorical_fields)) {
    stop(
      "compact_v1 obs manifest _categoricalFields must be a list.",
      call. = FALSE
    )
  }

  identity_fields <- list()
  manifest_keys <- character(0)
  # Both arrays write into obs/, so their payload indices share one space.
  payload_indices <- list()
  for (field in continuous_fields) {
    if (
      !is.list(field) ||
        !(length(field) %in% c(2L, 4L)) ||
        !is.character(field[[2]]) ||
        length(field[[2]]) != 1L ||
        is.na(field[[2]]) ||
        !nzchar(field[[2]])
    ) {
      stop(
        "compact_v1 continuous observation fields must be exact ",
        "[index, key] or [index, key, minValue, maxValue] tuples.",
        call. = FALSE
      )
    }
    payload_indices[[length(payload_indices) + 1L]] <- field[[1]]
    key <- field[[2]]
    manifest_keys <- c(manifest_keys, key)
    identity_fields[[length(identity_fields) + 1L]] <- list(
      key = key,
      kind = "continuous"
    )
  }

  for (field in categorical_fields) {
    if (
      !is.list(field) ||
        !(length(field) %in% c(6L, 8L)) ||
        !is.character(field[[2]]) ||
        length(field[[2]]) != 1L ||
        is.na(field[[2]]) ||
        !nzchar(field[[2]]) ||
        !(is.atomic(field[[3]]) || is.list(field[[3]])) ||
        !is.null(dim(field[[3]]))
    ) {
      stop(
        "compact_v1 categorical observation fields must be exact ",
        "six- or eight-member tuples with a category array.",
        call. = FALSE
      )
    }
    payload_indices[[length(payload_indices) + 1L]] <- field[[1]]
    key <- field[[2]]
    manifest_keys <- c(manifest_keys, key)
    identity_fields[[length(identity_fields) + 1L]] <- list(
      key = key,
      kind = "category",
      n_categories = as.integer(length(field[[3]]))
    )
  }

  .require_unique_identifiers(
    manifest_keys,
    what = "Observation field"
  )
  .require_dense_payload_indices(payload_indices, axis = "Observation")
  identity_fields
}
