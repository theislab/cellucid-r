# Reading a .cellucid-session file (session format 3), the file the viewer
# saves and cellucid_session_result() collects: the web app's framing
# (app/session/bundle/format.js) and cellucid-python's reader
# (session_bundle.py), check for check and message for message.
#
#   CELLUCID_SESSION/3\n          the signature
#   u32 LE, manifest JSON         createdAt, datasetFingerprint, source,
#                                 requires, chunks (each chunk's metadata),
#                                 cells (a subsample), genes (a gene sample)
#   [u32 LE, chunk bytes] ...     one frame per chunk, in the manifest's order
#
# The manifest is read and checked in full when the file is opened; a chunk's
# bytes are read, inflated and parsed only when asked for, so a session that
# holds analysis results of hundreds of megabytes costs only what is used.

.session_format_version <- 3L
.session_signature_prefix <- "CELLUCID_SESSION/"
.session_magic <- charToRaw(paste0(.session_signature_prefix, .session_format_version, "\n"))
.not_a_session_file_message <- paste0(
  "This file is not a Cellucid session file of format ", .session_format_version,
  " (it does not start with CELLUCID_SESSION/", .session_format_version, ")."
)

# The web app's sentence for a session file of another format version
# (bundle/format.js otherSessionFormatMessage, cellucid-python's
# other_session_format_message): saved by an earlier Cellucid, whose cell
# check this version cannot verify, or by a newer one.
.other_session_format_message <- function(version) {
  if (version < .session_format_version) {
    paste0(
      "This session was saved by an earlier version of Cellucid (session format ", version, "): ",
      "the check that it marks the same cells cannot be verified, so it was not opened. ",
      "Re-create it in this version of Cellucid."
    )
  } else {
    paste0(
      "This session was saved by a newer version of Cellucid (session format ", version, "), ",
      "which this version cannot read. Open it in that version of Cellucid."
    )
  }
}

# Why a file whose first bytes are `magic` (as many as this format's
# signature holds, fewer for a shorter file) is not read: another session
# format's signature, or none (bundle/format.js sessionSignatureRefusal).
.session_signature_refusal <- function(magic) {
  text <- rawToChar(magic[magic != as.raw(0L)])
  prefix <- .session_signature_prefix
  rest <- if (startsWith(text, prefix)) substring(text, nchar(prefix) + 1L) else NA_character_
  digits <- if (is.na(rest)) NA_character_ else regmatches(rest, regexpr("^[0-9]+(?=\n|$)", rest, perl = TRUE))
  version <- if (length(digits) == 1L && !is.na(digits)) as.numeric(digits) else NA_real_
  if (is.na(version) || version == 0 || version == .session_format_version) {
    .not_a_session_file_message
  } else {
    .other_session_format_message(version)
  }
}
.session_max_frame_bytes <- 4294967295
.session_max_safe_integer <- 9007199254740991

.session_manifest_error <- function(...) {
  .abort("Invalid session manifest (", ..., ")")
}

.read_u32 <- function(con) {
  bytes <- readBin(con, "raw", 4L)
  if (length(bytes) != 4L) .abort("Unexpected EOF while reading session bundle")
  sum(as.numeric(bytes) * c(1, 256, 65536, 16777216))
}

.require_session_integer <- function(value, label, maximum) {
  if (!.is_json_integer(value)) .session_manifest_error(label, " must be an integer")
  if (value < 0) .session_manifest_error(label, " must be non-negative")
  if (value > maximum) .session_manifest_error(label, " exceeds the ", format(maximum, scientific = FALSE), "-byte limit")
  as.numeric(value)
}

.session_key_report <- function(present, expected) {
  missing <- sort(setdiff(expected, present), method = "radix")
  unknown <- sort(setdiff(present, expected), method = "radix")
  paste(c(
    if (length(missing) > 0L) paste0("missing ", paste(missing, collapse = ", ")),
    if (length(unknown) > 0L) paste0("unknown ", paste(unknown, collapse = ", "))
  ), collapse = "; ")
}

.is_trimmed_text <- function(value) {
  # Python's str.strip(): Unicode white space, and the separators \x1c-\x1f.
  .is_json_string(value) && nzchar(value) &&
    !grepl("(*UCP)^[\\s\\x{1c}-\\x{1f}]|[\\s\\x{1c}-\\x{1f}]$", value, perl = TRUE)
}

# The fingerprint that ties a session to one dataset and data mapping;
# returns it with its mapping in canonical form.
.validate_session_fingerprint <- function(value) {
  if (!.is_json_object(value)) .session_manifest_error("datasetFingerprint must be an object")
  expected <- c("datasetId", "name", "cellCount", "varCount", "cellOrder", "dataMapping")
  if (!setequal(names(value), expected)) {
    .session_manifest_error("datasetFingerprint fields: ", .session_key_report(names(value), expected))
  }
  for (key in c("datasetId", "name")) {
    item <- value[[key]]
    if (!is.null(item) && !.is_trimmed_text(item)) {
      .session_manifest_error("datasetFingerprint.", key, " must be a non-empty trimmed string or null")
    }
  }
  for (key in c("cellCount", "varCount")) {
    .require_session_integer(value[[key]], paste0("datasetFingerprint.", key), 4294967295)
  }
  mapping <- tryCatch(
    .validate_data_mapping(value[["dataMapping"]], "datasetFingerprint.dataMapping"),
    cellucid_data_mapping_error = function(error) .session_manifest_error(conditionMessage(error))
  )
  order <- value[["cellOrder"]]
  if (!.is_json_object(order)) .session_manifest_error("datasetFingerprint.cellOrder must be an object")
  if (!setequal(names(order), c("embedding", "dimension", "digest"))) {
    .session_manifest_error("datasetFingerprint.cellOrder fields: ", .session_key_report(names(order), c("embedding", "dimension", "digest")))
  }
  dimension <- order[["dimension"]]
  if (!is.integer(dimension) || length(dimension) != 1L || is.na(dimension) || !dimension %in% 1:3) {
    .session_manifest_error("datasetFingerprint.cellOrder.dimension must be exactly 1, 2, or 3")
  }
  digest <- order[["digest"]]
  if (!.is_json_string(digest) || !grepl("^[0-9a-f]{16}$", digest)) {
    .session_manifest_error("datasetFingerprint.cellOrder.digest must be 16 lowercase hexadecimal characters")
  }
  embedding <- tryCatch(
    .validate_embedding_ref(order[["embedding"]], dimension, "datasetFingerprint.cellOrder.embedding"),
    cellucid_data_mapping_error = function(error) .session_manifest_error(conditionMessage(error))
  )
  value[["cellOrder"]] <- list(embedding = embedding, dimension = dimension, digest = digest)
  value[["dataMapping"]] <- mapping
  value
}

.validate_session_text <- function(value, label) {
  if (!.is_trimmed_text(value) || grepl("[\001-\037\177]", value)) {
    .session_manifest_error(label, " must be a non-empty, unpadded string")
  }
  value
}

.validate_session_source <- function(value) {
  if (is.null(value)) return(NULL)
  keys <- list(
    sample = c("kind", "dataset"), link = c("kind", "locator", "dataset"),
    remote = c("kind", "url", "dataset"), file = c("kind", "name", "dataset")
  )
  kind <- if (.is_json_object(value)) value[["kind"]] else NULL
  if (!.is_json_string(kind) || !kind %in% names(keys)) {
    .session_manifest_error('source kind must be "sample", "link", "remote" or "file"')
  }
  expected <- keys[[kind]]
  if (!setequal(names(value), expected) || anyDuplicated(names(value))) {
    .session_manifest_error("source fields must be exactly ", paste(sort(expected, method = "radix"), collapse = ", "))
  }
  for (key in sort(setdiff(expected, "kind"), method = "radix")) {
    .validate_session_text(value[[key]], paste0("source.", key))
  }
  if (identical(kind, "remote") && !grepl("^https?://", value[["url"]])) {
    .session_manifest_error("source.url must be an http(s) address")
  }
  value
}

# A string's UTF-16 code units: JavaScript's string order compares these.
.utf16_units <- function(text) {
  points <- .code_points(text)
  astral <- points > 0xFFFF
  if (!any(astral)) return(points)
  units <- integer()
  for (point in points) {
    if (point > 0xFFFF) {
      offset <- point - 0x10000
      units <- c(units, 0xD800 + offset %/% 1024, 0xDC00 + offset %% 1024)
    } else {
      units <- c(units, point)
    }
  }
  units
}

.utf16_less <- function(a, b) {
  x <- .utf16_units(a)
  y <- .utf16_units(b)
  shared <- min(length(x), length(y))
  if (shared > 0L) {
    differ <- which(x[seq_len(shared)] != y[seq_len(shared)])
    if (length(differ) > 0L) return(x[[differ[[1L]]]] < y[[differ[[1L]]]])
  }
  length(x) < length(y)
}

# The dataset fields the session uses: obs columns and genes, each list
# sorted as JavaScript sorts strings and without repeats.
.validate_session_requires <- function(value) {
  if (!.is_json_object(value) || !setequal(names(value), c("obs", "var")) || length(value) != 2L) {
    .session_manifest_error("requires must have exactly obs and var")
  }
  out <- list()
  for (source in c("obs", "var")) {
    keys <- value[[source]]
    if (!.is_json_array(keys)) .session_manifest_error("requires.", source, " must be an array")
    for (index in seq_along(keys)) {
      key <- keys[[index]]
      if (!.is_trimmed_text(key)) {
        .session_manifest_error("requires.", source, "[", index - 1L, "] must be a non-empty trimmed string")
      }
      if (index > 1L && !.utf16_less(keys[[index - 1L]], key)) {
        .session_manifest_error("requires.", source, " must be sorted and without duplicates")
      }
    }
    out[[source]] <- vapply(keys, identity, character(1L))
  }
  out
}

.session_chunk_keys <- c(
  "id", "contributorId", "priority", "kind", "codec", "label", "datasetDependent",
  "storedBytes", "uncompressedBytes"
)

.validate_session_manifest <- function(manifest) {
  if (!.is_json_object(manifest)) .session_manifest_error("expected an object")
  # `cells` is present only for a session saved on a subsample (which cells
  # were open: R/subsample.R, .session_cells()); `genes` only for one saved
  # while genes were read for a sample of the open cells (.session_genes()).
  expected <- c("createdAt", "datasetFingerprint", "source", "requires", "chunks",
                if ("cells" %in% names(manifest)) "cells",
                if ("genes" %in% names(manifest)) "genes")
  if (!setequal(names(manifest), expected) || anyDuplicated(names(manifest))) {
    .abort("Invalid session manifest fields (", .session_key_report(names(manifest), expected), ")")
  }
  created <- manifest[["createdAt"]]
  if (!.is_json_string(created) || !grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}\\.\\d{3}Z$", created, perl = TRUE)) {
    .session_manifest_error("createdAt must be an exact UTC millisecond timestamp")
  }
  manifest[["datasetFingerprint"]] <- .validate_session_fingerprint(manifest[["datasetFingerprint"]])
  .validate_session_source(manifest[["source"]])
  .validate_session_requires(manifest[["requires"]])
  if ("cells" %in% names(manifest)) {
    .session_cells(manifest[["cells"]], manifest[["datasetFingerprint"]][["cellCount"]])
  }
  if ("genes" %in% names(manifest)) {
    .session_genes(manifest[["genes"]], manifest[["datasetFingerprint"]][["cellCount"]])
  }
  chunks <- manifest[["chunks"]]
  if (!.is_json_array(chunks)) .session_manifest_error("chunks must be an array")
  seen <- character()
  lazy_seen <- FALSE
  for (position in seq_along(chunks)) {
    raw <- chunks[[position]]
    if (!.is_json_object(raw)) .session_manifest_error("chunk #", position - 1L, " must be an object")
    if (!setequal(names(raw), .session_chunk_keys) || anyDuplicated(names(raw))) {
      .session_manifest_error("chunk #", position - 1L, " fields: ", .session_key_report(names(raw), .session_chunk_keys))
    }
    id <- raw[["id"]]
    if (!.is_trimmed_text(id)) .session_manifest_error("chunk #", position - 1L, " has an invalid id")
    if (id %in% seen) .session_manifest_error("duplicate chunk id '", id, "'")
    seen <- c(seen, id)
    if (!.is_trimmed_text(raw[["contributorId"]])) .session_manifest_error("chunk '", id, "' contributorId")
    priority <- raw[["priority"]]
    if (!.is_json_string(priority) || !priority %in% c("eager", "lazy")) .session_manifest_error("chunk '", id, "' priority")
    if (identical(priority, "lazy")) {
      lazy_seen <- TRUE
    } else if (lazy_seen) {
      .session_manifest_error("all eager chunks must precede lazy chunks")
    }
    if (!.is_json_string(raw[["kind"]]) || !raw[["kind"]] %in% c("json", "binary")) .session_manifest_error("chunk '", id, "' kind")
    if (!.is_json_string(raw[["codec"]]) || !raw[["codec"]] %in% c("none", "gzip")) .session_manifest_error("chunk '", id, "' codec")
    if (!.is_json_string(raw[["label"]]) || !nzchar(raw[["label"]])) .session_manifest_error("chunk '", id, "' label")
    flag <- raw[["datasetDependent"]]
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) .session_manifest_error("chunk '", id, "' datasetDependent")
    stored <- .require_session_integer(raw[["storedBytes"]], paste0("chunk '", id, "' storedBytes"), .session_max_frame_bytes)
    plain <- .require_session_integer(raw[["uncompressedBytes"]], paste0("chunk '", id, "' uncompressedBytes"), .session_max_safe_integer)
    if (identical(raw[["codec"]], "none") && stored != plain) {
      .session_manifest_error("chunk '", id, "' storedBytes and uncompressedBytes must match for codec 'none'")
    }
  }
  manifest
}

# Open one .cellucid-session file: the manifest, checked in full, and where
# each chunk's bytes are. The chunks themselves are read by
# .session_chunk() when used.
.read_session_bundle <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .abort("A session is the path of a .cellucid-session file.")
  }
  path <- path.expand(path)
  if (!file.exists(path) || dir.exists(path)) .abort("No session file at ", path, ".")
  size <- file.size(path)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", length(.session_magic))
  if (!identical(magic, .session_magic)) .abort(.session_signature_refusal(magic))
  manifest_length <- .read_u32(con)
  if (manifest_length <= 0) .abort("Invalid session manifest length: ", manifest_length)
  if (length(.session_magic) + 4 + manifest_length > size) {
    .abort("Invalid session manifest length: ", format(manifest_length, scientific = FALSE),
           " ", .plural(manifest_length, "byte"), " is more than the session file holds")
  }
  manifest_bytes <- readBin(con, "raw", manifest_length)
  if (length(manifest_bytes) != manifest_length) .abort("Unexpected EOF while reading session bundle")
  manifest <- .parse_session_json(manifest_bytes, "Invalid session manifest JSON")
  manifest <- .validate_session_manifest(manifest)
  offsets <- numeric()
  for (meta in manifest[["chunks"]]) {
    stored <- .read_u32(con)
    if (stored != meta[["storedBytes"]]) {
      .abort("Invalid session chunk '", meta[["id"]], "': storedBytes declares ", meta[["storedBytes"]],
             ", frame contains ", stored)
    }
    offset <- seek(con)
    if (offset + stored > size) .abort("Invalid session chunk '", meta[["id"]], "': frame exceeds file size")
    seek(con, offset + stored)
    offsets <- c(offsets, offset)
  }
  if (seek(con) != size) .abort("Invalid session bundle (trailing bytes after final chunk)")
  ids <- vapply(manifest[["chunks"]], function(meta) meta[["id"]], character(1L))
  structure(
    list(
      path = normalizePath(path, winslash = "/"), manifest = manifest,
      chunks = stats::setNames(manifest[["chunks"]], ids), offsets = stats::setNames(offsets, ids)
    ),
    class = "cellucid_session_bundle"
  )
}

.parse_session_json <- function(bytes, what) {
  text <- .utf8_text(bytes)
  if (is.na(text)) .abort(what, ": not UTF-8 text")
  tryCatch(
    jsonlite::parse_json(text, simplifyVector = FALSE),
    error = function(error) .abort(what, ": ", conditionMessage(error))
  )
}

# One chunk's payload: inflated to exactly its declared size, then parsed
# (a json chunk) or returned as bytes (a binary chunk).
.session_chunk <- function(bundle, id) {
  meta <- bundle$chunks[[id]]
  if (is.null(meta)) .abort("The session has no chunk '", id, "'.")
  con <- file(bundle$path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, bundle$offsets[[id]])
  stored <- readBin(con, "raw", meta[["storedBytes"]])
  if (length(stored) != meta[["storedBytes"]]) .abort("Unexpected EOF while reading session bundle")
  payload <- if (identical(meta[["codec"]], "gzip")) {
    .native_inflate(stored, TRUE, meta[["uncompressedBytes"]])
  } else {
    stored
  }
  if (identical(meta[["kind"]], "binary")) payload else .parse_session_json(payload, "Invalid JSON chunk")
}
