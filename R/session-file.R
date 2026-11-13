# Reading a .cellucid-session file (session format 4):
# the history of what was done in the viewer, the file the viewer's Save
# session writes and cellucid_session_result() collects. The web app's
# app/session/README.md (what the manifest and chunks say), bundle/format.js
# and bundle/reader.js (the framing), codecs/cell-set.js and
# app/history/value-parts.js / value-pool.js (how values are stored and
# hashed), and app/history/model.js (every record, R/history-rules.R), check
# for check: a file the viewer refuses is refused here, and a file a newer
# Cellucid wrote, whose records this package does not know, is said to be
# newer, never damaged.
#
#   CELLUCID_SESSION/4\n          the signature
#   u32 LE, manifest JSON         createdAt, app, actions, scope,
#                                 datasetFingerprint, source, requires,
#                                 chrome (page chrome outside the history),
#                                 cells? (a subsample), genes? (a gene
#                                 sample), activeTracker, trackers (each the
#                                 model's Tracker, its steps, journal and
#                                 conversations as chunk ids), chunks
#   [u32 LE, chunk bytes] ...     one frame per chunk, in the manifest's
#                                 order: every eager one before any lazy one,
#                                 and within a priority the kinds in the
#                                 order below
#
# Chunks:
# - steps (steps/<n>/<k>, eager): a run of one tracker's Step records, in id
#   order; a tracker lists its steps chunks, read as one list;
# - parts (parts/<priority>/<k>): name -> skeleton, sub-values stored once;
# - array (array/<digest>): one typed array, stored once by its digest
#   <type>:<length>:<sha256 of its little-endian bytes>: its elements, or a
#   cell set for a strictly increasing Uint32Array;
# - values (values/<priority>/<k>): hash -> the value's skeleton; eager: the
#   values the active head shows (with the parts and arrays they name), lazy:
#   every other one; at least one eager chunk, none in a file shared as its
#   steps only (scope.content "steps");
# - journal (journal/<n>, lazy): one tracker's JournalEntry records;
# - chat (chat/<n>/<m>, lazy): a conversation.
# steps, parts and values are JSON cut into pieces of about 16 MiB of text;
# each chunk is stored gzip'd or as it is (codec "none", storedBytes equal to
# uncompressedBytes). A skeleton is JSON with {"$array": digest},
# {"$number": "NaN"|"Infinity"|"-Infinity"|"-0"} for what JSON cannot hold, and
# {"$part": <name>} for an object or list stored once as a part (a part names
# the parts inside it). A part is storage only: a value is checked against its
# own hash once joined, never by a part's name. No picture of a step is ever
# stored.
#
# The file is untrusted input, and is read whole and checked before anything
# is returned: every chunk inflated to exactly its declared size, every array
# against its digest, every value against its hash, every record against the
# history model and every tree whole, what the share left out absent, the
# values the active head shows in the eager chunks with what they name, no
# part that names itself (through others) and no value whose text would pass
# 2^29 - 24 characters (both refused before anything is expanded), no value,
# part or array nothing refers to, none a step needs that the file lacks, no
# trailing bytes.

.session_format_version <- 4L
.session_signature_prefix <- "CELLUCID_SESSION/"
.session_magic <- charToRaw(paste0(.session_signature_prefix, .session_format_version, "\n"))
.not_a_session_file_message <- paste0(
  "This file is not a Cellucid session file of format ", .session_format_version,
  " (it does not start with CELLUCID_SESSION/", .session_format_version, ")."
)

# A file holding a value nested more deeply than R's own stacks follow: no
# Cellucid writes one (cellucid-python's TOO_DEEP_SESSION_MESSAGE, the same
# sentence).
.session_too_deep_message <- paste0(
  "This session file holds a value nested more deeply than this package reads, ",
  "so it was not opened. Cellucid does not write such files: was it made by another program?"
)

# The web app's sentence for a session file of another format version
# (bundle/format.js otherSessionFormatMessage): saved by an earlier Cellucid
# (one state, where this version keeps the history of every step), or by a
# newer one.
.other_session_format_message <- function(version) {
  if (version < .session_format_version) {
    paste0(
      "This session was saved by an earlier version of Cellucid (session format ", version, "), ",
      "which this version no longer reads: a session file now holds the history of every step ",
      "(session format ", .session_format_version, "). It was not opened. Open it in the version ",
      "of Cellucid that saved it, or re-create the session in this version."
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
  digits <- if (is.na(rest)) character() else regmatches(rest, regexpr("^[0-9]+(?=\n|$)", rest, perl = TRUE))
  version <- if (length(digits) == 1L) as.numeric(digits) else NA_real_
  if (is.na(version) || version == 0 || version == .session_format_version) {
    .not_a_session_file_message
  } else {
    .other_session_format_message(version)
  }
}

# A frame's length is a u32: the largest manifest or stored chunk the format
# can frame.
.session_max_frame_bytes <- 4294967295
# The largest integer the web app reads exactly from JSON (2^53 - 1).
.session_max_safe_integer <- 9007199254740991

# A file that breaks the format, in the web app's words: one sentence, ending
# with its full stop (as cellucid-python says it).
.session_damaged <- function(...) {
  text <- paste0(...)
  .abort("This session file is damaged: ", text, if (!grepl("[.?]$", text)) ".")
}

.session_manifest_error <- function(...) .session_damaged("its manifest ", ...)

.is_safe_integer <- function(value) {
  is.numeric(value) && length(value) == 1L && !is.na(value) && is.finite(value) &&
    value == floor(value) && abs(value) <= .session_max_safe_integer
}

.is_finite_number <- function(value) is.numeric(value) && length(value) == 1L && !is.na(value) && is.finite(value)

# .is_text (agent-link.R): one nonempty string.
.require_session_text <- function(value, label) {
  if (!.is_text(value)) .session_damaged(label, " must be a non-empty text.")
  value
}

.require_nullable_text <- function(value, label) {
  if (!is.null(value)) .require_session_text(value, label)
  value
}

.session_key_report <- function(present, expected) {
  missing <- sort(setdiff(expected, present), method = "radix")
  unknown <- sort(setdiff(present, expected), method = "radix")
  paste(c(
    if (length(missing) > 0L) paste0("missing ", paste(missing, collapse = ", ")),
    if (length(unknown) > 0L) paste0("unknown ", paste(unknown, collapse = ", "))
  ), collapse = "; ")
}

.require_session_keys <- function(value, expected, label) {
  if (!.is_json_object(value)) .session_damaged(label, " must be a record.")
  if (anyDuplicated(names(value)) || !setequal(names(value), expected) || length(value) != length(expected)) {
    .session_damaged(label, " must hold exactly ", paste(expected, collapse = ", "), " (",
                     .session_key_report(names(value), expected), ").")
  }
  value
}

.is_trimmed_text <- function(value) {
  # Python's str.strip(): Unicode white space, and the separators \x1c-\x1f.
  .is_json_string(value) && nzchar(value) &&
    !grepl("(*UCP)^[\\s\\x{1c}-\\x{1f}]|[\\s\\x{1c}-\\x{1f}]$", value, perl = TRUE)
}

# A JSON value of finite numbers, objects with keys named once, and lists (a
# manifest's dataset record, a conversation). jsonlite keeps a repeated key,
# where JavaScript keeps its last value: a file that repeats one is refused.
.require_session_json <- function(value, label) {
  if (is.null(value) || is.character(value) || is.logical(value)) return(invisible(value))
  if (is.numeric(value)) {
    if (!all(is.finite(value))) .session_damaged(label, " must be JSON.")
    return(invisible(value))
  }
  if (!is.list(value)) .session_damaged(label, " must be JSON.")
  if (!is.null(names(value)) && anyDuplicated(names(value))) .session_damaged(label, " names a key twice.")
  for (item in value) .require_session_json(item, label)
  invisible(value)
}

# --- The dataset part of the manifest (app/session/session-context.js) ---------

.require_session_integer <- function(value, label, maximum) {
  if (!.is_safe_integer(value)) .session_manifest_error(label, " must be an integer")
  if (value < 0) .session_manifest_error(label, " must be non-negative")
  if (value > maximum) .session_manifest_error(label, " exceeds ", format(maximum, scientific = FALSE))
  as.numeric(value)
}

# The fingerprint that ties a session to one dataset and data mapping;
# returns it with its mapping in canonical form.
.validate_session_fingerprint <- function(value) {
  if (!.is_json_object(value)) .session_manifest_error("datasetFingerprint must be an object")
  expected <- c("datasetId", "name", "cellCount", "varCount", "cellOrder", "dataMapping")
  if (!setequal(names(value), expected) || anyDuplicated(names(value))) {
    .session_manifest_error("datasetFingerprint fields: ", .session_key_report(names(value), expected))
  }
  for (key in c("datasetId", "name")) {
    item <- value[[key]]
    if (!is.null(item) && !.is_trimmed_text(item)) {
      .session_manifest_error("datasetFingerprint.", key, " must be a non-empty trimmed string or null")
    }
  }
  for (key in c("cellCount", "varCount")) {
    .require_session_integer(value[[key]], paste0("datasetFingerprint.", key), .session_max_safe_integer)
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
    if (!.is_text(value[[key]])) .session_manifest_error("source.", key, " must be a non-empty text")
  }
  if (identical(kind, "remote") && !grepl("^https?://[^/?#]", value[["url"]])) {
    .session_manifest_error("source.url must be an http(s) address")
  }
  value
}

# The dataset fields the shown state uses: obs columns and genes, each list
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
      if (!.is_text(keys[[index]])) {
        .session_manifest_error("requires.", source, "[", index - 1L, "] must be a non-empty string")
      }
    }
    keys <- vapply(keys, identity, character(1L))
    if (length(keys) > 1L && (!identical(.utf16_order(keys), seq_along(keys)) || anyDuplicated(keys))) {
      .session_manifest_error("requires.", source, " must be sorted and without duplicates")
    }
    out[[source]] <- keys
  }
  out
}

# The page chrome a session carries outside its history (history-file.js
# assertPageChrome): list(pressed = list(<group id> = <button id>)), the
# pressed button of each of the web app's pressed-button groups (the selection
# tool), which no step records; empty from a file written without a page.
.is_stable_id <- function(value) {
  .is_text(value) && nzchar(value) && identical(value, trimws(value, whitespace = "[\\h\\v]"))
}

.validate_session_chrome <- function(chrome) {
  if (!.is_json_object(chrome) || !identical(names(chrome), "pressed")) .session_damaged("its page chrome must be exactly {pressed}.")
  pressed <- chrome$pressed
  if (!.is_json_object(pressed)) .session_damaged("its page chrome must map pressed-button groups to their pressed button.")
  for (at in seq_along(pressed)) {
    group <- names(pressed)[[at]]
    # The web app's rule for an id it puts back (utils/stable-id.js): nonempty, no white space around it.
    if (!.is_stable_id(group) || !.is_stable_id(pressed[[at]])) {
      .session_damaged("its page chrome must name the pressed button of the group \"", group, "\" by its id.")
    }
  }
  invisible(chrome)
}

# createdAt: the canonical ISO time JavaScript's toISOString() writes.
.is_canonical_time <- function(value) {
  if (!.is_json_string(value) || !grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}\\.\\d{3}Z$", value, perl = TRUE)) {
    return(FALSE)
  }
  parsed <- as.POSIXct(substring(value, 1L, 19L), format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
  !is.na(parsed) && identical(format(parsed, "%Y-%m-%dT%H:%M:%S", tz = "UTC"), substring(value, 1L, 19L))
}

# The dataset part of a manifest, checked; returns the fingerprint with its
# mapping in canonical form (the manifest itself is kept as the file says).
.validate_session_header <- function(manifest) {
  if (!.is_canonical_time(manifest[["createdAt"]])) .session_manifest_error("createdAt must be a canonical ISO time.")
  fingerprint <- .validate_session_fingerprint(manifest[["datasetFingerprint"]])
  .validate_session_source(manifest[["source"]])
  .validate_session_requires(manifest[["requires"]])
  .validate_session_chrome(manifest[["chrome"]])
  if ("cells" %in% names(manifest)) {
    .session_cells(manifest[["cells"]], manifest[["datasetFingerprint"]][["cellCount"]])
  }
  if ("genes" %in% names(manifest)) {
    .session_genes(manifest[["genes"]], manifest[["datasetFingerprint"]][["cellCount"]])
  }
  fingerprint
}

# --- Newer files --------------------------------------------------------------------

# A Cellucid version's numbers (0.9.7, 0.9.7.dev1 -> c(0, 9, 7)), or NULL.
.session_version_numbers <- function(text) {
  if (!.is_json_string(text)) return(NULL)
  parts <- regmatches(text, regexec("^(0|[1-9][0-9]*)\\.(0|[1-9][0-9]*)\\.(0|[1-9][0-9]*)", text))[[1L]]
  if (length(parts) != 4L) NULL else as.numeric(parts[2:4])
}

.session_package_version <- function() as.character(utils::packageVersion("cellucid"))

# Why a session file written by a newer Cellucid than this package is not
# read: the same sentence in cellucid-python and cellucid-r. A file whose
# records this package does not know, written by a newer web app, is newer,
# never "damaged".
.newer_session_message <- function(version, package) {
  paste0("This session was written by a newer Cellucid (", version, ") than this package reads (", package,
         "); update the cellucid package to read it.")
}

# The writing version when the manifest names one newer than this package.
.session_newer_than_package <- function(manifest) {
  version <- if (.is_json_object(manifest) && .is_json_object(manifest$app)) manifest$app$version else NULL
  written <- .session_version_numbers(version)
  current <- .session_version_numbers(.session_package_version())
  if (is.null(written) || is.null(current)) return(NULL)
  for (at in 1:3) {
    if (written[[at]] != current[[at]]) return(if (written[[at]] > current[[at]]) version else NULL)
  }
  NULL
}

# --- The history part of the manifest (app/session/README.md) ----------------------

.session_chunk_base_keys <- c("id", "kind", "priority", "codec", "storedBytes", "uncompressedBytes")
.session_chunk_keys <- list(
  steps = c(.session_chunk_base_keys, "tracker"),
  parts = .session_chunk_base_keys,
  array = c(.session_chunk_base_keys, "digest", "encoding"),
  values = .session_chunk_base_keys,
  journal = c(.session_chunk_base_keys, "tracker"),
  chat = c(.session_chunk_base_keys, "tracker", "conversation")
)
# The kinds in the order a file holds them within a priority (history-file.js KIND_ORDER).
.session_kind_order <- c("steps", "parts", "array", "values", "journal", "chat")
# The longest skeleton text, in UTF-16 code units, one value may expand to
# (value-parts.js MAX_VALUE_TEXT: V8's longest string).
.session_max_value_text <- 2^29 - 24
# A manifest tracker: the model's Tracker, its steps, journal and
# conversations as chunk ids.
.session_tracker_keys <- vapply(.history_fields$Tracker, function(field) field$key, "")

# The typed arrays a value may hold, by the name their digest starts with,
# and the bytes of one element (app/history/value-parts.js ARRAY_TYPES).
.session_array_types <- c(
  Int8Array = 1L, Uint8Array = 1L, Uint8ClampedArray = 1L, Int16Array = 2L, Uint16Array = 2L, Int32Array = 4L,
  Uint32Array = 4L, Float32Array = 4L, Float64Array = 8L, BigInt64Array = 8L, BigUint64Array = 8L, Float16Array = 2L
)

# What an array digest says: its type, length and the SHA-256 of its bytes.
.parse_array_digest <- function(digest) {
  parts <- if (.is_json_string(digest)) {
    regmatches(digest, regexec("^([A-Za-z0-9]+):(0|[1-9][0-9]*):([0-9a-f]{64})$", digest))[[1L]]
  } else {
    character()
  }
  if (length(parts) != 4L || !parts[[2L]] %in% names(.session_array_types) ||
      as.numeric(parts[[3L]]) > .session_max_safe_integer) {
    .session_damaged("\"", substr(as.character(digest)[1L], 1L, 80L), "\" is not the digest of an array a state value may hold.")
  }
  list(type = parts[[2L]], length = as.numeric(parts[[3L]]), sha256 = parts[[4L]])
}

.validate_tracker_record <- function(record, index) {
  label <- paste0("tracker ", index)
  .require_session_keys(record, .session_tracker_keys, paste0("its manifest ", label))
  .require_session_text(record$id, paste(label, "id"))
  if (!.is_json_array(record$steps) || length(record$steps) == 0L) .session_damaged(label, " steps must list its steps chunks.")
  for (at in seq_along(record$steps)) .require_session_text(record$steps[[at]], paste0(label, " steps[", at - 1L, "]"))
  if (!is.null(record$journal)) .require_session_text(record$journal, paste(label, "journal"))
  if (!.is_json_array(record$chats)) .session_damaged(label, " chats must be a list of chunk ids.")
  for (at in seq_along(record$chats)) .require_session_text(record$chats[[at]], paste0(label, " chats[", at - 1L, "]"))
}

.validate_chunk_record <- function(entry, index) {
  label <- paste0("chunk ", index)
  kind <- if (.is_json_object(entry)) entry$kind else NULL
  if (!.is_json_string(kind) || is.null(.session_chunk_keys[[kind]])) {
    .session_damaged(label, " kind must be one of ", paste(.session_kind_order, collapse = ", "), ".")
  }
  .require_session_keys(entry, .session_chunk_keys[[kind]], paste0("its manifest ", label))
  .require_session_text(entry$id, paste(label, "id"))
  if (!.is_json_string(entry$priority) || !entry$priority %in% c("eager", "lazy")) {
    .session_damaged(label, " priority must be eager or lazy.")
  }
  if (!.is_json_string(entry$codec) || !entry$codec %in% c("gzip", "none")) .session_damaged(label, " codec must be gzip or none.")
  if (!.is_safe_integer(entry$storedBytes) || entry$storedBytes < 0 || entry$storedBytes > .session_max_frame_bytes) {
    .session_damaged(label, " storedBytes must be a frame length.")
  }
  if (!.is_safe_integer(entry$uncompressedBytes) || entry$uncompressedBytes < 0) {
    .session_damaged(label, " uncompressedBytes must be a byte count.")
  }
  if (identical(entry$codec, "none") && entry$storedBytes != entry$uncompressedBytes) {
    .session_damaged(label, " is stored as it is, at another size.")
  }
  if (kind %in% c("steps", "chat", "journal")) .require_session_text(entry$tracker, paste(label, "tracker"))
  if (identical(kind, "steps") && !identical(entry$priority, "eager")) .session_damaged(label, ": steps are read first.")
  if (kind %in% c("chat", "journal") && !identical(entry$priority, "lazy")) .session_damaged(label, ": a ", kind, " is read last.")
  if (identical(kind, "chat")) .require_session_text(entry$conversation, paste(label, "conversation"))
  if (identical(kind, "array")) {
    digest <- .parse_array_digest(entry$digest)
    if (!.is_json_string(entry$encoding) || !entry$encoding %in% c("raw", "cell-set") ||
        (identical(entry$encoding, "cell-set") && !identical(digest$type, "Uint32Array"))) {
      .session_damaged(label, " encoding must be raw, or cell-set for a Uint32Array.")
    }
  }
}

.session_model_record <- function(name, value, label) {
  tryCatch(.history_check_record(name, value, label),
           cellucid_history_model_error = function(error) .session_damaged(conditionMessage(error), "."))
  value
}

# A whole format-4 manifest, checked before any chunk is read; returns the
# dataset's fingerprint, its mapping in canonical form.
.validate_history_manifest <- function(manifest) {
  if (!.is_json_object(manifest)) .session_damaged("its manifest must be a record.")
  header <- c("createdAt", "app", "actions", "scope", "datasetFingerprint", "source", "requires", "chrome",
              if ("cells" %in% names(manifest)) "cells", if ("genes" %in% names(manifest)) "genes")
  .require_session_keys(manifest, c(header, "activeTracker", "trackers", "chunks"), "its manifest")
  .session_model_record("AppStamp", manifest$app, "session.app")
  if (!.is_json_array(manifest$actions)) .session_damaged("session.actions must be a list.")
  for (at in seq_along(manifest$actions)) .session_model_record("ActionStamp", manifest$actions[[at]], paste0("session.actions[", at - 1L, "]"))
  scope <- .session_model_record("Scope", manifest$scope, "session.scope")
  fingerprint <- .validate_session_header(manifest)
  if (identical(scope$source, "name") && !is.null(manifest$source)) {
    .session_damaged("it was shared with the dataset's name only (scope.source), and names where it opens from.")
  }
  if (!.is_json_array(manifest$trackers) || length(manifest$trackers) == 0L) .session_damaged("it holds no tracker.")
  if (!.is_json_array(manifest$chunks)) .session_damaged("its manifest has no chunk list.")
  for (at in seq_along(manifest$trackers)) .validate_tracker_record(manifest$trackers[[at]], at - 1L)
  for (at in seq_along(manifest$chunks)) .validate_chunk_record(manifest$chunks[[at]], at - 1L)
  tracker_ids <- vapply(manifest$trackers, function(tracker) tracker$id, "")
  if (anyDuplicated(tracker_ids)) .session_damaged("it holds the tracker \"", tracker_ids[anyDuplicated(tracker_ids)], "\" twice.")
  for (tracker in manifest$trackers) {
    if (!is.null(tracker$journal) && !isTRUE(scope$journal)) {
      .session_damaged("tracker \"", tracker$id, "\" keeps its journal, which the share left out (scope.journal).")
    }
    if (length(tracker$chats) > 0L && !isTRUE(scope$chats)) {
      .session_damaged("tracker \"", tracker$id, "\" keeps conversations, which the share left out (scope.chats).")
    }
  }
  if (!.is_json_string(manifest$activeTracker) || !manifest$activeTracker %in% tracker_ids) {
    .session_damaged("its active tracker is not one it holds.")
  }
  chunk_ids <- vapply(manifest$chunks, function(entry) entry$id, "")
  if (anyDuplicated(chunk_ids)) .session_damaged("it lists the chunk \"", chunk_ids[anyDuplicated(chunk_ids)], "\" twice.")
  priorities <- vapply(manifest$chunks, function(entry) entry$priority, "")
  kinds <- vapply(manifest$chunks, function(entry) entry$kind, "")
  # Eager before lazy, and within a priority the kinds in their order.
  rank <- (priorities == "lazy") * length(.session_kind_order) + match(kinds, .session_kind_order)
  back <- which(diff(rank) < 0)
  if (length(back) > 0L) {
    .session_damaged("its chunk \"", chunk_ids[[back[[1L]] + 1L]], "\" is out of order (eager before lazy; ",
                     paste(.session_kind_order, collapse = ", "), ").")
  }
  digests <- vapply(manifest$chunks[kinds == "array"], function(entry) entry$digest, "")
  if (anyDuplicated(digests)) .session_damaged("it stores one array twice.")
  if (identical(scope$content, "steps")) {
    if (any(kinds == "values")) .session_damaged("it was shared as its steps only (scope.content), and holds values.")
  } else if (!any(kinds == "values" & priorities == "eager")) {
    .session_damaged("it holds no values of the step it shows.")
  }
  chunk_of <- function(id) {
    at <- match(id, chunk_ids)
    if (is.na(at)) NULL else manifest$chunks[[at]]
  }
  named <- character()
  for (tracker in manifest$trackers) {
    for (id in unlist(tracker$steps)) {
      steps <- chunk_of(id)
      if (!identical(steps$kind, "steps") || !identical(steps$tracker, tracker$id) || id %in% named) {
        .session_damaged("tracker \"", tracker$id, "\" names a steps chunk that is not its own.")
      }
      named <- c(named, id)
    }
    # A tracker's steps are read as one list, in the order it names its chunks: the file holds them in that order.
    if (is.unsorted(match(unlist(tracker$steps), chunk_ids))) {
      .session_damaged("tracker \"", tracker$id, "\" names its steps chunks in another order than the file holds them.")
    }
    if (!is.null(tracker$journal)) {
      journal <- chunk_of(tracker$journal)
      if (!identical(journal$kind, "journal") || !identical(journal$tracker, tracker$id) || tracker$journal %in% named) {
        .session_damaged("tracker \"", tracker$id, "\" names a journal chunk that is not its own.")
      }
      named <- c(named, tracker$journal)
    }
    conversations <- character()
    for (id in unlist(tracker$chats)) {
      chat <- chunk_of(id)
      if (!identical(chat$kind, "chat") || !identical(chat$tracker, tracker$id) || chat$conversation %in% conversations ||
          id %in% named) {
        .session_damaged("tracker \"", tracker$id, "\" names a conversation chunk that is not its own.")
      }
      conversations <- c(conversations, chat$conversation)
      named <- c(named, id)
    }
  }
  for (id in chunk_ids[kinds %in% c("steps", "chat", "journal")]) {
    if (!id %in% named) .session_damaged("its chunk \"", id, "\" belongs to no tracker.")
  }
  fingerprint
}

# --- Arrays and values -----------------------------------------------------------

# A canonical uvarint at the start of `bytes` (codecs/varint.js decodeUvarint):
# its value and the bytes it took.
.read_uvarint <- function(bytes, label) {
  value <- 0
  for (at in seq_len(min(10L, length(bytes)))) {
    byte <- as.integer(bytes[[at]])
    group <- bitwAnd(byte, 127L)
    value <- value + group * 2^(7 * (at - 1L))
    if (value > .session_max_safe_integer) .session_damaged(label, ": a varint exceeds the exact integer range.")
    if (bitwAnd(byte, 128L) == 0L) {
      if (at > 1L && group == 0L) .session_damaged(label, ": a varint is not canonical.")
      return(list(value = value, length = at))
    }
  }
  .session_damaged(label, ": a varint is truncated or too long.")
}

# The values of `count` little-endian elements of a typed array, as R holds
# them: integers for the 8-, 16- and signed 32-bit types, doubles for the
# others (a Uint32Array's values pass R's integers; a 64-bit integer is exact
# up to 2^53).
.decode_typed_values <- function(bytes, type, count) {
  if (count == 0) return(if (type %in% c("Int8Array", "Uint8Array", "Uint8ClampedArray", "Int16Array", "Uint16Array", "Int32Array")) integer(0) else double(0))
  switch(type,
    Int8Array = readBin(bytes, "integer", n = count, size = 1L, signed = TRUE),
    Uint8Array = , Uint8ClampedArray = readBin(bytes, "integer", n = count, size = 1L, signed = FALSE),
    Int16Array = readBin(bytes, "integer", n = count, size = 2L, signed = TRUE, endian = "little"),
    Uint16Array = readBin(bytes, "integer", n = count, size = 2L, signed = FALSE, endian = "little"),
    Int32Array = readBin(bytes, "integer", n = count, size = 4L, endian = "little"),
    Uint32Array = .decode_wire_values(bytes, "uint32", count),
    Float32Array = readBin(bytes, "double", n = count, size = 4L, endian = "little"),
    Float64Array = readBin(bytes, "double", n = count, size = 8L, endian = "little"),
    Float16Array = .decode_float16(readBin(bytes, "integer", n = count, size = 2L, signed = FALSE, endian = "little")),
    BigInt64Array = , BigUint64Array = {
      halves <- .decode_wire_values(bytes, "uint32", 2 * count)
      low <- halves[c(TRUE, FALSE)]
      high <- halves[c(FALSE, TRUE)]
      if (identical(type, "BigInt64Array")) high[high >= 2^31] <- high[high >= 2^31] - 2^32
      high * 2^32 + low
    }
  )
}

.decode_float16 <- function(bits) {
  sign <- ifelse(bitwAnd(bits, 0x8000L) != 0L, -1, 1)
  exponent <- bitwAnd(bitwShiftR(bits, 10L), 0x1FL)
  fraction <- bitwAnd(bits, 0x3FFL)
  value <- ifelse(exponent == 0L, fraction * 2^-24,
                  ifelse(exponent == 31L, ifelse(fraction == 0L, Inf, NaN), (1 + fraction / 1024) * 2^(exponent - 15)))
  sign * value
}

# One array chunk: its elements (raw) or its cell set, checked against its
# digest (type, length, SHA-256 of the little-endian bytes). The array as a
# state value holds it: an R vector whose attribute "cellucid_array" names its
# JavaScript type.
.decode_session_array <- function(entry, bytes) {
  digest <- .parse_array_digest(entry$digest)
  size <- .session_array_types[[digest$type]]
  if (identical(entry$encoding, "cell-set")) {
    bound <- .read_uvarint(bytes, paste0("the array chunk \"", entry$id, "\""))
    cells <- tryCatch(
      .native_decode_cell_set(bytes[-seq_len(bound$length)], bound$value),
      error = function(error) .session_damaged("the array chunk \"", entry$id, "\": ", conditionMessage(error))
    )
    if (length(cells) != digest$length) .session_damaged("the array chunk \"", entry$id, "\" does not match its digest.")
    raw_bytes <- .encode_uint32(cells)
    values <- cells
  } else {
    if (length(bytes) != digest$length * size) .session_damaged("the array chunk \"", entry$id, "\" has the wrong length.")
    raw_bytes <- bytes
    values <- .decode_typed_values(bytes, digest$type, digest$length)
  }
  if (!identical(.native_sha256(raw_bytes), digest$sha256)) {
    .session_damaged("the array chunk \"", entry$id, "\" does not match its digest.")
  }
  attr(values, "cellucid_array") <- digest$type
  values
}

# The numbers JSON cannot hold, by their name in a skeleton.
.session_special_numbers <- c(`NaN` = NaN, Infinity = Inf, `-Infinity` = -Inf, `-0` = -0)

# The placeholders of a skeleton (value-parts.js).
.session_placeholder_keys <- c("$array", "$number", "$part")

# JSON texts' lengths as JavaScript counts them (UTF-16 code units: a
# character past U+FFFF is two).
.utf16_lengths <- function(texts) {
  out <- nchar(texts, type = "chars")
  wide <- which(nchar(texts, type = "bytes") > out)
  for (at in wide) out[[at]] <- out[[at]] + sum(utf8ToInt(texts[[at]]) > 0xFFFF)
  out
}

# The lengths of a list's scalar items written as JSON.stringify writes them
# (each a length-one vector, as jsonlite reads JSON); a non-finite number is
# refused.
.session_scalar_lengths <- function(items, label) {
  out <- numeric(length(items))
  kinds <- vapply(items, typeof, "", USE.NAMES = FALSE)
  for (kind in unique(kinds)) {
    at <- which(kinds == kind)
    values <- unlist(items[at], use.names = FALSE)
    out[at] <- switch(kind,
      logical = ifelse(values, 4, 5),
      character = .utf16_lengths(.js_strings(values)),
      integer = , double = {
        if (!all(is.finite(values))) .session_damaged(label, ": a stored value has a non-finite number.")
        nchar(.js_numbers(values), type = "chars")
      },
      .session_damaged(label, ": a stored value is not JSON.")
    )
  }
  out
}

# The value a skeleton describes (value-parts.js joinValue), its arrays from
# `arrays` (digest -> array) and its parts from `part_of(name)` (joined once).
# The skeleton was checked by the reader's walk before it is joined.
.join_value <- function(node, arrays, part_of) {
  if (!is.list(node)) return(node)
  keys <- names(node)
  if (is.null(keys)) return(lapply(node, .join_value, arrays = arrays, part_of = part_of))
  if (identical(keys, "$array")) return(get(node[[1L]], envir = arrays))
  if (identical(keys, "$number")) return(.session_special_numbers[[node[[1L]]]])
  if (identical(keys, "$part")) return(part_of(node[[1L]]))
  out <- lapply(node, .join_value, arrays = arrays, part_of = part_of)
  if (length(out) == 0L) out <- structure(list(), names = character())
  out
}

# --- Reading ----------------------------------------------------------------------

.parse_session_json <- function(bytes, what) {
  text <- .utf8_text(bytes)
  if (is.na(text)) .session_damaged(what, " is not JSON.")
  tryCatch(
    jsonlite::parse_json(text, simplifyVector = FALSE),
    error = function(error) .session_damaged(what, " is not JSON.")
  )
}

.read_le_u32 <- function(bytes, at) {
  sum(as.numeric(bytes[at + 0:3]) * c(1, 256, 65536, 16777216))
}

# Every time field of a record and the records in it.
.session_walk_times <- function(name, value, visit, path) {
  if (is.null(value)) return(invisible())
  for (field in .history_fields[[name]]) {
    if (!field$key %in% names(value)) next
    item <- value[[field$key]]
    spec <- field$spec
    where <- paste0(path, ".", field$key)
    if (identical(spec$type, "time")) {
      visit(item, where)
    } else if (identical(spec$type, "record") && !is.null(item) && !spec$name %in% .history_map_records) {
      .session_walk_times(spec$name, item, visit, where)
    } else if (identical(spec$type, "array") && identical(spec$item$type, "record") && !is.null(item)) {
      for (at in seq_along(item)) .session_walk_times(spec$item$name, item[[at]], visit, paste0(where, "[", at - 1L, "]"))
    }
  }
}

# What the share left out is absent; its times are as it chose them.
.session_check_share <- function(trackers, scope) {
  for (tracker in trackers) {
    label <- paste0("tracker \"", tracker$id, "\"")
    visit <- function(time, path) {
      if (identical(scope$times, "none") && !is.null(time)) {
        .session_damaged(label, " holds a time at ", path, ", and the share left times out (scope.times).")
      }
      # A note never edited has no edit time, whatever the share kept.
      if (!identical(scope$times, "none") && is.null(time) && !endsWith(path, ".editedAt")) {
        .session_damaged(label, " has no time at ", path, ", and the share kept times (scope.times).")
      }
      if (identical(scope$times, "dates") && !is.null(time) && time %% 86400000 != 0) {
        .session_damaged(label, " holds a time of day at ", path, ", and the share kept dates only (scope.times).")
      }
    }
    .session_walk_times("Tracker", tracker, visit, "tracker")
  }
}

cellucid_read_session <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .abort("A session is the path of a .cellucid-session file.")
  }
  file <- path.expand(path)
  if (!file.exists(file) || dir.exists(file)) .abort("No session file at ", path, ".")
  size <- file.size(file)
  bytes <- readBin(file, "raw", n = size)
  if (length(bytes) != size) .abort("Could not read the whole session file ", path, ".")

  # The framing (bundle/reader.js).
  magic_length <- length(.session_magic)
  if (size < magic_length) .abort(.not_a_session_file_message)
  if (!identical(bytes[seq_len(magic_length)], .session_magic)) .abort(.session_signature_refusal(bytes[seq_len(magic_length)]))
  at <- magic_length + 1
  if (size < at + 3) .abort("Unexpected end of file.")
  manifest_length <- .read_le_u32(bytes, at)
  at <- at + 4
  if (manifest_length <= 0) .abort("Invalid manifest length: 0 bytes.")
  if (manifest_length > size - at + 1) {
    .abort("Invalid manifest length: ", format(manifest_length, scientific = FALSE), " ", .plural(manifest_length, "byte"),
           " is more than the session file holds (session file truncated?).")
  }
  manifest_text <- .utf8_text(bytes[at:(at + manifest_length - 1)])
  manifest <- if (is.na(manifest_text)) NULL else tryCatch(jsonlite::parse_json(manifest_text, simplifyVector = FALSE), error = function(error) NULL)
  if (is.null(manifest)) .abort("Invalid session file (manifest JSON parse failed).")
  newer <- .session_newer_than_package(manifest)
  withCallingHandlers(
    # R's own limit on what the reader follows (its expression and C stacks,
    # R >= 4.4 names their overflow), refused in a sentence as every other
    # file it cannot read is (cellucid-python's TOO_DEEP_SESSION_MESSAGE).
    tryCatch(
      .read_session_body(file, bytes, size, at + manifest_length, manifest),
      stackOverflowError = function(error) .abort(.session_too_deep_message)
    ),
    error = function(error) {
      if (!is.null(newer)) .abort(.newer_session_message(newer, .session_package_version()))
    }
  )
}

.read_session_body <- function(file, bytes, size, at, manifest) {
  if (!.is_json_object(manifest)) .abort("Invalid session file (manifest must be a JSON object).")
  if (!.is_json_array(manifest$chunks)) .abort("Invalid session file (manifest.chunks must be an array).")
  fingerprint <- .validate_history_manifest(manifest)
  scope <- manifest$scope

  # Every chunk's frame, in the manifest's order, and nothing after them.
  payloads <- vector("list", length(manifest$chunks))
  for (index in seq_along(manifest$chunks)) {
    entry <- manifest$chunks[[index]]
    if (size - at + 1 < 4) .abort("Unexpected end of file.")
    stored <- .read_le_u32(bytes, at)
    at <- at + 4
    if (stored != entry$storedBytes) {
      .abort("Chunk length mismatch for ", entry$id, ": header=", format(stored, scientific = FALSE),
             ", manifest=", format(entry$storedBytes, scientific = FALSE))
    }
    if (stored > size - at + 1) {
      .abort("Invalid chunk length for ", entry$id, ": ", format(stored, scientific = FALSE), " > remaining ",
             format(size - at + 1, scientific = FALSE), " (session file truncated?)")
    }
    frame <- if (stored == 0) raw(0) else bytes[at:(at + stored - 1)]
    at <- at + stored
    payload <- if (identical(entry$codec, "gzip")) {
      tryCatch(.native_inflate(frame, TRUE, entry$uncompressedBytes), error = function(error) {
        .session_damaged("the chunk \"", entry$id, "\" holds another size than its manifest says (", conditionMessage(error), ").")
      })
    } else {
      frame
    }
    if (length(payload) != entry$uncompressedBytes) {
      .session_damaged("the chunk \"", entry$id, "\" holds another size than its manifest says.")
    }
    payloads[[index]] <- payload
  }
  if (at <= size) .abort("Invalid session file: trailing bytes after declared chunks.")
  sha256 <- .native_sha256(bytes)
  rm(bytes)

  # Every chunk read, each value, part and array once (environments: no copy
  # per value, and lookups by name that do not grow with the file).
  arrays <- new.env(parent = emptyenv())
  lazy_arrays <- new.env(parent = emptyenv())
  referenced_arrays <- new.env(parent = emptyenv())
  parts <- new.env(parent = emptyenv())
  lazy_parts <- new.env(parent = emptyenv())
  skeletons <- new.env(parent = emptyenv())
  eager_hashes <- character()
  lazy_hashes <- character()
  # A tracker's steps chunks are pieces of one list, in the file's order.
  step_pieces <- stats::setNames(lapply(manifest$trackers, function(tracker) list()), vapply(manifest$trackers, function(tracker) tracker$id, ""))
  journal_of <- list()
  chats <- list()
  take_map <- function(target, record, label, lazy) {
    if (!.is_json_object(record)) .session_damaged(label, " must map names to skeletons.")
    keys <- names(record)
    if (anyDuplicated(keys) || !all(grepl("^[0-9a-f]{64}$", keys))) .session_damaged(label, " holds a name twice or malformed.")
    for (at in seq_along(record)) {
      if (exists(keys[[at]], envir = target, inherits = FALSE)) .session_damaged(label, " holds a name twice or malformed.")
      assign(keys[[at]], record[at], envir = target)
      if (!is.null(lazy)) assign(keys[[at]], TRUE, envir = lazy)
    }
    keys
  }
  for (index in seq_along(manifest$chunks)) {
    entry <- manifest$chunks[[index]]
    payload <- payloads[[index]]
    lazy <- identical(entry$priority, "lazy")
    switch(entry$kind,
      steps = {
        label <- paste0("the steps of tracker \"", entry$tracker, "\"")
        piece <- .parse_session_json(payload, label)
        if (!.is_json_array(piece)) .session_damaged(label, " must be a list.")
        # Pieces of one list, in the file's order (the order the tracker names them).
        step_pieces[[entry$tracker]][[length(step_pieces[[entry$tracker]]) + 1L]] <- piece
      },
      parts = take_map(parts, .parse_session_json(payload, "its parts"), "its parts", if (lazy) lazy_parts),
      journal = journal_of[[entry$tracker]] <- .parse_session_json(payload, paste0("the journal of tracker \"", entry$tracker, "\"")),
      array = {
        assign(entry$digest, .decode_session_array(entry, payload), envir = arrays)
        if (lazy) assign(entry$digest, TRUE, envir = lazy_arrays)
      },
      values = {
        keys <- take_map(skeletons, .parse_session_json(payload, "its values"), "its values", NULL)
        if (lazy) lazy_hashes <- c(lazy_hashes, keys) else eager_hashes <- c(eager_hashes, keys)
      },
      chat = {
        label <- paste0("the conversation \"", entry$conversation, "\"")
        chats[[entry$id]] <- list(conversation = entry$conversation, transcript = .require_session_json(.parse_session_json(payload, label), label))
      }
    )
    payloads[index] <- list(NULL)
  }

  # Every skeleton is walked once before it is joined (value-parts.js
  # partJoiner): its placeholders checked, the arrays it names referred to,
  # and its text length counted with each part's counted once from its own
  # text (a part that names itself, through others, is refused); a value that
  # would pass .session_max_value_text is refused before it is expanded.
  # While the values the active head shows are walked, a part or an array
  # only the lazy chunks hold is refused: the web app shows the head from the
  # eager chunks alone.
  eager_only <- FALSE
  part_lengths <- new.env(parent = emptyenv())
  visiting <- new.env(parent = emptyenv())
  refer_array <- function(digest) {
    if (!nzchar(digest) || !exists(digest, envir = arrays, inherits = FALSE)) {
      .session_damaged("a value names an array the file does not hold (", substr(digest, 1L, 40L), "...).")
    }
    if (eager_only && exists(digest, envir = lazy_arrays, inherits = FALSE)) {
      .session_damaged("a value of the step it shows names an array that is not among its first chunks (", substr(digest, 1L, 40L), "...).")
    }
    assign(digest, TRUE, envir = referenced_arrays)
  }
  part_length <- function(name, label) {
    if (!exists(name, envir = parts, inherits = FALSE)) {
      .session_damaged("a value names a part the file does not hold (", substr(name, 1L, 12L), "...).")
    }
    known <- part_lengths[[name]]
    if (!is.null(known)) return(known)
    if (eager_only && exists(name, envir = lazy_parts, inherits = FALSE)) {
      .session_damaged("a value of the step it shows names a part that is not among its first chunks (", substr(name, 1L, 12L), "...).")
    }
    if (exists(name, envir = visiting, inherits = FALSE)) .session_damaged(label, ": a stored value names a part inside itself.")
    assign(name, TRUE, envir = visiting)
    size <- walk(get(name, envir = parts)[[1L]], label)
    rm(list = name, envir = visiting)
    assign(name, size, envir = part_lengths)
    size
  }
  walk <- function(node, label) {
    if (is.null(node)) return(4)
    if (!is.list(node)) return(.session_scalar_lengths(list(node), label))
    keys <- names(node)
    count <- length(node)
    size <- 2 + max(count - 1, 0)
    if (!is.null(keys)) {
      if (anyDuplicated(keys)) .session_damaged(label, ": a stored value names a key twice.")
      placeholder <- keys %in% .session_placeholder_keys
      if (any(placeholder)) {
        if (count != 1L) .session_damaged(label, ": a stored value has a malformed placeholder.")
        name <- node[[1L]]
        if (identical(keys, "$part")) {
          # The placeholder's shape is checked before the part is looked for.
          if (!.is_json_string(name) || !grepl("^[0-9a-f]{64}$", name)) .session_damaged(label, ": a stored value names a part wrongly.")
          return(part_length(name, label))
        }
        if (identical(keys, "$array")) {
          if (!.is_json_string(name)) .session_damaged(label, ": a stored value names an array wrongly.")
          refer_array(name)
        } else if (!.is_json_string(name) || !name %in% names(.session_special_numbers)) {
          .session_damaged(label, ": a stored value has an unknown number.")
        }
      }
      size <- size + sum(.utf16_lengths(.js_strings(keys))) + count
    }
    if (count == 0L) return(size)
    empty <- vapply(node, is.null, TRUE, USE.NAMES = FALSE)
    nested <- vapply(node, is.list, TRUE, USE.NAMES = FALSE)
    scalar <- !empty & !nested
    size <- size + 4 * sum(empty)
    if (any(scalar)) size <- size + sum(.session_scalar_lengths(node[scalar], label))
    for (child in node[nested]) size <- size + walk(child, label)
    size
  }
  measure <- function(skeleton, label) {
    size <- walk(skeleton, label)
    if (size > .session_max_value_text) {
      .session_damaged(label, ": a stored value expands to ", format(size, scientific = FALSE), " ", .plural(size, "character"),
                       " of text, more than a page can hold (", format(.session_max_value_text, scientific = FALSE), ").")
    }
    size
  }
  # A part's canonical text and its joined value, each made once.
  part_texts <- new.env(parent = emptyenv())
  part_text <- function(name) {
    known <- part_texts[[name]]
    if (!is.null(known)) return(known)
    text <- .canonical_text(get(name, envir = parts)[[1L]], part_text)
    assign(name, text, envir = part_texts)
    text
  }
  joined_parts <- new.env(parent = emptyenv())
  part_of <- function(name) {
    if (exists(name, envir = joined_parts, inherits = FALSE)) return(get(name, envir = joined_parts))
    value <- .join_value(get(name, envir = parts)[[1L]], arrays, part_of)
    assign(name, value, envir = joined_parts)
    value
  }
  check_values <- function(hashes, label) {
    for (hash in hashes) {
      skeleton <- get(hash, envir = skeletons)[[1L]]
      measure(skeleton, label)
      if (!identical(.value_hash(skeleton, part_text), hash)) .session_damaged(label, " holds a value that does not match its hash.")
    }
  }
  values_shared <- !identical(scope$content, "steps")
  eager_only <- TRUE
  check_values(eager_hashes, "its eager values")
  eager_only <- FALSE
  check_values(lazy_hashes, "its other values")

  joined_step <- function(step, label) {
    if (!.is_json_object(step)) .session_damaged(label, " must be a Step record.")
    for (key in intersect(c("params", "outcome"), names(step))) {
      measure(step[[key]], label)
      step[key] <- list(.join_value(step[[key]], arrays, part_of))
    }
    step
  }
  joined_patch <- function(patch, label) {
    if (!.is_json_object(patch) || !.is_json_array(patch$steps)) return(patch)
    patch$steps <- lapply(seq_along(patch$steps), function(at) joined_step(patch$steps[[at]], paste0(label, ".steps[", at - 1L, "]")))
    patch
  }
  trackers <- lapply(manifest$trackers, function(record) {
    stored <- do.call(c, step_pieces[[record$id]])
    record$steps <- lapply(seq_along(stored), function(at) joined_step(stored[[at]], paste0("tracker \"", record$id, "\" step ", at - 1L)))
    journal <- list()
    if (!is.null(record$journal)) {
      entries <- journal_of[[record$id]]
      if (!.is_json_array(entries)) .session_damaged("the journal of tracker \"", record$id, "\" must be a list.")
      journal <- lapply(seq_along(entries), function(at) {
        entry <- entries[[at]]
        label <- paste0("tracker \"", record$id, "\" journal entry ", at - 1L)
        if (!.is_json_object(entry)) .session_damaged(label, " must be a JournalEntry record.")
        for (part in intersect(c("before", "after"), names(entry))) entry[part] <- list(joined_patch(entry[[part]], paste0(label, ".", part)))
        entry
      })
    }
    record$journal <- journal
    record$chats <- lapply(unlist(record$chats), function(id) list(chats[[id]]$conversation, chats[[id]]$transcript))
    record
  })
  session_record <- list(
    createdAt = manifest$createdAt, app = manifest$app, actions = manifest$actions, scope = scope,
    datasetFingerprint = manifest$datasetFingerprint, source = manifest$source, requires = manifest$requires,
    cells = manifest$cells, genes = manifest$genes, activeTracker = manifest$activeTracker, trackers = trackers
  )
  # cells and genes are present (null when absent), as the model reads a session.
  session_record["cells"] <- list(manifest$cells)
  session_record["genes"] <- list(manifest$genes)
  session_record["source"] <- list(manifest$source)
  tryCatch(.history_validate_session(session_record),
           cellucid_history_model_error = function(error) .session_damaged(conditionMessage(error), "."))
  .session_check_share(trackers, scope)

  referenced <- new.env(parent = emptyenv())
  refer_step <- function(step) {
    for (hash in c(unlist(step$before), unlist(step$after))) assign(hash, TRUE, envir = referenced)
  }
  for (tracker in trackers) {
    for (step in tracker$steps) refer_step(step)
    for (entry in tracker$journal) {
      if (isTRUE(entry$purged)) next
      for (step in c(entry$before$steps, entry$after$steps)) refer_step(step)
    }
  }
  if (values_shared) {
    active <- Filter(function(tracker) identical(tracker$id, manifest$activeTracker), trackers)[[1L]]
    by_id <- new.env(parent = emptyenv(), size = length(active$steps))
    for (step in active$steps) assign(step$id, step, envir = by_id)
    # The values at the head: each domain's after-value nearest the head on its path from the root.
    head <- list()
    id <- active$headId
    while (!is.null(id)) {
      after <- by_id[[id]]$after
      head <- c(head, after[setdiff(names(after), names(head))])
      id <- by_id[[id]]$parent
    }
    if (!all(unlist(head) %in% eager_hashes)) .session_damaged("the values of the step it shows are not among its first chunks.")
    held <- ls(skeletons, all.names = TRUE)
    wanted <- ls(referenced, all.names = TRUE)
    if (!all(wanted %in% held)) .session_damaged("a step refers to a value the file does not hold.")
    if (!all(held %in% wanted)) .session_damaged("it holds a value no step refers to.")
  }
  if (!all(ls(parts, all.names = TRUE) %in% ls(part_lengths, all.names = TRUE))) .session_damaged("it holds a part nothing refers to.")
  if (!all(ls(arrays, all.names = TRUE) %in% ls(referenced_arrays, all.names = TRUE))) .session_damaged("it holds an array nothing refers to.")

  cells <- if (is.null(manifest$cells)) NULL else .session_cells(manifest$cells, fingerprint$cellCount)
  genes <- if (is.null(manifest$genes)) NULL else .session_genes(manifest$genes, fingerprint$cellCount)
  facts <- list(
    name = fingerprint$name,
    cells = if (is.null(cells)) fingerprint$cellCount else cells$spec$nObs,
    genes = fingerprint$varCount,
    open_cells = if (is.null(cells)) NULL else list(open = cells$spec$k, of = cells$spec$nObs, method = cells$spec$method, seed = cells$spec$seed),
    version = manifest$app$version
  )
  reading <- list(
    value = function(hash) .join_value(get(hash, envir = skeletons)[[1L]], arrays, part_of),
    values_included = values_shared,
    dataset = facts
  )
  built <- lapply(trackers, .history_tracker, reading = reading)
  names(built) <- vapply(trackers, function(tracker) tracker$id, "")
  structure(
    list(
      path = normalizePath(file, winslash = "/"),
      sha256 = sha256,
      manifest = manifest,
      created_at = manifest$createdAt,
      app = .history_build("AppStamp", manifest$app),
      actions = lapply(manifest$actions, function(stamp) .history_build("ActionStamp", stamp)),
      scope = .history_build("Scope", scope),
      dataset_fingerprint = fingerprint,
      source = manifest$source,
      requires = lapply(manifest$requires, function(keys) as.character(unlist(keys))),
      # Page chrome outside the history (the pressed selection tool): read, never applied.
      chrome = manifest$chrome,
      cells = cells,
      genes = genes,
      active_tracker = manifest$activeTracker,
      trackers = built
    ),
    class = "cellucid_session"
  )
}

# --- What a session holds ----------------------------------------------------------

.require_session <- function(session) {
  if (!inherits(session, "cellucid_session")) {
    .abort("session must be what cellucid_read_session() or cellucid_session_result() returned.")
  }
  session
}

cellucid_conversations <- function(session, tracker = NULL) {
  found <- if (inherits(tracker, "cellucid_tracker")) tracker else cellucid_tracker(session, tracker)
  lapply(found$chats, function(pair) list(conversation = pair[[1L]], transcript = pair[[2L]]))
}

print.cellucid_session <- function(x, ...) {
  fingerprint <- x$dataset_fingerprint
  active <- x$trackers[[x$active_tracker]]
  head <- cellucid_step(active, active$head_id)
  cat(
    "<Cellucid session> ", fingerprint$name %||% fingerprint$datasetId %||% "a dataset", " (",
    fingerprint$datasetId %||% "no id", "), ", format(fingerprint$cellCount, scientific = FALSE), " ",
    .plural(fingerprint$cellCount, "cell"), if (!is.null(x$cells)) paste0(" of ", x$cells$spec$nObs), "\n",
    "  saved:    ", x$created_at, " by Cellucid ", x$app$version, "\n",
    "  trackers: ", paste(vapply(x$trackers, function(item) {
      paste0(item$name, " (", length(item$steps), " ", .plural(length(item$steps), "step"), if (!is.null(item$removed)) ", removed", ")")
    }, ""), collapse = ", "), "\n",
    "  shown:    ", active$name, ", step ", head$id, ": ", head$label, "\n",
    sep = ""
  )
  invisible(x)
}
