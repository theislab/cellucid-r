# The viewer channel: what this package and the Cellucid viewer it opens say
# to each other. The web app's data/viewer-protocol.js is the source of truth;
# its document (the shared fixture tests/testthat/fixtures/viewer_protocol.json,
# kept byte for byte in the web app and in cellucid-python, SHA-256 pinned in
# each) is what tests/testthat/test-viewer-protocol.R checks these
# declarations and checks against.
#
# One WebSocket per viewer at <server>/_cellucid/ws. The viewer's first frame
# is hello {type, viewerId, viewerToken}; the host answers welcome, or closes
# with 4401. From then on:
#
# - commands (R -> viewer): text frames, one JSON object each, with exactly
#   the keys .viewer_command_fields names; no credentials. A command that
#   names cells or carries a column (a highlight, a visibility change, a
#   pushed field) is a binary frame whose header has exactly the keys
#   .viewer_binary_command_fields names, padded with spaces so its cells
#   start on a 4-byte boundary;
# - events (viewer -> R): text frames with exactly the keys
#   .viewer_event_fields names, and binary frames for a selection, a result
#   and a session bundle.
#
# Every cell either side names is an original position: the 0-based row of
# the dataset's obs, whatever subsample the viewer opened. A result is one
# thing a person made in the viewer (a highlight group, a categorical or
# continuous field, a gene set, the visible or the open cells): its body is
# its cells (uint32, when its header counts them) then its values (one per
# named cell, or per cell of the dataset when it names none).
#
# A binary frame, either way, is a little-endian uint32 header length, that
# many bytes of UTF-8 JSON header, then the body; a body of cells is
# little-endian uint32 cell indices, strictly ascending.
#
# Every frame is checked exactly -- exact key sets, exact types -- because a
# page in a browser is not a trusted caller.

.viewer_channel_path <- "_cellucid/ws"

.viewer_close_codes <- c(
  normal = 1000L, goingAway = 1001L, protocolViolation = 1008L,
  unauthorized = 4401L, replaced = 4409L
)

.viewer_command_fields <- list(
  ping = c("type", "requestId"),
  debug_snapshot = c("type", "requestId"),
  requestSessionBundle = c("type", "requestId"),
  clearHighlights = "type",
  setColorBy = c("type", "field"),
  resetCamera = "type",
  freeze = "type",
  setDataMapping = c("type", "requestId", "mapping"),
  setGeneSet = c("type", "requestId", "name", "description", "genes")
)

# The header of each binary command; its body is the cells. setVisibility of
# every cell is all = TRUE and no cells.
.viewer_binary_command_fields <- list(
  highlight = c("type", "color", "label"),
  setVisibility = c("type", "visible", "all"),
  setField = c("type", "requestId", "key", "kind", "categories", "nObs", "cells", "values")
)

# The keys a command is given as here, its cells under `cells` (NULL for
# setVisibility of every cell, and for a setField whose values cover every
# cell); a setField's values under `values`, in the type `valueType` names.
.viewer_binary_command_arguments <- list(
  highlight = c("type", "cells", "color", "label"),
  setVisibility = c("type", "cells", "visible"),
  setField = c("type", "requestId", "key", "kind", "categories", "nObs", "cells", "values", "valueType")
)

.viewer_command_types <- c(names(.viewer_command_fields), names(.viewer_binary_command_fields))

.viewer_event_fields <- list(
  hover = c("type", "cell", "position"),
  click = c("type", "cell", "button", "shift", "ctrl"),
  ready = c("type", "n_cells", "n_obs", "dimensions"),
  pong = c("type", "requestId", "t"),
  debug_snapshot = c("type", "requestId", "ts", "locationHref", "origin", "serverUrl", "connected", "userAgent"),
  command_error = c("type", "command", "reason"),
  data_mapping = c("type", "requestId", "origin", "mapping"),
  color_by = c("type", "field", "source", "symbol", "cells_read"),
  result_removed = c("type", "epoch", "revision", "kind", "id"),
  results = c("type", "epoch", "dataset", "revision", "nObs", "items")
)

.viewer_binary_event_fields <- list(
  selection = c("type", "source"),
  session_bundle = c("type", "requestId"),
  result = c("type", "requestId", "epoch", "revision", "kind", "id", "nObs", "meta", "cells", "values")
)

# One page's account of one dataset it opened (viewer-protocol.json
# results$epochPattern): every result event carries it, and a host keeps the
# results of one epoch.
.viewer_result_epoch_pattern <- "^[0-9a-f]{16}$"
.is_epoch <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && grepl(.viewer_result_epoch_pattern, value)
}

# Each result kind, with the exact keys of its meta.
.viewer_result_meta_fields <- list(
  highlight = c("label", "tool", "pageId", "pageName", "color", "enabled"),
  categorical = c("key", "categories", "origin"),
  continuous = c("key", "origin"),
  geneSet = c("name", "description", "genes", "symbols", "origin"),
  visibility = character(0),
  subset = "digest"
)
.viewer_highlight_tools <- c("annotation", "category", "combined", "knn", "lasso", "proximity", "range")
.viewer_result_origins <- c("analysis", "markers", "notebook", "user")
.viewer_field_kinds <- c("categorical", "continuous")
# Each value type of a body, and its size in bytes.
.viewer_value_types <- c(uint8 = 1L, uint16 = 2L, uint32 = 4L, float32 = 4L, float64 = 8L)

.viewer_data_mapping_origins <- c("default", "best_guess", "session", "user", "notebook")

# The longest command_error reason this host accepts (the viewer sends at most
# 240 characters; the margin keeps a change of the viewer's own limit from
# breaking the channel).
.viewer_reason_limit <- 500L

# The body of GET /_cellucid/protocol.
.viewer_protocol_document <- function() {
  list(
    channel = .viewer_channel_path,
    commands = I(sort(.viewer_command_types, method = "radix")),
    binaryCommands = I(sort(names(.viewer_binary_command_fields), method = "radix")),
    events = I(sort(c(names(.viewer_event_fields), names(.viewer_binary_event_fields)), method = "radix"))
  )
}

.protocol_error <- function(message) {
  structure(
    class = c("cellucid_protocol_error", "error", "condition"),
    list(message = message, call = NULL)
  )
}

.refuse_frame <- function(...) stop(.protocol_error(paste0(...)))

.is_exact_text <- function(value, whitespace = FALSE) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value) &&
    identical(value, trimws(value)) && !grepl("[\001-\037\177]", value) &&
    (whitespace || !grepl("[[:space:]]", value))
}

.is_count <- function(value) {
  is.numeric(value) && length(value) == 1L && !is.na(value) && is.finite(value) &&
    value >= 0 && value == floor(value)
}

.require_exact_keys <- function(value, keys, label) {
  if (!is.list(value) || is.null(names(value)) || anyDuplicated(names(value)) ||
      !setequal(names(value), keys) || length(value) != length(keys)) {
    .refuse_frame(label, " must carry exactly ", paste(keys, collapse = ", "), ".")
  }
  invisible(value)
}

# A text frame from the browser: httpuv hands its UTF-8 bytes over as native
# text, which is read as UTF-8 (.as_utf8()), so a name the browser sends
# reads the same in every R session.
.parse_json_frame <- function(text, label) {
  text <- .as_utf8(text)
  if (is.na(text)) .refuse_frame(label, " must be UTF-8.")
  value <- tryCatch(jsonlite::parse_json(text, simplifyVector = FALSE), error = function(error) NULL)
  if (!is.list(value) || is.null(names(value))) .refuse_frame(label, " must be one JSON object.")
  value
}

# The hello frame: the viewer id and token it carries.
.parse_viewer_hello <- function(text) {
  hello <- .parse_json_frame(text, "The first frame")
  .require_exact_keys(hello, c("type", "viewerId", "viewerToken"), "The hello frame")
  if (!identical(hello$type, "hello") || !.is_exact_text(hello$viewerId) || !.is_exact_text(hello$viewerToken)) {
    .refuse_frame("The first frame must be hello with a viewer id and token.")
  }
  hello
}

# One text event, checked: its type and its fields without `type`.
.parse_viewer_event <- function(text) {
  event <- .parse_json_frame(text, "An event")
  type <- event$type
  if (!.is_exact_text(type) || is.null(.viewer_event_fields[[type]])) {
    .refuse_frame("Unknown event type ", jsonlite::toJSON(type, auto_unbox = TRUE), ".")
  }
  .require_exact_keys(event, .viewer_event_fields[[type]], paste0("A ", type, " event"))
  bad <- function(what) .refuse_frame(type, " ", what, ".")
  switch(type,
    hover = {
      if (!is.null(event$cell) && !.is_count(event$cell)) bad("cell must be null or a cell index")
      position <- event$position
      if (!is.null(position)) {
        .require_exact_keys(position, c("x", "y", "z"), "A hover position")
        for (axis in c("x", "y", "z")) {
          value <- position[[axis]]
          if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) bad("position must be finite numbers")
        }
      }
    },
    click = {
      if (!.is_count(event$cell)) bad("cell must be a cell index")
      if (!.is_count(event$button) || !event$button %in% 0:2) bad("button must be 0, 1 or 2")
      if (!isTRUE(event$shift) && !isFALSE(event$shift)) bad("shift must be true or false")
      if (!isTRUE(event$ctrl) && !isFALSE(event$ctrl)) bad("ctrl must be true or false")
    },
    color_by = {
      if (is.null(event$field) != is.null(event$source)) bad("names a field and its source, or neither")
      if (!is.null(event$field)) {
        if (!.is_exact_text(event$field, whitespace = TRUE)) bad("field must be exact text")
        if (!event$source %in% c("obs", "var")) bad("source must be obs or var")
      }
      if (!is.null(event$symbol)) {
        if (!identical(event$source, "var")) bad("names a symbol only for a gene")
        if (!.is_exact_text(event$symbol, whitespace = TRUE)) bad("symbol must be exact text")
      }
      # A gene read for a sample of the open cells: the first `cells_read` of them.
      if (!is.null(event$cells_read)) {
        if (!identical(event$source, "var")) bad("names the cells read only for a gene")
        if (!.is_count(event$cells_read)) bad("cells_read must be a count")
      }
    },
    result_removed = {
      if (!.is_epoch(event$epoch)) bad("epoch must be 16 lowercase hexadecimal characters")
      if (!.is_count(event$revision) || event$revision < 1) bad("revision must be a positive whole number")
      if (!.is_exact_text(event$kind) || is.null(.viewer_result_meta_fields[[event$kind]])) bad("kind must be a result kind")
      if (!.is_exact_text(event$id, whitespace = TRUE)) bad("id must be exact text")
    },
    results = {
      if (!.is_epoch(event$epoch)) bad("epoch must be 16 lowercase hexadecimal characters")
      if (!is.null(event$dataset) && !.is_exact_text(event$dataset)) bad("dataset must be a dataset id or null")
      if (!.is_count(event$revision)) bad("revision must be a count")
      if (!.is_count(event$nObs)) bad("nObs must be a count")
      if (!is.list(event$items) || !is.null(names(event$items))) bad("items must be a list")
      if (is.null(event$dataset) && length(event$items) > 0L) bad("of a dataset the host did not open must list nothing")
      for (item in event$items) {
        .require_exact_keys(item, c("kind", "id", "revision"), "A results item")
        if (!.is_exact_text(item$kind) || is.null(.viewer_result_meta_fields[[item$kind]])) bad("item kind must be a result kind")
        if (!.is_exact_text(item$id, whitespace = TRUE) || !.is_count(item$revision) || item$revision < 1) {
          bad("item must name a result and its revision")
        }
      }
    },
    ready = {
      if (!.is_count(event$n_cells)) bad("n_cells must be a count")
      if (!.is_count(event$n_obs) || event$n_obs < event$n_cells) bad("n_obs must be a count of at least n_cells")
      if (!.is_count(event$dimensions) || !event$dimensions %in% 1:3) bad("dimensions must be 1, 2 or 3")
    },
    pong = {
      if (!.is_exact_text(event$requestId)) bad("requestId must be exact text")
      if (!.is_count(event$t)) bad("t must be a non-negative whole number")
    },
    debug_snapshot = {
      for (key in c("requestId", "ts", "locationHref", "origin", "serverUrl")) {
        if (!.is_exact_text(event[[key]], whitespace = TRUE)) bad(paste(key, "must be exact text"))
      }
      if (!isTRUE(event$connected) && !isFALSE(event$connected)) bad("connected must be true or false")
      if (!is.null(event$userAgent) && !.is_exact_text(event$userAgent, whitespace = TRUE)) {
        bad("userAgent must be null or text")
      }
    },
    command_error = {
      if (!.is_exact_text(event$command) || !event$command %in% .viewer_command_types) {
        bad("command must name a command of this protocol")
      }
      if (!.is_exact_text(event$reason, whitespace = TRUE) || nchar(event$reason) > .viewer_reason_limit) {
        bad(paste0("reason must be text of at most ", .viewer_reason_limit, " ", .plural(.viewer_reason_limit, "character")))
      }
    },
    data_mapping = {
      if (!.is_exact_text(event$origin) || !event$origin %in% .viewer_data_mapping_origins) {
        bad(paste0("origin must be one of ", paste(.viewer_data_mapping_origins, collapse = ", ")))
      }
      if (identical(event$origin, "notebook")) {
        if (!.is_exact_text(event$requestId)) bad("requestId must name the setDataMapping it answers")
      } else if (!is.null(event$requestId)) {
        bad("requestId must be null unless the event answers setDataMapping")
      }
      .require_data_mapping_object(event$mapping, "data_mapping mapping")
    }
  )
  list(type = type, payload = event[setdiff(names(event), "type")])
}

# A data mapping as it crosses the wire: one JSON object whose kind is
# "anndata". The viewer applies the full data-mapping rules (and answers a
# command it refuses with command_error); the host checks the shape.
.require_data_mapping_object <- function(mapping, label) {
  if (!is.list(mapping) || is.null(names(mapping)) || !identical(mapping$kind, "anndata")) {
    .refuse_frame(label, " must be a data mapping: a JSON object whose kind is \"anndata\".")
  }
  invisible(mapping)
}

# One binary event: header length, JSON header, body. A selection's body is
# little-endian uint32 cell indices, strictly ascending, read in one
# vectorised pass however many cells it holds.
.parse_viewer_binary_event <- function(frame) {
  if (!is.raw(frame) || length(frame) < 4L) .refuse_frame("A binary event is shorter than its header length.")
  header_length <- .read_uint32(frame[1:4])
  if (header_length < 2 || 4 + header_length > length(frame)) {
    .refuse_frame("A binary event's header length is outside the frame.")
  }
  header <- .utf8_text(frame[4 + seq_len(header_length)])
  if (is.na(header)) .refuse_frame("A binary event header must be UTF-8.")
  header <- .parse_json_frame(header, "A binary event header")
  type <- header$type
  if (!.is_exact_text(type) || is.null(.viewer_binary_event_fields[[type]])) {
    .refuse_frame("Unknown binary event type.")
  }
  .require_exact_keys(header, .viewer_binary_event_fields[[type]], paste0("A ", type, " header"))
  body <- if (length(frame) > 4 + header_length) frame[(5 + header_length):length(frame)] else raw(0)
  if (identical(type, "selection")) {
    if (!.is_exact_text(header$source)) .refuse_frame("selection source must be exact text.")
    list(type = type, payload = list(cells = .decode_selection_cells(body), source = header$source))
  } else if (identical(type, "result")) {
    list(type = type, payload = .decode_result(header, body))
  } else {
    if (!.is_exact_text(header$requestId)) .refuse_frame("session_bundle requestId must be exact text.")
    list(type = type, payload = list(requestId = header$requestId), body = body)
  }
}

.read_uint32 <- function(bytes) {
  value <- readBin(bytes, "integer", n = 1L, size = 4L, endian = "little")
  if (value < 0L) as.double(value) + 4294967296 else as.double(value)
}

# 0-based cell indices, as doubles when one exceeds R's integer range.
.decode_selection_cells <- function(body) {
  n <- length(body) %/% 4L
  if (n == 0L || length(body) %% 4L != 0L) {
    .refuse_frame("A selection body must be one or more little-endian uint32 cell indices.")
  }
  cells <- readBin(body, "integer", n = n, size = 4L, endian = "little")
  if (anyNA(cells) || any(cells < 0L)) {
    cells <- as.double(cells)
    negative <- is.na(cells) | cells < 0
    # NA_integer_ is the bit pattern 0x80000000, 2^31 as a uint32.
    cells[is.na(cells)] <- 2147483648
    cells[!is.na(cells) & cells < 0] <- cells[!is.na(cells) & cells < 0] + 4294967296
  }
  if (n > 1L && !all(cells[-1L] > cells[-n])) {
    .refuse_frame("A selection must list its cells strictly ascending, each once.")
  }
  cells
}

.encode_uint32 <- function(values) {
  values <- as.double(values)
  values[values >= 2147483648] <- values[values >= 2147483648] - 4294967296
  # -2^31 has no R integer: it becomes NA_integer_, whose bit pattern is
  # exactly 0x80000000.
  writeBin(suppressWarnings(as.integer(values)), raw(), size = 4L, endian = "little")
}

# One command, checked and serialized as the frame the viewer reads: a
# command that names cells as one binary frame (a raw vector), every other one
# as JSON text.
.viewer_command_frame <- function(command) {
  type <- command$type
  if (!.is_exact_text(type) || !type %in% .viewer_command_types) {
    .abort("Unknown viewer command ", .format_value(type), ".")
  }
  binary <- type %in% names(.viewer_binary_command_fields)
  keys <- if (binary) .viewer_binary_command_arguments[[type]] else .viewer_command_fields[[type]]
  if (!setequal(names(command), keys) || length(command) != length(keys)) {
    .abort("A ", type, " command carries exactly ", paste(keys, collapse = ", "), ".")
  }
  if ("requestId" %in% keys && !.is_exact_text(command$requestId)) {
    .abort(type, " requestId must be exact text.")
  }
  if (identical(type, "highlight")) {
    if (!is.character(command$color) || length(command$color) != 1L || !grepl("^#[0-9a-fA-F]{6}$", command$color)) {
      .abort("highlight color must be a six-digit hex colour such as \"#e41a1c\".")
    }
    if (!is.null(command$label) && !.is_exact_text(command$label, whitespace = TRUE)) {
      .abort("highlight label must be exact text, or NULL.")
    }
    return(.binary_command_frame(list(type = type, color = command$color, label = command$label), command$cells))
  }
  if (identical(type, "setField")) return(.set_field_frame(command))
  if (identical(type, "setVisibility")) {
    if (!isTRUE(command$visible) && !isFALSE(command$visible)) {
      .abort("setVisibility visible must be TRUE or FALSE.")
    }
    every <- is.null(command$cells)
    return(.binary_command_frame(list(type = type, visible = command$visible, all = every), command$cells, every))
  }
  if (identical(type, "setGeneSet")) {
    if (!.is_exact_text(command$name, whitespace = TRUE)) .abort("setGeneSet name must be exact text.")
    if (!is.character(command$description) || length(command$description) != 1L || is.na(command$description)) {
      .abort("setGeneSet description must be one string.")
    }
    genes <- unlist(command$genes)
    if (length(genes) > 0L && (!is.character(genes) || anyNA(genes) || !all(vapply(genes, .is_exact_text, logical(1))) || anyDuplicated(genes))) {
      .abort("setGeneSet genes must be gene names, each once, without spaces.")
    }
    command$genes <- I(as.character(genes))
  } else if (identical(type, "setColorBy")) {
    if (!.is_exact_text(command$field, whitespace = TRUE)) .abort("setColorBy field must be exact text.")
  } else if (identical(type, "setDataMapping")) {
    tryCatch(
      .require_data_mapping_object(command$mapping, "setDataMapping mapping"),
      cellucid_protocol_error = function(error) .abort(conditionMessage(error))
    )
  }
  as.character(jsonlite::toJSON(command, auto_unbox = TRUE, null = "null", digits = NA))
}

# A binary command frame: the header length, the JSON header padded with
# spaces to a multiple of 4 bytes (the cells then start on a 4-byte
# boundary, where the viewer reads them in place), and the cells -- 0-based
# cell indices, strictly ascending (each once), at least one, each below 2^32
# -- as little-endian uint32, checked and written in one pass (src/chunks.c).
# `every`: a setVisibility of every cell, which names none.
.binary_command_frame <- function(header, cells, every = FALSE) {
  text <- as.character(jsonlite::toJSON(header, auto_unbox = TRUE, null = "null"))
  text <- paste0(text, strrep(" ", (4L - nchar(text, type = "bytes") %% 4L) %% 4L))
  frame <- .native_cells_frame(charToRaw(text), if (!every) cells, allow_empty = every)
  if (is.null(frame)) {
    .abort(header$type, " cells must be cell indices, strictly ascending, at least one.")
  }
  frame
}

# --- Results and pushed fields --------------------------------------------------

# Values of one value type as little-endian bytes, and back.
.encode_wire_values <- function(values, type) {
  switch(type,
    uint8 = writeBin(as.integer(values), raw(), size = 1L),
    uint16 = writeBin(as.integer(values), raw(), size = 2L, endian = "little"),
    uint32 = .encode_uint32(values),
    float32 = writeBin(as.double(values), raw(), size = 4L, endian = "little"),
    float64 = writeBin(as.double(values), raw(), size = 8L, endian = "little")
  )
}

.decode_wire_values <- function(bytes, type, n) {
  if (n == 0L) return(if (type %in% c("float32", "float64")) double(0) else integer(0))
  switch(type,
    uint8 = readBin(bytes, "integer", n = n, size = 1L, signed = FALSE),
    uint16 = readBin(bytes, "integer", n = n, size = 2L, signed = FALSE, endian = "little"),
    uint32 = {
      values <- as.double(readBin(bytes, "integer", n = n, size = 4L, endian = "little"))
      values[is.na(values)] <- 2147483648
      values[values < 0] <- values[values < 0] + 4294967296
      values
    },
    float32 = readBin(bytes, "double", n = n, size = 4L, endian = "little"),
    float64 = readBin(bytes, "double", n = n, size = 8L, endian = "little")
  )
}

.is_text_list <- function(value, whitespace = TRUE, nonempty = FALSE) {
  values <- unlist(value)
  is.list(value) && (!nonempty || length(value) > 0L) && length(values) == length(value) &&
    (length(values) == 0L || (is.character(values) && !anyNA(values) && !anyDuplicated(values) &&
      (whitespace || all(vapply(values, .is_exact_text, logical(1))))))
}

# Check one result's header against the rules of its kind.
.require_result_header <- function(header) {
  bad <- function(what) .refuse_frame("A result ", what, ".")
  kind <- header$kind
  if (!.is_exact_text(kind) || is.null(.viewer_result_meta_fields[[kind]])) bad("kind must be a result kind")
  if (!is.null(header$requestId) && !.is_exact_text(header$requestId)) bad("requestId must be exact text or null")
  if (!.is_epoch(header$epoch)) bad("epoch must be 16 lowercase hexadecimal characters")
  if (!.is_count(header$revision) || header$revision < 1) bad("revision must be a positive whole number")
  if (!.is_exact_text(header$id, whitespace = TRUE)) bad("id must be exact text")
  if (!.is_count(header$nObs)) bad("nObs must be a cell count")
  meta <- header$meta
  fields <- .viewer_result_meta_fields[[kind]]
  if (length(fields) == 0L) {
    if (!is.list(meta) || length(meta) != 0L) bad("meta of this kind must be empty")
  } else {
    .require_exact_keys(meta, fields, paste0("A ", kind, " result's meta"))
  }
  if ("origin" %in% fields && !identical(length(meta$origin), 1L)) bad("must name its origin")
  if ("origin" %in% fields && !meta$origin %in% .viewer_result_origins) bad("origin must be one of the protocol's")
  cells <- header$cells
  values <- header$values
  if (!is.null(cells) && !.is_count(cells)) bad("cells must be a count or null")
  if (!is.null(values) && (!.is_exact_text(values) || !values %in% names(.viewer_value_types))) bad("values must be a value type or null")
  switch(kind,
    highlight = {
      if (!.is_exact_text(meta$label, whitespace = TRUE) || !meta$tool %in% .viewer_highlight_tools) bad("names its label and tool")
      if (!.is_exact_text(meta$pageId) || !.is_exact_text(meta$pageName, whitespace = TRUE)) bad("names its page")
      if (!is.character(meta$color) || !grepl("^#[0-9a-fA-F]{6}$", meta$color) || !(isTRUE(meta$enabled) || isFALSE(meta$enabled))) {
        bad("carries its page colour and whether it is enabled")
      }
      if (is.null(cells) || cells < 1 || !is.null(values)) bad("of a highlight names its cells and no values")
    },
    categorical = {
      if (!.is_exact_text(meta$key, whitespace = TRUE) || !.is_text_list(meta$categories, nonempty = TRUE)) bad("names its key and categories")
      if (!values %in% c("uint8", "uint16", "uint32") || identical(cells, 0L)) bad("of a categorical carries codes")
    },
    continuous = {
      if (!.is_exact_text(meta$key, whitespace = TRUE)) bad("names its key")
      if (!values %in% c("float32", "float64") || identical(cells, 0L)) bad("of a continuous field carries floats")
    },
    geneSet = {
      if (!.is_exact_text(meta$name, whitespace = TRUE) || !is.character(meta$description)) bad("names its gene set")
      if (!.is_text_list(meta$genes, whitespace = FALSE)) bad("lists its genes, each once")
      if (!is.null(meta$symbols) && (!is.list(meta$symbols) || length(meta$symbols) != length(meta$genes) ||
          !all(vapply(meta$symbols, function(symbol) is.character(symbol) && length(symbol) == 1L, TRUE)))) {
        bad("names a symbol for each of its genes, or none")
      }
      if (!identical(as.integer(cells), 0L) || !is.null(values)) bad("of a gene set names no cells and no values")
    },
    visibility = {
      if (!identical(header$id, "visible") || !is.null(values)) bad("of the visible cells is \"visible\" with no values")
    },
    subset = {
      if (!identical(header$id, "open") || !is.null(values) || identical(as.integer(cells), 0L)) bad("of the open cells is \"open\"")
      if (!is.character(meta$digest) || !grepl("^[0-9a-f]{16}$", meta$digest)) bad("carries the open set's digest")
    }
  )
  invisible(header)
}

# One result: its header checked, its cells (0-based, ascending, or NULL for
# every cell) and values read.
.decode_result <- function(header, body) {
  .require_result_header(header)
  n_cells <- header$cells
  type <- header$values
  cell_bytes <- if (is.null(n_cells)) 0 else 4 * n_cells
  n_values <- if (is.null(type)) 0 else if (is.null(n_cells)) header$nObs else n_cells
  value_bytes <- if (is.null(type)) 0 else n_values * .viewer_value_types[[type]]
  if (length(body) != cell_bytes + value_bytes) {
    .refuse_frame("The ", header$kind, " result body must be ", cell_bytes + value_bytes, " ", .plural(cell_bytes + value_bytes, "byte"), ", not ", length(body), ".")
  }
  cells <- NULL
  if (!is.null(n_cells)) {
    cells <- if (n_cells == 0) numeric(0) else .decode_selection_cells(body[seq_len(cell_bytes)])
    if (length(cells) > 0L && cells[[length(cells)]] >= header$nObs) .refuse_frame("A result names a cell beyond its dataset.")
  }
  values <- NULL
  if (!is.null(type)) {
    values <- .decode_wire_values(if (value_bytes == 0) raw(0) else body[cell_bytes + seq_len(value_bytes)], type, n_values)
    if (identical(header$kind, "categorical")) {
      missing <- 2^(8 * .viewer_value_types[[type]]) - 1
      if (any(values >= length(header$meta$categories) & values != missing)) .refuse_frame("A categorical result code names no category.")
    } else if (any(is.infinite(values))) {
      .refuse_frame("A continuous result holds an infinity.")
    }
  }
  list(
    requestId = header$requestId, epoch = header$epoch, revision = header$revision, kind = header$kind, id = header$id,
    nObs = header$nObs, meta = header$meta, cells = cells, values = values, value_type = type
  )
}

# A setField command: a column for the viewer to add as a field. `cells` are
# 0-based ascending positions or NULL (values for every cell), `values` codes
# (`valueType` uint8/16/32, the type's largest value "no category") or
# numbers (float32/64, NA or NaN "no value").
.set_field_frame <- function(command) {
  if (!.is_exact_text(command$requestId)) .abort("setField requestId must be exact text.")
  if (!.is_exact_text(command$key, whitespace = TRUE)) .abort("setField key must be exact text.")
  if (!.is_exact_text(command$kind) || !command$kind %in% .viewer_field_kinds) .abort("setField kind must be categorical or continuous.")
  n_obs <- command$nObs
  if (!.is_count(n_obs) || n_obs < 1) .abort("setField nObs must be the dataset's cell count.")
  type <- command$valueType
  categorical <- identical(command$kind, "categorical")
  if (!.is_exact_text(type) || !type %in% if (categorical) c("uint8", "uint16", "uint32") else c("float32", "float64")) {
    .abort("setField values must be ", if (categorical) "uint8, uint16 or uint32 codes" else "float32 or float64 numbers", ".")
  }
  cells <- command$cells
  if (!is.null(cells) && (length(cells) == 0L || anyNA(cells) || any(cells < 0) || any(cells >= n_obs) ||
                          is.unsorted(cells, strictly = TRUE))) {
    .abort("setField cells must be cell positions, strictly ascending, inside the dataset.")
  }
  values <- command$values
  count <- if (is.null(cells)) n_obs else length(cells)
  if (length(values) != count) .abort("setField holds ", length(values), " ", .plural(length(values), "value"), " for ", count, " ", .plural(count, "cell"), ".")
  if (categorical) {
    categories <- as.character(unlist(command$categories))
    missing <- 2^(8 * .viewer_value_types[[type]]) - 1
    if (length(categories) == 0L || anyNA(categories) || anyDuplicated(categories) || length(categories) > missing) {
      .abort("setField categories must be text, each once, fewer than the code type holds.")
    }
    if (anyNA(values) || any(values < 0 | (values >= length(categories) & values != missing))) {
      .abort("setField codes must name a category, or be ", missing, " for none.")
    }
    categories <- I(categories)
  } else {
    categories <- NULL
    values <- as.double(values)
    values[is.na(values)] <- NaN
    if (any(is.infinite(values))) .abort("setField values must be finite numbers, or NA for none.")
  }
  header <- list(
    type = "setField", requestId = command$requestId, key = command$key, kind = command$kind,
    categories = categories, nObs = n_obs, cells = if (is.null(cells)) NULL else length(cells), values = type
  )
  text <- as.character(jsonlite::toJSON(header, auto_unbox = TRUE, null = "null", digits = NA))
  text <- paste0(text, strrep(" ", (4L - nchar(text, type = "bytes") %% 4L) %% 4L))
  header_bytes <- charToRaw(text)
  c(.encode_uint32(length(header_bytes)), header_bytes,
    if (!is.null(cells)) .encode_uint32(cells), .encode_wire_values(values, type))
}
