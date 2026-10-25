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
#   names cells (a highlight, a visibility change) is a binary frame whose
#   header has exactly the keys .viewer_binary_command_fields names, padded
#   with spaces so its cells start on a 4-byte boundary;
# - events (viewer -> R): text frames with exactly the keys
#   .viewer_event_fields names, and binary frames for a selection and a
#   session bundle.
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
  setDataMapping = c("type", "requestId", "mapping")
)

# The header of each binary command; its body is the cells. setVisibility of
# every cell is all = TRUE and no cells.
.viewer_binary_command_fields <- list(
  highlight = c("type", "color"),
  setVisibility = c("type", "visible", "all")
)

# The keys a command is given as here, its cells under `cells` (NULL for
# setVisibility of every cell).
.viewer_binary_command_arguments <- list(
  highlight = c("type", "cells", "color"),
  setVisibility = c("type", "cells", "visible")
)

.viewer_command_types <- c(names(.viewer_command_fields), names(.viewer_binary_command_fields))

.viewer_event_fields <- list(
  hover = c("type", "cell", "position"),
  click = c("type", "cell", "button", "shift", "ctrl"),
  ready = c("type", "n_cells", "dimensions"),
  pong = c("type", "requestId", "t"),
  debug_snapshot = c("type", "requestId", "ts", "locationHref", "origin", "serverUrl", "connected", "userAgent"),
  command_error = c("type", "command", "reason"),
  data_mapping = c("type", "requestId", "origin", "mapping")
)

.viewer_binary_event_fields <- list(
  selection = c("type", "source"),
  session_bundle = c("type", "requestId")
)

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
    ready = {
      if (!.is_count(event$n_cells)) bad("n_cells must be a count")
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
    return(.binary_command_frame(list(type = type, color = command$color), command$cells))
  }
  if (identical(type, "setVisibility")) {
    if (!isTRUE(command$visible) && !isFALSE(command$visible)) {
      .abort("setVisibility visible must be TRUE or FALSE.")
    }
    every <- is.null(command$cells)
    return(.binary_command_frame(list(type = type, visible = command$visible, all = every), command$cells, every))
  }
  if (identical(type, "setColorBy")) {
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
  text <- as.character(jsonlite::toJSON(header, auto_unbox = TRUE))
  text <- paste0(text, strrep(" ", (4L - nchar(text, type = "bytes") %% 4L) %% 4L))
  frame <- .native_cells_frame(charToRaw(text), if (!every) cells, allow_empty = every)
  if (is.null(frame)) {
    .abort(header$type, " cells must be cell indices, strictly ascending, at least one.")
  }
  frame
}
