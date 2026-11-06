# The viewer channel: what this package and the Cellucid viewer it opens say
# to each other. The web app's data/viewer-protocol.js is the source of truth;
# its document (the shared fixture tests/testthat/fixtures/viewer_protocol.json,
# kept byte for byte in the web app and in cellucid-python, SHA-256 pinned in
# each) is what tests/testthat/test-viewer-protocol.R checks these
# declarations and the frames made here against.
#
# One WebSocket per viewer at <server>/_cellucid/ws. The viewer's first frame
# is hello {type, viewerId, viewerToken}; the host answers welcome, or closes
# with 4401. From then on:
#
# - commands (R -> viewer): text frames, one JSON object each, with exactly
#   the keys .viewer_command_fields names; no credentials. Everything R
#   changes or reads in the viewer is one `action` of the catalogue
#   (inst/action_catalog.json, R/actions.R): {type, requestId, action, params,
#   turn, turnLabel} -- the notebook turn it belongs to and what the turn was
#   asked to do (cellucid_turn()), both null outside one. An action whose
#   parameters hold typed arrays (cells, a column's values) is a binary frame
#   whose header has exactly the keys .viewer_binary_command_fields names --
#   `arrays` lists where each array goes in params, its type and length --
#   padded with spaces so its body starts on a 4-byte boundary; the body is
#   the arrays' elements, little-endian, in order (.encode_viewer_action()).
#   A session file pushed into the viewer (cellucid_load_session()) is one
#   binary load_session frame: {type, requestId, name}, then the file;
# - agent commands (R -> viewer): an agent app's tool list, call, stop and
#   answer (agent_list, agent_call, agent_cancel, agent_answer), carried for a
#   local MCP server (cellucid_mcp(), R/mcp.R, and R/agent-link.R); the
#   viewer answers them with agent events, which reach the viewer's agent
#   port and never the hooks;
# - events (viewer -> R): text frames with exactly the keys
#   .viewer_event_fields names, and binary frames for a selection, a result,
#   an action's answer (action_result: its outcome as UTF-8 JSON), a
#   session bundle, and an agent's tools and a call's result (agent_tools,
#   agent_result: UTF-8 JSON). An action is answered by one action_result, or by one
#   command_error naming its requestId; a load_session by session_loaded (the
#   tracker shown and its head) or a command_error. Every change of the
#   viewer's history is a history_changed event: its kind
#   (.viewer_history_change_kinds), the tracker in use, the step shown, the
#   step recorded (a "step" change) and the active branch.
#
# Every cell either side names is an original position: the 0-based row of
# the dataset's obs, whatever subsample the viewer opened. A result is one
# thing a person made in the viewer (a highlight group, a categorical or
# continuous field, a gene set, the visible or the open cells): its body is
# its cells (uint32, when its header counts them) then its values (one per
# named cell, or per cell of the dataset when it names none).
#
# A binary frame, either way, is a little-endian uint32 header length, that
# many bytes of UTF-8 JSON header, then the body.
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
  freeze = "type",
  action = c("type", "requestId", "action", "params", "turn", "turnLabel"),
  # An agent app on this computer driving the viewer through a local MCP
  # server (cellucid_mcp(), R/mcp.R; cellucid-python's cellucid mcp): its
  # tool list, a call, a stop, and the person's answer given in the app. The
  # viewer runs each call as the agent through a local connection, under the
  # person's permissions, and answers with the agent events below
  # (R/agent-link.R carries them).
  agent_list = c("type", "requestId", "connection"),
  agent_call = c("type", "requestId", "tool", "args", "agent", "connection"),
  agent_cancel = c("type", "requestId"),
  agent_answer = c("type", "requestId", "allowed")
)

# The keys of an agent command's connection: the local connection's id (one
# MCP session of one app) and the app's own name, or null; the pattern of the
# id; an agent's tool name as the tool list names it ("color_by"). An agent's
# name and an app's are exact text of any length.
.viewer_agent_connection_fields <- c("id", "client")
.viewer_agent_connection_id_pattern <- "^[A-Za-z0-9_-]{1,64}$"
.viewer_agent_tool_pattern <- "^[a-z][a-z0-9_]*$"

# The header of each binary command; its body is the arrays `arrays` lists.
.viewer_binary_command_fields <- list(
  action = c("type", "requestId", "action", "params", "turn", "turnLabel", "arrays"),
  load_session = c("type", "requestId", "name")
)

# The keys of one entry of a binary action's `arrays`, and of an
# action_result's `step`.
.viewer_action_array_fields <- c("path", "type", "length")
.viewer_action_step_fields <- c("id", "label")

# A dotted action type, as the catalogue names one ("color.by").
.viewer_action_type_pattern <- "^[a-z][a-z0-9]*(?:-[a-z0-9]+)*(?:\\.[a-z][a-z0-9]*(?:-[a-z0-9]+)*)+$"

.viewer_command_types <- sort(unique(c(names(.viewer_command_fields), names(.viewer_binary_command_fields))),
                              method = "radix")

.viewer_event_fields <- list(
  hover = c("type", "cell", "position"),
  click = c("type", "cell", "button", "shift", "ctrl"),
  ready = c("type", "n_cells", "n_obs", "dimensions"),
  pong = c("type", "requestId", "t"),
  debug_snapshot = c("type", "requestId", "ts", "locationHref", "origin", "serverUrl", "connected", "userAgent"),
  command_error = c("type", "command", "requestId", "reason"),
  data_mapping = c("type", "origin", "mapping"),
  color_by = c("type", "field", "source", "symbol", "cells_read"),
  result_removed = c("type", "epoch", "revision", "kind", "id"),
  results = c("type", "epoch", "dataset", "revision", "nObs", "items"),
  history_changed = c("type", "kind", "tracker", "head", "step", "branch"),
  session_loaded = c("type", "requestId", "tracker", "head"),
  # The person is asked about an agent's call (the viewer's confirm card),
  # and answered it; the tools an agent may use changed (the person's
  # permission settings).
  agent_asking = c("type", "requestId", "message", "detail"),
  agent_answered = c("type", "requestId", "allowed"),
  agent_tools_changed = "type"
)

# Every agent event: handed to the viewer's agent port (R/agent-link.R),
# never to the hooks, viewer$last or a printed notice.
.viewer_agent_events <- c("agent_tools", "agent_result", "agent_asking", "agent_answered", "agent_tools_changed")

# What a history_changed says changed: a step recorded, the step shown moved
# (undo, redo, a jump, a branch switched to), an edit of the history's
# structure (a note, a bookmark, a branch renamed, locked, archived or
# deleted), a tracker started, used, renamed or removed, a session loaded.
.viewer_history_change_kinds <- c("step", "moved", "edited", "tracker", "loaded")
# The keys of session_loaded's tracker.
.viewer_session_loaded_tracker_fields <- c("id", "name")

.viewer_binary_event_fields <- list(
  selection = c("type", "source"),
  session_bundle = c("type", "requestId"),
  result = c("type", "epoch", "revision", "kind", "id", "nObs", "meta", "cells", "values"),
  action_result = c("type", "requestId", "step"),
  # An agent's tools (a JSON array) and a call's MCP tool result
  # ({content, isError, ...}), each as the body's UTF-8 JSON.
  agent_tools = c("type", "requestId"),
  agent_result = c("type", "requestId")
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
# Each value type of a body, and its size in bytes.
.viewer_value_types <- c(uint8 = 1L, uint16 = 2L, uint32 = 4L, float32 = 4L, float64 = 8L)

.viewer_data_mapping_origins <- c("default", "best_guess", "session", "user", "notebook")

# The longest command_error reason this host accepts (the viewer sends at most
# the protocol's commandErrorReasonLimit, 1000 characters, room for an
# action's refusal naming what it expected; the margin keeps a change of the
# viewer's own limit from breaking the channel), as cellucid-python accepts.
.viewer_reason_limit <- 4000L

# The body of GET /_cellucid/protocol.
.viewer_protocol_document <- function() {
  list(
    channel = .viewer_channel_path,
    commands = I(.viewer_command_types),
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
      # The request the refused command named: null for one that names none.
      if (!is.null(event$requestId) && !.is_exact_text(event$requestId)) bad("requestId must be exact text or null")
      if (!.is_exact_text(event$reason, whitespace = TRUE) || nchar(event$reason) > .viewer_reason_limit) {
        bad(paste0("reason must be text of at most ", .viewer_reason_limit, " ", .plural(.viewer_reason_limit, "character")))
      }
    },
    data_mapping = {
      if (!.is_exact_text(event$origin) || !event$origin %in% .viewer_data_mapping_origins) {
        bad(paste0("origin must be one of ", paste(.viewer_data_mapping_origins, collapse = ", ")))
      }
      .require_data_mapping_object(event$mapping, "data_mapping mapping")
    },
    history_changed = {
      if (!.is_exact_text(event$kind) || !event$kind %in% .viewer_history_change_kinds) {
        bad(paste0("kind must be one of ", paste(.viewer_history_change_kinds, collapse = ", ")))
      }
      for (key in c("tracker", "head")) if (!.is_exact_text(event[[key]])) bad(paste(key, "must be exact text"))
      for (key in c("step", "branch")) {
        if (!is.null(event[[key]]) && !.is_exact_text(event[[key]])) bad(paste(key, "must be exact text or null"))
      }
    },
    session_loaded = {
      if (!.is_exact_text(event$requestId)) bad("requestId must be exact text")
      .require_exact_keys(event$tracker, .viewer_session_loaded_tracker_fields, "A session_loaded tracker")
      if (!.is_exact_text(event$tracker$id)) bad("tracker id must be exact text")
      if (!.is_exact_text(event$tracker$name, whitespace = TRUE)) bad("tracker name must be exact text")
      if (!.is_exact_text(event$head)) bad("head must be exact text")
    },
    agent_asking = {
      if (!.is_exact_text(event$requestId)) bad("requestId must be exact text")
      if (!.is_exact_text(event$message, whitespace = TRUE)) bad("message must be exact text")
      if (!is.null(event$detail) && !.is_exact_text(event$detail, whitespace = TRUE)) bad("detail must be exact text or null")
    },
    agent_answered = {
      if (!.is_exact_text(event$requestId)) bad("requestId must be exact text")
      if (!isTRUE(event$allowed) && !isFALSE(event$allowed)) bad("allowed must be true or false")
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
  switch(type,
    selection = {
      if (!.is_exact_text(header$source)) .refuse_frame("selection source must be exact text.")
      list(type = type, payload = list(cells = .decode_selection_cells(body), source = header$source))
    },
    result = list(type = type, payload = .decode_result(header, body)),
    action_result = list(type = type, payload = .decode_action_result(header, body)),
    agent_tools = list(type = type, payload = .decode_agent_json(header, body, "tools")),
    agent_result = list(type = type, payload = .decode_agent_json(header, body, "result")),
    session_bundle = {
      if (!.is_exact_text(header$requestId)) .refuse_frame("session_bundle requestId must be exact text.")
      list(type = type, payload = list(requestId = header$requestId), body = body)
    }
  )
}

# The answer to one action: the request it answers, the history step it made
# ({id, label}, or NULL: a read, or a change that changed nothing), and its
# outcome, the body's UTF-8 JSON (NULL when the action had nothing to report;
# a read's result), as jsonlite reads JSON (simplifyVector = FALSE).
.decode_action_result <- function(header, body) {
  if (!.is_exact_text(header$requestId)) .refuse_frame("action_result requestId must be exact text.")
  step <- header$step
  if (!is.null(step)) {
    .require_exact_keys(step, .viewer_action_step_fields, "An action_result step")
    if (!.is_exact_text(step$id) || !.is_exact_text(step$label, whitespace = TRUE)) {
      .refuse_frame("An action_result step names its id and label.")
    }
  }
  text <- .utf8_text(body)
  if (is.na(text)) .refuse_frame("An action_result body must be UTF-8 JSON.")
  outcome <- tryCatch(
    list(value = jsonlite::parse_json(text, simplifyVector = FALSE)),
    error = function(error) NULL
  )
  if (is.null(outcome)) .refuse_frame("An action_result body must be one JSON value.")
  list(requestId = header$requestId, step = step, outcome = outcome$value)
}

# One binary agent event (agent_tools, agent_result): its request, and its
# body's UTF-8 JSON as `key` (as jsonlite reads JSON, simplifyVector =
# FALSE): the tools, each an object with its name, or an MCP tool result.
.decode_agent_json <- function(header, body, key) {
  type <- header$type
  if (!.is_exact_text(header$requestId)) .refuse_frame(type, " requestId must be exact text.")
  text <- .utf8_text(body)
  value <- if (is.na(text)) NULL else tryCatch(
    list(value = jsonlite::parse_json(text, simplifyVector = FALSE)),
    error = function(error) NULL
  )
  if (is.null(value)) .refuse_frame("The body of ", type, " ", header$requestId, " is not UTF-8 JSON.")
  value <- value$value
  if (identical(type, "agent_tools")) {
    named <- function(tool) is.list(tool) && !is.null(names(tool)) && is.character(tool[["name"]]) && length(tool[["name"]]) == 1L
    if (!is.list(value) || !is.null(names(value)) || !all(vapply(value, named, TRUE))) {
      .refuse_frame("agent_tools body must be a list of tools, each with its name.")
    }
  } else if (!is.list(value) || is.null(names(value)) || !is.list(value[["content"]]) || !is.null(names(value[["content"]])) ||
             !(isTRUE(value[["isError"]]) || isFALSE(value[["isError"]]))) {
    .refuse_frame("agent_result body must be an MCP tool result: {content, isError, ...}.")
  }
  stats::setNames(list(header$requestId, value), c("requestId", key))
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

# A command other than an action, checked and serialized as the text frame
# the viewer reads: ping, debug_snapshot and requestSessionBundle name a
# request, freeze none.
.viewer_command_frame <- function(command) {
  type <- command$type
  if (!.is_exact_text(type) || !type %in% setdiff(names(.viewer_command_fields), "action")) {
    .abort("Unknown viewer command ", .format_value(type), ".")
  }
  keys <- .viewer_command_fields[[type]]
  if (!setequal(names(command), keys) || length(command) != length(keys)) {
    .abort("A ", type, " command carries exactly ", paste(keys, collapse = ", "), ".")
  }
  if ("requestId" %in% keys && !.is_exact_text(command$requestId)) {
    .abort(type, " requestId must be exact text.")
  }
  if (startsWith(type, "agent_")) .require_agent_command(command)
  .js_json(command[keys])
}

# An agent command's fields beyond its keys, checked as the viewer checks
# them (requireViewerAgentCommand) and as cellucid-python does
# (_require_agent_command).
.require_agent_command <- function(command) {
  type <- command$type
  if (type %in% c("agent_list", "agent_call")) .require_agent_connection(command$connection, type)
  if (identical(type, "agent_call")) {
    tool <- command$tool
    if (!is.character(tool) || length(tool) != 1L || is.na(tool) || !grepl(.viewer_agent_tool_pattern, tool)) {
      .abort("agent_call names its tool as the tool list does, such as 'color_by'.")
    }
    args <- command$args
    if (!is.list(args) || is.null(names(args))) .abort("agent_call args must be one object.")
    if (!.is_exact_text(command$agent, whitespace = TRUE)) .abort("agent_call agent must be exact text.")
  }
  if (identical(type, "agent_answer") && !isTRUE(command$allowed) && !isFALSE(command$allowed)) {
    .abort("agent_answer allowed must be TRUE or FALSE.")
  }
  invisible(command)
}

.require_agent_connection <- function(connection, label) {
  fields <- .viewer_agent_connection_fields
  if (!is.list(connection) || is.null(names(connection)) || anyDuplicated(names(connection)) ||
      !setequal(names(connection), fields) || length(connection) != length(fields)) {
    .abort(label, " connection must be {id, client}.")
  }
  id <- connection$id
  if (!is.character(id) || length(id) != 1L || is.na(id) || !grepl(.viewer_agent_connection_id_pattern, id)) {
    .abort(label, " connection id must be 1 to 64 letters, digits, \"_\" or \"-\".")
  }
  client <- connection$client
  if (!is.null(client) && !.is_exact_text(client, whitespace = TRUE)) {
    .abort(label, " connection client must be exact text or null.")
  }
  invisible(connection)
}

# --- Actions ------------------------------------------------------------------

# A typed array of an action's parameters: `values` sent as `type` (one of
# .viewer_value_types) in a binary frame's body, never as JSON. Cells are
# original positions as uint32 (an x-cells parameter); a column's values any
# value type (x-values).
.typed_array <- function(values, type) {
  if (!.is_exact_text(type) || !type %in% names(.viewer_value_types)) {
    .abort("A typed array's type is one of ", paste(names(.viewer_value_types), collapse = ", "), ".")
  }
  if (!is.numeric(values) && !is.logical(values)) .abort("A typed array holds numbers.")
  values <- as.vector(values)
  if (type %in% c("float32", "float64")) {
    values <- as.double(values)
    values[is.na(values)] <- NaN
  } else {
    largest <- 2^(8 * .viewer_value_types[[type]]) - 1
    extent <- if (length(values) > 0L) .native_extent(as.double(values), whole = TRUE) else c(0, 0, 0, 0, 1)
    if (extent[[3L]] > 0 || extent[[1L]] < 0 || extent[[2L]] > largest || extent[[5L]] != 1) {
      .abort("A ", type, " array holds whole numbers from 0 to ", format(largest, scientific = FALSE), ".")
    }
  }
  structure(list(values = values, type = type), class = "cellucid_typed_array")
}

.typed_array_bytes <- function(array) .encode_wire_values(array$values, array$type)

# The typed arrays of an action's parameters taken out: the JSON rest, and
# each array with the path to it (keys and 0-based list indices, ending in a
# key), in the order met (data/viewer-protocol.js splitActionParams). An array
# goes under an object's key, never directly in a list.
.split_action_params <- function(params) {
  arrays <- list()
  dropped <- structure(list(), class = "cellucid_dropped")
  visit <- function(node, path, in_list) {
    if (inherits(node, "cellucid_typed_array")) {
      if (in_list) {
        .abort("An action parameter list holds a typed array at ", paste(unlist(path), collapse = "."), "; put it under a key.")
      }
      arrays[[length(arrays) + 1L]] <<- list(path = path, array = node)
      return(dropped)
    }
    if (is.list(node) && is.null(names(node))) {
      return(lapply(seq_along(node), function(at) visit(node[[at]], c(path, list(at - 1L)), TRUE)))
    }
    if (is.list(node)) {
      out <- structure(list(), names = character())
      for (key in names(node)) {
        kept <- visit(node[[key]], c(path, list(key)), FALSE)
        if (!inherits(kept, "cellucid_dropped")) out[key] <- list(kept)
      }
      return(out)
    }
    node
  }
  list(params = visit(params, list(), FALSE), arrays = arrays)
}

# One `action` command as R sends it: a JSON text frame when its parameters
# hold no typed array, else a binary frame (a raw vector) whose header lists
# the arrays and whose body holds them, the header padded with spaces to a
# multiple of 4 bytes (the web app's encodeViewerAction makes the same bytes
# from the same parameters: the fixture's binaryCommands). `turn` and
# `turn_label` name the notebook turn it belongs to (cellucid_turn()), both
# NULL outside one.
.encode_viewer_action <- function(request_id, action, params, turn = NULL, turn_label = NULL) {
  if (!.is_exact_text(request_id)) .abort("An action command's requestId must be exact text.")
  if (!.is_exact_text(action) || !grepl(.viewer_action_type_pattern, action, perl = TRUE)) {
    .abort("An action is named by its dotted type, such as \"color.by\"; not ", .format_value(action), ".")
  }
  if (!is.list(params) || is.null(names(params))) .abort("An action's parameters are one named list.")
  if (is.null(turn) != is.null(turn_label)) .abort("An action names its turn and the turn's label together, or neither.")
  if (!is.null(turn) && (!.is_exact_text(turn) || !.is_exact_text(turn_label, whitespace = TRUE))) {
    .abort("An action's turn and turn label must be exact text.")
  }
  split <- .split_action_params(params)
  header <- list(type = "action", requestId = request_id, action = action, params = split$params, turn = turn,
                 turnLabel = turn_label)
  if (length(split$arrays) == 0L) return(.js_json(header))
  header$arrays <- lapply(split$arrays, function(entry) {
    list(path = entry$path, type = entry$array$type, length = length(entry$array$values))
  })
  .binary_command_frame(header, unlist(lapply(split$arrays, function(entry) .typed_array_bytes(entry$array)), use.names = FALSE))
}

# A binary command: the header's compact JSON padded with spaces to a
# multiple of 4 bytes, its length first, then the body.
.binary_command_frame <- function(header, body) {
  bytes <- charToRaw(enc2utf8(.js_json(header)))
  bytes <- c(bytes, rep(charToRaw(" "), (4L - length(bytes) %% 4L) %% 4L))
  c(.encode_uint32(length(bytes)), bytes, body)
}

# One `load_session` command: a .cellucid-session file pushed into the viewer
# (cellucid_load_session()), named by its file name without a folder.
.encode_viewer_load_session <- function(request_id, name, file) {
  if (!.is_exact_text(request_id)) .abort("A load_session command's requestId must be exact text.")
  if (!.is_exact_text(name, whitespace = TRUE) || grepl("[/\\\\]", name)) {
    .abort("A load_session command names the file by its name, without a folder.")
  }
  if (!is.raw(file) || length(file) == 0L) .abort("A load_session command carries the session file's bytes.")
  .binary_command_frame(list(type = "load_session", requestId = request_id, name = name), file)
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
    epoch = header$epoch, revision = header$revision, kind = header$kind, id = header$id,
    nObs = header$nObs, meta = header$meta, cells = cells, values = values, value_type = type
  )
}
