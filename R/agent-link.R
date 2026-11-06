# Local agent connections: how a local MCP server -- cellucid_mcp() here
# (R/mcp.R), cellucid mcp in cellucid-python -- reaches a Cellucid viewer.
#
# An agent app on this computer (Claude Code, Claude Desktop, Cursor) runs the
# server over stdio. Each of its tool calls is carried to a viewer page over
# that viewer's channel as an agent command (agent_list, agent_call,
# agent_cancel, agent_answer; R/viewer-protocol.R), and the page runs it
# exactly as it runs a connector agent's call: as the agent, through a
# connection of kind "local", under the person's permissions, asked in the
# confirm card, one turn the person can Stop and undo. The page answers with
# agent events (agent_tools, agent_result, agent_asking, agent_answered,
# agent_tools_changed), which reach the viewer's agent port below and never
# its hooks (cellucid_on()), viewer$last or a printed notice.
#
# Two ways to a viewer, one interface (a link, R/mcp.R):
#
# - in the same process (cellucid_mcp(open = <file>) serves the viewer
#   itself, or a viewer of the R session running cellucid_mcp()): the
#   viewer's port sends the commands on its channel and hands each answer to
#   whoever waits for it;
# - in another process (an R session, a notebook, a Python script): every
#   viewer with an HTTP server is published to agents on this computer. An R
#   viewer's own server (httpuv) answers its endpoint, a route beside its
#   data routes, and its record goes into the agents folder (.agents_dir()):
#   where the endpoint is and the secret token it takes, readable by this
#   user only. A notebook widget's viewer (IRkernel, viewer-widget.R) has no
#   server and is not published: an agent reaches it only through a viewer
#   with one. A viewer bound to one non-loopback interface is not published
#   either: the endpoint answers loopback peers only, and such a server has
#   no loopback address. The server reads the records (.agent_discover()) and
#   drives one over its endpoint (.remote_viewer_link()); cellucid-python
#   publishes its viewers the same way (cellucid/_agent_link.py), so either
#   package's server drives either package's viewer.
#
# The endpoint (<endpoint> is http://127.0.0.1:<port>/_cellucid/agent/<viewer id>):
#
# - GET <endpoint>: the viewer, {schema, viewer: {viewer, runtime, dataset,
#   url, connected}};
# - POST <endpoint>/command: one agent command, exactly as the viewer
#   protocol has it; 202 once sent, 409 with the reason when no page is
#   connected to the viewer, 400 for a command that is not one;
# - GET <endpoint>/events?consumer=<connection id>&wait=<seconds>: the agent
#   events for that consumer (the answers to its commands, and every
#   agent_tools_changed), oldest first, as {"events": [...]}; it waits up to
#   `wait` seconds (at most .agent_event_wait_limit_seconds) for one, without
#   holding R (the answer is a promise httpuv sends when it settles). A
#   command the page refused is {"type": "agent_failed", requestId, reason},
#   the viewer stopping {"type": "viewer_stopped", reason}.
#
# Every request carries Authorization: Bearer <token>; one from a browser
# (any Origin header), from another host, or naming another Host is refused
# (403), before its body is read. The token is compared in constant time.
# tests/testthat/fixtures/local_agents.json, shared byte for byte with
# cellucid-python, pins the folder, the record, the endpoint and the
# sentences.

# The version of a viewer record and of the endpoint's answers.
.agent_record_schema <- 2L

# The route every viewer's endpoint is under.
.agent_route <- "_cellucid/agent"

# The longest a GET .../events waits for an event before answering none.
.agent_event_wait_limit_seconds <- 25

# How long a consumer that stopped asking for its events keeps them: its
# server ended without saying so (killed, the app closed).
.agent_consumer_idle_seconds <- 120

# The largest command body the endpoint reads: one tool call written by a
# model (a model writes at most a few hundred kilobytes in one answer).
.agent_command_body_limit_bytes <- 4 * 1024 * 1024

# How long the endpoint lets a page list its tools.
.agent_list_seconds <- 60

# The keys of a viewer record, exactly (cellucid-python writes the same).
.agent_record_fields <- c("dataset", "endpoint", "host", "pid", "runtime", "schema", "started", "token", "url", "viewer")

# The sentences an agent reads (the shared fixture's): a call whose page went
# before it was answered; a call the page refused to read ({reason} is the
# page's); a call whose viewer stopped; an agents folder shared with others.
.agent_page_gone <- paste0(
  "The Cellucid viewer's page went (it was reloaded or closed) before this call ended, so whether it ",
  "was done is not known. Read the history (history_get) before trying again."
)
.agent_not_taken <- "The Cellucid viewer page did not take the call: {reason} Ask the person to reload the viewer's page."
.agent_stopped <- "The Cellucid viewer stopped ({reason}); nothing more was done."
.agent_folder_unsafe <- paste0(
  "The agents folder {folder} is not this user's alone ({why}), so no viewer is published there or ",
  "read from it. Make it yours and not writable by others (chmod 700), or set CELLUCID_AGENTS_DIR."
)

.agent_sentence <- function(template, ...) {
  values <- list(...)
  for (key in names(values)) template <- sub(paste0("{", key, "}"), values[[key]], template, fixed = TRUE)
  template
}

# CELLUCID_AGENTS=off (or 0, false, no; any case) publishes no viewer to
# agents on this computer.
.agent_opt_out_variable <- "CELLUCID_AGENTS"
.agent_opt_out_values <- c("off", "0", "false", "no")

.agents_disabled <- function() {
  tolower(trimws(Sys.getenv(.agent_opt_out_variable))) %in% .agent_opt_out_values
}

# This computer's name, as records carry it (cellucid-python's
# socket.gethostname()): a home folder shared between machines holds other
# computers' records, which a reader skips and never removes.
.agent_this_host <- function() Sys.info()[["nodename"]]

# The program a viewer of this package names itself by (the protocol's runtimes).
.agent_runtime <- "r"

.agent_clock <- function() proc.time()[["elapsed"]]

# No page shows a viewer now, or the viewer stopped: one sentence that says
# which (cellucid-python's ViewerUnavailableError).
.viewer_unavailable <- function(...) {
  error <- structure(
    class = c("cellucid_viewer_unavailable", "error", "condition"),
    list(message = .condition_text(...), call = NULL)
  )
  stop(error)
}

# What an agent reads when no page shows a viewer now: open its page again
# (the address its program printed: it carries the viewer's secret, which an
# agent is never given), or (a notebook widget, which has none) show its cell
# again.
.not_open_sentence <- function(description) {
  name <- .json_field(description[["dataset"]], "name") %||% "the dataset"
  url <- description[["url"]]
  if (!is.null(url)) {
    return(paste0(
      "No page shows the Cellucid viewer of ", name, " (served at ", url, ") now. Ask the person to open ",
      "its page again (the viewer address its program printed, or reload the page), then try again."
    ))
  }
  paste0(
    "No page shows the Cellucid viewer of ", name, " now. Ask the person to show it again ",
    "(run the notebook cell that shows it), then try again."
  )
}

# A field of a JSON object read from outside, by its exact name: NULL when
# `value` is no object or lacks it (`$` would match a prefix).
.json_field <- function(value, key) {
  if (is.list(value) && !is.null(names(value))) value[[key]] else NULL
}

.is_single_text <- function(value) is.character(value) && length(value) == 1L && !is.na(value)

# =============================================================================
# One viewer's end, in its own process: its agent port
# =============================================================================

# The viewer's agent port, made with the viewer (.start_viewer()): the calls
# waiting for their events, the tool lists asked for, who hears that the
# tools changed, why it stopped, and its publication.
.agent_port <- function(viewer) {
  port <- viewer$agent
  if (is.null(port)) {
    port <- new.env(parent = emptyenv())
    port$calls <- new.env(parent = emptyenv())
    port$lists <- new.env(parent = emptyenv())
    port$tools_listeners <- new.env(parent = emptyenv())
    port$listener_count <- 0L
    port$stopped <- NULL
    port$published <- NULL
    viewer$agent <- port
  }
  port
}

.agent_hear <- function(listener, kind, data) {
  tryCatch(listener(kind, data), error = function(error) {
    .warn("A Cellucid agent call's listener failed: ", conditionMessage(error))
  })
  invisible(NULL)
}

# One agent event from the page (`payload` without its type).
.agent_port_deliver <- function(viewer, type, payload) {
  port <- .agent_port(viewer)
  if (identical(type, "agent_tools_changed")) {
    for (key in ls(port$tools_listeners)) {
      listener <- get0(key, envir = port$tools_listeners, inherits = FALSE)
      if (is.function(listener)) {
        tryCatch(listener(), error = function(error) {
          .warn("A Cellucid agent tools listener failed: ", conditionMessage(error))
        })
      }
    }
    return(invisible(NULL))
  }
  id <- payload$requestId
  if (identical(type, "agent_tools")) {
    .agent_port_settle_list(port, id, "tools", payload$tools)
    return(invisible(NULL))
  }
  listener <- get0(id, envir = port$calls, inherits = FALSE)
  if (identical(type, "agent_result") && !is.null(listener)) rm(list = id, envir = port$calls)
  if (is.null(listener)) return(invisible(NULL))
  switch(type,
    agent_result = .agent_hear(listener, "result", list(result = payload$result)),
    agent_asking = .agent_hear(listener, "asking", list(message = payload$message, detail = payload$detail)),
    agent_answered = .agent_hear(listener, "answered", list(allowed = payload$allowed))
  )
  invisible(NULL)
}

# The page refused an agent command (its command_error): an old page, or a
# command it could not read.
.agent_port_command_failed <- function(viewer, payload) {
  port <- .agent_port(viewer)
  id <- payload$requestId
  if (is.null(id)) return(invisible(NULL))
  reason <- .agent_sentence(.agent_not_taken, reason = payload$reason)
  if (identical(payload$command, "agent_list")) {
    .agent_port_settle_list(port, id, "failed", reason)
    return(invisible(NULL))
  }
  listener <- get0(id, envir = port$calls, inherits = FALSE)
  if (!is.null(listener)) {
    rm(list = id, envir = port$calls)
    .agent_hear(listener, "failed", list(reason = reason))
  }
  invisible(NULL)
}

# The page holding the channel went (closed, reloaded, replaced by another):
# no answer of a call it was running will come, so each fails now with the
# sentence an agent reads, and so does every tool list asked of it.
.agent_port_page_gone <- function(viewer) {
  port <- .agent_port(viewer)
  ids <- ls(port$calls, all.names = TRUE)
  listeners <- mget(ids, envir = port$calls)
  rm(list = ids, envir = port$calls)
  for (listener in listeners) .agent_hear(listener, "failed", list(reason = .agent_page_gone))
  for (id in ls(port$lists, all.names = TRUE)) .agent_port_settle_list(port, id, "failed", .agent_page_gone)
  invisible(NULL)
}

# A page holds the channel again: its tool list may say other things (its
# settings), so whoever lists the tools is told.
.agent_port_page_back <- function(viewer) {
  .agent_port_deliver(viewer, "agent_tools_changed", list())
}

# The viewer stopped: everything waiting ends with `reason`.
.agent_port_stopped <- function(viewer, reason) {
  port <- .agent_port(viewer)
  port$stopped <- reason
  ids <- ls(port$calls, all.names = TRUE)
  listeners <- mget(ids, envir = port$calls)
  rm(list = ids, envir = port$calls)
  for (listener in listeners) .agent_hear(listener, "stopped", list(reason = reason))
  for (id in ls(port$lists, all.names = TRUE)) {
    .agent_port_settle_list(port, id, "stopped", paste0("The Cellucid viewer stopped (", reason, ")."))
  }
  invisible(NULL)
}

.agent_port_connected <- function(viewer) {
  port <- .agent_port(viewer)
  is.null(port$stopped) && isTRUE(viewer$running) && !is.null(viewer$socket)
}

# The viewer as an agent connection names it (cellucid_status): its id, the
# program serving it, its dataset, where it is served (scheme, host and port,
# never its page address, which carries the viewer's token an agent is not
# given; NULL for a notebook widget, shown in a cell), and whether a page
# shows it now.
.agent_port_describe <- function(viewer) {
  list(
    viewer = viewer$viewer_id,
    runtime = .agent_runtime,
    dataset = list(id = viewer$dataset_id, name = viewer$dataset_name),
    url = if (is.null(viewer$widget)) viewer$url,
    connected = .agent_port_connected(viewer)
  )
}

.agent_port_require_open <- function(viewer) {
  port <- .agent_port(viewer)
  if (!is.null(port$stopped)) .viewer_unavailable("The Cellucid viewer stopped (", port$stopped, ").")
  if (!.agent_port_connected(viewer)) .viewer_unavailable(.not_open_sentence(.agent_port_describe(viewer)))
  invisible(port)
}

# Ask the page for the tools as `connection` may use them; `done(kind, value)`
# hears "tools" and the list, or "failed" and the sentence that says why
# (the page refused, the viewer stopped, `timeout` seconds passed).
.agent_port_list <- function(viewer, connection, timeout, done) {
  port <- .agent_port_require_open(viewer)
  id <- paste0(connection$id, "-list-", .random_hex(4L))
  slot <- new.env(parent = emptyenv())
  slot$done <- done
  assign(id, slot, envir = port$lists)
  slot$cancel <- later::later(function() {
    .agent_port_settle_list(port, id, "failed", "The Cellucid viewer did not list its tools in time.")
  }, timeout)
  tryCatch(
    .viewer_send(viewer, list(type = "agent_list", requestId = id, connection = connection)),
    error = function(error) {
      if (exists(id, envir = port$lists, inherits = FALSE)) rm(list = id, envir = port$lists)
      slot$cancel()
      stop(error)
    }
  )
  invisible(id)
}

.agent_port_settle_list <- function(port, id, kind, value) {
  slot <- get0(id, envir = port$lists, inherits = FALSE)
  if (is.null(slot)) return(invisible(NULL))
  rm(list = id, envir = port$lists)
  if (is.function(slot$cancel)) slot$cancel()
  # Anything but a list of tools is a sentence (a stop is one too).
  if (!identical(kind, "tools")) kind <- "failed"
  tryCatch(slot$done(kind, value), error = function(error) {
    .warn("A Cellucid agent tool list's listener failed: ", conditionMessage(error))
  })
  invisible(NULL)
}

# Send one call; `listener(kind, data)` hears its events: "asking" and
# "answered" while the person is asked (any number of times), then exactly
# one of "result" (the MCP tool result), "failed" (the page refused the
# command: reason) or "stopped" (the viewer stopped: reason).
.agent_port_start_call <- function(viewer, request_id, tool, args, agent, connection, listener) {
  port <- .agent_port_require_open(viewer)
  assign(request_id, listener, envir = port$calls)
  tryCatch(
    .viewer_send(viewer, list(
      type = "agent_call", requestId = request_id, tool = tool, args = args, agent = agent, connection = connection
    )),
    error = function(error) {
      if (exists(request_id, envir = port$calls, inherits = FALSE)) rm(list = request_id, envir = port$calls)
      stop(error)
    }
  )
  invisible(NULL)
}

.agent_port_cancel <- function(viewer, request_id) {
  port <- .agent_port(viewer)
  waiting <- exists(request_id, envir = port$calls, inherits = FALSE)
  if (waiting) rm(list = request_id, envir = port$calls)
  if (waiting && .agent_port_connected(viewer)) {
    try(.viewer_send(viewer, list(type = "agent_cancel", requestId = request_id)), silent = TRUE)
  }
  invisible(NULL)
}

.agent_port_answer <- function(viewer, request_id, allowed) {
  if (.agent_port_connected(viewer)) {
    .viewer_send(viewer, list(type = "agent_answer", requestId = request_id, allowed = isTRUE(allowed)))
  }
  invisible(NULL)
}

# Hear when the page says its tools changed; returns the function that stops it.
.agent_port_add_tools_listener <- function(viewer, listener) {
  port <- .agent_port(viewer)
  port$listener_count <- port$listener_count + 1L
  key <- sprintf("l%08d", port$listener_count)
  assign(key, listener, envir = port$tools_listeners)
  function() {
    if (exists(key, envir = port$tools_listeners, inherits = FALSE)) rm(list = key, envir = port$tools_listeners)
    invisible(NULL)
  }
}

# =============================================================================
# The agents folder: one record per published viewer
# =============================================================================

.is_absolute_path <- function(path) {
  if (.Platform$OS.type == "windows") return(grepl("^([A-Za-z]:[/\\\\]|[/\\\\]{2})", path))
  startsWith(path, "/")
}

# The agents folder refused (not this user's alone, or named relatively):
# one sentence that says why (cellucid-python's AgentsFolderError).
.agents_folder_error <- function(message) {
  error <- structure(class = c("cellucid_agents_folder", "error", "condition"), list(message = message, call = NULL))
  stop(error)
}

# The folder of this user's published viewers: $CELLUCID_AGENTS_DIR when set
# (an absolute path), else <user cache>/cellucid/agents: %LOCALAPPDATA% on
# Windows, ~/Library/Caches on macOS, $XDG_CACHE_HOME or ~/.cache elsewhere --
# the folder cellucid-python reads and writes (not tools::R_user_dir(), which
# names R's own), so either package's server finds either package's viewers.
.agents_dir <- function() {
  override <- Sys.getenv("CELLUCID_AGENTS_DIR")
  if (nzchar(override)) {
    if (!.is_absolute_path(override)) .agents_folder_error("CELLUCID_AGENTS_DIR must be an absolute path")
    return(override)
  }
  if (.Platform$OS.type == "windows") {
    local <- Sys.getenv("LOCALAPPDATA")
    home <- Sys.getenv("USERPROFILE")
    if (!nzchar(home)) home <- path.expand("~")
    base <- if (nzchar(local) && .is_absolute_path(local)) local else file.path(home, "AppData", "Local")
  } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
    base <- file.path(path.expand("~"), "Library", "Caches")
  } else {
    xdg <- Sys.getenv("XDG_CACHE_HOME")
    base <- if (nzchar(xdg) && .is_absolute_path(xdg)) xdg else file.path(path.expand("~"), ".cache")
  }
  file.path(base, "cellucid", "agents")
}

# On Linux and macOS the agents folder must be this user's and writable by no
# one else: a record names an endpoint and its secret, and a folder others
# can write could hand the server a record of theirs. A home folder shared
# between machines is fine (records name their host).
.agents_require_private <- function(folder) {
  if (.Platform$OS.type != "unix") return(invisible(folder))
  info <- file.info(folder, extra_cols = TRUE)
  why <- if (!identical(info$uname, Sys.info()[["effective_user"]])) {
    "another user owns it"
  } else if (bitwAnd(as.integer(info$mode), strtoi("022", 8L)) != 0L) {
    "others may write in it"
  }
  if (!is.null(why)) .agents_folder_error(.agent_sentence(.agent_folder_unsafe, folder = folder, why = why))
  invisible(folder)
}

# The agents folder, made (0700) when missing, and checked.
.agents_private_folder <- function() {
  folder <- .agents_dir()
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE, mode = "0700")
    if (!dir.exists(folder)) .abort("cannot create the agents folder ", folder)
    if (.Platform$OS.type == "unix" && !isTRUE(Sys.chmod(folder, "0700", use_umask = FALSE))) {
      .abort("cannot make the agents folder ", folder, " readable by this user only")
    }
  }
  .agents_require_private(folder)
  folder
}

# Write one record atomically (a temporary file in the folder, made 0600 by
# the umask it is created under, renamed over the record).
.agent_write_record <- function(record) {
  folder <- .agents_private_folder()
  target <- file.path(folder, paste0(record$viewer, ".json"))
  temporary <- tempfile(".", tmpdir = folder, fileext = ".tmp")
  mask <- Sys.umask("077")
  on.exit(Sys.umask(mask), add = TRUE)
  written <- tryCatch({
    connection <- file(temporary, "wb")
    tryCatch(writeBin(charToRaw(enc2utf8(.js_json(record))), connection), finally = close(connection))
    if (.Platform$OS.type == "unix" && !isTRUE(Sys.chmod(temporary, "0600", use_umask = FALSE))) .abort("chmod")
    file.rename(temporary, target)
  }, error = function(error) FALSE, warning = function(warning) FALSE)
  if (!isTRUE(written)) {
    unlink(temporary)
    .abort("cannot write ", target)
  }
  target
}

# Whether an endpoint is http on this computer's loopback (127.0.0.1,
# localhost, ::1) with a port: a record naming any other is not read.
.agent_is_loopback_endpoint <- function(endpoint) {
  match <- regmatches(endpoint, regexec("^http://(127\\.0\\.0\\.1|localhost|\\[::1\\]):([0-9]{1,5})(/[^?#]*)?$", endpoint,
                                        ignore.case = TRUE))[[1L]]
  length(match) > 0L && as.integer(match[[3L]]) %in% 1:65535
}

.is_text <- function(value) .is_single_text(value) && nzchar(value)

# A viewer record, or NULL for anything that is not exactly one: every field
# of its type (the shared fixture's record.types), its endpoint on this
# computer's loopback.
.agent_read_record <- function(path) {
  bytes <- tryCatch(readBin(path, "raw", file.size(path)), error = function(error) NULL, warning = function(warning) NULL)
  if (is.null(bytes)) return(NULL)
  text <- .utf8_text(bytes)
  if (is.na(text)) return(NULL)
  record <- tryCatch(jsonlite::parse_json(text, simplifyVector = FALSE), error = function(error) NULL)
  if (!is.list(record) || is.null(names(record)) || anyDuplicated(names(record)) ||
      !setequal(names(record), .agent_record_fields) || length(record) != length(.agent_record_fields)) {
    return(NULL)
  }
  pid <- record[["pid"]]
  dataset <- record[["dataset"]]
  url <- record[["url"]]
  valid <- identical(record[["schema"]], .agent_record_schema) &&
    all(vapply(c("viewer", "runtime", "host", "started", "endpoint", "token"), function(key) .is_text(record[[key]]), TRUE)) &&
    is.numeric(pid) && length(pid) == 1L && !is.na(pid) && pid == floor(pid) &&
    (is.null(url) || .is_text(url)) &&
    (is.null(dataset) || (is.list(dataset) && !is.null(names(dataset)) && !anyDuplicated(names(dataset)) &&
                            setequal(names(dataset), c("id", "name")) && length(dataset) == 2L &&
                            .is_text(dataset[["id"]]) && .is_text(dataset[["name"]]))) &&
    .agent_is_loopback_endpoint(record[["endpoint"]])
  if (!isTRUE(valid)) return(NULL)
  record
}

# Every published viewer of this user on this computer whose process still
# runs, newest first. Records of this computer's processes that ended are
# removed; another computer's (a home folder shared between machines) are
# left alone. A folder that is not this user's alone is refused (a
# cellucid_agents_folder error with the sentence).
.agent_discover <- function() {
  folder <- .agents_dir()
  if (!dir.exists(folder)) return(list())
  .agents_require_private(folder)
  paths <- sort(list.files(folder, pattern = "\\.json$", full.names = TRUE), method = "radix")
  host <- .agent_this_host()
  records <- list()
  for (path in paths) {
    record <- .agent_read_record(path)
    if (is.null(record) || !identical(record[["host"]], host)) next
    if (!isTRUE(.native_process_alive(record[["pid"]]))) {
      unlink(path)
      next
    }
    records[[length(records) + 1L]] <- record
  }
  if (length(records) == 0L) return(list())
  started <- vapply(records, function(record) record[["started"]], "")
  records[order(started, decreasing = TRUE, method = "radix")]
}

# =============================================================================
# A viewer's endpoint: its own server, reachable by agents on this computer
# =============================================================================

# The authority agents reach the viewer's server at: a loopback address it
# listens on, or NULL when it listens on one non-loopback interface only.
.agent_endpoint_authority <- function(viewer) {
  host <- tolower(viewer$host)
  address <- if (host %in% c("0.0.0.0", "", "localhost")) {
    "127.0.0.1"
  } else if (identical(host, "::")) {
    "[::1]"
  } else if (.is_loopback_host(host)) {
    if (grepl(":", host, fixed = TRUE)) paste0("[", host, "]") else host
  } else {
    return(NULL)
  }
  paste0(address, ":", viewer$port)
}

# Publish a viewer to agents on this computer (.start_viewer(), for a viewer
# with a server): its record in the agents folder; its commands and events
# answer on its own server. CELLUCID_AGENTS=off publishes nothing. A folder
# that cannot be written, or that is not this user's alone, leaves the
# viewer unpublished, said in a warning (as cellucid-python logs it): the
# viewer works on, and cellucid_mcp(open =) still reaches its own.
.agent_publish <- function(viewer) {
  port <- .agent_port(viewer)
  if (.agents_disabled()) {
    port$unpublished <- paste0(.agent_opt_out_variable, " turns publication off")
    return(invisible(NULL))
  }
  authority <- .agent_endpoint_authority(viewer)
  if (is.null(authority)) {
    port$unpublished <- "it listens on no loopback address"
    return(invisible(NULL))
  }
  published <- new.env(parent = emptyenv())
  published$viewer <- viewer
  published$token <- .random_hex(32L)
  published$consumers <- new.env(parent = emptyenv())
  published$stopped <- NULL
  published$record_path <- NULL
  published$endpoint <- paste0("http://", authority, "/", .agent_route, "/", viewer$viewer_id)
  failure <- tryCatch({
    description <- .agent_port_describe(viewer)
    published$record_path <- .agent_write_record(list(
      dataset = description$dataset,
      endpoint = published$endpoint,
      host = .agent_this_host(),
      pid = Sys.getpid(),
      runtime = .agent_runtime,
      schema = .agent_record_schema,
      started = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
      token = published$token,
      url = description$url,
      viewer = viewer$viewer_id
    ))
    NULL
  }, error = conditionMessage)
  if (!is.null(failure)) {
    port$unpublished <- failure
    .warn("This Cellucid viewer is not published to agents on this computer: ", failure)
    return(invisible(NULL))
  }
  published$remove_tools_listener <- .agent_port_add_tools_listener(viewer, function() {
    for (name in ls(published$consumers, all.names = TRUE)) {
      .agent_push(published, name, list(type = "agent_tools_changed"))
    }
  })
  port$published <- published
  port$unpublished <- NULL
  .agent_register_exit_hook()
  .agent_schedule_sweep(published)
  invisible(published)
}

# Consumers that went silent are swept on a clock (every quarter of
# .agent_consumer_idle_seconds), not only when another one asks, while the
# viewer is published.
.agent_schedule_sweep <- function(published) {
  published$cancel_sweep <- later::later(function() {
    if (!is.null(published$stopped)) return(invisible(NULL))
    .agent_sweep(published)
    .agent_schedule_sweep(published)
  }, .agent_consumer_idle_seconds / 4)
  invisible(NULL)
}

# Withdraw a published viewer: its record goes, its consumers hear it
# stopped (a consumer waiting for events now, at once). Returns what was
# published, or NULL.
.agent_withdraw <- function(viewer, reason = "the viewer was stopped") {
  port <- viewer$agent
  published <- if (is.null(port)) NULL else port$published
  if (is.null(published)) return(invisible(NULL))
  port$published <- NULL
  port$withdrawn <- published
  if (is.function(published$cancel_sweep)) published$cancel_sweep()
  published$remove_tools_listener()
  published$stopped <- reason
  for (name in ls(published$consumers, all.names = TRUE)) {
    entry <- get(name, envir = published$consumers)
    entry$events[[length(entry$events) + 1L]] <- list(type = "viewer_stopped", reason = reason)
    for (key in ls(entry$waiters, all.names = TRUE)) .agent_resolve_waiter(published, entry, key)
  }
  if (!is.null(published$record_path)) unlink(published$record_path)
  invisible(published)
}

# After .agent_withdraw(), before the server stops: R's event loop runs until
# every consumer that asked for events lately has been told the viewer
# stopped -- its waiting request answered with viewer_stopped, or its next
# request with 404 -- or `seconds` pass. A consumer that never hears it would
# keep asking a server that is gone.
.agent_tell_withdrawn <- function(published, seconds) {
  if (is.null(published)) return(invisible(NULL))
  recent <- function() {
    now <- .agent_clock()
    names <- ls(published$consumers, all.names = TRUE)
    names[vapply(names, function(name) {
      entry <- get(name, envir = published$consumers)
      !isTRUE(entry$told) && now - entry$seen <= .agent_event_wait_limit_seconds + 10
    }, TRUE)]
  }
  deadline <- .agent_clock() + seconds
  answered <- isTRUE(published$answered)
  while (length(recent()) > 0L && .agent_clock() < deadline) later::run_now(0.01)
  # httpuv writes an answer on its own thread once R handed it over.
  if (answered || isTRUE(published$answered)) {
    settle <- .agent_clock() + 0.05
    while (.agent_clock() < settle) later::run_now(0.01)
  }
  invisible(NULL)
}

# Every record this process published goes when R exits (a process killed
# leaves its records behind; .agent_discover() removes them).
.agent_exit_state <- new.env(parent = emptyenv())

.agent_register_exit_hook <- function() {
  if (isTRUE(.agent_exit_state$registered)) return(invisible(NULL))
  .agent_exit_state$registered <- TRUE
  reg.finalizer(.agent_exit_state, function(state) .agent_remove_records(), onexit = TRUE)
  invisible(NULL)
}

.agent_remove_records <- function() {
  for (id in ls(.viewer_registry, all.names = TRUE)) {
    viewer <- get0(id, envir = .viewer_registry, inherits = FALSE)
    path <- if (is.environment(viewer) && !is.null(viewer$agent)) viewer$agent$published$record_path
    if (!is.null(path)) unlink(path)
  }
  invisible(NULL)
}

# --- Its consumers ---------------------------------------------------------------

# One consumer (one MCP connection) of a published viewer's events; a
# consumer idle longer than .agent_consumer_idle_seconds goes, its calls
# stopped in the page.
# Drop the consumers silent for .agent_consumer_idle_seconds (their server
# ended without saying so): their calls are stopped in the page (a question
# withdrawn, a change rolled back) and their queued events let go.
.agent_sweep <- function(published, keep = NULL) {
  now <- .agent_clock()
  for (other in ls(published$consumers, all.names = TRUE)) {
    if (identical(other, keep)) next
    entry <- get(other, envir = published$consumers)
    if (now - entry$seen > .agent_consumer_idle_seconds && length(ls(entry$waiters)) == 0L) {
      for (request_id in entry$requests) try(.agent_port_cancel(published$viewer, request_id), silent = TRUE)
      entry$events <- list()
      rm(list = other, envir = published$consumers)
    }
  }
  invisible(NULL)
}

.agent_consumer <- function(published, name) {
  now <- .agent_clock()
  .agent_sweep(published, keep = name)
  entry <- get0(name, envir = published$consumers, inherits = FALSE)
  if (is.null(entry)) {
    entry <- new.env(parent = emptyenv())
    entry$events <- list()
    entry$requests <- character()
    entry$waiters <- new.env(parent = emptyenv())
    entry$waiter_count <- 0L
    entry$told <- FALSE
    assign(name, entry, envir = published$consumers)
  }
  entry$seen <- now
  entry
}

.agent_push <- function(published, name, event) {
  entry <- get0(name, envir = published$consumers, inherits = FALSE)
  if (is.null(entry)) return(invisible(NULL))
  entry$events[[length(entry$events) + 1L]] <- event
  if (event$type %in% c("agent_result", "agent_failed")) entry$requests <- setdiff(entry$requests, event$requestId)
  waiting <- sort(ls(entry$waiters, all.names = TRUE), method = "radix")
  if (length(waiting) > 0L) .agent_resolve_waiter(published, entry, waiting[[1L]])
  invisible(NULL)
}

# A consumer's events, oldest first, taken: with viewer_stopped once the
# viewer stopped.
.agent_take_events <- function(published, entry) {
  events <- entry$events
  entry$events <- list()
  if (!is.null(published$stopped)) {
    if (!any(vapply(events, function(event) identical(event$type, "viewer_stopped"), TRUE))) {
      events[[length(events) + 1L]] <- list(type = "viewer_stopped", reason = published$stopped)
    }
    entry$told <- TRUE
  }
  events
}

.agent_resolve_waiter <- function(published, entry, key) {
  waiter <- get0(key, envir = entry$waiters, inherits = FALSE)
  if (is.null(waiter)) return(invisible(NULL))
  rm(list = key, envir = entry$waiters)
  waiter$cancel()
  published$answered <- TRUE
  waiter$resolve(.agent_answer(200L, list(events = .agent_take_events(published, entry))))
  invisible(NULL)
}

# --- Its requests ----------------------------------------------------------------

.agent_answer <- function(status, value) {
  list(
    status = as.integer(status),
    headers = list(
      `Content-Type` = "application/json; charset=utf-8", `Cache-Control` = "no-store",
      # Never gzipped by httpuv: a loopback answer, sent as it is.
      `Content-Encoding` = "identity"
    ),
    body = charToRaw(enc2utf8(.js_json(value)))
  )
}

.is_agent_request <- function(req) {
  path <- req$PATH_INFO %||% ""
  path <- sub("^/+", "/", path)
  identical(path, paste0("/", .agent_route)) || startsWith(path, paste0("/", .agent_route, "/"))
}

# The checks every request to the endpoint passes before anything else, run
# from onHeaders (so a refused request's body is never read) and again with
# the request: a browser, another host or another Host is refused, and a
# command's declared body must be one this endpoint reads. NULL to go on.
.agent_request_refusal <- function(req) {
  if (!is.null(.request_header(req, "Origin"))) return(.agent_answer(403L, list(error = "Browsers cannot use this endpoint.")))
  host <- .parse_authority(.request_header(req, "Host"))
  if (is.null(host) || !.is_loopback_host(host$host) || !.is_loopback_host(req$REMOTE_ADDR %||% "")) {
    return(.agent_answer(403L, list(error = "Only this computer can use this endpoint.")))
  }
  if (identical(req$REQUEST_METHOD, "POST")) {
    # A whole number of bytes, or the command is not read.
    declared <- trimws(req$CONTENT_LENGTH %||% .request_header(req, "Content-Length") %||% "")
    length <- if (grepl("^[0-9]{1,15}$", declared)) as.numeric(declared) else NA
    if (is.na(length) || length <= 0 || length > .agent_command_body_limit_bytes) {
      return(.agent_answer(400L, list(error = "A command is one JSON object.")))
    }
  }
  NULL
}

.agent_query <- function(text) {
  text <- sub("^\\?", "", text %||% "")
  out <- list()
  for (pair in strsplit(text, "&", fixed = TRUE)[[1L]]) {
    if (!nzchar(pair)) next
    parts <- strsplit(pair, "=", fixed = TRUE)[[1L]]
    key <- .as_utf8(utils::URLdecode(gsub("+", " ", parts[[1L]], fixed = TRUE)))
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else ""
    value <- .as_utf8(utils::URLdecode(gsub("+", " ", value, fixed = TRUE)))
    if (!is.na(key) && !is.na(value) && is.null(out[[key]])) out[[key]] <- value
  }
  out
}

# One request to the endpoint (.on_call()): its answer, or a promise of it
# for a GET .../events that waits.
.agent_endpoint_call <- function(viewer, req) {
  refusal <- .agent_request_refusal(req)
  if (!is.null(refusal)) return(refusal)
  segments <- strsplit(gsub("^/+|/+$", "", req$PATH_INFO %||% ""), "/", fixed = TRUE)[[1L]]
  if (length(segments) < 3L || !identical(paste(segments[1:2], collapse = "/"), .agent_route)) {
    return(.agent_answer(404L, list(error = "Not found.")))
  }
  port <- .agent_port(viewer)
  published <- port$published
  query <- .agent_query(req$QUERY_STRING)
  if (is.null(published) || !identical(segments[[3L]], viewer$viewer_id)) {
    # A consumer of a viewer just withdrawn is told so (.agent_tell_withdrawn()).
    withdrawn <- port$withdrawn
    if (!is.null(withdrawn) && identical(segments[[3L]], viewer$viewer_id) && !is.null(query$consumer)) {
      entry <- get0(query$consumer, envir = withdrawn$consumers, inherits = FALSE)
      if (!is.null(entry)) entry$told <- TRUE
    }
    return(.agent_answer(404L, list(error = "This viewer is not published (it stopped).")))
  }
  given <- sub("^Bearer ", "", .request_header(req, "Authorization") %||% "")
  if (!.same_secret(enc2utf8(given), published$token)) {
    return(.agent_answer(401L, list(error = "The token does not match this viewer.")))
  }
  rest <- paste(segments[-(1:3)], collapse = "/")
  method <- req$REQUEST_METHOD
  if (identical(method, "GET")) {
    if (identical(rest, "")) {
      return(.agent_answer(200L, list(schema = .agent_record_schema, viewer = .agent_port_describe(viewer))))
    }
    if (identical(rest, "events")) return(.agent_events_answer(published, query))
    return(.agent_answer(404L, list(error = "Not found.")))
  }
  if (identical(method, "POST")) {
    if (!identical(rest, "command")) return(.agent_answer(404L, list(error = "Not found.")))
    return(.agent_command_answer(viewer, published, req))
  }
  .agent_answer(405L, list(error = paste0(method, " is not supported here.")))
}

.agent_events_answer <- function(published, query) {
  consumer <- query$consumer %||% ""
  wait <- suppressWarnings(as.numeric(query$wait %||% "0"))
  if (length(wait) != 1L || is.na(wait)) wait <- 0
  wait <- min(max(wait, 0), .agent_event_wait_limit_seconds)
  if (!nzchar(consumer)) return(.agent_answer(400L, list(error = "Name the consumer (?consumer=<connection id>).")))
  entry <- .agent_consumer(published, consumer)
  if (length(entry$events) > 0L || !is.null(published$stopped) || wait <= 0) {
    return(.agent_answer(200L, list(events = .agent_take_events(published, entry))))
  }
  promises::promise(function(resolve, reject) {
    entry$waiter_count <- entry$waiter_count + 1L
    key <- sprintf("w%08d", entry$waiter_count)
    waiter <- new.env(parent = emptyenv())
    waiter$resolve <- resolve
    waiter$cancel <- later::later(function() .agent_resolve_waiter(published, entry, key), wait)
    assign(key, waiter, envir = entry$waiters)
  })
}

.agent_command_answer <- function(viewer, published, req) {
  body <- tryCatch(req$rook.input$read(), error = function(error) raw(0))
  refused <- .agent_answer(400L, list(error = "A command is one JSON object."))
  if (length(body) == 0L || length(body) > .agent_command_body_limit_bytes) return(refused)
  text <- .utf8_text(body)
  if (is.na(text)) return(refused)
  command <- tryCatch(list(value = jsonlite::parse_json(text, simplifyVector = FALSE)), error = function(error) NULL)
  if (is.null(command)) return(refused)
  tryCatch({
    .agent_run_command(viewer, published, command$value)
    .agent_answer(202L, structure(list(), names = character()))
  },
  cellucid_viewer_unavailable = function(error) .agent_answer(409L, list(error = conditionMessage(error))),
  error = function(error) .agent_answer(400L, list(error = conditionMessage(error))))
}

.agent_run_command <- function(viewer, published, command) {
  type <- .json_field(command, "type")
  if (!.is_single_text(type)) .abort("A command is one JSON object with its type.")
  request_id <- command[["requestId"]]
  if (!.is_single_text(request_id)) .abort("A command names its requestId.")
  if (identical(type, "agent_list")) {
    connection <- command[["connection"]]
    .require_agent_connection(connection, type)
    name <- connection$id
    .agent_consumer(published, name)
    failed <- function(reason) .agent_push(published, name, list(type = "agent_failed", requestId = request_id, reason = reason))
    tryCatch(
      .agent_port_list(viewer, connection, .agent_list_seconds, function(kind, value) {
        if (identical(kind, "tools")) {
          .agent_push(published, name, list(type = "agent_tools", requestId = request_id, tools = value))
        } else {
          failed(value)
        }
      }),
      cellucid_viewer_unavailable = function(error) failed(conditionMessage(error))
    )
    return(invisible(NULL))
  }
  if (identical(type, "agent_call")) {
    connection <- command[["connection"]]
    .require_agent_connection(connection, type)
    name <- connection$id
    entry <- .agent_consumer(published, name)
    entry$requests <- union(entry$requests, request_id)
    listener <- function(kind, data) {
      event <- switch(kind,
        result = list(type = "agent_result", requestId = request_id, result = data$result),
        asking = list(type = "agent_asking", requestId = request_id, message = data$message, detail = data$detail),
        answered = list(type = "agent_answered", requestId = request_id, allowed = data$allowed),
        stopped = list(type = "viewer_stopped", reason = data$reason),
        list(type = "agent_failed", requestId = request_id, reason = data$reason)
      )
      .agent_push(published, name, event)
    }
    tryCatch(
      .agent_port_start_call(viewer, request_id, command[["tool"]], command[["args"]], command[["agent"]], connection, listener),
      error = function(error) {
        entry$requests <- setdiff(entry$requests, request_id)
        stop(error)
      }
    )
    return(invisible(NULL))
  }
  if (identical(type, "agent_cancel")) return(.agent_port_cancel(viewer, request_id))
  if (identical(type, "agent_answer")) return(.agent_port_answer(viewer, request_id, isTRUE(command[["allowed"]])))
  .abort("'", type, "' is no agent command.")
}

# =============================================================================
# A viewer in another process
# =============================================================================

# A curl handle for one request to a published viewer's endpoint: never
# through a proxy (the endpoint is on this computer), with its token.
.agent_handle <- function(record, method, body = NULL, timeout = 10) {
  handle <- curl::new_handle(customrequest = method, noproxy = "*", timeout = timeout, connecttimeout = min(timeout, 5))
  headers <- list(Authorization = paste0("Bearer ", record[["token"]]))
  if (!is.null(body)) {
    curl::handle_setopt(handle, postfields = charToRaw(enc2utf8(.js_json(body))))
    headers$`Content-Type` <- "application/json"
  }
  curl::handle_setheaders(handle, .list = headers)
  handle
}

.agent_payload <- function(content) {
  text <- .utf8_text(content)
  value <- if (is.na(text) || !nzchar(text)) NULL else tryCatch(jsonlite::parse_json(text, simplifyVector = FALSE), error = function(error) NULL)
  if (is.list(value) && !is.null(names(value))) value else list()
}

# One request to a published viewer's endpoint, waited for: list(status,
# payload). An error when its process cannot be reached. It is made through
# curl's multi interface while R's event loop runs (later::run_now()), so the
# viewers this process serves -- and their agents' requests -- keep being
# answered while it waits; the server's standard input waits in its pipe
# meanwhile, for at most `timeout` seconds. A viewer of this R session is
# reached through its port, never over HTTP (.agent_link_of()).
.agent_request <- function(record, method, path, body = NULL, timeout = 10) {
  result <- NULL
  failure <- NULL
  pool <- curl::new_pool()
  curl::curl_fetch_multi(paste0(record[["endpoint"]], path), done = function(response) result <<- response,
                         fail = function(message) failure <<- message, pool = pool,
                         handle = .agent_handle(record, method, body, timeout))
  while (is.null(result) && is.null(failure)) {
    curl::multi_run(timeout = 0.01, poll = TRUE, pool = pool)
    later::run_now(0)
  }
  if (!is.null(failure)) .abort(failure)
  list(status = result$status_code, payload = .agent_payload(result$content))
}

# A published viewer as its process describes it now; from its record, and
# not connected, when its process does not answer.
.agent_describe_record <- function(record) {
  viewer <- .agent_own_viewer(record)
  if (!is.null(viewer)) return(.agent_port_describe(viewer))
  answer <- tryCatch(.agent_request(record, "GET", "", timeout = 5), error = function(error) NULL)
  described <- if (!is.null(answer) && identical(answer$status, 200L)) .json_field(answer$payload, "viewer")
  if (is.list(described) && !is.null(names(described))) return(described)
  list(viewer = record[["viewer"]], runtime = record[["runtime"]], dataset = record[["dataset"]], url = record[["url"]],
       connected = FALSE)
}

# The viewer of this R session a record names, or NULL.
.agent_own_viewer <- function(record) {
  if (!identical(as.integer(record[["pid"]]), Sys.getpid())) return(NULL)
  viewer <- get0(record[["viewer"]], envir = .viewer_registry, inherits = FALSE)
  if (is.environment(viewer) && isTRUE(viewer$running)) viewer else NULL
}

# --- One interface to either ------------------------------------------------------

# A link to one viewer an agent connection drives (cellucid-python's
# ViewerLink): describe(), connected(), list_tools(connection, timeout,
# done), start_call(request_id, tool, args, agent, connection, listener),
# cancel(request_id), answer(request_id, allowed), add_tools_listener(fn)
# (returns its remover), stopped() (why it stopped, or NULL), pump() (its
# share of the server's loop) and close().
.agent_link_of <- function(record, consumer) {
  viewer <- .agent_own_viewer(record)
  if (!is.null(viewer)) .local_viewer_link(viewer) else .remote_viewer_link(record, consumer)
}

# A viewer of this process, through its port.
.local_viewer_link <- function(viewer) {
  link <- new.env(parent = emptyenv())
  link$viewer_id <- viewer$viewer_id
  link$describe <- function() .agent_port_describe(viewer)
  link$connected <- function() .agent_port_connected(viewer)
  link$list_tools <- function(connection, timeout, done) .agent_port_list(viewer, connection, timeout, done)
  link$start_call <- function(request_id, tool, args, agent, connection, listener) {
    .agent_port_start_call(viewer, request_id, tool, args, agent, connection, listener)
  }
  link$cancel <- function(request_id) .agent_port_cancel(viewer, request_id)
  link$answer <- function(request_id, allowed) .agent_port_answer(viewer, request_id, allowed)
  link$add_tools_listener <- function(listener) .agent_port_add_tools_listener(viewer, listener)
  link$stopped <- function() .agent_port(viewer)$stopped
  link$pump <- function() invisible(NULL)
  link$close <- function() invisible(NULL)
  link
}

# A viewer another process published, driven over its endpoint. Its events
# are asked for without holding R (curl's multi interface, run by pump());
# an answer with no events is asked again after a short pause (an endpoint
# may answer at once), one with events at once.
.remote_viewer_link <- function(record, consumer) {
  link <- new.env(parent = emptyenv())
  state <- new.env(parent = emptyenv())
  state$calls <- new.env(parent = emptyenv())
  state$lists <- new.env(parent = emptyenv())
  state$tools_listeners <- list()
  state$listener_count <- 0L
  state$stopped <- NULL
  state$closed <- FALSE
  state$pool <- curl::new_pool()
  state$asking <- FALSE
  state$next_ask <- 0
  state$failures <- 0L
  link$viewer_id <- record[["viewer"]]
  link$record <- record
  link$state <- state

  command <- function(value) {
    if (!is.null(state$stopped)) .viewer_unavailable("The Cellucid viewer stopped (", state$stopped, ").")
    answer <- tryCatch(.agent_request(record, "POST", "/command", value), error = function(error) {
      .viewer_unavailable("The Cellucid viewer's process cannot be reached (", conditionMessage(error), ").")
    })
    if (identical(answer$status, 202L)) return(invisible(NULL))
    reason <- .json_field(answer$payload, "error")
    if (!.is_single_text(reason)) reason <- paste0("its process answered ", answer$status)
    if (identical(answer$status, 409L)) .viewer_unavailable(reason)
    if (answer$status %in% c(401L, 404L)) .viewer_unavailable("The Cellucid viewer is gone (", reason, ").")
    .abort(reason)
  }

  end <- function(reason) {
    state$stopped <- reason
    ids <- ls(state$calls, all.names = TRUE)
    listeners <- mget(ids, envir = state$calls)
    rm(list = ids, envir = state$calls)
    for (listener in listeners) .agent_hear(listener, "stopped", list(reason = reason))
    for (id in ls(state$lists, all.names = TRUE)) settle(id, "failed", paste0("The Cellucid viewer stopped (", reason, ")."))
    state$closed <- TRUE
  }

  settle <- function(id, kind, value) {
    slot <- get0(id, envir = state$lists, inherits = FALSE)
    if (is.null(slot)) return(invisible(NULL))
    rm(list = id, envir = state$lists)
    slot$cancel()
    tryCatch(slot$done(kind, value), error = function(error) {
      .warn("A Cellucid agent tool list's listener failed: ", conditionMessage(error))
    })
  }

  route <- function(event) {
    kind <- .json_field(event, "type")
    if (identical(kind, "agent_tools_changed")) {
      for (listener in state$tools_listeners) try(listener(), silent = TRUE)
      return(invisible(NULL))
    }
    if (identical(kind, "viewer_stopped")) {
      reason <- .json_field(event, "reason")
      end(if (.is_single_text(reason)) reason else "the viewer was stopped")
      return(invisible(NULL))
    }
    id <- .json_field(event, "requestId")
    if (!.is_single_text(id)) return(invisible(NULL))
    if (exists(id, envir = state$lists, inherits = FALSE)) {
      if (identical(kind, "agent_tools")) {
        settle(id, "tools", .json_field(event, "tools") %||% list())
      } else {
        reason <- .json_field(event, "reason")
        settle(id, "failed", if (.is_single_text(reason)) reason else "the viewer did not list its tools")
      }
      return(invisible(NULL))
    }
    listener <- get0(id, envir = state$calls, inherits = FALSE)
    if (is.null(listener)) return(invisible(NULL))
    if (kind %in% c("agent_result", "agent_failed")) rm(list = id, envir = state$calls)
    switch(kind,
      agent_result = .agent_hear(listener, "result", list(result = .json_field(event, "result"))),
      agent_failed = .agent_hear(listener, "failed", list(reason = paste(.json_field(event, "reason")))),
      agent_asking = .agent_hear(listener, "asking", list(message = .json_field(event, "message"), detail = .json_field(event, "detail"))),
      agent_answered = .agent_hear(listener, "answered", list(allowed = isTRUE(.json_field(event, "allowed"))))
    )
    invisible(NULL)
  }

  failed <- function() {
    state$failures <- state$failures + 1L
    if (!isTRUE(.native_process_alive(record[["pid"]]))) return(end("its process ended"))
    # An R viewer's endpoint is its own server, which stops with the viewer:
    # its record gone, the viewer was stopped.
    path <- file.path(tryCatch(.agents_dir(), error = function(error) ""), paste0(record[["viewer"]], ".json"))
    if (!file.exists(path)) return(end("the viewer was stopped"))
    # Back off while its process is busy or restarting its endpoint.
    state$next_ask <- .agent_clock() + min(5, 0.2 * state$failures)
  }

  ask <- function() {
    state$asking <- TRUE
    query <- paste0("?consumer=", utils::URLencode(consumer, reserved = TRUE), "&wait=", .agent_event_wait_limit_seconds)
    curl::curl_fetch_multi(
      paste0(record[["endpoint"]], "/events", query),
      done = function(response) {
        state$asking <- FALSE
        if (state$closed) return(invisible(NULL))
        payload <- .agent_payload(response$content)
        events <- .json_field(payload, "events")
        if (identical(response$status_code, 200L) && is.list(events) && is.null(names(events))) {
          state$failures <- 0L
          state$next_ask <- if (length(events) > 0L) 0 else .agent_clock() + 0.1
          for (event in events) {
            if (state$closed) break
            route(event)
          }
        } else if (response$status_code %in% c(401L, 404L)) {
          end("the viewer was stopped")
        } else {
          failed()
        }
      },
      fail = function(message) {
        state$asking <- FALSE
        if (!state$closed) failed()
      },
      pool = state$pool,
      handle = .agent_handle(record, "GET", timeout = .agent_event_wait_limit_seconds + 10)
    )
  }

  link$pump <- function() {
    if (state$closed) return(invisible(NULL))
    if (!state$asking && .agent_clock() >= state$next_ask) ask()
    curl::multi_run(timeout = 0, poll = TRUE, pool = state$pool)
    invisible(NULL)
  }
  link$describe <- function() {
    description <- .agent_describe_record(record)
    if (!is.null(state$stopped)) description$connected <- FALSE
    description
  }
  link$connected <- function() is.null(state$stopped) && isTRUE(link$describe()$connected)
  link$list_tools <- function(connection, timeout, done) {
    id <- paste0(connection$id, "-list-", .random_hex(4L))
    slot <- new.env(parent = emptyenv())
    slot$done <- done
    slot$cancel <- later::later(function() settle(id, "failed", "The Cellucid viewer did not list its tools in time."), timeout)
    assign(id, slot, envir = state$lists)
    tryCatch(command(list(type = "agent_list", requestId = id, connection = connection)), error = function(error) {
      if (exists(id, envir = state$lists, inherits = FALSE)) rm(list = id, envir = state$lists)
      slot$cancel()
      stop(error)
    })
    invisible(id)
  }
  link$start_call <- function(request_id, tool, args, agent, connection, listener) {
    assign(request_id, listener, envir = state$calls)
    tryCatch(
      command(list(type = "agent_call", requestId = request_id, tool = tool, args = args, agent = agent, connection = connection)),
      error = function(error) {
        if (exists(request_id, envir = state$calls, inherits = FALSE)) rm(list = request_id, envir = state$calls)
        stop(error)
      }
    )
    invisible(NULL)
  }
  link$cancel <- function(request_id) {
    waiting <- exists(request_id, envir = state$calls, inherits = FALSE)
    if (waiting) {
      rm(list = request_id, envir = state$calls)
      tryCatch(command(list(type = "agent_cancel", requestId = request_id)), error = function(error) NULL)
    }
    invisible(NULL)
  }
  link$answer <- function(request_id, allowed) {
    tryCatch(command(list(type = "agent_answer", requestId = request_id, allowed = isTRUE(allowed))), error = function(error) NULL)
    invisible(NULL)
  }
  link$add_tools_listener <- function(listener) {
    state$listener_count <- state$listener_count + 1L
    key <- sprintf("l%08d", state$listener_count)
    state$tools_listeners[[key]] <- listener
    function() {
      state$tools_listeners[[key]] <- NULL
      invisible(NULL)
    }
  }
  link$stopped <- function() state$stopped
  # Stop reading this viewer's events. A call still waiting can get no answer
  # through this link now: it fails with the page-gone sentence (whether it
  # was done is not known).
  link$close <- function() {
    state$closed <- TRUE
    ids <- ls(state$calls, all.names = TRUE)
    listeners <- mget(ids, envir = state$calls)
    rm(list = ids, envir = state$calls)
    for (listener in listeners) .agent_hear(listener, "failed", list(reason = .agent_page_gone))
    for (id in ls(state$lists, all.names = TRUE)) settle(id, "failed", .agent_page_gone)
    invisible(NULL)
  }
  link
}
