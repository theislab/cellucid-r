# cellucid_mcp(): a local MCP server that lets an agent app on this computer
# -- Claude Code, Claude Desktop, Cursor -- work in a Cellucid viewer. The
# help page is hand-written: man/cellucid_mcp.Rd.
#
# It speaks MCP (JSON-RPC 2.0, one message per line) on this process's
# standard input and output and drives a viewer the cellucid R or Python
# package serves on this computer, with no relay and no deployed site in
# between: cellucid_mcp(open = <file>) serves the file and opens its viewer
# itself, and stops it when the app ends the server; cellucid_mcp() drives a
# viewer opened elsewhere (cellucid_view() in an R session, show_anndata() in
# Python), found among those published on this computer (R/agent-link.R).
#
# cellucid-python's cellucid mcp is built on the official MCP SDK
# (cellucid/mcp.py, cellucid/_mcp_core.py); R has no official one. The
# mcptools package on CRAN serves R functions as tools defined in R (through
# ellmer), where these tools are the web app's own -- generated into
# inst/mcp_tools.json, and listed by the open page itself, the same whatever
# the person's settings -- and every call is run by the page, so the protocol is
# written here, the part this server needs: the handshake of protocol
# versions 2024-11-05 to 2025-11-25 (initialize, notifications/initialized),
# ping, tools/list, tools/call, notifications/cancelled, and the server's own
# requests elicitation/create and notifications/progress; and the per-request protocol 2026-07-28 (its version and the
# client's capabilities in each request's _meta, server/discover), whose
# results say resultType "complete" and carry the server's info, as the web
# app's relay answers it (relay/mcp.js). A per-request client is not asked
# in the app (an elicitation is a request of the handshake's session): the
# viewer's confirm card asks the person.
#
# R runs one thread: standard input is read without blocking (src/stdio.c:
# poll() on Linux and macOS, PeekNamedPipe() for an app's pipe on Windows)
# while R's event loop runs the viewer this process serves (its httpuv
# server, the page's channel) and the requests to viewers of other processes
# (curl's multi interface), so a call waiting for the person never stops the
# server from reading the app's next message: an answer to its question, a
# cancellation, another call. Nothing but protocol messages reaches standard
# output: R's own output is sunk to standard error while the server runs, and
# the messages are written past it, by the output's descriptor.

.mcp_tools_schema <- 1L
.mcp_server_name <- "cellucid"
.mcp_server_title <- "Cellucid"

# How long a page may take to list its tools.
.mcp_list_seconds <- 30

# The protocol versions of the initialize handshake this server speaks, the
# newest first (the answer to a client that asks for another), and of the
# per-request protocol.
.mcp_handshake_versions <- c("2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05")
.mcp_modern_versions <- "2026-07-28"
.mcp_meta_version <- "io.modelcontextprotocol/protocolVersion"
.mcp_meta_capabilities <- "io.modelcontextprotocol/clientCapabilities"

# How long a per-request client may keep server/discover and the tool list
# (relay/mcp.js's TOOL_LIST_TTL_MS).
.mcp_ttl_ms <- 60000L

.mcp_errors <- list(parse = -32700L, invalid_request = -32600L, method_not_found = -32601L,
                    invalid_params = -32602L, internal = -32603L, unsupported_version = -32022L)

# What a call is answered with when no viewer can take it, and when the
# viewer `viewer =` names is not found (the shared fixture's sentences).
.mcp_no_viewer <- paste0(
  "No Cellucid viewer is open on this computer. Ask the person to open their dataset in a viewer ",
  "(cellucid mcp --open <file>, or cellucid serve <file>; cellucid.show_anndata(...) in Python; ",
  "cellucid::cellucid_view(...) in R), ",
  "then try again; cellucid_status says which viewers are found."
)
.mcp_not_found <- "No Cellucid viewer on this computer is {choice}. cellucid_status lists the viewers found."

# What the person is asked in the agent's app: allow (accept) or deny (decline).
.mcp_elicit_note <- " (You can also answer in the Cellucid viewer.)"

.mcp_cache <- new.env(parent = emptyenv())

# The generated document (inst/mcp_tools.json, the web connector's tool list
# for this package's own action catalog, written by the web repository's
# scripts/generate-action-catalog.mjs): the tools, the instructions, the
# status tool and the app names, checked against the catalog it was made
# from; each tool's JSON text kept.
.mcp_tools_document <- function() {
  document <- .mcp_cache$document
  if (!is.null(document)) return(document)
  path <- system.file("mcp_tools.json", package = "cellucid")
  catalog <- system.file("action_catalog.json", package = "cellucid")
  if (!nzchar(path) || !nzchar(catalog)) {
    .abort("This cellucid installation lacks its agent tools: inst/mcp_tools.json or inst/action_catalog.json is missing.")
  }
  document <- jsonlite::parse_json(.utf8_text(readBin(path, "raw", file.size(path))), simplifyVector = FALSE)
  if (!is.list(document) || !identical(document[["schema"]], .mcp_tools_schema)) {
    .abort("mcp_tools.json is not of the schema this cellucid reads")
  }
  if (!identical(document[["catalogDigest"]], .native_sha256(readBin(catalog, "raw", file.size(catalog))))) {
    .abort("mcp_tools.json was made from another action catalog than this cellucid's: reinstall cellucid")
  }
  document$tools_json <- vapply(document$tools, .js_json, "", USE.NAMES = FALSE)
  document$status_json <- .js_json(document$statusTool)
  .mcp_cache$document <- document
  document
}

# Text on one line: runs of white space as one space, none at either end.
.mcp_one_line <- function(text) {
  text <- gsub("(*UCP)^\\s+|\\s+$", "", enc2utf8(text), perl = TRUE)
  gsub("(*UCP)\\s+", " ", text, perl = TRUE)
}

# The agent an app's clientInfo names, as the viewer shows it: its title when
# it gives one, else the name Cellucid knows it by (claude-code is Claude
# Code), else its name as given, else "Local agent". On one line, whole.
.mcp_agent_name <- function(name, title = NULL) {
  document <- .mcp_tools_document()
  known <- if (.is_single_text(name) && nzchar(name)) document$agentNames[[name]]
  for (candidate in list(title, known, name)) {
    if (.is_single_text(candidate)) {
      text <- .mcp_one_line(candidate)
      if (nzchar(text)) return(text)
    }
  }
  document$agentFallback
}

# JSON as a result's data, as the Cellucid tab writes every result
# (relay-protocol.js resultDataText): between the data tags of mcp_tools.json
# (resultData$tags), every "<" written as the six characters \u003c, so a
# name written by others in it is data an agent never takes for instructions
# and nothing in it can close the tags.
.mcp_result_data <- function(json) {
  tags <- .mcp_tools_document()$resultData$tags
  paste0(tags$open, gsub("<", "\\u003c", enc2utf8(json), fixed = TRUE), tags$close)
}

# A name written by others marked as data where the sentence is built, as the
# tab marks every name it quotes (targets/names.js quoteName): the data tags
# around the name as one JSON string, every "<" in it written \u003c, so no
# name can end its quotes, its line or its mark.
.mcp_quote_name <- function(name) {
  tags <- .mcp_tools_document()$resultData$tags
  paste0(tags$open, gsub("<", "\\u003c", .js_strings(name), fixed = TRUE), tags$close)
}

# A name .mcp_quote_name marked: the tags around one JSON string with no raw
# "<", newline or unescaped quote in it (targets/names.js MARKED_NAME).
.mcp_marked_name <- "<cellucid-data>(\"(?:[^\"\\\\\n<]|\\\\.)*\")</cellucid-data>"

# A sentence of this package's as an agent reads it, as the tab sends its
# refusals (mcp-content.js refusalText): the names in it were marked where it
# was built (.mcp_quote_name); every other "<" is escaped, so nothing else can
# open or close a mark. Nothing is marked here.
.mcp_refusal_text <- function(text) {
  text <- enc2utf8(text)
  escape <- function(part) gsub("<", "\\u003c", part, fixed = TRUE)
  marks <- gregexpr(.mcp_marked_name, text, perl = TRUE)[[1L]]
  if (marks[[1L]] == -1L) return(escape(text))
  sizes <- attr(marks, "match.length")
  out <- character()
  at <- 1L
  for (index in seq_along(marks)) {
    out <- c(out, escape(substr(text, at, marks[[index]] - 1L)), substr(text, marks[[index]], marks[[index]] + sizes[[index]] - 1L))
    at <- marks[[index]] + sizes[[index]]
  }
  paste(c(out, escape(substr(text, at, nchar(text)))), collapse = "")
}

# An agent's tool result; an error is a sentence whose quoted names are marked as data.
.mcp_tool_result <- function(text, error) {
  if (isTRUE(error)) text <- .mcp_refusal_text(text)
  .js_json(list(content = list(list(type = "text", text = text)), isError = error))
}

# =============================================================================
# One agent connection: which viewer it drives, and what it answers
# =============================================================================

# One agent connection (one cellucid_mcp() process and its app). `choice`
# names the viewer to drive among those published on this computer (its
# viewer id, dataset id or dataset name), or NULL for the most recently
# opened one whose page is open (else the most recent); `discover` lists the
# published viewers (tests give their own).
.mcp_session <- function(choice = NULL, discover = .agent_discover) {
  session <- new.env(parent = emptyenv())
  session$connection_id <- paste0("local-", .random_hex(6L))
  session$client_name <- NULL
  session$agent <- .mcp_tools_document()$agentFallback
  session$opened <- NULL
  session$opening <- NULL
  session$choice <- choice
  session$discover <- discover
  session$link <- NULL
  session$remote_viewer <- NULL
  # Why the agents folder was refused, when it was (its sentence).
  session$folder_error <- NULL
  session$number <- 0L
  session
}

# The app this connection serves (its clientInfo).
.mcp_adopt_client <- function(session, name, title = NULL) {
  client <- if (.is_single_text(name)) .mcp_one_line(name)
  session$client_name <- if (!is.null(client) && nzchar(client)) client
  session$agent <- .mcp_agent_name(name, title)
  invisible(session)
}

# The connection as agent commands carry it ({id, client}).
.mcp_connection <- function(session) list(id = session$connection_id, client = session$client_name)

.mcp_drop_link <- function(session) {
  if (!is.null(session$link) && !identical(session$link, session$opened)) session$link$close()
  session$link <- NULL
  session$remote_viewer <- NULL
  invisible(NULL)
}

# The viewer this process serves itself, once it started (open =); `opening`
# says it is still starting (its file), or why it failed.
.mcp_set_opened <- function(session, link, opening = NULL) {
  session$opened <- link
  session$opening <- opening
  session$link <- link
  invisible(session)
}

.mcp_matches <- function(record, choice) {
  dataset <- record[["dataset"]]
  any(vapply(list(record[["viewer"]], .json_field(dataset, "id"), .json_field(dataset, "name")),
             function(value) identical(value, choice), TRUE))
}

# The newest viewer whose page is open; else the newest.
.mcp_choose <- function(records) {
  if (length(records) == 0L) return(NULL)
  for (record in records) if (isTRUE(.agent_describe_record(record)$connected)) return(record)
  records[[1L]]
}

# The published viewers, or none with the sentence of a refused agents folder
# kept for what the agent reads.
.mcp_discover <- function(session) {
  tryCatch({
    records <- session$discover()
    session$folder_error <- NULL
    records
  }, cellucid_agents_folder = function(error) {
    session$folder_error <- conditionMessage(error)
    list()
  })
}

# The viewer this connection drives now: the one it opened, else one
# published on this computer, found again when the one it had stopped, or,
# when none was named, when its page is not open and another viewer's is (a
# link let go fails its waiting calls with the page-gone sentence).
.mcp_link <- function(session) {
  if (!is.null(session$opened)) return(session$opened)
  current <- session$link
  if (!is.null(current) && is.null(current$stopped()) && (!is.null(session$choice) || isTRUE(current$connected()))) {
    return(current)
  }
  records <- .mcp_discover(session)
  if (!is.null(session$choice)) records <- Filter(function(record) .mcp_matches(record, session$choice), records)
  chosen <- .mcp_choose(records)
  if (is.null(chosen)) {
    .mcp_drop_link(session)
    return(NULL)
  }
  if (identical(session$remote_viewer, chosen[["viewer"]]) && !is.null(session$link)) return(session$link)
  .mcp_drop_link(session)
  link <- .agent_link_of(chosen, session$connection_id)
  session$link <- link
  session$remote_viewer <- chosen[["viewer"]]
  link
}

# `quote`: how the viewer asked for is named (marked as data in a refusal;
# inside the status's data already, quoted as JSON).
.mcp_unavailable_sentence <- function(session, quote = .mcp_quote_name) {
  if (!is.null(session$opening)) return(session$opening)
  if (!is.null(session$folder_error)) return(session$folder_error)
  if (!is.null(session$choice)) return(sub("{choice}", quote(session$choice), .mcp_not_found, fixed = TRUE))
  .mcp_no_viewer
}

# cellucid_status: the viewer driven, every viewer found, and what to tell the person.
.mcp_status <- function(session) {
  link <- .mcp_link(session)
  driven <- if (!is.null(link)) link$describe()
  found <- lapply(.mcp_discover(session), function(record) {
    list(
      viewer = record[["viewer"]], runtime = record[["runtime"]], dataset = record[["dataset"]], url = record[["url"]],
      driven = !is.null(driven) && identical(driven[["viewer"]], record[["viewer"]])
    )
  })
  advice <- if (is.null(driven)) {
    .mcp_unavailable_sentence(session, quote = .js_strings)
  } else if (!isTRUE(driven[["connected"]])) {
    .not_open_sentence(driven)
  }
  list(
    agent = session$agent,
    connected = !is.null(driven) && isTRUE(driven[["connected"]]),
    viewer = driven,
    viewers = found,
    advice = advice
  )
}

.mcp_next_request_id <- function(session) {
  session$number <- session$number + 1L
  paste0(session$connection_id, "-", session$number)
}

# Send one call to the viewer this connection drives; a
# cellucid_viewer_unavailable error with the sentence an agent reads when
# none can take it.
.mcp_start_call <- function(session, request_id, tool, args, listener) {
  link <- .mcp_link(session)
  if (is.null(link)) .viewer_unavailable(.mcp_unavailable_sentence(session))
  link$start_call(request_id, tool, args, session$agent, .mcp_connection(session), listener)
  link
}

# The tools, handed to `done` as JSON texts: the open page's list, else the
# generated one (the same tools when the page is of this package's Cellucid
# version), and the status tool. Neither depends on the person's settings
# (the page tells the agent those after its tool results), and the server
# never says the list changed: an agent app lists it once, and a changed list
# makes it rewrite its cached prompt.
.mcp_list_tools <- function(session, done) {
  document <- .mcp_tools_document()
  generated <- function() done(c(document$tools_json, document$status_json))
  link <- .mcp_link(session)
  if (is.null(link) || !isTRUE(link$connected())) return(generated())
  tryCatch(
    link$list_tools(.mcp_connection(session), .mcp_list_seconds, function(kind, value) {
      if (identical(kind, "tools")) {
        done(c(vapply(value, .js_json, "", USE.NAMES = FALSE), document$status_json))
      } else {
        generated()
      }
    }),
    error = function(error) generated()
  )
  invisible(NULL)
}

.mcp_session_close <- function(session) {
  .mcp_drop_link(session)
  invisible(NULL)
}

# =============================================================================
# The protocol: JSON-RPC 2.0 messages, one per line
# =============================================================================

.mcp_server <- function(session, write) {
  server <- new.env(parent = emptyenv())
  server$session <- session
  server$write <- write
  server$closed <- FALSE
  server$initialized <- FALSE
  server$answered_initialize <- FALSE
  server$protocol <- NULL
  server$elicitation <- FALSE
  server$calls <- new.env(parent = emptyenv())
  server$elicitations <- new.env(parent = emptyenv())
  server$out_id <- 0L
  server
}

.mcp_send <- function(server, text) {
  if (server$closed) return(invisible(FALSE))
  sent <- tryCatch(isTRUE(server$write(text)), error = function(error) FALSE)
  # The app closed its end: nothing more can be said.
  if (!sent) server$closed <- TRUE
  invisible(sent)
}

.mcp_server_info <- function() {
  list(name = .mcp_server_name, title = .mcp_server_title, version = as.character(utils::packageVersion("cellucid")))
}

.mcp_capabilities <- function() list(tools = list(listChanged = FALSE))

# A result (its JSON text) as the request's protocol answers it: a
# per-request client's says it is complete and names the server.
.mcp_complete <- function(result, modern) {
  if (!modern) return(result)
  inner <- substr(result, 2L, nchar(result) - 1L)
  meta <- paste0("\"_meta\":", .js_json(list(`io.modelcontextprotocol/serverInfo` = .mcp_server_info())))
  paste0("{\"resultType\":\"complete\",", if (nzchar(inner)) paste0(inner, ","), meta, "}")
}

.mcp_answer <- function(server, id, result, modern = FALSE) {
  .mcp_send(server, paste0("{\"jsonrpc\":\"2.0\",\"id\":", id, ",\"result\":", .mcp_complete(result, modern), "}"))
}

.mcp_fail <- function(server, id, code, message, data = NULL) {
  error <- list(code = code, message = message)
  if (!is.null(data)) error$data <- data
  .mcp_send(server, paste0("{\"jsonrpc\":\"2.0\",\"id\":", id %||% "null", ",\"error\":", .js_json(error), "}"))
}

.mcp_notify <- function(server, method, params) {
  .mcp_send(server, .js_json(list(jsonrpc = "2.0", method = method, params = params)))
}

.mcp_failure <- function(code, message, data = NULL) {
  error <- structure(class = c("cellucid_mcp_failure", "error", "condition"),
                     list(message = message, code = code, data = data, call = NULL))
  stop(error)
}

# One line from the app.
.mcp_handle_line <- function(server, bytes) {
  if (length(bytes) > 0L && bytes[[length(bytes)]] == as.raw(13L)) bytes <- bytes[-length(bytes)]
  if (length(bytes) == 0L) return(invisible(NULL))
  text <- .utf8_text(bytes)
  if (!is.na(text) && !grepl("\\S", text, perl = TRUE)) return(invisible(NULL))
  message <- if (is.na(text)) NULL else tryCatch(list(value = jsonlite::parse_json(text, simplifyVector = FALSE)),
                                                 error = function(error) NULL)
  if (is.null(message)) return(.mcp_fail(server, NULL, .mcp_errors$parse, "Parse error"))
  .mcp_handle_message(server, message$value)
}

.mcp_handle_message <- function(server, message) {
  if (!is.list(message) || is.null(names(message))) {
    return(.mcp_fail(server, NULL, .mcp_errors$invalid_request, "Invalid Request: one JSON-RPC message, an object, per line"))
  }
  method <- message[["method"]]
  has_id <- "id" %in% names(message)
  id <- message[["id"]]
  valid_id <- (.is_single_text(id)) || (is.numeric(id) && length(id) == 1L && !is.na(id) && id == floor(id))
  id_text <- if (valid_id) .js_json(id)
  if (!identical(message[["jsonrpc"]], "2.0")) {
    return(.mcp_fail(server, id_text, .mcp_errors$invalid_request, "Invalid Request: jsonrpc must be \"2.0\""))
  }
  if (!.is_single_text(method)) {
    # An answer to a request of this server's (an elicitation).
    if (valid_id && ("result" %in% names(message) || "error" %in% names(message))) {
      return(.mcp_heard_answer(server, id_text, message[["result"]]))
    }
    return(.mcp_fail(server, id_text, .mcp_errors$invalid_request, "Invalid Request: a request names its method"))
  }
  if (!has_id) return(.mcp_handle_notification(server, method, message[["params"]]))
  if (!valid_id) return(.mcp_fail(server, NULL, .mcp_errors$invalid_request, "Invalid Request: a request id is a string or an integer"))
  tryCatch(
    .mcp_handle_request(server, id_text, method, message),
    cellucid_mcp_failure = function(failure) .mcp_fail(server, id_text, failure$code, conditionMessage(failure), failure$data),
    error = function(error) .mcp_fail(server, id_text, .mcp_errors$internal, conditionMessage(error))
  )
}

.mcp_handle_notification <- function(server, method, params) {
  if (identical(method, "notifications/initialized")) {
    server$initialized <- TRUE
  } else if (identical(method, "notifications/cancelled")) {
    # The app gave up on a call: it is stopped in the viewer (a change rolled
    # back, a question withdrawn), and it is not answered.
    id <- .json_field(params, "requestId")
    if (.is_single_text(id) || (is.numeric(id) && length(id) == 1L && !is.na(id))) .mcp_cancel_call(server, .js_json(id))
  }
  invisible(NULL)
}

# Which protocol a request speaks: the per-request one when its _meta names
# a version (checked as relay/mcp.js checks it), else the handshake's.
.mcp_is_modern <- function(method, params) {
  meta <- .json_field(params, "_meta")
  version <- .json_field(meta, .mcp_meta_version)
  if (identical(method, "initialize") || is.null(version)) return(FALSE)
  if (!.is_single_text(version)) {
    .mcp_failure(.mcp_errors$invalid_params, paste0("The request has no _meta \"", .mcp_meta_version, "\""))
  }
  capabilities <- .json_field(meta, .mcp_meta_capabilities)
  if (!is.list(capabilities) || is.null(names(capabilities))) {
    .mcp_failure(.mcp_errors$invalid_params, paste0("The request has no _meta \"", .mcp_meta_capabilities, "\""))
  }
  if (!version %in% .mcp_modern_versions) {
    .mcp_failure(.mcp_errors$unsupported_version, "Unsupported protocol version",
                 list(supported = I(c(.mcp_modern_versions, .mcp_handshake_versions)), requested = version))
  }
  TRUE
}

.mcp_handle_request <- function(server, id, method, message) {
  params <- message[["params"]]
  if (!is.null(params) && !(is.list(params) && !is.null(names(params)))) {
    .mcp_failure(.mcp_errors$invalid_params, "params must be an object")
  }
  modern <- .mcp_is_modern(method, params)
  session <- server$session
  document <- .mcp_tools_document()
  switch(method,
    initialize = {
      requested <- .json_field(params, "protocolVersion")
      version <- if (.is_single_text(requested) && requested %in% .mcp_handshake_versions) requested else .mcp_handshake_versions[[1L]]
      server$protocol <- version
      # A form elicitation goes to an app that declared elicitation, and not
      # for URLs alone (an empty declaration means forms, as before modes).
      elicitation <- .json_field(.json_field(params, "capabilities"), "elicitation")
      server$elicitation <- !is.null(elicitation) &&
        (!is.null(.json_field(elicitation, "form")) || is.null(.json_field(elicitation, "url")))
      info <- .json_field(params, "clientInfo")
      .mcp_adopt_client(session, .json_field(info, "name"), .json_field(info, "title"))
      .mcp_answer(server, id, .js_json(list(
        protocolVersion = version, capabilities = .mcp_capabilities(), serverInfo = .mcp_server_info(),
        instructions = document$instructions
      )))
      server$answered_initialize <- TRUE
    },
    ping = .mcp_answer(server, id, "{}", modern),
    `server/discover` = {
      if (!modern) .mcp_failure(.mcp_errors$method_not_found, "server/discover is part of protocol 2026-07-28")
      .mcp_answer(server, id, .js_json(list(
        supportedVersions = I(c(.mcp_modern_versions, .mcp_handshake_versions)), capabilities = .mcp_capabilities(),
        instructions = document$instructions, ttlMs = .mcp_ttl_ms, cacheScope = "public"
      )), modern)
    },
    `tools/list` = {
      .mcp_list_tools(session, function(tools) {
        extra <- if (modern) paste0(",\"ttlMs\":", .mcp_ttl_ms, ",\"cacheScope\":\"private\"") else ""
        .mcp_answer(server, id, paste0("{\"tools\":[", paste(tools, collapse = ","), "]", extra, "}"), modern)
      })
    },
    `tools/call` = .mcp_call_tool(server, id, params, modern),
    .mcp_failure(.mcp_errors$method_not_found, paste0("Method not found: ", method))
  )
  invisible(NULL)
}

# --- A call ---------------------------------------------------------------------

.mcp_call_tool <- function(server, id, params, modern) {
  session <- server$session
  name <- .json_field(params, "name")
  if (!.is_single_text(name)) .mcp_failure(.mcp_errors$invalid_params, "tools/call names its tool (params.name)")
  args <- .json_field(params, "arguments")
  if (is.null(args)) args <- structure(list(), names = character())
  if (!is.list(args) || (length(args) > 0L && is.null(names(args)))) {
    .mcp_failure(.mcp_errors$invalid_params, "tools/call arguments must be an object")
  }
  if (is.null(names(args))) args <- structure(list(), names = character())
  document <- .mcp_tools_document()
  if (identical(name, document$statusTool$name)) {
    return(.mcp_answer(server, id, .mcp_tool_result(.mcp_result_data(.js_json(.mcp_status(session))), FALSE), modern))
  }
  call <- new.env(parent = emptyenv())
  call$id <- id
  call$modern <- modern
  call$progress <- .json_field(.json_field(params, "_meta"), "progressToken")
  call$request_id <- .mcp_next_request_id(session)
  call$asking <- NULL
  call$finished <- FALSE
  call$link <- NULL
  # A call's progress notifications count 1, 2, ... (they must grow).
  call$step <- 0L
  finish <- function(result) {
    if (call$finished) return(invisible(NULL))
    call$finished <- TRUE
    .mcp_forget_call(server, call)
    .mcp_answer(server, id, result, modern)
  }
  listener <- function(kind, data) {
    if (call$finished) return(invisible(NULL))
    switch(kind,
      asking = {
        message <- data$message
        if (!is.null(call$progress)) {
          call$step <- call$step + 1L
          .mcp_notify(server, "notifications/progress", list(progressToken = call$progress, progress = call$step, message = message))
        }
        if (is.null(call$asking) && server$elicitation && !modern) .mcp_elicit(server, call, message)
      },
      answered = .mcp_forget_elicitation(server, call),
      result = finish(.js_json(data$result)),
      # The page refused it, or went before it answered: the sentence says which.
      failed = finish(.mcp_tool_result(data$reason, TRUE)),
      stopped = finish(.mcp_tool_result(.agent_sentence(.agent_stopped, reason = data$reason), TRUE))
    )
    invisible(NULL)
  }
  assign(id, call, envir = server$calls)
  call$link <- tryCatch(
    .mcp_start_call(session, call$request_id, name, args, listener),
    error = function(error) {
      finish(.mcp_tool_result(conditionMessage(error), TRUE))
      NULL
    }
  )
  invisible(NULL)
}

.mcp_forget_call <- function(server, call) {
  if (exists(call$id, envir = server$calls, inherits = FALSE)) rm(list = call$id, envir = server$calls)
  .mcp_forget_elicitation(server, call)
}

# The question asked in the app is withdrawn there (the page answered it,
# or the call ended): notifications/cancelled for the elicitation request,
# and its late answer is not heard.
.mcp_forget_elicitation <- function(server, call) {
  asking <- call$asking
  call$asking <- NULL
  if (!is.null(asking) && exists(asking, envir = server$elicitations, inherits = FALSE)) {
    rm(list = asking, envir = server$elicitations)
    .mcp_send(server, paste0("{\"jsonrpc\":\"2.0\",\"method\":\"notifications/cancelled\",\"params\":{\"requestId\":", asking,
                             ",\"reason\":\"The question was answered in the Cellucid viewer, or the call ended.\"}}"))
  }
  invisible(NULL)
}

# The person's question, asked in the agent's app too (an MCP form
# elicitation): Accept allows, Decline denies; dismissing it leaves the
# viewer's card to answer. The model never answers it: the app shows it.
.mcp_elicit <- function(server, call, message) {
  server$out_id <- server$out_id + 1L
  id <- .js_json(server$out_id)
  params <- list(message = paste0(message, .mcp_elicit_note),
                 requestedSchema = list(type = "object", properties = structure(list(), names = character())))
  # The form mode is named from 2025-11-25 on; before it, every elicitation is one.
  if (identical(server$protocol, "2025-11-25")) params <- c(list(mode = "form"), params)
  call$asking <- id
  assign(id, call, envir = server$elicitations)
  .mcp_send(server, paste0("{\"jsonrpc\":\"2.0\",\"id\":", id, ",\"method\":\"elicitation/create\",\"params\":", .js_json(params), "}"))
}

.mcp_heard_answer <- function(server, id, result) {
  call <- get0(id, envir = server$elicitations, inherits = FALSE)
  if (is.null(call)) return(invisible(NULL))
  rm(list = id, envir = server$elicitations)
  if (call$finished || !identical(call$asking, id)) return(invisible(NULL))
  # Answered in the app: nothing to withdraw.
  call$asking <- NULL
  action <- .json_field(result, "action")
  if (.is_single_text(action) && action %in% c("accept", "decline") && !is.null(call$link)) {
    try(call$link$answer(call$request_id, identical(action, "accept")), silent = TRUE)
  }
  invisible(NULL)
}

.mcp_cancel_call <- function(server, id) {
  call <- get0(id, envir = server$calls, inherits = FALSE)
  if (is.null(call) || call$finished) return(invisible(NULL))
  call$finished <- TRUE
  .mcp_forget_call(server, call)
  if (!is.null(call$link)) try(call$link$cancel(call$request_id), silent = TRUE)
  invisible(NULL)
}

# =============================================================================
# cellucid_mcp()
# =============================================================================

# Serve `path` in a viewer of this process (open =) and open its page; NULL,
# said to the agent and on standard error, when it cannot open.
.mcp_open_viewer <- function(session, path, port, browser, web_source_url, web_cache_dir) {
  viewer <- tryCatch(
    .start_viewer(
      path.expand(path), list(), label = path, port = port %||% .default_port, port_given = !is.null(port),
      host = .default_host, allowed_hosts = NULL, serve_web_ui = TRUE, web_cache_dir = web_cache_dir,
      web_source_url = web_source_url, cache_bytes = NULL, quiet = FALSE,
      # The viewer says what opened it: this server, which the agent app runs.
      runtime = "r-mcp"
    ),
    error = function(error) error
  )
  if (inherits(viewer, "error")) {
    sentence <- paste0("The Cellucid viewer of ", path, " could not open: ", conditionMessage(viewer))
    .inform(sentence)
    .mcp_set_opened(session, NULL, opening = sentence)
    return(NULL)
  }
  .inform("Cellucid viewer: ", viewer$viewer_url)
  if (browser) {
    opened <- tryCatch({
      utils::browseURL(viewer$viewer_url)
      TRUE
    }, error = function(error) FALSE)
    if (!opened) .inform("Open the address above in a browser: no browser could be opened from here.")
  }
  .mcp_set_opened(session, .local_viewer_link(viewer))
  viewer
}

cellucid_mcp <- function(
    open = NULL,
    viewer = NULL,
    port = NULL,
    browser = TRUE,
    web_source_url = "https://www.cellucid.com",
    web_cache_dir = NULL
) {
  # Standard output is the protocol's from the first line: what R printed
  # before goes out now, and whatever it prints while the server runs goes to
  # standard error; the messages are written past both, by the output's
  # descriptor.
  flush(stdout())
  sinks <- sink.number()
  sink(stderr())
  on.exit(while (sink.number() > sinks) sink(), add = TRUE)
  if (!is.null(open)) open <- .validate_required_string(open, "open")
  if (!is.null(viewer)) viewer <- .validate_required_string(viewer, "viewer")
  if (!is.null(open) && !is.null(viewer)) .abort("Give open or viewer, not both: open drives the viewer it opens.")
  if (!is.null(port)) port <- .validate_port(port)
  browser <- .validate_flag(browser, "browser")
  if (isTRUE(.native_stdin_console())) {
    .abort(
      "cellucid_mcp() talks to an agent app over a pipe, and this is a console: add it to the app's ",
      "settings (Rscript -e \"cellucid::cellucid_mcp()\"), which starts it."
    )
  }
  .mcp_tools_document()
  session <- .mcp_session(choice = viewer)
  server <- .mcp_server(session, function(text) .native_stdout_write(charToRaw(paste0(enc2utf8(text), "\n"))))
  served <- NULL
  # Run before the sink is removed: what stopping the viewer prints still
  # goes to standard error.
  on.exit({
    .mcp_session_close(session)
    if (!is.null(served)) try(cellucid_stop(served), silent = TRUE)
  }, add = TRUE, after = FALSE)
  # A viewer to open is opened once the app's handshake is answered (or after
  # two seconds without one): its server, and at first use the web app's
  # download, take a while, and the app waits for its first answer.
  opening <- !is.null(open)
  if (opening) {
    .mcp_set_opened(session, NULL, opening = paste0("The Cellucid viewer of ", open, " is still opening; try again in a moment."))
  } else if (is.null(.mcp_link(session))) {
    .inform("cellucid_mcp(): ", .mcp_status(session)$advice, " The server waits; its tools work once one is open.")
  }
  started <- .agent_clock()
  buffer <- raw(0)
  repeat {
    chunk <- .native_stdin_read(0.02)
    ended <- !is.null(chunk) && length(chunk) == 0L
    if (length(chunk) > 0L) {
      buffer <- c(buffer, chunk)
      while (!server$closed && length(ends <- which(buffer == as.raw(10L))) > 0L) {
        line <- buffer[seq_len(ends[[1L]] - 1L)]
        buffer <- buffer[-seq_len(ends[[1L]])]
        .mcp_handle_line(server, line)
      }
    }
    if (ended) {
      # The input ended: a last line without its newline is read too.
      if (length(buffer) > 0L && !server$closed) .mcp_handle_line(server, buffer)
      break
    }
    if (server$closed) break
    if (opening && (server$answered_initialize || .agent_clock() - started > 2)) {
      opening <- FALSE
      served <- .mcp_open_viewer(session, open, port, browser, web_source_url, web_cache_dir)
    }
    later::run_now(0)
    link <- session$link
    if (!is.null(link)) link$pump()
  }
  invisible(NULL)
}
