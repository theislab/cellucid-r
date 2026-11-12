# cellucid_mcp(): the local MCP server for agent apps on this computer.
#
# - inst/mcp_tools.json is the web connector's tool list for this package's
#   own inst/action_catalog.json (written by the web repository's
#   node scripts/generate-action-catalog.mjs, which pins the same digest);
#   cellucid-python pins the same file.
# - tests/testthat/fixtures/local_agents.json is shared with cellucid-python
#   byte for byte: the agents folder, the record, the endpoint, the app names
#   and the sentences an agent reads, so either package's server drives
#   either package's viewer and says the same things.
# - The agent commands and events of the viewer channel, a viewer's agent
#   port, its publication on this computer and its endpoint, and the MCP
#   server itself over stdio (initialize, tools/list, tools/call, errors,
#   cancellation, the person asked in the app, a tool list that never
#   changes, a viewer it opens itself) against a fake page; and cellucid-python's server
#   driving an R viewer, where that Python is found.

.mcp_tools_sha256 <- "65b1e61533966676c2fb60f5c307cc00b43b0db952f044dc72fa8839297ffaa0"
.local_agents_sha256 <- "f5a5f265380c9566ccfa5f5bb67165fb95261d4da8f4ac792c179b4b4d75ad00"

agents_fixture <- function() cellucid_read_fixture("local_agents.json")

package_file <- function(name) system.file(name, package = "cellucid")

empty_object <- function() structure(list(), names = character())

# A fresh agents folder for the calling test, the previous one back after it.
local_agents_dir <- function(env = parent.frame()) {
  folder <- normalizePath(tempfile("cellucid-agents-"), winslash = "/", mustWork = FALSE)
  previous <- Sys.getenv("CELLUCID_AGENTS_DIR")
  Sys.setenv(CELLUCID_AGENTS_DIR = folder)
  do.call(base::on.exit, list(bquote(Sys.setenv(CELLUCID_AGENTS_DIR = .(previous))), TRUE, FALSE), envir = env)
  folder
}

# A served viewer of two cells named "Agents", no web build, stopped when the
# calling test ends.
agents_viewer <- function(env = parent.frame()) {
  viewer <- cellucid_serve(
    NULL, embeddings = list(X_umap = matrix(c(0, 0, 1, 1), 2, 2, byrow = TRUE)),
    dataset_id = "agents", dataset_name = "Agents", port = 0L, serve_web_ui = FALSE, quiet = TRUE
  )
  do.call(base::on.exit, list(bquote(cellucid_stop(.(viewer))), TRUE, FALSE), envir = env)
  viewer
}

# Run R's event loop, and every pump, until `done()` holds.
pump_until <- function(done, pumps = list(), timeout = 30, what = "the condition") {
  deadline <- Sys.time() + timeout
  while (!done()) {
    if (Sys.time() > deadline) stop("timed out waiting for ", what, call. = FALSE)
    later::run_now(0.01)
    for (pump in pumps) pump()
  }
  invisible(TRUE)
}

# --- A page on the viewer channel that answers agent commands as the web app does ---

# One whole frame from the socket's bytes, without waiting: list(opcode,
# payload), or NULL.
ws_take_frame <- function(client) {
  client$buffer <- c(client$buffer, readBin(client$connection, "raw", 65536L))
  buffer <- client$buffer
  if (length(buffer) < 2L) return(NULL)
  n <- as.integer(buffer[[2L]]) %% 128L
  at <- 3L
  if (n == 126L) {
    if (length(buffer) < 4L) return(NULL)
    n <- as.integer(buffer[[3L]]) * 256L + as.integer(buffer[[4L]])
    at <- 5L
  } else if (n == 127L) {
    if (length(buffer) < 10L) return(NULL)
    n <- sum(as.numeric(buffer[3:10]) * 256^(7:0))
    at <- 11L
  }
  if (length(buffer) < at + n - 1L) return(NULL)
  payload <- if (n > 0) buffer[at:(at + n - 1L)] else raw(0)
  client$buffer <- if (length(buffer) > at + n - 1L) buffer[(at + n):length(buffer)] else raw(0)
  list(opcode = as.integer(buffer[[1L]]) %% 16L, payload = payload)
}

agent_json <- function(value) as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = NA))

# A binary agent event: header length, its JSON header, the body's UTF-8 JSON.
agent_binary_event <- function(type, request_id, value) {
  header <- charToRaw(agent_json(list(type = type, requestId = request_id)))
  c(cellucid:::.encode_uint32(length(header)), header, charToRaw(enc2utf8(agent_json(value))))
}

text_result <- function(text, error = FALSE) list(content = list(list(type = "text", text = text)), isError = error)

# The status result's JSON: its data, between the data tags (names written by others are marked as data).
status_of <- function(text) {
  tags <- cellucid:::.mcp_tools_document()$resultData$tags
  stopifnot(startsWith(text, tags$open), endsWith(text, tags$close))
  jsonlite::fromJSON(substr(text, nchar(tags$open) + 1L, nchar(text) - nchar(tags$close)), simplifyVector = FALSE)
}

fake_page <- function(viewer, tools = NULL) {
  page <- new.env()
  page$client <- ws_welcomed(viewer)
  page$tools <- tools %||% list(list(
    name = "color_by", title = "Color", description = "From the page.", inputSchema = list(type = "object"),
    annotations = empty_object()
  ))
  page$commands <- list()
  page$ask <- character()
  page$hold <- character()
  # Tools whose question the person answers in the viewer's own card (the test does, with page_card_answer()).
  page$card <- character()
  page$open <- TRUE
  page
}

page_result <- function(page, request_id, result) {
  ws_send_binary(page$client, agent_binary_event("agent_result", request_id, result))
}

# The person answers in the viewer's card.
page_card_answer <- function(page, request_id, allowed) {
  ws_send_json(page$client, list(type = "agent_answered", requestId = request_id, allowed = allowed))
  page_result(page, request_id, text_result(if (allowed) "saved from the card" else "declined in the card", !allowed))
}

# Read and answer every command the page holds now.
page_pump <- function(page) {
  if (!page$open) return(invisible(NULL))
  while (!is.null(frame <- ws_take_frame(page$client))) {
    if (frame$opcode != 1L) next
    command <- jsonlite::fromJSON(rawToChar(frame$payload), simplifyVector = FALSE)
    page$commands[[length(page$commands) + 1L]] <- command
    id <- command$requestId
    switch(command$type,
      agent_list = ws_send_binary(page$client, agent_binary_event("agent_tools", id, page$tools)),
      agent_call = {
        if (command$tool %in% page$card) {
          for (step in c("asks to save a copy.", "asks again to save a copy.")) {
            ws_send_json(page$client, list(type = "agent_asking", requestId = id, message = paste("Claude Code (local)", step), detail = NULL))
          }
        } else if (command$tool %in% page$ask) {
          ws_send_json(page$client, list(
            type = "agent_asking", requestId = id,
            message = "Claude Code (local) asks to save a copy. It writes a file.", detail = NULL
          ))
        } else if (!command$tool %in% page$hold) {
          page_result(page, id, text_result(paste0("ran ", command$tool, " as ", command$agent)))
        }
      },
      agent_answer = {
        ws_send_json(page$client, list(type = "agent_answered", requestId = id, allowed = command$allowed))
        page_result(page, id, text_result(if (isTRUE(command$allowed)) "saved" else "The person declined: save a copy.",
                                          !isTRUE(command$allowed)))
      },
      agent_cancel = page_result(page, id, text_result("Stopped before it finished: the agent stopped waiting", TRUE))
    )
  }
  invisible(NULL)
}

# The next command of `kind` the page received (taken).
page_next <- function(page, kind, pumps = list(), timeout = 30) {
  found <- NULL
  pump_until(function() {
    page_pump(page)
    at <- which(vapply(page$commands, function(command) identical(command$type, kind), TRUE))
    if (length(at) == 0L) return(FALSE)
    found <<- page$commands[[at[[1L]]]]
    page$commands[[at[[1L]]]] <- NULL
    TRUE
  }, pumps, timeout, paste0("the page to receive ", kind))
  found
}

# The page goes as a browser's does: its close frame, then its socket let go.
# (Any test with a server child needs the frame: ws_close() in helper-viewer.R.)
page_close <- function(page) {
  page$open <- FALSE
  ws_send_close(page$client)
  ws_close(page$client)
}

# The page is lost with no word (its browser crashed): only for a test that
# started no child while the page was open (ws_close() in helper-viewer.R).
page_lost <- function(page) {
  page$open <- FALSE
  ws_close(page$client)
}

listen <- function() {
  heard <- new.env()
  heard$events <- list()
  heard$listener <- function(kind, data) heard$events[[length(heard$events) + 1L]] <- list(kind = kind, data = data)
  heard
}

heard_next <- function(heard, pumps = list(), timeout = 30) {
  pump_until(function() length(heard$events) > 0L, pumps, timeout, "an agent event")
  event <- heard$events[[1L]]
  heard$events[[1L]] <- NULL
  event
}

# One HTTP request to this process's own server (curl's multi interface,
# R's event loop answering it): list(status, json).
agent_http <- function(url, method = "GET", headers = list(), body = NULL, timeout = 30) {
  handle <- curl::new_handle(customrequest = method, noproxy = "*")
  if (!is.null(body)) curl::handle_setopt(handle, postfields = charToRaw(enc2utf8(body)))
  if (length(headers) > 0L) curl::handle_setheaders(handle, .list = headers)
  result <- NULL
  failure <- NULL
  pool <- curl::new_pool()
  curl::curl_fetch_multi(url, done = function(response) result <<- response,
                         fail = function(message) failure <<- message, pool = pool, handle = handle)
  pump_until(function() {
    curl::multi_run(timeout = 0, poll = TRUE, pool = pool)
    !is.null(result) || !is.null(failure)
  }, timeout = timeout, what = url)
  if (!is.null(failure)) stop(failure)
  text <- rawToChar(result$content)
  list(status = result$status_code, json = if (nzchar(text)) jsonlite::fromJSON(text, simplifyVector = FALSE))
}

# =============================================================================
# The generated tools and the shared fixture
# =============================================================================

test_that("mcp_tools.json is the list the web app's source makes now", {
  path <- package_file("mcp_tools.json")
  expect_identical(cellucid_rendered_by_web("scripts/generate-action-catalog.mjs", "renderMcpTools"),
                   readChar(path, file.size(path), useBytes = TRUE),
                   info = "Run node scripts/generate-action-catalog.mjs in the web app.")
})

test_that("mcp_tools.json is the web connector's list for this package's catalog", {
  path <- package_file("mcp_tools.json")
  expect_identical(cellucid_file_sha256(path), .mcp_tools_sha256,
                   info = "regenerate it: node scripts/generate-action-catalog.mjs in the web repository; pin the digest it prints")
  document <- cellucid:::.mcp_tools_document()
  expect_identical(document$catalogDigest, cellucid_file_sha256(package_file("action_catalog.json")))
  catalog <- jsonlite::read_json(package_file("action_catalog.json"), simplifyVector = FALSE)
  # One tool per action, named, titled and annotated from it, in the catalog's order.
  expect_identical(vapply(document$tools, `[[`, "", "name"), vapply(catalog$actions, `[[`, "", "tool"))
  for (index in seq_along(document$tools)) {
    tool <- document$tools[[index]]
    action <- catalog$actions[[index]]
    expect_identical(tool$title, action$title)
    expect_identical(tool$annotations$title, action$title)
    expect_true(startsWith(tool$description, action$summary), info = tool$name)
    expect_identical(tool$annotations$readOnlyHint, identical(action$effect, "read"), info = tool$name)
  }
  expect_identical(document$statusTool$name, agents_fixture()$server$statusTool)
  expect_false(document$statusTool$name %in% vapply(document$tools, `[[`, "", "name"))
  expect_match(document$instructions, "a Cellucid viewer on this computer", fixed = TRUE)
  names <- c(vapply(document$tools, `[[`, "", "name"), document$statusTool$name)
  used <- regmatches(document$instructions, gregexpr("`[a-z][a-z0-9_]*`", document$instructions))[[1L]]
  for (name in gsub("`", "", used)) expect_true(name %in% names, info = name)
})

test_that("the tools document is cellucid-python's, byte for byte", {
  expect_identical(cellucid_file_sha256(cellucid_sibling_file("cellucid-python", "src", "cellucid", "mcp_tools.json")), .mcp_tools_sha256)
})

test_that("the local instructions name the R calls that bring an agent's work back, and they exist", {
  # A local agent told the person to run a function that does not exist (live-8 L26-6).
  text <- cellucid:::.mcp_tools_document()$instructions
  for (call in c("plan <- cellucid_apply_session(object, cellucid_request_session(viewer))", "cellucid_apply_plan(plan)")) {
    expect_true(grepl(call, text, fixed = TRUE), info = call)
  }
  for (name in c("cellucid_apply_session", "cellucid_request_session", "cellucid_apply_plan")) {
    expect_true(name %in% getNamespaceExports("cellucid"), info = name)
  }
})

test_that("a tools document made from another catalog is refused by name", {
  saved <- cellucid:::.mcp_tools_document()
  on.exit(assign("document", saved, envir = cellucid:::.mcp_cache))
  rm("document", envir = cellucid:::.mcp_cache)
  local_mocked_bindings(.native_sha256 = function(bytes) "0000", .package = "cellucid")
  expect_error(cellucid:::.mcp_tools_document(), "another action catalog")
})

test_that("the shared local agents fixture is cellucid-python's, byte for byte", {
  expect_identical(cellucid_file_sha256(cellucid_sibling_file("cellucid-python", "tests", "fixtures", "local_agents.json")), .local_agents_sha256)
})

test_that("the shared local agents fixture is this package's behavior", {
  expect_pinned_fixture("local_agents.json", .local_agents_sha256)
  fixture <- agents_fixture()
  expect_identical(fixture$server$name, cellucid:::.mcp_server_name)
  expect_identical(fixture$server$title, cellucid:::.mcp_server_title)
  expect_identical(fixture$record$schema, cellucid:::.agent_record_schema)
  expect_identical(as.character(unlist(fixture$record$fields)), cellucid:::.agent_record_fields)
  endpoint <- fixture$endpoint
  expect_identical(endpoint$route, cellucid:::.agent_route)
  expect_equal(endpoint$eventWaitLimitSeconds, cellucid:::.agent_event_wait_limit_seconds)
  expect_equal(endpoint$consumerIdleSeconds, cellucid:::.agent_consumer_idle_seconds)
  expect_equal(endpoint$commandBodyLimitBytes, cellucid:::.agent_command_body_limit_bytes)
  expect_identical(fixture$agentsDir$override, "CELLUCID_AGENTS_DIR")
  for (case in fixture$agentNames) {
    expect_identical(cellucid:::.mcp_agent_name(case$name, case$title), case$agent, info = paste(case$name, case$title))
  }
  sentences <- fixture$sentences
  expect_identical(sentences$noViewer, cellucid:::.mcp_no_viewer)
  expect_identical(sentences$notFound, cellucid:::.mcp_not_found)
  expect_identical(sentences$notOpen, cellucid:::.not_open_sentence(list(dataset = list(name = "PBMC"), url = "http://127.0.0.1:8765")))
  expect_identical(sentences$notOpenWidget, cellucid:::.not_open_sentence(list(dataset = list(name = "PBMC"), url = NULL)))
  expect_identical(sentences$pageGone, cellucid:::.agent_page_gone)
  expect_identical(sentences$notTaken, cellucid:::.agent_not_taken)
  expect_identical(sentences$stopped, cellucid:::.agent_stopped)
  expect_identical(sentences$folderUnsafe, cellucid:::.agent_folder_unsafe)
  expect_identical(fixture$optOut$variable, cellucid:::.agent_opt_out_variable)
  expect_identical(as.character(unlist(fixture$optOut$values)), cellucid:::.agent_opt_out_values)
  expect_identical(fixture$endpoint$statuses$badLength, 400L)
  session <- cellucid:::.mcp_session(discover = list)
  expect_match(session$connection_id, fixture$server$connectionIdPattern)
  # The protocol's own pattern of a connection id takes it.
  expect_match(session$connection_id, cellucid:::.viewer_agent_connection_id_pattern)
})

test_that("the agents folder follows the user cache and its override", {
  previous <- Sys.getenv("CELLUCID_AGENTS_DIR")
  on.exit(Sys.setenv(CELLUCID_AGENTS_DIR = previous))
  folder <- normalizePath(tempdir(), winslash = "/")
  Sys.setenv(CELLUCID_AGENTS_DIR = folder)
  expect_identical(cellucid:::.agents_dir(), folder)
  Sys.setenv(CELLUCID_AGENTS_DIR = "relative/folder")
  expect_error(cellucid:::.agents_dir(), "absolute")
  # Nothing is published, and nothing is found, where the folder is refused: the agent reads why.
  expect_error(cellucid:::.agent_discover(), "absolute", class = "cellucid_agents_folder")
  expect_identical(cellucid:::.mcp_status(cellucid:::.mcp_session())$advice, "CELLUCID_AGENTS_DIR must be an absolute path")
  Sys.unsetenv("CELLUCID_AGENTS_DIR")
  default <- cellucid:::.agents_dir()
  expect_identical(basename(dirname(default)), "cellucid")
  expect_identical(basename(default), "agents")
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    expect_identical(default, file.path(path.expand("~"), "Library", "Caches", "cellucid", "agents"))
  } else if (.Platform$OS.type != "windows") {
    xdg <- Sys.getenv("XDG_CACHE_HOME")
    base <- if (nzchar(xdg) && startsWith(xdg, "/")) xdg else file.path(path.expand("~"), ".cache")
    expect_identical(default, file.path(base, "cellucid", "agents"))
  }
})

# =============================================================================
# The agent commands and events on the viewer channel
# =============================================================================

test_that("agent commands are checked and serialized as the viewer protocol has them", {
  protocol <- cellucid_read_fixture("viewer_protocol.json")
  for (example in protocol$examples$commands) {
    if (!startsWith(example$type, "agent_")) next
    frame <- cellucid:::.viewer_command_frame(example)
    expect_identical(jsonlite::fromJSON(frame, simplifyVector = FALSE), example, info = example$type)
  }
  frame <- cellucid:::.viewer_command_frame
  connection <- list(id = "local-1", client = NULL)
  expect_error(frame(list(type = "agent_list", requestId = "a", connection = list(id = "has space", client = NULL))), "connection id")
  expect_error(frame(list(type = "agent_list", requestId = "a", connection = list(id = "x"))), "\\{id, client\\}")
  expect_error(frame(list(type = "agent_list", requestId = "a", connection = list(id = "x", client = " padded"))), "exact text")
  expect_error(frame(list(type = "agent_call", requestId = "a", tool = "color_by", args = list(), agent = "A", connection = connection)), "args must be one object")
  expect_error(frame(list(type = "agent_call", requestId = "a", tool = "Color By", args = empty_object(), agent = "A", connection = connection)), "names its tool")
  expect_error(frame(list(type = "agent_call", requestId = "a", tool = "color_by", args = empty_object(), agent = "two\nlines", connection = connection)), "exact text")
  # Names of any length: no cap of this package's own.
  long <- strrep("a", 500)
  expect_match(frame(list(type = "agent_call", requestId = "a", tool = "color_by", args = empty_object(), agent = long,
                          connection = list(id = "x", client = long))), long, fixed = TRUE)
  expect_error(frame(list(type = "agent_answer", requestId = "a", allowed = "yes")), "allowed")
  expect_error(frame(list(type = "agent_cancel", requestId = "")), "requestId")
  # An empty object stays an object.
  expect_match(frame(list(type = "agent_call", requestId = "a", tool = "history_get", args = empty_object(), agent = "A", connection = connection)),
               "\"args\":{}", fixed = TRUE)
})

test_that("agent events are read as the viewer protocol has them", {
  protocol <- cellucid_read_fixture("viewer_protocol.json")
  for (example in protocol$examples$events) {
    if (!startsWith(example$type, "agent_")) next
    parsed <- cellucid:::.parse_viewer_event(agent_json(example))
    expect_identical(c(list(type = parsed$type), parsed$payload), example)
  }
  binary <- function(header, body) {
    bytes <- charToRaw(header)
    c(cellucid:::.encode_uint32(length(bytes)), bytes, body)
  }
  refused <- function(frame, pattern) {
    expect_error(cellucid:::.parse_viewer_binary_event(frame), pattern, class = "cellucid_protocol_error")
  }
  refused(binary("{\"type\":\"agent_result\",\"requestId\":\"a\"}", charToRaw("{\"content\":\"no\",\"isError\":false}")), "MCP tool result")
  refused(binary("{\"type\":\"agent_result\",\"requestId\":\"a\"}", charToRaw("{\"contentX\":[],\"isError\":false}")), "MCP tool result")
  refused(binary("{\"type\":\"agent_tools\",\"requestId\":\"a\"}", as.raw(0xff)), "not UTF-8 JSON")
  refused(binary("{\"type\":\"agent_tools\",\"requestId\":\"a\"}", charToRaw("[{\"title\":\"x\"}]")), "each with its name")
  refused(binary("{\"type\":\"agent_tools\",\"requestId\":\"\"}", charToRaw("[]")), "requestId")
  expect_error(cellucid:::.parse_viewer_event("{\"type\":\"agent_answered\",\"requestId\":\"a\",\"allowed\":\"yes\"}"), "allowed",
               class = "cellucid_protocol_error")
  expect_error(cellucid:::.parse_viewer_event("{\"type\":\"agent_asking\",\"requestId\":\"a\",\"message\":\"\",\"detail\":null}"), "message",
               class = "cellucid_protocol_error")
})

# =============================================================================
# A viewer's agent port, its publication and its endpoint, with a fake page
# =============================================================================

test_that("a viewer is published, and its port carries calls both ways", {
  folder <- local_agents_dir()
  viewer <- agents_viewer()
  expect_false(cellucid:::.agent_port_connected(viewer))
  expect_error(
    cellucid:::.agent_port_start_call(viewer, "r1", "color_by", empty_object(), "Claude Code", list(id = "local-1", client = NULL), function(...) NULL),
    paste0("No page shows the Cellucid viewer of Agents \\(served at http://127\\.0\\.0\\.1:", viewer$port, "\\) now\\."),
    class = "cellucid_viewer_unavailable"
  )
  records <- cellucid:::.agent_discover()
  expect_length(records, 1L)
  record <- records[[1L]]
  expect_identical(sort(names(record)), cellucid:::.agent_record_fields)
  expect_identical(record$viewer, viewer$viewer_id)
  expect_identical(record$dataset, list(id = "agents", name = "Agents"))
  expect_identical(record$pid, Sys.getpid())
  expect_identical(record$runtime, "r")
  # Where it is served, never its page address (which carries the viewer's token).
  expect_identical(record$url, viewer$url)
  expect_false(grepl(viewer$viewer_token, paste(unlist(record), collapse = " "), fixed = TRUE))
  expect_identical(record$host, Sys.info()[["nodename"]])
  expect_identical(record$schema, 2L)
  expect_identical(record$endpoint, paste0("http://127.0.0.1:", viewer$port, "/_cellucid/agent/", viewer$viewer_id))
  expect_match(record$started, "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\\.[0-9]{6}Z$")
  expect_match(record$token, "^[0-9a-f]{64}$")
  path <- file.path(folder, paste0(viewer$viewer_id, ".json"))
  if (.Platform$OS.type != "windows") {
    expect_identical(as.character(file.info(path)$mode), "600")
    expect_identical(as.character(file.info(folder)$mode), "700")
  }
  page <- fake_page(viewer)
  expect_true(cellucid:::.agent_port_connected(viewer))
  connection <- list(id = "local-abc", client = "claude-code")
  listed <- NULL
  cellucid:::.agent_port_list(viewer, connection, 10, function(kind, value) listed <<- list(kind = kind, value = value))
  pump_until(function() !is.null(listed), list(function() page_pump(page)))
  expect_identical(listed$kind, "tools")
  expect_equal(listed$value, page$tools)
  heard <- listen()
  cellucid:::.agent_port_start_call(viewer, "local-abc-1", "color_by", list(field = "x"), "Claude Code", connection, heard$listener)
  command <- page_next(page, "agent_call")
  expect_identical(command, list(type = "agent_call", requestId = "local-abc-1", tool = "color_by", args = list(field = "x"),
                                 agent = "Claude Code", connection = connection))
  expect_identical(heard_next(heard), list(kind = "result", data = list(result = text_result("ran color_by as Claude Code"))))
  # The person is asked: heard, answered from the app, settled.
  page$ask <- "save_copy"
  heard <- listen()
  cellucid:::.agent_port_start_call(viewer, "local-abc-2", "save_copy", empty_object(), "Claude Code", connection, heard$listener)
  asking <- heard_next(heard, list(function() page_pump(page)))
  expect_identical(asking$kind, "asking")
  expect_match(asking$data$message, "^Claude Code \\(local\\) asks")
  expect_null(asking$data$detail)
  cellucid:::.agent_port_answer(viewer, "local-abc-2", FALSE)
  expect_identical(heard_next(heard, list(function() page_pump(page))), list(kind = "answered", data = list(allowed = FALSE)))
  settled <- heard_next(heard, list(function() page_pump(page)))
  expect_identical(settled$kind, "result")
  expect_true(settled$data$result$isError)
  # No agent event reaches the hooks or the viewer's last events.
  expect_null(viewer$last$agent_result)
  expect_null(viewer$last$agent_asking)
  # Stopping the viewer ends what waits, and withdraws it.
  page$hold <- "slow"
  page$commands <- list()
  heard <- listen()
  cellucid:::.agent_port_start_call(viewer, "local-abc-3", "slow", empty_object(), "Claude Code", connection, heard$listener)
  expect_identical(page_next(page, "agent_call")$requestId, "local-abc-3")
  page_close(page)
  cellucid_stop(viewer)
  expect_identical(heard_next(heard), list(kind = "stopped", data = list(reason = "the viewer was stopped")))
  expect_false(file.exists(path))
  expect_identical(cellucid:::.agent_discover(), list())
  expect_error(
    cellucid:::.agent_port_start_call(viewer, "r9", "color_by", empty_object(), "A", connection, function(...) NULL),
    "The Cellucid viewer stopped \\(the viewer was stopped\\)\\.", class = "cellucid_viewer_unavailable"
  )
})

test_that("a refused agent command reaches its caller, never a hook or a warning", {
  local_agents_dir()
  viewer <- agents_viewer()
  hooked <- FALSE
  cellucid_on(viewer, "command_error", function(event) hooked <<- TRUE)
  heard <- listen()
  assign("r1", heard$listener, envir = cellucid:::.agent_port(viewer)$calls)
  expect_no_warning(cellucid:::.deliver_event(viewer, list(type = "command_error", payload = list(
    command = "agent_call", requestId = "r1", reason = "This viewer does not know that host command; reload the page."
  ))))
  expect_identical(heard_next(heard), list(kind = "failed", data = list(reason = paste(
    "The Cellucid viewer page did not take the call: This viewer does not know that host command; reload the page.",
    "Ask the person to reload the viewer's page."
  ))))
  expect_false(hooked)
  expect_null(viewer$last$command_error)
  # A refused tool list is said to whoever asked for it.
  listed <- NULL
  slot <- new.env()
  slot$done <- function(kind, value) listed <<- list(kind, value)
  slot$cancel <- function() NULL
  assign("l1", slot, envir = cellucid:::.agent_port(viewer)$lists)
  cellucid:::.deliver_event(viewer, list(type = "command_error", payload = list(command = "agent_list", requestId = "l1", reason = "Old page.")))
  expect_identical(listed, list("failed", "The Cellucid viewer page did not take the call: Old page. Ask the person to reload the viewer's page."))
})

test_that("records are read only when whole, of this computer, and on its loopback", {
  folder <- local_agents_dir()
  dir.create(folder, recursive = TRUE, mode = "0700")
  host <- Sys.info()[["nodename"]]
  # A process id no system gives (2^31 - 2): a process that ended.
  stale <- list(dataset = list(id = "d", name = "D"), endpoint = "http://127.0.0.1:1/_cellucid/agent/v", host = host,
                pid = 2147483646L, runtime = "python", schema = 2L, started = "2026-01-01T00:00:00.000000Z", token = "t",
                url = NULL, viewer = "v")
  write_record <- function(record, name) writeLines(cellucid:::.js_json(record), file.path(folder, name))
  write_record(stale, "v.json")
  writeLines("{not json", file.path(folder, "broken.json"))
  live <- stale
  live$pid <- Sys.getpid()
  # Every field of its type, its endpoint on this computer's loopback.
  broken <- list(
    c(live[setdiff(names(live), "token")]),
    replace(live, "schema", list(1L)),
    replace(live, "token", list("")),
    replace(live, "pid", list("1")),
    replace(live, "url", list(7L)),
    replace(live, "dataset", list(list(id = "d"))),
    replace(live, "dataset", list(list(id = "", name = "D"))),
    replace(live, "endpoint", list("http://192.168.1.4:1/_cellucid/agent/v")),
    replace(live, "endpoint", list("https://127.0.0.1:1/_cellucid/agent/v")),
    replace(live, "endpoint", list("http://127.0.0.1/_cellucid/agent/v")),
    replace(live, "endpoint", list("http://127.0.0.1.evil.example:1/_cellucid/agent/v"))
  )
  for (index in seq_along(broken)) write_record(broken[[index]], paste0("b", index, ".json"))
  for (index in seq_along(broken)) expect_null(cellucid:::.agent_read_record(file.path(folder, paste0("b", index, ".json"))), info = index)
  # Another computer's record (a home folder shared between machines): skipped, never removed, even with its process gone.
  elsewhere <- replace(stale, "host", list(paste0(host, "-elsewhere")))
  write_record(elsewhere, "x.json")
  ok <- replace(live, "endpoint", list("http://[::1]:8/_cellucid/agent/v"))
  write_record(ok, "ok.json")
  found <- cellucid:::.agent_discover()
  expect_identical(vapply(found, `[[`, "", "endpoint"), "http://[::1]:8/_cellucid/agent/v")
  expect_false(file.exists(file.path(folder, "v.json")))
  expect_true(file.exists(file.path(folder, "x.json")))
  expect_true(file.exists(file.path(folder, "broken.json")))
})

test_that("an agents folder others may write in is refused, said to the agent and on publication", {
  skip_on_os("windows")
  folder <- local_agents_dir()
  dir.create(folder, recursive = TRUE)
  Sys.chmod(folder, "0777", use_umask = FALSE)
  sentence <- sub("{why}", "others may write in it", sub("{folder}", folder, cellucid:::.agent_folder_unsafe, fixed = TRUE), fixed = TRUE)
  expect_error(cellucid:::.agent_discover(), sentence, fixed = TRUE, class = "cellucid_agents_folder")
  expect_warning(viewer <- agents_viewer(), "not published to agents on this computer", fixed = TRUE)
  expect_null(cellucid:::.agent_port(viewer)$published)
  expect_length(list.files(folder), 0L)
  session <- cellucid:::.mcp_session()
  status <- cellucid:::.mcp_status(session)
  expect_identical(status$advice, sentence)
  expect_identical(status$viewers, list())
  expect_error(cellucid:::.mcp_start_call(session, "r", "color_by", empty_object(), function(...) NULL), sentence, fixed = TRUE)
})

test_that("the agents folder's owner is this process's user id, not a user name", {
  skip_on_os("windows")
  folder <- local_agents_dir()
  dir.create(folder, recursive = TRUE, mode = "0700")
  # The process's own uid owns the folder it made, whatever name (if any) it has.
  expect_identical(cellucid:::.agent_process_uid(), file.info(folder, extra_cols = TRUE)$uid)
  expect_identical(cellucid:::.agents_require_private(folder), folder)
  # Another uid's folder is refused, as cellucid-python refuses it (st_uid against os.getuid()).
  local_mocked_bindings(.agent_process_uid = function() file.info(folder, extra_cols = TRUE)$uid + 1L, .package = "cellucid")
  sentence <- sub("{why}", "another user owns it", sub("{folder}", folder, cellucid:::.agent_folder_unsafe, fixed = TRUE), fixed = TRUE)
  expect_error(cellucid:::.agent_discover(), sentence, fixed = TRUE, class = "cellucid_agents_folder")
})

test_that("CELLUCID_AGENTS=off publishes nothing", {
  folder <- local_agents_dir()
  previous <- Sys.getenv("CELLUCID_AGENTS")
  on.exit(Sys.setenv(CELLUCID_AGENTS = previous), add = TRUE)
  for (value in c(" Off ", "0", "FALSE", "no")) {
    Sys.setenv(CELLUCID_AGENTS = value)
    expect_true(cellucid:::.agents_disabled(), info = value)
  }
  Sys.setenv(CELLUCID_AGENTS = "off")
  viewer <- agents_viewer()
  expect_null(cellucid:::.agent_port(viewer)$published)
  expect_false(dir.exists(folder))
  Sys.setenv(CELLUCID_AGENTS = "on")
  expect_false(cellucid:::.agents_disabled())
})

test_that("a page that goes fails the calls it held with the page-gone sentence", {
  local_agents_dir()
  viewer <- agents_viewer()
  page <- fake_page(viewer)
  page$hold <- "slow"
  connection <- list(id = "local-abc", client = NULL)
  heard <- listen()
  cellucid:::.agent_port_start_call(viewer, "local-abc-1", "slow", empty_object(), "A", connection, heard$listener)
  page_next(page, "agent_call")
  # A reload: the page goes, and a new one takes the channel.
  page_close(page)
  expect_identical(heard_next(heard), list(kind = "failed", data = list(reason = cellucid:::.agent_page_gone)))
  page <- fake_page(viewer)
  on.exit(page_close(page), add = TRUE)
  # A page opened elsewhere replaces the one that held the channel: its calls fail too.
  page$hold <- "slow"
  heard <- listen()
  cellucid:::.agent_port_start_call(viewer, "local-abc-2", "slow", empty_object(), "A", connection, heard$listener)
  page_next(page, "agent_call")
  listed <- NULL
  cellucid:::.agent_port_list(viewer, connection, 30, function(kind, value) listed <<- list(kind, value))
  other <- ws_welcomed(viewer)
  on.exit(ws_close(other), add = TRUE)
  expect_identical(heard_next(heard), list(kind = "failed", data = list(reason = cellucid:::.agent_page_gone)))
  expect_identical(listed, list("failed", cellucid:::.agent_page_gone))
})

test_that("a page's close frame is heard with its socket still open, and a page lost without one is heard too", {
  local_agents_dir()
  viewer <- agents_viewer()
  connection <- list(id = "local-abc", client = NULL)
  gone <- list(kind = "failed", data = list(reason = cellucid:::.agent_page_gone))
  # Its socket stays open here, as it stays open on Windows when a child
  # process started meanwhile inherited a handle on it: the frame alone ends
  # the channel, and the server closes the connection in answer (httpuv may
  # close before its answering frame is written, so that frame is not awaited).
  page <- fake_page(viewer)
  page$hold <- "slow"
  heard <- listen()
  cellucid:::.agent_port_start_call(viewer, "local-abc-1", "slow", empty_object(), "A", connection, heard$listener)
  page_next(page, "agent_call")
  ws_send_close(page$client)
  expect_identical(heard_next(heard), gone)
  expect_false(cellucid:::.agent_port_connected(viewer))
  ws_close(page$client)
  # No frame at all: the connection's end is the only word.
  page <- fake_page(viewer)
  page$hold <- "slow"
  heard <- listen()
  cellucid:::.agent_port_start_call(viewer, "local-abc-2", "slow", empty_object(), "A", connection, heard$listener)
  page_next(page, "agent_call")
  page_lost(page)
  expect_identical(heard_next(heard), gone)
  expect_false(cellucid:::.agent_port_connected(viewer))
})

test_that("silent consumers are swept on a clock: their calls stopped in the page, their events let go", {
  local_agents_dir()
  viewer <- agents_viewer()
  page <- fake_page(viewer)
  on.exit(page_close(page), add = TRUE)
  page$hold <- "slow"
  published <- cellucid:::.agent_port(viewer)$published
  expect_true(is.function(published$cancel_sweep))
  record <- cellucid:::.agent_discover()[[1L]]
  good <- list(Authorization = paste0("Bearer ", record$token), `Content-Type` = "application/json")
  sent <- agent_http(paste0(record$endpoint, "/command"), "POST", good, agent_json(list(
    type = "agent_call", requestId = "local-gone-1", tool = "slow", args = empty_object(), agent = "A",
    connection = list(id = "local-gone", client = NULL)
  )))
  expect_identical(sent$status, 202L)
  page_next(page, "agent_call")
  entry <- get("local-gone", envir = published$consumers)
  entry$events <- list(list(type = "agent_failed", requestId = "local-gone-0", reason = "not read"))
  entry$seen <- entry$seen - cellucid:::.agent_consumer_idle_seconds - 1
  cellucid:::.agent_sweep(published)
  expect_identical(page_next(page, "agent_cancel")$requestId, "local-gone-1")
  expect_false(exists("local-gone", envir = published$consumers, inherits = FALSE))
  expect_identical(entry$events, list())
})

test_that("the endpoint refuses browsers, other hosts and wrong tokens, and holds a request for events", {
  local_agents_dir()
  viewer <- agents_viewer()
  record <- cellucid:::.agent_discover()[[1L]]
  endpoint <- record$endpoint
  good <- list(Authorization = paste0("Bearer ", record$token))
  ask <- function(headers = good, path = "", method = "GET", body = NULL) agent_http(paste0(endpoint, path), method, headers, body)
  described <- ask()
  expect_identical(described$status, 200L)
  expect_identical(described$json$schema, 2L)
  expect_identical(described$json$viewer, list(viewer = viewer$viewer_id, runtime = "r", dataset = list(id = "agents", name = "Agents"),
                                               url = viewer$url, connected = FALSE))
  expect_identical(ask(list(Authorization = "Bearer nope"))$status, 401L)
  expect_identical(ask(list())$status, 401L)
  expect_identical(ask(c(good, list(Origin = "https://evil.example")))$status, 403L)
  expect_identical(ask(c(good, list(Host = "evil.example")))$status, 403L)
  expect_identical(ask(c(good, list(Host = paste0("evil.example:", viewer$port))))$status, 403L)
  expect_identical(ask(path = "/nothing")$status, 404L)
  expect_identical(ask(path = "/events")$status, 400L)
  expect_identical(agent_http(sub(viewer$viewer_id, "0000000000000000", endpoint), headers = good)$status, 404L)
  call <- agent_json(list(type = "agent_call", requestId = "r", tool = "color_by", args = empty_object(), agent = "A",
                          connection = list(id = "local-x", client = NULL)))
  closed <- ask(c(good, list(`Content-Type` = "application/json")), "/command", "POST", call)
  expect_identical(closed$status, 409L)
  expect_match(closed$json$error, paste0("No page shows the Cellucid viewer of Agents (served at ", viewer$url, ") now."), fixed = TRUE)
  expect_identical(ask(good, "/command", "POST", "{\"type\":\"nothing\",\"requestId\":\"r\"}")$status, 400L)
  expect_identical(ask(good, "/command", "POST", "not json")$status, 400L)
  expect_identical(ask(good, "/command", "POST", agent_json(list(type = "agent_call", requestId = "r")))$status, 400L)
  # A Content-Length that is not a whole number of bytes: the command is not read -- answered 400, or
  # (httpuv's own parser refuses the request before R sees it) the connection closed with no answer.
  for (length in c("1.5", "1e3", "-4", "twelve")) {
    status <- tryCatch(ask(c(good, list(`Content-Length` = length)), "/command", "POST", "{}")$status, error = function(error) "closed")
    expect_true(status %in% list(400L, "closed"), info = paste(length, status))
  }
  expect_null(cellucid:::.agent_request_refusal(list(REQUEST_METHOD = "GET", HTTP_HOST = "127.0.0.1:1", REMOTE_ADDR = "127.0.0.1")))
  for (length in c("1.5", "1e3", "-4", "twelve", "", "0")) {
    refused <- cellucid:::.agent_request_refusal(list(REQUEST_METHOD = "POST", HTTP_HOST = "127.0.0.1:1", REMOTE_ADDR = "127.0.0.1",
                                                      CONTENT_LENGTH = length))
    expect_identical(refused$status, 400L, info = length)
  }
  # A request for events waits for one (a promise httpuv answers when it settles).
  page <- fake_page(viewer)
  on.exit(page_close(page), add = TRUE)
  # A page taking the channel says nothing to the consumers: its tools are the same whatever its settings.
  expect_identical(ask(path = "/events?consumer=local-x&wait=0")$json, list(events = list()))
  started <- Sys.time()
  expect_identical(ask(path = "/events?consumer=local-x&wait=1")$json, list(events = list()))
  expect_gte(as.numeric(Sys.time() - started, units = "secs"), 0.9)
  held <- NULL
  pool <- curl::new_pool()
  handle <- curl::new_handle(noproxy = "*")
  curl::handle_setheaders(handle, .list = good)
  curl::curl_fetch_multi(paste0(endpoint, "/events?consumer=local-x&wait=20"), done = function(response) held <<- response,
                         pool = pool, handle = handle)
  pump_until(function() {
    curl::multi_run(timeout = 0, poll = TRUE, pool = pool)
    length(ls(get("local-x", envir = cellucid:::.agent_port(viewer)$published$consumers)$waiters)) > 0L
  })
  # A command sent through the endpoint reaches the page; its answer, the
  # consumer's event, settles the request held for one.
  started <- Sys.time()
  listed <- ask(c(good, list(`Content-Type` = "application/json")), "/command", "POST",
                agent_json(list(type = "agent_list", requestId = "l1", connection = list(id = "local-x", client = NULL))))
  expect_identical(listed$status, 202L)
  page_next(page, "agent_list")
  pump_until(function() {
    page_pump(page)
    curl::multi_run(timeout = 0, poll = TRUE, pool = pool)
    !is.null(held)
  })
  expect_lt(as.numeric(Sys.time() - started, units = "secs"), 10)
  events <- jsonlite::fromJSON(rawToChar(held$content), simplifyVector = FALSE)$events
  expect_length(events, 1L)
  expect_identical(events[[1L]]$type, "agent_tools")
  expect_identical(events[[1L]]$requestId, "l1")
  expect_equal(events[[1L]]$tools, page$tools)
})

test_that("the session chooses the newest viewer whose page is open and names its agent", {
  session <- cellucid:::.mcp_session(discover = list)
  expect_null(cellucid:::.mcp_link(session))
  expect_identical(cellucid:::.mcp_status(session), list(agent = "Local agent", connected = FALSE, viewer = NULL, viewers = list(),
                                                         advice = cellucid:::.mcp_no_viewer))
  cellucid:::.mcp_adopt_client(session, "claude-code")
  expect_identical(session$agent, "Claude Code")
  expect_identical(cellucid:::.mcp_connection(session), list(id = session$connection_id, client = "claude-code"))
  tools <- NULL
  cellucid:::.mcp_list_tools(session, function(listed) tools <<- listed)
  names <- vapply(tools, function(text) jsonlite::fromJSON(text, simplifyVector = FALSE)$name, "", USE.NAMES = FALSE)
  document <- cellucid:::.mcp_tools_document()
  expect_identical(names, c(vapply(document$tools, `[[`, "", "name"), "cellucid_status"))
  expect_error(cellucid:::.mcp_start_call(session, "r", "color_by", empty_object(), function(...) NULL),
               "No Cellucid viewer is open on this computer", class = "cellucid_viewer_unavailable")
  chosen <- cellucid:::.mcp_session(choice = "pbmc", discover = list)
  expect_identical(cellucid:::.mcp_status(chosen)$advice, sub("{choice}", "\"pbmc\"", cellucid:::.mcp_not_found, fixed = TRUE))

  local_agents_dir()
  older <- agents_viewer()
  viewer <- agents_viewer()
  page <- fake_page(viewer)
  on.exit(page_close(page), add = TRUE)
  session <- cellucid:::.mcp_session()
  cellucid:::.mcp_adopt_client(session, "cursor-vscode")
  tools <- NULL
  cellucid:::.mcp_list_tools(session, function(listed) tools <<- listed)
  pump_until(function() !is.null(tools), list(function() page_pump(page)))
  expect_equal(lapply(tools, jsonlite::fromJSON, simplifyVector = FALSE), c(page$tools, list(document$statusTool)))
  expect_identical(page_next(page, "agent_list")$connection, list(id = session$connection_id, client = "cursor-vscode"))
  status <- cellucid:::.mcp_status(session)
  expect_true(status$connected)
  expect_identical(status$viewer$viewer, viewer$viewer_id)
  expect_identical(vapply(status$viewers, function(found) found$driven, TRUE)[vapply(status$viewers, `[[`, "", "viewer") == viewer$viewer_id], TRUE)
  expect_length(status$viewers, 2L)
  heard <- listen()
  cellucid:::.mcp_start_call(session, cellucid:::.mcp_next_request_id(session), "color_by", empty_object(), heard$listener)
  expect_identical(page_next(page, "agent_call")$agent, "Cursor")
  expect_identical(heard_next(heard, list(function() page_pump(page)))$kind, "result")
  # Named by its dataset, the older one is chosen.
  chosen <- cellucid:::.mcp_session(choice = older$viewer_id)
  expect_identical(cellucid:::.mcp_status(chosen)$viewer$viewer, older$viewer_id)
  expect_identical(cellucid:::.mcp_status(chosen)$advice, cellucid:::.not_open_sentence(cellucid:::.agent_port_describe(older)))
  cellucid:::.mcp_session_close(session)
})

# =============================================================================
# The MCP server over stdio
# =============================================================================

# The R expression a child R runs to serve cellucid_mcp(): the package as the
# suite loaded it (the source tree under test_local(), the installed package
# under R CMD check).
mcp_expression <- function(arguments = "") {
  path <- getNamespaceInfo("cellucid", "path")
  if (file.exists(file.path(path, "R", "mcp.R"))) {
    # The source tree is loaded by pkgload (test_local(), load_all()), which the child loads it with too.
    if (!isNamespaceLoaded("pkgload")) testthat::skip("the source tree is loaded without pkgload")
    return(sprintf("pkgload::load_all(%s, quiet = TRUE, helpers = FALSE, export_all = FALSE); cellucid::cellucid_mcp(%s)",
                   deparse(path), arguments))
  }
  sprintf("cellucid::cellucid_mcp(%s)", arguments)
}

# An MCP client app: JSON-RPC 2.0 over a child process's stdin and stdout.
# Every wait ends on what it waits for or on the server ending (mcp_next()
# fails at once, with its stderr). Starting and ending the server are waited
# for on those events alone, never on a clock: they are R's start and
# package load, and its teardown, as long as the computer makes them (an
# antivirus starved by CPU load stretched a start past a minute on a
# developer's Mac), and the server answers or ends either way. The other
# waits are bounded by what the test sends and the fake page answers; their
# deadline only turns a hang into a failure.
mcp_client <- function(command = file.path(R.home("bin"), "Rscript"), args = c("--vanilla", "-e", mcp_expression()), env = character()) {
  testthat::skip_if_not_installed("processx")
  client <- new.env()
  # MCP's stdio is UTF-8 whatever the session's native encoding (Windows before
  # R 4.2: a code page), so its output is not read in the session's encoding.
  client$process <- processx::process$new(
    command, args, stdin = "|", stdout = "|", stderr = "|", encoding = "UTF-8",
    env = c("current", R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep), env)
  )
  client$held <- list()
  client$stderr <- character()
  client$noise <- character()
  client$partial <- ""
  client
}

mcp_pump <- function(client) {
  process <- client$process
  out <- process$read_output()
  if (nzchar(out)) {
    text <- paste0(client$partial, out)
    lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
    if (!endsWith(text, "\n")) {
      client$partial <- lines[[length(lines)]]
      lines <- lines[-length(lines)]
    } else {
      client$partial <- ""
    }
    for (line in lines[nzchar(lines)]) {
      message <- tryCatch(jsonlite::parse_json(line, simplifyVector = FALSE), error = function(error) NULL)
      if (is.null(message)) client$noise <- c(client$noise, line) else client$held[[length(client$held) + 1L]] <- message
    }
  }
  errors <- process$read_error_lines()
  if (length(errors) > 0L) client$stderr <- c(client$stderr, errors)
  invisible(NULL)
}

# Bytes: given text, write_input() reads it as native whatever its mark.
mcp_write <- function(client, text) {
  client$process$write_input(charToRaw(enc2utf8(text)))
}

mcp_send <- function(client, message) {
  mcp_write(client, paste0(agent_json(c(list(jsonrpc = "2.0"), message)), "\n"))
}

mcp_next <- function(client, match, pumps = list(), timeout = 60) {
  found <- NULL
  pump_until(function() {
    mcp_pump(client)
    for (index in seq_along(client$held)) {
      if (isTRUE(match(client$held[[index]]))) {
        found <<- client$held[[index]]
        client$held[[index]] <- NULL
        return(TRUE)
      }
    }
    if (!client$process$is_alive()) {
      mcp_pump(client)
      stop("the server ended: ", paste(client$stderr, collapse = "\n"), call. = FALSE)
    }
    FALSE
  }, pumps, timeout, "the server's message")
  found
}

mcp_reply <- function(client, id, pumps = list(), timeout = 60) {
  mcp_next(client, function(message) identical(message$id, id) && is.null(message$method), pumps, timeout)
}

mcp_initialize <- function(client, capabilities = empty_object(), name = "claude-code", version = "2025-06-18", pumps = list()) {
  mcp_send(client, list(id = 1L, method = "initialize", params = list(
    protocolVersion = version, capabilities = capabilities, clientInfo = list(name = name, version = "1")
  )))
  # The server's start: until it answers or ends.
  answer <- mcp_reply(client, 1L, pumps, timeout = Inf)
  mcp_send(client, list(method = "notifications/initialized"))
  answer
}

# Close the server's input; its exit status (it ends by itself).
mcp_close <- function(client, pumps = list()) {
  try(close(client$process$get_input_connection()), silent = TRUE)
  pump_until(function() {
    mcp_pump(client)
    !client$process$is_alive()
  }, pumps, Inf, "the server to end")
  mcp_pump(client)
  client$process$get_exit_status()
}

expect_clean_streams <- function(client) {
  expect_identical(client$noise, character(), info = "standard output carries protocol messages only")
  expect_false(any(startsWith(client$stderr, "{\"jsonrpc\"")))
}

test_that("the server without a viewer lists the generated tools and says what is missing", {
  testthat::skip_on_cran()
  local_agents_dir()
  client <- mcp_client()
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  document <- cellucid:::.mcp_tools_document()
  answer <- mcp_initialize(client)$result
  expect_identical(answer$protocolVersion, "2025-06-18")
  expect_identical(answer$serverInfo[c("name", "title")], list(name = "cellucid", title = "Cellucid"))
  expect_identical(answer$serverInfo$version, as.character(utils::packageVersion("cellucid")))
  # The list never changes while the server runs, and it says so.
  expect_identical(answer$capabilities$tools, list(listChanged = FALSE))
  expect_identical(answer$instructions, document$instructions)
  mcp_send(client, list(id = 2L, method = "tools/list"))
  tools <- mcp_reply(client, 2L)$result$tools
  expect_equal(tools, c(document$tools[seq_along(document$tools)], list(document$statusTool)))
  mcp_send(client, list(id = 3L, method = "tools/call", params = list(name = "color_by", arguments = list(field = list(source = "obs", key = "x")))))
  expect_identical(mcp_reply(client, 3L)$result, text_result(cellucid:::.mcp_refusal_text(cellucid:::.mcp_no_viewer), TRUE))
  mcp_send(client, list(id = 4L, method = "tools/call", params = list(name = "cellucid_status", arguments = empty_object())))
  status_text <- mcp_reply(client, 4L)$result$content[[1L]]$text
  tags <- document$resultData$tags
  expect_true(startsWith(status_text, tags$open) && endsWith(status_text, tags$close))
  status <- jsonlite::fromJSON(substr(status_text, nchar(tags$open) + 1L, nchar(status_text) - nchar(tags$close)), simplifyVector = FALSE)
  expect_identical(status, list(agent = "Claude Code", connected = FALSE, viewer = NULL, viewers = list(), advice = cellucid:::.mcp_no_viewer))
  mcp_send(client, list(id = 5L, method = "no/such/method"))
  expect_identical(mcp_reply(client, 5L)$error$code, -32601L)
  mcp_send(client, list(id = 6L, method = "ping"))
  expect_identical(mcp_reply(client, 6L)$result, empty_object())
  mcp_send(client, list(id = "seven", method = "tools/call", params = list(arguments = empty_object())))
  expect_identical(mcp_reply(client, "seven")$error$code, -32602L)
  mcp_write(client, "{not json\n")
  expect_identical(mcp_next(client, function(message) identical(message$error$code, -32700L))$id, NULL)
  # The per-request protocol: its version and capabilities in each request's _meta.
  meta <- list(`io.modelcontextprotocol/protocolVersion` = "2026-07-28", `io.modelcontextprotocol/clientCapabilities` = empty_object())
  mcp_send(client, list(id = 8L, method = "server/discover", params = list(`_meta` = meta)))
  discovered <- mcp_reply(client, 8L)$result
  expect_identical(discovered$resultType, "complete")
  expect_identical(unlist(discovered$supportedVersions), c("2026-07-28", "2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05"))
  expect_identical(discovered$`_meta`$`io.modelcontextprotocol/serverInfo`$name, "cellucid")
  mcp_send(client, list(id = 9L, method = "tools/list", params = list(`_meta` = meta)))
  listed <- mcp_reply(client, 9L)$result
  expect_identical(listed$resultType, "complete")
  expect_length(listed$tools, length(document$tools) + 1L)
  mcp_send(client, list(id = 10L, method = "ping", params = list(`_meta` = list(`io.modelcontextprotocol/protocolVersion` = "2099-01-01",
                                                                               `io.modelcontextprotocol/clientCapabilities` = empty_object()))))
  refused <- mcp_reply(client, 10L)$error
  expect_identical(refused$code, -32022L)
  expect_identical(refused$data$requested, "2099-01-01")
  mcp_send(client, list(id = 11L, method = "server/discover"))
  expect_identical(mcp_reply(client, 11L)$error$code, -32601L)
  expect_identical(mcp_close(client), 0L)
  expect_true(any(grepl("No Cellucid viewer is open on this computer", client$stderr, fixed = TRUE)))
  expect_clean_streams(client)
})

test_that("an unknown protocol version is answered with the newest, and a last line without its newline is read", {
  testthat::skip_on_cran()
  local_agents_dir()
  client <- mcp_client()
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  expect_identical(mcp_initialize(client, version = "1999-01-01")$result$protocolVersion, "2025-11-25")
  mcp_write(client, "{\"jsonrpc\":\"2.0\",\"id\":2,\"method\":\"ping\"}")
  try(close(client$process$get_input_connection()), silent = TRUE)
  expect_identical(mcp_reply(client, 2L)$result, empty_object())
  expect_identical(mcp_close(client), 0L)
  expect_clean_streams(client)
})

test_that("the tool list is the same whatever the page does, and is never said to change", {
  # The person's settings are told to the agent after its tool results, never
  # in its tool list (L26-9): a list change makes the app rewrite its cached
  # prompt and drop the earlier reasoning. So the list is the generated one
  # before a page opens, the same once a page (whatever its settings) took the
  # channel, reloaded or was replaced, and notifications/tools/list_changed is
  # never sent.
  testthat::skip_on_cran()
  local_agents_dir()
  viewer <- agents_viewer()
  document <- cellucid:::.mcp_tools_document()
  generated <- c(document$tools[seq_along(document$tools)], list(document$statusTool))
  client <- mcp_client()
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  expect_identical(mcp_initialize(client)$result$capabilities$tools, list(listChanged = FALSE))
  mcp_send(client, list(id = 2L, method = "tools/list"))
  expect_equal(mcp_reply(client, 2L)$result$tools, generated)
  page <- fake_page(viewer, tools = document$tools)
  page_close(page)
  page <- fake_page(viewer, tools = document$tools)
  page_close(page)
  page <- fake_page(viewer, tools = document$tools)
  on.exit(page_close(page), add = TRUE)
  pumps <- list(function() page_pump(page))
  mcp_send(client, list(id = 3L, method = "tools/list"))
  expect_equal(mcp_reply(client, 3L, pumps)$result$tools, generated)
  page_next(page, "agent_list", pumps)
  is_list_changed <- function(message) identical(message$method, "notifications/tools/list_changed")
  expect_error(mcp_next(client, is_list_changed, pumps, timeout = 2), "timed out")
  expect_identical(mcp_close(client, pumps), 0L)
  expect_clean_streams(client)
})

test_that("the server drives a published viewer, cancels, and asks in the app", {
  testthat::skip_on_cran()
  local_agents_dir()
  viewer <- agents_viewer()
  page <- fake_page(viewer)
  pumps <- list(function() page_pump(page))
  client <- mcp_client()
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  mcp_initialize(client, capabilities = list(elicitation = empty_object()), name = "claude-code", pumps = pumps)
  mcp_send(client, list(id = 2L, method = "tools/list"))
  expect_identical(vapply(mcp_reply(client, 2L, pumps)$result$tools, `[[`, "", "name"), c("color_by", "cellucid_status"))
  connection <- page_next(page, "agent_list")$connection
  expect_identical(connection$client, "claude-code")
  expect_match(connection$id, agents_fixture()$server$connectionIdPattern)
  # A call runs in the page as the app's agent.
  mcp_send(client, list(id = 3L, method = "tools/call", params = list(name = "color_by", arguments = list(field = list(source = "obs", key = "leiden")))))
  expect_identical(mcp_reply(client, 3L, pumps)$result, text_result("ran color_by as Claude Code"))
  call <- page_next(page, "agent_call")
  expect_identical(call$args, list(field = list(source = "obs", key = "leiden")))
  expect_identical(call$connection, connection)
  expect_identical(call$requestId, paste0(connection$id, "-1"))
  # The app gives up on a call: it is stopped in the page, and not answered.
  page$hold <- "slow"
  mcp_send(client, list(id = 4L, method = "tools/call", params = list(name = "slow", arguments = empty_object(), `_meta` = list(progressToken = "p4"))))
  request_id <- page_next(page, "agent_call", pumps)$requestId
  mcp_send(client, list(method = "notifications/cancelled", params = list(requestId = 4L, reason = "The person pressed Esc")))
  expect_identical(page_next(page, "agent_cancel", pumps)$requestId, request_id)
  # The person is asked: in the viewer's card, and in the app (an elicitation, and the progress the call asked for).
  page$ask <- "save_copy"
  mcp_send(client, list(id = 5L, method = "tools/call", params = list(name = "save_copy", arguments = empty_object(), `_meta` = list(progressToken = 55L))))
  progress <- mcp_next(client, function(message) identical(message$method, "notifications/progress"), pumps)
  expect_identical(progress$params$progressToken, 55L)
  expect_match(progress$params$message, "^Claude Code \\(local\\) asks to save a copy\\.")
  asked <- mcp_next(client, function(message) identical(message$method, "elicitation/create"), pumps)
  expect_identical(asked$params$message, "Claude Code (local) asks to save a copy. It writes a file. (You can also answer in the Cellucid viewer.)")
  expect_identical(asked$params$requestedSchema, list(type = "object", properties = empty_object()))
  expect_null(asked$params$mode)
  mcp_send(client, list(id = asked$id, result = list(action = "accept", content = empty_object())))
  expect_true(page_next(page, "agent_answer", pumps)$allowed)
  expect_identical(mcp_reply(client, 5L, pumps)$result, text_result("saved"))
  # Declined in the app: denied in the page.
  mcp_send(client, list(id = 6L, method = "tools/call", params = list(name = "save_copy", arguments = empty_object())))
  asked <- mcp_next(client, function(message) identical(message$method, "elicitation/create"), pumps)
  mcp_send(client, list(id = asked$id, result = list(action = "decline")))
  expect_false(page_next(page, "agent_answer", pumps)$allowed)
  expect_identical(mcp_reply(client, 6L, pumps)$result, text_result("The person declined: save a copy.", TRUE))
  # The page goes: a call says so in one sentence.
  page_close(page)
  pump_until(function() !cellucid:::.agent_port_connected(viewer))
  mcp_send(client, list(id = 7L, method = "tools/call", params = list(name = "color_by", arguments = empty_object())))
  answer <- mcp_reply(client, 7L)$result
  expect_true(answer$isError)
  expect_identical(answer$content[[1L]]$text, cellucid:::.not_open_sentence(cellucid:::.agent_port_describe(viewer)))
  # Its status names the viewer, not connected, and what to tell the person.
  mcp_send(client, list(id = 8L, method = "tools/call", params = list(name = "cellucid_status", arguments = empty_object())))
  status <- status_of(mcp_reply(client, 8L)$result$content[[1L]]$text)
  expect_false(status$connected)
  expect_identical(status$viewer$viewer, viewer$viewer_id)
  expect_identical(status$viewers, list(list(viewer = viewer$viewer_id, runtime = "r", dataset = list(id = "agents", name = "Agents"),
                                             url = viewer$url, driven = TRUE)))
  expect_identical(status$advice, cellucid:::.not_open_sentence(cellucid:::.agent_port_describe(viewer)))
  expect_identical(mcp_close(client), 0L)
  expect_clean_streams(client)
})

test_that("a page reloaded mid-call fails the call with the page-gone sentence, and the next call reaches the new page", {
  testthat::skip_on_cran()
  local_agents_dir()
  viewer <- agents_viewer()
  page <- fake_page(viewer)
  page$hold <- "slow"
  pumps <- list(function() page_pump(page))
  client <- mcp_client()
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  mcp_initialize(client, pumps = pumps)
  mcp_send(client, list(id = 2L, method = "tools/call", params = list(name = "slow", arguments = empty_object())))
  page_next(page, "agent_call", pumps)
  started <- Sys.time()
  page_close(page)
  expect_identical(mcp_reply(client, 2L)$result, text_result(cellucid:::.mcp_refusal_text(cellucid:::.agent_page_gone), TRUE))
  expect_lt(as.numeric(Sys.time() - started, units = "secs"), 10)
  page <- fake_page(viewer)
  on.exit(page_close(page), add = TRUE)
  mcp_send(client, list(id = 3L, method = "tools/call", params = list(name = "color_by", arguments = empty_object())))
  expect_identical(mcp_reply(client, 3L, list(function() page_pump(page)))$result, text_result("ran color_by as Claude Code"))
  expect_identical(mcp_close(client), 0L)
  expect_clean_streams(client)
})

test_that("a question the person answers in the viewer is withdrawn in the app, and progress grows", {
  testthat::skip_on_cran()
  local_agents_dir()
  viewer <- agents_viewer()
  page <- fake_page(viewer)
  on.exit(page_close(page), add = TRUE)
  page$card <- "save_copy"
  pumps <- list(function() page_pump(page))
  client <- mcp_client()
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  mcp_initialize(client, capabilities = list(elicitation = list(form = empty_object())), version = "2025-11-25", pumps = pumps)
  mcp_send(client, list(id = 2L, method = "tools/call", params = list(name = "save_copy", arguments = empty_object(), `_meta` = list(progressToken = "t2"))))
  call <- page_next(page, "agent_call", pumps)
  asked <- mcp_next(client, function(message) identical(message$method, "elicitation/create"), pumps)
  expect_identical(asked$params$mode, "form")
  first <- mcp_next(client, function(message) identical(message$method, "notifications/progress"), pumps)
  second <- mcp_next(client, function(message) identical(message$method, "notifications/progress"), pumps)
  expect_identical(c(first$params$progress, second$params$progress), c(1L, 2L))
  # One question in the app for the call, however often the card asks.
  page_card_answer(page, call$requestId, TRUE)
  withdrawn <- mcp_next(client, function(message) identical(message$method, "notifications/cancelled"), pumps)
  expect_identical(withdrawn$params$requestId, asked$id)
  expect_identical(mcp_reply(client, 2L, pumps)$result, text_result("saved from the card"))
  expect_length(Filter(function(message) identical(message$method, "elicitation/create"), client$held), 0L)
  # Its late answer is not heard.
  mcp_send(client, list(id = asked$id, result = list(action = "decline")))
  mcp_send(client, list(id = 3L, method = "ping"))
  mcp_reply(client, 3L, pumps)
  expect_false(any(vapply(page$commands, function(command) identical(command$type, "agent_answer"), TRUE)))
  expect_identical(mcp_close(client, pumps), 0L)
  expect_clean_streams(client)
})

test_that("an app that takes URL elicitations only is never sent a form", {
  testthat::skip_on_cran()
  local_agents_dir()
  viewer <- agents_viewer()
  page <- fake_page(viewer)
  on.exit(page_close(page), add = TRUE)
  page$card <- "save_copy"
  pumps <- list(function() page_pump(page))
  client <- mcp_client()
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  mcp_initialize(client, capabilities = list(elicitation = list(url = empty_object())), version = "2025-11-25", pumps = pumps)
  mcp_send(client, list(id = 2L, method = "tools/call", params = list(name = "save_copy", arguments = empty_object(), `_meta` = list(progressToken = 9L))))
  call <- page_next(page, "agent_call", pumps)
  mcp_next(client, function(message) identical(message$method, "notifications/progress") && identical(message$params$progress, 2L), pumps)
  page_card_answer(page, call$requestId, FALSE)
  expect_identical(mcp_reply(client, 2L, pumps)$result, text_result("declined in the card", TRUE))
  expect_false(any(vapply(client$held, function(message) identical(message$method, "elicitation/create"), TRUE)))
  expect_identical(mcp_close(client, pumps), 0L)
})

test_that("with no viewer named, the server leaves a viewer whose page went for one whose page is open", {
  testthat::skip_on_cran()
  local_agents_dir()
  older <- agents_viewer()
  newer <- agents_viewer()
  page <- fake_page(newer)
  client <- mcp_client()
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  mcp_initialize(client, pumps = list(function() page_pump(page)))
  mcp_send(client, list(id = 2L, method = "tools/call", params = list(name = "color_by", arguments = empty_object())))
  expect_identical(mcp_reply(client, 2L, list(function() page_pump(page)))$result, text_result("ran color_by as Claude Code"))
  page_close(page)
  pump_until(function() !cellucid:::.agent_port_connected(newer))
  other <- fake_page(older)
  on.exit(page_close(other), add = TRUE)
  pumps <- list(function() page_pump(other))
  mcp_send(client, list(id = 3L, method = "tools/call", params = list(name = "color_by", arguments = empty_object())))
  expect_identical(mcp_reply(client, 3L, pumps)$result, text_result("ran color_by as Claude Code"))
  expect_identical(page_next(other, "agent_call", pumps)$tool, "color_by")
  mcp_send(client, list(id = 4L, method = "tools/call", params = list(name = "cellucid_status", arguments = empty_object())))
  status <- status_of(mcp_reply(client, 4L, pumps)$result$content[[1L]]$text)
  expect_identical(status$viewer$viewer, older$viewer_id)
  expect_identical(mcp_close(client, pumps), 0L)
})

test_that("the server reports a viewer that stops while a call waits", {
  testthat::skip_on_cran()
  local_agents_dir()
  viewer <- agents_viewer()
  page <- fake_page(viewer)
  page$hold <- "slow"
  pumps <- list(function() page_pump(page))
  client <- mcp_client()
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  mcp_initialize(client, pumps = pumps)
  mcp_send(client, list(id = 2L, method = "tools/call", params = list(name = "slow", arguments = empty_object())))
  page_next(page, "agent_call", pumps)
  # The server is waiting for the call's events when the viewer stops.
  consumers <- cellucid:::.agent_port(viewer)$published$consumers
  pump_until(function() any(vapply(ls(consumers), function(name) length(ls(get(name, envir = consumers)$waiters)) > 0L, TRUE)), pumps)
  page_close(page)
  cellucid_stop(viewer)
  expect_identical(mcp_reply(client, 2L)$result,
                   text_result(cellucid:::.mcp_refusal_text("The Cellucid viewer stopped (the viewer was stopped); nothing more was done."), TRUE))
  expect_identical(mcp_close(client), 0L)
  expect_clean_streams(client)
})

test_that("open serves a file in the server's own viewer, and stops it when the app ends the server", {
  testthat::skip_on_cran()
  folder <- local_agents_dir()
  store <- file.path(tempfile("cellucid-mcp-open-"), "agents.zarr")
  suppressMessages(cellucid_prepare(store, "zarr", embeddings = list(X_umap = matrix(c(0, 0, 1, 1), 2, 2, byrow = TRUE)),
                                    dataset_id = "opened", dataset_name = "Opened"))
  arguments <- sprintf("open = %s, port = 0L, browser = FALSE, web_source_url = NULL, web_cache_dir = %s",
                       deparse(store), deparse(viewer_fixture_web_cache()))
  client <- mcp_client(args = c("--vanilla", "-e", mcp_expression(arguments)))
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  mcp_initialize(client)
  # The viewer's start, its package code and the web app's copy included:
  # until it says it started or could not, or the server ends.
  pump_until(function() {
    mcp_pump(client)
    any(startsWith(client$stderr, "Cellucid viewer: ") | grepl(" could not open: ", client$stderr, fixed = TRUE)) ||
      !client$process$is_alive()
  }, timeout = Inf, what = "the viewer to open")
  mcp_pump(client)
  expect_true(any(startsWith(client$stderr, "Cellucid viewer: ")), info = paste(client$stderr, collapse = "\n"))
  address <- sub("^Cellucid viewer: ", "", grep("^Cellucid viewer: ", client$stderr, value = TRUE)[[1L]])
  query <- function(key) sub(paste0(".*[?&]", key, "=([^&]+).*"), "\\1", address)
  served <- list(port = as.integer(sub("^http://127\\.0\\.0\\.1:([0-9]+)/.*", "\\1", address)),
                 viewer_id = query("viewerId"), viewer_token = query("viewerToken"))
  served$url <- paste0("http://127.0.0.1:", served$port)
  # Published too: its record names this child's viewer.
  record <- cellucid:::.agent_read_record(file.path(folder, paste0(served$viewer_id, ".json")))
  expect_identical(record$dataset, list(id = "opened", name = "Opened"))
  expect_identical(record$pid, client$process$get_pid())
  page <- fake_page(served)
  pumps <- list(function() page_pump(page))
  mcp_send(client, list(id = 2L, method = "tools/call", params = list(name = "color_by", arguments = empty_object())))
  expect_identical(mcp_reply(client, 2L, pumps)$result, text_result("ran color_by as Claude Code"))
  mcp_send(client, list(id = 3L, method = "tools/call", params = list(name = "cellucid_status", arguments = empty_object())))
  status <- status_of(mcp_reply(client, 3L, pumps)$result$content[[1L]]$text)
  expect_true(status$connected)
  expect_identical(status$viewer$dataset, list(id = "opened", name = "Opened"))
  expect_identical(mcp_close(client, pumps), 0L)
  # Its viewer stopped with it: the page was told to freeze, and the record is gone.
  pump_until(function() {
    page_pump(page)
    any(vapply(page$commands, function(command) identical(command$type, "freeze"), TRUE))
  }, timeout = 10, what = "the freeze")
  page_close(page)
  expect_false(file.exists(file.path(folder, paste0(served$viewer_id, ".json"))))
  expect_clean_streams(client)
})

test_that("cellucid_mcp() refuses open with viewer, and a bad argument, by name", {
  expect_error(cellucid_mcp(open = "data.h5ad", viewer = "pbmc"), "not both")
  expect_error(cellucid_mcp(port = -1), "port")
  expect_error(cellucid_mcp(browser = NA), "browser")
})

# =============================================================================
# cellucid-python's server driving an R viewer
# =============================================================================

test_that("cellucid-python's cellucid mcp drives an R viewer through its endpoint", {
  python <- require_python(
    "import mcp.server.lowlevel.server, cellucid.mcp, cellucid._agent_link",
    c("CELLUCID_MCP_PYTHON", "CELLUCID_PYTHON"), "cellucid mcp (cellucid[mcp], cellucid-python of the matching commit)"
  )
  testthat::skip_if_not_installed("processx")
  local_agents_dir()
  viewer <- agents_viewer()
  page <- fake_page(viewer)
  pumps <- list(function() page_pump(page))
  client <- mcp_client(command = python, args = c("-m", "cellucid.cli", "mcp"))
  on.exit(if (client$process$is_alive()) client$process$kill(), add = TRUE)
  mcp_initialize(client, capabilities = list(elicitation = empty_object()), pumps = pumps)
  mcp_send(client, list(id = 2L, method = "tools/list"))
  expect_identical(vapply(mcp_reply(client, 2L, pumps)$result$tools, `[[`, "", "name"), c("color_by", "cellucid_status"))
  # The same result (the SDK writes its keys in its own order).
  same_result <- function(actual, expected) {
    expect_identical(lapply(actual$content, function(item) item[c("type", "text")]), expected$content)
    expect_identical(actual$isError, expected$isError)
  }
  mcp_send(client, list(id = 3L, method = "tools/call", params = list(name = "color_by", arguments = empty_object())))
  same_result(mcp_reply(client, 3L, pumps)$result, text_result("ran color_by as Claude Code"))
  page$ask <- "save_copy"
  mcp_send(client, list(id = 4L, method = "tools/call", params = list(name = "save_copy", arguments = empty_object())))
  asked <- mcp_next(client, function(message) identical(message$method, "elicitation/create"), pumps)
  mcp_send(client, list(id = asked$id, result = list(action = "accept", content = empty_object())))
  expect_true(page_next(page, "agent_answer", pumps)$allowed)
  same_result(mcp_reply(client, 4L, pumps)$result, text_result("saved"))
  page$hold <- "slow"
  page$commands <- list()
  mcp_send(client, list(id = 5L, method = "tools/call", params = list(name = "slow", arguments = empty_object())))
  expect_identical(page_next(page, "agent_call", pumps)$tool, "slow")
  page_close(page)
  cellucid_stop(viewer)
  answer <- mcp_reply(client, 5L)$result
  expect_identical(answer$content[[1L]]$text, "The Cellucid viewer stopped (the viewer was stopped); nothing more was done.")
  expect_identical(mcp_close(client), 0L)
})

test_that("names written by others are marked as the tab marks them, byte for byte (one rule in three languages)", {
  # The examples mcp_tools.json carries are made by the web app's own functions, with a name meant to close the tags.
  rule <- cellucid:::.mcp_tools_document()$resultData
  expect_identical(rule$tags[c("open", "close")], list(open = "<cellucid-data>", close = "</cellucid-data>"))
  expect_identical(cellucid:::.mcp_result_data(rule$examples$data$json), rule$examples$data$text)
  expect_identical(cellucid:::.mcp_result_data(cellucid:::.js_json(jsonlite::parse_json(rule$examples$data$json))), rule$examples$data$text)
  expect_identical(cellucid:::.mcp_quote_name(rule$examples$name$name), rule$examples$name$text)
  expect_identical(sub("{choice}", cellucid:::.mcp_quote_name(rule$examples$name$name), cellucid:::.mcp_not_found, fixed = TRUE),
                   rule$examples$sentence$built)
  expect_identical(cellucid:::.mcp_refusal_text(rule$examples$sentence$built), rule$examples$sentence$text)
  # Nothing in a name can close a mark: one closing tag per marked name, the real one.
  expect_identical(lengths(regmatches(rule$examples$sentence$text, gregexpr(rule$tags$close, rule$examples$sentence$text, fixed = TRUE))), 1L)
  expect_identical(cellucid:::.mcp_tool_result(rule$examples$sentence$built, TRUE),
                   cellucid:::.js_json(list(content = list(list(type = "text", text = rule$examples$sentence$text)), isError = TRUE)))
  # A quoted name nobody marked gains no mark: its tags are escaped.
  expect_identical(cellucid:::.mcp_refusal_text(rule$examples$unmarked$built), rule$examples$unmarked$text)
  expect_false(grepl(rule$tags$close, rule$examples$unmarked$text, fixed = TRUE))
})
