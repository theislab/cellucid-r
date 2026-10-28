# The viewer in a Jupyter notebook: cellucid_view() in an IRkernel session
# shows the viewer in the cell's output as a notebook widget, the same one
# cellucid-python's show_anndata() shows (inst/widget/widget.js, byte for byte
# the Python package's cellucid/jupyter/static/widget.js). The widget frames
# the Cellucid app (`<web_source_url>/?viewer=embedded`) and relays between it
# and this R session over the kernel's own connection, a Jupyter comm:
#
# - the page's requests (its documents and the store's bytes, by Range) are
#   answered by .route(), the routes of the viewer's HTTP server, without one;
# - the viewer channel's frames go through .on_ws_open() over a socket that
#   is the comm (.widget_socket()), so the hello, the commands, the events
#   and the hooks are those of a viewer in a browser tab.
#
# No port is opened, so nothing depends on how the notebook's browser reaches
# this machine. IRkernel runs a comm's messages between cells, on R's one
# thread: while a cell runs, the viewer waits, and a cell cannot wait for the
# viewer's answer (.wait_for() says so at once). An httpuv server is no help
# here: IRkernel never runs R's event loop between cells, so it would answer
# nothing at all. IRkernel's comms carry no binary buffers, so bodies travel
# as base64 text.

# The widget protocol the notebook's widget manager speaks (ipywidgets 8).
.widget_protocol_version <- "2.1.0"

# The anywidget front-end module the widget is a model of: any version the
# notebook has, as the module the Python package's widget names.
.widget_module_version <- "*"

# A binary frame from the page arrives in pieces of at most this size (the
# front end's PIECE_BYTES); answers are sent in messages of at most this many
# body bytes, a message's shape and not a limit on what is read.
.widget_answer_message_bytes <- 8388608

.in_jupyter_kernel <- function() isTRUE(getOption("jupyter.in_kernel"))

# Google Colab's runtime (its Jupyter server sets COLAB_RELEASE_TAG, and
# every kernel it starts, R's included, inherits it).
.in_colab <- function() nzchar(Sys.getenv("COLAB_RELEASE_TAG"))

# Colab shows a widget only when its display names the widget manager to
# load, in the display's metadata: what google.colab.output's
# enable_custom_widget_manager() makes the Python kernel add to every widget
# display (google/colab/output/_widgets.py, version '5.0.0a'), and what
# anywidget adds itself. The manager loads the widget's module from jsDelivr
# and talks to the kernel over its comms, whatever the kernel's language.
.colab_widget_manager_url <- paste0(
  "https://ssl.gstatic.com/colaboratory-static/widgets/colab-cdn-widget-manager/",
  "2b70e893a8ba7c0f/manager.min.js"
)

# The metadata of the widget's display: Colab's widget manager there, nothing
# elsewhere.
.widget_display_metadata <- function() {
  if (!.in_colab()) return(NULL)
  list(`application/vnd.jupyter.widget-view+json` = list(
    colab = list(custom_widget_manager = list(url = .colab_widget_manager_url))
  ))
}

# The kernel's comms: IRkernel's comm manager, or a test's stand-in with the
# same three functions.
.kernel_comms <- function() {
  for (package in c("IRkernel", "IRdisplay")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      .abort("Showing the viewer in a Jupyter notebook needs the ", package, " package: install.packages(\"", package, "\").")
    }
  }
  manager <- IRkernel::comm_manager()
  # The widget protocol's version travels in the comm_open message's own
  # metadata, which the notebook's widget manager checks; IRkernel's
  # send_open() puts metadata in the content instead, so the message is made
  # here with the kernel's own reply and wire functions.
  kernel <- environment(manager$send_response)$.self
  send_multipart <- get("zmq.send.multipart", envir = asNamespace("IRkernel"))
  list(
    open = function(comm_id, data, on_msg) {
      comm <- manager$new_comm("jupyter.widget", comm_id)
      manager$register_comm(comm)
      comm$on_msg(on_msg)
      # Opened by the cell that runs: its request is the message's parent.
      parent <- kernel$executor$current_request %||% manager$parent_request
      message <- kernel$new_reply("comm_open", parent)
      message$content <- list(comm_id = comm_id, target_name = "jupyter.widget", data = data)
      message$metadata <- list(version = .widget_protocol_version)
      send_multipart(kernel$sockets$iopub, kernel$msg_to_wire(message), serialize = FALSE)
      invisible(comm)
    },
    send = function(comm_id, data) manager$send_msg(comm_id, "jupyter.widget", data),
    display = function(bundle, metadata) IRdisplay::publish_mimebundle(bundle, metadata)
  )
}

.widget_file <- function(name) {
  path <- system.file("widget", name, package = "cellucid")
  if (!nzchar(path)) .abort("The installed cellucid package has no widget/", name, ": reinstall it.")
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

# The address the browser loads the app from: an exact http(s) origin and
# optional path, without a query, fragment, credentials or trailing slash.
.validate_widget_source_url <- function(value) {
  ok <- is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^https?://[^/?#@[:space:]]+(/[^?#[:space:]]*)?$", value) && !endsWith(value, "/")
  if (!ok) {
    .abort(
      "In a Jupyter notebook the browser loads the viewer from web_source_url, which must be an ",
      "exact http(s) address without a query, fragment, credentials or trailing slash, such as ",
      "\"https://www.cellucid.com\"."
    )
  }
  value
}

# Show `viewer` in the cell's output, and route the widget's messages to it.
.widget_show <- function(viewer, web_source_url, height = 600L, comms = .kernel_comms()) {
  web_source_url <- .validate_widget_source_url(web_source_url)
  layout_id <- .random_hex(16L)
  model_id <- .random_hex(16L)
  widget <- new.env(parent = emptyenv())
  widget$model_id <- model_id
  widget$sockets <- list()
  widget$send <- function(content) comms$send(model_id, list(method = "custom", content = content))
  viewer$widget <- widget
  state <- list(
    `_anywidget_id` = "cellucid.Widget",
    `_dom_classes` = list(),
    `_esm` = .widget_file("widget.js"),
    `_css` = .widget_file("widget.css"),
    `_model_module` = "anywidget",
    `_model_module_version` = .widget_module_version,
    `_model_name` = "AnyModel",
    `_view_module` = "anywidget",
    `_view_module_version` = .widget_module_version,
    `_view_name` = "AnyView",
    `_binary` = "base64",
    height = as.integer(height),
    web_source_url = web_source_url,
    layout = paste0("IPY_MODEL_", layout_id)
  )
  comms$open(layout_id, list(state = list(
    `_model_module` = "@jupyter-widgets/base", `_model_module_version` = "2.0.0",
    `_model_name` = "LayoutModel", `_view_module` = "@jupyter-widgets/base",
    `_view_module_version` = "2.0.0", `_view_name` = "LayoutView"
  ), buffer_paths = list()), function(data) invisible(NULL))
  comms$open(model_id, list(state = state, buffer_paths = list()), function(data) {
    tryCatch(.widget_receive(viewer, state, comms, data), error = function(error) {
      .warn("The Cellucid widget could not handle a message: ", conditionMessage(error))
    })
  })
  comms$display(list(
    `application/vnd.jupyter.widget-view+json` = list(model_id = model_id, version_major = 2L, version_minor = 1L),
    `text/plain` = paste0("Cellucid viewer of ", viewer$dataset_name, " (a notebook widget)")
  ), .widget_display_metadata())
  invisible(viewer)
}

# One message of the widget's front end (a comm_msg's data).
.widget_receive <- function(viewer, state, comms, data) {
  widget <- viewer$widget
  if (identical(data$method, "request_state")) {
    # A front end that reloaded asks for the model again.
    comms$send(widget$model_id, list(method = "update", state = state, buffer_paths = list()))
    return(invisible(NULL))
  }
  if (!identical(data$method, "custom")) return(invisible(NULL))
  content <- data$content
  view <- content$v
  if (!is.character(view) || length(view) != 1L) return(invisible(NULL))
  socket_key <- function() paste0(view, "/", format(content$s, scientific = FALSE))
  switch(content$t %||% "",
    attach = widget$send(list(
      t = "credentials", v = view, viewerId = viewer$viewer_id, viewerToken = viewer$viewer_token
    )),
    fetch = .widget_answer(viewer, view, content$requests),
    open = {
      key <- socket_key()
      previous <- widget$sockets[[key]]
      if (!is.null(previous)) .widget_peer_closed(previous)
      ws <- .widget_socket(widget, view, content$s)
      widget$sockets[[key]] <- ws
      .on_ws_open(viewer, ws)
    },
    text = {
      ws <- widget$sockets[[socket_key()]]
      if (!is.null(ws) && !ws$closed && !is.null(ws$message_handler)) ws$message_handler(FALSE, content$text)
    },
    binary = {
      key <- socket_key()
      ws <- widget$sockets[[key]]
      if (is.null(ws) || ws$closed) return(invisible(NULL))
      ws$pieces[[length(ws$pieces) + 1L]] <- jsonlite::base64_dec(content$b64 %||% "")
      if (isTRUE(content$last)) {
        frame <- unlist(ws$pieces, use.names = FALSE)
        ws$pieces <- list()
        if (!is.null(ws$message_handler)) ws$message_handler(TRUE, if (is.null(frame)) raw(0) else frame)
      }
    },
    close = {
      key <- socket_key()
      ws <- widget$sockets[[key]]
      widget$sockets[[key]] <- NULL
      if (!is.null(ws)) .widget_peer_closed(ws)
    },
    detach = {
      # The output was removed (cleared, re-rendered, the page left).
      for (key in names(widget$sockets)) {
        if (startsWith(key, paste0(view, "/"))) {
          ws <- widget$sockets[[key]]
          widget$sockets[[key]] <- NULL
          .widget_peer_closed(ws)
        }
      }
    }
  )
  invisible(NULL)
}

# A socket of the viewer channel whose frames are the widget's messages: the
# parts of httpuv's WebSocket that .on_ws_open() and .viewer_send() use.
.widget_socket <- function(widget, view, socket_id) {
  ws <- new.env(parent = emptyenv())
  ws$request <- list(PATH_INFO = paste0("/", .viewer_channel_path))
  ws$closed <- FALSE
  ws$pieces <- list()
  ws$message_handler <- NULL
  ws$close_handler <- NULL
  ws$onMessage <- function(handler) ws$message_handler <- handler
  ws$onClose <- function(handler) ws$close_handler <- handler
  ws$send <- function(frame) {
    if (ws$closed) .abort("The viewer's channel is closed.")
    widget$send(if (is.raw(frame)) {
      list(t = "binary", v = view, s = socket_id, b64 = jsonlite::base64_enc(frame))
    } else {
      list(t = "text", v = view, s = socket_id, text = frame)
    })
  }
  ws$close <- function(code = 1000L, reason = "") {
    if (ws$closed) return(invisible(NULL))
    ws$closed <- TRUE
    widget$send(list(t = "close", v = view, s = socket_id, code = as.integer(code), reason = reason))
    if (!is.null(ws$close_handler)) ws$close_handler()
    invisible(NULL)
  }
  ws
}

# The page closed the socket (or its output went away): nothing is sent back.
.widget_peer_closed <- function(ws) {
  if (ws$closed) return(invisible(NULL))
  ws$closed <- TRUE
  if (!is.null(ws$close_handler)) ws$close_handler()
  invisible(NULL)
}

# The page's requests, as IRkernel hands them over: a data frame (its JSON
# reader simplifies an array of objects) or a list of records.
.widget_requests <- function(requests) {
  if (is.data.frame(requests)) {
    return(lapply(seq_len(nrow(requests)), function(row) as.list(requests[row, , drop = FALSE])))
  }
  if (is.list(requests)) requests else list()
}

# Answer the page's requests from the routes of the viewer's server.
.widget_answer <- function(viewer, view, requests) {
  widget <- viewer$widget
  answers <- list()
  size <- 0
  flush <- function() {
    if (length(answers) > 0L) widget$send(list(t = "answers", v = view, answers = answers))
    answers <<- list()
    size <<- 0
  }
  for (request in .widget_requests(requests)) {
    response <- .widget_route(viewer, request)
    body <- response$body
    if (length(body) > 0L && size > 0 && size + length(body) > .widget_answer_message_bytes) flush()
    answers[[length(answers) + 1L]] <- list(
      id = request$id, status = as.integer(response$status),
      headers = .widget_headers(response$headers), b64 = jsonlite::base64_enc(body)
    )
    size <- size + length(body)
  }
  flush()
}

.widget_headers <- function(headers) {
  out <- list()
  for (name in names(headers)) out[[name]] <- as.character(headers[[name]])
  out
}

# One request answered as the HTTP server answers it, with the body read whole.
.widget_route <- function(viewer, request) {
  if (!isTRUE(viewer$running)) return(.text_response(410L, "The viewer was stopped"))
  method <- request$method
  path <- request$path
  range <- request$range
  if (!is.character(method) || length(method) != 1L || !is.character(path) || length(path) != 1L) {
    return(.text_response(400L, "Malformed request"))
  }
  parts <- strsplit(path, "?", fixed = TRUE)[[1L]]
  req <- list(
    REQUEST_METHOD = method,
    PATH_INFO = paste0("/", parts[[1L]] %||% ""),
    QUERY_STRING = if (length(parts) > 1L) parts[[2L]] else "",
    HTTP_RANGE = if (is.character(range) && length(range) == 1L && !is.na(range)) range
  )
  response <- tryCatch(.route(viewer, req), error = function(error) {
    .inform("Cellucid viewer: a request for ", req$PATH_INFO, " failed: ", conditionMessage(error))
    .text_response(500L, "Internal server error")
  })
  body <- response$body
  if (is.list(body)) {
    # A file the HTTP server would stream from disk (removed after, when owned).
    file <- body$file
    bytes <- readBin(file, "raw", n = file.info(file)$size)
    if (isTRUE(body$owned)) unlink(file)
    body <- bytes
  } else if (is.character(body)) {
    body <- charToRaw(enc2utf8(paste(body, collapse = "")))
  }
  if (identical(method, "HEAD")) body <- raw(0)
  response$body <- body
  response
}
