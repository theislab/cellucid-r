# The viewer's server: one httpuv server per viewer, mirroring the routes of
# cellucid-python's AnnData server (anndata_server.py, _server_base.py), so
# the web app's served-AnnData source reads it unchanged:
#
#   GET  /, /index.html, /<asset>           the verified web build (web-cache.R)
#   GET  /_cellucid/health                  {status, type, version, format, n_cells, n_genes}
#   GET  /_cellucid/info                    the same, with host and port
#   GET  /_cellucid/datasets                the dataset, its store and launch defaults
#   GET  /_cellucid/protocol                the viewer channel (viewer-protocol.R)
#   GET  /_cellucid/anndata/<name>[/<key>]  the store's bytes, with Range,
#                                           If-Range, If-None-Match and ETag
#   WS   /_cellucid/ws                      the viewer channel
#
# An R object is served as a Cellucid-ready AnnData (or MuData) Zarr store
# held in memory: the plan cellucid_prepare() writes (dataset-plan.R),
# uncompressed (the browser reads it over loopback), every chunk made when the
# browser asks for it, and the chunks made by a transform kept in a
# least-recently-used cache bounded by `cache_bytes` (zarr-store.R). A file or
# a Zarr folder given by path is served from disk. The object is read when the
# viewer starts: changing it afterwards does not change what the viewer shows
# (R copies on change), so view it again to see the change.
#
# httpuv does its network work on a background thread and runs these
# handlers on R's main thread whenever R is idle (an RStudio or Positron
# console waiting at its prompt), so the console stays usable while a viewer
# is open. Every request is checked first (onHeaders): a Host header that
# does not name this server is answered 421 (DNS rebinding), and a WebSocket
# upgrade from another origin 403.

.viewer_registry <- new.env(parent = emptyenv())

.json_response <- function(value, status = 200L) {
  list(
    status = status,
    headers = list(`Content-Type` = "application/json"),
    body = .json_text(value)
  )
}

.text_response <- function(status, message) {
  list(status = status, headers = list(`Content-Type` = "text/plain"), body = message)
}

# A body httpuv streams from a file, `owned` when it removes the file after.
# httpuv opens the file by the bytes of its name read as UTF-8 (on Windows
# with CreateFileW() after converting them from UTF-8), so the name is
# handed over in UTF-8: in an R session that is not UTF-8 (Windows R before
# 4.2) a native name with a character such as an e-acute -- a store under a
# user's folder, the web app's cache in it, the session's temporary folder
# -- would name no file.
.file_body <- function(path, owned = FALSE) {
  list(file = enc2utf8(path), owned = owned)
}

# --- Requests ----------------------------------------------------------------

.request_header <- function(req, name) {
  value <- req[[paste0("HTTP_", toupper(gsub("-", "_", name, fixed = TRUE)))]]
  if (is.null(value) || !nzchar(value)) NULL else value
}

# The host and port a Host (or Origin) authority names.
.parse_authority <- function(text) {
  if (is.null(text)) return(NULL)
  match <- regmatches(text, regexec("^(\\[[0-9A-Fa-f:.]+\\]|[^:/\\[\\]]+)(?::([0-9]{1,5}))?$", text, perl = TRUE))[[1L]]
  if (length(match) == 0L) return(NULL)
  list(host = tolower(gsub("^\\[|\\]$", "", match[[2L]])), port = if (nzchar(match[[3L]])) as.integer(match[[3L]]) else NULL)
}

# Whether a bind address is an IP literal (IPv4 dotted quad, or IPv6).
.is_ip_literal <- function(host) {
  host <- gsub("^\\[|\\]$", "", host)
  ipv4 <- regmatches(host, regexec("^([0-9]{1,3})\\.([0-9]{1,3})\\.([0-9]{1,3})\\.([0-9]{1,3})$", host))[[1L]]
  if (length(ipv4) > 0L) return(all(as.integer(ipv4[-1L]) <= 255L))
  grepl(":", host, fixed = TRUE) && grepl("^[0-9A-Fa-f:.]+$", host)
}

# The Host check is skipped only when the server was positively bound to a
# non-loopback IP address (0.0.0.0, ::, or one interface's address) and no
# allowed_hosts were named: an explicit request for network exposure, whose
# legitimate names cannot be known here. A host name bound, or anything
# undetermined, is checked (cellucid-python's _requires_loopback_host_header).
.host_check_applies <- function(viewer) {
  length(viewer$allowed_hosts) > 0L || !.is_ip_literal(viewer$host) || .is_loopback_host(viewer$host)
}

.host_is_authorized <- function(viewer, req) {
  allowed <- viewer$allowed_hosts
  if (!.host_check_applies(viewer)) return(TRUE)
  authority <- .parse_authority(.request_header(req, "Host"))
  if (is.null(authority)) return(FALSE)
  if (authority$host %in% allowed) return(TRUE)
  port <- if (is.null(authority$port)) 80L else authority$port
  identical(port, as.integer(viewer$port)) && .is_loopback_host(authority$host)
}

# A WebSocket upgrade is accepted from a page on this server's own host (the
# request's Host), on a loopback name, or on a host named in allowed_hosts;
# a handshake with no Origin is refused (cellucid-python's
# origin_is_allowed: a browser always sends one).
.origin_is_allowed <- function(viewer, req) {
  origin <- .request_header(req, "Origin")
  if (is.null(origin)) return(FALSE)
  parts <- regmatches(origin, regexec("^(https?)://([^/]+)$", origin))[[1L]]
  if (length(parts) == 0L) return(FALSE)
  authority <- .parse_authority(parts[[3L]])
  if (is.null(authority)) return(FALSE)
  if (.is_loopback_host(authority$host) || authority$host %in% viewer$allowed_hosts) return(TRUE)
  own <- .parse_authority(.request_header(req, "Host"))
  !is.null(own) && identical(own$host, authority$host)
}

# The CORS headers of every answer, as cellucid-python sends them: the
# requesting origin when it is local or Cellucid's own, the methods and
# headers a reader uses, and the headers a reader acts on exposed.
.cors_headers <- function(req) {
  headers <- list(
    `Access-Control-Allow-Methods` = "GET, HEAD, OPTIONS",
    `Access-Control-Allow-Headers` = "Content-Type, Range",
    `Access-Control-Expose-Headers` = "Accept-Ranges, Content-Encoding, Content-Length, Content-Range, ETag",
    `Cache-Control` = "no-cache"
  )
  origin <- .request_header(req, "Origin")
  if (!is.null(origin)) {
    parts <- regmatches(origin, regexec("^(https?)://([^/:]+)(:[0-9]+)?$", origin))[[1L]]
    allowed <- length(parts) > 0L && (
      parts[[3L]] %in% c("localhost", "127.0.0.1") ||
        identical(origin, .cellucid_web_url) ||
        (identical(parts[[2L]], "https") && (parts[[3L]] == "cellucid.com" || endsWith(parts[[3L]], ".cellucid.com")))
    )
    if (allowed) {
      headers <- c(list(`Access-Control-Allow-Origin` = origin, Vary = "Origin"), headers)
    }
  }
  headers
}

.with_cors <- function(response, req) {
  response$headers <- c(response$headers, .cors_headers(req))
  response
}

.on_headers <- function(viewer, req) {
  if (!.host_is_authorized(viewer, req)) {
    return(.with_cors(.text_response(421L, "Host header does not name this server"), req))
  }
  upgrade <- .request_header(req, "Upgrade")
  if (!is.null(upgrade) && identical(tolower(upgrade), "websocket") && !.origin_is_allowed(viewer, req)) {
    return(.with_cors(.text_response(403L, "WebSocket origin is not allowed"), req))
  }
  NULL
}

# One Range header against `size` bytes: NULL for none, c(start, end)
# (0-based, inclusive), or FALSE when it is not one satisfiable interval
# (cellucid-python's parse_byte_range).
.parse_byte_range <- function(value, size) {
  if (is.null(value)) return(NULL)
  match <- regmatches(value, regexec("^bytes=([0-9]{0,20})-([0-9]{0,20})$", value))[[1L]]
  if (length(match) == 0L || (!nzchar(match[[2L]]) && !nzchar(match[[3L]])) || size == 0) return(FALSE)
  if (nzchar(match[[2L]])) {
    start <- as.numeric(match[[2L]])
    end <- if (nzchar(match[[3L]])) as.numeric(match[[3L]]) else size - 1
    if (start >= size || end < start) return(FALSE)
    return(c(start, min(end, size - 1)))
  }
  suffix <- as.numeric(match[[3L]])
  if (suffix <= 0) return(FALSE)
  c(max(0, size - suffix), size - 1)
}

# The interval to answer: NULL for the whole resource, FALSE for 416. A
# Range sent with an If-Range naming another version is answered whole.
.range_interval <- function(req, size, etag) {
  range <- .request_header(req, "Range")
  condition <- .request_header(req, "If-Range")
  if (!is.null(range) && !is.null(condition) && !identical(trimws(condition), etag)) return(NULL)
  .parse_byte_range(range, size)
}

.not_modified <- function(req, etag) {
  condition <- .request_header(req, "If-None-Match")
  if (is.null(condition)) return(FALSE)
  tags <- sub("^W/", "", trimws(strsplit(condition, ",", fixed = TRUE)[[1L]]))
  etag %in% tags || "*" %in% tags
}

.range_not_satisfiable <- function(size) {
  list(
    status = 416L,
    headers = list(
      `Content-Type` = "text/plain", `Content-Range` = paste0("bytes */", format(size, scientific = FALSE)),
      `Accept-Ranges` = "bytes"
    ),
    body = "Requested byte range is not satisfiable"
  )
}

# The headers of the store's bytes. `Content-Encoding: identity` keeps httpuv
# from gzipping the answer for a browser that accepts gzip, which it does by
# default, even for a 206: gzip would run on R's main thread (the cost an
# uncompressed served store avoids; a file's chunks are compressed already),
# and a gzipped answer is sent chunked, without the Content-Length a reader
# sizes a key from (a HEAD then names no size at all) -- cellucid-python
# sends them as they are.
.bytes_headers <- function(size, interval, etag) {
  headers <- list(
    `Content-Type` = "application/octet-stream", `Content-Encoding` = "identity",
    `Accept-Ranges` = "bytes", ETag = etag
  )
  if (!is.null(interval)) {
    headers$`Content-Range` <- paste0(
      "bytes ", format(interval[[1L]], scientific = FALSE), "-",
      format(interval[[2L]], scientific = FALSE), "/", format(size, scientific = FALSE)
    )
  }
  headers
}

# Bytes held in memory, or a range of them: exactly the bytes asked for.
.bytes_response <- function(req, bytes, etag, head_only) {
  if (.not_modified(req, etag)) return(list(status = 304L, headers = list(ETag = etag), body = raw(0)))
  size <- length(bytes)
  interval <- .range_interval(req, size, etag)
  if (isFALSE(interval)) return(.range_not_satisfiable(size))
  body <- if (is.null(interval)) bytes else bytes[(interval[[1L]] + 1):(interval[[2L]] + 1)]
  # A HEAD answer is given its body too: httpuv sends a HEAD's headers
  # alone, with the body's Content-Length, which is what a reader sizes the
  # resource from.
  list(
    status = if (is.null(interval)) 200L else 206L,
    headers = .bytes_headers(size, interval, etag),
    body = body
  )
}

# The bytes read from disk at once for a range; a longer range is copied
# through a temporary file in blocks of this size, so its bytes are never
# held in R's memory whole (a scheduling choice, not a limit on ranges).
.range_block_bytes <- 8388608

# A file on disk, or a range of it. A whole file is streamed by httpuv from
# disk; a range is read with one seek and one read.
.file_response <- function(req, path, head_only) {
  info <- file.info(path, extra_cols = FALSE)
  if (is.na(info$size) || isTRUE(info$isdir)) return(.text_response(404L, paste0(basename(path), " is no longer readable")))
  size <- info$size
  etag <- paste0("\"", as.hexmode(as.integer(size %% 2^31)), "-", format(as.numeric(info$mtime) * 1e6, scientific = FALSE), "\"")
  if (.not_modified(req, etag)) return(list(status = 304L, headers = list(ETag = etag), body = raw(0)))
  interval <- .range_interval(req, size, etag)
  if (isFALSE(interval)) return(.range_not_satisfiable(size))
  headers <- .bytes_headers(size, interval, etag)
  status <- if (is.null(interval)) 200L else 206L
  # httpuv streams a file body from disk on its own thread, and answers a
  # HEAD with the body's Content-Length and no body.
  if (size == 0) {
    return(list(status = status, headers = headers, body = raw(0)))
  }
  if (is.null(interval) || (interval[[1L]] == 0 && interval[[2L]] == size - 1)) {
    return(list(status = status, headers = headers, body = .file_body(path)))
  }
  length <- interval[[2L]] - interval[[1L]] + 1
  connection <- file(path, "rb")
  on.exit(close(connection))
  seek(connection, interval[[1L]])
  if (length <= .range_block_bytes) {
    return(list(status = status, headers = headers, body = readBin(connection, "raw", n = length)))
  }
  # A long range is copied block by block into a file httpuv sends and then
  # removes (owned), so memory holds one block at a time.
  copy <- tempfile("cellucid-range-")
  output <- file(copy, "wb")
  left <- length
  tryCatch({
    while (left > 0) {
      block <- readBin(connection, "raw", n = min(left, .range_block_bytes))
      if (length(block) == 0L) .abort(basename(path), " ended before the range it was asked for")
      writeBin(block, output)
      left <- left - length(block)
    }
    close(output)
  }, error = function(error) {
    try(close(output), silent = TRUE)
    unlink(copy)
    stop(error)
  })
  list(status = status, headers = headers, body = .file_body(copy, owned = TRUE))
}

# --- The served store ------------------------------------------------------------

# A key of the store as the browser names it: a plain relative key, or NULL.
.canonical_store_key <- function(text) {
  if (!nzchar(text) || startsWith(text, "/") || grepl("[\\\\\001-\037]", text)) return(NULL)
  parts <- strsplit(text, "/", fixed = TRUE)[[1L]]
  if (endsWith(text, "/") || any(parts %in% c("", ".", ".."))) return(NULL)
  text
}

.store_response <- function(viewer, req, route, head_only) {
  served <- viewer$served
  not_found <- function() .text_response(404L, paste0("Not found: _cellucid/anndata/", route))
  if (!identical(served$layout, "folder")) {
    if (!identical(route, served$route_name)) return(not_found())
    return(.file_response(req, served$path, head_only))
  }
  prefix <- paste0(served$route_name, "/")
  if (!startsWith(route, prefix)) return(not_found())
  key <- .canonical_store_key(substring(route, nchar(prefix) + 1L))
  if (is.null(key)) return(not_found())
  if (!is.null(viewer$store)) {
    bytes <- .store_get(viewer$store, key)
    if (is.null(bytes)) return(.text_response(404L, paste0("No key '", key, "' in the store")))
    return(.bytes_response(req, bytes, viewer$etag, head_only))
  }
  generated <- served$generated
  if (!is.null(generated) && identical(key, generated$key)) {
    return(.bytes_response(req, generated$bytes, viewer$etag, head_only))
  }
  path <- file.path(served$path, key)
  if (!file.exists(path) || dir.exists(path)) return(.text_response(404L, paste0("No key '", key, "' in the store")))
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!startsWith(resolved, paste0(served$resolved, "/"))) {
    return(.text_response(404L, paste0("No key '", key, "' in the store")))
  }
  .file_response(req, resolved, head_only)
}

# --- HTTP routes -----------------------------------------------------------------

.on_call <- function(viewer, req) {
  response <- tryCatch(.route(viewer, req), error = function(error) {
    # The error is the R session's to read, not the page's: a request is
    # answered without the internal message.
    .inform("Cellucid viewer: a request for ", req$PATH_INFO %||% "/", " failed: ", conditionMessage(error))
    .text_response(500L, "Internal server error")
  })
  .with_cors(response, req)
}

.route <- function(viewer, req) {
  method <- req$REQUEST_METHOD
  if (identical(method, "OPTIONS")) return(list(status = 204L, headers = list(), body = raw(0)))
  if (!method %in% c("GET", "HEAD")) {
    return(.text_response(405L, paste0(method, " is not supported; the viewer channel is the WebSocket at /", .viewer_channel_path)))
  }
  head_only <- identical(method, "HEAD")
  # A path that is not valid percent-encoded UTF-8 (a %00, a lone %ff) names
  # nothing this server holds.
  path <- tryCatch({
    raw_path <- req$PATH_INFO %||% "/"
    # URLdecode() would cut the path at an encoded NUL.
    if (grepl("%00", raw_path, fixed = TRUE)) .abort("an encoded NUL")
    # URLdecode() returns the bytes as native text; they are UTF-8.
    decoded <- .as_utf8(utils::URLdecode(raw_path))
    if (is.na(decoded)) .abort("not UTF-8")
    decoded
  }, error = function(error) NULL)
  if (is.null(path)) return(.text_response(400L, "The request path is not valid percent-encoded UTF-8"))
  if (!is.null(viewer$web)) {
    asset <- if (path %in% c("/", "/index.html")) "index.html" else substring(path, 2L)
    if (startsWith(path, "/") && !identical(asset, .web_inventory_file) && !grepl("[\\\\?#[:space:]]", asset)) {
      found <- .read_web_asset(viewer$web, asset)
      if (!is.null(found)) {
        return(list(
          status = 200L, headers = list(`Content-Type` = found$content_type),
          # httpuv drops a HEAD's body and keeps its Content-Length.
          body = .file_body(found$file)
        ))
      }
    }
    if (path == paste0("/", .web_inventory_file) || path == "/assets" || startsWith(path, "/assets/")) {
      return(.text_response(404L, "Web asset is not declared by the active build"))
    }
  }
  route <- sub("^/+", "", path)
  if (route %in% c("", "index.html")) return(.text_response(503L, "Cellucid viewer UI unavailable"))
  switch(route,
    `_cellucid/health` = return(.json_response(viewer$health)),
    `_cellucid/info` = return(.json_response(viewer$info)),
    `_cellucid/datasets` = return(.json_response(viewer$datasets)),
    `_cellucid/protocol` = return(.json_response(.viewer_protocol_document()))
  )
  if (startsWith(route, "_cellucid/anndata/")) {
    return(.store_response(viewer, req, substring(route, nchar("_cellucid/anndata/") + 1L), head_only))
  }
  .text_response(404L, paste0("Not found: ", route))
}

# --- The WebSocket channel ---------------------------------------------------------

.same_secret <- function(given, expected) {
  a <- charToRaw(given)
  b <- charToRaw(expected)
  if (length(a) != length(b)) return(FALSE)
  sum(as.integer(xor(a, b))) == 0L
}

# How long a socket may stay open without saying hello, and how often an
# open channel is pinged (cellucid-python's KEEPALIVE_SECONDS): a proxy
# between the page and this server (RStudio Server, Workbench, a notebook
# proxy) closes a connection that carries nothing for a while.
.viewer_hello_seconds <- 60
.viewer_keepalive_seconds <- 20

# The requestId prefix of the keep-alive pings: their pongs keep the
# connection busy and are not events of the caller's.
.keepalive_prefix <- "keepalive-"

.on_ws_open <- function(viewer, ws) {
  path <- ws$request$PATH_INFO %||% ""
  if (!identical(sub("^/+", "", path), .viewer_channel_path)) {
    ws$close(.viewer_close_codes[["protocolViolation"]], "Not the viewer channel")
    return(invisible(NULL))
  }
  state <- new.env(parent = emptyenv())
  state$welcomed <- FALSE
  state$closed <- FALSE
  close_socket <- function(code, reason) {
    state$closed <- TRUE
    try(ws$close(code, enc2utf8(substr(reason, 1L, 120L))), silent = TRUE)
    invisible(NULL)
  }
  violation <- function(reason) close_socket(.viewer_close_codes[["protocolViolation"]], reason)
  # A socket that never says hello is closed, so it holds nothing open.
  later::later(function() {
    if (!state$welcomed && !state$closed) {
      close_socket(.viewer_close_codes[["protocolViolation"]], "No hello frame was received")
    }
  }, .viewer_hello_seconds)
  ws$onMessage(function(binary, message) {
    # Nothing a socket sends after it was closed for breaking the rules is
    # read: the frames httpuv already held are dropped here.
    if (state$closed) return(invisible(NULL))
    if (!state$welcomed) {
      hello <- tryCatch(
        if (binary) .refuse_frame("The first frame must be the text frame hello.") else .parse_viewer_hello(message),
        cellucid_protocol_error = function(error) error
      )
      if (inherits(hello, "error")) return(violation(conditionMessage(hello)))
      if (!identical(hello$viewerId, viewer$viewer_id) || !.same_secret(hello$viewerToken, viewer$viewer_token)) {
        close_socket(.viewer_close_codes[["unauthorized"]], "Invalid viewer credentials")
        return(invisible(NULL))
      }
      previous <- viewer$socket
      viewer$socket <- ws
      state$welcomed <- TRUE
      if (!is.null(previous)) {
        try(previous$close(.viewer_close_codes[["replaced"]], "Replaced by a newer connection"), silent = TRUE)
      }
      ws$send(.json_text(list(type = "welcome")))
      queued <- viewer$queue
      viewer$queue <- list()
      for (frame in queued) ws$send(frame)
      .keep_alive(viewer, ws, state)
      return(invisible(NULL))
    }
    event <- tryCatch(
      .check_event_cells(viewer, if (binary) .parse_viewer_binary_event(message) else .parse_viewer_event(message)),
      cellucid_protocol_error = function(error) error
    )
    if (inherits(event, "error")) return(violation(conditionMessage(event)))
    .deliver_event(viewer, event)
  })
  ws$onClose(function() {
    state$closed <- TRUE
    if (identical(viewer$socket, ws)) viewer$socket <- NULL
  })
  invisible(NULL)
}

# A ping command every .viewer_keepalive_seconds while this socket is the
# viewer's channel; the viewer answers each with a pong, which is not
# delivered (.deliver_event()).
.keep_alive <- function(viewer, ws, state) {
  tick <- function() {
    if (state$closed || !isTRUE(viewer$running) || !identical(viewer$socket, ws)) return(invisible(NULL))
    sent <- tryCatch({
      ws$send(.viewer_command_frame(list(type = "ping", requestId = paste0(.keepalive_prefix, .random_hex(4L)))))
      TRUE
    }, error = function(error) FALSE)
    if (sent) later::later(tick, .viewer_keepalive_seconds)
    invisible(NULL)
  }
  later::later(tick, .viewer_keepalive_seconds)
}

# Every cell index an event names is one of the dataset's cells (when the
# server knows how many it has): a page in a browser is not a trusted caller.
.check_event_cells <- function(viewer, event) {
  n <- viewer$n_cells
  if (is.null(n)) return(event)
  cells <- switch(event$type,
    selection = event$payload$cells,
    hover = event$payload$cell,
    click = event$payload$cell,
    NULL
  )
  if (length(cells) > 0L && max(cells) >= n) {
    .refuse_frame(event$type, " names cell ", format(max(cells), scientific = FALSE),
                  ", and the dataset has ", n, " ", .plural(n, "cell"), ".")
  }
  event
}

# Send one checked command, or keep it until the viewer's channel is open.
.viewer_send <- function(viewer, command) {
  frame <- .viewer_command_frame(command)
  if (!isTRUE(viewer$running)) {
    .abort("This viewer is stopped; open a new one with cellucid_view().")
  }
  if (is.null(viewer$socket)) {
    viewer$queue[[length(viewer$queue) + 1L]] <- frame
  } else {
    viewer$socket$send(frame)
  }
  invisible(viewer)
}

# An event from the viewer: kept as the latest of its type, and handed to
# every hook registered for it (cellucid_on()).
.deliver_event <- function(viewer, event) {
  type <- event$type
  payload <- event$payload
  # Cells are R's 1-based positions in every event, as cellucid_highlight()
  # and the other commands take them; the viewer's wire is 0-based.
  position <- function(index) {
    index <- index + 1
    if (all(index <= .Machine$integer.max)) as.integer(index) else index
  }
  if (identical(type, "pong") && startsWith(payload$requestId, .keepalive_prefix)) {
    return(invisible(NULL))
  }
  if (identical(type, "selection")) {
    indices <- position(payload$cells)
    payload <- list(
      cells = indices,
      cell_names = if (!is.null(viewer$cell_names)) viewer$cell_names[indices],
      source = payload$source
    )
    viewer$selection <- payload
  } else if (type %in% c("hover", "click") && !is.null(payload$cell)) {
    payload$cell <- position(payload$cell)
    if (!is.null(viewer$cell_names)) payload$cell_name <- viewer$cell_names[[payload$cell]]
  } else if (identical(type, "session_bundle")) {
    request <- viewer$session_requests[[payload$requestId]]
    if (is.null(request)) {
      # Nobody waits for it any more (the wait timed out): dropped, and the
      # channel stays open.
      return(invisible(NULL))
    }
    viewer$session_requests[[payload$requestId]] <- NULL
    if (!.is_session_bundle(event$body)) {
      # The one waiting is told now, in the form its wait accepts.
      return(.deliver_event(viewer, list(type = "command_error", payload = list(
        command = "requestSessionBundle",
        reason = paste0("The viewer sent an invalid session bundle: ", .not_a_session_file)
      ))))
    }
    path <- tempfile("cellucid-session-", fileext = ".cellucid-session")
    writeBin(event$body, path)
    payload <- list(requestId = payload$requestId, status = "ok", bytes = length(event$body), path = path)
  } else if (identical(type, "data_mapping")) {
    viewer$data_mapping <- payload
  } else if (identical(type, "ready")) {
    viewer$ready <- payload
  }
  viewer$last[[type]] <- payload
  for (hook in viewer$hooks[[type]]) {
    tryCatch(hook(payload), error = function(error) {
      .warn("A cellucid_on(\"", type, "\") hook failed: ", conditionMessage(error))
    })
  }
  invisible(NULL)
}

# A .cellucid-session file starts with its format's signature (the web app's
# bundle/format.js, cellucid-python's SESSION_BUNDLE_MAGIC).
.session_format_version <- 2L
.session_bundle_magic <- charToRaw(paste0("CELLUCID_SESSION/", .session_format_version, "\n"))
.not_a_session_file <- paste0(
  "This file is not a Cellucid session file of format ", .session_format_version,
  " (it does not start with CELLUCID_SESSION/", .session_format_version, ")."
)

.is_session_bundle <- function(body) {
  length(body) >= length(.session_bundle_magic) &&
    identical(body[seq_along(.session_bundle_magic)], .session_bundle_magic)
}
