# Talking to a viewer's server from the test process itself. httpuv runs the
# server's handlers on this R thread whenever later's event loop runs, so a
# client here must never block: HTTP requests go through curl's multi
# interface and the WebSocket through a non-blocking socket, and both turn
# the event loop between polls. Nothing is downloaded: the web app served is
# the small fixture build in fixtures/web (its own inventory and digests).

viewer_fixture_web_cache <- function() {
  cache <- tempfile("cellucid_web_")
  dir.create(cache)
  source <- testthat::test_path("fixtures", "web")
  for (file in list.files(source, recursive = TRUE, all.files = TRUE)) {
    dir.create(dirname(file.path(cache, file)), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(source, file), file.path(cache, file))
  }
  writeLines(
    "{\"version\":1,\"source_url\":\"https://www.cellucid.com\"}",
    file.path(cache, ".cellucid-web-source.json"),
    sep = ""
  )
  cache
}

# A viewer on a free port with the fixture build, stopped when the calling
# test ends.
test_viewer <- function(data, ..., env = parent.frame()) {
  viewer <- cellucid_serve(
    data, ...,
    port = 0L, web_cache_dir = viewer_fixture_web_cache(), web_source_url = NULL, quiet = TRUE
  )
  # Stopped when the calling test's frame exits (on.exit() registered in
  # that frame, as withr::defer() does).
  do.call(base::on.exit, list(bquote(cellucid_stop(.(viewer))), TRUE, FALSE), envir = env)
  # Its expression matrices laid out (in the background), so every key answers.
  if (!is.null(viewer$store)) cellucid:::.store_laid_out(viewer$store, timeout = 300)
  viewer
}

run_loop_until <- function(done, timeout = 10) {
  deadline <- Sys.time() + timeout
  while (!done()) {
    if (Sys.time() > deadline) stop("timed out")
    later::run_now(0.01)
  }
  invisible(TRUE)
}

# One HTTP request: list(status, headers (lower-case names), body (raw)).
http_request <- function(url, headers = list(), method = "GET") {
  handle <- curl::new_handle(customrequest = method, nobody = identical(method, "HEAD"))
  if (length(headers) > 0L) do.call(curl::handle_setheaders, c(list(handle), headers))
  result <- NULL
  failure <- NULL
  pool <- curl::new_pool()
  curl::curl_fetch_multi(url, done = function(response) result <<- response,
                         fail = function(message) failure <<- message, pool = pool, handle = handle)
  run_loop_until(function() {
    curl::multi_run(timeout = 0, poll = TRUE, pool = pool)
    !is.null(result) || !is.null(failure)
  })
  if (!is.null(failure)) stop(failure)
  list(
    status = result$status_code,
    headers = curl::parse_headers_list(result$headers),
    body = result$content
  )
}

http_json <- function(url) {
  response <- http_request(url)
  jsonlite::fromJSON(rawToChar(response$body), simplifyVector = FALSE)
}

# --- A minimal WebSocket client (RFC 6455) -----------------------------------

# The position of the last byte of the first CR LF CR LF in `bytes`, or 0.
handshake_end <- function(bytes) {
  n <- length(bytes)
  if (n < 4L) return(0L)
  at <- which(
    bytes[1:(n - 3L)] == as.raw(13L) & bytes[2:(n - 2L)] == as.raw(10L) &
      bytes[3:(n - 1L)] == as.raw(13L) & bytes[4:n] == as.raw(10L)
  )
  if (length(at) == 0L) 0L else at[[1L]] + 3L
}

ws_open <- function(viewer, path = "/_cellucid/ws", origin = viewer$url, host = NULL) {
  connection <- socketConnection("127.0.0.1", viewer$port, blocking = FALSE, open = "r+b", timeout = 5)
  if (is.null(host)) host <- paste0("127.0.0.1:", viewer$port)
  request <- paste0(
    "GET ", path, " HTTP/1.1\r\nHost: ", host, "\r\nUpgrade: websocket\r\n",
    "Connection: Upgrade\r\nSec-WebSocket-Key: dGhlIHNhbXBsZSBub25jZQ==\r\n",
    "Sec-WebSocket-Version: 13\r\n",
    if (!is.null(origin)) paste0("Origin: ", origin, "\r\n"), "\r\n"
  )
  writeBin(charToRaw(request), connection)
  client <- new.env()
  client$connection <- connection
  client$buffer <- raw(0)
  # The handshake answer, up to its blank line, found in the bytes: frames the
  # server sends at once (the close of a refused channel) can arrive in the
  # same read, and they are not text.
  run_loop_until(function() {
    client$buffer <- c(client$buffer, readBin(connection, "raw", 65536L))
    handshake_end(client$buffer) > 0L
  })
  end <- handshake_end(client$buffer)
  head <- rawToChar(client$buffer[seq_len(end)])
  client$status <- as.integer(sub("^HTTP/1\\.1 ([0-9]+).*", "\\1", strsplit(head, "\r\n")[[1L]][[1L]]))
  client$head <- head
  client$buffer <- if (length(client$buffer) > end) client$buffer[(end + 1L):length(client$buffer)] else raw(0)
  client
}

ws_frame <- function(opcode, payload) {
  mask <- as.raw(c(1, 2, 3, 4))
  n <- length(payload)
  header <- as.raw(0x80 + opcode)
  length_bytes <- if (n < 126) {
    as.raw(0x80 + n)
  } else if (n < 65536) {
    c(as.raw(0x80 + 126), as.raw(c(n %/% 256, n %% 256)))
  } else {
    c(as.raw(0x80 + 127), as.raw(rev((n %/% 256^(0:7)) %% 256)))
  }
  masked <- xor(payload, rep_len(mask, n))
  c(header, length_bytes, mask, masked)
}

ws_send_text <- function(client, text) {
  writeBin(ws_frame(1, charToRaw(enc2utf8(text))), client$connection)
}

ws_send_binary <- function(client, bytes) {
  writeBin(ws_frame(2, bytes), client$connection)
}

ws_send_json <- function(client, value) {
  ws_send_text(client, as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null")))
}

# The next frame from the server: list(opcode, payload), opcode 8 a close
# whose `code` is set.
ws_receive <- function(client, timeout = 10) {
  frame <- NULL
  # A frame that never comes says what did: the bytes held, if any.
  received <- function(error) {
    held <- client$buffer
    stop("No whole frame from the server (", conditionMessage(error), " after ", timeout, " s); ", length(held), " bytes held",
         if (length(held) > 0L) paste0(": ", paste(utils::head(as.character(held), 48L), collapse = " ")), ".",
         call. = FALSE)
  }
  tryCatch(run_loop_until(function() {
    client$buffer <- c(client$buffer, readBin(client$connection, "raw", 65536L))
    buffer <- client$buffer
    if (length(buffer) < 2L) return(FALSE)
    n <- as.integer(buffer[[2L]]) %% 128L
    at <- 3L
    if (n == 126L) {
      if (length(buffer) < 4L) return(FALSE)
      n <- as.integer(buffer[[3L]]) * 256L + as.integer(buffer[[4L]])
      at <- 5L
    } else if (n == 127L) {
      if (length(buffer) < 10L) return(FALSE)
      n <- sum(as.numeric(buffer[3:10]) * 256^(7:0))
      at <- 11L
    }
    if (length(buffer) < at + n - 1L) return(FALSE)
    payload <- if (n > 0) buffer[at:(at + n - 1L)] else raw(0)
    frame <<- list(opcode = as.integer(buffer[[1L]]) %% 16L, payload = payload)
    client$buffer <- if (length(buffer) > at + n - 1L) buffer[(at + n):length(buffer)] else raw(0)
    TRUE
  }, timeout = timeout), error = received)
  if (frame$opcode == 8L && length(frame$payload) >= 2L) {
    frame$code <- as.integer(frame$payload[[1L]]) * 256L + as.integer(frame$payload[[2L]])
  }
  if (frame$opcode == 1L) frame$json <- jsonlite::fromJSON(rawToChar(frame$payload), simplifyVector = FALSE)
  if (frame$opcode == 2L) frame$json <- binary_command(frame$payload)
  frame
}

# The close frame the server sends a socket it refused (421, 403): httpuv
# upgrades the connection after writing the refusal all the same -- its 101
# head, and the refusal's body before or after it -- and the server closes
# the socket as it opens (.on_server_ws_open()). Waiting for that close also
# lets the open run before the client goes: a socket closed first has its
# open skipped by httpuv but its close delivered, which warns "object '...'
# not found" from httpuv's own handler in whichever test turns the event
# loop next.
ws_refused_close <- function(client) {
  body <- as.integer(sub("(?s).*\r\ncontent-length: *([0-9]+).*", "\\1", tolower(client$head), perl = TRUE))
  frames <- NULL
  run_loop_until(function() {
    client$buffer <- c(client$buffer, readBin(client$connection, "raw", 65536L))
    at <- grepRaw("HTTP/1.1 101 ", client$buffer, fixed = TRUE)
    if (length(at) == 0L) return(FALSE)
    end <- handshake_end(client$buffer[at:length(client$buffer)])
    if (end == 0L) return(FALSE)
    after <- client$buffer[-seq_len(at + end - 1L)]
    rest <- body - (at - 1L)
    if (length(after) < rest) return(FALSE)
    frames <<- if (rest > 0L) after[-seq_len(rest)] else after
    TRUE
  })
  client$buffer <- frames
  ws_receive(client)
}

# `params` with `value` placed at `path` (keys, and 0-based list indices),
# its last key appended: where a binary action's array goes
# (data/viewer-protocol.js placeArray).
place_at <- function(params, path, value) {
  key <- path[[1L]]
  if (length(path) == 1L) {
    params[key] <- list(value)
    return(params)
  }
  at <- if (is.numeric(key)) key + 1L else key
  params[[at]] <- place_at(params[[at]], path[-1L], value)
  params
}

# A binary command read as the viewer reads it (decodeViewerBinaryCommand):
# its header, padded to a multiple of 4 bytes, and each array of an action's
# body placed in its params at its path, as numbers; a load_session's file.
binary_command <- function(payload) {
  header_length <- readBin(payload[1:4], "integer", size = 4L, endian = "little")
  testthat::expect_identical(header_length %% 4L, 0L)
  header <- jsonlite::fromJSON(rawToChar(payload[4L + seq_len(header_length)]), simplifyVector = FALSE)
  at <- 4L + header_length
  # A session file pushed in: its bytes follow the header.
  if (identical(header$type, "load_session")) return(c(header, list(file = payload[(at + 1L):length(payload)])))
  for (array in header$arrays) {
    size <- cellucid:::.viewer_value_types[[array$type]]
    bytes <- if (array$length == 0L) raw(0) else payload[at + seq_len(array$length * size)]
    values <- cellucid:::.decode_wire_values(bytes, array$type, array$length)
    header$params <- place_at(header$params, array$path, as.numeric(values))
    at <- at + array$length * size
  }
  testthat::expect_identical(at, length(payload))
  header
}

# The command a frame R sent holds, text or binary, as the viewer reads it.
sent_command <- function(frame) {
  if (is.raw(frame)) binary_command(frame) else jsonlite::fromJSON(frame, simplifyVector = FALSE)
}

# Let go of the socket, saying nothing: to the server the connection is lost
# (a browser that crashed) only once the last handle on it is closed. On
# Windows that may be long after this close(): R's client socket is made
# without WSA_FLAG_NO_HANDLE_INHERIT (src/modules/internet/Rsock.c; only its
# listening sockets have it) and a child started while the socket is open
# inherits a handle on it (processx and system2() both start children with
# bInheritHandles), so no FIN leaves until that child ends too. So a test that
# waits for the server to notice a socket going uses ws_send_close() when a
# child runs (page_close()), and this alone only when none was started while
# the socket was open.
ws_close <- function(client) {
  try(close(client$connection), silent = TRUE)
}

# What a browser sends when a page goes (1001, going away): the server answers
# and closes the connection itself, whoever else holds a handle on this end.
ws_send_close <- function(client, code = 1001L) {
  frame <- ws_frame(8, as.raw(c(code %/% 256L, code %% 256L)))
  # The server may have closed first (a stopped viewer): nothing to tell then.
  suppressWarnings(try(writeBin(frame, client$connection), silent = TRUE))
  invisible(NULL)
}

# A welcomed channel for `viewer`.
ws_welcomed <- function(viewer) {
  client <- ws_open(viewer)
  ws_send_json(client, list(type = "hello", viewerId = viewer$viewer_id, viewerToken = viewer$viewer_token))
  welcome <- ws_receive(client)
  testthat::expect_identical(welcome$json, list(type = "welcome"))
  client
}
