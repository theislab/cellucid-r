# The viewer's server and channel under conditions the happy path does not
# reach: HEAD answers, long ranges, stale answers, cell indices out of range,
# sockets that never say hello or keep sending after a violation, the
# package unloading, and the random credentials.

served_url <- function(viewer, key) {
  paste0(viewer$url, "/_cellucid/anndata/", viewer$served$route_name, "/", key)
}

test_that("viewer credentials come from the system and leave .Random.seed alone", {
  set.seed(42)
  first <- cellucid:::.random_hex(16L)
  after_first <- runif(1L)
  set.seed(42)
  second <- cellucid:::.random_hex(16L)
  set.seed(42)
  expected <- runif(1L)
  expect_match(first, "^[0-9a-f]{32}$")
  # set.seed() does not make a token predictable ...
  expect_false(identical(first, second))
  # ... and drawing one does not move the caller's random stream.
  expect_identical(after_first, expected)
})

test_that("a HEAD answer carries the Content-Length of the GET", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  for (key in c(".zmetadata", "obsm/X_UMAP/0.0")) {
    get <- http_request(served_url(viewer, key))
    head <- http_request(served_url(viewer, key), method = "HEAD")
    expect_identical(head$status, 200L, info = key)
    expect_identical(as.numeric(head$headers[["content-length"]]), as.numeric(length(get$body)), info = key)
    expect_length(head$body, 0L)
    # httpuv gzips any answer to a client that accepts gzip, unless told the
    # body is sent as it is: then there would be no Content-Length at all.
    expect_identical(head$headers[["content-encoding"]], "identity", info = key)
    expect_identical(get$headers[["content-encoding"]], "identity", info = key)
  }
  ranged <- http_request(served_url(viewer, ".zmetadata"), list(Range = "bytes=2-9"), method = "HEAD")
  expect_identical(ranged$status, 206L)
  expect_identical(ranged$headers[["content-length"]], "8")
  plain <- list(`Accept-Encoding` = "identity")
  asset <- http_request(paste0(viewer$url, "/index.html"), plain)
  head <- http_request(paste0(viewer$url, "/index.html"), plain, method = "HEAD")
  expect_identical(as.numeric(head$headers[["content-length"]]), as.numeric(length(asset$body)))
})

test_that("a file's range is read block by block, never whole, and a HEAD sizes it", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(file.path(tempfile(), "f.zarr.zip"))
  args$format <- "zarr.zip"
  zip <- suppressMessages(do.call(cellucid_prepare, args))
  viewer <- test_viewer(zip)
  url <- paste0(viewer$url, "/_cellucid/anndata/", viewer$served$route_name)
  bytes <- readBin(zip, "raw", file.size(zip))
  head <- http_request(url, method = "HEAD")
  expect_identical(as.numeric(head$headers[["content-length"]]), as.numeric(length(bytes)))
  # A range longer than one block goes through a file httpuv sends and
  # removes; blocks of 16 bytes make every range here such a range.
  local_mocked_bindings(.range_block_bytes = 16, .package = "cellucid")
  before <- list.files(tempdir(), pattern = "^cellucid-range-")
  part <- http_request(url, list(Range = "bytes=100-1099"))
  expect_identical(part$status, 206L)
  expect_identical(part$body, bytes[101:1100])
  expect_identical(part$headers[["content-range"]], paste0("bytes 100-1099/", length(bytes)))
  whole <- http_request(url, list(Range = "bytes=0-"))
  expect_identical(whole$status, 206L)
  expect_identical(whole$body, bytes)
  run_loop_until(function() {
    length(setdiff(list.files(tempdir(), pattern = "^cellucid-range-"), before)) == 0L
  })
  expect_length(setdiff(list.files(tempdir(), pattern = "^cellucid-range-"), before), 0L)
})

test_that("a damaged request path is a 400, and a failed request names nothing internal", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  for (path in c("/_cellucid/%00", "/_cellucid/%ff")) {
    response <- http_request(paste0(viewer$url, path))
    expect_identical(response$status, 400L, info = path)
    expect_false(grepl("Error|error in|URLdecode", rawToChar(response$body)), info = path)
  }
  # A web asset changed on disk after it was served is refused (it no longer
  # matches its digest), and the page is told only that the request failed.
  first <- http_request(paste0(viewer$url, "/index.html"))
  expect_identical(first$status, 200L)
  index <- file.path(viewer$web$dir, "index.html")
  writeBin(c(readBin(index, "raw", file.size(index)), charToRaw("<!-- changed -->")), index)
  failed <- suppressMessages(http_request(paste0(viewer$url, "/index.html")))
  expect_identical(failed$status, 500L)
  expect_identical(rawToChar(failed$body), "Internal server error")
})

test_that("a verified web asset is served from disk, not re-read and re-hashed per request", {
  cache <- viewer_fixture_web_cache()
  build <- cellucid:::.establish_web_build(cache, NULL, quiet = TRUE)
  found <- cellucid:::.read_web_asset(build, "index.html")
  expect_identical(found$file, file.path(cache, "index.html"))
  expect_true(exists("index.html", envir = build$verified, inherits = FALSE))
  calls <- 0L
  local_mocked_bindings(
    .verify_asset_bytes = function(asset, bytes) calls <<- calls + 1L,
    .package = "cellucid"
  )
  for (i in 1:3) cellucid:::.read_web_asset(build, "index.html")
  expect_identical(calls, 0L)
  # Rewritten in place, the same size, its modification time put back: it
  # is checked again.
  index <- file.path(cache, "index.html")
  modified <- file.mtime(index)
  bytes <- readBin(index, "raw", file.size(index))
  Sys.sleep(0.05)
  writeBin(rev(bytes), index)
  Sys.setFileTime(index, modified)
  cellucid:::.read_web_asset(build, "index.html")
  expect_identical(calls, 1L)
})

test_that("the Host check applies unless the bind is a non-loopback IP address", {
  applies <- function(host, allowed = character()) {
    cellucid:::.host_check_applies(list(host = host, allowed_hosts = allowed))
  }
  expect_true(applies("127.0.0.1"))
  expect_true(applies("::1"))
  expect_true(applies("localhost"))
  # A host name bound is checked (cellucid-python's rule) ...
  expect_true(applies("my-laptop.local"))
  # ... an explicit non-loopback IP bind is a request for network exposure ...
  expect_false(applies("0.0.0.0"))
  expect_false(applies("::"))
  expect_false(applies("192.168.1.20"))
  # ... unless allowed_hosts names the proxy's host.
  expect_true(applies("0.0.0.0", "hub.example.org"))
})

test_that("hover and click name cells by 1-based position, and out-of-range cells are refused", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  viewer <- test_viewer(object, dataset_id = "s", dataset_name = "S")
  seen <- list()
  cellucid_on(viewer, "click", function(event) seen$click <<- event)
  client <- ws_welcomed(viewer)
  ws_send_json(client, list(type = "click", cell = 0, button = 0, shift = FALSE, ctrl = FALSE))
  run_loop_until(function() !is.null(seen$click))
  expect_identical(seen$click$cell, 1L)
  expect_identical(seen$click$cell_name, colnames(object)[[1L]])
  # The selection of the same cell names the same position.
  header <- charToRaw("{\"type\":\"selection\",\"source\":\"lasso\"}")
  ws_send_binary(client, c(cellucid:::.encode_uint32(length(header)), header, cellucid:::.encode_uint32(0)))
  run_loop_until(function() !is.null(viewer$selection))
  expect_identical(viewer$selection$cells, seen$click$cell)
  # A cell the dataset does not have closes the channel (1008).
  ws_send_binary(client, c(cellucid:::.encode_uint32(length(header)), header, cellucid:::.encode_uint32(c(1, 12))))
  closed <- ws_receive(client)
  expect_identical(closed$opcode, 8L)
  expect_identical(closed$code, 1008L)
  expect_identical(viewer$selection$cells, 1L)
  ws_close(client)
})

test_that("nothing a socket sends after a violation is delivered", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  clicks <- 0L
  cellucid_on(viewer, "click", function(event) clicks <<- clicks + 1L)
  client <- ws_welcomed(viewer)
  # Two frames written at once: a frame of no known type, then a click.
  bad <- cellucid:::.json_text(list(type = "zoom", level = 2))
  good <- cellucid:::.json_text(list(type = "click", cell = 1, button = 0, shift = FALSE, ctrl = FALSE))
  writeBin(c(ws_frame(1, charToRaw(bad)), ws_frame(1, charToRaw(good))), client$connection)
  closed <- ws_receive(client)
  expect_identical(closed$code, 1008L)
  for (i in 1:20) later::run_now(0.02)
  expect_identical(clicks, 0L)
  ws_close(client)
})

test_that("a socket that never says hello is closed, and an open channel is kept alive", {
  skip_without_sce()
  local_mocked_bindings(.viewer_hello_seconds = 0.3, .viewer_keepalive_seconds = 0.3, .package = "cellucid")
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  silent <- ws_open(viewer)
  closed <- ws_receive(silent, timeout = 5)
  expect_identical(closed$opcode, 8L)
  expect_identical(closed$code, 1008L)
  ws_close(silent)

  pongs <- 0L
  cellucid_on(viewer, "pong", function(event) pongs <<- pongs + 1L)
  client <- ws_welcomed(viewer)
  ping <- ws_receive(client, timeout = 5)$json
  expect_identical(ping$type, "ping")
  expect_match(ping$requestId, "^keepalive-")
  ws_send_json(client, list(type = "pong", requestId = ping$requestId, t = 1))
  for (i in 1:20) later::run_now(0.02)
  # The keep-alive's answer is not an event of the caller's.
  expect_identical(pongs, 0L)
  expect_null(viewer$last$pong)
  ws_close(client)
})

test_that("cellucid_debug_connection reports the routes, the store and the channel", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  unconnected <- cellucid_debug_connection(viewer, timeout = 5)
  expect_identical(unconnected$server_health$status, "ok")
  expect_identical(unconnected$server_info$port, viewer$port)
  expect_identical(unconnected$server_datasets$datasets[[1L]]$id, "s")
  expect_identical(unconnected$anndata_store_probe$head$status, 200L)
  expect_identical(unconnected$anndata_store_probe$range$status, 206L)
  expect_false(unconnected$channel_connected)
  expect_false(unconnected$frontend_roundtrip$ok)

  client <- ws_welcomed(viewer)
  # The page's side, answering each command as it arrives without blocking
  # R's event loop (the report waits on that loop for its HTTP probes too).
  answering <- TRUE
  respond <- function() {
    if (!answering) return(invisible(NULL))
    client$buffer <- c(client$buffer, readBin(client$connection, "raw", 65536L))
    buffer <- client$buffer
    n <- if (length(buffer) >= 2L) as.integer(buffer[[2L]]) %% 128L else NA
    if (!is.na(n) && n < 126L && length(buffer) >= 2L + n) {
      command <- jsonlite::fromJSON(rawToChar(buffer[2L + seq_len(n)]), simplifyVector = FALSE)
      client$buffer <- if (length(buffer) > 2L + n) buffer[(3L + n):length(buffer)] else raw(0)
      if (identical(command$type, "ping")) {
        ws_send_json(client, list(type = "pong", requestId = command$requestId, t = 5))
      } else if (identical(command$type, "debug_snapshot")) {
        ws_send_json(client, list(
          type = "debug_snapshot", requestId = command$requestId, ts = "now",
          locationHref = viewer$viewer_url, origin = viewer$url, serverUrl = viewer$url,
          connected = TRUE, userAgent = "test"
        ))
      }
    }
    later::later(respond, 0.01)
  }
  later::later(respond, 0.01)
  report <- cellucid_debug_connection(viewer, timeout = 10)
  expect_true(report$channel_connected)
  expect_true(report$frontend_roundtrip$ok)
  expect_identical(report$frontend_roundtrip$answer$t, 5L)
  expect_true(report$frontend_debug_snapshot$ok)
  expect_identical(report$frontend_debug_snapshot$answer$userAgent, "test")
  answering <- FALSE
  ws_close(client)
})

test_that("a logical vector of cells with NA is refused", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  cells <- rep(FALSE, 12L)
  cells[[1L]] <- TRUE
  cells[[3L]] <- NA
  expect_error(cellucid_highlight(viewer, cells), "holds NA at cell 3")
})

test_that("quiet = TRUE silences the notes reading the data writes", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  SummarizedExperiment::colData(object)$when <- Sys.Date() + seq_len(12L)
  expect_silent(viewer <- cellucid_serve(
    object, dataset_id = "s", dataset_name = "S", port = 0L, serve_web_ui = FALSE, quiet = TRUE
  ))
  cellucid_stop(viewer)
})

test_that("a Zarr store whose identity cannot be read is said to be damaged", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(file.path(tempfile(), "f.zarr"))
  folder <- suppressMessages(do.call(cellucid_prepare, args))
  writeLines("{ not json", file.path(folder, "uns/cellucid/id/.zarray"))
  expect_error(
    cellucid_serve(folder, port = 0L, serve_web_ui = FALSE, quiet = TRUE),
    "holds uns/cellucid/id and it cannot be read"
  )
  # Named by the caller, the store's own identity is not needed.
  viewer <- cellucid_serve(folder, dataset_id = "f", dataset_name = "F", port = 0L, serve_web_ui = FALSE, quiet = TRUE)
  expect_identical(viewer$dataset_id, "f")
  cellucid_stop(viewer)
})

test_that("unloading the package stops every viewer first", {
  skip_without_sce()
  viewer <- cellucid_serve(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S",
                           port = 0L, serve_web_ui = FALSE, quiet = TRUE)
  expect_true(viewer$running)
  # .onUnload() stops every viewer before it unloads the DLL; the stop is
  # called here alone, so the suite's copy of the DLL stays loaded.
  expect_match(paste(deparse(body(cellucid:::.onUnload)), collapse = "\n"), "^\\{\n\\s*\\.stop_viewers_on_unload\\(\\)")
  cellucid:::.stop_viewers_on_unload()
  expect_false(viewer$running)
  expect_length(cellucid_viewers(), 0L)
})

test_that("DESCRIPTION says the package views data and downloads the web app", {
  description <- read.dcf(file.path(cellucid_repository_root(), "DESCRIPTION"))
  text <- gsub("[[:space:]]+", " ", paste(description[, c("Title", "Description")], collapse = " "))
  expect_match(text, "View", fixed = TRUE)
  expect_match(text, "cellucid_view()", fixed = TRUE)
  expect_match(text, "cellucid_serve()", fixed = TRUE)
  expect_match(text, "downloaded from <https://www.cellucid.com> at first use", fixed = TRUE)
})
