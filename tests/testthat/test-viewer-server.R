# The viewer's server: the routes cellucid-python's server answers, the
# served store with exact Range answers, the Host and Origin checks, and the
# viewer channel (WebSocket) end to end. Every server listens on a free
# loopback port and is stopped when its test ends; nothing is downloaded (the
# fixture web build) and no browser is opened.

test_that("the data routes describe the dataset as cellucid-python's server does", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "sce", dataset_name = "SCE")
  health <- http_json(paste0(viewer$url, "/_cellucid/health"))
  expect_identical(health, list(status = "ok", type = "anndata", version = as.character(utils::packageVersion("cellucid")),
                                format = "zarr", n_cells = 12L, n_genes = 30L))
  info <- http_json(paste0(viewer$url, "/_cellucid/info"))
  expect_identical(info$port, viewer$port)
  expect_identical(info[c("format", "n_cells", "n_genes")], health[c("format", "n_cells", "n_genes")])
  expect_identical(
    names(info), c("version", "type", "format", "host", "port", "n_cells", "n_genes", "runtime")
  )
  # The program serving the viewer, by a name the shared protocol declares.
  expect_identical(info$runtime, "r")
  expect_true(info$runtime %in% names(cellucid_read_fixture("viewer_protocol.json")$runtimes))
  # The page address names a host program's viewer by the protocol's key.
  page <- cellucid_read_fixture("viewer_protocol.json")$page
  query <- strsplit(sub("^[^?]*\\?", "", viewer$viewer_url), "&", fixed = TRUE)[[1]]
  keys <- vapply(strsplit(query, "=", fixed = TRUE), `[[`, "", 1L)
  expect_identical(keys, as.character(unlist(page$parameters)))
  expect_identical(query[[1]], paste0(page$mode, "=true"))
  datasets <- http_json(paste0(viewer$url, "/_cellucid/datasets"))
  expect_identical(datasets, list(datasets = list(list(
    id = "sce", name = "SCE",
    anndata = list(file = list(url = "_cellucid/anndata/sce.zarr/", kind = "zarr", layout = "folder"), defaults = structure(list(), names = character()))
  ))))
  protocol <- http_json(paste0(viewer$url, "/_cellucid/protocol"))
  expect_identical(protocol, cellucid_read_fixture("viewer_protocol.json")$protocolDocument)
  expect_identical(http_request(paste0(viewer$url, "/_cellucid/nothing"))$status, 404L)
  expect_identical(http_request(paste0(viewer$url, "/_cellucid/health"), method = "POST")$status, 405L)
})

test_that("the web build is served from the verified cache, and nothing else of it", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "sce", dataset_name = "SCE")
  index <- http_request(paste0(viewer$url, "/"))
  expect_identical(index$status, 200L)
  expect_identical(index$headers[["content-type"]], "text/html; charset=utf-8")
  expect_match(rawToChar(index$body), "fixture-build-1", fixed = TRUE)
  script <- http_request(paste0(viewer$url, "/assets/app.js"))
  expect_identical(script$headers[["content-type"]], "text/javascript; charset=utf-8")
  expect_identical(http_request(paste0(viewer$url, "/assets/other.js"))$status, 404L)
  expect_identical(http_request(paste0(viewer$url, "/cellucid-web-assets.json"))$status, 404L)
  # A file changed on disk after the check is refused, not served.
  writeLines("tampered", file.path(viewer$web$dir, "assets", "app.js"))
  expect_identical(http_request(paste0(viewer$url, "/assets/app.js"))$status, 500L)
})

test_that("without the web build the data routes still answer", {
  skip_without_sce()
  viewer <- cellucid_serve(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S", port = 0L, serve_web_ui = FALSE, quiet = TRUE)
  on.exit(cellucid_stop(viewer))
  expect_identical(http_request(paste0(viewer$url, "/"))$status, 503L)
  expect_identical(http_request(paste0(viewer$url, "/_cellucid/health"))$status, 200L)
})

test_that("every key of the served store answers with the prepared store's arrays, with exact ranges", {
  skip_without_sce()
  object <- small_sce()
  # Served in priority order under the served seed, as cellucid_prepare() writes it by default.
  folder <- suppressMessages(cellucid_prepare(file.path(tempfile(), "s.zarr"), "zarr", data = object, dataset_id = "s", dataset_name = "S"))
  viewer <- test_viewer(object, dataset_id = "s", dataset_name = "S")
  served <- tempfile("cellucid_answered_")
  for (key in list.files(folder, recursive = TRUE, all.files = TRUE)) {
    response <- http_request(served_url(viewer, key))
    expect_identical(response$status, 200L, info = key)
    expect_identical(response$headers[["content-length"]], as.character(length(response$body)), info = key)
    dir.create(dirname(file.path(served, key)), recursive = TRUE, showWarnings = FALSE)
    writeBin(response$body, file.path(served, key))
  }
  expect_same_arrays(folder, served)
  key <- "mod/main/X/data/0"
  whole <- http_request(served_url(viewer, key))$body
  size <- length(whole)
  part <- http_request(served_url(viewer, key), list(Range = "bytes=3-9"))
  expect_identical(part$status, 206L)
  expect_identical(part$body, whole[4:10])
  expect_identical(part$headers[["content-range"]], paste0("bytes 3-9/", size))
  tail <- http_request(served_url(viewer, key), list(Range = "bytes=-4"))
  expect_identical(tail$body, whole[(size - 3):size])
  open <- http_request(served_url(viewer, key), list(Range = paste0("bytes=", size - 2, "-")))
  expect_identical(open$body, whole[(size - 1):size])
  expect_identical(http_request(served_url(viewer, key), list(Range = paste0("bytes=", size, "-")))$status, 416L)
  expect_identical(http_request(served_url(viewer, key), list(Range = "bytes=1-2,4-5"))$status, 416L)
  etag <- part$headers[["etag"]]
  expect_identical(http_request(served_url(viewer, key), list(`If-None-Match` = etag))$status, 304L)
  # A range of another version is answered whole.
  expect_identical(http_request(served_url(viewer, key), list(Range = "bytes=0-1", `If-Range` = "\"other\""))$status, 200L)
  expect_identical(http_request(served_url(viewer, "no/such/key"))$status, 404L)
  expect_identical(http_request(served_url(viewer, "../escape"))$status, 404L)
  expect_identical(http_request(paste0(viewer$url, "/_cellucid/anndata/other.zarr/.zmetadata"))$status, 404L)
  head <- http_request(served_url(viewer, ".zmetadata"), method = "HEAD")
  expect_identical(head$status, 200L)
})

test_that("chunks made by a transform are kept within the cache budget", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S", cache_bytes = 2000)
  keys <- grep("/[0-9.]+$", viewer$store$keys, value = TRUE)
  for (key in keys) http_request(served_url(viewer, key))
  cache <- viewer$store$cache
  expect_lte(cache$used, 2000)
  expect_gt(length(cache$sizes), 0L)
  expect_identical(cache$used, sum(cache$sizes))
  # Chunks of text and of a matrix are kept; a chunk of a numeric vector, a
  # slice of it, never is, nor X's gene-major arrays (laid out once, outside
  # the cache: .laid_out_arrays()).
  transformed <- vapply(keys, function(key) {
    meta <- jsonlite::fromJSON(rawToChar(cellucid:::.store_get(viewer$store, sub("/[^/]+$", "/.zarray", key))))
    identical(meta$dtype, "|O") || length(meta$shape) == 2L
  }, logical(1L))
  expect_true(any(!transformed))
  expect_true(all(names(cache$sizes) %in% keys[transformed]))
})

# The layouts' scratch files in R's session folder. Linux and macOS remove a
# scratch file's name as soon as it is mapped; Windows keeps the name of an
# open file and removes it when the file is closed (FILE_FLAG_DELETE_ON_CLOSE,
# src/chunks.c), which the viewer stopping does -- and the system does
# however R ends.
layout_scratch <- function() list.files(tempdir(), pattern = "^cellucid-layout-")

test_that("a served expression matrix is laid out once in the background, its chunks made from the matrix until then", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S", cache_bytes = 2000)
  store <- viewer$store
  # Laid out (test_viewer waits): X's chunks are slices of the layout, the
  # cache never holds them, and no scratch file is named in the folder.
  expect_identical(cellucid:::.native_layout_state(store$layouts, 0L)[[1L]], 2)
  if (.Platform$OS.type != "windows") expect_identical(layout_scratch(), character())
  laid <- lapply(c("X/data/0", "X/indices/0"), function(key) http_request(served_url(viewer, key)))
  expect_identical(vapply(laid, `[[`, 0L, "status"), c(200L, 200L))
  expect_false(any(startsWith(as.character(names(store$cache$sizes)), "X")))
  # While a matrix is laid out, a chunk (or a range of one) is made from the
  # matrix itself: the same bytes, never a wait.
  testthat::local_mocked_bindings(.native_layout_state = function(layouts, index) c(1, 250, 1000, 3), .package = "cellucid")
  during <- lapply(c("X/data/0", "X/indices/0"), function(key) http_request(served_url(viewer, key)))
  expect_identical(vapply(during, `[[`, 0L, "status"), c(200L, 200L))
  expect_identical(lapply(during, `[[`, "body"), lapply(laid, `[[`, "body"))
  ranged <- http_request(served_url(viewer, "X/indices/0"), list(Range = "bytes=4-11"))
  expect_identical(ranged$status, 206L)
  expect_identical(ranged$body, laid[[2L]]$body[5:12])
})

test_that("stopping a viewer stops its layouts and lets their files go", {
  skip_without_sce()
  before <- layout_scratch()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  layouts <- viewer$store$layouts
  # Its scratch files (two a matrix laid out): named only on Windows, while
  # they are open.
  named <- setdiff(layout_scratch(), before)
  expect_gt(length(viewer$store$layout_labels), 0L)
  expect_length(named, if (.Platform$OS.type == "windows") 2L * length(viewer$store$layout_labels) else 0L)
  cellucid_stop(viewer)
  expect_error(cellucid:::.native_layout_state(layouts, 0L), "stopped")
  expect_identical(intersect(layout_scratch(), named), character())
})

test_that("a Host that does not name the server is refused, and so is a foreign WebSocket", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  expect_identical(http_request(paste0(viewer$url, "/_cellucid/health"), list(Host = "attacker.example"))$status, 421L)
  expect_identical(http_request(paste0(viewer$url, "/_cellucid/health"), list(Host = "127.0.0.1:1"))$status, 421L)
  expect_identical(http_request(paste0(viewer$url, "/_cellucid/health"), list(Host = paste0("localhost:", viewer$port)))$status, 200L)
  # Each refused socket is closed by the server as httpuv opens it, never
  # left to say hello (ws_refused_close(), helper-viewer.R).
  foreign <- ws_open(viewer, origin = "https://attacker.example")
  expect_identical(foreign$status, 403L)
  expect_identical(ws_refused_close(foreign)$code, 1008L)
  ws_close(foreign)
  anonymous <- ws_open(viewer, origin = NULL)
  expect_identical(anonymous$status, 403L)
  expect_identical(ws_refused_close(anonymous)$code, 1008L)
  ws_close(anonymous)
  rebound <- ws_open(viewer, host = "attacker.example")
  expect_identical(rebound$status, 421L)
  expect_identical(ws_refused_close(rebound)$code, 1008L)
  ws_close(rebound)
  expect_null(viewer$socket)
})

test_that("cross-origin reads are open to the hosted app and loopback pages only, as in Python", {
  allowed <- function(origin) {
    headers <- cellucid:::.cors_headers(list(HTTP_ORIGIN = origin))
    identical(headers[["Access-Control-Allow-Origin"]], origin)
  }
  for (origin in c("https://www.cellucid.com", "https://cellucid.com", "http://localhost:3000", "http://127.0.0.1:8000")) {
    expect_true(allowed(origin), info = origin)
  }
  # No other page on a cellucid.com name: the documentation site runs its
  # host's scripts and must not read a served dataset.
  for (origin in c("https://docs.cellucid.com", "https://preview.cellucid.com", "http://www.cellucid.com",
                   "https://www.cellucid.com.example.org", "https://theislab.github.io", "null")) {
    expect_false(allowed(origin), info = origin)
  }
})

test_that("allowed_hosts lets a proxy's name through", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S", allowed_hosts = "rstudio.example.org")
  expect_identical(http_request(paste0(viewer$url, "/_cellucid/health"), list(Host = "rstudio.example.org:443"))$status, 200L)
  expect_error(cellucid_serve(small_sce(), port = 0L, allowed_hosts = "https://x.org", quiet = TRUE), "bare host names")
})

test_that("the channel welcomes the viewer's credentials and refuses others", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  wrong <- ws_open(viewer)
  ws_send_json(wrong, list(type = "hello", viewerId = viewer$viewer_id, viewerToken = strrep("0", 32)))
  expect_identical(ws_receive(wrong)$code, 4401L)
  ws_close(wrong)
  rude <- ws_open(viewer)
  ws_send_json(rude, list(type = "ready", n_cells = 1, n_obs = 1, dimensions = 2))
  expect_identical(ws_receive(rude)$code, 1008L)
  ws_close(rude)
  elsewhere <- ws_open(viewer, path = "/elsewhere")
  expect_identical(ws_receive(elsewhere)$code, 1008L)
  ws_close(elsewhere)

  first <- ws_welcomed(viewer)
  second <- ws_welcomed(viewer)
  expect_identical(ws_receive(first)$code, 4409L)
  ws_close(first)
  # A frame that breaks the protocol closes the channel with 1008.
  ws_send_text(second, "{\"type\":\"ready\"}")
  expect_identical(ws_receive(second)$code, 1008L)
  ws_close(second)
})

test_that("actions wait for the channel and go out in order, each answered or refused by its request", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  viewer <- test_viewer(object, dataset_id = "s", dataset_name = "S")
  highlight <- cellucid_highlight_cells(viewer, c("cell3", "cell1"), color = "#123456")
  color <- cellucid_color_by(viewer, field = list(source = "obs", key = "cluster"))
  hide <- cellucid_cells_visibility_set(viewer, c(TRUE, FALSE, rep(TRUE, 10)), visible = FALSE)
  reset <- cellucid_camera_reset(viewer, "live")
  expect_identical(highlight$state, "waiting for the viewer")
  client <- ws_welcomed(viewer)
  frames <- lapply(1:4, function(i) ws_receive(client))
  received <- lapply(frames, function(frame) frame$json)
  # An action that names cells is a binary frame; every other one is text.
  expect_identical(vapply(frames, function(frame) frame$opcode, 0L), c(2L, 1L, 2L, 1L))
  expect_identical(vapply(received, function(x) x$type, ""), c("action", "action", "action", "action"))
  expect_identical(vapply(received[1:4], function(x) x$action, ""), c("highlight.cells", "color.by", "cells.visibility.set", "camera.reset"))
  # Cells by original position, 0-based, in the order given.
  expect_identical(received[[1L]]$params, list(color = "#123456", cells = c(2, 0)))
  expect_identical(received[[2L]]$params, list(field = list(source = "obs", key = "cluster")))
  expect_identical(received[[3L]]$params, list(visible = FALSE, cells = c(0, 2:11)))
  expect_identical(received[[4L]]$params, list(view = "live"))
  expect_identical(highlight$state, "sent")
  expect_identical(received[[1L]]$requestId, highlight$id)
  # Once welcomed, an action goes out at once.
  everyone <- cellucid_cells_visibility_set(viewer, NULL, visible = TRUE)
  expect_identical(ws_receive(client)$json$params, list(cells = NULL, visible = TRUE))
  single <- cellucid_highlight_cells(viewer, 12)
  expect_identical(ws_receive(client)$json$params$cells, 11)
  expect_error(cellucid_highlight_cells(viewer, 13), "1-based positions from 1 to 12")
  expect_error(cellucid_highlight_cells(viewer, "nope"), "names cells the dataset does not have")
  expect_error(cellucid_highlight_cells(viewer, c(2, 2)), "names the cell 2 twice")

  # The viewer answers each: an outcome and the step it made, or a refusal
  # naming its request.
  answer <- function(id, outcome, step) {
    header <- charToRaw(as.character(jsonlite::toJSON(list(type = "action_result", requestId = id, step = step, unrecorded = NULL), auto_unbox = TRUE, null = "null")))
    body <- charToRaw(as.character(jsonlite::toJSON(outcome, auto_unbox = TRUE, null = "null")))
    ws_send_binary(client, c(cellucid:::.encode_uint32(length(header)), header, body))
  }
  answer(highlight$id, list(group = "highlight_1", page = "Page 1", cells = 2L), list(id = "s1", label = "Highlighted 2 cells"))
  answer(color$id, NULL, list(id = "s2", label = "Colored by cluster"))
  ws_send_json(client, list(type = "command_error", command = "action", requestId = hide$id,
                            reason = "There is no field \"nope\"."))
  seen <- list()
  cellucid_on(viewer, "action_result", function(event) seen[[length(seen) + 1L]] <<- event)
  expect_warning(
    run_loop_until(function() hide$done),
    "The viewer did not take cellucid_cells_visibility_set\\(\\): There is no field"
  )
  expect_identical(cellucid_action_result(highlight, timeout = 5), list(group = "highlight_1", page = "Page 1", cells = 2L))
  expect_identical(highlight$step, list(id = "s1", label = "Highlighted 2 cells"))
  expect_null(cellucid_action_result(color, timeout = 5))
  expect_identical(color$state, "answered")
  expect_error(cellucid_action_result(hide, timeout = 5), "The viewer did not take cellucid_cells_visibility_set\\(\\): There is no field")
  expect_identical(hide$state, "refused")
  expect_output(print(highlight), "highlight.cells .*answered, step s1: Highlighted 2 cells")
  # A refusal the call that asked waits for is said by that call alone.
  later::later(function() {
    ws_send_json(client, list(type = "command_error", command = "action", requestId = reset$id, reason = "No view."))
  }, 0.1)
  expect_no_warning(expect_error(cellucid_action_result(reset, timeout = 10), "did not take cellucid_camera_reset\\(\\): No view"))
  expect_length(seen, 2L)
  # The answer is kept while the request is; an answer to a request nobody
  # holds any more is let go as it arrives.
  orphan <- everyone$id
  rm(everyone)
  invisible(gc())
  expect_true(viewer$action_requests[[orphan]]$orphan)
  answer(orphan, NULL, NULL)
  answer(single$id, list(group = "highlight_2"), list(id = "s3", label = "Highlighted 1 cell"))
  run_loop_until(function() is.null(viewer$action_requests[[orphan]]) && single$done)
  expect_identical(cellucid_action_result(single), list(group = "highlight_2"))
  expect_identical(cellucid_action_result(single), list(group = "highlight_2"))
  ws_close(client)
})

test_that("events reach the hooks, a selection by index and by name", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  viewer <- test_viewer(object, dataset_id = "s", dataset_name = "S")
  seen <- list()
  cellucid_on(viewer, "selection", function(event) seen$selection <<- event)
  cellucid_on(viewer, "hover", function(event) seen$hover <<- event)
  cellucid_on(viewer, "click", function(event) stop("a broken hook"))
  client <- ws_welcomed(viewer)
  ws_send_json(client, list(type = "ready", n_cells = 12, n_obs = 12, dimensions = 2))
  ws_send_json(client, list(type = "hover", cell = 3, position = list(x = 0.5, y = -1, z = 0)))
  header <- charToRaw("{\"type\":\"selection\",\"source\":\"lasso\"}")
  ws_send_binary(client, c(cellucid:::.encode_uint32(length(header)), header, cellucid:::.encode_uint32(c(0, 4, 11))))
  run_loop_until(function() !is.null(seen$selection))
  expect_identical(seen$selection$cells, c(1L, 5L, 12L))
  expect_identical(seen$selection$cell_names, colnames(object)[c(1, 5, 12)])
  expect_identical(cellucid_selection(viewer)$source, "lasso")
  # Cells are R's 1-based positions in every event (the wire's 3 is cell 4).
  expect_identical(seen$hover$cell, 4L)
  expect_identical(seen$hover$cell_name, colnames(object)[[4L]])
  expect_identical(cellucid_wait_ready(viewer, timeout = 5), list(n_cells = 12L, n_obs = 12L, dimensions = 2L))
  ws_send_json(client, list(type = "click", cell = 1, button = 0, shift = FALSE, ctrl = FALSE))
  expect_warning(run_loop_until(function() !is.null(viewer$last$click)), "a broken hook")
  cellucid_off(viewer, "selection")
  expect_null(viewer$hooks$selection)
  ws_close(client)
})

test_that("a large selection is decoded at once", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  # A viewer that knows neither its cells' names nor their count (a file
  # only the browser reads), so a million indices are all in range.
  viewer$cell_names <- NULL
  viewer$n_cells <- NULL
  client <- ws_welcomed(viewer)
  header <- charToRaw("{\"type\":\"selection\",\"source\":\"lasso\"}")
  cells <- seq(0, 1e6 - 1)
  ws_send_binary(client, c(cellucid:::.encode_uint32(length(header)), header, cellucid:::.encode_uint32(cells)))
  run_loop_until(function() !is.null(viewer$selection), timeout = 30)
  expect_identical(length(viewer$selection$cells), 1000000L)
  expect_identical(viewer$selection$cells[[1000000L]], 1000000L)
  ws_close(client)
})

test_that("a session bundle arrives as the session file it asked for, read and checked", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  client <- ws_welcomed(viewer)
  # Unasked bundles are not kept.
  unasked <- charToRaw("{\"type\":\"session_bundle\",\"requestId\":\"nobody\"}")
  ws_send_binary(client, c(cellucid:::.encode_uint32(length(unasked)), unasked, charToRaw("CELLUCID")))
  suppressWarnings(for (i in 1:20) later::run_now(0.02))
  expect_null(viewer$last$session_bundle)
  fixture <- testthat::test_path("fixtures", "history_golden", "history_golden.cellucid-session")
  send <- function(body) {
    function() {
      command <- ws_receive(client)$json
      header <- charToRaw(paste0("{\"type\":\"session_bundle\",\"requestId\":\"", command$requestId, "\"}"))
      ws_send_binary(client, c(cellucid:::.encode_uint32(length(header)), header, body))
    }
  }
  later::later(send(readBin(fixture, "raw", file.size(fixture))), 0.2)
  target <- tempfile(fileext = ".cellucid-session")
  session <- cellucid_session_result(cellucid_request_session(viewer), file = target, timeout = 10)
  expect_s3_class(session, "cellucid_session")
  expect_identical(session$path, normalizePath(target, winslash = "/"))
  expect_identical(cellucid_file_sha256(target), cellucid_file_sha256(fixture))
  expect_identical(cellucid_trackers(session)$name, c("Annotation", "Annotation (replay)", "Annotation"))
  expect_length(viewer$session_requests, 0L)
  # A file of another session format is reported to the one waiting for it,
  # in the viewer's own words.
  later::later(send(charToRaw("CELLUCID_SESSION/3\nbody")), 0.2)
  expect_no_warning(expect_error(
    cellucid_session_result(cellucid_request_session(viewer), timeout = 10),
    "invalid session file: This session was saved by an earlier version of Cellucid \\(session format 3\\).*Ask again with cellucid_request_session\\(\\)"
  ))
  # One nobody waits for is said as it arrives (as every refused command is).
  request <- cellucid_request_session(viewer)
  send(charToRaw("not a session"))()
  expect_warning(
    run_loop_until(function() !is.null(viewer$session_requests[[request$id]]$failed)),
    "The viewer did not take cellucid_request_session\\(\\): The viewer sent an invalid session file: This file is not a Cellucid session file of format 4"
  )
  expect_error(cellucid_session_result(request, timeout = 1), "not a Cellucid session file of format 4")
  ws_close(client)
})

test_that("waiting on a viewer that never answers stops with what to do", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  request <- cellucid_request_session(viewer)
  error <- expect_error(cellucid_session_result(request, timeout = 0.3), "did not answer within 0.3 seconds.*No viewer is connected: open the viewer's address \\(viewer\\$viewer_url\\)\\.$")
  # An error is logged and pasted: it never carries the viewer's token.
  expect_false(grepl(viewer$viewer_token, conditionMessage(error), fixed = TRUE))
  # The request stays open: collected later, once a viewer answers it.
  expect_length(viewer$session_requests, 1L)
})

test_that("stopping freezes the page, closes the channel and frees the port for a new viewer", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  viewer <- cellucid_serve(object, dataset_id = "s", dataset_name = "S", port = 0L, web_cache_dir = viewer_fixture_web_cache(), web_source_url = NULL, quiet = TRUE)
  port <- viewer$port
  client <- ws_welcomed(viewer)
  expect_true(viewer$viewer_id %in% vapply(cellucid_viewers(), function(v) v$viewer_id, ""))
  cellucid_stop(viewer)
  # The open page is told to freeze (it keeps its last view), as Python's stop() tells it.
  expect_identical(ws_receive(client)$json$type, "freeze")
  expect_identical(ws_receive(client)$code, 1001L)
  ws_close(client)
  expect_false(viewer$running)
  expect_error(cellucid_highlight_cells(viewer, 1), "This viewer is stopped")
  cellucid_stop(viewer)
  again <- cellucid_serve(object, dataset_id = "s", dataset_name = "S", port = port, web_cache_dir = viewer_fixture_web_cache(), web_source_url = NULL, quiet = TRUE)
  on.exit(cellucid_stop(again))
  expect_identical(again$port, port)
  expect_identical(http_request(paste0(again$url, "/_cellucid/health"))$status, 200L)
})

test_that("several viewers run at once, a taken port is refused or passed over", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  first <- test_viewer(object, dataset_id = "one", dataset_name = "One")
  second <- test_viewer(object, dataset_id = "two", dataset_name = "Two")
  expect_false(identical(first$port, second$port))
  expect_identical(http_json(paste0(first$url, "/_cellucid/datasets"))$datasets[[1L]]$id, "one")
  expect_identical(http_json(paste0(second$url, "/_cellucid/datasets"))$datasets[[1L]]$id, "two")
  expect_error(
    cellucid_serve(object, port = first$port, web_cache_dir = viewer_fixture_web_cache(), web_source_url = NULL, quiet = TRUE),
    paste0(
      "^Port ", first$port, " is already in use on 127.0.0.1, so Cellucid could not start its server there\\. ",
      ".*Serve on a different port with port = 9000, or let the operating system pick a free one with port = 0"
    )
  )
  # The default port, taken, is passed over for the next free one.
  blocker <- tryCatch(httpuv::startServer("127.0.0.1", 8765L, list(call = function(req) list(status = 200L, headers = list(), body = ""))),
                      error = function(error) NULL)
  if (!is.null(blocker)) on.exit(blocker$stop(), add = TRUE)
  # Every message captured: the address it says carries the viewer's token,
  # which test output (a CI log) never shows.
  said <- capture_messages(
    third <- cellucid_serve(object, dataset_id = "three", dataset_name = "Three", web_cache_dir = viewer_fixture_web_cache(), web_source_url = NULL)
  )
  on.exit(cellucid_stop(third), add = TRUE)
  expect_match(said[[1L]], "Port 8765 is in use; serving on [0-9]+")
  # The person is told the address to open, once, token and all.
  expect_identical(sum(grepl(third$viewer_token, said, fixed = TRUE)), 1L)
  expect_gt(third$port, 8765L)
  expect_identical(cellucid_stop_all() >= 3L, TRUE)
  expect_length(cellucid_viewers(), 0L)
})

test_that("a file rewritten while it is served is a new version, in place or replaced", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(file.path(tempfile(), "f.zarr.zip"))
  args$format <- "zarr.zip"
  zip <- suppressMessages(do.call(cellucid_prepare, args))
  viewer <- test_viewer(zip)
  url <- paste0(viewer$url, "/_cellucid/anndata/", viewer$served$route_name)
  listed <- function() http_json(paste0(viewer$url, "/_cellucid/datasets"))$datasets[[1L]]$anndata$file
  first <- listed()
  expect_identical(http_request(url, list(`If-None-Match` = first$etag))$status, 304L)
  expect_identical(listed()$etag, first$etag)

  # Rewritten in place, the same size, its modification time put back (as a
  # copy that keeps it would): still a new version.
  bytes <- readBin(zip, "raw", file.size(zip))
  modified <- file.mtime(zip)
  Sys.sleep(0.05)
  connection <- file(zip, "r+b")
  seek(connection, 0L, rw = "write")
  writeBin(rev(bytes[1:8]), connection)
  close(connection)
  Sys.setFileTime(zip, modified)
  second <- listed()
  expect_identical(second$size, first$size)
  expect_false(identical(second$etag, first$etag))
  stale <- http_request(url, list(`If-None-Match` = first$etag))
  expect_identical(stale$status, 200L)
  expect_identical(stale$headers[["etag"]], second$etag)
  # A range asked of the old version is answered whole, as the new one.
  ranged <- http_request(url, list(Range = "bytes=0-7", `If-Range` = first$etag))
  expect_identical(ranged$status, 200L)
  expect_identical(ranged$body[1:8], rev(bytes[1:8]))

  # Replaced by another file moved over it.
  Sys.sleep(0.05)
  other <- tempfile(tmpdir = dirname(zip))
  writeBin(bytes, other)
  expect_true(file.rename(other, zip))
  expect_false(listed()$etag %in% c(first$etag, second$etag))

  # A file that is gone is not listed as if it were there.
  unlink(zip)
  gone <- http_request(paste0(viewer$url, "/_cellucid/datasets"))
  expect_identical(gone$status, 404L)
  expect_identical(rawToChar(gone$body), paste0(basename(zip), " is no longer readable"))
  expect_identical(http_request(url)$status, 404L)
})

test_that("a served file's validator is read from the file opened, never from its path", {
  # A name outside every single-byte code page: on Windows the file is opened by its UTF-16 name.
  path <- file.path(tempfile("cellucid_version_"), "café α.h5ad")
  dir.create(dirname(path))
  writeBin(as.raw(0:255), path)
  version <- cellucid:::.file_version(path)
  expect_identical(version$size, 256)
  # The size, modification time, status-change time and identity, as hex.
  expect_match(version$etag, "^\"100-[0-9a-f]+-[0-9a-f]+-[0-9a-f]+\"$")
  expect_identical(cellucid:::.file_version(path), version)
  # What a path's file.info() reports is not the validator (on Windows its
  # ctime is the creation time, which a rewrite in place keeps): times it
  # reports otherwise change nothing.
  file_info <- base::file.info
  testthat::local_mocked_bindings(
    file.info = function(...) {
      info <- file_info(...)
      info$mtime <- info$mtime - 1
      info$ctime <- info$ctime - 1
      info
    },
    .package = "base"
  )
  expect_identical(cellucid:::.file_version(path), version)
  # A copy with the same bytes, its modification time copied (cp -p), is
  # another file.
  copy <- file.path(dirname(path), "copy.h5ad")
  expect_true(file.copy(path, copy, copy.date = TRUE))
  expect_false(identical(cellucid:::.file_version(copy)$etag, version$etag))
  # A folder, or a file that is gone, has no version.
  expect_null(cellucid:::.file_version(dirname(path)))
  unlink(copy)
  expect_null(cellucid:::.file_version(copy))
})

test_that("a prepared store or file given by path is served from disk", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(file.path(tempfile(), "f.zarr"))
  folder <- suppressMessages(do.call(cellucid_prepare, args))
  args$out_path <- sub("\\.zarr$", ".zarr.zip", args$out_path)
  args$format <- "zarr.zip"
  zip <- suppressMessages(do.call(cellucid_prepare, args))

  viewer <- test_viewer(zip)
  expect_identical(viewer$served[c("kind", "layout")], list(kind = "zarr", layout = "zip"))
  url <- paste0(viewer$url, "/_cellucid/anndata/", viewer$served$route_name)
  datasets <- http_json(paste0(viewer$url, "/_cellucid/datasets"))$datasets[[1L]]
  expect_identical(datasets$anndata$file$size, as.integer(file.size(zip)))
  # The store's own identity names it.
  expect_identical(datasets[c("id", "name")], list(id = "fixture", name = "Fixture"))
  expect_identical(datasets$anndata$file$url, "_cellucid/anndata/fixture.zarr.zip")
  whole <- http_request(url)
  expect_identical(whole$body, readBin(zip, "raw", file.size(zip)))
  part <- http_request(url, list(Range = "bytes=100-163"))
  expect_identical(part$body, whole$body[101:164])
  # The list states the validator the bytes are answered with: the version
  # the viewer files what it keeps of the file under.
  expect_match(datasets$anndata$file$etag, '^"[0-9a-f]+-[0-9a-f]+-[0-9a-f]+-[0-9a-f]+"$')
  expect_identical(whole$headers[["etag"]], datasets$anndata$file$etag)

  viewer <- test_viewer(folder, latent_key = "X_pca", obs_keys = c("cluster"))
  health <- http_json(paste0(viewer$url, "/_cellucid/health"))
  expect_identical(health[c("format", "n_cells", "n_genes")], list(format = "zarr", n_cells = 6L, n_genes = 4L))
  expect_identical(
    http_json(paste0(viewer$url, "/_cellucid/datasets"))$datasets[[1L]]$anndata$defaults,
    list(latentKey = "X_pca", obsKeys = list("cluster"))
  )
  # The neighbor graph and the expression matrix a file opens on, in the
  # viewer's spelling, as cellucid-python's serve_anndata(knn_key=,
  # expression_source=) states them.
  viewer <- test_viewer(folder, knn_key = "connectivities", expression_source = "layer:counts")
  expect_identical(
    http_json(paste0(viewer$url, "/_cellucid/datasets"))$datasets[[1L]]$anndata$defaults,
    list(knn = "connectivities", expression = "layers/counts")
  )
  for (case in list(c("X", "X"), c("raw", "raw/X"), c("layer:counts", "layers/counts"))) {
    expect_identical(cellucid:::.served_defaults(list(expression_source = case[[1L]]))$expression, case[[2L]])
  }
  for (bad in list("raw.X", "layer:", "layer:a/b", "layers/counts", NA_character_, 1)) {
    expect_error(cellucid:::.served_defaults(list(expression_source = bad)), "^expression_source must be \"X\", \"raw\" or \"layer:<name>\" \\(a layer's name, without \"/\"\\); got ")
  }
  expect_error(cellucid:::.served_defaults(list(knn_key = "")), "^knn_key must name a key")
  expect_error(cellucid:::.served_defaults(list(knn_key = c("a", "b"))), "^knn_key must be exactly one string")
  key <- http_request(served_url(viewer, "X/indptr/0"))
  expect_identical(key$body, readBin(file.path(folder, "X/indptr/0"), "raw", file.size(file.path(folder, "X/indptr/0"))))
  expect_identical(http_request(served_url(viewer, "../f.zarr.zip"))$status, 404L)

  # A folder without consolidated metadata gets it built in memory.
  unlink(file.path(folder, ".zmetadata"))
  viewer <- test_viewer(folder)
  consolidated <- jsonlite::fromJSON(rawToChar(http_request(served_url(viewer, ".zmetadata"))$body), simplifyVector = FALSE)
  expect_identical(consolidated$zarr_consolidated_format, 1L)
  expect_true(all(c(".zgroup", "obs/.zattrs", "X/data/.zarray") %in% names(consolidated$metadata)))
  expect_error(cellucid_serve(folder, latent = "x", quiet = TRUE), "For a path, cellucid_view\\(\\) takes")
})

test_that("arguments that do not choose what is shown are refused by name", {
  skip_without_sce()
  expect_error(cellucid_serve(small_sce(), zarr_format = 3, quiet = TRUE), "has no argument c\\(\"zarr_format\"\\)")
  expect_error(cellucid_serve(small_sce(), "positional", quiet = TRUE), "must be named")
  expect_error(cellucid_view(), "needs data")
  expect_error(cellucid_serve(small_sce(), port = 70000, quiet = TRUE), "port must be one whole number")
  # The seed of a served object's priority order is fixed (0, a prefix for the default subsample).
  expect_error(cellucid_serve(small_sce(), priority_seed = 3, quiet = TRUE), "has no argument c\\(\"priority_seed\"\\)")
})

test_that("a served object takes total_counts, the column of its cells' totals, as cellucid_prepare() does", {
  skip_without_sce()
  sce <- small_sce()
  adt <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, "ADT"))
  SummarizedExperiment::colData(sce)$adt_total <- as.numeric(Matrix::colSums(adt)) + 1
  viewer <- test_viewer(sce, total_counts = c(ADT = "adt_total"))
  expect_true("total_counts" %in% cellucid:::.view_plan_arguments)
  expect_error(cellucid_serve(sce, total_counts = c(ADT = "nowhere"), quiet = TRUE), "which obs does not hold")
})

test_that("a file served as it is says what a subsample of it reads, unless it is in priority order", {
  skip_if_not_installed("Matrix")
  # An HDF5 file (its signature), which R cannot read to put in priority order.
  loom <- tempfile(fileext = ".loom")
  writeBin(c(as.raw(c(0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a)), raw(4000L)), loom)
  note <- .path_order_note(.classify_path(loom))
  expect_match(note, "is served as it is. Its cells are in their own order:", fixed = TRUE)
  expect_match(note, "cellucid_view(object)", fixed = TRUE)
  expect_match(note, "cellucid serve <path>", fixed = TRUE)
  expect_no_match(note, "unless cellucid_prepare", fixed = TRUE)
  h5ad <- tempfile(fileext = ".h5ad")
  file.copy(loom, h5ad)
  expect_match(.path_order_note(.classify_path(h5ad)), "unless cellucid_prepare() or prepare() wrote it", fixed = TRUE)
  # A store cellucid_prepare() wrote in priority order: nothing to say, folder or zip.
  args <- prepare_fixture(file.path(tempfile(), "f.zarr"))
  args$priority_seed <- 0
  folder <- suppressMessages(do.call(cellucid_prepare, args))
  expect_null(.path_order_note(.classify_path(folder)))
  args$out_path <- sub("\\.zarr$", ".zarr.zip", args$out_path)
  args$format <- "zarr.zip"
  expect_null(.path_order_note(.classify_path(suppressMessages(do.call(cellucid_prepare, args)))))
  # One in its own order (no uns['cellucid']['priority']) is told.
  unlink(file.path(folder, "uns", "cellucid", "priority"), recursive = TRUE)
  expect_match(.path_order_note(.classify_path(folder)), "no reader of Zarr chunks", fixed = TRUE)
  # Said when it is served (not when quiet).
  messages <- character()
  viewer <- withCallingHandlers(
    cellucid_serve(loom, dataset_id = "l", dataset_name = "L", port = 0L, serve_web_ui = FALSE),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  )
  on.exit(cellucid_stop(viewer), add = TRUE)
  expect_true(any(grepl("Its cells are in their own order", messages, fixed = TRUE)))
})

test_that("a host this machine cannot serve from is refused saying so, not as a port taken", {
  skip_without_sce()
  # cellucid-python's bind_refusal() sentence, with R's spelling.
  for (host in c("no-such-host.invalid", "192.0.2.1")) {
    expect_error(
      cellucid_serve(small_sce(alt = FALSE), dataset_id = "h", dataset_name = "H", host = host,
                     port = 0L, web_cache_dir = viewer_fixture_web_cache(), web_source_url = NULL, quiet = TRUE),
      paste0("Host '", host, "' is not an address this machine can serve from. Use host = \"127.0.0.1\" to serve ",
             "only this computer, which is the default, or host = \"0.0.0.0\" to accept connections from other machines."),
      fixed = TRUE
    )
  }
})

test_that("zarr-python read every pair of stores this file compared", {
  expect_zarr_python_compared()
})
