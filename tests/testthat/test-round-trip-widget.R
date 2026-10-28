# The round trip in a Jupyter notebook (viewer-widget.R) against stand-in
# comms. IRkernel hands the widget the viewer's messages between cells, so
# nothing in the round trip may ask the viewer and wait: the page sends its
# snapshot (every result, then the inventory) and each change, and the next
# cell pulls, subsets, lists gene sets and pushes at once. A session is asked
# in one cell and collected in a later one.

round_trip_comms <- function() {
  comms <- new.env(parent = emptyenv())
  comms$handlers <- list()
  comms$sent <- list()
  wire <- function(value) jsonlite::fromJSON(as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null")))
  comms$api <- list(
    open = function(comm_id, data, on_msg) comms$handlers[[comm_id]] <- on_msg,
    send = function(comm_id, data) comms$sent[[length(comms$sent) + 1L]] <- c(list(comm_id = comm_id), wire(data)),
    display = function(bundle, metadata) invisible(NULL)
  )
  comms$post <- function(model_id, content) comms$handlers[[model_id]](wire(list(method = "custom", content = content)))
  comms$to <- function(t) Filter(function(m) identical(m$method, "custom") && identical(m$content$t, t), comms$sent)
  comms
}

round_trip_widget <- function(comms, env = parent.frame()) {
  args <- prepare_fixture(file.path(tempfile(), "f.zarr"))
  folder <- suppressMessages(do.call(cellucid_prepare, args))
  viewer <- cellucid:::.start_viewer(
    folder, list(), label = "folder",
    port = 0L, port_given = TRUE, host = "127.0.0.1", allowed_hosts = NULL, serve_web_ui = FALSE,
    web_cache_dir = NULL, web_source_url = NULL, cache_bytes = 268435456, quiet = TRUE, bind = FALSE
  )
  do.call(base::on.exit, list(bquote(cellucid_stop(.(viewer))), TRUE, FALSE), envir = env)
  cellucid:::.widget_show(viewer, "https://www.cellucid.com", comms = comms$api)
  viewer
}

page_text <- function(comms, viewer, value) {
  comms$post(viewer$widget$model_id, list(t = "text", v = "out-1", s = 0L,
    text = as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = NA))))
}

page_binary <- function(comms, viewer, frame) {
  comms$post(viewer$widget$model_id, list(t = "binary", v = "out-1", s = 0L, last = TRUE, b64 = jsonlite::base64_enc(frame)))
}

widget_result <- function(kind, id, meta, cells, values = NULL, type = NULL, revision) {
  header <- list(type = "result", requestId = NULL, epoch = "0123456789abcdef", revision = revision, kind = kind, id = id, nObs = 6, meta = meta,
                 cells = if (is.null(cells)) NULL else length(cells), values = type)
  text <- as.character(jsonlite::toJSON(header, auto_unbox = TRUE, null = "null", digits = NA))
  bytes <- charToRaw(text)
  c(cellucid:::.encode_uint32(length(bytes)), bytes,
    if (length(cells) > 0L) cellucid:::.encode_uint32(cells),
    if (!is.null(values)) cellucid:::.encode_wire_values(values, type))
}

test_that("a lasso, then cellucid_pull() in the next cell: nothing waits", {
  skip_if_not_installed("Matrix")
  comms <- round_trip_comms()
  viewer <- round_trip_widget(comms)
  cells <- data.frame(n = 1:6, row.names = paste0("cell-", 1:6))
  # Before the snapshot: what to do, at once.
  expect_error(cellucid_pull(cells, viewer), "has not sent what it holds yet")
  comms$post(viewer$widget$model_id, list(t = "open", v = "out-1", s = 0L))
  page_text(comms, viewer, list(type = "hello", viewerId = viewer$viewer_id, viewerToken = viewer$viewer_token))
  # Between cells: the page opened its dataset and sends its snapshot.
  page_text(comms, viewer, list(type = "ready", n_cells = 6, n_obs = 6, dimensions = 2))
  page_binary(comms, viewer, widget_result("highlight", "highlight_1",
    list(label = "Lasso", tool = "lasso", pageId = "page_1", pageName = "Page 1", color = "#e41a1c", enabled = TRUE),
    c(1, 4), revision = 1))
  page_binary(comms, viewer, widget_result("categorical", "user_cat_1",
    list(key = "annotation", categories = I(c("B", "T")), origin = "user"),
    NULL, c(0, 1, 255, 0, 1, 0), "uint8", revision = 2))
  page_binary(comms, viewer, widget_result("geneSet", "gene_set_1",
    list(name = "T", description = "", genes = I("CD3E"), symbols = I("CD3E"), origin = "markers"), numeric(0), revision = 3))
  page_text(comms, viewer, list(type = "results", epoch = "0123456789abcdef", dataset = viewer$dataset_id, revision = 3, nObs = 6, items = list(
    list(kind = "highlight", id = "highlight_1", revision = 1),
    list(kind = "categorical", id = "user_cat_1", revision = 2),
    list(kind = "geneSet", id = "gene_set_1", revision = 3)
  )))

  # The next cell: everything reads what arrived.
  elapsed <- system.time(pulled <- cellucid_pull(cells, viewer))[["elapsed"]]
  expect_lt(elapsed, 1)
  expect_identical(pulled$cellucid_highlight__highlight_1, c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_identical(as.character(pulled$annotation), c("B", "T", NA, "B", "T", "B"))
  expect_identical(cellucid_gene_sets(viewer), list(T = "CD3E"))
  expect_identical(rownames(cellucid_subset(cells, viewer, "Lasso")), c("cell-2", "cell-5"))

  # Pushing sends at once.
  cellucid_push(viewer, c(1.5, NA, 2, 3, 4, 5), key = "score")
  sent <- comms$to("binary")
  frame <- jsonlite::base64_dec(sent[[length(sent)]]$content$b64)
  header_length <- readBin(frame[1:4], "integer", size = 4L, endian = "little")
  expect_identical(jsonlite::fromJSON(rawToChar(frame[4L + seq_len(header_length)]))$key, "score")

  # A session: asked in one cell, collected in a later one.
  request <- cellucid_request_session(viewer)
  expect_error(cellucid_session_result(request), "Run cellucid_session_result\\(\\) again in a later cell")
  asked <- Filter(function(m) grepl("requestSessionBundle", m$content$text %||% ""), comms$to("text"))
  id <- jsonlite::fromJSON(asked[[1L]]$content$text)$requestId
  bundle <- c(charToRaw("CELLUCID_SESSION/2\n"), as.raw(c(2, 0, 0, 0)), charToRaw("{}"))
  head <- charToRaw(as.character(jsonlite::toJSON(list(type = "session_bundle", requestId = id), auto_unbox = TRUE)))
  page_binary(comms, viewer, c(cellucid:::.encode_uint32(length(head)), head, bundle))
  path <- cellucid_session_result(request)
  expect_identical(readBin(path, "raw", 1000L), bundle)
  expect_error(cellucid_session_result(request), "collected or withdrawn already")
})
