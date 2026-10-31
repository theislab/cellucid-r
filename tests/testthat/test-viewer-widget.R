# The viewer in a Jupyter notebook (viewer-widget.R) against stand-in comms:
# this test plays the kernel's comm manager and the widget's front end
# (inst/widget/widget.js), sending what the front end relays from the page
# and reading what R sends back, every message through the JSON IRkernel
# reads and writes. The routes and the channel behind them are the viewer's
# own (.route(), .on_ws_open()); no server is started and no port is bound.

# The kernel's comms, recording what R sends, and the front end's side.
fake_comms <- function() {
  comms <- new.env(parent = emptyenv())
  comms$opened <- list()
  comms$handlers <- list()
  comms$sent <- list()
  comms$displayed <- list()
  comms$display_metadata <- list()
  # IRkernel writes with jsonlite::toJSON(auto_unbox = TRUE) and reads with
  # jsonlite::fromJSON() (vectors simplified): every message goes through both.
  wire <- function(value) jsonlite::fromJSON(as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", json_verbatim = TRUE)))
  comms$api <- list(
    open = function(comm_id, data, on_msg) {
      comms$opened[[comm_id]] <- wire(data)
      comms$handlers[[comm_id]] <- on_msg
    },
    # The message's binary buffers travel beside its JSON, as raw vectors.
    send = function(comm_id, data, buffers = list()) {
      stopifnot(is.list(buffers), all(vapply(buffers, is.raw, TRUE)))
      comms$sent[[length(comms$sent) + 1L]] <- c(list(comm_id = comm_id), wire(data), list(buffers = buffers))
    },
    display = function(bundle, metadata) {
      comms$displayed[[length(comms$displayed) + 1L]] <- wire(bundle)
      comms$display_metadata[[length(comms$displayed)]] <- if (is.null(metadata)) NULL else wire(metadata)
    }
  )
  comms$post <- function(model_id, content) {
    comms$handlers[[model_id]](wire(list(method = "custom", content = content)))
  }
  # The custom messages R sent to output `view`, of kind `t`.
  comms$to <- function(view, t = NULL) {
    Filter(function(message) {
      identical(message$method, "custom") && identical(message$content$v, view) &&
        (is.null(t) || identical(message$content$t, t))
    }, comms$sent)
  }
  comms
}

# A viewer of a prepared store (6 cells, 4 genes) shown as a notebook
# widget, stopped when the calling test ends.
widget_viewer <- function(comms, env = parent.frame()) {
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

model_of <- function(viewer) viewer$widget$model_id

connect <- function(comms, viewer, view = "out-1", socket = 0L) {
  comms$post(model_of(viewer), list(t = "open", v = view, s = socket))
  comms$post(model_of(viewer), list(
    t = "text", v = view, s = socket,
    text = as.character(jsonlite::toJSON(list(type = "hello", viewerId = viewer$viewer_id, viewerToken = viewer$viewer_token), auto_unbox = TRUE))
  ))
  texts <- comms$to(view, "text")
  expect_identical(jsonlite::fromJSON(texts[[1L]]$content$text), list(type = "welcome"))
}

example_event <- function(type) {
  for (event in cellucid_read_fixture("viewer_protocol.json")$examples$events) {
    if (identical(event$type, type)) return(event)
  }
  stop("no example ", type)
}

test_that("the widget is an anywidget model of the shared front end, with no server and no credentials in it", {
  skip_if_not_installed("Matrix")
  comms <- fake_comms()
  viewer <- widget_viewer(comms)
  expect_null(viewer$server)
  expect_null(viewer$port)
  state <- comms$opened[[model_of(viewer)]]$state
  expect_identical(state$`_model_name`, "AnyModel")
  expect_identical(state$`_model_module`, "anywidget")
  expect_identical(state$`_binary`, "base64")
  expect_identical(state$web_source_url, "https://www.cellucid.com")
  expect_identical(state$`_esm`, paste(readLines(system.file("widget", "widget.js", package = "cellucid"), warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  layout <- sub("^IPY_MODEL_", "", state$layout)
  expect_identical(comms$opened[[layout]]$state$`_model_name`, "LayoutModel")
  view <- comms$displayed[[1L]]$`application/vnd.jupyter.widget-view+json`
  expect_identical(view$model_id, model_of(viewer))
  everything <- as.character(jsonlite::toJSON(list(comms$opened, comms$displayed), auto_unbox = TRUE))
  expect_false(grepl(viewer$viewer_token, everything, fixed = TRUE))
  # The credentials go only to the output that asks, in a message.
  comms$post(model_of(viewer), list(t = "attach", v = "out-1"))
  credentials <- comms$to("out-1", "credentials")[[1L]]$content
  expect_identical(credentials$viewerToken, viewer$viewer_token)
  expect_identical(credentials$viewerId, viewer$viewer_id)
  # A front end that reloaded is sent the model again.
  comms$handlers[[model_of(viewer)]](list(method = "request_state"))
  expect_identical(comms$sent[[length(comms$sent)]]$state$`_model_name`, "AnyModel")
})

test_that("the page reads the documents and the store through the widget, by range", {
  skip_if_not_installed("Matrix")
  comms <- fake_comms()
  viewer <- widget_viewer(comms)
  comms$post(model_of(viewer), list(t = "fetch", v = "out-1", requests = list(
    list(id = 1L, method = "GET", path = "_cellucid/health", range = NULL),
    list(id = 2L, method = "GET", path = "_cellucid/info", range = NULL),
    list(id = 3L, method = "GET", path = "_cellucid/datasets", range = NULL),
    list(id = 4L, method = "GET", path = "_cellucid/nothing", range = NULL)
  )))
  message <- comms$to("out-1", "answers")[[1L]]
  answers <- message$content$answers
  expect_identical(answers$status, c(200L, 200L, 200L, 404L))
  # Each body is a binary buffer of the message, named by its index; none is base64.
  expect_null(answers$b64)
  expect_identical(answers$buffer, 0:3)
  body <- function(row) jsonlite::fromJSON(rawToChar(message$buffers[[answers$buffer[[row]] + 1L]]), simplifyVector = FALSE)
  expect_identical(body(1)$n_cells, 6L)
  info <- body(2)
  expect_identical(info$runtime, "r")
  expect_null(info$port)
  expect_null(info$host)
  url <- body(3)$datasets[[1L]]$anndata$file$url
  expect_identical(url, "_cellucid/anndata/fixture.zarr/")
  root_key <- if (file.exists(file.path(viewer$served$path, "zarr.json"))) "zarr.json" else ".zmetadata"

  comms$post(model_of(viewer), list(t = "fetch", v = "out-1", requests = list(
    list(id = 5L, method = "GET", path = paste0(url, root_key), range = NULL),
    list(id = 6L, method = "GET", path = paste0(url, root_key), range = "bytes=0-9"),
    list(id = 7L, method = "HEAD", path = paste0(url, root_key), range = NULL)
  )))
  message <- comms$to("out-1", "answers")[[2L]]
  answers <- message$content$answers
  whole <- message$buffers[[1L]]
  expect_identical(answers$status, c(200L, 206L, 200L))
  expect_identical(answers$buffer, c(0L, 1L, -1L))
  expect_identical(message$buffers[[2L]], whole[1:10])
  expect_identical(answers$headers$`Content-Range`[[2L]], paste0("bytes 0-9/", length(whole)))
  # A HEAD has no body, and no buffer.
  expect_length(message$buffers, 2L)
})

test_that("the channel runs over the comm: commands first, then events and hooks", {
  skip_if_not_installed("Matrix")
  comms <- fake_comms()
  viewer <- widget_viewer(comms)
  cellucid_color_by(viewer, "cluster")
  seen <- list()
  cellucid_on(viewer, "selection", function(event) seen$selection <<- event)
  connect(comms, viewer)
  texts <- comms$to("out-1", "text")
  expect_identical(jsonlite::fromJSON(texts[[2L]]$content$text), list(type = "setColorBy", field = "cluster"))

  ready <- example_event("ready")
  comms$post(model_of(viewer), list(t = "text", v = "out-1", s = 0L, text = as.character(jsonlite::toJSON(ready, auto_unbox = TRUE))))
  expect_false(is.null(viewer$ready))
  # A binary frame in two pieces, as the front end splits one.
  header <- charToRaw("{\"type\":\"selection\",\"source\":\"lasso\"}")
  frame <- c(cellucid:::.encode_uint32(length(header)), header, cellucid:::.encode_uint32(c(0, 4, 5)))
  comms$post(model_of(viewer), list(t = "binary", v = "out-1", s = 0L, last = FALSE, b64 = jsonlite::base64_enc(frame[1:7])))
  comms$post(model_of(viewer), list(t = "binary", v = "out-1", s = 0L, last = TRUE, b64 = jsonlite::base64_enc(frame[8:length(frame)])))
  expect_identical(seen$selection$cells, c(1L, 5L, 6L))
  # A command that names cells goes as a binary frame: the message's buffer.
  cellucid_highlight(viewer, c(2L, 3L))
  binary <- comms$to("out-1", "binary")
  expect_length(binary, 1L)
  expect_null(binary[[1L]]$content$b64)
  expect_identical(binary_command(binary[[1L]]$buffers[[1L]])$cells, c(1L, 2L))
})

test_that("a cell cannot wait for the viewer, and says what to do instead", {
  skip_if_not_installed("Matrix")
  comms <- fake_comms()
  viewer <- widget_viewer(comms)
  expect_error(cellucid_wait_ready(viewer, timeout = 30), "only between cells")
  connect(comms, viewer)
  ready <- example_event("ready")
  comms$post(model_of(viewer), list(t = "text", v = "out-1", s = 0L, text = as.character(jsonlite::toJSON(ready, auto_unbox = TRUE))))
  # Once the viewer said it is ready (between cells), the answer is there.
  expect_false(is.null(cellucid_wait_ready(viewer, timeout = 30)))
})

test_that("a second output replaces the first, and stopping closes the channel and the reads", {
  skip_if_not_installed("Matrix")
  comms <- fake_comms()
  viewer <- widget_viewer(comms)
  connect(comms, viewer, view = "out-1", socket = 1L)
  connect(comms, viewer, view = "out-2", socket = 0L)
  closed <- comms$to("out-1", "close")
  expect_identical(closed[[1L]]$content$code, 4409L)
  cellucid_stop(viewer)
  expect_identical(comms$to("out-2", "close")[[1L]]$content$code, 1001L)
  comms$post(model_of(viewer), list(t = "fetch", v = "out-2", requests = list(list(id = 9L, method = "GET", path = "_cellucid/health", range = NULL))))
  expect_identical(comms$to("out-2", "answers")[[1L]]$content$answers$status, 410L)
})

test_that("a comm message with buffers is the kernel's own four parts, then each buffer as a frame", {
  # A stand-in for IRkernel's Kernel: its reply and signing functions.
  kernel <- list(
    new_reply = function(msg_type, parent) list(header = list(msg_type = msg_type), parent_header = parent$header, identities = list(charToRaw("id")), metadata = structure(list(), names = character())),
    sign_msg = function(parts) paste(vapply(parts, length, 0L), collapse = "-")
  )
  bodies <- list(as.raw(0:255), raw(3))
  for (data in list(
    list(method = "custom", content = list(t = "answers")),
    structure("{\"method\":\"custom\",\"content\":{\"t\":\"answers\"}}", class = "json")
  )) {
    parts <- cellucid:::.comm_msg_parts(kernel, list(header = list(msg_id = "p")), "comm-1", data, bodies)
    expect_length(parts, 9L)
    expect_identical(parts[[1L]], charToRaw("id"))
    expect_identical(parts[[2L]], charToRaw("<IDS|MSG>"))
    # Signed over header, parent header, metadata and content, as msg_to_wire signs.
    expect_identical(rawToChar(parts[[3L]]), paste(vapply(parts[4:7], length, 0L), collapse = "-"))
    expect_identical(parts[8:9], bodies)
    expect_identical(jsonlite::fromJSON(rawToChar(parts[[4L]]))$msg_type, "comm_msg")
    expect_identical(jsonlite::fromJSON(rawToChar(parts[[5L]]))$msg_id, "p")
    expect_identical(rawToChar(parts[[6L]]), "{}")
    content <- jsonlite::fromJSON(rawToChar(parts[[7L]]))
    expect_identical(content$comm_id, "comm-1")
    expect_identical(content$data$content$t, "answers")
  }
})

test_that("the answers' JSON is jsonlite's, every string escaped", {
  headers <- list(`Content-Type` = "application/json", ETag = "\"v-1\"", `X-Odd` = "back\\slash\ttab", `Content-Length` = 12L)
  json <- cellucid:::.widget_answer_json(17, 206L, headers, 3L)
  expect_identical(jsonlite::fromJSON(json, simplifyVector = FALSE), list(
    id = 17L, status = 206L,
    headers = list(`Content-Type` = "application/json", ETag = "\"v-1\"", `X-Odd` = "back\\slash\ttab", `Content-Length` = "12"),
    buffer = 3L
  ))
  expect_identical(jsonlite::fromJSON(cellucid:::.widget_answer_json(2^40, 404L, list(), -1L), simplifyVector = FALSE),
                   list(id = 2^40, status = 404L, headers = structure(list(), names = character()), buffer = -1L))
})

test_that("the widget has IRkernel take two fetch messages of the page at once", {
  comms <- fake_comms()
  viewer <- widget_viewer(comms)
  expect_identical(comms$opened[[viewer$widget$model_id]]$state$`_fetch_messages`, 2L)
})

test_that("the widget address is exact", {
  for (url in c("https://www.cellucid.com/", "ftp://x.org", "https://x.org/?a=1", "https://u:p@x.org", "x")) {
    expect_error(cellucid:::.validate_widget_source_url(url), "web_source_url")
  }
  expect_identical(cellucid:::.validate_widget_source_url("https://cellucid.example.org/app"), "https://cellucid.example.org/app")
})

# The widget's front end, byte for byte the Python package's: pinned here and
# in cellucid-python's tests/test_notebook_widget_contract.py.
.widget_js_sha256 <- "75258ef90ac240b3902dded04695a7b51efbd7755606076d598b6a64c49a466b"
.widget_css_sha256 <- "3b9acad3ebc04ece34acf580446a72ed7f73f3497b747539a11bc33fe281d82b"

test_that("the widget's front end is the Python package's, byte for byte", {
  expect_identical(cellucid_file_sha256(system.file("widget", "widget.js", package = "cellucid")), .widget_js_sha256)
  expect_identical(cellucid_file_sha256(system.file("widget", "widget.css", package = "cellucid")), .widget_css_sha256)
})

test_that("in a Jupyter notebook no server is started that the kernel could never run", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(file.path(tempfile(), "f.zarr"))
  folder <- suppressMessages(do.call(cellucid_prepare, args))
  withr_option <- options(jupyter.in_kernel = TRUE)
  on.exit(options(withr_option), add = TRUE)
  expect_error(cellucid_serve(folder, port = 0L, serve_web_ui = FALSE, quiet = TRUE), "cannot run a server between cells")
  expect_error(cellucid_view(folder, pane = "browser", quiet = TRUE), "the cell's output")
  expect_length(cellucid_viewers(), 0L)
  options(jupyter.in_kernel = FALSE)
  expect_error(cellucid_view(folder, pane = "notebook", quiet = TRUE), "IRkernel session")
  expect_length(cellucid_viewers(), 0L)
})

test_that("in Google Colab the display names Colab's widget manager, and nowhere else", {
  skip_if_not_installed("Matrix")
  # Outside Colab the display carries no metadata.
  withr_env <- Sys.getenv("COLAB_RELEASE_TAG", unset = NA)
  on.exit(if (is.na(withr_env)) Sys.unsetenv("COLAB_RELEASE_TAG") else Sys.setenv(COLAB_RELEASE_TAG = withr_env), add = TRUE)
  Sys.unsetenv("COLAB_RELEASE_TAG")
  comms <- fake_comms()
  widget_viewer(comms)
  expect_length(comms$display_metadata, 0L)
  # In Colab: exactly what google.colab's enable_custom_widget_manager() adds
  # to a Python widget's display (google/colab/output/_widgets.py).
  Sys.setenv(COLAB_RELEASE_TAG = "release-colab-test")
  comms <- fake_comms()
  widget_viewer(comms)
  expect_identical(comms$display_metadata[[1L]], list(
    `application/vnd.jupyter.widget-view+json` = list(colab = list(custom_widget_manager = list(
      url = "https://ssl.gstatic.com/colaboratory-static/widgets/colab-cdn-widget-manager/2b70e893a8ba7c0f/manager.min.js"
    )))
  ))
})
