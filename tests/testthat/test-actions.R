# The viewer's actions from R: the catalog the package carries (pinned in
# the three repositories), the functions and help pages generated from it
# (tools/generate-actions.R), and how R values become an action's parameters.
# The viewer is played here by the frames R queues before its page connects
# (R/viewer-server.R .viewer_send_frame()), read as the viewer reads them.

.action_catalog_sha256 <- "725d32e5357cebbcff332be5a599667cf42ebb63ce55bd73e39ce4a1469556c4"

catalog_path <- function() file.path(cellucid_repository_root(), "inst", "action_catalog.json")
read_catalog <- function() jsonlite::read_json(catalog_path(), simplifyVector = FALSE)

# A viewer that is running but whose page has not connected: what R sends
# waits in its queue, in order.
queued_viewer <- function(n_cells = 6L, cell_names = paste0("c", seq_len(n_cells))) {
  viewer <- new.env(parent = emptyenv())
  class(viewer) <- "cellucid_viewer"
  viewer$running <- TRUE
  viewer$queue <- list()
  viewer$action_requests <- list()
  viewer$session_requests <- list()
  viewer$hooks <- list()
  viewer$last <- list()
  viewer$n_cells <- n_cells
  viewer$cell_names <- cell_names
  viewer
}

last_sent <- function(viewer) sent_command(viewer$queue[[length(viewer$queue)]]$frame)

test_that("the three repositories carry the same pinned action catalog", {
  expect_identical(cellucid_file_sha256(catalog_path()), .action_catalog_sha256)
  catalog <- read_catalog()
  expect_identical(catalog$schema, 1L)
  expect_gt(length(catalog$actions), 100L)
  for (copy in list(c("cellucid", "tests", "fixtures", "action-catalog.json"), c("cellucid-python", "src", "cellucid", "action_catalog.json"))) {
    path <- do.call(cellucid_sibling_file, as.list(copy))
    expect_identical(cellucid_file_sha256(path), .action_catalog_sha256, info = path)
  }
  # The installed package reads the same file.
  expect_identical(names(cellucid:::.action_catalog()), vapply(catalog$actions, function(action) action$type, ""))
})

test_that("the generated functions, help pages and exports are the catalog's, fresh", {
  root <- cellucid_repository_root()
  catalog <- read_catalog()
  names <- vapply(catalog$actions, function(action) paste0("cellucid_", action$tool), "")
  namespace <- getNamespaceExports("cellucid")
  # One exported function per action, its first argument the viewer.
  for (index in seq_along(names)) {
    expect_true(names[[index]] %in% namespace, info = names[[index]])
    expect_identical(names(formals(get(names[[index]], envir = asNamespace("cellucid"))))[[1L]], "viewer", info = names[[index]])
  }
  if (!cellucid_is_source_checkout()) {
    # The generator is not in the built package (.Rbuildignore); the checkout
    # held the files below to it.
    expect_false(file.exists(file.path(root, "tools", "generate-actions.R")))
    return(invisible(NULL))
  }
  generator <- new.env()
  sys.source(file.path(root, "tools", "generate-actions.R"), envir = generator)
  files <- generator$render_actions(root)
  expect_setequal(names(files), c("R/actions-generated.R", "NAMESPACE", "_pkgdown.yml", file.path("man", paste0(names, ".Rd"))))
  stale <- character()
  for (path in names(files)) {
    on_disk <- readLines(file.path(root, path), warn = FALSE, encoding = "UTF-8")
    if (!identical(on_disk, enc2utf8(files[[path]]))) stale <- c(stale, path)
  }
  expect_identical(stale, character(), info = "run `Rscript tools/generate-actions.R` from the package root")
  # No page is left over from an action the catalog no longer has.
  expect_identical(generator$stale_pages(root, files), character())
})

test_that("a generated name never takes a hand-written export's, and none takes an R reserved word", {
  root <- cellucid_repository_root()
  skip_if_not(cellucid_is_source_checkout(), "the generator is in the checkout only")
  generator <- new.env()
  sys.source(file.path(root, "tools", "generate-actions.R"), envir = generator)
  scratch <- tempfile("cellucid_generate_")
  dir.create(file.path(scratch, "inst"), recursive = TRUE)
  file.copy(catalog_path(), file.path(scratch, "inst", "action_catalog.json"))
  file.copy(file.path(root, "_pkgdown.yml"), file.path(scratch, "_pkgdown.yml"))
  namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  writeLines(c("export(cellucid_color_by)", namespace), file.path(scratch, "NAMESPACE"))
  expect_error(generator$render_actions(scratch), "takes a hand-written export's name: cellucid_color_by")
  writeLines(namespace, file.path(scratch, "NAMESPACE"))
  expect_silent(generator$render_actions(scratch))
  # A parameter named as an R reserved word is the argument <name>_.
  expect_identical(generator$argument_name("in"), "in_")
  expect_identical(names(formals(cellucid_analysis_run))[1:3], c("viewer", "analysis", "in_"))
  expect_identical(cellucid:::.action_param_name("in_"), "in")
  expect_identical(cellucid:::.action_param_name("field_"), "field_")
})

test_that("a generated function sends its action: required parameters, others only when given", {
  viewer <- queued_viewer()
  request <- cellucid_color_by(viewer, field = list(source = "obs", key = "cell_type"))
  expect_s3_class(request, "cellucid_action_request")
  command <- last_sent(viewer)
  expect_identical(command[c("type", "action")], list(type = "action", action = "color.by"))
  expect_identical(command$requestId, request$id)
  expect_identical(command$params, list(field = list(source = "obs", key = "cell_type")))
  expect_identical(request$state, "waiting for the viewer")
  expect_false(request$done)
  expect_null(request$step)
  expect_output(print(request), "<cellucid action color.by \\(a[0-9a-f]{16}\\): waiting for the viewer>")
  # Optional parameters are left out unless given; a default is the viewer's own.
  cellucid_history_get(viewer)
  expect_identical(last_sent(viewer)$params, structure(list(), names = character()))
  cellucid_history_get(viewer, all_branches = TRUE)
  expect_identical(last_sent(viewer)$params, list(all_branches = TRUE))
  # An array parameter is an array whatever its length; NULL is JSON null.
  cellucid_gene_set_save(viewer, "T cells", "CD3E")
  expect_identical(last_sent(viewer)$params, list(name = "T cells", genes = list("CD3E")))
  cellucid_gene_set_save(viewer, "T cells", c("CD3E", "CD3D"), description = "from R")
  expect_identical(last_sent(viewer)$params, list(name = "T cells", genes = list("CD3E", "CD3D"), description = "from R"))
  visibility <- cellucid_cells_visibility_set(viewer, NULL, TRUE)
  expect_identical(viewer$queue[[length(viewer$queue)]]$frame, paste0(
    "{\"type\":\"action\",\"requestId\":\"", visibility$id,
    "\",\"action\":\"cells.visibility.set\",\"params\":{\"cells\":null,\"visible\":true},\"turn\":null,\"turnLabel\":null}"
  ))
  # A reserved word's argument is sent under the parameter's own name.
  cellucid_analysis_run(viewer, "markers", in_ = 2)
  expect_identical(last_sent(viewer)$params, list(analysis = "markers", `in` = 2L))
  # A required parameter not given says which.
  expect_error(cellucid_color_by(viewer), "cellucid_color_by\\(\\) needs field")
  expect_error(cellucid_category_rename(viewer, list(source = "obs", key = "x"), "a"), "needs name")
})

test_that("cells are R's: 1-based positions, TRUE/FALSE or names, sent as 0-based original positions in order", {
  viewer <- queued_viewer()
  cellucid_highlight_cells(viewer, c(3, 1), label = "Two")
  command <- last_sent(viewer)
  expect_identical(command$params, list(label = "Two", cells = c(2, 0)))
  expect_identical(vapply(command$arrays, function(array) array$type, ""), "uint32")
  cellucid_highlight_cells(viewer, c("c6", "c2"))
  expect_identical(last_sent(viewer)$params$cells, c(5, 1))
  cellucid_highlight_cells(viewer, c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
  expect_identical(last_sent(viewer)$params$cells, c(1, 2, 5))
  cellucid_highlight_cells(viewer, factor("c4"))
  expect_identical(last_sent(viewer)$params$cells, 3)
  # In a oneOf (cells, or NULL for every cell) the cells take the x-cells form.
  cellucid_cells_visibility_set(viewer, 2:3, visible = FALSE)
  expect_identical(last_sent(viewer)$params, list(visible = FALSE, cells = c(1, 2)))
  count <- length(viewer$queue)
  expect_error(cellucid_highlight_cells(viewer, c(2, 2)), "names the cell 2 twice")
  expect_error(cellucid_highlight_cells(viewer, c("c1", "c1")), "names the cell 1 twice")
  expect_error(cellucid_highlight_cells(viewer, 7), "1-based positions from 1 to 6")
  expect_error(cellucid_highlight_cells(viewer, 0), "1-based positions from 1 to 6")
  expect_error(cellucid_highlight_cells(viewer, 1.5), "1-based positions")
  expect_error(cellucid_highlight_cells(viewer, NA_real_), "1-based positions")
  expect_error(cellucid_highlight_cells(viewer, c(TRUE, FALSE)), "one value per cell \\(6\\)")
  expect_error(cellucid_highlight_cells(viewer, c(TRUE, NA, FALSE, FALSE, FALSE, FALSE)), "holds NA at cell 2")
  expect_error(cellucid_highlight_cells(viewer, "nope"), "names cells the dataset does not have: c\\(\"nope\"\\)")
  expect_error(cellucid_highlight_cells(viewer, list(1, "a")), "1-based positions")
  # Nothing refused was sent.
  expect_length(viewer$queue, count)
  # Names need the cell names of an object; positions do for a file.
  file_viewer <- queued_viewer(n_cells = NULL, cell_names = NULL)
  expect_error(cellucid_highlight_cells(file_viewer, "c1"), "by name needs a viewer of an R object")
  cellucid_highlight_cells(file_viewer, 4294967296)
  expect_identical(last_sent(file_viewer)$params$cells, 4294967295)
  # Repeated names in the object cannot name one cell.
  twice <- queued_viewer(n_cells = 3L, cell_names = c("a", "a", "b"))
  expect_error(cellucid_highlight_cells(twice, "a"), "holds more than once")
  cellucid_highlight_cells(twice, "b")
  expect_identical(last_sent(twice)$params$cells, 2)
})

test_that("a column takes its value type: codes of a factor, 0/1, whole numbers, numbers", {
  viewer <- queued_viewer()
  push <- function(values, kind = "continuous", categories = NULL) {
    cellucid_field_push(viewer, key = "k", kind = kind, categories = categories, n_obs = 6L, cells = NULL, values = values)
    command <- last_sent(viewer)
    list(values = command$params$values, type = command$arrays[[1L]]$type, params = command$params)
  }
  sent <- push(factor(c("a", NA, "b", "a", "b", "b")), "categorical", c("a", "b"))
  expect_identical(sent$type, "uint8")
  expect_identical(sent$values, c(0, 255, 1, 0, 1, 1))
  expect_identical(sent$params$categories, list("a", "b"))
  expect_null(sent$params$cells)
  sent <- push(factor(c(paste0("x", 1:300), NA))[1:6], "categorical", paste0("x", 1:300))
  expect_identical(sent$type, "uint16")
  sent <- push(c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_identical(sent[c("values", "type")], list(values = c(1, 0, 1, 1, 0, 0), type = "uint8"))
  sent <- push(c(0L, 3L, 255L, 0L, 1L, 2L))
  expect_identical(sent$type, "uint8")
  sent <- push(c(0L, 3L, 256L, 0L, 1L, 2L))
  expect_identical(sent$type, "uint16")
  sent <- push(c(0L, 3L, 70000L, 0L, 1L, 2L))
  expect_identical(sent$type, "uint32")
  sent <- push(c(-1L, 3L, 7L, 0L, 1L, 2L))
  expect_identical(sent$type, "float64")
  sent <- push(c(0.5, NA, -2, 1e10, 0, 1))
  expect_identical(sent$type, "float64")
  expect_true(is.nan(sent$values[[2L]]))
  expect_identical(sent$values[-2L], c(0.5, -2, 1e10, 0, 1))
  expect_error(push(c(1L, NA, 2L, 3L, 4L, 5L)), "holds NA: give a factor")
  expect_error(push(c(TRUE, NA, TRUE, TRUE, FALSE, FALSE)), "holds NA: give a factor")
  expect_error(push(c(1, Inf, 2, 3, 4, 5)), "holds an infinity")
  expect_error(push(c("a", "b")), "takes numbers, TRUE/FALSE, or a factor")
  # Cells of a column: x-cells too.
  cellucid_field_push(viewer, key = "k", kind = "continuous", categories = NULL, n_obs = 6L, cells = c("c2", "c5"),
                      values = c(1.5, 2.5))
  expect_identical(last_sent(viewer)$params$cells, c(1, 4))
})

test_that("cellucid_act() sends any action by its type, a newer viewer's too", {
  viewer <- queued_viewer()
  request <- cellucid_act(viewer, "highlight.cells", cells = 2:3, page = "Notebook")
  expect_identical(last_sent(viewer)$params, list(page = "Notebook", cells = c(1, 2)))
  expect_identical(request$action, "highlight.cells")
  # A parameter named like the start of cellucid_act()'s own arguments is a
  # parameter (`view` is not partly `viewer`).
  cellucid_act(viewer, "view.focus", view = 2)
  expect_identical(last_sent(viewer)$params, list(view = 2L))
  cellucid_act(viewer, "analysis.run", analysis = "markers", in_ = 1)
  expect_identical(last_sent(viewer)$params, list(analysis = "markers", `in` = 1L))
  # An action the catalog does not describe goes as JSON of its shape.
  cellucid_act(viewer, "future.thing", a = c(1, 2), b = list(c = "d"), e = NULL, f = I("one"))
  expect_identical(last_sent(viewer)$params, list(a = list(1L, 2L), b = list(c = "d"), e = NULL, f = list("one")))
  expect_error(cellucid_act(viewer, "colorby"), "dotted type")
  expect_error(cellucid_act(viewer, "color.by", list(source = "obs")), "Name every parameter")
  expect_error(cellucid_act(viewer, "x.y", a = 1, a = 2), "given twice")
  expect_error(cellucid_act(viewer, "x.y", a = NA), "a is NA")
  expect_error(cellucid_act(viewer, "x.y", a = list(b = Inf)), "a\\$b must be a finite number")
  expect_error(cellucid_act(viewer, "x.y", a = mean), "must be text, numbers")
  expect_error(cellucid_act("not a viewer", "x.y"), "Pass the viewer")
  viewer$running <- FALSE
  expect_error(cellucid_act(viewer, "camera.reset", view = "live"), "This viewer is stopped")
})

test_that("an answer reaches its request; a refusal stops its collection, and is said when nobody waits", {
  viewer <- queued_viewer()
  request <- cellucid_history_get(viewer)
  deliver <- function(type, payload) cellucid:::.deliver_event(viewer, list(type = type, payload = payload))
  outcome <- list(head = "s4", steps = list(list(id = "s0", label = "Opened")))
  deliver("action_result", list(requestId = request$id, step = NULL, unrecorded = NULL, outcome = outcome))
  expect_true(request$done)
  expect_identical(request$state, "answered")
  expect_identical(cellucid_action_result(request, timeout = 0), outcome)
  expect_identical(viewer$last$action_result$outcome, outcome)
  change <- cellucid_camera_reset(viewer, "live")
  deliver("action_result", list(requestId = change$id, step = list(id = "s5", label = "Reset the camera"), unrecorded = NULL, outcome = NULL))
  # A step, or a read, warns nothing.
  expect_no_warning(expect_null(cellucid_action_result(change)))
  expect_null(change$unrecorded)
  expect_null(request$unrecorded)
  expect_identical(change$step, list(id = "s5", label = "Reset the camera"))
  expect_output(print(change), "answered, step s5: Reset the camera")
  refused <- cellucid_category_rename(viewer, list(source = "obs", key = "cell_type"), "T", "T cells")
  expect_warning(
    deliver("command_error", list(command = "action", requestId = refused$id, reason = "There is no category \"T\".")),
    "The viewer did not take cellucid_category_rename\\(\\): There is no category \"T\"."
  )
  expect_identical(refused$state, "refused")
  expect_error(cellucid_action_result(refused), "The viewer did not take cellucid_category_rename\\(\\): There is no category")
  expect_output(print(refused), "refused: There is no category")
  # A handler of command_error takes the refusal instead of the warning.
  caught <- NULL
  cellucid_on(viewer, "command_error", function(event) caught <<- event)
  other <- cellucid_act(viewer, "x.y")
  expect_no_warning(deliver("command_error", list(command = "action", requestId = other$id, reason = "Unknown action x.y.")))
  expect_identical(caught$requestId, other$id)
  # Stopping the viewer answers every open request: none will come.
  open <- cellucid_camera_reset(viewer, "live")
  cellucid:::.stop_action_requests(viewer)
  expect_error(cellucid_action_result(open), "stopped before it answered cellucid_camera_reset\\(\\)")
  # In a notebook a cell cannot wait: it says what to do instead.
  notebook <- queued_viewer()
  notebook$widget <- new.env()
  pending <- cellucid_camera_reset(notebook, "live")
  expect_error(cellucid_action_result(pending), "in a Jupyter notebook it answers between cells. Run cellucid_action_result\\(\\) again in a later cell")
  expect_error(cellucid_action_result(list()), "request must be what cellucid_act\\(\\)")
  expect_error(cellucid_action_result(pending, timeout = -1), "timeout must be")
})

test_that("a change that made no step says why: kept on its request, warned when collected", {
  viewer <- queued_viewer()
  deliver <- function(type, payload) cellucid:::.deliver_event(viewer, list(type = type, payload = payload))
  sentences <- cellucid_read_fixture("viewer_protocol.json")$action$unrecorded
  expect_identical(as.list(cellucid:::.viewer_action_unrecorded), sentences)
  for (reason in names(sentences)) {
    change <- cellucid_color_by(viewer, list(source = "obs", key = "cell_type"))
    expect_null(change$unrecorded)
    deliver("action_result", list(requestId = change$id, step = NULL, unrecorded = reason, outcome = list(field = "cell_type")))
    expect_null(change$step)
    expect_identical(change$unrecorded, reason)
    expect_output(print(change), paste0("answered, no step: ", reason), fixed = TRUE)
    caught <- NULL
    outcome <- withCallingHandlers(
      cellucid_action_result(change, timeout = 0),
      cellucid_action_not_recorded = function(condition) {
        caught <<- condition
        invokeRestart("muffleWarning")
      }
    )
    expect_identical(outcome, list(field = "cell_type"))
    expect_s3_class(caught, c("cellucid_action_not_recorded", "warning", "condition"))
    expect_identical(conditionMessage(caught), paste0("cellucid_color_by(): ", sentences[[reason]]))
    expect_identical(caught$unrecorded, reason)
    expect_null(conditionCall(caught))
  }
})

test_that("the catalog is the one the web app's action definitions make now", {
  # Against the definitions themselves, not only the web app's committed copy.
  expect_identical(cellucid_rendered_by_web("scripts/generate-action-catalog.mjs", "renderActionCatalog"),
                   readChar(catalog_path(), file.size(catalog_path()), useBytes = TRUE),
                   info = "Run node scripts/generate-action-catalog.mjs in the web app, then Rscript tools/generate-actions.R here.")
})
