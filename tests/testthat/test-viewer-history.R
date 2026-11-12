# The viewer's history, live, over the channel (R/viewer-api.R "The history,
# live"): notebook turns carried by every action sent inside one, a session
# file pushed into the viewer and answered by session_loaded, every
# history_changed kept and handed to the hooks, and the history's actions as
# the generated cellucid_history_*() functions -- cellucid-python's
# viewer.turn(), viewer.load_session() and viewer.history.

history_session_file <- function() testthat::test_path("fixtures", "history_golden", "history_golden.cellucid-session")

test_that("actions inside a turn carry it, and turns do not nest", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  client <- ws_welcomed(viewer)
  inside <- list()
  turn <- cellucid_turn(viewer, "  Normalize   and cluster ", {
    inside$color <- cellucid_color_by(viewer, field = list(source = "obs", key = "cluster"))
    inside$highlight <- cellucid_highlight_cells(viewer, c("cell1", "cell2"), color = "#00ff00")
    expect_error(cellucid_turn(viewer, "Inner", NULL), 'A turn \\("Normalize and cluster"\\) is open already; turns do not nest.')
  })
  expect_match(turn, "^r-[0-9a-f]{12}$")
  after <- cellucid_history_undo(viewer)
  first <- ws_receive(client)$json
  expect_identical(list(first$requestId, first$turn, first$turnLabel), list(inside$color$id, turn, "Normalize and cluster"))
  second <- ws_receive(client)$json
  expect_identical(list(second$requestId, second$turn, second$turnLabel), list(inside$highlight$id, turn, "Normalize and cluster"))
  third <- ws_receive(client)$json
  expect_identical(list(third$action, third$turn, third$turnLabel, third$requestId), list("history.undo", NULL, NULL, after$id))
  # A turn that stops on an error closes all the same.
  expect_error(cellucid_turn(viewer, "Broken", stop("no")), "no")
  expect_null(viewer$turn)
  expect_error(cellucid_turn(viewer, "   ", NULL), "A turn needs a label")
  expect_error(cellucid_turn(viewer, "QC"), "runs code inside the turn")
  ws_close(client)
})

test_that("a session pushed in is one binary frame, answered by session_loaded", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  client <- ws_welcomed(viewer)
  path <- history_session_file()
  request <- cellucid_load_session(viewer, path)
  frame <- ws_receive(client)
  expect_identical(frame$opcode, 2L)
  command <- frame$json
  expect_identical(command[c("type", "requestId", "name")],
                   list(type = "load_session", requestId = request$id, name = basename(path)))
  expect_identical(command$file, readBin(path, "raw", file.size(path)))
  expect_output(print(request), "waiting")
  ws_send_json(client, list(type = "session_loaded", requestId = request$id,
                            tracker = list(id = "t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10", name = "Annotation"), head = "s13"))
  expect_identical(cellucid_load_session_result(request, timeout = 5),
                   list(tracker = list(id = "t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10", name = "Annotation"), head = "s13"))
  expect_error(cellucid_load_session_result(request), "collected already")
  refused <- cellucid_load_session(viewer, path)
  ws_receive(client)
  ws_send_json(client, list(type = "command_error", command = "load_session", requestId = refused$id,
                            reason = "This session was saved on another dataset."))
  expect_error(cellucid_load_session_result(refused, timeout = 5),
               "The viewer did not load history_golden.cellucid-session: This session was saved on another dataset.", fixed = TRUE)
  expect_error(cellucid_load_session(viewer, file.path(tempdir(), "absent.cellucid-session")), "No session file at")
  ws_close(client)
})

test_that("every history change is kept and reaches the hooks", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  client <- ws_welcomed(viewer)
  seen <- list()
  cellucid_on(viewer, "history_changed", function(change) seen[[length(seen) + 1L]] <<- change)
  expect_null(cellucid_last_history_change(viewer))
  change <- list(type = "history_changed", kind = "step", tracker = "t-1", head = "s4", step = "s4", branch = "b0")
  ws_send_json(client, change)
  moved <- utils::modifyList(change, list(kind = "moved", head = "s2"))
  moved["step"] <- list(NULL)
  moved["branch"] <- list(NULL)
  ws_send_json(client, moved)
  run_loop_until(function() length(seen) == 2L)
  expect_identical(seen[[1L]], change[-1L])
  expect_identical(seen[[2L]], moved[-1L])
  expect_identical(cellucid_last_history_change(viewer), moved[-1L])
  # A kind the protocol does not have breaks it.
  ws_send_json(client, utils::modifyList(change, list(kind = "renamed")))
  expect_identical(ws_receive(client)$code, 1008L)
  ws_close(client)
})

test_that("the history's actions are the generated functions", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  client <- ws_welcomed(viewer)
  requests <- list(
    list(cellucid_history_jump(viewer, step = "s3"), "history.jump", list(step = "s3")),
    list(cellucid_history_note_add(viewer, target = list(step = "s3"), text = "Looks right."), "history.note.add",
         list(target = list(step = "s3"), text = "Looks right.")),
    list(cellucid_history_turn_undo(viewer, turn = "r-1", only_own = TRUE), "history.turn.undo", list(turn = "r-1", only_own = TRUE))
  )
  for (entry in requests) {
    command <- ws_receive(client)$json
    expect_identical(list(command$requestId, command$action, command$params), list(entry[[1L]]$id, entry[[2L]], entry[[3L]]))
  }
  ws_close(client)
})
