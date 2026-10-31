# A notebook cell run again replaces its output: the viewer cellucid_view()
# showed in the cell's earlier run, and no other cell shows, is stopped (its
# object's matrices and chunk cache released), as cellucid-python stops the
# widget a re-run cell replaced (round-2 audit MB-9). Where no cell id is
# known (a terminal, a console, a front end that sends none), showing a
# dataset again stops the earlier viewer of the same dataset id, and says so.

fake_viewer <- function(dataset_id) {
  viewer <- new.env(parent = emptyenv())
  class(viewer) <- "cellucid_viewer"
  viewer$viewer_id <- .random_hex(8L)
  viewer$dataset_id <- dataset_id
  viewer$running <- TRUE
  viewer$store <- list(matrices = "held")
  assign(viewer$viewer_id, viewer, envir = .viewer_registry)
  viewer
}

test_that("a cell run again stops the viewer it showed before, unless another cell shows it", {
  on.exit({
    cellucid_stop_all()
    rm(list = ls(.shown_by_cell), envir = .shown_by_cell)
  }, add = TRUE)
  run <- list(cell = "cell-a", count = 1L)
  local_mocked_bindings(.notebook_run = function() run, .package = "cellucid")
  first <- fake_viewer("demo")
  .release_replaced_viewers(first, quiet = TRUE)
  second <- fake_viewer("demo")
  .release_replaced_viewers(second, quiet = TRUE)
  expect_true(first$running && second$running)  # one run's viewers are kept
  run <- list(cell = "cell-a", count = 2L)
  third <- fake_viewer("demo")
  .release_replaced_viewers(third, quiet = TRUE)
  expect_false(first$running)
  expect_false(second$running)
  expect_null(first$store)
  expect_true(third$running)
  # Shown by another cell too: kept when its first cell runs again.
  run <- list(cell = "cell-b", count = 3L)
  .release_replaced_viewers(third, quiet = TRUE)
  run <- list(cell = "cell-a", count = 4L)
  fourth <- fake_viewer("demo")
  .release_replaced_viewers(fourth, quiet = TRUE)
  expect_true(third$running)
})

test_that("where no cell id is known, showing a dataset again stops its earlier viewer and says so", {
  on.exit(cellucid_stop_all(), add = TRUE)
  local_mocked_bindings(.notebook_run = function() NULL, .package = "cellucid")
  first <- fake_viewer("demo")
  other <- fake_viewer("other")
  again <- fake_viewer("demo")
  expect_message(.release_replaced_viewers(again, quiet = FALSE), "Stopped the earlier viewer of 'demo'")
  expect_false(first$running)
  expect_true(other$running)
  expect_true(again$running)
})

test_that("the running notebook cell is read from IRkernel's executor on the call stack", {
  local_mocked_bindings(.in_jupyter_kernel = function() TRUE, .package = "cellucid")
  Executor <- setRefClass("Executor", fields = list(current_request = "list", execution_count = "integer"),
    methods = list(execute = function(code) eval(code)))
  executor <- Executor$new(current_request = list(metadata = list(cellId = "cell-7")), execution_count = 5L)
  expect_equal(executor$execute(quote(cellucid:::.notebook_run())), list(cell = "cell-7", count = 5L))
  silent <- Executor$new(current_request = list(metadata = list()), execution_count = 6L)
  expect_null(silent$execute(quote(cellucid:::.notebook_run())))
  expect_null(.notebook_run())
})

test_that("a served store keeps a 64th of the machine's memory by default, as cellucid-python's", {
  machine <- .native_physical_memory()
  expect_true(is.numeric(machine) && length(machine) == 1L)
  expect_true(is.na(machine) || machine > 2^28)
  expected <- if (is.na(machine)) 8 * 2^30 / 64 else floor(machine / 64)
  expect_equal(.default_cache_bytes(), expected)
})
