# What opened a viewer, as each entry point names it. The viewer names its
# host in every sentence it says from the runtime of _cellucid/info (the
# shared protocol's runtimes): named as an R session, an agent app's MCP
# server that ends would freeze the page with "The notebook or R session that
# opened this viewer stopped". An R session (cellucid_view(),
# cellucid_serve()) names itself "r", the agent app's cellucid_mcp(open = )
# "r-mcp"; the agents record carries the same.

published_runtime <- function(viewer) {
  records <- Filter(function(record) identical(record$viewer, viewer$viewer_id), cellucid:::.agent_discover())
  expect_length(records, 1L)
  records[[1L]]$runtime
}

test_that("the protocol declares both runtimes this package names, as R's", {
  runtimes <- cellucid_read_fixture("viewer_protocol.json")$runtimes
  expect_setequal(names(runtimes), c("jupyter", "python", "cellucid-serve", "cellucid-mcp", "r", "r-mcp"))
  for (runtime in c("r", "r-mcp")) expect_identical(runtimes[[runtime]]$language, "r")
})

test_that("cellucid_serve() and cellucid_view() name the R session", {
  skip_without_sce()
  served <- test_viewer(small_sce(alt = FALSE), dataset_id = "served", dataset_name = "Served")
  expect_identical(http_json(paste0(served$url, "/_cellucid/info"))$runtime, "r")
  expect_identical(published_runtime(served), "r")
  expect_identical(cellucid:::.agent_port_describe(served)$runtime, "r")

  viewed <- cellucid_view(
    small_sce(alt = FALSE), dataset_id = "viewed", pane = "none", port = 0L,
    web_cache_dir = viewer_fixture_web_cache(), web_source_url = NULL, quiet = TRUE
  )
  on.exit(cellucid_stop(viewed), add = TRUE)
  expect_identical(http_json(paste0(viewed$url, "/_cellucid/info"))$runtime, "r")
  expect_identical(published_runtime(viewed), "r")
})

test_that("cellucid_mcp(open = ) names the agent app's server", {
  folder <- suppressMessages(do.call(cellucid_prepare, prepare_fixture(file.path(tempfile(), "f.zarr"))))
  session <- cellucid:::.mcp_session()
  # Its address goes to standard error as a message (never shown here: it
  # carries the viewer's token).
  viewer <- suppressMessages(cellucid:::.mcp_open_viewer(session, folder, 0L, FALSE, NULL, viewer_fixture_web_cache()))
  on.exit(cellucid_stop(viewer), add = TRUE)
  expect_identical(viewer$runtime, "r-mcp")
  expect_identical(http_json(paste0(viewer$url, "/_cellucid/info"))$runtime, "r-mcp")
  expect_identical(published_runtime(viewer), "r-mcp")
  expect_identical(cellucid:::.agent_port_describe(viewer)$runtime, "r-mcp")
})
