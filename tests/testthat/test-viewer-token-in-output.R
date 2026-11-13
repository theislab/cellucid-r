# The viewer token never reaches what is logged or pasted. The address that
# carries it is said once, to the person, where they are told to open it;
# an error, a debug report, a printed viewer and the test suite's own output
# say the address without it (.public_viewer_url(), as cellucid-python's
# public_viewer_url() and tests/test_viewer_token_in_output_contract.py).

test_that(".public_viewer_url() drops the token wherever it stands, as Python's", {
  cases <- c(
    "http://127.0.0.1:8765/?viewer=true&viewerId=ab&viewerToken=cd&anndata=true" =
      "http://127.0.0.1:8765/?viewer=true&viewerId=ab&anndata=true",
    "a?viewerToken=X&b=1" = "a?b=1",
    "a?b=1&viewerToken=X" = "a?b=1",
    "a?viewerToken=X" = "a",
    "a?b=1" = "a?b=1",
    "x?viewerToken=1#frag" = "x#frag",
    "GET /?viewer=true&viewerToken=ff&anndata=true HTTP/1.1" = "GET /?viewer=true&anndata=true HTTP/1.1",
    "open a?viewerToken=1 or b?c=2&viewerToken=3" = "open a or b?c=2"
  )
  expect_identical(unname(cellucid:::.public_viewer_url(names(cases))), unname(cases))
})

test_that("a printed viewer and its errors say the address without the token", {
  skip_without_sce()
  viewer <- test_viewer(small_sce(alt = FALSE), dataset_id = "s", dataset_name = "S")
  printed <- paste(capture.output(print(viewer)), collapse = "\n")
  expect_false(grepl(viewer$viewer_token, printed, fixed = TRUE))
  expect_match(printed, paste0("address:  ", cellucid:::.public_viewer_url(viewer$viewer_url), " (opened by viewer$viewer_url)"), fixed = TRUE)
  expect_match(viewer$viewer_url, paste0("&viewerToken=", viewer$viewer_token, "&"), fixed = TRUE)
})

test_that("cellucid_view(pane = \"none\") says the address to open, once, with its token", {
  skip_without_sce()
  said <- capture_messages(viewer <- cellucid_view(
    small_sce(alt = FALSE), dataset_id = "none", dataset_name = "None", pane = "none", port = 0L,
    web_cache_dir = viewer_fixture_web_cache(), web_source_url = NULL
  ))
  on.exit(cellucid_stop(viewer), add = TRUE)
  expect_identical(grep(viewer$viewer_token, said, fixed = TRUE, value = TRUE), paste0("Open ", viewer$viewer_url, " in a browser.\n"))
})

test_that("cellucid_mcp(open = ) says the address to open once, on standard error, and no other line says its token", {
  folder <- suppressMessages(do.call(cellucid_prepare, prepare_fixture(file.path(tempfile(), "f.zarr"))))
  session <- cellucid:::.mcp_session()
  said <- capture_messages(
    viewer <- cellucid:::.mcp_open_viewer(session, folder, 0L, FALSE, NULL, viewer_fixture_web_cache())
  )
  on.exit(cellucid_stop(viewer), add = TRUE)
  expect_identical(grep(viewer$viewer_token, said, fixed = TRUE, value = TRUE), paste0("Cellucid viewer: ", viewer$viewer_url, "\n"))
})
