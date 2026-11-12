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

test_that("only the sentences that tell the person what to open say the tokened address", {
  # Every line of R/ naming a viewer's address, outside comments: the address
  # opened, stripped, or said to the person to open. A new one is classified
  # here, so an error or a printed object never carries the token unseen.
  said_to_open <- c(
    # cellucid_serve(): the address to open, with how to reach it over SSH.
    "\"Open \", viewer$viewer_url, \" in a browser\",",
    # cellucid_view(pane = "none"), or outside an interactive session
    # (Rscript): no browser opened.
    "if (!quiet) .inform(\"Open \", viewer$viewer_url, \" in a browser.\")",
    # cellucid_mcp(open = ): on stderr, the agent app's console.
    ".inform(\"Cellucid viewer: \", viewer$viewer_url)"
  )
  root <- test_path("..", "..", "R")
  skip_if_not(dir.exists(root), "the package sources are not beside the tests")
  files <- sort(list.files(root, pattern = "[.]R$", full.names = TRUE))
  unclassified <- character()
  for (file in files) {
    lines <- trimws(readLines(file, warn = FALSE))
    for (index in grep("viewer_url", lines, fixed = TRUE)) {
      line <- lines[[index]]
      if (startsWith(line, "#")) next
      # Made, opened, or stripped of its token; or named, not written out.
      if (grepl("^viewer[$]viewer_url <- |browseURL[(]viewer[$]viewer_url[)]|viewer_option[(]viewer[$]viewer_url[)]|[.]public_viewer_url[(]|^[.]public_viewer_url <- |is[.]null[(](x|viewer)[$]viewer_url[)]|^viewer_url = public_url,$|[(]viewer[$]viewer_url[)][.,]", line)) next
      if (line %in% said_to_open) next
      unclassified <- c(unclassified, paste0(basename(file), ":", index, ": ", line))
    }
  }
  expect_identical(unclassified, character())
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
