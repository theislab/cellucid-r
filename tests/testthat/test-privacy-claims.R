# What the package says about privacy matches what Cellucid does: community
# annotation signs in with GitHub on www.cellucid.com (port 8765 is Google
# Drive's local address, not community annotation's); a server on 127.0.0.1 is
# readable by every user of a shared machine; a session applied to an object
# records no `source`.

test_that("no R surface repeats a privacy claim the code contradicts", {
  repository_root <- cellucid_repository_root()
  files <- unlist(lapply(file.path(repository_root, c("R", "man", "vignettes")), list.files,
                         pattern = "\\.(R|Rd|Rmd)$", recursive = TRUE, full.names = TRUE))
  expect_true(file.path(repository_root, "man", "cellucid_view.Rd") %in% files)
  stale <- c(
    "Google sign-in \\(community annotation\\)",
    "\\(community annotation\\) works only on port",
    "only your machine can"
  )
  found <- character()
  for (file in files) {
    lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
    for (pattern in stale) {
      hits <- grep(pattern, lines, perl = TRUE)
      if (length(hits) > 0L) found <- c(found, sprintf("%s:%d", basename(file), hits))
    }
  }
  expect_identical(found, character())
  view <- paste(readLines(file.path(repository_root, "man", "cellucid_view.Rd"), warn = FALSE), collapse = "\n")
  expect_match(view, "is every user of it", fixed = TRUE)
  apply <- paste(readLines(file.path(repository_root, "man", "cellucid_apply_session.Rd"), warn = FALSE), collapse = "\n")
  expect_match(apply, "without\n\\code{source}", fixed = TRUE)
})
