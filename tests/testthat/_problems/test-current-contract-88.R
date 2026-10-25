# Extracted from test-current-contract.R:88

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "cellucid", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
repository_root <- cellucid_repository_root()
roots <- c("R", "man", "vignettes")
scalars <- c("README.md", "DESCRIPTION", "NAMESPACE")
if (cellucid_is_source_checkout()) {
    roots <- c(roots, ".github")
    scalars <- c(scalars, "_pkgdown.yml", "publishing.md")
  }
files <- unlist(
    lapply(
      file.path(repository_root, roots),
      list.files,
      recursive = TRUE,
      full.names = TRUE,
      all.files = TRUE
    ),
    use.names = FALSE
  )
files <- c(files, file.path(repository_root, scalars))
text_suffix <- "\\.(R|Rd|Rmd|md|ya?ml)$"
files <- files[
    !grepl("\\.", basename(files), fixed = TRUE) |
      grepl(text_suffix, files, ignore.case = TRUE)
  ]
expect_true(all(
    file.path(repository_root, scalars) %in% files
  ))
expect_true(all(
    file.path(
      repository_root,
      c("R/cellucid_prepare.R", "man/cellucid_prepare.Rd")
    ) %in% files
  ))
retired_language <- paste0(
    "\\bfallback\\b|fall(?:ing)?[[:space:]-]+back|",
    "best[[:space:]_-]*effort|\\blegacy\\b|\\bdeprecated\\b|",
    "deprecation|backwards?[[:space:]_-]*compat|\\bshim\\b"
  )
violations <- character()
for (file in files) {
    lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
    matches <- grep(retired_language, lines, ignore.case = TRUE, perl = TRUE)
    if (length(matches) > 0L) {
      relative <- substring(
        normalizePath(file, mustWork = TRUE),
        nchar(repository_root) + 2L
      )
      violations <- c(
        violations,
        sprintf("%s:%d: %s", relative, matches, trimws(lines[matches]))
      )
    }
  }
expect_identical(violations, character())
