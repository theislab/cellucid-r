# Every size or count limit in the package is a real one, and says why.
#
# A package-level constant whose name has max, limit, cap, capacity, budget or
# ceiling in it (R/), and every `MAX`/`LIMIT` bound in the C sources (src/),
# has to appear below with its class and the reason it exists. A new one fails
# this test until someone classifies it. The classes are FORMAT (fixed by a
# file format), PLAT (a real platform limit), HOST (a bound on untrusted input),
# PERF (decides what is kept, never what is refused) and MSG (how much a
# message quotes); a limit made up for convenience has no class to take.

size_limit_allowlist <- list(
  "R/matrix-input.R:.float32_max_magnitude" = c(
    "FORMAT",
    "The largest finite IEEE 754 binary32 value; a larger one is not a float32."
  ),
  "R/matrix-input.R:.float32_exact_integer_limit" = c(
    "FORMAT",
    "The largest integer float32 holds exactly (2^24); counts above it are stored as float64."
  ),
  "R/zip.R:.zip_u16_max" = c(
    "FORMAT",
    "The largest classic ZIP entry count; above it the archive carries ZIP64 records."
  ),
  "R/zip.R:.zip_u32_max" = c(
    "FORMAT",
    "The largest classic ZIP size or offset; above it the archive carries ZIP64 fields."
  ),
  "R/error-messages.R:.non_finite_example_limit" = c(
    "MSG",
    "How many offending values an error message quotes; every value is checked."
  )
)

size_limit_classes <- c("FORMAT", "PLAT", "HOST", "PERF", "MSG")

package_root <- function() {
  candidates <- c(
    testthat::test_path("..", ".."),
    file.path(testthat::test_path("..", ".."), "00_pkg_src", "cellucid")
  )
  for (candidate in candidates) {
    if (dir.exists(file.path(candidate, "R"))) {
      return(normalizePath(candidate))
    }
  }
  NULL
}

scan_size_limits <- function(root) {
  found <- character()
  name_pattern <- "(^|[._])(max|maximum|limit|limits|cap|capacity|budget|ceiling)([._]|$)"
  for (path in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) {
    lines <- readLines(path, warn = FALSE)
    declared <- regmatches(lines, regexpr("^[.A-Za-z][.A-Za-z0-9_]*(?=\\s*<-)", lines, perl = TRUE))
    for (name in declared) {
      if (grepl(name_pattern, tolower(name))) {
        found <- c(found, paste0("R/", basename(path), ":", name))
      }
    }
  }
  for (path in list.files(file.path(root, "src"), pattern = "[.](c|h)$", full.names = TRUE)) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    names <- regmatches(text, gregexpr("\\b[A-Z][A-Z0-9_]*(MAX|LIMIT)[A-Z0-9_]*\\b", text))[[1]]
    for (name in unique(names)) {
      found <- c(found, paste0("src/", basename(path), ":", name))
    }
  }
  sort(unique(found))
}

test_that("every size or count limit is classified with a reason", {
  root <- package_root()
  skip_if(is.null(root), "package sources are not available")
  found <- scan_size_limits(root)
  expect_identical(setdiff(found, names(size_limit_allowlist)), character())
  if (dir.exists(file.path(root, "src"))) {
    expect_identical(setdiff(names(size_limit_allowlist), found), character())
  }
  for (key in names(size_limit_allowlist)) {
    entry <- size_limit_allowlist[[key]]
    expect_true(entry[[1]] %in% size_limit_classes, info = key)
    expect_true(nzchar(trimws(entry[[2]])), info = key)
  }
})
