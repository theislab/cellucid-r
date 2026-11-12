# One spelling: US, in every word a person or an agent reads from this package.
#
# The web app held its strings to US spelling ("color", "neighbor", "canceled",
# "catalog") while this package's DESCRIPTION, help, vignettes and code said
# "neighbour", "catalogue", "labelled" and "centre" (round 9, R9N-1): the
# contract lived only in the web repository. The British word list
# (fixtures/spelling.json) is shared byte for byte with the web app
# (tests/failure-wording-contract.test.mjs) and cellucid-python
# (tests/test_spelling_contract.py), at the digest pinned here.
#
# Read here: every string of R/ that is prose (words with a space between, or a
# one-word sentence such as "Canceled."), the help pages (man/, outside \code{}
# and its kin), the vignettes, the README, DESCRIPTION and the unreleased section
# of NEWS.md. A word in backticks is a name, and so is a single word in quotes
# (a colormap id "greys", a parameter a message echoes), as is a word that is
# part of an identifier, a path or a flag.

.spelling_sha256 <- "f195abcd84ce15a4d39575042ac1ad205a0b4e1b6b8b439d361d4c0dfbbfa77c"

spelling_pattern <- function() {
  patterns <- vapply(cellucid_read_fixture("spelling.json")$british, function(entry) entry$pattern, character(1))
  paste0("(?i)(?<![\\w.$/\\\\]|--)(?:", paste(patterns, collapse = "|"), ")(?![\\w$/\\\\]|\\.\\w|\\s*=(?!=)|\\()")
}

british_words <- function(text, pattern = spelling_pattern()) {
  text <- gsub("`+[^`\n]*?`+", " ", text, perl = TRUE)
  text <- gsub("\\\\?([\"'])[\\w-]+\\\\?\\1", " ", text, perl = TRUE)
  found <- regmatches(text, gregexpr(pattern, text, perl = TRUE))
  unlist(found, use.names = FALSE)
}

text_offenders <- function(lines, where, pattern) {
  offenders <- character()
  for (index in seq_along(lines)) {
    words <- british_words(lines[[index]], pattern)
    if (length(words) > 0L) offenders <- c(offenders, paste0(where, ":", index, ": ", paste(words, collapse = ", ")))
  }
  offenders
}

string_offenders <- function(file, pattern) {
  data <- utils::getParseData(parse(file, keep.source = TRUE), includeText = TRUE)
  strings <- data[data$token == "STR_CONST", ]
  offenders <- character()
  for (index in seq_len(nrow(strings))) {
    body <- substr(strings$text[[index]], 2L, nchar(strings$text[[index]]) - 1L)
    if (!grepl("\\s", body) && !grepl("^[A-Za-z]+[.!?]+$", body)) next
    words <- british_words(body, pattern)
    if (length(words) > 0L) {
      offenders <- c(offenders, paste0(basename(dirname(file)), "/", basename(file), ":", strings$line1[[index]], ": ", paste(words, collapse = ", ")))
    }
  }
  offenders
}

test_that("the word list is the one shared with the web app and cellucid-python, at the pinned digest", {
  expect_pinned_fixture("spelling.json", .spelling_sha256)
  expect_identical(cellucid_read_fixture("spelling.json")$version, 1L)
})

test_that("the check finds British words in prose and never in names", {
  pattern <- spelling_pattern()
  pairs <- list(
    c("neighbour", "neighbor"), c("Colours", "Colors"), c("greyed", "grayed"), c("cancelled", "canceled"),
    c("catalogue", "catalog"), c("unrecognised", "unrecognized"), c("reorganisation", "reorganization"),
    c("learnt", "learned"), c("analysed", "analyzed"), c("relabelled", "relabeled"), c("centred", "centered")
  )
  for (pair in pairs) {
    expect_identical(british_words(paste("the", pair[[1]], "one"), pattern), pair[[1]])
    expect_identical(british_words(paste("the", pair[[2]], "one"), pattern), character())
  }
  expect_identical(british_words("x.cancelled notifications/cancelled --colour colour= normalise( `cancelled` \"greys\"", pattern), character())
  expect_identical(british_words("cancellation, analyses, organism, summaries, precise, raised", pattern), character())
})

test_that("every message, help page, vignette and page of the package is US spelling", {
  root <- cellucid_repository_root()
  pattern <- spelling_pattern()
  offenders <- unlist(lapply(sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)), string_offenders, pattern = pattern))
  code <- "\\\\(?:code|verb|link|linkS4class|pkg|file|env|option|var|samp|kbd|href|url|method|S3method|alias|name|keyword|usage|item\\{[^{}]*\\}(?=\\{))\\{(?:[^{}]|\\{[^{}]*\\})*\\}"
  for (page in sort(list.files(file.path(root, "man"), pattern = "\\.Rd$", full.names = TRUE))) {
    lines <- gsub(code, " ", readLines(page, warn = FALSE, encoding = "UTF-8"), perl = TRUE)
    offenders <- c(offenders, text_offenders(lines, paste0("man/", basename(page)), pattern))
  }
  for (vignette in sort(list.files(file.path(root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE))) {
    offenders <- c(offenders, text_offenders(readLines(vignette, warn = FALSE, encoding = "UTF-8"), paste0("vignettes/", basename(vignette)), pattern))
  }
  for (name in c("README.md", "DESCRIPTION")) {
    offenders <- c(offenders, text_offenders(readLines(file.path(root, name), warn = FALSE, encoding = "UTF-8"), name, pattern))
  }
  news <- readLines(file.path(root, "NEWS.md"), warn = FALSE, encoding = "UTF-8")
  released <- grep("^# cellucid [0-9]", news)
  expect_gte(length(released), 2L)
  offenders <- c(offenders, text_offenders(news[seq_len(released[[2L]] - 1L)], "NEWS.md", pattern))
  expect_identical(offenders, character())
})

test_that("the tests pin the words the package says", {
  # A test that expects British text forces the source back to it, or fails the day the source is
  # right. Every prose string of the suite is read: an expectation, and an input the package
  # echoes back ("Normalize and cluster", a turn's label). A misspelled name a test feeds in on
  # purpose (an unknown key "colour") is one word, so no prose. This file's own examples are the
  # word list's.
  pattern <- spelling_pattern()
  files <- setdiff(sort(list.files(test_path(), pattern = "\\.R$", full.names = TRUE)), test_path("test-spelling.R"))
  offenders <- unlist(lapply(files, string_offenders, pattern = pattern))
  expect_identical(offenders, character())
})
