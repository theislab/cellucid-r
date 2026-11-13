# Comments, test descriptions and workflow comments hold no history: no finding, audit or plan id, no
# round or hunt story, no CI run address and no pointer into the private plans.
# The rule is "Comments say what the code does" in CONTRIBUTING.md; each pattern
# is an exact shape, so a hit is a defect.

history_patterns <- c(
  "CI run" = "actions/runs/[0-9]+|\\brun [0-9]{8,}\\b|\\b4[0-4][0-9]{9}\\b",
  "finding id" = paste0(
    "(?<!kbd )(?<!noqa: )(?<!Shift\\+)(?<!Ctrl\\+)(?<!Alt\\+)(?<!press )(?<!key )(?<!`)",
    "(?:\\bF[0-9]{2,3}\\b|\\bL[0-9]{1,2}-[0-9]{1,2}\\b|\\bR[0-9]{1,2}[A-Za-z]{0,2}-[0-9]{1,3}\\b",
    "|\\bRC[A-Z]-[0-9]{1,2}\\b|\\bCEL-(?:AUDIT-)?[0-9]{2,4}\\b|\\bCEL-[A-Z]+-[0-9]+\\b)"
  ),
  "bug hunt" = "(?i)bug[- ]?hunt|\\bround[- ][0-9]+\\b|\\blive[- ]pass\\b|\\blive[- ][0-9]+\\b",
  "plan pointer" = "History v2 ?§|cellucid-private"
)

find_history <- function(text) {
  hits <- character()
  for (class in names(history_patterns)) {
    found <- regmatches(text, gregexpr(history_patterns[[class]], text, perl = TRUE))[[1L]]
    hits <- c(hits, if (length(found)) paste0(class, ": ", found))
  }
  hits
}

# The `#` comments of R source: COMMENT tokens of the parse data, never a `#` in a string.
r_comments <- function(lines, name = "x") {
  data <- utils::getParseData(parse(text = lines, keep.source = TRUE, srcfile = srcfilecopy(name, lines)))
  data <- data[data$token == "COMMENT", ]
  list(line = data$line1, text = data$text)
}

# The comments of C source: a small scan that skips string and character literals.
c_comments <- function(lines) {
  chars <- strsplit(paste(lines, collapse = "\n"), "")[[1L]]
  n <- length(chars)
  line <- 1L
  i <- 1L
  out_line <- integer()
  out_text <- character()
  while (i <= n) {
    ch <- chars[[i]]
    nxt <- if (i < n) chars[[i + 1L]] else ""
    if (ch == "\n") {
      line <- line + 1L
      i <- i + 1L
    } else if (ch == "\"" || ch == "'") {
      i <- i + 1L
      while (i <= n && chars[[i]] != ch && chars[[i]] != "\n") i <- i + if (chars[[i]] == "\\") 2L else 1L
      i <- i + 1L
    } else if (ch == "/" && nxt == "/") {
      start <- i
      while (i <= n && chars[[i]] != "\n") i <- i + 1L
      out_line <- c(out_line, line)
      out_text <- c(out_text, paste(chars[start:(i - 1L)], collapse = ""))
    } else if (ch == "/" && nxt == "*") {
      start <- i
      first <- line
      i <- i + 2L
      while (i < n && !(chars[[i]] == "*" && chars[[i + 1L]] == "/")) {
        if (chars[[i]] == "\n") line <- line + 1L
        i <- i + 1L
      }
      i <- i + 2L
      out_line <- c(out_line, first)
      out_text <- c(out_text, paste(chars[start:min(i - 1L, n)], collapse = ""))
    } else {
      i <- i + 1L
    }
  }
  list(line = out_line, text = out_text)
}

# The `#` comments of a YAML file: a `#` that opens a line or follows white space outside
# a quoted string.
yaml_comments <- function(lines) {
  out_line <- integer()
  out_text <- character()
  for (number in seq_along(lines)) {
    chars <- strsplit(lines[[number]], "")[[1L]]
    quote <- ""
    i <- 1L
    while (i <= length(chars)) {
      ch <- chars[[i]]
      before <- if (i > 1L) chars[[i - 1L]] else ""
      if (nzchar(quote)) {
        if (ch == "\\" && quote == "\"") i <- i + 1L else if (ch == quote) quote <- ""
      } else if ((ch == "\"" || ch == "'") && (i == 1L || before %in% c(" ", "\t", ":", "[", "{", ",", "-"))) {
        quote <- ch
      } else if (ch == "#" && (i == 1L || before %in% c(" ", "\t"))) {
        out_line <- c(out_line, number)
        out_text <- c(out_text, paste(chars[(i + 1L):max(i + 1L, length(chars))], collapse = ""))
        break
      }
      i <- i + 1L
    }
  }
  list(line = out_line, text = out_text)
}

# The first argument of each test_that() call.
test_descriptions <- function(lines, name = "x") {
  data <- utils::getParseData(parse(text = lines, keep.source = TRUE, srcfile = srcfilecopy(name, lines)))
  calls <- data[data$token == "SYMBOL_FUNCTION_CALL" & data$text == "test_that", ]
  line <- integer()
  text <- character()
  for (index in seq_len(nrow(calls))) {
    call <- str2lang(utils::getParseText(data, data$parent[data$id == calls$parent[[index]]]))
    if (length(call) > 1L && is.character(call[[2L]])) {
      line <- c(line, calls$line1[[index]])
      text <- c(text, call[[2L]])
    }
  }
  list(line = line, text = text)
}

source_files <- function(root) {
  c(
    list.files(file.path(root, c("R", "tests", "tools")), "\\.R$", recursive = TRUE, full.names = TRUE),
    list.files(file.path(root, "src"), "\\.[ch]$", full.names = TRUE)
  )
}

workflow_files <- function(root) {
  list.files(file.path(root, ".github", "workflows"), "\\.ya?ml$", full.names = TRUE)
}

history_in_file <- function(file) {
  lines <- readLines(file, encoding = "UTF-8", warn = FALSE)
  parts <- if (grepl("\\.R$", file)) {
    list(r_comments(lines, basename(file)), test_descriptions(lines, basename(file)))
  } else if (grepl("\\.ya?ml$", file)) {
    list(yaml_comments(lines))
  } else {
    list(c_comments(lines))
  }
  found <- character()
  for (part in parts) {
    for (index in seq_along(part$text)) {
      hits <- find_history(part$text[[index]])
      if (length(hits)) found <- c(found, paste0(file, ":", part$line[[index]], ": ", hits))
    }
  }
  found
}

test_that("the checker finds history planted in R, C and test_that text", {
  planted <- c("F123 was a bug", "see L4-12", "R9a-3: order", "CEL-AUDIT-42", "found in round 9",
               "the bug hunt", "the live pass", "https://x/actions/runs/123456", "run 12345678901",
               "cellucid-private/plans/open-items.md", "History v2 §3.8")
  for (text in planted) expect_gt(length(find_history(text)), 0L, label = text)
  expect_length(find_history("the key F5 and a 4-byte id"), 0L)

  r <- c("x <- \"F123 # in a string\"", "y <- 1  # F123 in a comment")
  comments <- r_comments(r)
  expect_equal(comments$line, 2L)
  expect_length(find_history(comments$text), 1L)

  c_source <- c("char *s = \"/* F123 */\";", "int a; // round 9", "/* line one", "   F123 */", "int b;")
  comments <- c_comments(c_source)
  expect_equal(comments$line, c(2L, 3L))
  expect_length(find_history(comments$text[[1L]]), 1L)
  expect_length(find_history(comments$text[[2L]]), 1L)

  workflow <- yaml_comments(c("name: a", "# run 12345678901 failed", "  - run: x # round 9",
                              "A: \"a # round 9\"", "U: https://x.test/p#round-9", "B: '#round 9'"))
  expect_equal(workflow$line, 2:3)
  expect_length(find_history(workflow$text[[1L]]), 1L)
  expect_length(find_history(workflow$text[[2L]]), 1L)

  described <- test_descriptions(c("test_that(\"F123 breaks\", {", "  expect_true(TRUE)", "})"))
  expect_equal(described$line, 1L)
  expect_length(find_history(described$text), 1L)
})

test_that("workflow comments of the repository hold no history", {
  skip_if_not(cellucid_is_source_checkout(), "the workflows are in the checkout only")
  root <- cellucid_repository_root()
  expect_gt(length(workflow_files(root)), 0L)
  problems <- unlist(lapply(workflow_files(root), history_in_file))
  if (length(problems)) fail(paste(c("history in workflow comments:", problems), collapse = "\n"))
  expect_length(problems, 0L)
})

test_that("comments and test descriptions of the package hold no history", {
  # The root is the unpacked sources under R CMD check as well
  # (helper-repository.R), so the scan runs there; tools/ and the workflows
  # are build-ignored and are scanned in the checkout only.
  root <- cellucid_repository_root()
  expect_true(dir.exists(file.path(root, "R")) && dir.exists(file.path(root, "src")))
  problems <- unlist(lapply(source_files(root), history_in_file))
  expect_gt(length(source_files(root)), 100L)
  expect_length(problems, 0L)
  if (length(problems)) fail(paste(c("history in comments:", problems), collapse = "\n"))
})
