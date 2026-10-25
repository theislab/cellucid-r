# A count and its noun agree: "1 cell", "2 cells", never "1 cells".
#
# Every message that writes a count before a noun takes the noun from
# .plural(); paste0(n, " cells") reads "1 cells" the day n is one. The
# package's sources are parsed, and an argument " cells" (a plural noun after
# a space) that follows a count in the same call is refused. What is not a
# count is not refused: a label before a plural noun (paste0(label, " fields:")).

plural_nouns <- c(
  "cells", "points", "genes", "values", "entries", "categories", "options", "columns", "rows",
  "pages", "items", "fields", "characters", "edges", "datasets", "views", "codes", "groups",
  "pixels", "indices", "stores", "positions", "frames", "samples", "records", "files",
  "workers", "results", "arrays", "layers", "features", "suggestions", "requests", "steps",
  "objects", "modalities", "embeddings", "annotators", "voters", "votes", "chunks", "parts",
  "matrices", "sets", "nodes", "users", "labels", "lines", "errors", "warnings", "retries",
  "attempts", "bytes", "digits", "dimensions", "axes", "keys", "components", "messages",
  "times", "reads", "assays", "reductions", "graphs"
)
noun_argument <- paste0("^[\"'] (", paste(plural_nouns, collapse = "|"), ")\\b")
label_like <- "^(where|what|label|name|key|type|kind|path|prefix|item\\$label|header\\$type|field|context|position - 1L)$"

count_noun_offenders <- function(file) {
  parsed <- parse(file, keep.source = TRUE)
  data <- utils::getParseData(parsed, includeText = TRUE)
  data <- data[order(data$line1, data$col1), ]
  strings <- data[data$token == "STR_CONST" & grepl(noun_argument, data$text), ]
  offenders <- character()
  for (index in seq_len(nrow(strings))) {
    holder <- strings$parent[[index]]
    call_id <- data$parent[data$id == holder]
    siblings <- data[data$parent == call_id, ]
    siblings <- siblings[order(siblings$line1, siblings$col1), ]
    at <- which(siblings$id == holder)
    if (at < 3L || siblings$token[[at - 1L]] != "','") next
    previous <- siblings[at - 2L, ]
    if (previous$token != "expr") next
    text <- utils::getParseText(data, previous$id)
    if (grepl("^[\"']", text) || grepl(label_like, text)) next
    offenders <- c(offenders, paste0(basename(file), ":", strings$line1[[index]], ": ", text, ", ", strings$text[[index]]))
  }
  offenders
}

test_that(".plural takes the singular for one", {
  expect_identical(.plural(1, "cell"), "cell")
  expect_identical(.plural(1L, "cell"), "cell")
  expect_identical(.plural(0, "cell"), "cells")
  expect_identical(.plural(2, "entry", "entries"), "entries")
  expect_identical(.plural(NA, "cell"), "cells")
  expect_identical(.plural(NULL, "cell"), "cells")
})

test_that("no count is written before a plural noun it may not agree with", {
  root <- test_path("..", "..", "R")
  skip_if_not(dir.exists(root), "the package sources are not beside the tests")
  offenders <- unlist(lapply(sort(list.files(root, pattern = "\\.R$", full.names = TRUE)), count_noun_offenders))
  expect_identical(offenders, character(), info = "use .plural() for the noun")
})

test_that("the check catches a hard-coded plural", {
  probe <- tempfile(fileext = ".R")
  on.exit(unlink(probe), add = TRUE)
  writeLines(c(
    'f <- function(n, label) stop("has ", n, " cells and ", length(n), " genes.")',
    'g <- function(label) paste0(label, " fields: none")'
  ), probe)
  expect_length(count_noun_offenders(probe), 2L)
})
