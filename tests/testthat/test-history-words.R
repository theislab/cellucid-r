# A history in words: the methods paragraph and the pseudocode, by the web
# app's rules (app/history/exports/{methods,pseudocode}.js) -- cellucid-python's
# tests/test_history_words.py, string for string.
#
# The golden history pins the text of three branches of its annotation tracker
# byte for byte (test-sessions.R). These tests cover the wording rules the
# golden branches do not reach, on the golden tracker's own steps changed: a
# squashed step of each cause (recording paused, steps merged, the history
# cleared, a state shared) listed in the methods and kept as a pseudocode
# comment, the sentence that says a squashed step holds only its result,
# everyone who acted while recording was paused, and the "Before"/"After"
# headings around a group.

words_tracker <- function() {
  cellucid_tracker(cellucid_read_session(testthat::test_path("fixtures", "history_golden", "history_golden.cellucid-session")), "Annotation")
}
# The Scratch branch up to its paused period (s14).
words_path <- function(tracker) .history_path(tracker, "s14")
words_dataset <- list(name = "History golden", cells = 600L, genes = 4L)
words_opening <- paste0(
  "Single-cell data were explored interactively in Cellucid. The dataset “History golden” (600 cells and 4 genes) was",
  " opened in full."
)
words_complete <- "The complete history, with every parameter, can be reopened from the accompanying Cellucid session file."

test_that("squash phrases", {
  phrase <- function(cause, count) .history_squash_phrase(list(cause = cause, count = count, by = list(), carried = list(notes = list(), bookmarks = list())))
  expect_identical(phrase("paused", NULL), "Changes made with recording paused")
  expect_identical(phrase("steps", 3L), "3 steps merged into one")
  expect_identical(phrase("clear", 9L), "A fresh start from the state reached before the history was cleared")
  expect_identical(phrase("share", NULL), "A shared state, received without its history")
})

test_that("a paused period that changed the data is listed and said plainly", {
  tracker <- words_tracker()
  scratch <- words_path(tracker)
  # Recording was paused while the filters and the colours changed (the person and a notebook acted).
  scratch[[length(scratch)]]$touches <- list("filters", "colors")
  text <- .history_methods_paragraph(scratch, words_dataset, NULL, list(), tracker$groups)
  expect_identical(text, paste(c(
    paste0(words_opening, " The 6 steps applied are listed below, in order.",
           " For step 6, only the state it reached was kept, not the operations that led to it.",
           " 1 further step only changed the camera, the panel layout, colors or display settings and is not listed. ",
           words_complete),
    "Before QC\nThe following 1 step was applied: (1) renamed NK cell to NK in cell_type.",
    paste0("QC\nThe following 2 steps were applied, in this order: (2) highlighted 40 cells as \"Doublets\" on page \"QC\";",
           " (3) added the continuous column score. Steps 2–3 were taken by a Python notebook through the Cellucid viewer API."),
    paste0("After QC\nThe following 3 steps were applied, in this order: (4) ran Marker genes; (5) kept the gene set T cell",
           " markers (2 genes); (6) changes made with recording paused (filters, colors). Steps 4–5 were taken by an AI",
           " assistant (Claude, claude-sonnet-4-6) at the user's request, and step 6 by a Python notebook through the",
           " Cellucid viewer API.")
  ), collapse = "\n\n"))
})

test_that("each squash cause in the methods and the pseudocode", {
  tracker <- words_tracker()
  scratch <- words_path(tracker)
  root <- scratch[[1L]]
  paused <- scratch[[length(scratch)]]
  cases <- list(
    steps = list(3L, "3 steps merged into one"),
    clear = list(5L, "a fresh start from the state reached before the history was cleared"),
    share = list(NULL, "a shared state, received without its history")
  )
  for (cause in names(cases)) {
    said <- cases[[cause]][[2L]]
    step <- paused
    step$parent <- root$id
    step["group"] <- list(NULL)
    step["bookmark"] <- list(NULL)
    step$touches <- list("field-overlays", "highlights")
    step$squashed <- list(cause = cause, count = cases[[cause]][[1L]], by = if (identical(cause, "share")) list() else list(paused$by),
                         carried = list(notes = list(), bookmarks = list()))
    expect_identical(.history_methods_paragraph(list(root, step), words_dataset, NULL, list()), paste0(
      words_opening, " The following 1 step was applied: (1) ", said, " (fields and labels, highlights).",
      " For step 1, only the state it reached was kept, not the operations that led to it. ", words_complete
    ))
    expect_identical(.history_pseudocode_text(list(root, step), NULL, NULL, "T", "B"), paste0(
      "# Cellucid\n# Tracker \"T\", branch \"B\": 1 step\n",
      "# ", toupper(substring(said, 1L, 1L)), substring(said, 2L), ": Fields and labels, Highlights (values only)\n"
    ))
  }
})

test_that("a squashed step that only presented is counted", {
  tracker <- words_tracker()
  path <- words_path(tracker)
  text <- .history_methods_paragraph(path, words_dataset, NULL, list(), tracker$groups)
  expect_true(grepl("2 further steps only changed the camera, the panel layout, colors or display settings and are not listed.", text, fixed = TRUE))
  expect_false(grepl("recording paused", text, fixed = TRUE))
  expect_true(endsWith(.history_pseudocode_text(path, NULL, NULL, "T", "B", tracker$groups),
                       "# Changes made with recording paused: Colors, Settings (values only)\n"))
})
