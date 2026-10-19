# The shared data-mapping conformance fixture, run against this package.
#
# tests/testthat/fixtures/data_mapping_conformance.json is a byte-identical
# copy of cellucid-python's tests/fixtures/data_mapping_conformance.json, which
# the web app also carries; all three repositories pin the same SHA-256.
#
# A writer applies two of its rules: making gene names unique
# (makeUniqueNames) and the one-line note that reports what that changed
# (uniqueNamesNote). Both are run here case for case. The other sections pin
# rules this package has no counterpart for: the best-guess embedding and
# expression choices, the representations enabled when a file opens, the
# default vector bindings, the "nothing to plot" text, and DataMapping
# validation are what a reader does with an AnnData file; exportDataMapping is
# what the web app shows an export with, and its cases include an export
# declared the way cellucid_prepare() declares it (caller-chosen ids, one 2D
# default), which the reader treats exactly like a Python export; and
# embeddingColumns and canonicalEmbeddingId describe the 0-based, typed column
# selection of cellucid-python's guessing API. cellucid_prepare() takes
# caller-named embedding ids and 1-based columns instead.

.conformance_sha256 <- "dd08f84297c179b37e21a09dd820781cbfa21bf33f063bd79f925ec44a026660"

test_that("both repositories hold the same pinned conformance fixture", {
  expect_pinned_fixture("data_mapping_conformance.json", .conformance_sha256)
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  expect_identical(fixture$version, 2L)
  expect_identical(
    names(fixture),
    c(
      "_about", "version", "makeUniqueNames", "uniqueNamesNote",
      "bestGuessTier", "embeddingSlots", "enabledEmbeddings",
      "embeddingColumns", "canonicalEmbeddingId", "defaultExpression",
      "defaultVectorBindings", "nothingPlottable", "dataMapping",
      "exportDataMapping"
    )
  )
})

test_that("make-unique matches every case of the shared conformance fixture", {
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  for (case in fixture$makeUniqueNames) {
    index <- if (is.null(case$index)) {
      NULL
    } else {
      vapply(case$index, function(value) if (is.null(value)) NA_character_ else value, character(1))
    }
    run <- function() {
      cellucid:::.make_unique_names(case$names, index = index, join = case$join)
    }
    if (!is.null(case$expected$error)) {
      condition <- tryCatch(run(), cellucid_data_mapping_error = identity)
      expect_s3_class(condition, "cellucid_data_mapping_error")
      expect_identical(condition$code, case$expected$error, info = case$name)
      next
    }
    result <- run()
    expect_identical(result$names, as.character(unlist(case$expected$names)), info = case$name)
    expect_identical(anyDuplicated(result$names), 0L, info = case$name)
    expect_identical(
      lapply(result$renamed, as.character),
      lapply(case$expected$renamed, function(pair) as.character(unlist(pair))),
      info = case$name
    )
    expect_identical(
      as.integer(result$filled_from_index),
      as.integer(case$expected$filledFromIndex),
      info = case$name
    )
  }
})

test_that("the make-unique note matches the shared conformance fixture word for word", {
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  for (case in fixture$uniqueNamesNote) {
    renamed <- lapply(case$renamed, function(pair) as.character(unlist(pair)))
    note <- cellucid:::.unique_names_note(case$filledFromIndex, renamed, case$source)
    if (is.null(case$expected)) {
      expect_null(note)
    } else {
      expect_identical(note, enc2utf8(case$expected))
    }
  }
})
