# The shared data-mapping conformance fixture, run against this package.
#
# tests/testthat/fixtures/data_mapping_conformance.json is a byte-identical
# copy of cellucid-python's tests/fixtures/data_mapping_conformance.json, which
# the web app also carries; all three repositories pin the same SHA-256.
#
# A writer applies three of its rules: making gene names unique
# (makeUniqueNames), the one-line note that reports what that changed
# (uniqueNamesNote), and the vector-field key rule (vectorFieldKeys), by which
# cellucid_prepare() names each vector field '<field>_<basis>' and checks that
# the viewer reads it back as that field on 'X_<basis>'. All three are run here
# case for case. The other sections pin rules this package has no counterpart
# for: the best-guess embedding and expression choices, the representations
# enabled when a file opens, the default vector bindings, the "nothing to
# plot" text, DataMapping validation, and the check of a mapping's modality
# (mappingModality) are what a reader does with an AnnData or MuData file; and
# embeddingColumns and canonicalEmbeddingId describe the
# 0-based, typed column selection of cellucid-python's guessing API.
# cellucid_prepare() takes caller-named embedding ids and 1-based columns
# instead.

.conformance_sha256 <- "cfc84c8c332f54824c6d8a438ee0c8ae014852674d7ae6bdc19dff1a8c84089e"

test_that("both repositories hold the same pinned conformance fixture", {
  expect_pinned_fixture("data_mapping_conformance.json", .conformance_sha256)
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  expect_identical(fixture$version, 6L)
  expect_identical(
    names(fixture),
    c(
      "_about", "version", "makeUniqueNames", "uniqueNamesNote",
      "bestGuessTier", "embeddingSlots", "enabledEmbeddings",
      "embeddingColumns", "canonicalEmbeddingId", "defaultExpression",
      "vectorFieldKeys", "defaultVectorBindings", "nothingPlottable",
      "dataMapping", "defaultModality", "mappingModality"
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

test_that("the vector-field key rule matches every case of the shared conformance fixture", {
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  for (case in fixture$vectorFieldKeys) {
    keys <- vapply(case$keys, enc2utf8, character(1))
    found <- list()
    for (key in keys) {
      parsed <- cellucid:::.parse_vector_field_key(key, keys)
      if (!is.null(parsed)) {
        expect_identical(paste0(parsed$field, "_", parsed$basis), key, info = case$name)
        found[[length(found) + 1L]] <- list(key = key, field = parsed$field, basis = parsed$basis)
      }
    }
    expected <- lapply(case$expected, function(entry) {
      list(key = enc2utf8(entry$key), field = enc2utf8(entry$field), basis = enc2utf8(entry$basis))
    })
    expect_identical(found, expected, info = case$name)
  }
})

test_that("the default modality matches every case of the shared conformance fixture", {
  # The modality a multimodal dataset opens with; cellucid_prepare() and the
  # viewer apply it to choose the modality the gene arguments describe.
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  expect_gt(length(fixture$defaultModality), 0L)
  for (case in fixture$defaultModality) {
    names <- as.character(unlist(case$names))
    expect_identical(cellucid:::.default_modality(names), case$expected, info = paste(names, collapse = ","))
  }
})
