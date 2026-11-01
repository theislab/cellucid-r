# The shared data-mapping conformance fixture, run against this package.
#
# tests/testthat/fixtures/data_mapping_conformance.json is a byte-identical
# copy of cellucid-python's tests/fixtures/data_mapping_conformance.json, which
# the web app also carries; all three repositories pin the same SHA-256.
#
# A writer applies four of its rules: making gene names unique
# (makeUniqueNames), the one-line note that reports what that changed
# (uniqueNamesNote), the vector-field key rule (vectorFieldKeys), by which
# cellucid_prepare() names each vector field '<field>_<basis>' and checks that
# the viewer reads it back as that field on 'X_<basis>' (in its modality's
# namespace), and the default neighbour graph (defaultGraph), by which it
# chooses the graph it writes. All four are run here case for case. The other sections pin rules this package has no counterpart
# for: the best-guess embedding and expression choices, the representations
# enabled when a file opens, the default vector bindings, the "nothing to
# plot" text, DataMapping validation, and the check of a mapping's modality
# (mappingModality) are what a reader does with an AnnData or MuData file; and
# embeddingColumns and canonicalEmbeddingId describe the
# 0-based, typed column selection of cellucid-python's guessing API.
# cellucid_prepare() takes caller-named embedding ids and 1-based columns
# instead.

.conformance_sha256 <- "6f7b7ff0e4392ea00c7e21c7a868ec68d3ec83a63558227f9f7d2c300a666c32"

test_that("both repositories hold the same pinned conformance fixture", {
  expect_pinned_fixture("data_mapping_conformance.json", .conformance_sha256)
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  expect_identical(fixture$version, 7L)
  expect_identical(
    names(fixture),
    c(
      "_about", "version", "makeUniqueNames", "uniqueNamesNote",
      "bestGuessTier", "embeddingSlots", "enabledEmbeddings",
      "embeddingColumns", "canonicalEmbeddingId", "defaultExpression",
      "vectorFieldKeys", "defaultVectorBindings", "nothingPlottable",
      "dataMapping", "defaultModality", "mappingModality", "defaultGraph"
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
        found[[length(found) + 1L]] <- list(
          key = key, field = parsed$field, basis = parsed$basis, embeddingKey = parsed$embedding
        )
      }
    }
    expected <- lapply(case$expected, function(entry) {
      list(
        key = enc2utf8(entry$key), field = enc2utf8(entry$field), basis = enc2utf8(entry$basis),
        embeddingKey = enc2utf8(entry$embeddingKey)
      )
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

test_that("the default neighbour graph matches every case of the shared conformance fixture", {
  # cellucid_prepare() writes the graph this rule chooses (Seurat graphs, SCE
  # colPairs) as obsp['connectivities'], the one the viewer opens with.
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  expect_gt(length(fixture$defaultGraph), 0L)
  for (case in fixture$defaultGraph) {
    keys <- as.character(unlist(case$keys))
    expect_identical(cellucid:::.default_graph(keys, case$modality), case$expected, info = case$name)
  }
})

test_that("a self-loop is no edge: a graph given with a unit diagonal (Seurat SNN) is kept without it", {
  # The web app and cellucid-python leave the diagonal out the same way.
  ring <- Matrix::sparseMatrix(i = c(1:5, 2:6), j = c(2:6, 1:5), x = 0.5, dims = c(6L, 6L))
  graph <- cellucid:::.resolve_connectivities(ring + Matrix::Diagonal(6), 6L)
  expect_identical(sum(Matrix::diag(graph)), 0)
  expect_identical(as.matrix(graph), as.matrix(methods::as(ring, "generalMatrix")))
  expect_error(
    cellucid:::.resolve_connectivities(Matrix::sparseMatrix(i = 1, j = 2, x = 1, dims = c(2L, 2L)), 2L),
    "exactly symmetric: a directed neighbour graph"
  )
})
