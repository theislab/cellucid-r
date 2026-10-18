# One export, written by both writers, byte for byte.
#
# tests/testthat/fixtures/writer_parity_v3.json is a byte-identical copy of
# cellucid-python's tests/fixtures/writer_parity_v3.json: the inputs of one
# small export -- two cells, five embeddings including a wide array read at two
# of its columns and a plain vector, two vector fields bound to embedding ids,
# a neighbour graph, and a var whose default and alternative gene names both
# repeat and have empty values -- together with every payload byte and every
# manifest cellucid-python wrote for them. This writer exports the same inputs
# and must produce exactly those bytes and manifests. Both repositories pin the
# file's SHA-256, so neither copy can change alone.

.parity_v3_sha256 <- "44779fb1fce0955d6bdbdfc85df9458c8f7edff945cec07a4c9c65fbc7bb9aad"
.conformance_sha256 <- "b5b35f7d86f4472d4fc8117c49ada6341aaf9834998bff778ef8d4a07bc4fb43"

.fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

.file_sha256 <- function(path) {
  # tools::md5sum is the only digest base R ships, so the SHA-256 pin is
  # computed by the system's shasum/sha256sum when one exists.
  for (tool in c("shasum", "sha256sum")) {
    binary <- Sys.which(tool)
    if (nzchar(binary)) {
      arguments <- if (identical(tool, "shasum")) c("-a", "256", shQuote(path)) else shQuote(path)
      output <- system2(binary, arguments, stdout = TRUE)
      return(strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]])
    }
  }
  NA_character_
}

.read_fixture <- function(name) {
  jsonlite::fromJSON(.fixture_path(name), simplifyVector = FALSE)
}

.as_number_matrix <- function(rows) {
  matrix(
    as.numeric(unlist(rows)),
    nrow = length(rows),
    byrow = TRUE
  )
}

.parity_v3_arguments <- function(inputs, out) {
  obs <- data.frame(row.names = seq_along(inputs$obs[[1L]]$values))
  for (key in names(inputs$obs)) {
    field <- inputs$obs[[key]]
    values <- unlist(field$values)
    obs[[key]] <- if (identical(field$kind, "category")) factor(values) else as.numeric(values)
  }
  var_columns <- lapply(inputs$var$columns, function(values) {
    vapply(values, function(value) if (is.null(value)) NA_character_ else value, character(1))
  })
  var <- as.data.frame(var_columns, stringsAsFactors = FALSE)
  rownames(var) <- unlist(inputs$var$index)

  embeddings <- list()
  for (key in names(inputs$embeddings)) {
    item <- inputs$embeddings[[key]]
    values <- if (is.list(item$values[[1L]])) {
      .as_number_matrix(item$values)
    } else {
      as.numeric(unlist(item$values))
    }
    embeddings[[key]] <- cellucid_embedding(
      values,
      columns = if (is.null(item$columns)) NULL else as.integer(unlist(item$columns)) + 1L,
      label = item$label,
      source_key = item$source_key
    )
  }
  vector_fields <- lapply(inputs$vector_fields, function(bindings) {
    lapply(bindings, .as_number_matrix)
  })
  list(
    latent_space = .as_number_matrix(inputs$latent_space),
    obs = obs,
    var = var,
    gene_expression = .as_number_matrix(inputs$gene_expression),
    var_gene_id_column = inputs$var_gene_id_column,
    var_gene_name_columns = unlist(inputs$var_gene_name_columns),
    connectivities = .as_number_matrix(inputs$connectivities),
    out_dir = out,
    centroid_min_points = inputs$centroid_min_points,
    force = TRUE,
    obs_categorical_dtype = inputs$obs_categorical_dtype,
    dataset_name = inputs$dataset_name,
    dataset_id = inputs$dataset_id,
    created_at = inputs$created_at,
    embeddings = embeddings,
    embedding_defaults = unlist(inputs$embedding_defaults),
    vector_fields = vector_fields,
    vector_field_default = inputs$vector_field_default
  )
}

.export_relative_files <- function(root) {
  sort(list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE))
}

.raw_hex <- function(path) {
  bytes <- readBin(path, "raw", file.size(path))
  paste(as.character(bytes), collapse = "")
}

test_that("both repositories hold the same pinned fixtures", {
  for (pin in list(
    list("writer_parity_v3.json", .parity_v3_sha256),
    list("data_mapping_conformance.json", .conformance_sha256)
  )) {
    digest <- .file_sha256(.fixture_path(pin[[1L]]))
    skip_if(is.na(digest), "no SHA-256 tool on this system")
    expect_identical(digest, pin[[2L]])
  }
})

test_that("the v3 parity fixture exports to exactly cellucid-python's bytes and manifests", {
  fixture <- .read_fixture("writer_parity_v3.json")
  out <- file.path(tempdir(), "cellucid_parity_v3")
  unlink(out, recursive = TRUE, force = TRUE)
  arguments <- .parity_v3_arguments(fixture$inputs, out)
  # The fixture's default and alternative gene names both repeat and both
  # hold an empty value, which is reported once per name set.
  expect_message(
    expect_message(do.call(cellucid_prepare, arguments), "var column 'gene_name'"),
    "var column 'symbol'"
  )

  expected_payloads <- fixture$expected$payloads
  expected_manifests <- fixture$expected$manifests
  expect_identical(
    .export_relative_files(out),
    sort(c(names(expected_payloads), names(expected_manifests)))
  )
  for (relative in names(expected_payloads)) {
    expect_identical(
      .raw_hex(file.path(out, relative)),
      expected_payloads[[relative]],
      info = relative
    )
  }
  for (relative in names(expected_manifests)) {
    written <- jsonlite::fromJSON(file.path(out, relative), simplifyVector = FALSE)
    expected <- expected_manifests[[relative]]
    if (identical(relative, "dataset_identity.json")) {
      expect_identical(
        written$cellucid_data_version,
        as.character(utils::packageVersion("cellucid"))
      )
      written$cellucid_data_version <- NULL
    }
    expect_null(
      cellucid:::.first_json_difference(expected, written),
      info = relative
    )
  }
})

test_that("make-unique matches every case of the shared conformance fixture", {
  fixture <- .read_fixture("data_mapping_conformance.json")
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
  fixture <- .read_fixture("data_mapping_conformance.json")
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
