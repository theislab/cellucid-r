# cellucid_apply_session() against the session file the web app's own
# encoders wrote: tests/testthat/fixtures/session_format3/ is a byte-identical
# copy of cellucid-python's tests/fixtures/session_format3/ (the same file its
# test_web_app_format3_file_reads_exactly applies to an AnnData), so both
# packages are held to one session and one expected result.

.session_fixture_sha256 <- c(
  "session.cellucid-session" = "d3bc22cef6310ca01adc5021895ac909776ad0cf2ef62decd9a9c10d918218b1",
  "expected.json" = "7b9a08e2e59500e35bab5d4c3b9eedb62cef8e5e9176cb4b3939c2b9811de7db"
)

session_fixture <- function(name = "session.cellucid-session") {
  testthat::test_path("fixtures", "session_format3", name)
}

session_expected <- function() jsonlite::read_json(session_fixture("expected.json"))

# The fixture dataset's cell names, cell_0 ... cell_<n - 1>.
format3_cells <- function(n = length(session_expected()$ctype)) paste0("cell_", seq_len(n) - 1L)

# The dataset the fixture session was saved on, as a SingleCellExperiment:
# 1024 cells, the category ctype, one gene, a 2D UMAP.
format3_sce <- function(ctype = unlist(session_expected()$ctype)) {
  n <- length(ctype)
  cells <- format3_cells(n)
  SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = Matrix::Matrix(0, 1L, n, sparse = TRUE, dimnames = list("gene_0", cells))),
    colData = S4Vectors::DataFrame(ctype = factor(ctype), row.names = cells),
    reducedDims = list(umap = matrix(0, n, 2L, dimnames = list(cells, NULL)))
  )
}

test_that("both repositories hold the same pinned session fixture", {
  for (name in names(.session_fixture_sha256)) {
    expect_identical(
      cellucid:::.native_sha256(readBin(session_fixture(name), "raw", file.size(session_fixture(name)))),
      .session_fixture_sha256[[name]]
    )
  }
  python <- file.path(cellucid_repository_root(), "..", "cellucid-python", "tests", "fixtures", "session_format3")
  skip_if_not(dir.exists(python), "cellucid-python is not checked out beside this repository")
  for (name in names(.session_fixture_sha256)) {
    expect_identical(unname(tools::md5sum(file.path(python, name))), unname(tools::md5sum(session_fixture(name))))
  }
})

test_that("a session of another format is refused naming its version, in the web app's words", {
  # A format-2 file's cell-order digest was taken over display coordinates,
  # which no later version can check (bundle/format.js otherSessionFormatMessage).
  expect_identical(cellucid:::.session_format_version, 3L)
  path <- tempfile(fileext = ".cellucid-session")
  for (case in list(
    list(2L, paste0("This session was saved by an earlier version of Cellucid (session format 2): ",
                    "the check that it marks the same cells cannot be verified, so it was not opened. ",
                    "Re-create it in this version of Cellucid.")),
    list(4L, paste0("This session was saved by a newer version of Cellucid (session format 4), ",
                    "which this version cannot read. Open it in that version of Cellucid."))
  )) {
    writeBin(c(charToRaw(paste0("CELLUCID_SESSION/", case[[1L]], "\n")), cellucid:::.encode_uint32(2L), charToRaw("{}")), path)
    expect_error(cellucid_apply_session(format3_sce(), path, expected_dataset_id = "format3"), case[[2L]], fixed = TRUE)
  }
  expect_identical(cellucid:::.session_signature_refusal(charToRaw("CELLUCID_SESSION/x\n")), cellucid:::.not_a_session_file_message)
})

test_that("the session's codecs read the web app's file exactly", {
  expected <- session_expected()
  bundle <- cellucid:::.read_session_bundle(session_fixture())
  expect_identical(bundle$manifest$datasetFingerprint$name, "Format 3 fixture")
  expect_identical(unlist(bundle$manifest$requires$obs), "ctype")
  expect_equal(cellucid:::.native_digest_cell_set(0:99), expected$digests$A)
  expect_equal(cellucid:::.native_digest_cell_set(100:899), expected$digests$B)
  expect_equal(
    cellucid:::.native_digest_cell_set(c(0, 1, 65535, 65536, 561946, 4294967295)),
    expected$digests$large
  )
  for (group in names(expected$encodings)) {
    payload <- cellucid:::.session_chunk(bundle, group)
    expect_identical(as.integer(payload[[1L]]), expected$encodings[[group]], info = group)
  }
})

test_that("a saved session's highlights and fields become the object's columns", {
  skip_without_sce()
  expected <- session_expected()
  object <- format3_sce()
  out <- cellucid_apply_session(object, session_fixture(), expected_dataset_id = "format3")
  added <- SummarizedExperiment::colData(out)
  for (group in names(expected$groups)) {
    expect_identical(which(added[[paste0("cellucid_highlight__", group)]]) - 1L, as.integer(unlist(expected$groups[[group]])), info = group)
  }
  regions <- c("R1", "R2", "rest")
  expect_identical(as.character(added$regions), regions[unlist(expected$regionsCodes) + 1L])
  expect_identical(levels(added$regions), regions)
  record <- S4Vectors::metadata(out)$cellucid$session
  expect_identical(record$applied, list(expected_dataset_id = "format3", contract = "exact"))
  # Where the viewer opened the dataset (a link, a server, a file name) never
  # travels with the object.
  expect_false("source" %in% names(record$manifest))
  expect_true(all(c("createdAt", "datasetFingerprint", "requires", "chunks") %in% names(record$manifest)))
  expect_setequal(names(record$materialized_fields), c(paste0("cellucid_highlight__", names(expected$groups)), "regions"))
  # The caller's object is not changed.
  expect_false("regions" %in% names(SummarizedExperiment::colData(object)))
  summary <- cellucid_apply_session(object, session_fixture(), expected_dataset_id = "format3", return_summary = TRUE)
  expect_setequal(summary$added_columns, names(record$materialized_fields))
})

test_that("a Seurat object gets the same columns in meta.data", {
  skip_without_seurat()
  expected <- session_expected()
  counts <- Matrix::Matrix(0, 1L, length(expected$ctype), sparse = TRUE, dimnames = list("gene_0", format3_cells()))
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  object$ctype <- factor(unlist(expected$ctype))
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    matrix(0, ncol(object), 2L, dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))), key = "UMAP_", assay = "RNA"
  )
  # The fixture was saved on a store without the object's own metadata
  # columns: the session asks only for ctype.
  out <- cellucid_apply_session(object, session_fixture(), expected_dataset_id = "format3")
  meta <- methods::slot(out, "meta.data")
  expect_identical(which(meta$cellucid_highlight__highlight_3) - 1L, as.integer(unlist(expected$groups$highlight_3)))
  expect_identical(as.character(meta$regions), c("R1", "R2", "rest")[unlist(expected$regionsCodes) + 1L])
  expect_identical(methods::slot(out, "misc")$cellucid$session$applied$expected_dataset_id, "format3")
})

test_that("a session that does not match the object is refused before anything is written", {
  skip_without_sce()
  expected <- session_expected()
  expect_error(
    cellucid_apply_session(format3_sce(), session_fixture()),
    "expected_dataset_id is required"
  )
  expect_error(
    cellucid_apply_session(format3_sce(), session_fixture(), expected_dataset_id = "other"),
    "datasetId 'format3' != expected_dataset_id 'other'"
  )
  moved <- unlist(expected$ctype)
  moved[[1L]] <- "B"
  expect_error(
    cellucid_apply_session(format3_sce(moved), session_fixture(), expected_dataset_id = "format3"),
    "no longer holds the same cells"
  )
  missing <- format3_sce()
  SummarizedExperiment::colData(missing)$ctype <- NULL
  expect_error(
    cellucid_apply_session(missing, session_fixture(), expected_dataset_id = "format3"),
    "This session uses the obs column 'ctype', which this dataset does not have."
  )
  taken <- format3_sce()
  SummarizedExperiment::colData(taken)$regions <- 1
  expect_error(
    cellucid_apply_session(taken, session_fixture(), expected_dataset_id = "format3"),
    "Column already exists: regions"
  )
  fewer <- format3_sce()[, 1:10]
  expect_error(
    cellucid_apply_session(fewer, session_fixture(), expected_dataset_id = "format3"),
    "cellCount 1024 != the object's 10 cells"
  )
})

test_that("a file that is not a session is refused with the viewer's sentence", {
  path <- tempfile(fileext = ".cellucid-session")
  on.exit(unlink(path))
  writeBin(charToRaw("NOT_A_SESSION_FILE_AT_ALL\n"), path)
  expect_error(cellucid:::.read_session_bundle(path), cellucid:::.not_a_session_file_message, fixed = TRUE)
})

test_that("a gzip stream is inflated exactly, or refused", {
  inflate <- cellucid:::.native_inflate
  member <- memCompress(charToRaw(strrep("A", 100)), "gzip")
  # memCompress writes zlib; a gzip member is made by gzcon.
  gz <- function(text) {
    path <- tempfile()
    con <- gzfile(path, "wb")
    writeBin(charToRaw(text), con)
    close(con)
    on.exit(unlink(path))
    readBin(path, "raw", file.size(path))
  }
  one <- gz(strrep("A", 100))
  expect_identical(rawToChar(inflate(one, TRUE, 100)), strrep("A", 100))
  expect_identical(rawToChar(inflate(member, FALSE)), strrep("A", 100))
  expect_error(inflate(c(one, gz("B")), TRUE, 100), "trailing or concatenated")
  expect_error(inflate(c(one, as.raw(1:3)), TRUE, 100), "trailing or concatenated")
  expect_error(inflate(one[seq_len(length(one) - 6L)], TRUE, 100), "truncated|Invalid gzip")
  expect_error(inflate(one, TRUE, 99), "smaller than the gzip payload")
  expect_error(inflate(one, TRUE, 101), "does not match the decoded payload")
})

test_that("the shared conformance fixture's mapping rules hold in R", {
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  for (case in fixture$dataMapping) {
    if (!is.null(case$expected$error)) {
      condition <- tryCatch(cellucid:::.validate_data_mapping(case$value), cellucid_data_mapping_error = identity)
      expect_s3_class(condition, "cellucid_data_mapping_error")
      expect_identical(condition$code, case$expected$error, info = case$name)
      next
    }
    canonical <- cellucid:::.validate_data_mapping(case$value)
    expect_identical(names(canonical), names(case$expected), info = case$name)
  }
  for (case in fixture$mappingModality) {
    modalities <- if (is.null(case$modalities)) NULL else unlist(case$modalities)
    mapping <- if (length(case$mapping) == 0L) structure(list(), names = character()) else case$mapping
    if (is.null(case$expected)) {
      expect_silent(cellucid:::.require_mapping_modality(mapping, modalities))
      next
    }
    condition <- tryCatch(cellucid:::.require_mapping_modality(mapping, modalities), cellucid_data_mapping_error = identity)
    expect_identical(condition$code, case$expected$error, info = case$name)
  }
})

test_that("a session saved on a subsample applies at the original positions, the others NA", {
  skip_without_seurat()
  expected <- session_expected()
  k <- length(expected$ctype)
  n <- k + 76L
  # The fixture session, as if saved on the first k of n cells.
  bytes <- readBin(session_fixture(), "raw", file.size(session_fixture()))
  magic <- length(cellucid:::.session_magic)
  manifest_length <- readBin(bytes[magic + 1:4], "integer", size = 4L, endian = "little")
  manifest <- jsonlite::fromJSON(rawToChar(bytes[magic + 4L + seq_len(manifest_length)]), simplifyVector = FALSE)
  kept <- seq_len(k) - 1
  manifest$cells <- list(
    spec = list(method = "uniform", nObs = n, k = k, seed = 0L),
    kept = gsub("\n", "", jsonlite::base64_enc(cellucid:::.subsample_encode_cell_set(kept)), fixed = TRUE),
    digest = cellucid:::.subsample_digest(n, kept)
  )
  encoded <- charToRaw(as.character(jsonlite::toJSON(manifest, auto_unbox = TRUE, null = "null", digits = NA)))
  path <- tempfile(fileext = ".cellucid-session")
  writeBin(c(bytes[seq_len(magic)], cellucid:::.encode_uint32(length(encoded)), encoded,
             bytes[(magic + 4L + manifest_length + 1L):length(bytes)]), path)

  cells <- c(format3_cells(), paste0("extra_", seq_len(76L)))
  counts <- Matrix::Matrix(0, 1L, n, sparse = TRUE, dimnames = list("gene_0", cells))
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  object$ctype <- factor(c(unlist(expected$ctype), rep(unlist(expected$ctype)[[1L]], 76L)))
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    matrix(0, n, 2L, dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))), key = "UMAP_", assay = "RNA"
  )
  out <- cellucid_apply_session(object, path, expected_dataset_id = "format3")
  meta <- methods::slot(out, "meta.data")
  # The session's cell i is original position kept[i]: per-cell chunks name
  # cells by original position, ascending, whatever order they were open in.
  expect_identical(sort(which(meta$cellucid_highlight__highlight_3) - 1L),
                   sort(as.integer(kept[unlist(expected$groups$highlight_3) + 1L])))
  expect_true(all(is.na(meta$cellucid_highlight__highlight_3[(k + 1L):n])))
  expect_identical(as.character(meta$regions[kept + 1L]), c("R1", "R2", "rest")[unlist(expected$regionsCodes) + 1L])
  expect_true(all(is.na(meta$regions[(k + 1L):n])))

  manifest$cells$digest <- strrep("0", 16L)
  encoded <- charToRaw(as.character(jsonlite::toJSON(manifest, auto_unbox = TRUE, null = "null", digits = NA)))
  writeBin(c(bytes[seq_len(magic)], cellucid:::.encode_uint32(length(encoded)), encoded,
             bytes[(magic + 4L + manifest_length + 1L):length(bytes)]), path)
  expect_error(cellucid_apply_session(object, path, expected_dataset_id = "format3"), "digest")
})
