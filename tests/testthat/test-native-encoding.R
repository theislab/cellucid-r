# Text is the same in every R session, whatever its native encoding
# (R/encoding.R). Windows R before 4.2 runs in a code page (CP1252 in Western
# Europe and the Americas), and R makes text native wherever the C level takes
# it as plain bytes. A session like that is made here by switching LC_CTYPE
# to "C" for the length of a block: only ASCII is native there, so every
# other character crosses those borders.

store_path <- function(name = "store.zarr") file.path(tempfile("cellucid_encoding_"), name)

with_non_utf8_session <- function(code) {
  old <- Sys.getlocale("LC_CTYPE")
  on.exit(Sys.setlocale("LC_CTYPE", old), add = TRUE)
  switched <- suppressWarnings(Sys.setlocale("LC_CTYPE", "C"))
  if (!nzchar(switched) || isTRUE(l10n_info()[["UTF-8"]])) {
    testthat::skip("this R cannot switch to a locale whose native encoding is not UTF-8")
  }
  force(code)
}

# A string as Windows R before 4.2 holds text from a CP1252 file or console.
latin1 <- function(text) iconv(text, "UTF-8", "latin1")

# Every message a call emits, muffled, in order.
messages_of <- function(code) {
  said <- character()
  value <- withCallingHandlers(code, message = function(condition) {
    said <<- c(said, conditionMessage(condition))
    invokeRestart("muffleMessage")
  })
  list(value = value, messages = said)
}

test_that("messages, warnings and errors carry their text exactly in a session that is not UTF-8", {
  text <- "e.g. 'A' \u2192 'A-1', a gene named \u03b1"
  with_non_utf8_session({
    expect_identical(tryCatch(cellucid:::.abort("x ", text), error = conditionMessage), paste0("x ", text))
    expect_identical(tryCatch(cellucid:::.warn("x ", text), warning = conditionMessage), paste0("x ", text))
    expect_identical(tryCatch(cellucid:::.inform("x ", text), message = conditionMessage), paste0("x ", text, "\n"))
  })
  expect_silent(suppressMessages(cellucid:::.inform("quiet")))
  # Arguments are joined as message(), warning() and stop() join them.
  expect_identical(tryCatch(cellucid:::.abort("a", c("b", "c"), 1L, NULL), error = conditionMessage), "abc1")
  error <- tryCatch(cellucid:::.abort("x"), error = identity)
  expect_s3_class(error, "simpleError")
  expect_null(conditionCall(error))
  expect_null(conditionCall(tryCatch(cellucid:::.warn("x"), warning = identity)))
})

test_that("every message, warning and error of the package is raised as a condition object", {
  root <- test_path("..", "..", "R")
  skip_if_not(dir.exists(root), "the package sources are not beside the tests")
  constructors <- c("simpleError", "simpleWarning", "simpleMessage", ".data_mapping_error", ".protocol_error")
  offenders <- character()
  for (file in sort(list.files(root, pattern = "\\.R$", full.names = TRUE))) {
    data <- utils::getParseData(parse(file, keep.source = TRUE), includeText = TRUE)
    calls <- data[data$token == "SYMBOL_FUNCTION_CALL" & data$text %in% c("stop", "warning", "message"), ]
    for (index in seq_len(nrow(calls))) {
      head <- calls$parent[[index]]
      # x$stop(): a method of x, not base's.
      if (any(data$parent == head & data$token %in% c("'$'", "'@'"))) next
      call <- str2lang(utils::getParseText(data, data$parent[data$id == head]))
      arguments <- as.list(call)[-1L]
      condition <- length(arguments) == 1L && is.null(names(arguments)) && (
        identical(arguments[[1L]], as.name("error")) ||
          (is.call(arguments[[1L]]) && as.character(arguments[[1L]][[1L]])[[1L]] %in% constructors)
      )
      if (!condition) offenders <- c(offenders, paste0(basename(file), ":", calls$line1[[index]]))
    }
  }
  expect_identical(offenders, character(), info = "raise it with .abort(), .warn() or .inform() (R/encoding.R)")
})

test_that("a store is written byte for byte the same, with the same notes, in a session that is not UTF-8", {
  skip_if_not_installed("Matrix")
  arguments <- function(out_path) {
    args <- prepare_fixture(out_path)
    args$dataset_name <- "Caf\u00e9 \u03b1tlas"
    rownames(args$gene_expression) <- c("z\u00e9ro", "\u03b1", "\u03b1", "c4", "c5", "c6")
    args$var$symbol <- c("CD3\u03b5", "CD8\u03b1", "MS4A1", "CD3\u03b5")
    args$var_gene_id_column <- "symbol"
    args$obs$kind <- c("T \u03b1", "B \u00e9", "T \u03b1", "NK \u2192", "B \u00e9", "T \u03b1")
    args
  }
  here <- messages_of(do.call(cellucid_prepare, arguments(store_path())))
  there <- with_non_utf8_session(messages_of(do.call(cellucid_prepare, arguments(store_path()))))
  expect_identical(here$messages, c(
    "Gene names from var column 'symbol': 1 duplicate name was made unique, e.g. 'CD3\u03b5' \u2192 'CD3\u03b5-1'.\n",
    "Cell names: 1 duplicate name was made unique, e.g. '\u03b1' \u2192 '\u03b1-1'.\n"
  ))
  expect_identical(there$messages, here$messages)
  files <- sort(list.files(here$value, recursive = TRUE, all.files = TRUE))
  expect_identical(sort(list.files(there$value, recursive = TRUE, all.files = TRUE)), files)
  expect_identical(
    unname(tools::md5sum(file.path(there$value, files))),
    unname(tools::md5sum(file.path(here$value, files)))
  )
  expect_identical(store_array(here$value, "obs/_index")[1:3], c("z\u00e9ro", "\u03b1", "\u03b1-1"))
  expect_identical(store_array(here$value, "obs/kind/categories"), c("B \u00e9", "NK \u2192", "T \u03b1"))
})

test_that("categories are in code-point order whatever each string's encoding", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  # One latin1 'e-acute' (the byte 0xE9) among UTF-8 strings: a byte-wise sort
  # put it after the Greek alpha (0xCE 0xB1), though U+00E9 comes first.
  args$obs$kind <- c(latin1("\u00e9"), "\u03b1", "a", latin1("\u00e9"), "\u03b1", "a")
  store <- suppressMessages(do.call(cellucid_prepare, args))
  expect_identical(store_array(store, "obs/kind/categories"), c("a", "\u00e9", "\u03b1"))
  expect_identical(store_array(store, "obs/kind/codes"), c(1L, 2L, 0L, 1L, 2L, 0L))
})

test_that("repeated names are made unique by code point in a session that is not UTF-8", {
  fixture <- cellucid_read_fixture("data_mapping_conformance.json")
  with_non_utf8_session({
    # A name that is the escape R writes for another name stays apart from it.
    expect_no_warning(result <- cellucid:::.make_unique_names(c("\u03b1", "<U+03B1>", "\u03b1")))
    for (case in fixture$makeUniqueNames) {
      if (!is.null(case$expected$error)) next
      index <- if (is.null(case$index)) NULL else vapply(case$index, function(value) {
        if (is.null(value)) NA_character_ else value
      }, character(1))
      expect_no_warning(shared <- cellucid:::.make_unique_names(case$names, index = index, join = case$join))
      expect_identical(shared$names, as.character(unlist(case$expected$names)), info = case$name)
    }
  })
  expect_identical(result$names, c("\u03b1", "<U+03B1>", "\u03b1-1"))
  expect_identical(result$renamed, list(c("\u03b1", "\u03b1-1")))
})

test_that("JSON and frames from the browser are read as UTF-8 in a session that is not UTF-8", {
  document <- tempfile(fileext = ".json")
  writeBin(cellucid:::.utf8_bytes("{\"name\":\"Caf\u00e9 \u03b1\"}"), document)
  header <- cellucid:::.utf8_bytes("{\"type\":\"selection\",\"source\":\"lasso-\u03b1\"}")
  binary <- c(writeBin(length(header), raw(), size = 4L, endian = "little"), header, writeBin(c(1L, 3L), raw(), size = 4L, endian = "little"))
  with_non_utf8_session({
    expect_identical(cellucid:::.read_json_file(document)$name, "Caf\u00e9 \u03b1")
    # httpuv hands a text frame's UTF-8 bytes over as native text.
    text <- rawToChar(cellucid:::.utf8_bytes("{\"type\":\"hover\",\"cell\":\"\u03b1\",\"position\":null}"))
    expect_identical(cellucid:::.parse_json_frame(text, "An event")$cell, "\u03b1")
    event <- cellucid:::.parse_viewer_binary_event(binary)
    expect_identical(event$payload$source, "lasso-\u03b1")
    expect_identical(event$payload$cells, c(1L, 3L))
  })
  expect_identical(cellucid:::.utf8_text(as.raw(c(0x61, 0x00))), NA_character_)
  expect_identical(cellucid:::.utf8_text(as.raw(0xE9)), NA_character_)
  expect_identical(cellucid:::.as_utf8(rawToChar(as.raw(0xE9))), NA_character_)
})

test_that("a store folder's consolidated metadata is the same bytes in a session that is not UTF-8", {
  skip_if_not_installed("Matrix")
  folder <- suppressMessages(do.call(cellucid_prepare, prepare_fixture(store_path())))
  unlink(file.path(folder, ".zmetadata"))
  attrs <- file.path(folder, "uns", ".zattrs")
  writeBin(cellucid:::.utf8_bytes("{\"note\":\"Caf\u00e9 \u03b1 \u2192\"}"), attrs)
  here <- cellucid:::.consolidated_zarr_metadata(folder, 2L)
  there <- with_non_utf8_session(cellucid:::.consolidated_zarr_metadata(folder, 2L))
  expect_identical(there, here)
  expect_identical(
    jsonlite::parse_json(cellucid:::.utf8_text(here$bytes))$metadata[["uns/.zattrs"]]$note,
    "Caf\u00e9 \u03b1 \u2192"
  )
})

test_that("a memory store and its cache find a key the browser asks for in a session that is not UTF-8", {
  with_non_utf8_session({
    store <- cellucid:::.store_open_memory(2L)
    expect_no_warning(cellucid:::.store_put(store, "obs/\u03b1/.zattrs", as.raw(1:3)))
    # The request path, percent-decoded as the server decodes it.
    asked <- cellucid:::.as_utf8(utils::URLdecode("obs/%CE%B1/.zattrs"))
    expect_no_warning(expect_identical(cellucid:::.store_get(store, asked), as.raw(1:3)))
    expect_null(cellucid:::.store_get(store, "obs/<U+03B1>/.zattrs"))
    expect_identical(store$keys, "obs/\u03b1/.zattrs")
    made <- 0L
    make <- function() {
      made <<- made + 1L
      as.raw(1:4)
    }
    cache <- cellucid:::.chunk_cache(4)
    expect_no_warning(cellucid:::.cached(cache, "obs/\u03b1/0", make))
    expect_identical(cellucid:::.cached(cache, cellucid:::.as_utf8(utils::URLdecode("obs/%CE%B1/0")), make), as.raw(1:4))
    expect_identical(made, 1L)
    # Over budget, the oldest key is dropped by its own name.
    cellucid:::.cached(cache, "obs/\u03b2/0", make)
    expect_identical(made, 2L)
    cellucid:::.cached(cache, "obs/\u03b1/0", make)
    expect_identical(made, 3L)
  })
})

test_that("an object's column named in Greek is served from memory in a session that is not UTF-8", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  SummarizedExperiment::colData(object)[["\u03b1"]] <- factor(rep(c("x", "y"), length.out = ncol(object)))
  with_non_utf8_session({
    viewer <- test_viewer(object, dataset_id = "greek", dataset_name = "Greek")
    url <- paste0(viewer$url, "/_cellucid/anndata/", viewer$served$route_name, "/obs/%CE%B1/.zattrs")
    response <- http_request(url)
    expect_identical(response$status, 200L)
    expect_identical(jsonlite::parse_json(cellucid:::.utf8_text(response$body))[["encoding-type"]], "categorical")
  })
})

test_that("a file httpuv streams is named in UTF-8, whatever the name's encoding", {
  # httpuv opens a body file by its name's bytes read as UTF-8.
  path <- latin1(file.path(tempfile("cellucid_"), "caf\u00e9.bin"))
  expect_identical(Encoding(path), "latin1")
  body <- cellucid:::.file_body(path, owned = TRUE)
  expect_identical(charToRaw(body$file), charToRaw(enc2utf8(path)))
  expect_true(body$owned)
  # A latin1 session on a UTF-8 file system (macOS) cannot name the file.
  dir.create(dirname(path))
  made <- tryCatch({
    writeBin(as.raw(1:10), path)
    TRUE
  }, error = function(error) FALSE, warning = function(warning) FALSE)
  skip_if_not(made, "this session cannot name a file with an e-acute")
  response <- cellucid:::.file_response(list(), path, FALSE)
  expect_identical(response$status, 200L)
  expect_true(validUTF8(response$body$file))
  expect_identical(charToRaw(response$body$file), charToRaw(enc2utf8(path)))
  expect_false(response$body$owned)
})

test_that("code points are read from a string whatever its encoding", {
  expect_identical(cellucid:::.code_points(latin1("caf\u00e9")), c(99L, 97L, 102L, 233L))
  expect_identical(cellucid:::.compare_code_points(latin1("\u00e9"), "\u03b1"), -1L)
  expect_identical(cellucid:::.compare_code_points("\u03b1", latin1("\u00e9")), 1L)
  expect_error(
    cellucid:::.validate_display_text(latin1("Caf\u00e9 "), "dataset_name"),
    "glyph: 'Caf\u00e9 ' ends with U+0020 SPACE.",
    fixed = TRUE
  )
  expect_identical(cellucid:::.sort_code_points(c("\u03b1", latin1("\u00e9"), "a")), c("a", "\u00e9", "\u03b1"))
})

test_that("a stream is inflated whole however long its output, and a declared size R cannot hold is refused", {
  member <- memCompress(raw(1e6), "gzip")
  # Unknown size: the output grows from four times the input.
  expect_identical(cellucid:::.native_inflate(member, FALSE), raw(1e6))
  expect_identical(cellucid:::.native_inflate(member, FALSE, NA), raw(1e6))
  expect_identical(cellucid:::.native_inflate(member, FALSE, 1e6), raw(1e6))
  expect_error(cellucid:::.native_inflate(member, FALSE, 10.5), "whole number of bytes")
  expect_error(cellucid:::.native_inflate(member, FALSE, 2^53), "more than one R raw vector holds")
})

test_that("a chunk longer than one R raw vector is refused before its byte count can wrap", {
  # A length R_xlen_t holds whose bytes it does not: 2^30 four-byte values on
  # 32-bit R (an int R_xlen_t), 2^51 on 64-bit (R_XLEN_T_MAX is 2^52).
  n_values <- if (.Machine$sizeof.pointer == 4L) 2^30 else 2^51
  expect_error(cellucid:::.native_encode_values(1, 1, 1, n_values, "<f4"), "more than one R raw vector holds")
  expect_error(cellucid:::.native_encode_values(1, 1, 1, 2^60, "<f4"), "whole number from 0 to")
  expect_error(cellucid:::.native_raw_slice(as.raw(1:4), 0, 2, -1), "whole number from 0 to")
})
