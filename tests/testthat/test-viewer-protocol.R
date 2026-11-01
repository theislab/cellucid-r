# The viewer channel's declarations, against the shared protocol fixture.
#
# tests/testthat/fixtures/viewer_protocol.json is a byte-identical copy of the
# web app's tests/fixtures/viewer-protocol.json (written by
# assets/js/data/viewer-protocol.js, the source of truth) and of
# cellucid-python's tests/fixtures/viewer_protocol.json; all three
# repositories pin the same SHA-256. R receives events and sends commands, so
# it must accept every example event and binary event, and produce every
# example command.

.viewer_protocol_sha256 <- "24b8606e4ce533f5320b4adb88b6aaac314091943bbb5f3c6590f55b7adfd32a"

protocol_fixture <- function() cellucid_read_fixture("viewer_protocol.json")

hex_raw <- function(hex) {
  as.raw(strtoi(substring(hex, seq(1, nchar(hex), 2), seq(2, nchar(hex), 2)), 16L))
}

test_that("the three repositories hold the same pinned viewer protocol", {
  expect_pinned_fixture("viewer_protocol.json", .viewer_protocol_sha256)
  expect_identical(protocol_fixture()$schema, 2L)
})

test_that("R declares exactly the fixture's channel, commands and events", {
  fixture <- protocol_fixture()
  strings <- function(x) lapply(x, function(keys) as.character(unlist(keys)))
  expect_identical(cellucid:::.viewer_channel_path, fixture$channel$path)
  expect_identical(as.list(cellucid:::.viewer_close_codes), fixture$channel$closeCodes)
  expect_identical(cellucid:::.viewer_command_fields, strings(fixture$commands))
  expect_identical(cellucid:::.viewer_binary_command_fields, strings(fixture$binaryCommands))
  expect_identical(cellucid:::.viewer_event_fields, strings(fixture$events))
  expect_identical(cellucid:::.viewer_binary_event_fields, strings(fixture$binaryEvents))
  expect_identical(cellucid:::.viewer_data_mapping_origins, as.character(unlist(fixture$dataMappingOrigins)))
  expect_identical(cellucid:::.viewer_result_meta_fields, lapply(fixture$results$kinds, function(keys) as.character(unlist(keys))))
  expect_identical(cellucid:::.viewer_highlight_tools, as.character(unlist(fixture$results$highlightTools)))
  expect_identical(cellucid:::.viewer_result_origins, as.character(unlist(fixture$results$origins)))
  expect_identical(cellucid:::.viewer_field_kinds, as.character(unlist(fixture$results$fieldKinds)))
  expect_identical(as.list(cellucid:::.viewer_value_types), fixture$results$valueTypes)
  expect_identical(cellucid:::.viewer_result_epoch_pattern, fixture$results$epochPattern)
  expect_gte(cellucid:::.viewer_reason_limit, fixture$commandErrorReasonLimit)
  document <- jsonlite::fromJSON(cellucid:::.json_text(cellucid:::.viewer_protocol_document()), simplifyVector = FALSE)
  expect_identical(document, fixture$protocolDocument)
})

test_that("every example event of the fixture is accepted", {
  fixture <- protocol_fixture()
  for (event in fixture$examples$events) {
    text <- as.character(jsonlite::toJSON(event, auto_unbox = TRUE, null = "null", digits = NA))
    parsed <- cellucid:::.parse_viewer_event(text)
    expect_identical(parsed$type, event$type)
    expect_setequal(names(parsed$payload), setdiff(names(event), "type"))
  }
  hello <- fixture$examples$hello
  parsed <- cellucid:::.parse_viewer_hello(as.character(jsonlite::toJSON(hello, auto_unbox = TRUE)))
  expect_identical(parsed$viewerToken, hello$viewerToken)
})

test_that("every example binary event of the fixture is decoded", {
  for (example in protocol_fixture()$examples$binaryEvents) {
    parsed <- cellucid:::.parse_viewer_binary_event(hex_raw(example$hex))
    expect_identical(parsed$type, example$header$type)
    if (identical(parsed$type, "selection")) {
      expect_equal(parsed$payload$cells, as.numeric(unlist(example$cells)))
      expect_identical(parsed$payload$source, example$header$source)
    } else if (identical(parsed$type, "result")) {
      expected_cells <- if (is.null(example$cells)) NULL else as.numeric(unlist(example$cells))
      expect_equal(parsed$payload$cells, expected_cells, info = example$header$kind)
      expect_equal(parsed$payload$values, if (is.null(example$values)) NULL else as.numeric(unlist(example$values)), info = example$header$kind)
      expect_identical(parsed$payload$kind, example$header$kind)
    } else {
      expect_identical(rawToChar(parsed$body), example$bodyText)
      expect_identical(parsed$payload$requestId, example$header$requestId)
    }
  }
})

test_that("every example command of the fixture is produced as the same frame", {
  for (command in protocol_fixture()$examples$commands) {
    frame <- cellucid:::.viewer_command_frame(command)
    expect_true(is.character(frame), info = command$type)
    expect_identical(jsonlite::fromJSON(frame, simplifyVector = FALSE), command, info = command$type)
  }
  # A command that names cells is the example's binary frame, byte for byte.
  types <- character()
  for (example in protocol_fixture()$examples$binaryCommands) {
    header <- example$header
    if (identical(header$type, "setField")) {
      # A pushed column: its cells (NULL: every cell), its values and their type.
      given <- header[setdiff(names(header), c("cells", "values"))]
      given["cells"] <- list(if (!is.null(example$cells)) as.numeric(unlist(example$cells)))
      given$values <- as.numeric(unlist(example$values))
      given$valueType <- header$values
      frame <- cellucid:::.viewer_command_frame(given)
      expect_identical(frame, hex_raw(example$hex), info = example$hex)
      types <- c(types, header$type)
      next
    }
    every <- isTRUE(header$all)
    given <- header[setdiff(names(header), "all")]
    given["cells"] <- list(if (!every) as.numeric(unlist(example$cells)))
    frame <- cellucid:::.viewer_command_frame(given)
    expect_identical(frame, hex_raw(example$hex), info = example$hex)
    types <- c(types, header$type)
  }
  expect_setequal(types, names(cellucid:::.viewer_binary_command_fields))
})

test_that("frames that break the protocol are refused with the reason", {
  refused <- function(text) {
    expect_error(cellucid:::.parse_viewer_event(text), class = "cellucid_protocol_error")
  }
  refused("not json")
  refused("[1, 2]")
  refused("{\"type\":\"mystery\"}")
  refused("{\"type\":\"ready\",\"n_cells\":10,\"dimensions\":2}")
  refused("{\"type\":\"ready\",\"n_cells\":10,\"n_obs\":10,\"dimensions\":4}")
  refused("{\"type\":\"ready\",\"n_cells\":-1,\"n_obs\":10,\"dimensions\":2}")
  refused("{\"type\":\"ready\",\"n_cells\":11,\"n_obs\":10,\"dimensions\":2}")
  refused("{\"type\":\"color_by\",\"field\":\"x\",\"source\":null,\"symbol\":null}")
  refused("{\"type\":\"color_by\",\"field\":\"x\",\"source\":\"obs\",\"symbol\":\"X\"}")
  refused("{\"type\":\"color_by\",\"field\":\"x\",\"source\":\"var\",\"label\":null}")
  refused("{\"type\":\"result_removed\",\"epoch\":\"0123456789abcdef\",\"revision\":0,\"kind\":\"highlight\",\"id\":\"h\"}")
  refused("{\"type\":\"result_removed\",\"epoch\":\"0123456789ABCDEF\",\"revision\":1,\"kind\":\"highlight\",\"id\":\"h\"}")
  refused("{\"type\":\"results\",\"requestId\":null,\"revision\":1,\"nObs\":3,\"items\":[{\"kind\":\"nope\",\"id\":\"a\",\"revision\":1}]}")
  refused("{\"type\":\"results\",\"epoch\":\"0123456789abcdef\",\"dataset\":null,\"revision\":1,\"nObs\":3,\"items\":[{\"kind\":\"highlight\",\"id\":\"a\",\"revision\":1}]}")
  refused("{\"type\":\"click\",\"cell\":1,\"button\":3,\"shift\":false,\"ctrl\":false}")
  refused("{\"type\":\"hover\",\"cell\":1,\"position\":{\"x\":1,\"y\":2}}")
  refused("{\"type\":\"command_error\",\"command\":\"nope\",\"reason\":\"x\"}")
  refused(paste0("{\"type\":\"command_error\",\"command\":\"ping\",\"reason\":\"", strrep("x", 501), "\"}"))
  refused("{\"type\":\"data_mapping\",\"requestId\":\"r1\",\"origin\":\"user\",\"mapping\":{\"kind\":\"anndata\"}}")
  refused("{\"type\":\"data_mapping\",\"requestId\":null,\"origin\":\"notebook\",\"mapping\":{\"kind\":\"anndata\"}}")
  refused("{\"type\":\"pong\",\"requestId\":\"r1\",\"t\":1,\"viewerId\":\"v\"}")
  expect_error(
    cellucid:::.parse_viewer_hello("{\"type\":\"hello\",\"viewerId\":\"a\"}"),
    class = "cellucid_protocol_error"
  )

  header <- charToRaw("{\"type\":\"selection\",\"source\":\"lasso\"}")
  frame <- function(cells) c(cellucid:::.encode_uint32(length(header)), header, cellucid:::.encode_uint32(cells))
  expect_error(cellucid:::.parse_viewer_binary_event(frame(c(3, 2))), "strictly ascending", class = "cellucid_protocol_error")
  expect_error(cellucid:::.parse_viewer_binary_event(frame(c(1, 1))), "strictly ascending", class = "cellucid_protocol_error")
  expect_error(cellucid:::.parse_viewer_binary_event(frame(numeric(0))), "one or more", class = "cellucid_protocol_error")
  expect_error(
    cellucid:::.parse_viewer_binary_event(c(frame(1), as.raw(0))),
    "one or more", class = "cellucid_protocol_error"
  )
  expect_error(cellucid:::.parse_viewer_binary_event(as.raw(c(200, 0, 0, 0, 1))), "outside the frame", class = "cellucid_protocol_error")
})

test_that("a selection of millions of cells and of indices past 2^31 decodes exactly", {
  cells <- c(seq(0, 2e6 - 1, by = 1), 2147483647, 2147483648, 4294967295)
  header <- charToRaw("{\"type\":\"selection\",\"source\":\"lasso\"}")
  frame <- c(cellucid:::.encode_uint32(length(header)), header, cellucid:::.encode_uint32(cells))
  elapsed <- system.time(parsed <- cellucid:::.parse_viewer_binary_event(frame))[["elapsed"]]
  expect_identical(parsed$payload$cells, cells)
  # Vectorised: two million indices decode in well under a second.
  expect_lt(elapsed, 2)
})

test_that("commands with wrong values are refused before they are sent", {
  frame <- cellucid:::.viewer_command_frame
  expect_error(frame(list(type = "highlight", cells = numeric(0), color = "#ff0000", label = NULL)), "at least one")
  expect_error(frame(list(type = "highlight", cells = c(1, 1), color = "#ff0000", label = NULL)), "strictly ascending")
  expect_error(frame(list(type = "highlight", cells = c(2, 1), color = "#ff0000", label = NULL)), "strictly ascending")
  expect_error(frame(list(type = "highlight", cells = c(1, 2^32), color = "#ff0000", label = NULL)), "strictly ascending")
  expect_error(frame(list(type = "highlight", cells = c(-1, 2), color = "#ff0000", label = NULL)), "strictly ascending")
  expect_error(frame(list(type = "highlight", cells = c(1.5, 2), color = "#ff0000", label = NULL)), "strictly ascending")
  expect_error(frame(list(type = "highlight", cells = 1, color = "red", label = NULL)), "six-digit")
  expect_error(frame(list(type = "setVisibility", cells = NULL, visible = NA)), "TRUE or FALSE")
  expect_error(frame(list(type = "highlight", color = "#ff0000", label = NULL)), "carries exactly type, cells, color, label")
  expect_error(frame(list(type = "ping")), "carries exactly type, requestId")
  expect_error(frame(list(type = "explode")), "Unknown viewer command")
  expect_error(frame(list(type = "setDataMapping", requestId = "r", mapping = list(kind = "seurat"))), "data mapping")
})

test_that("a highlight of millions of cells is one binary frame, made in one pass", {
  cells <- seq(0, 1e7 - 1, by = 1)
  elapsed <- system.time(frame <- cellucid:::.viewer_command_frame(list(type = "highlight", cells = cells, color = "#00ff00", label = NULL)))[["elapsed"]]
  header_length <- readBin(frame[1:4], "integer", size = 4L, endian = "little")
  expect_identical(length(frame), 4L + header_length + 40000000L)
  body <- frame[(5L + header_length):length(frame)]
  expect_identical(readBin(body, "integer", n = 1e7, size = 4L, endian = "little"), as.integer(cells))
  expect_lt(elapsed, 2)
})
