# The viewer channel's declarations, against the shared protocol fixture.
#
# tests/testthat/fixtures/viewer_protocol.json is a byte-identical copy of the
# web app's tests/fixtures/viewer-protocol.json (written by
# assets/js/data/viewer-protocol.js, the source of truth) and of
# cellucid-python's tests/fixtures/viewer_protocol.json; all three
# repositories pin the same SHA-256. R receives events and sends commands, so
# it must accept every example event and binary event, and produce every
# example command -- a binary action byte for byte.

.viewer_protocol_sha256 <- "c117822a5d87df03439fdbed9af460617764b59d9e560ba7c0e519ec5f708dae"

protocol_fixture <- function() cellucid_read_fixture("viewer_protocol.json")

hex_raw <- function(hex) {
  as.raw(strtoi(substring(hex, seq(1, nchar(hex), 2), seq(2, nchar(hex), 2)), 16L))
}

raw_hex <- function(bytes) paste(as.character(bytes), collapse = "")

test_that("the three repositories hold the same pinned viewer protocol", {
  expect_pinned_fixture("viewer_protocol.json", .viewer_protocol_sha256)
  expect_identical(protocol_fixture()$schema, 5L)
  for (copy in list(c("cellucid", "tests", "fixtures", "viewer-protocol.json"), c("cellucid-python", "tests", "fixtures", "viewer_protocol.json"))) {
    path <- do.call(cellucid_sibling_file, as.list(copy))
    expect_identical(cellucid_file_sha256(path), .viewer_protocol_sha256, info = path)
  }
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
  expect_identical(cellucid:::.viewer_action_type_pattern, fixture$action$typePattern)
  expect_identical(cellucid:::.viewer_action_array_fields, as.character(unlist(fixture$action$arrayFields)))
  expect_identical(cellucid:::.viewer_action_step_fields, as.character(unlist(fixture$action$stepFields)))
  expect_identical(as.list(cellucid:::.viewer_action_unrecorded), fixture$action$unrecorded)
  expect_identical(cellucid:::.viewer_history_change_kinds, as.character(unlist(fixture$history$changeKinds)))
  expect_identical(cellucid:::.viewer_session_loaded_tracker_fields, as.character(unlist(fixture$history$sessionLoadedTrackerFields)))
  expect_identical(cellucid:::.viewer_data_mapping_origins, as.character(unlist(fixture$dataMappingOrigins)))
  expect_identical(cellucid:::.viewer_result_meta_fields, lapply(fixture$results$kinds, function(keys) as.character(unlist(keys))))
  expect_identical(cellucid:::.viewer_highlight_tools, as.character(unlist(fixture$results$highlightTools)))
  expect_identical(cellucid:::.viewer_result_origins, as.character(unlist(fixture$results$origins)))
  expect_identical(as.list(cellucid:::.viewer_value_types), fixture$results$valueTypes)
  expect_identical(cellucid:::.viewer_result_epoch_pattern, fixture$results$epochPattern)
  expect_gte(cellucid:::.viewer_reason_limit, fixture$commandErrorReasonLimit)
  # An agent command's connection and its id; names of any length.
  expect_identical(cellucid:::.viewer_agent_connection_fields, as.character(unlist(fixture$agent$connectionFields)))
  expect_identical(cellucid:::.viewer_agent_connection_id_pattern, fixture$agent$connectionIdPattern)
  expect_null(fixture$agent$nameLimit)
  expect_setequal(cellucid:::.viewer_agent_events, grep("^agent_", c(names(fixture$events), names(fixture$binaryEvents)), value = TRUE))
  # GET /_cellucid/protocol: every command R sends, the binary ones, and
  # every event R accepts (an action's answer included).
  document <- jsonlite::fromJSON(cellucid:::.json_text(cellucid:::.viewer_protocol_document()), simplifyVector = FALSE)
  expect_identical(document, fixture$protocolDocument)
  expect_identical(unlist(document$commands), c(
    "action", "agent_answer", "agent_call", "agent_cancel", "agent_list", "debug_snapshot", "freeze", "load_session",
    "ping", "requestSessionBundle"
  ))
  expect_identical(unlist(document$binaryCommands), c("action", "load_session"))
  expect_true(all(c("action_result", "history_changed", "session_loaded", "agent_tools", "agent_result", "agent_asking",
                    "agent_answered") %in% unlist(document$events)))
})

test_that("a handler takes every event of the fixture but an agent app's, and a misspelt one is refused", {
  # cellucid-python's hooks take the same names (cellucid.jupyter.HOOK_EVENTS),
  # pinned to the same fixture. A misspelt event used to be partially matched
  # (match.arg()) or refused with R's generic words.
  fixture <- protocol_fixture()
  events <- c(names(fixture$events), names(fixture$binaryEvents))
  expect_setequal(cellucid:::.hook_events(), c(grep("^agent_", events, value = TRUE, invert = TRUE), "message"))
  viewer <- new.env()
  class(viewer) <- "cellucid_viewer"
  for (event in list("selections", "sel", "agent_result", 1, NA_character_, c("click", "hover"))) {
    expect_error(cellucid_on(viewer, event, identity), "^event must be one of \"action_result\", .*\"message\"; got ")
    expect_error(cellucid_off(viewer, event), "^event must be one of ")
  }
  expect_error(cellucid_on(viewer, "selections", identity), "; got c\\(\"selections\"\\)\\.$")
  expect_error(cellucid_on(viewer, 1, identity), "; got a numeric\\.$")
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
  types <- character()
  for (example in protocol_fixture()$examples$binaryEvents) {
    parsed <- cellucid:::.parse_viewer_binary_event(hex_raw(example$hex))
    expect_identical(parsed$type, example$header$type)
    types <- c(types, parsed$type)
    if (identical(parsed$type, "selection")) {
      expect_equal(parsed$payload$cells, as.numeric(unlist(example$cells)))
      expect_identical(parsed$payload$source, example$header$source)
    } else if (identical(parsed$type, "result")) {
      expected_cells <- if (is.null(example$cells)) NULL else as.numeric(unlist(example$cells))
      expect_equal(parsed$payload$cells, expected_cells, info = example$header$kind)
      expect_equal(parsed$payload$values, if (is.null(example$values)) NULL else as.numeric(unlist(example$values)), info = example$header$kind)
      expect_identical(parsed$payload$kind, example$header$kind)
    } else if (identical(parsed$type, "action_result")) {
      # The outcome is the body's UTF-8 JSON, as jsonlite reads JSON.
      expect_identical(parsed$payload$requestId, example$header$requestId)
      expect_identical(parsed$payload$step, example$header$step)
      expect_identical(parsed$payload$unrecorded, example$header$unrecorded)
      expect_identical(parsed$payload$outcome, example$outcome)
    } else if (parsed$type %in% c("agent_tools", "agent_result")) {
      # The body is UTF-8 JSON: the tools, or the call's MCP tool result.
      key <- if (identical(parsed$type, "agent_tools")) "tools" else "result"
      expect_identical(parsed$payload$requestId, example$header$requestId)
      expect_identical(parsed$payload[[key]], example$body)
    } else {
      expect_identical(rawToChar(parsed$body), example$bodyText)
      expect_identical(parsed$payload$requestId, example$header$requestId)
    }
  }
  expect_setequal(types, names(cellucid:::.viewer_binary_event_fields))
})

test_that("every example command is produced, a binary action byte for byte", {
  for (command in protocol_fixture()$examples$commands) {
    frame <- if (identical(command$type, "action")) {
      cellucid:::.encode_viewer_action(command$requestId, command$action, command$params, command$turn, command$turnLabel)
    } else {
      cellucid:::.viewer_command_frame(command)
    }
    expect_true(is.character(frame), info = command$type)
    expect_identical(jsonlite::fromJSON(frame, simplifyVector = FALSE), command, info = frame)
  }
  # A binary action: the header's params with each array placed at its path
  # (its key appended, as JavaScript appends one), encoded as the web app's
  # encodeViewerAction encodes it -- compact JSON in key order, padded with
  # spaces to 4 bytes, then each array little-endian.
  examples <- protocol_fixture()$examples$binaryCommands
  expect_length(examples, 7L)
  for (example in examples) {
    if (identical(example$header$type, "load_session")) {
      # A session file pushed in: the header, then the file's bytes.
      frame <- cellucid:::.encode_viewer_load_session(example$header$requestId, example$header$name, charToRaw(example$bodyText))
      expect_identical(raw_hex(frame), example$hex)
      next
    }
    params <- example$header$params
    for (array in example$arrays) {
      params <- place_at(params, array$path, cellucid:::.typed_array(as.numeric(unlist(array$elements)), array$type))
    }
    frame <- cellucid:::.encode_viewer_action(example$header$requestId, example$header$action, params, example$header$turn,
                                              example$header$turnLabel)
    expect_identical(raw_hex(frame), example$hex, info = example$header$requestId)
    header_length <- readBin(frame[1:4], "integer", size = 4L, endian = "little")
    expect_identical(header_length %% 4L, 0L)
  }
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
  # A command_error names its command, the request it refused (or null) and why.
  refused("{\"type\":\"command_error\",\"command\":\"nope\",\"requestId\":null,\"reason\":\"x\"}")
  refused("{\"type\":\"command_error\",\"command\":\"action\",\"reason\":\"x\"}")
  refused("{\"type\":\"command_error\",\"command\":\"action\",\"requestId\":\"\",\"reason\":\"x\"}")
  refused("{\"type\":\"command_error\",\"command\":\"action\",\"requestId\":7,\"reason\":\"x\"}")
  refused(paste0("{\"type\":\"command_error\",\"command\":\"ping\",\"requestId\":null,\"reason\":\"",
                 strrep("x", cellucid:::.viewer_reason_limit + 1L), "\"}"))
  # The commands schema 2 had are no longer the protocol's.
  refused("{\"type\":\"command_error\",\"command\":\"setColorBy\",\"requestId\":null,\"reason\":\"x\"}")
  refused("{\"type\":\"data_mapping\",\"requestId\":\"r1\",\"origin\":\"notebook\",\"mapping\":{\"kind\":\"anndata\"}}")
  refused("{\"type\":\"data_mapping\",\"origin\":\"nowhere\",\"mapping\":{\"kind\":\"anndata\"}}")
  refused("{\"type\":\"pong\",\"requestId\":\"r1\",\"t\":1,\"viewerId\":\"v\"}")
  # A history change names a kind of the protocol, the tracker and the step shown.
  change <- "{\"type\":\"history_changed\",\"kind\":\"%s\",\"tracker\":%s,\"head\":\"s2\",\"step\":null,\"branch\":%s}"
  expect_identical(cellucid:::.parse_viewer_event(sprintf(change, "moved", "\"t-1\"", "null"))$payload$branch, NULL)
  refused(sprintf(change, "renamed", "\"t-1\"", "null"))
  refused(sprintf(change, "moved", "null", "null"))
  refused(sprintf(change, "moved", "\"t-1\"", "7"))
  loaded <- "{\"type\":\"session_loaded\",\"requestId\":\"l1\",\"tracker\":%s,\"head\":\"s9\"}"
  expect_identical(cellucid:::.parse_viewer_event(sprintf(loaded, "{\"id\":\"t-1\",\"name\":\"T cell analysis\"}"))$payload$head, "s9")
  refused(sprintf(loaded, "{\"id\":\"t-1\"}"))
  refused(sprintf(loaded, "\"t-1\""))
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

  # An action's answer: exact header keys, a step of exactly {id, label} or
  # null, a body of one UTF-8 JSON value.
  answer <- function(header, body) {
    bytes <- charToRaw(header)
    c(cellucid:::.encode_uint32(length(bytes)), bytes, body)
  }
  ok <- answer("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":null,\"unrecorded\":null}", charToRaw("null"))
  expect_null(cellucid:::.parse_viewer_binary_event(ok)$payload$outcome)
  # Why a change made no step: one of the closed set, never beside a step.
  for (reason in names(cellucid:::.viewer_action_unrecorded)) {
    said <- answer(paste0("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":null,\"unrecorded\":\"", reason, "\"}"), charToRaw("null"))
    expect_identical(cellucid:::.parse_viewer_binary_event(said)$payload$unrecorded, reason)
  }
  refused_binary <- function(frame, pattern) {
    expect_error(cellucid:::.parse_viewer_binary_event(frame), pattern, class = "cellucid_protocol_error")
  }
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":null,\"unrecorded\":null}", raw(0)), "one JSON value")
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":null,\"unrecorded\":null}", charToRaw("{nope")), "one JSON value")
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":null,\"unrecorded\":null}", as.raw(c(0x22, 0xff, 0x22))), "UTF-8")
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":{\"id\":\"s1\"},\"unrecorded\":null}", charToRaw("null")), "exactly")
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"\",\"step\":null,\"unrecorded\":null}", charToRaw("null")), "requestId")
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"a1\"}", charToRaw("null")), "exactly")
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":null}", charToRaw("null")), "exactly")
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":null,\"unrecorded\":\"forgotten\"}", charToRaw("null")), "unrecorded must be one of unchanged, retracted")
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":null,\"unrecorded\":3}", charToRaw("null")), "unrecorded must be one of")
  refused_binary(answer("{\"type\":\"action_result\",\"requestId\":\"a1\",\"step\":{\"id\":\"s1\",\"label\":\"Colored\"},\"unrecorded\":\"unchanged\"}", charToRaw("null")), "unrecorded must be null when the action made a step")
  refused_binary(answer("{\"type\":\"result\",\"requestId\":null,\"epoch\":\"0123456789abcdef\",\"revision\":1,\"kind\":\"visibility\",\"id\":\"visible\",\"nObs\":1,\"meta\":{},\"cells\":null,\"values\":null}", raw(0)), "exactly")
})

test_that("a categorical result carries the colors it is shown in, one a category", {
  # The colors of the view in focus (RCC-1): one lowercase #rrggbb a
  # category, in the categories' order; a pull writes them with the column.
  categorical <- function(colors) {
    header <- list(type = "result", epoch = "0123456789abcdef", revision = 1, kind = "categorical", id = "user_cat_1", nObs = 3,
                   meta = list(key = "stage", categories = list("early", "late"), colors = colors, origin = "user"),
                   cells = NULL, values = "uint8")
    bytes <- charToRaw(as.character(jsonlite::toJSON(header, auto_unbox = TRUE, null = "null", digits = NA)))
    c(cellucid:::.encode_uint32(length(bytes)), bytes, as.raw(c(0, 1, 255)))
  }
  decoded <- cellucid:::.parse_viewer_binary_event(categorical(list("#808080", "#008080")))$payload
  expect_identical(unlist(decoded$meta$colors), c("#808080", "#008080"))
  expect_identical(as.double(decoded$values), c(0, 1, 255))
  for (colors in list(list("#808080"), list("#808080", "#008080", "#000000"), list("#808080", "#00808G"),
                      list("#808080", "#00FF00"), list("#808080", "#0f0"), list("#808080", 7), list())) {
    expect_error(cellucid:::.parse_viewer_binary_event(categorical(colors)), "one #rrggbb color (lowercase) a category",
                 fixed = TRUE, class = "cellucid_protocol_error")
  }
  # A categorical without its colors is not the protocol's.
  header <- list(type = "result", epoch = "0123456789abcdef", revision = 1, kind = "categorical", id = "user_cat_1", nObs = 3,
                 meta = list(key = "stage", categories = list("early", "late"), origin = "user"), cells = NULL, values = "uint8")
  bytes <- charToRaw(as.character(jsonlite::toJSON(header, auto_unbox = TRUE, null = "null", digits = NA)))
  expect_error(cellucid:::.parse_viewer_binary_event(c(cellucid:::.encode_uint32(length(bytes)), bytes, as.raw(c(0, 1, 255)))),
               class = "cellucid_protocol_error")
})

test_that("a selection of millions of cells and of indices past 2^31 decodes exactly", {
  cells <- c(seq(0, 2e6 - 1, by = 1), 2147483647, 2147483648, 4294967295)
  header <- charToRaw("{\"type\":\"selection\",\"source\":\"lasso\"}")
  frame <- c(cellucid:::.encode_uint32(length(header)), header, cellucid:::.encode_uint32(cells))
  elapsed <- system.time(parsed <- cellucid:::.parse_viewer_binary_event(frame))[["elapsed"]]
  expect_identical(parsed$payload$cells, cells)
  # Vectorized: two million indices decode in well under a second.
  expect_lt(elapsed, 2)
})

test_that("commands and actions with wrong fields are refused before they are sent", {
  frame <- cellucid:::.viewer_command_frame
  expect_error(frame(list(type = "ping")), "carries exactly type, requestId")
  expect_error(frame(list(type = "explode")), "Unknown viewer command")
  # Every change is an action: the commands of schema 2 are gone.
  for (gone in c("highlight", "setColorBy", "setVisibility", "setDataMapping", "setGeneSet", "setField", "resetCamera", "clearHighlights")) {
    expect_error(frame(list(type = gone)), "Unknown viewer command", info = gone)
  }
  expect_error(frame(list(type = "action", requestId = "r", action = "color.by", params = list())), "Unknown viewer command")
  encode <- cellucid:::.encode_viewer_action
  expect_error(encode("", "color.by", list(a = 1)), "requestId")
  expect_error(encode("r1", "color", list(a = 1)), "dotted type")
  expect_error(encode("r1", "Color.By", list(a = 1)), "dotted type")
  expect_error(encode("r1", "color.by", list(1)), "one named list")
  # A turn is named with its label, or neither.
  expect_error(encode("r1", "color.by", list(a = 1), "r-1", NULL), "turn and the turn's label together")
  expect_error(encode("r1", "color.by", list(a = 1), "", "QC"), "exact text")
  expect_match(encode("r1", "color.by", list(a = 1)), "\"turn\":null,\"turnLabel\":null}$")
  load <- cellucid:::.encode_viewer_load_session
  expect_error(load("r1", "a/b.cellucid-session", as.raw(1)), "without a folder")
  expect_error(load("r1", "b.cellucid-session", raw(0)), "session file's bytes")
  # A typed array goes under a key, never directly in a list.
  expect_error(encode("r1", "x.y", list(items = list(cellucid:::.typed_array(1, "uint32")))), "put it under a key")
  expect_error(encode("r1", "x.y", list(a = NA)), "NA")
  expect_error(encode("r1", "x.y", list(a = Inf)), "finite")
  expect_error(cellucid:::.typed_array(c(1, 256), "uint8"), "from 0 to 255")
  expect_error(cellucid:::.typed_array(c(-1, 2), "uint32"), "whole numbers")
  expect_error(cellucid:::.typed_array(1.5, "uint16"), "whole numbers")
  expect_error(cellucid:::.typed_array(1, "int64"), "type is one of")
  # Floats: NA is NaN, no value.
  bytes <- encode("r1", "x.y", list(values = cellucid:::.typed_array(c(1, NA), "float64")))
  header_length <- readBin(bytes[1:4], "integer", size = 4L, endian = "little")
  expect_true(is.nan(readBin(bytes[(4L + header_length + 9L):length(bytes)], "double", size = 8L, endian = "little")))
})

test_that("an action of millions of cells is one binary frame, made in seconds", {
  cells <- seq(0, 1e7 - 1, by = 1)
  elapsed <- system.time(
    frame <- cellucid:::.encode_viewer_action("r1", "highlight.cells", list(cells = cellucid:::.typed_array(cells, "uint32")))
  )[["elapsed"]]
  header_length <- readBin(frame[1:4], "integer", size = 4L, endian = "little")
  expect_identical(length(frame), 4L + header_length + 40000000L)
  body <- frame[(5L + header_length):length(frame)]
  expect_identical(readBin(body, "integer", n = 1e7, size = 4L, endian = "little"), as.integer(cells))
  expect_lt(elapsed, 5)
})
