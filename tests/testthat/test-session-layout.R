# Session format 4's chunk layout, read as the web app writes it
# (app/session/history-file.js; app/session/README.md, "The file"): steps,
# parts and values cut into pieces, each chunk gzip'd or stored as it is,
# every sub-value of at least 512 characters of skeleton text once as a part,
# and the arrays and parts only step parameters name in the lazy chunks. The
# files are written by the test helper the same way (helper-session-files.R)
# and then taken apart -- cellucid-python's tests/test_session_layout.py, test
# for test.

layout_genes <- as.list(sprintf("GENE%04d", 0:119))
layout_other_genes <- as.list(sprintf("OTHER%04d", 0:119))

layout_gene_sets <- function(...) {
  list(sets = lapply(c(...), function(name) list(name = name, genes = if (identical(name, "Third")) layout_other_genes else layout_genes)))
}

# Three steps whose gene-sets values share the first set's genes, and a step
# whose parameters hold a cell array no shown value holds.
layout_session <- function(...) {
  steps <- list(
    list(set = list(`gene-sets` = layout_gene_sets("Markers", "More")), type = "gene-set.add",
         params = list(name = "More", genes = layout_genes)),
    list(set = list(highlights = session_highlights(session_page("p0", "Page 1", session_group("g0", "Lasso", c(0, 1, 2))))),
         type = "highlight.add", params = list(cells = session_cells(c(0, 1, 2, 3)))),
    list(set = list(`gene-sets` = layout_gene_sets("Markers", "More", "Third")), type = "gene-set.add", params = list(name = "Third"))
  )
  session_file(tempfile(fileext = ".cellucid-session"), session_fingerprint(6L, length(layout_genes), dataset_id = "layout"),
               list(list(id = "t1", name = "Tracker 1", root = session_state(`gene-sets` = layout_gene_sets("Markers")), steps = steps)), ...)
}

layout_rewrite <- function(path, edit) {
  file <- edit(session_unframe(readBin(path, "raw", file.size(path))))
  out <- tempfile(fileext = ".cellucid-session")
  writeBin(session_frame(file$manifest, file$chunks), out)
  out
}

layout_index <- function(file, id) which(vapply(file$manifest$chunks, function(entry) entry$id, "") == id)

layout_json <- function(file, id) {
  at <- layout_index(file, id)
  jsonlite::parse_json(rawToChar(session_payload(file$manifest$chunks[[at]], file$chunks[[at]])), simplifyVector = FALSE)
}

layout_put_json <- function(file, id, payload) {
  at <- layout_index(file, id)
  entry <- file$manifest$chunks[[at]]
  stored <- session_store(entry[setdiff(names(entry), c("codec", "storedBytes", "uncompressedBytes"))], session_json(payload), entry$codec)
  file$manifest$chunks[[at]] <- stored$meta
  file$chunks[[at]] <- stored$bytes
  file
}

# A chunk put where the order puts it (eager before lazy, then by kind).
layout_insert <- function(file, meta, payload) {
  rank <- function(entry) (entry$priority == "lazy") * 10 + match(entry$kind, cellucid:::.session_kind_order)
  ranks <- vapply(file$manifest$chunks, rank, 0)
  at <- sum(ranks <= rank(meta))
  stored <- session_store(meta, payload)
  file$manifest$chunks <- append(file$manifest$chunks, list(stored$meta), after = at)
  file$chunks <- append(file$chunks, list(stored$bytes), after = at)
  file
}

# An eager value's hash and the first part its skeleton names.
layout_eager_part <- function(file) {
  values <- layout_json(file, "values/eager/0")
  for (hash in names(values)) {
    names <- session_part_names(values[[hash]])
    if (length(names) > 0L) return(list(hash = hash, name = names[[1L]]))
  }
  stop("no eager value names a part")
}

layout_refused <- function(path) {
  error <- tryCatch({
    cellucid_read_session(path)
    NULL
  }, error = function(error) conditionMessage(error))
  expect_false(is.null(error))
  expect_true(startsWith(error, "This session file is damaged: "), info = error)
  error
}

layout_json_values <- function(value) jsonlite::toJSON(value, auto_unbox = TRUE, digits = NA, null = "null")

layout_reading <- function(path) {
  session <- cellucid_read_session(path)
  tracker <- cellucid_tracker(session)
  list(
    steps = cellucid_steps(tracker),
    states = lapply(c("s0", "s1", "s2", "s3"), function(step) layout_json_values(cellucid_step_state(tracker, step))),
    cells = as.numeric(cellucid_step(tracker, "s2")$params$cells)
  )
}

# --- Parts ------------------------------------------------------------------------

test_that("parts are stored once and joined once", {
  path <- layout_session()
  file <- session_unframe(readBin(path, "raw", file.size(path)))
  kinds <- vapply(file$manifest$chunks, function(entry) entry$kind, "")
  expect_true("parts" %in% kinds)
  parts <- list()
  for (at in which(kinds == "parts")) parts <- c(parts, layout_json(file, file$manifest$chunks[[at]]$id))
  # The genes of "Markers" and "More" in three values and s1's parameters: one part; "Third"'s another.
  expect_identical(sort(vapply(parts, function(part) part[[1L]], "", USE.NAMES = FALSE)), c("GENE0000", "OTHER0000"))
  session <- cellucid_read_session(path)
  tracker <- cellucid_tracker(session)
  expect_identical(layout_json_values(cellucid_step_state(tracker, "s0")[["gene-sets"]]), layout_json_values(layout_gene_sets("Markers")))
  expect_identical(layout_json_values(cellucid_step_state(tracker, "s3")[["gene-sets"]]),
                   layout_json_values(layout_gene_sets("Markers", "More", "Third")))
  expect_identical(unlist(cellucid_step(tracker, "s1")$params$genes), unlist(layout_genes))
})

test_that("a part that names itself through others is refused", {
  path <- layout_rewrite(layout_session(), function(file) {
    name <- layout_eager_part(file)$name
    parts <- layout_json(file, "parts/eager/0")
    other <- strrep("f", 64L)
    parts[[other]] <- list(list(`$part` = name))
    parts[[name]] <- list(loop = list(`$part` = other))
    layout_put_json(file, "parts/eager/0", parts)
  })
  expect_match(layout_refused(path), "its eager values: a stored value names a part inside itself.", fixed = TRUE)
})

test_that("a value that would expand past a page is refused before it is expanded", {
  path <- layout_rewrite(layout_session(), function(file) {
    hash <- layout_eager_part(file)$hash
    parts <- layout_json(file, "parts/eager/0")
    names <- sprintf("%064x", 0:23)
    parts[[names[[1L]]]] <- list(strrep("x", 1000L))
    for (level in 2:24) parts[[names[[level]]]] <- list(list(`$part` = names[[level - 1L]]), list(`$part` = names[[level - 1L]]))
    file <- layout_put_json(file, "parts/eager/0", parts)
    values <- layout_json(file, "values/eager/0")
    values[[hash]] <- list(sets = list(list(`$part` = names[[24L]])))
    layout_put_json(file, "values/eager/0", values)
  })
  started <- proc.time()[["elapsed"]]
  error <- layout_refused(path)
  # 2^23 copies of 1000 characters: counted, never made.
  expect_lt(proc.time()[["elapsed"]] - started, 10)
  expect_match(error, "a stored value expands to", fixed = TRUE)
  expect_match(error, "more than a page can hold (536870888)", fixed = TRUE)
  expect_identical(cellucid:::.session_max_value_text, 2^29 - 24)
})

test_that("a part nothing refers to is refused", {
  path <- layout_rewrite(layout_session(), function(file) {
    parts <- layout_json(file, "parts/eager/0")
    parts[[strrep("0", 64L)]] <- list("unused")
    layout_put_json(file, "parts/eager/0", parts)
  })
  expect_match(layout_refused(path), "it holds a part nothing refers to.$")
})

test_that("a part's name is never trusted", {
  original <- layout_session()
  file <- session_unframe(readBin(original, "raw", file.size(original)))
  name <- layout_eager_part(file)$name
  renamed <- layout_rewrite(original, function(file) {
    for (at in seq_along(file$chunks)) {
      entry <- file$manifest$chunks[[at]]
      if (!entry$kind %in% c("steps", "parts", "values")) next
      text <- gsub(name, strrep("e", 64L), rawToChar(session_payload(entry, file$chunks[[at]])), fixed = TRUE)
      stored <- session_store(entry[setdiff(names(entry), c("codec", "storedBytes", "uncompressedBytes"))], charToRaw(text), entry$codec)
      file$manifest$chunks[[at]] <- stored$meta
      file$chunks[[at]] <- stored$bytes
    }
    file
  })
  expect_identical(layout_reading(renamed)$states, layout_reading(original)$states)
  tampered <- layout_rewrite(original, function(file) {
    parts <- layout_json(file, "parts/eager/0")
    parts[[name]][[1L]] <- "Tampered"
    layout_put_json(file, "parts/eager/0", parts)
  })
  expect_match(layout_refused(tampered), "its eager values holds a value that does not match its hash.$")
})

test_that("a shown value naming a lazy part is refused", {
  path <- layout_rewrite(layout_session(), function(file) {
    name <- layout_eager_part(file)$name
    parts <- layout_json(file, "parts/eager/0")
    moved <- parts[name]
    parts[[name]] <- NULL
    file <- layout_put_json(file, "parts/eager/0", parts)
    layout_insert(file, list(id = "parts/lazy/9", kind = "parts", priority = "lazy"), session_json(moved))
  })
  expect_match(layout_refused(path), "a value of the step it shows names a part that is not among its first chunks", fixed = TRUE)
})

test_that("a part the file lacks is refused", {
  path <- layout_rewrite(layout_session(), function(file) {
    name <- layout_eager_part(file)$name
    parts <- layout_json(file, "parts/eager/0")
    parts[[name]] <- NULL
    layout_put_json(file, "parts/eager/0", parts)
  })
  expect_match(layout_refused(path), "a value names a part the file does not hold", fixed = TRUE)
})

test_that("a malformed part reference is refused", {
  for (case in list(list(list(`$part` = 5L), "a stored value names a part wrongly."),
                    # The placeholder's shape is checked before the part is looked for.
                    list(list(`$part` = "not-a-hash"), "a stored value names a part wrongly."),
                    list(list(`$part` = ""), "a stored value names a part wrongly."),
                    list(list(`$part` = strrep("f", 64L)), "a value names a part the file does not hold (ffffffffffff"),
                    list(list(`$array` = ""), "a value names an array the file does not hold ("),
                    list(list(`$part` = strrep("a", 64L), x = 1L), "a stored value has a malformed placeholder."))) {
    path <- layout_rewrite(layout_session(), function(file) {
      hash <- layout_eager_part(file)$hash
      values <- layout_json(file, "values/eager/0")
      values[[hash]] <- list(sets = list(case[[1L]]))
      layout_put_json(file, "values/eager/0", values)
    })
    expect_match(layout_refused(path), case[[2L]], fixed = TRUE)
  }
})

# --- Pieces, codecs and lazy parameters -------------------------------------------

test_that("steps and values in many pieces read as one", {
  whole <- layout_session()
  cut <- layout_session(piece_chars = 600)
  file <- session_unframe(readBin(cut, "raw", file.size(cut)))
  ids <- vapply(file$manifest$chunks, function(entry) entry$id, "")
  expect_gt(length(file$manifest$trackers[[1L]]$steps), 1L)
  expect_identical(unlist(file$manifest$trackers[[1L]]$steps), ids[startsWith(ids, "steps/")])
  expect_gt(sum(startsWith(ids, "values/")), 2L)
  expect_gt(sum(startsWith(ids, "parts/")), 1L)
  expect_identical(layout_reading(cut), layout_reading(whole))
})

test_that("chunks stored as they are read the same", {
  plain <- layout_session(codec = "none")
  file <- session_unframe(readBin(plain, "raw", file.size(plain)))
  expect_identical(unique(vapply(file$manifest$chunks, function(entry) entry$codec, "")), "none")
  expect_true(all(vapply(file$manifest$chunks, function(entry) entry$storedBytes == entry$uncompressedBytes, TRUE)))
  expect_identical(layout_reading(plain), layout_reading(layout_session(codec = "gzip")))
  dishonest <- layout_rewrite(plain, function(file) {
    file$manifest$chunks[[1L]]$uncompressedBytes <- file$manifest$chunks[[1L]]$uncompressedBytes + 1
    file
  })
  expect_match(layout_refused(dishonest), "chunk 0 is stored as it is, at another size.$")
})

test_that("parameters naming a lazy array are joined", {
  path <- layout_session()
  file <- session_unframe(readBin(path, "raw", file.size(path)))
  digest <- session_cells(c(0, 1, 2, 3))$digest
  entry <- Filter(function(entry) identical(entry$kind, "array") && identical(entry$digest, digest), file$manifest$chunks)[[1L]]
  expect_identical(entry$priority, "lazy")
  cells <- cellucid_step(cellucid_tracker(cellucid_read_session(path)), "s2")$params$cells
  expect_identical(as.numeric(cells), c(0, 1, 2, 3))
  expect_identical(attr(cells, "cellucid_array"), "Uint32Array")
})

test_that("steps chunks named in another order than the file holds them are refused", {
  path <- layout_rewrite(layout_session(piece_chars = 600), function(file) {
    file$manifest$trackers[[1L]]$steps[1:2] <- file$manifest$trackers[[1L]]$steps[2:1]
    file
  })
  expect_match(layout_refused(path), paste0('tracker "', session_tracker_id("t1"), '" names its steps chunks in another order than the file holds them.'),
               fixed = TRUE)
})

test_that("a steps-only file holds the parts and arrays its steps name", {
  path <- layout_session(content = "steps")
  file <- session_unframe(readBin(path, "raw", file.size(path)))
  kinds <- vapply(file$manifest$chunks, function(entry) entry$kind, "")
  expect_false("values" %in% kinds)
  expect_true(all(c("parts", "array") %in% kinds))
  tracker <- cellucid_tracker(cellucid_read_session(path))
  expect_identical(unlist(cellucid_step(tracker, "s1")$params$genes), unlist(layout_genes))
  expect_identical(as.numeric(cellucid_step(tracker, "s2")$params$cells), c(0, 1, 2, 3))
})

test_that("a broken layout is refused", {
  cases <- list(
    list(function(file) { file$manifest$trackers[[1L]]$steps <- list(); file }, "tracker 0 steps must list its steps chunks."),
    list(function(file) { file$manifest$trackers[[1L]]$steps <- "steps/0/0"; file }, "tracker 0 steps must list its steps chunks."),
    list(function(file) {
      file$manifest$trackers[[1L]]$steps <- c(file$manifest$trackers[[1L]]$steps, list("values/eager/0"))
      file
    }, "names a steps chunk that is not its own."),
    list(function(file) {
      file$manifest$trackers[[1L]]$steps <- c(file$manifest$trackers[[1L]]$steps, file$manifest$trackers[[1L]]$steps[1L])
      file
    }, "names a steps chunk that is not its own."),
    list(function(file) {
      at <- layout_index(file, "values/eager/0")
      file$manifest$chunks <- file$manifest$chunks[-at]
      file$chunks <- file$chunks[-at]
      file
    }, "it holds no values of the step it shows."),
    list(function(file) {
      at <- layout_index(file, "values/eager/0")
      order <- seq_along(file$chunks)
      order[c(at - 1L, at)] <- order[c(at, at - 1L)]
      file$manifest$chunks <- file$manifest$chunks[order]
      file$chunks <- file$chunks[order]
      file
    }, "is out of order (eager before lazy; steps, parts, array, values, journal, chat).")
  )
  for (case in cases) expect_match(layout_refused(layout_rewrite(layout_session(), case[[1L]])), case[[2L]], fixed = TRUE)
})

test_that("the layout constants are the web app's", {
  writer <- cellucid_sibling_file("cellucid", "assets", "js", "app", "session", "history-file.js")
  parts <- cellucid_sibling_file("cellucid", "assets", "js", "app", "history", "value-parts.js")
  expect_true(grepl("const KIND_ORDER = Object.freeze(['steps', 'parts', 'array', 'values', 'journal', 'chat']);",
                    paste(readLines(writer, warn = FALSE), collapse = "\n"), fixed = TRUE))
  expect_identical(cellucid:::.session_kind_order, c("steps", "parts", "array", "values", "journal", "chat"))
  text <- paste(readLines(parts, warn = FALSE), collapse = "\n")
  expect_true(grepl("export const MAX_VALUE_TEXT = 2 ** 29 - 24;", text, fixed = TRUE))
  expect_true(grepl("export const PART_MIN_CHARS = 512;", text, fixed = TRUE))
  expect_identical(session_part_min_chars, 512L)
})
