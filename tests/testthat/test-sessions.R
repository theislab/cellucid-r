# Session format 4 (History v2) in R: the golden history read exactly, and
# applied -- cellucid-python's tests/test_sessions.py, case for case.
#
# The shared golden fixture (fixtures/history_golden, written by the web
# repository's scripts/generate-history-golden.mjs with the history model's
# own constructors and pinned by digest in the web app, cellucid-python and
# cellucid-r) holds every kind of record the model has: trackers (a replay,
# one removed, one paused), branches (a nested fork, a locked and an archived
# one), turns (a notebook's, an agent's with its summary), a group, a squashed
# step of what changed while recording was paused, notes, a bookmark, a deleted
# step and its journal entry, a flagged step, a stale result, each step's
# parameters version. history_golden.expected.json says what reading it must give, and
# the apply plans of a sequence of applies on one object: these tests compare
# field by field. Synthesized files (helper-session-files.R) test every other
# rule: what the share left out, damaged and newer files, and each apply rule
# on its own.

golden_path <- function(name = "history_golden.cellucid-session") testthat::test_path("fixtures", "history_golden", name)
golden_expected <- function() jsonlite::read_json(golden_path("history_golden.expected.json"), simplifyVector = FALSE)
# The golden files' digests; the web app (tests/history-golden-fixture.test.mjs) and cellucid-python pin the same.
.golden_sha256 <- c(
  "history_golden.cellucid-session" = "9711bbc4dfd19fd91809a2ab3967cab54b6982c2a3907c73432200e72414ca4a",
  "history_golden.expected.json" = "8684353157ec936017b777c395e1d60342dcd4f393ef2a0648fe4a7891a67c34",
  "history_golden_state.cellucid-session" = "f9f0a6025679dacb91cbd00b0abd1cd6380ab1e829fa9ea529cbc59c1801e3dd"
)
golden <- function() cellucid_read_session(golden_path())

# The dataset the golden session was saved on, as expected.json describes it:
# a SingleCellExperiment (its reducedDim umap is obsm X_umap).
golden_sce <- function() {
  dataset <- golden_expected()$dataset
  cells <- unlist(dataset$obs_names)
  genes <- unlist(dataset$var_names)
  obs <- S4Vectors::DataFrame(row.names = cells)
  for (name in names(dataset$obs)) {
    column <- dataset$obs[[name]]
    categories <- unlist(column$categories)
    obs[[name]] <- factor(categories[unlist(column$codes) + 1L], levels = categories)
  }
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = Matrix::Matrix(0, length(genes), length(cells), sparse = TRUE, dimnames = list(genes, cells))),
    colData = obs
  )
  SingleCellExperiment::reducedDim(object, "umap") <- do.call(rbind, lapply(dataset$obsm$X_umap, function(row) as.double(unlist(row))))
  object
}

# A plan's rows as the shared expected plans write them.
plan_rows <- function(plan) {
  rows <- plan$rows
  lapply(seq_len(nrow(rows)), function(at) list(
    kind = rows$kind[[at]], domain = rows$domain[[at]], name = rows$name[[at]],
    column = if (is.na(rows$column[[at]])) NULL else rows$column[[at]], from = as.list(rows$from[[at]]),
    cells = if (is.na(rows$cells[[at]])) NULL else rows$cells[[at]], status = rows$status[[at]],
    note = if (is.na(rows$note[[at]])) NULL else rows$note[[at]]
  ))
}

as_json_value <- function(value) jsonlite::parse_json(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = NA),
                                                       simplifyVector = FALSE)

# =============================================================================
# The golden fixture, read
# =============================================================================

test_that("the golden files are pinned, the web app's and cellucid-python's", {
  for (name in names(.golden_sha256)) expect_identical(cellucid_file_sha256(golden_path(name)), .golden_sha256[[name]], info = name)
  expect_identical(golden_expected()$sha256, .golden_sha256[["history_golden.cellucid-session"]])
  expect_identical(cellucid_file_sha256(system.file("extdata", "history_golden.cellucid-session", package = "cellucid")),
                   .golden_sha256[["history_golden.cellucid-session"]])
  copies <- list(
    list(c("cellucid", "tests", "fixtures", "history-golden.cellucid-session"), c("cellucid", "tests", "fixtures", "history-golden.expected.json"),
         c("cellucid", "tests", "fixtures", "history-golden-state.cellucid-session")),
    lapply(names(.golden_sha256), function(name) c("cellucid-python", "tests", "fixtures", "history_golden", name))
  )
  for (pair in copies) {
    paths <- vapply(pair, function(copy) do.call(cellucid_sibling_file, as.list(copy)), "")
    expect_identical(unname(vapply(paths, cellucid_file_sha256, "")), unname(.golden_sha256), info = paths[[1L]])
  }
})

test_that("the web app's current state only reads as one state", {
  # Share > Current state only: one tracker of the opening and one squashed step, nobody named, no time, no note.
  state <- cellucid_read_session(golden_path("history_golden_state.cellucid-session"))
  expect_identical(list(state$scope$content, state$scope$times, state$scope$notes, state$scope$chats), list("state", "none", FALSE, FALSE))
  expect_identical(state$created_at, "1970-01-01T00:00:00.000Z")
  expect_length(state$trackers, 1L)
  tracker <- cellucid_tracker(state)
  expect_null(tracker$created_at)
  expect_length(tracker$steps, 2L)
  root <- tracker$steps[[1L]]
  shared <- tracker$steps[[2L]]
  expect_identical(list(root$kind, root$by$origin, shared$kind, shared$squashed$cause, shared$by$origin),
                   list("root", "startup", "squashed", "share", "user"))
  expect_identical(tracker$head_id, shared$id)
  # What the golden shows, value for value; a user field's creation time is left out as every time is.
  session <- golden()
  was_tracker <- cellucid_tracker(session)
  shown <- cellucid:::.history_hashes(tracker, shared$id)
  was <- cellucid:::.history_hashes(was_tracker, was_tracker$head_id)
  expect_setequal(names(shown), names(was))
  for (domain in setdiff(names(shown), "field-overlays")) expect_identical(shown[[domain]], was[[domain]], info = domain)
  fields <- cellucid_step_state(tracker, shared$id)[["field-overlays"]]$userFields
  expect_true(all(vapply(fields, function(field) is.null(field$createdAt), TRUE)))
})

test_that("the header reads as written", {
  expected <- golden_expected()$session
  session <- golden()
  expect_identical(session$created_at, expected$created_at)
  expect_identical(c(session$app$version, session$app$catalog), c(expected$app$version, expected$app$catalog))
  expect_identical(vapply(session$actions, function(stamp) stamp$type, ""), unlist(expected$actions))
  expect_true(all(grepl("^[0-9a-f]{64}$", vapply(session$actions, function(stamp) stamp$params, ""))))
  for (key in c("content", "results", "notes", "chats", "journal", "times", "source")) {
    expect_identical(session$scope[[key]], expected$scope[[key]], info = key)
  }
  expect_identical(unlist(session$scope$trackers), unlist(expected$scope$trackers))
  expect_identical(as_json_value(lapply(session$requires, as.list)), expected$requires)
  # Page chrome outside the history (the pressed selection tool): read, never applied.
  expect_identical(session$chrome, expected$chrome)
  expect_identical(session$active_tracker, expected$active_tracker)
  expect_identical(session$sha256, .golden_sha256[["history_golden.cellucid-session"]])
  expect_identical(session$dataset_fingerprint$cellOrder$digest, golden_expected()$dataset$cell_order_digest)
  expect_identical(session$source, list(kind = "sample", dataset = golden_expected()$dataset$id))
  expect_null(session$cells)
  expect_null(session$genes)
  expect_output(print(session), "Cellucid session")
})

test_that("every tracker reads as expected", {
  session <- golden()
  for (expected in golden_expected()$trackers) {
    tracker <- cellucid_tracker(session, expected$id)
    expect_s3_class(tracker, "cellucid_tracker")
    expect_identical(list(tracker$id, tracker$name, !is.null(tracker$removed), tracker$start$kind),
                     list(expected$id, expected$name, expected$removed, expected$start))
    expect_identical(c(tracker$head_id, tracker$active_branch, tracker$main_branch),
                     c(expected$head, expected$active_branch, expected$main_branch))
    paused <- tracker$paused
    expect_identical(as_json_value(if (is.null(paused)) NULL else list(at = paused$at, origins = lapply(paused$by, function(by) by$origin),
                                                                           draft = paused$draft)),
                     expected$paused, info = expected$name)
    groups <- cellucid_groups(tracker)
    expect_identical(
      as_json_value(lapply(seq_len(nrow(groups)), function(at) list(
        id = groups$id[[at]], name = groups$name[[at]], branch = groups$branch[[at]], first = groups$first[[at]],
        last = groups$last[[at]], steps = as.list(cellucid_group(tracker, groups$id[[at]])$steps), n_notes = groups$n_notes[[at]]
      ))),
      expected$groups, info = expected$name
    )
    squashed <- Filter(function(step) identical(step$kind, "squashed"), tracker$steps)
    expect_identical(
      as_json_value(lapply(squashed, function(step) list(
        step = step$id, cause = step$squashed$cause, count = step$squashed$count,
        origins = lapply(step$squashed$by, function(by) by$origin), touches = step$touches,
        carried = list(notes = step$squashed$carried$notes, bookmarks = step$squashed$carried$bookmarks)
      ))),
      expected$squashed, info = expected$name
    )
    branches <- cellucid_branches(tracker)
    steps <- cellucid_steps(tracker)
    got <- lapply(seq_len(nrow(branches)), function(at) list(
      id = branches$id[[at]], name = branches$name[[at]], base = if (is.na(branches$base[[at]])) NULL else branches$base[[at]],
      tip = branches$tip[[at]], main = branches$main[[at]], locked = branches$locked[[at]], archived = branches$archived[[at]],
      n_notes = branches$n_notes[[at]], path = as.list(cellucid_steps(tracker, branches$id[[at]])$id),
      own = as.list(steps$id[steps$branch == branches$id[[at]]])
    ))
    expect_identical(as_json_value(got), expected$branches, info = expected$name)
    expect_identical(names(steps), unlist(expected$steps$columns))
    rows <- lapply(seq_len(nrow(steps)), function(at) {
      row <- as.list(steps[at, , drop = FALSE])
      row$at <- if (is.na(row$at)) NA else round(as.numeric(row$at) * 1000)
      lapply(row, function(value) if (length(value) == 1L && is.na(value)) NULL else value)
    })
    expect_identical(as_json_value(lapply(rows, unname)), expected$steps$rows, info = expected$name)
    expect_s3_class(steps$at, "POSIXct")
    turns <- cellucid_turns(tracker)
    expect_identical(
      as_json_value(lapply(seq_len(nrow(turns)), function(at) list(
        id = turns$id[[at]], label = turns$label[[at]], origin = turns$origin[[at]], branch = turns$branch[[at]],
        steps = as.list(steps$id[!is.na(steps$turn) & steps$turn == turns$id[[at]]]),
        summary = if (is.na(turns$summary[[at]])) NULL else turns$summary[[at]]
      ))),
      expected$turns
    )
    notes <- cellucid_notes(tracker)
    expect_identical(names(notes), cellucid:::.history_note_columns)
    expect_identical(
      as_json_value(lapply(seq_len(nrow(notes)), function(at) list(
        id = notes$id[[at]], on = stats::setNames(list(notes$target[[at]]), notes$on[[at]]), text = notes$text[[at]],
        origin = notes$origin[[at]]
      ))),
      expected$notes
    )
    journal <- cellucid_journal(tracker)
    expect_identical(
      as_json_value(lapply(seq_len(nrow(journal)), function(at) list(id = journal$id[[at]], op = journal$op[[at]],
                                                                    label = journal$label[[at]], purged = journal$purged[[at]]))),
      expected$journal
    )
    stale <- lapply(cellucid_stale(tracker), function(entry) list(step = entry$step, cause = entry$cause, domains = as.list(entry$domains)))
    expect_identical(as_json_value(stale), expected$stale)
    expect_identical(as_json_value(lapply(cellucid_conversations(session, tracker), function(entry) entry$conversation)),
                     expected$conversations)
  }
  expect_identical(cellucid_trackers(session)$name, c("Annotation", "Annotation (replay)", "Annotation"))
  expect_output(print(cellucid_tracker(session)), "Annotation")
})

test_that("the state at a step is the one written", {
  session <- golden()
  expected <- golden_expected()
  for (state in expected$states) {
    tracker <- cellucid_tracker(session, state$tracker)
    expect_identical(as.list(cellucid:::.history_hashes(tracker, state$step)), state$values)
    expect_setequal(names(cellucid_step_state(tracker, state$step)), names(state$values))
  }
  head <- cellucid_step_state(cellucid_tracker(session), "s13")
  fields <- head[["field-overlays"]]$userFields
  score <- Filter(function(field) identical(field$key, "score"), fields)[[1L]]
  expect_identical(as.numeric(score$values$cells), as.numeric(unlist(expected$values$score$cells)))
  expect_identical(attr(score$values$values, "cellucid_array"), "Float32Array")
  written <- vapply(expected$values$score$values, function(value) if (is.null(value)) NaN else value, 0)
  expect_equal(as.numeric(score$values$values), written, tolerance = 1e-6)
  merged <- Filter(function(field) identical(field$key, "cell_type (merged)"), fields)[[1L]]
  expect_identical(merged$categories, expected$values$merged$categories)
  expect_identical(as.numeric(merged$codes$values), as.numeric(unlist(expected$values$merged$codes)))
  expect_identical(as.numeric(head$highlights$pages[[1L]]$groups[[1L]]$cells), as.numeric(unlist(expected$values$doublets)))
  # Numbers JSON cannot hold come back as they were written.
  expect_identical(1 / head$camera$views[[1L]]$orbit$theta, -Inf)
  expect_true(is.nan(head$settings$controls$gamma$value))
})

test_that("diffs say what differs, domain by domain", {
  session <- golden()
  for (entry in golden_expected()$diffs) {
    a <- cellucid_tracker(session, entry$a$tracker)
    b <- cellucid_tracker(session, entry$b$tracker)
    diff <- cellucid_step_diff(a, entry$a$step, entry$b$step, other = b)
    expect_identical(names(diff), c("domain", "name", "line"))
    expect_identical(as_json_value(lapply(seq_len(nrow(diff)), function(at) as.list(diff[at, ]))), entry$rows)
    expect_identical(is.null(attr(diff, "common")), !identical(entry$a$tracker, entry$b$tracker))
  }
  tracker <- cellucid_tracker(session)
  expect_identical(attr(cellucid_step_diff(tracker, "s5", "s13"), "common"), "s5")
  expect_identical(attr(cellucid_step_diff(tracker, "s13", "s9"), "common"), "s2")
})

test_that("the methods paragraph and pseudocode are the web app's words", {
  session <- golden()
  expected <- golden_expected()
  for (entry in expected$methods) expect_identical(cellucid_methods(cellucid_tracker(session, entry$tracker), entry$branch), entry$text)
  for (entry in expected$pseudocode) expect_identical(cellucid_pseudocode(cellucid_tracker(session, entry$tracker), entry$branch), entry$text)
  # Main when no branch is named.
  expect_identical(cellucid_methods(cellucid_tracker(session)), expected$methods[[1L]]$text)
})

test_that("records are found by name or id and refused naming the others", {
  session <- golden()
  ids <- vapply(golden_expected()$trackers, function(tracker) tracker$id, "")
  # The removed tracker shares the active one's name: names are among the trackers not removed.
  expect_identical(cellucid_tracker(session, "Annotation")$id, ids[[1L]])
  expect_identical(cellucid_tracker(session, " annotation ")$id, ids[[1L]])
  expect_false(is.null(cellucid_tracker(session, ids[[3L]])$removed))
  expect_error(cellucid_tracker(session, "Nothing"),
               'There is no tracker "Nothing"; the trackers are "Annotation" and "Annotation (replay)".', fixed = TRUE)
  tracker <- cellucid_tracker(session)
  expect_identical(cellucid_branch(tracker, "without SCORE")$id, "b1")
  expect_identical(cellucid_branch(tracker, "b2")$name, "Clusters")
  expect_identical(cellucid_branch(tracker)$id, tracker$active_branch)
  expect_error(cellucid_branch(tracker, "Gone"),
               'There is no branch "Gone"; the branches are "Main", "Without score", "Clusters", "Scratch" and "Tidy fork".', fixed = TRUE)
  expect_error(cellucid_step(tracker, "s11"),
               'There is no step "s11" in the tracker "Annotation"; cellucid_steps(tracker) lists its steps.', fixed = TRUE)
})

test_that("steps carry every field of the model", {
  session <- golden()
  tracker <- cellucid_tracker(session)
  agent <- cellucid_step(tracker, "s5")
  expect_identical(agent$by[c("origin", "agent", "model", "turn")],
                   list(origin = "agent", agent = "Claude", model = "claude-sonnet-4-6", turn = "claude-1"))
  expect_identical(agent$bookmark, "Markers saved")
  expect_identical(vapply(agent$notes, function(note) note$text, ""), "Top hits match the paper.")
  expect_identical(list(agent$kind, agent$v, agent$group, agent$squashed), list("action", 1L, NULL, NULL))
  expect_false(any(c("picture", "picture_png") %in% names(agent)))
  expect_false(exists("cellucid_step_picture", envir = asNamespace("cellucid")))
  expect_identical(list(cellucid_step(tracker, "s0")$kind, cellucid_step(tracker, "s0")$v), list("root", NULL))
  # The notebook's whole turn, grouped.
  qc <- cellucid_group(tracker, "qc")
  expect_identical(list(qc$id, qc$name, qc$steps), list("g0", "QC", c("s2", "s3")))
  expect_identical(cellucid_group(tracker, "g0")$name, "QC")
  expect_identical(vapply(qc$notes, function(note) note$text, ""), "Doublets marked and the score pushed from the notebook.")
  expect_identical(cellucid_step(tracker, "s2")$group, "g0")
  expect_error(cellucid_group(tracker, "Clustering"), 'There is no group "Clustering"; the groups are "QC".', fixed = TRUE)
  # What changed while recording was paused: one squashed step with no action of its own.
  paused <- cellucid_step(tracker, "s14")
  expect_identical(list(paused$kind, paused$type, paused$params, paused$v, paused$code), list("squashed", NULL, NULL, NULL, NULL))
  expect_identical(list(paused$squashed$cause, paused$squashed$count, vapply(paused$squashed$by, function(by) by$origin, "")),
                   list("paused", NULL, c("user", "host")))
  expect_identical(paused$label, "While paused: Colors, Settings")
  expect_identical(cellucid_step_state(tracker, "s14")$settings$controls$`point-size`$value, 5L)
  paused_tracker <- cellucid_tracker(session, "Annotation (replay)")$paused
  expect_identical(list(vapply(paused_tracker$by, function(by) by$origin, ""), paused_tracker$draft), list("user", NULL))
  expect_false(is.na(cellucid_trackers(session)$paused_at[[2L]]))
  expect_false(cellucid_trackers(session)$recording[[2L]])
  # A squash in mid-branch keeps the last replaced step's id; what the replaced steps carried is kept.
  tidied <- cellucid_step(tracker, "s16")
  expect_identical(list(tidied$kind, tidied$squashed$cause, tidied$squashed$count), list("squashed", "steps", 2L))
  expect_identical(unlist(tidied$squashed$carried$bookmarks), "Before tidying")
  expect_length(tidied$squashed$carried$notes, 1L)
  expect_identical(cellucid_step(tracker, "s17")$parent, "s16")
  expect_identical(cellucid_branch(tracker, "Tidy fork")$base, "s16")
  expect_error(cellucid_step(tracker, "s15"), "There is no step")
  expect_setequal(names(cellucid_step(tracker, "s4")$reads), c("cells", "mapping", "expression", "field-overlays", "views", "filters", "visibility", "highlights"))
  camera <- cellucid_step(tracker, "s10")
  expect_identical(list(camera$merged, camera$until - camera$at, camera$ms, camera$coalesce_key), list(3L, 2000, NULL, "camera:View 1"))
  flagged <- cellucid_step(tracker, "s12")
  expect_identical(flagged$status, "unapplicable")
  expect_identical(flagged$before, flagged$after)
  expect_identical(flagged$reason, 'There is no highlight "Gate" on page "QC"; its highlights are "Doublets".')
  replay <- cellucid_tracker(session, "Annotation (replay)")
  expect_identical(cellucid_step(replay, "s1")$replay_of[c("tracker", "step")], list(tracker = tracker$id, step = "s1"))
  expect_identical(replay$start$from$branch$name, "Main")
  expect_identical(vapply(replay$start$report, function(entry) entry$outcome, ""), c("applied", "applied", "unmatched"))
  expect_identical(replay$start$mapping$genes[[1L]]$from, "Cd3e")
  expect_identical(cellucid_journal(tracker)$op,
                   c("note.add", "bookmark.set", "group.create", "branch.lock", "branch.archive", "step.delete", "squash"))
  grouped <- tracker$journal[[3L]]$after
  # A group edit records which steps the group held (memberships), never the steps.
  expect_identical(vapply(grouped$groups, function(group) group$name, ""), "QC")
  expect_length(grouped$steps, 0L)
  expect_identical(grouped$memberships, list(list(group = "g0", from = "s2", to = "s3")))
  expect_identical(vapply(tracker$journal[[6L]]$before$steps, function(step) step$id, ""), c("s6", "s11", "s12"))
})

# =============================================================================
# Applying the golden fixture: the shared plans
# =============================================================================

test_that("the apply sequence gives the shared plans", {
  skip_without_sce()
  session <- golden()
  expected <- golden_expected()
  object <- golden_sce()
  expect_identical(cellucid:::.apply_cell_order_digest(SingleCellExperiment::reducedDim(object, "umap")),
                   expected$dataset$cell_order_digest)
  first <- NULL
  for (entry in expected$apply$sequence) {
    plan <- cellucid_apply_session(object, session, tracker = entry$tracker, branch = entry$branch, step = entry$step)
    expect_s3_class(plan, "cellucid_apply_plan")
    expect_identical(list(plan$tracker, plan$branch, plan$step), unname(entry$expect))
    expect_identical(plan$session_sha256, .golden_sha256[["history_golden.cellucid-session"]])
    expect_identical(as_json_value(lapply(plan$requires, as.list)), entry$requires)
    expect_identical(as_json_value(plan_rows(plan)), entry$rows, info = entry$name)
    expect_identical(names(plan$rows), c("kind", "domain", "name", "column", "from", "cells", "status", "note"))
    object <- cellucid_apply_plan(plan)
    report <- attr(object, "cellucid_apply_report")
    expect_identical(report$applied, entry$applied, info = entry$name)
    expect_identical(is.null(report$provenance), !entry$applied)
    # The gene sets' columns of the var the mapping reads (rowData), the web app's and cellucid-python's.
    written <- attr(report$plan, "cellucid_context")$writes$var_columns
    expect_identical(lapply(written, function(column) list(name = column$name, values = as.list(column$values))),
                     entry$var_columns, info = entry$name)
    for (column in entry$var_columns) {
      expect_identical(SummarizedExperiment::rowData(object)[[column$name]], unlist(column$values), info = entry$name)
    }
    if (is.null(first)) first <- report
  }
  provenance <- first$provenance
  expect_identical(as_json_value(provenance[setdiff(names(provenance), c("applied_at", "cellucid_version", "package"))]),
                   expected$apply$provenance)
  expect_identical(provenance$package, "r")
  expect_identical(provenance$cellucid_version, as.character(utils::packageVersion("cellucid")))
  expect_match(provenance$applied_at, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}\\.\\d{3}Z$")
  # What the object holds now.
  data <- SummarizedExperiment::colData(object)
  expect_identical(levels(data$cell_type), c("B cell", "T cell", "NK"))
  expect_true("cluster" %in% names(data))
  expect_false("leiden" %in% names(data))
  expect_identical(sum(data$cellucid_highlight__QC__Doublets), length(expected$values$doublets))
  expect_identical(sum(data$cellucid_highlight__Types__cell_type_T_cell), length(expected$values$t_cells))
  expect_type(data$score, "double")
  expect_identical(levels(data[["cell_type (merged)"]]), unlist(expected$values$merged$categories))
  gene_set <- S4Vectors::metadata(object)$cellucid$gene_sets[["T cell markers"]]
  expect_identical(list(gene_set$genes, gene_set$origin), list(c("CD3E", "NKG7"), "session"))
  expect_identical(names(SummarizedExperiment::rowData(object)), "cellucid_gene_set__T_cell_markers")
  records <- cellucid_history_applied(object)
  expect_identical(lapply(records, function(record) c(record$step, record$branch)), list(c("s13", "b0"), c("s8", "b1")))
})

test_that("applying again changes nothing, and dry_run = FALSE applies at once", {
  skip_without_sce()
  object <- cellucid_apply_session(golden_sce(), golden_path(), dry_run = FALSE)
  expect_true(attr(object, "cellucid_apply_report")$applied)
  again <- cellucid_apply_session(object, golden(), dry_run = FALSE)
  report <- attr(again, "cellucid_apply_report")
  expect_false(report$applied)
  expect_null(report$provenance)
  expect_setequal(unique(report$plan$rows$status), c("unchanged", "not-applied"))
  expect_length(cellucid_history_applied(again), 1L)
  expect_output(print(report$plan), "Cellucid apply plan")
})

test_that("keep_session keeps the file to rebuild the history", {
  skip_without_sce()
  object <- cellucid_apply_session(golden_sce(), golden(), dry_run = FALSE, keep_session = TRUE)
  provenance <- attr(object, "cellucid_apply_report")$provenance
  rebuilt <- tempfile(fileext = ".cellucid-session")
  writeBin(jsonlite::base64_dec(provenance$session), rebuilt)
  expect_identical(cellucid_file_sha256(rebuilt), .golden_sha256[["history_golden.cellucid-session"]])
  expect_identical(cellucid_steps(cellucid_tracker(cellucid_read_session(rebuilt)))$id, cellucid_steps(cellucid_tracker(golden()))$id)
})

test_that("a column the person made is never replaced", {
  skip_without_sce()
  object <- golden_sce()
  object$score <- as.double(seq_len(ncol(object)) - 1L)
  plan <- cellucid_apply_session(object, golden())
  row <- plan$rows[plan$rows$name == "score", ]
  expect_identical(c(row$status, row$note), c(
    "not-applied", 'The object already has an obs column "score" that Cellucid did not write; it is never replaced.'
  ))
  applied <- cellucid_apply_plan(plan)
  provenance <- attr(applied, "cellucid_apply_report")$provenance
  expect_identical(applied$score, as.double(seq_len(ncol(object)) - 1L))
  expect_false("score" %in% unlist(provenance$columns))
  expect_true(any(vapply(provenance$skipped, function(item) identical(item$kind, "user-field") && identical(item$name, "score"), TRUE)))
})

test_that("a plan refuses an object changed since it was made", {
  skip_without_sce()
  object <- golden_sce()
  plan <- cellucid_apply_session(object, golden())
  context <- attr(plan, "cellucid_context")
  context$data$score <- 1
  attr(plan, "cellucid_context") <- context
  expect_error(cellucid_apply_plan(plan), "The object changed since this plan was made")
})

test_that("the cells must be the session's, in its order", {
  skip_without_sce()
  object <- golden_sce()
  shuffled <- object[, c(2L, 1L, 3:ncol(object))]
  expect_error(cellucid_apply_session(shuffled, golden()),
               "This object's cells are not the cells the session was saved on, in the same order")
  # Other coordinates (an embedding computed again) are other cells too: the digest is over the stored values.
  moved <- object
  SingleCellExperiment::reducedDim(moved, "umap") <- SingleCellExperiment::reducedDim(moved, "umap") + 1
  expect_error(cellucid_apply_session(moved, golden()), "cell-order digest of obsm['X_umap'] (columns 0, 1)", fixed = TRUE)
  # The digest is the web app's on every cell count, 32-bit words exact.
  values <- matrix(c(0, -0, NaN, 1e30, -1.5, 3), ncol = 2L)
  expect_identical(cellucid:::.apply_cell_order_digest(values), cellucid:::.apply_cell_order_digest(values))
  expect_match(cellucid:::.apply_cell_order_digest(values[, 1L, drop = FALSE], c(7, 2^32 - 1, 0)), "^[0-9a-f]{16}$")
})

test_that("the dataset must be the session's", {
  skip_without_sce()
  object <- golden_sce()
  S4Vectors::metadata(object)$cellucid <- list(id = "another")
  expect_error(cellucid_apply_session(object, golden()), "saved on the dataset 'history_golden', and this object is 'another'")
  expect_error(cellucid_apply_session(golden_sce()[, 1:599], golden()), "saved on 600 cells, and this object has 599")
})

test_that("only what the chosen step uses is required", {
  skip_without_sce()
  object <- golden_sce()
  names(SummarizedExperiment::colData(object))[names(SummarizedExperiment::colData(object)) == "cell_type"] <- "kind"
  # The bookmark colors by cell_type; the head colors by a gene.
  expect_error(cellucid_apply_session(object, golden(), step = "s5"), "The state at step s5 uses the obs column 'cell_type'")
  head <- cellucid_apply_session(object, golden())
  renamed <- head$rows[head$rows$kind == "categories-renamed", ]
  expect_identical(c(renamed$status, renamed$note), c("not-applied", 'The object has no obs column "cell_type".'))
})

test_that("a step off the branch named is refused", {
  skip_without_sce()
  expect_error(cellucid_apply_session(golden_sce(), golden(), branch = "Clusters", step = "s13"),
               'Step s13 is not on the branch "Clusters"', fixed = TRUE)
  plan <- cellucid_apply_session(golden_sce(), golden(), branch = "Clusters", step = "s2")
  expect_identical(c(plan$branch, plan$step), c("b2", "s2"))
})

test_that("the options are exact", {
  skip_without_sce()
  for (option in c("dry_run", "overwrite", "keep_session")) {
    arguments <- list(golden_sce(), golden())
    arguments[[option]] <- 1
    expect_error(do.call(cellucid_apply_session, arguments), option)
  }
  expect_error(cellucid_apply_session(golden_sce(), golden(), step = 3), "step must be one string")
  expect_error(cellucid_apply_session(golden_sce(), golden(), nope = 1), "nope")
  expect_error(cellucid_apply_plan(list()), "plan must be what cellucid_apply_session")
})

test_that("a Seurat object and a data frame are applied to as a SingleCellExperiment is", {
  skip_without_seurat()
  skip_without_sce()
  sce <- golden_sce()
  counts <- SummarizedExperiment::assay(sce, "logcounts")
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  object$cell_type <- sce$cell_type
  object$leiden <- sce$leiden
  umap <- SingleCellExperiment::reducedDim(sce, "umap")
  dimnames(umap) <- list(colnames(object), c("UMAP_1", "UMAP_2"))
  object[["umap"]] <- SeuratObject::CreateDimReducObject(umap, key = "UMAP_", assay = "RNA")
  out <- cellucid_apply_session(object, golden(), dry_run = FALSE)
  meta <- methods::slot(out, "meta.data")
  expect_true(all(c("score", "cell_type (merged)", "cellucid_highlight__QC__Doublets") %in% names(meta)))
  expect_identical(levels(meta$cell_type), c("B cell", "T cell", "NK"))
  expect_identical(cellucid_history_applied(out)[[1L]]$step, "s13")
  expect_identical(methods::slot(out, "misc")$cellucid$history_applied[[1L]]$package, "r")
})

# =============================================================================
# Apply rules, one by one (synthesized files)
# =============================================================================

small_umap <- function(n = 6L) cbind(seq_len(n) - 1, (seq_len(n) - 1) * 2)

# A small SingleCellExperiment: categories, text and numbers in its colData,
# two genes, a 2D umap (obsm X_umap).
tiny_sce <- function(n = 6L) {
  cells <- paste0("cell_", seq_len(n) - 1L)
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = matrix(1, 2L, n, dimnames = list(c("G1", "G2"), cells))),
    colData = S4Vectors::DataFrame(
      cell_type = factor(c("B", "T", "T", "NK", "B", "T")[seq_len(n)], levels = c("B", "T", "NK")),
      leiden = factor(as.character((seq_len(n) - 1L) %% 3L)),
      label = paste0("x", (seq_len(n) - 1L) %% 2L),
      count = as.integer(seq_len(n) - 1L),
      row.names = cells
    )
  )
  SingleCellExperiment::reducedDim(object, "umap") <- small_umap(n)
  object
}

small_apply <- function(path, data = tiny_sce(), ..., dry_run = TRUE) cellucid_apply_session(data, path, ..., dry_run = dry_run)

small_fingerprint <- function(n = 6L, kept = NULL) {
  values <- small_umap(n)
  digest <- cellucid:::.apply_cell_order_digest(if (is.null(kept)) values else values[kept + 1, , drop = FALSE], kept)
  session_fingerprint(cell_count = if (is.null(kept)) n else length(kept), var_count = 2L, dataset_id = "exact", digest = digest)
}

small_session <- function(..., root = session_state(), n = 6L, content = "history") {
  path <- tempfile(fileext = ".cellucid-session")
  session_file(path, small_fingerprint(n), list(list(id = "t1", name = "Tracker 1", root = root, steps = list(...))), content = content)
}

test_that("labels rename columns and categories in place", {
  skip_without_sce()
  overlays <- session_overlays(
    renamed_fields = list(list(source = "obs", key = "leiden", name = "cluster"), list(source = "obs", key = "absent", name = "x"),
                          list(source = "var", key = "G1", name = "Gene one")),
    renamed_categories = list(
      list(source = "obs", field = "cell_type", category = "B", name = "B cell"),
      list(source = "obs", field = "cell_type", category = "NK", name = "T"),
      list(source = "obs", field = "label", category = "x1", name = "y"),
      list(source = "obs", field = "count", category = 1L, name = "one"),
      list(source = "obs", field = "leiden", category = 1L, name = "One")
    ),
    deleted_fields = list(list(source = "obs", key = "label", purged = FALSE))
  )
  plan <- small_apply(small_session(list(set = list(`field-overlays` = overlays))))
  rows <- plan$rows
  expect_identical(
    lapply(seq_len(nrow(rows)), function(at) c(rows$kind[[at]], rows$name[[at]], rows$column[[at]], paste(rows$from[[at]], collapse = ","),
                                                rows$status[[at]], rows$note[[at]])),
    list(
      c("column-renamed", "cluster", NA, "leiden", "write", NA),
      c("column-renamed", "x", NA, "absent", "not-applied", 'The object has no obs column "absent".'),
      c("column-renamed", "Gene one", NA, "G1", "not-applied", "A gene's display name is recorded, never applied to the object."),
      c("categories-renamed", "B cell", "cell_type", "B", "write", NA),
      c("categories-merged", "T", "cell_type", "NK", "write", NA),
      c("categories-renamed", "y", "label", "x1", "write", NA),
      c("categories-renamed", "one", "count", "1", "not-applied", 'The obs column "count" holds numbers, not categories.'),
      # leiden holds "1" (text), not the number 1: labels told apart as JavaScript tells them.
      c("categories-renamed", "One", "cluster", "1", "not-applied", 'The obs column "cluster" has no category "1".'),
      c("field-deleted", "label", NA, "", "not-applied", "A deletion is recorded, never applied to the object: the column stays.")
    )
  )
  expect_identical(rows$cells[startsWith(rows$kind, "categories-") & rows$status == "write"], c(2L, 1L, 3L))
  out <- cellucid_apply_plan(plan)
  expect_identical(names(SummarizedExperiment::colData(out)), c("cell_type", "cluster", "label", "count"))
  expect_identical(levels(out$cell_type), c("B cell", "T"))
  expect_identical(as.character(out$cell_type), c("B cell", "T", "T", "T", "B cell", "T"))
  expect_identical(out$label, c("x0", "y", "x0", "y", "x0", "y"))
})

test_that("a logical column's categories are false and true, as the viewer reads it", {
  skip_without_sce()
  overlays <- session_overlays(renamed_categories = list(
    list(source = "obs", field = "flag", category = TRUE, name = "kept")
  ))
  data <- tiny_sce()
  data$flag <- c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE)
  plan <- small_apply(small_session(list(set = list(`field-overlays` = overlays))), data = data)
  expect_identical(plan$rows$status, "write")
  expect_identical(plan$rows$from[[1]], "true")
  expect_identical(plan$rows$cells, 3L)
  out <- cellucid_apply_plan(plan)
  expect_identical(out$flag, c("kept", "false", "kept", "false", "false", "kept"))
})

test_that("swapped column names rename at once", {
  skip_without_sce()
  overlays <- session_overlays(renamed_fields = list(
    list(source = "obs", key = "leiden", name = "label"), list(source = "obs", key = "label", name = "leiden"),
    list(source = "obs", key = "count", name = "cell_type")
  ))
  out <- small_apply(small_session(list(set = list(`field-overlays` = overlays))), dry_run = FALSE)
  expect_identical(attr(out, "cellucid_apply_report")$plan$rows$status, c("write", "write", "not-applied"))
  expect_identical(out$label, tiny_sce()$leiden)
})

test_that("user fields apply by kind", {
  skip_without_sce()
  fields <- session_overlays(user_fields = list(
    session_user_category("user_cat_1", "annotation", c("a", "b"), c(0, 1, 255, 0, 1, 0)),
    session_user_values("user_obs_value_1", "pushed", c(0.5, NaN, 2), cells_of = c(0, 2, 4)),
    session_user_copy("user_obs_value_2", "copy of count", "obs", "count"),
    session_user_category("user_cat_2", "deleted one", "c", rep(0, 6), deleted = TRUE)
  ))
  plan <- small_apply(small_session(list(set = list(`field-overlays` = fields))))
  expect_identical(plan$rows$name, c("annotation", "pushed", "copy of count", "deleted one"))
  expect_identical(plan$rows$cells, c(5L, 2L, 6L, NA))
  expect_identical(plan$rows$status, c("write", "write", "write", "not-applied"))
  expect_identical(plan$rows$note, c(NA, NA, NA, "Deleted in the session."))
  out <- cellucid_apply_plan(plan)
  expect_identical(as.character(out$annotation), c("a", "b", NA, "a", "b", "a"))
  expect_true(all(is.na(out$pushed[c(2, 3, 4, 6)])))
  expect_identical(out$`copy of count`, as.double(0:5))
  expect_false("deleted one" %in% names(SummarizedExperiment::colData(out)))
})

test_that("a user category field's palette is read and checked", {
  skip_without_sce()
  # The colors a field made from highlight pages opens with (its pages'), one
  # #rrggbb a category, or null (F124).
  made <- session_user_category("user_cat_1", "annotation", c("a", "b"), c(0, 1, 255, 0, 1, 0))
  made$palette <- list("#dc2626", "#16a34a")
  plan <- small_apply(small_session(list(set = list(`field-overlays` = session_overlays(user_fields = list(made))))))
  expect_identical(plan$rows$status, "write")
  for (palette in list(list("#dc2626"), list("#DC2626", "#16a34a"), list("red", "#16a34a"), "#dc2626")) {
    bad <- made
    bad$palette <- palette
    expect_error(
      small_apply(small_session(list(set = list(`field-overlays` = session_overlays(user_fields = list(bad)))))),
      "palette must be null or one #rrggbb color a category"
    )
  }
})

test_that("a field of the intersections of 41 pages applies", {
  skip_without_sce()
  # Each intersection combines any number of pages (F269); the session keys
  # its label by its pages' positions ("0+1+...+40"), a record of labels; the
  # column is its codes.
  pages <- lapply(0:40, function(p) list(pageId = paste0("page_", p), label = paste0("P", p)))
  every <- paste(0:40, collapse = "+")
  categories <- c(paste0("P", 0:40), "First two", "Every page")
  for (labels in list(stats::setNames(list("First two", "Every page"), c("0+1", every)))) {
    made <- session_user_category("user_cat_1", "pages", categories, c(42, 41, 0, 1, 2, 40))
    made$sourcePages <- pages
    made$overlapStrategy <- "intersections"
    made$intersectionLabels <- labels
    plan <- small_apply(small_session(list(set = list(`field-overlays` = session_overlays(user_fields = list(made))))))
    expect_identical(plan$rows$status, "write")
    out <- cellucid_apply_plan(plan)
    expect_identical(as.character(out$pages), c("Every page", "First two", "P0", "P1", "P2", "P40"))
  }
})

test_that("a written column's colors are the ones it is shown in", {
  skip_without_sce()
  # A field's palette goes with its column as Scanpy keeps colors,
  # uns['<key>_colors'] (metadata(x) of a SingleCellExperiment), one a
  # category in the column's order -- every category, used or not (F139). A
  # field without a palette opens in Cellucid's own category colors, and is
  # written with them (F239); a column of no categories (a highlight) takes
  # no colors: those under its name described another column.
  made <- session_user_category("user_cat_1", "annotation", c("a", "b", "unused"), c(0, 1, 255, 0, 1, 0))
  made$palette <- list("#dc2626", "#16a34a", "#2563eb")
  plain <- session_user_category("user_cat_2", "plain", "x", rep(0, 6))
  groups <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", c(0, 1))))
  path <- small_session(list(set = list(`field-overlays` = session_overlays(user_fields = list(made, plain)), highlights = groups)))
  data <- tiny_sce()
  S4Vectors::metadata(data)$plain_colors <- c("#000000", "#111111", "#222222")
  S4Vectors::metadata(data)$cellucid_highlight__QC__Gate_colors <- "#000000"
  S4Vectors::metadata(data)$cell_type_colors <- c("#aaaaaa", "#bbbbbb", "#cccccc")
  plan <- small_apply(path, data)
  expect_null(attr(plan, "cellucid_context")$data@metadata$annotation_colors) # a dry run writes nothing
  out <- cellucid_apply_plan(plan)
  meta <- S4Vectors::metadata(out)
  expect_identical(levels(out$annotation), c("a", "b", "unused"))
  expect_identical(meta$annotation_colors, c("#dc2626", "#16a34a", "#2563eb"))
  expect_identical(meta$plain_colors, "#ff9900")
  expect_false("cellucid_highlight__QC__Gate_colors" %in% names(meta))
  expect_identical(meta$cell_type_colors, c("#aaaaaa", "#bbbbbb", "#cccccc")) # a column left as it is keeps its colors
  # Applied again: unchanged, colors as they are.
  S4Vectors::metadata(out)$annotation_colors <- c("#dc2626", "#16a34a", "#000000")
  again <- small_apply(path, out, dry_run = FALSE)
  expect_false(attr(again, "cellucid_apply_report")$applied)
  expect_identical(S4Vectors::metadata(again)$annotation_colors, c("#dc2626", "#16a34a", "#000000"))
  # Replaced (overwrite): the new palette replaces the colors.
  other <- made
  other$codes <- list(cells = NULL, values = typed("Uint8Array", rep(1, 6)))
  other$palette <- list("#111111", "#222222", "#333333")
  path <- small_session(list(set = list(`field-overlays` = session_overlays(user_fields = list(other, plain)), highlights = groups)))
  replaced <- small_apply(path, again, dry_run = FALSE, overwrite = TRUE)
  rows <- attr(replaced, "cellucid_apply_report")$plan$rows
  expect_identical(rows$status[rows$name == "annotation"], "replace")
  expect_identical(S4Vectors::metadata(replaced)$annotation_colors, c("#111111", "#222222", "#333333"))
})

test_that("colors follow a renamed column and merged categories", {
  skip_without_sce()
  overlays <- session_overlays(
    renamed_fields = list(list(source = "obs", key = "leiden", name = "cluster"), list(source = "obs", key = "count", name = "number")),
    renamed_categories = list(
      # B merges into T, in B's place: the category keeps B's color.
      list(source = "obs", field = "cell_type", category = "B", name = "T"),
      list(source = "obs", field = "label", category = "x1", name = "y")
    )
  )
  data <- tiny_sce()
  S4Vectors::metadata(data)$cell_type_colors <- c("#aa0000", "#00aa00", "#0000aa")
  S4Vectors::metadata(data)$leiden_colors <- c("#111111", "#222222", "#333333")
  S4Vectors::metadata(data)$cluster_colors <- "#999999" # stale: no column of that name
  S4Vectors::metadata(data)$count_colors <- matrix(0, 2, 4) # not Scanpy's form: not carried
  S4Vectors::metadata(data)$label_colors <- c("#444444", "#555555")
  out <- small_apply(small_session(list(set = list(`field-overlays` = overlays))), data, dry_run = FALSE)
  meta <- S4Vectors::metadata(out)
  expect_identical(levels(out$cell_type), c("T", "NK"))
  expect_identical(meta$cell_type_colors, c("#aa0000", "#0000aa"))
  expect_identical(meta$cluster_colors, c("#111111", "#222222", "#333333"))
  expect_false(any(c("leiden_colors", "count_colors", "number_colors") %in% names(meta)))
  # A text column relabeled: Scanpy colors it by its sorted labels, which the rename reorders.
  expect_false("label_colors" %in% names(meta))
})

test_that("merged categories drop colors that did not cover them", {
  skip_without_sce()
  overlays <- session_overlays(renamed_categories = list(
    list(source = "obs", field = "cell_type", category = "NK", name = "T"),
    list(source = "obs", field = "leiden", category = "0", name = "zero")
  ))
  data <- tiny_sce()
  S4Vectors::metadata(data)$cell_type_colors <- c("#aa0000", "#00aa00") # fewer than its 3 categories: Scanpy draws its own
  S4Vectors::metadata(data)$leiden_colors <- "#111111" # a rename alone: left as it is
  out <- small_apply(small_session(list(set = list(`field-overlays` = overlays))), data, dry_run = FALSE)
  expect_identical(levels(out$cell_type), c("B", "T"))
  expect_false("cell_type_colors" %in% names(S4Vectors::metadata(out)))
  expect_identical(S4Vectors::metadata(out)$leiden_colors, "#111111")
})

# One field's entry of a view's colors: its categories recolored.
recolor_entry <- function(field_id, ...) {
  pairs <- list(...)
  list(field = list(source = "obs", user = field_id),
       colors = lapply(pairs, function(pair) list(category = pair[[1L]], color = pair[[2L]])))
}

views_value <- function(focus, count) {
  view <- list(embedding = list(key = "X_umap", columns = list(0L, 1L)), followsDefault = TRUE,
               centroids = list(points = FALSE, labels = FALSE))
  list(layout = "grid", focus = focus, liveHidden = FALSE, camerasLinked = TRUE, views = rep(list(view), count))
}

test_that("a field's colors are the ones the view in focus shows", {
  skip_without_sce()
  # The colors a person or an agent gave a field's categories reach
  # <key>_colors (L26-7/8): in the view in focus at the step, each category
  # recolored there, the rest in the colors it opens with -- its palette,
  # else Cellucid's own colors by position. Another view's colors stay how
  # the data was shown.
  stage <- session_user_category("user_cat_1", "stage", c("early", "mid", "late"), c(0, 1, 2, 0, 1, 2))
  fate <- session_user_category("user_cat_2", "fate", c("a", "b"), c(0, 1, 0, 1, 0, 1))
  fate$palette <- list("#dc2626", "#16a34a")
  live <- list(recolor_entry("user_cat_1", list("early", "#111111")))
  kept <- list(
    list(field = list(source = "obs", key = "cell_type"), colors = list(list(category = "B", color = "#000000"))),
    recolor_entry("user_cat_1", list("late", "#008080"), list("early", "#808080"), list("absent", "#ffffff")),
    recolor_entry("user_cat_2", list("b", "#ffd700"))
  )
  made <- list(set = list(`field-overlays` = session_overlays(user_fields = list(stage, fate))))
  root <- session_state(colors = list(views = list(list(), list())), views = views_value(0L, 2L))
  recolored <- function(focus) list(set = list(colors = list(views = list(live, kept)), views = views_value(focus, 2L)))
  out <- small_apply(small_session(made, recolored(1L), root = root), dry_run = FALSE)
  meta <- S4Vectors::metadata(out)
  expect_identical(meta$stage_colors, c("#808080", "#0077c8", "#008080"))
  expect_identical(meta$fate_colors, c("#dc2626", "#ffd700"))
  # The dataset column recolored in the view in focus (RCC-1).
  expect_identical(meta$cell_type_colors, c("#000000", "#0077c8", "#ad50e0"))
  rows <- attr(out, "cellucid_apply_report")$plan$rows
  colors <- rows[rows$domain == "colors", ]
  expect_identical(colors$kind, c("column-colors", "presentation"))
  expect_identical(colors$name[[1L]], "cell_type")
  expect_identical(colors$status, c("write", "not-applied"))
  expect_identical(colors$note[[2L]], cellucid:::.apply_colors_note)
  expect_identical(rows$note[rows$domain == "views"], "Presentation: not applied to the object.")
  meta <- S4Vectors::metadata(small_apply(small_session(made, recolored(0L), root = root), dry_run = FALSE))
  expect_identical(meta$stage_colors, c("#111111", "#0077c8", "#ad50e0"))
  expect_identical(meta$fate_colors, c("#dc2626", "#16a34a"))
  expect_false("cell_type_colors" %in% names(meta)) # recolored in another view only
  # The colors at the step before (none recolored): the colors the fields open with.
  path <- small_session(made, recolored(1L), root = root)
  before <- small_apply(path, step = "s1", dry_run = FALSE)
  expect_identical(S4Vectors::metadata(before)$stage_colors, c("#ff9900", "#0077c8", "#ad50e0"))
  # At the head: the columns are unchanged (their colors as they are); the
  # dataset column recolored there takes its colors, once.
  again <- small_apply(path, before, dry_run = FALSE)
  rows <- attr(again, "cellucid_apply_report")$plan$rows
  expect_identical(rows$kind[rows$domain == "colors"], c("column-colors", "presentation"))
  expect_identical(rows$status[rows$domain == "colors"], c("write", "not-applied"))
  expect_identical(S4Vectors::metadata(again)$stage_colors, c("#ff9900", "#0077c8", "#ad50e0"))
  third <- small_apply(path, again, dry_run = FALSE)
  expect_false(attr(third, "cellucid_apply_report")$applied)
  rows <- attr(third, "cellucid_apply_report")$plan$rows
  expect_identical(rows$status[rows$domain == "colors"], c("unchanged", "not-applied"))
  expect_identical(rows$note[rows$domain == "colors" & rows$kind == "presentation"], cellucid:::.apply_colors_note)
})

# A dataset column's entry of a view's colors: its categories recolored.
dataset_recolor <- function(key, ...) {
  pairs <- list(...)
  list(field = list(source = "obs", key = key),
       colors = lapply(pairs, function(pair) list(category = pair[[1L]], color = pair[[2L]])))
}

test_that("a recolored dataset column is written with the colors it is shown in", {
  skip_without_sce()
  # A column that came with the dataset, recolored in the view in focus,
  # takes the colors it is shown in (RCC-1): the categories recolored, the
  # rest in the colors it opens with -- the colors held for it, as drawn
  # (#rgb widened, alpha dropped), else Cellucid's own by position. Applied
  # again: unchanged. A column it cannot be written to says why (a column the
  # object lacks is refused before, as every field the state uses).
  recolors <- list(
    dataset_recolor("cell_type", list("T", "#123456")),
    dataset_recolor("leiden", list("1", "#00ff00")),
    dataset_recolor("label", list("x0", "#ffffff")),
    dataset_recolor("count", list("absent", "#ffffff"))
  )
  path <- small_session(list(set = list(colors = list(views = list(recolors)), views = views_value(0L, 1L))),
                        root = session_state(colors = list(views = list(list())), views = views_value(0L, 1L)))
  data <- tiny_sce()
  S4Vectors::metadata(data)$cell_type_colors <- c("#AA0000", "#0a0", "#0000aaff")
  S4Vectors::metadata(data)$leiden_colors <- "#111111" # fewer than its 3 categories: the viewer opens it in Cellucid's colors
  plan <- small_apply(path, data)
  rows <- plan$rows[plan$rows$domain == "colors", ]
  expect_identical(rows$kind, c(rep("column-colors", 4L), "presentation"))
  expect_identical(rows$name, c("cell_type", "leiden", "label", "count", "Colors"))
  expect_identical(rows$status, c("write", "write", "not-applied", "not-applied", "not-applied"))
  expect_identical(rows$note, c(
    NA, NA,
    "The obs column \"label\" is not categorical; Scanpy keeps colors for a categorical column only.",
    "The obs column \"count\" is not categorical; Scanpy keeps colors for a categorical column only.",
    cellucid:::.apply_colors_note
  ))
  expect_identical(rows$column, rep(NA_character_, 5L))
  expect_identical(rows$cells, rep(NA_integer_, 5L))
  expect_identical(S4Vectors::metadata(attr(plan, "cellucid_context")$data)$cell_type_colors,
                   c("#AA0000", "#0a0", "#0000aaff")) # a dry run writes nothing
  out <- cellucid_apply_plan(plan)
  report <- attr(out, "cellucid_apply_report")
  expect_true(report$applied)
  meta <- S4Vectors::metadata(out)
  expect_identical(meta$cell_type_colors, c("#aa0000", "#123456", "#0000aa"))
  expect_identical(meta$leiden_colors, c("#ff9900", "#00ff00", "#ad50e0"))
  expect_false("label_colors" %in% names(meta))
  skipped <- Filter(function(entry) identical(entry$kind, "column-colors"), report$provenance$skipped)
  expect_identical(vapply(skipped, function(entry) entry$name, ""), c("label", "count"))
  again <- small_apply(path, out, dry_run = FALSE)
  expect_false(attr(again, "cellucid_apply_report")$applied)
  rows <- attr(again, "cellucid_apply_report")$plan$rows
  expect_identical(rows$status[rows$kind == "column-colors"], c("unchanged", "unchanged", "not-applied", "not-applied"))
})

test_that("a recolored dataset column renamed in the session takes its colors by its new labels", {
  skip_without_sce()
  # The colors value names a dataset column by the key the dataset gives it
  # and its categories by the labels shown at the step: the column renamed
  # and a category relabeled by the same apply are written under the new
  # names, the colors in the column's order; a label the column lacks is
  # said, the others written.
  overlays <- session_overlays(
    renamed_fields = list(list(source = "obs", key = "cell_type", name = "kind")),
    renamed_categories = list(list(source = "obs", field = "cell_type", category = "B", name = "Bcell"))
  )
  recolors <- list(dataset_recolor("cell_type", list("Bcell", "#101010"), list("ghost", "#202020")))
  root <- session_state(colors = list(views = list(list())), views = views_value(0L, 1L))
  path <- small_session(list(set = list(`field-overlays` = overlays, colors = list(views = list(recolors)), views = views_value(0L, 1L))),
                        root = root)
  data <- tiny_sce()
  S4Vectors::metadata(data)$cell_type_colors <- c("#aa0000", "#00aa00", "#0000aa")
  out <- small_apply(path, data, dry_run = FALSE)
  expect_identical(levels(SummarizedExperiment::colData(out)$kind), c("Bcell", "T", "NK"))
  meta <- S4Vectors::metadata(out)
  expect_identical(meta$kind_colors, c("#101010", "#00aa00", "#0000aa"))
  expect_false("cell_type_colors" %in% names(meta))
  rows <- attr(out, "cellucid_apply_report")$plan$rows
  row <- rows[rows$kind == "column-colors", ]
  expect_identical(c(row$name, row$status, row$note),
                   c("kind", "write", "Not written: the color of \"ghost\", which the obs column \"kind\" does not have."))
  # Only labels the column lacks: nothing to write.
  only <- list(dataset_recolor("leiden", list("9", "#202020"), list("8", "#303030")))
  plan <- small_apply(small_session(list(set = list(colors = list(views = list(only)), views = views_value(0L, 1L))), root = root))
  row <- plan$rows[plan$rows$kind == "column-colors", ]
  expect_identical(c(row$status, row$note),
                   c("not-applied", "Not written: the colors of \"9\", \"8\", which the obs column \"leiden\" does not have."))
})

test_that("colors the session holds are checked", {
  skip_without_sce()
  stage <- session_user_category("user_cat_1", "stage", "early", rep(0, 6))
  made <- list(set = list(`field-overlays` = session_overlays(user_fields = list(stage))))
  root <- session_state(colors = list(views = list(list())), views = views_value(0L, 1L))
  bad_focus <- views_value(0L, 1L)
  bad_focus$focus <- -1L
  cases <- list(
    list(colors = list(views = list(list(recolor_entry("user_cat_1", list("early", "#FFFFFF"))))), views = NULL,
         message = "color must be one #rrggbb color"),
    list(colors = list(views = list(list(list(field = list(source = "obs", user = ""))))), views = NULL,
         message = "field.user must be a non-empty text"),
    list(colors = list(views = list(list(list(field = list(source = "nope", key = "cell_type"))))), views = NULL,
         message = "field.source must be \"obs\" or \"var\""),
    list(colors = list(views = list(list())), views = views_value(1L, 2L),
         message = "colors value holds 1 view, and the view in focus is view 1"),
    list(colors = list(views = list(list())), views = bad_focus, message = "focus must be a non-negative whole number")
  )
  for (case in cases) {
    changes <- list(colors = case$colors)
    if (!is.null(case$views)) changes$views <- case$views
    expect_error(small_apply(small_session(made, list(set = changes), root = root)), case$message, fixed = TRUE)
  }
})

test_that("a session saved under the id an object was shown under applies when named", {
  skip_without_sce()
  # A viewer given dataset_id = serves the object under that id, and its
  # session is saved under it (L26-7): naming it applies the session to the
  # object, which declares its own id; another id is refused.
  record <- small_fingerprint()
  record$datasetId <- "pancreas-notebook"
  gate <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", c(0, 1))))
  path <- tempfile(fileext = ".cellucid-session")
  session_file(path, record, list(list(id = "t1", name = "Tracker 1", root = session_state(), steps = list(list(set = list(highlights = gate))))))
  data <- tiny_sce()
  S4Vectors::metadata(data)$cellucid <- list(id = "pancreas")
  expect_error(cellucid_apply_session(data, path),
               "saved on the dataset 'pancreas-notebook', and this object is 'pancreas' (its cellucid record's id). Apply it to the dataset it was saved on; when it is that dataset, shown under another id (dataset_id = when it was shown), name that id: dataset_id = 'pancreas-notebook'.",
               fixed = TRUE)
  expect_error(cellucid_apply_session(data, path, dataset_id = "pancreas"),
               "saved on the dataset 'pancreas-notebook', and this object was shown as 'pancreas'", fixed = TRUE)
  expect_error(cellucid_apply_session(data, path, dataset_id = "no/slash"), "dataset_id")
  out <- cellucid_apply_session(data, path, dataset_id = "pancreas-notebook", dry_run = FALSE)
  expect_true(attr(out, "cellucid_apply_report")$applied)
  expect_identical(as.logical(out$cellucid_highlight__QC__Gate), c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  # A plan made so is applied so (the plan is made again first).
  plan <- cellucid_apply_session(data, path, dataset_id = "pancreas-notebook")
  expect_true(attr(cellucid_apply_plan(plan), "cellucid_apply_report")$applied)
  # A pull of the session file takes the same id.
  expect_error(cellucid_pull(data, path), "this object is 'pancreas'")
  pulled <- cellucid_pull(data, path, dataset_id = "pancreas-notebook", return_summary = TRUE)
  expect_identical(names(pulled$summary$columns), "cellucid_highlight__highlight_1")
})

test_that("a viewer given another dataset id takes its own session back", {
  skip_without_sce()
  skip_if_not_installed("httpuv")
  # cellucid_view(object, dataset_id = ) serves the object under that id, and
  # the viewer saves its session under it; applied back through the viewer's
  # request, the session goes to the object it served, which declares its own
  # id (L26-7). An object declaring another id is refused, as before.
  data <- tiny_sce()
  S4Vectors::metadata(data)$cellucid <- list(id = "pancreas")
  viewer <- cellucid:::.start_viewer(
    data, list(dataset_id = "pancreas-notebook"), label = "data",
    port = 0L, port_given = TRUE, host = "127.0.0.1", allowed_hosts = NULL, serve_web_ui = FALSE,
    web_cache_dir = NULL, web_source_url = NULL, cache_bytes = 0, quiet = TRUE, bind = FALSE
  )
  on.exit(cellucid_stop(viewer), add = TRUE)
  expect_identical(c(viewer$dataset_id, viewer$declared_id), c("pancreas-notebook", "pancreas"))
  record <- small_fingerprint()
  record$datasetId <- "pancreas-notebook"
  stage <- session_user_category("user_cat_1", "stage", c("early", "late"), c(0, 1, 0, 1, 0, 1))
  answered <- function() {
    path <- tempfile(fileext = ".cellucid-session")
    session_file(path, record, list(list(id = "t1", name = "Tracker 1", root = session_state(),
                                         steps = list(list(set = list(`field-overlays` = session_overlays(user_fields = list(stage))))))))
    entry <- new.env(parent = emptyenv())
    entry$answer <- list(path = path)
    id <- paste0("s", length(viewer$session_requests) + 1L)
    viewer$session_requests[[id]] <- entry
    structure(list(viewer = viewer, id = id), class = "cellucid_session_request")
  }
  out <- cellucid_apply_session(data, answered(), dry_run = FALSE)
  expect_identical(as.character(out$stage), c("early", "late", "early", "late", "early", "late"))
  other <- data
  S4Vectors::metadata(other)$cellucid <- list(id = "garcia")
  expect_error(cellucid_apply_session(other, answered()), "saved on the dataset 'pancreas-notebook', and this object is 'garcia'")
})

test_that("a pull of a session file writes a field's palette under its column", {
  skip_without_sce()
  made <- session_user_category("user_cat_1", "annotation", c("a", "b"), c(0, 1, 255, 0, 1, 0))
  made$palette <- list("#dc2626", "#16a34a")
  path <- small_session(list(set = list(`field-overlays` = session_overlays(user_fields = list(made)))))
  data <- tiny_sce()
  S4Vectors::metadata(data)$field__annotation_colors <- "#000000"
  out <- cellucid_pull(data, path, field_prefix = "field__", return_summary = TRUE)
  expect_identical(names(out$summary$columns), "field__annotation")
  expect_identical(S4Vectors::metadata(out$data)$field__annotation_colors, c("#dc2626", "#16a34a"))
  # Written beside a column of that name (on_conflict = "suffix"): under the name written.
  data <- tiny_sce()
  data$annotation <- 1
  S4Vectors::metadata(data)$annotation_colors <- "#000000"
  out <- cellucid_pull(data, path, on_conflict = "suffix", return_summary = TRUE)
  written <- names(out$summary$columns)
  expect_length(written, 1L)
  expect_false(identical(written, "annotation"))
  expect_identical(S4Vectors::metadata(out$data)[[paste0(written, "_colors")]], c("#dc2626", "#16a34a"))
  expect_identical(S4Vectors::metadata(out$data)$annotation_colors, "#000000")
})

test_that("a column's colors live where the object keeps uns", {
  # Seurat misc, metadata(x), a data frame's attributes: the cellucid
  # record's places, under uns['<column>_colors']'s own name.
  frame <- data.frame(a = factor(c("x", "y")), row.names = c("c1", "c2"))
  attr(frame, "a_colors") <- c("#000000", "#111111")
  out <- cellucid:::.object_with_columns(frame, list(list(name = "a", values = factor(c("y", "y")), colors = c("#123456", "#abcdef"))))
  expect_identical(attr(out, "a_colors"), c("#123456", "#abcdef"))
  out <- cellucid:::.object_with_columns(out, list(), remove = "a")
  expect_null(attr(out, "a_colors"))
  # The cellucid record is read by its exact name: a highlight column's
  # colors (cellucid_highlight__..._colors) are not it.
  skip_without_sce()
  sce <- tiny_sce()
  S4Vectors::metadata(sce)$cellucid_highlight__QC__Gate_colors <- c("#000000", "#111111")
  expect_identical(cellucid:::.object_record(sce), list())
})

test_that("highlight groups of every page are named by page and label", {
  skip_without_sce()
  groups <- session_highlights(
    session_page("page_1", "QC / pass", session_group("highlight_1", "Gate", c(0, 1)), session_group("highlight_2", "Gate!", 2)),
    session_page("page_2", "Types", session_group("highlight_3", "cell_type: T", c(1, 2, 5), kind = "category",
                                                  field = list(source = "obs", key = "cell_type"), category = "T"))
  )
  out <- small_apply(small_session(list(set = list(highlights = groups))), dry_run = FALSE)
  plan <- attr(out, "cellucid_apply_report")$plan
  highlight <- plan$rows[plan$rows$kind == "highlight", ]
  expect_identical(highlight$name, c("cellucid_highlight__QC_pass__Gate", "cellucid_highlight__QC_pass__Gate_2",
                                     "cellucid_highlight__Types__cell_type_T"))
  expect_identical(highlight$cells, c(2L, 1L, 3L))
  expect_type(out$cellucid_highlight__QC_pass__Gate, "logical")
  expect_identical(plan$requires, list(obs = "cell_type", var = character()))
})

test_that("a subsample knows only its cells", {
  skip_without_sce()
  n <- 6L
  subsample <- session_subsample_record(n, 4L, seed = 3L)
  kept <- as.double(subsample$kept)
  groups <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", kept[[1L]])))
  path <- session_file(tempfile(fileext = ".cellucid-session"), small_fingerprint(n, kept),
                       list(list(id = "t1", name = "Tracker 1", root = session_state(), steps = list(list(set = list(highlights = groups))))),
                       cells_record = subsample$record)
  session <- cellucid_read_session(path)
  expect_identical(as.double(session$cells$kept), kept)
  out <- small_apply(session, dry_run = FALSE)
  column <- out$cellucid_highlight__QC__Gate
  unknown <- setdiff(seq_len(n) - 1, kept)
  expect_true(all(is.na(column[unknown + 1])))
  expect_false(anyNA(column[kept + 1]))
  expect_error(small_apply(session, data = tiny_sce(5L)), "saved on 4 of 6 cells, and this object has 5")
})

test_that("an earlier apply's column is replaced only with overwrite", {
  skip_without_sce()
  one <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", c(0, 1))))
  two <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", c(2, 3, 4))))
  path <- small_session(list(set = list(highlights = one)), list(set = list(highlights = two)))
  first <- small_apply(path, step = "s1", dry_run = FALSE)
  plan <- small_apply(path, data = first, step = "s2")
  row <- plan$rows[plan$rows$kind == "highlight", ]
  expect_identical(c(row$status, row$note), c(
    "not-applied", 'An earlier apply (step s1) wrote "cellucid_highlight__QC__Gate" with other values; overwrite replaces it.'
  ))
  replaced <- small_apply(path, data = first, step = "s2", overwrite = TRUE, dry_run = FALSE)
  expect_identical(attr(replaced, "cellucid_apply_report")$plan$rows$status, "replace")
  expect_identical(replaced$cellucid_highlight__QC__Gate, c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(vapply(cellucid_history_applied(replaced), function(record) record$step, ""), c("s1", "s2"))
})

test_that("gene sets are written once and never over the person's", {
  skip_without_sce()
  sets <- list(sets = list(list(name = "Markers", description = "d", genes = list("G1")),
                           list(name = "Mine", description = "", genes = list("G2"))))
  data <- tiny_sce()
  S4Vectors::metadata(data)$cellucid <- list(gene_sets = list(Mine = list(genes = c("G1", "G2"), description = "", origin = "user")))
  out <- small_apply(small_session(list(set = list(`gene-sets` = sets))), data = data, dry_run = FALSE)
  rows <- attr(out, "cellucid_apply_report")$plan$rows
  expect_identical(rows$status, c("write", "not-applied"))
  expect_identical(rows$note[[2L]], 'The object already has a gene set "Mine" that no session wrote; it is never replaced.')
  record <- S4Vectors::metadata(out)$cellucid
  expect_identical(record$gene_sets$Mine$genes, c("G1", "G2"))
  expect_identical(record$gene_sets$Markers$origin, "session")
})

# A session whose steps hold the gene sets `...` (each a list of sets), saved
# on `mapping` over `var_count` genes.
gene_set_session <- function(..., mapping = session_mapping(), var_count = 2L) {
  fingerprint <- small_fingerprint()
  fingerprint$varCount <- var_count
  fingerprint$dataMapping <- mapping
  steps <- lapply(list(...), function(sets) list(set = list(`gene-sets` = list(sets = sets))))
  session_file(tempfile(fileext = ".cellucid-session"), fingerprint,
               list(list(id = "t1", name = "Tracker 1", root = session_state(), steps = steps)))
}
gene_set <- function(name, ...) list(name = name, description = "", genes = list(...))

# The value of `expr`, and the warnings it gave.
with_warnings <- function(expr) {
  warned <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warned)
}

test_that("gene sets are columns of the var the data mapping reads", {
  skip_without_sce()
  data <- tiny_sce()
  SummarizedExperiment::rowData(data)$symbol <- c("g1", "g2")
  path <- gene_set_session(list(
    gene_set("T cell markers", "G2"), gene_set("T-cell markers", "G1", "G2"), gene_set("!!!", "G1"),
    gene_set("Wide", "G1", "X1", "X2", "X3", "X4"), gene_set("Few", "G2", "Y1")
  ))
  out <- with_warnings(small_apply(path, data = data, dry_run = FALSE))
  frame <- SummarizedExperiment::rowData(out$value)
  expect_identical(names(frame), c("symbol", "cellucid_gene_set__T_cell_markers", "cellucid_gene_set__T_cell_markers_2",
                                   "cellucid_gene_set__gene_set", "cellucid_gene_set__Wide", "cellucid_gene_set__Few"))
  expect_identical(frame$cellucid_gene_set__T_cell_markers, c(FALSE, TRUE))
  expect_identical(frame$cellucid_gene_set__T_cell_markers_2, c(TRUE, TRUE))
  expect_identical(frame$cellucid_gene_set__Wide, c(TRUE, FALSE))
  expect_identical(out$warnings, c(
    'The gene set "Wide" names 4 genes the dataset does not have (X1, X2, X3, …), left out of its var column.',
    'The gene set "Few" names 1 gene the dataset does not have (Y1), left out of its var column.'
  ))
  # The plan's rows are the shared model's: no var column in them, none in the provenance.
  report <- attr(out$value, "cellucid_apply_report")
  expect_identical(report$plan$rows$status, rep("write", 5L))
  expect_length(report$provenance$columns, 0L)
  # The mapping names genes by a var column: the column's names mark the rows.
  named <- gene_set_session(list(gene_set("By symbol", "g2")), mapping = session_mapping(geneNames = "column:symbol"))
  frame <- SummarizedExperiment::rowData(small_apply(named, data = data, dry_run = FALSE))
  expect_identical(frame$cellucid_gene_set__By_symbol, c(FALSE, TRUE))
})

test_that("a gene set an earlier apply wrote keeps its var column's place when replaced", {
  skip_without_sce()
  path <- gene_set_session(list(gene_set("Markers", "G1")), list(gene_set("Markers", "G2")))
  first <- small_apply(path, step = "s1", dry_run = FALSE)
  SummarizedExperiment::rowData(first)$later <- c(1, 2)
  kept <- small_apply(path, data = first, step = "s2", dry_run = FALSE)
  expect_identical(attr(kept, "cellucid_apply_report")$plan$rows$status, "not-applied")
  expect_identical(SummarizedExperiment::rowData(kept)$cellucid_gene_set__Markers, c(TRUE, FALSE))
  replaced <- small_apply(path, data = first, step = "s2", overwrite = TRUE, dry_run = FALSE)
  expect_identical(attr(replaced, "cellucid_apply_report")$plan$rows$status, "replace")
  frame <- SummarizedExperiment::rowData(replaced)
  expect_identical(names(frame), c("cellucid_gene_set__Markers", "later"))
  expect_identical(frame$cellucid_gene_set__Markers, c(FALSE, TRUE))
  expect_identical(S4Vectors::metadata(replaced)$cellucid$gene_sets$Markers$genes, "G2")
  # A set the step no longer holds has no row: its record and its column stay.
  gone <- gene_set_session(list(gene_set("Markers", "G1")), list(gene_set("Other", "G2")))
  later <- small_apply(gone, data = small_apply(gone, step = "s1", dry_run = FALSE), step = "s2", dry_run = FALSE)
  expect_identical(names(SummarizedExperiment::rowData(later)), c("cellucid_gene_set__Markers", "cellucid_gene_set__Other"))
  expect_identical(names(S4Vectors::metadata(later)$cellucid$gene_sets), c("Markers", "Other"))
})

test_that("a var column the dataset holds is its own, never a gene set's", {
  skip_without_sce()
  data <- tiny_sce()
  SummarizedExperiment::rowData(data)$cellucid_gene_set__Markers <- c("mine", "mine")
  out <- with_warnings(small_apply(gene_set_session(list(gene_set("Markers", "G1"))), data = data, dry_run = FALSE))
  expect_identical(attr(out$value, "cellucid_apply_report")$plan$rows$status, "write")
  expect_identical(SummarizedExperiment::rowData(out$value)$cellucid_gene_set__Markers, c("mine", "mine"))
  expect_identical(S4Vectors::metadata(out$value)$cellucid$gene_sets$Markers$genes, "G1")
  expect_identical(out$warnings,
                   'The gene set "Markers" has no var column: the var column "cellucid_gene_set__Markers" is the dataset\'s own.')
})

test_that("a mapping that reads no genes keeps the gene sets in the record only", {
  skip_without_sce()
  path <- gene_set_session(list(gene_set("Markers", "G1")), mapping = session_mapping(expression = NULL), var_count = 0L)
  out <- with_warnings(small_apply(path, dry_run = FALSE))
  expect_identical(names(SummarizedExperiment::rowData(out$value)), character())
  expect_identical(S4Vectors::metadata(out$value)$cellucid$gene_sets$Markers$genes, "G1")
  expect_identical(out$warnings, paste0(
    "The gene sets are kept in the cellucid record's gene_sets only: the data mapping reads its genes from another ",
    "matrix than var, so they have no var columns."
  ))
  # Nothing written, nothing said: applying again finds the set unchanged.
  again <- with_warnings(small_apply(path, data = out$value, dry_run = FALSE))
  expect_identical(again$warnings, character())
})

test_that("a gene set's column goes to the modality the mapping reads", {
  skip_without_sce()
  data <- tiny_sce()
  SingleCellExperiment::altExp(data, "ADT") <- SummarizedExperiment::SummarizedExperiment(
    list(counts = matrix(1, 3L, 6L, dimnames = list(c("P1", "P2", "P3"), colnames(data))))
  )
  path <- gene_set_session(list(gene_set("Proteins", "P2", "P3")), mapping = session_mapping(modality = "ADT"), var_count = 3L)
  out <- small_apply(path, data = data, dry_run = FALSE)
  expect_identical(SummarizedExperiment::rowData(SingleCellExperiment::altExp(out, "ADT"))$cellucid_gene_set__Proteins,
                   c(FALSE, TRUE, TRUE))
  expect_identical(names(SummarizedExperiment::rowData(out)), character())
  main <- gene_set_session(list(gene_set("Genes", "G1")), mapping = session_mapping(modality = "main"))
  out <- small_apply(main, data = data, dry_run = FALSE)
  expect_identical(SummarizedExperiment::rowData(out)$cellucid_gene_set__Genes, c(TRUE, FALSE))
  expect_identical(names(SummarizedExperiment::rowData(SingleCellExperiment::altExp(out, "ADT"))), character())
})

test_that("a Seurat assay's feature metadata takes a gene set's column, an Assay and an Assay5 alike", {
  skip_without_seurat()
  skip_without_sce()
  sce <- tiny_sce()
  rna <- Matrix::Matrix(SummarizedExperiment::assay(sce, "logcounts"), sparse = TRUE)
  adt <- Matrix::Matrix(matrix(1, 3L, 6L, dimnames = list(c("P1", "P2", "P3"), colnames(sce))), sparse = TRUE)
  creators <- list(Assay = SeuratObject::CreateAssayObject)
  if (seurat_v5()) creators$Assay5 <- SeuratObject::CreateAssay5Object
  for (class in names(creators)) {
    object <- suppressWarnings(SeuratObject::CreateSeuratObject(creators[[class]](counts = rna)))
    suppressWarnings(object[["ADT"]] <- creators[[class]](counts = adt))
    umap <- small_umap()
    dimnames(umap) <- list(colnames(object), c("UMAP_1", "UMAP_2"))
    object[["umap"]] <- SeuratObject::CreateDimReducObject(umap, key = "UMAP_", assay = "RNA")
    object[["ADT"]][["existing"]] <- stats::setNames(c(1, 2, 3), c("P1", "P2", "P3"))
    path <- gene_set_session(list(gene_set("Proteins", "P3")), list(gene_set("Proteins", "P1")),
                             mapping = session_mapping(modality = "ADT"), var_count = 3L)
    first <- small_apply(path, data = object, step = "s1", dry_run = FALSE)
    expect_s4_class(first[["ADT"]], class)
    expect_identical(names(first[["ADT"]][[]]), c("existing", "cellucid_gene_set__Proteins"), info = class)
    expect_identical(first[["ADT"]][[]]$cellucid_gene_set__Proteins, c(FALSE, FALSE, TRUE), info = class)
    expect_false("cellucid_gene_set__Proteins" %in% names(first[["RNA"]][[]]))
    replaced <- small_apply(path, data = first, step = "s2", overwrite = TRUE, dry_run = FALSE)
    expect_identical(names(replaced[["ADT"]][[]]), c("existing", "cellucid_gene_set__Proteins"), info = class)
    expect_identical(replaced[["ADT"]][[]]$cellucid_gene_set__Proteins, c(TRUE, FALSE, FALSE), info = class)
  }
})

test_that("an apply never writes, renames or replaces the priority-order column", {
  skip_without_sce()
  message <- paste0("obs['cellucid_original_position'] holds each row's original position in a file stored in priority order, ",
                    "so an apply never writes, renames or replaces it.")
  field <- session_overlays(user_fields = list(session_user_values("f1", "cellucid_original_position", c(0, 1, 2, 3, 4, 5))))
  expect_error(small_apply(small_session(list(set = list(`field-overlays` = field)))), message, fixed = TRUE)
  onto <- session_overlays(renamed_fields = list(list(source = "obs", key = "leiden", name = "cellucid_original_position")))
  expect_error(small_apply(small_session(list(set = list(`field-overlays` = onto)))), message, fixed = TRUE)
  # Its own column (no priority record beside it, so the cells are in their order) is not renamed away either.
  data <- tiny_sce()
  data$cellucid_original_position <- as.double(seq_len(ncol(data)) - 1L)
  away <- session_overlays(renamed_fields = list(list(source = "obs", key = "cellucid_original_position", name = "position")))
  expect_error(suppressWarnings(small_apply(small_session(list(set = list(`field-overlays` = away))), data = data, dry_run = FALSE)),
               message, fixed = TRUE)
})

test_that("presentation domains are reported, never applied", {
  skip_without_sce()
  camera <- list(views = list(list(zoom = 2, ratio = 1)))
  plan <- small_apply(small_session(list(set = list(camera = camera))))
  expect_identical(unlist(plan$rows[, c("kind", "domain", "name", "status", "note")], use.names = FALSE),
                   c("presentation", "camera", "Camera", "not-applied", "Presentation: not applied to the object."))
  expect_false(attr(cellucid_apply_plan(plan), "cellucid_apply_report")$applied)
})

test_that("the step's expression values are recorded in its apply's record, never applied", {
  skip_without_sce()
  gate <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", c(0, 1))))
  root <- session_state(expression = list(normalize = FALSE, log1p = FALSE))
  shown <- list(normalize = TRUE, log1p = TRUE)
  path <- small_session(list(set = list(expression = shown)), list(set = list(highlights = gate)), root = root)
  data <- tiny_sce()
  applied <- small_apply(path, data = data, dry_run = FALSE)
  report <- attr(applied, "cellucid_apply_report")
  expect_true(report$applied)
  presentation <- report$plan$rows[report$plan$rows$kind == "presentation", c("kind", "domain", "name", "status", "note")]
  expect_identical(unlist(presentation, use.names = FALSE),
                   c("presentation", "expression", "Expression values", "not-applied", "Presentation: not applied to the object."))
  expect_identical(report$provenance$expression, shown)
  expect_identical(lapply(cellucid_history_applied(applied), function(record) record$expression), list(shown))
  # Recorded, never applied: the object's expression values are as stored.
  expect_identical(SummarizedExperiment::assay(applied, "logcounts"), SummarizedExperiment::assay(data, "logcounts"))
  # The step before (only the expression values changed) writes nothing, so records nothing.
  expect_false(attr(small_apply(path, step = "s1", dry_run = FALSE), "cellucid_apply_report")$applied)
  # The opening's values (as stored) are recorded as they are.
  as_stored <- small_apply(small_session(list(set = list(highlights = gate)), root = root), dry_run = FALSE)
  expect_identical(cellucid_history_applied(as_stored)[[1L]]$expression, list(normalize = FALSE, log1p = FALSE))
})

test_that("a session without expression values records NULL, and another shape is refused", {
  skip_without_sce()
  gate <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", c(0, 1))))
  applied <- small_apply(small_session(list(set = list(highlights = gate))), dry_run = FALSE)
  record <- cellucid_history_applied(applied)[[1L]]
  expect_true("expression" %in% names(record))
  expect_null(record$expression)
  expect_identical(as_json_value(record["expression"]), list(expression = NULL))
  root <- session_state(expression = list(normalize = FALSE, log1p = FALSE))
  expect_error(small_apply(small_session(list(set = list(expression = list(normalize = TRUE), highlights = gate)), root = root)),
               "expression value at step s1 has other fields than this version reads (missing log1p)", fixed = TRUE)
  expect_error(small_apply(small_session(list(set = list(expression = list(normalize = TRUE, log1p = 1), highlights = gate)), root = root)),
               "expression value at step s1's log1p must be true or false.", fixed = TRUE)
})

test_that("a steps-only file reads and is not applied", {
  skip_without_sce()
  path <- small_session(list(set = list(camera = list(views = list(list(zoom = 2, ratio = 1))))), content = "steps")
  session <- cellucid_read_session(path)
  expect_identical(c(session$scope$content, session$scope$results), c("steps", "recompute"))
  expect_identical(cellucid_steps(cellucid_tracker(session))$id, c("s0", "s1"))
  expect_error(cellucid_step_state(cellucid_tracker(session), "s1"), "shared as its steps only")
  expect_error(small_apply(session), "shared as its steps only, without the values they leave behind, so it cannot be applied")
})

test_that("a state shared without its history names nobody", {
  groups <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", c(0, 1))))
  shared <- list(set = list(highlights = groups), label = "Shared state", agent = "Claude", squashed = list(cause = "share", count = NULL, by = list()))
  expect_error(cellucid_read_session(small_session(shared, content = "state")),
               "session.scope of a state names nobody: every by is the person or the opening", fixed = TRUE)
})

test_that("a state-only file reads and applies as any file", {
  skip_without_sce()
  groups <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", c(0, 1))))
  shared <- list(set = list(highlights = groups), label = "Shared state", squashed = list(cause = "share", count = NULL, by = list()))
  session <- cellucid_read_session(small_session(shared, content = "state"))
  expect_identical(list(session$scope$content, session$scope$notes, session$scope$chats, session$scope$journal),
                   list("state", FALSE, FALSE, FALSE))
  tracker <- cellucid_tracker(session)
  steps <- cellucid_steps(tracker)
  expect_identical(steps$kind, c("root", "squashed"))
  expect_identical(steps$type, c(NA_character_, NA_character_))
  expect_identical(steps$v, c(NA_integer_, NA_integer_))
  expect_identical(steps$label, c("Opened", "Shared state"))
  expect_identical(cellucid_step(tracker, "s1")$squashed$cause, "share")
  out <- small_apply(session, dry_run = FALSE)
  report <- attr(out, "cellucid_apply_report")
  expect_true(report$applied)
  expect_identical(report$plan$rows$name[report$plan$rows$kind == "highlight"], "cellucid_highlight__QC__Gate")
  expect_identical(out$cellucid_highlight__QC__Gate, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("a state file holds nothing but the state", {
  groups <- session_highlights(session_page("page_1", "QC", session_group("highlight_1", "Gate", c(0, 1))))
  path <- small_session(list(set = list(highlights = groups)), list(set = list(camera = list(views = list(list(zoom = 2, ratio = 1))))),
                        content = "state")
  expect_error(cellucid_read_session(path),
               "session.scope of a state holds the root and one squashed step (cause share), shown", fixed = TRUE)
})

# =============================================================================
# Damaged and newer files
# =============================================================================

golden_bytes <- function() readBin(golden_path(), "raw", file.size(golden_path()))

# The golden file with `edit(file)` applied to its manifest and chunks
# (list(manifest, chunks)), written again.
rewritten <- function(edit) {
  file <- edit(session_unframe(golden_bytes()))
  path <- tempfile(fileext = ".cellucid-session")
  writeBin(session_frame(file$manifest, file$chunks), path)
  path
}

chunk_index <- function(file, id) which(vapply(file$manifest$chunks, function(entry) entry$id, "") == id)

restore_sizes <- function(file) {
  for (at in seq_along(file$chunks)) {
    entry <- file$manifest$chunks[[at]]
    file$manifest$chunks[[at]]$storedBytes <- length(file$chunks[[at]])
    file$manifest$chunks[[at]]$uncompressedBytes <- if (identical(entry$codec, "gzip")) {
      length(cellucid:::.native_inflate(file$chunks[[at]], TRUE))
    } else {
      length(file$chunks[[at]])
    }
  }
  file
}

edit_json_chunk <- function(file, id, edit) {
  at <- chunk_index(file, id)
  entry <- file$manifest$chunks[[at]]
  payload <- jsonlite::parse_json(rawToChar(session_payload(entry, file$chunks[[at]])), simplifyVector = FALSE)
  stored <- session_store(entry[setdiff(names(entry), c("codec", "storedBytes", "uncompressedBytes"))], session_json(edit(payload)), entry$codec)
  file$manifest$chunks[[at]] <- stored$meta
  file$chunks[[at]] <- stored$bytes
  file
}

step_edit <- function(edit) function(file) edit_json_chunk(file, "steps/0/0", function(steps) {
  at <- which(vapply(steps, function(step) step$id, "") == "s1")
  steps[[at]] <- edit(steps[[at]])
  steps
})

manifest_edit <- function(edit) function(file) {
  file$manifest <- edit(file$manifest)
  file
}

damages <- list(
  `cut gzip` = list(function(file) {
    at <- chunk_index(file, "values/eager/0")
    file$chunks[[at]] <- file$chunks[[at]][seq_len(length(file$chunks[[at]]) %/% 2L)]
    file$manifest$chunks[[at]]$storedBytes <- length(file$chunks[[at]])
    file
  }, 'the chunk "values/eager/0"'),
  `flipped byte` = list(function(file) {
    at <- chunk_index(file, "values/eager/0")
    middle <- length(file$chunks[[at]]) %/% 2L
    file$chunks[[at]][[middle]] <- xor(file$chunks[[at]][[middle]], as.raw(0xff))
    file
  }, 'the chunk "values/eager/0"'),
  `wrong hash` = list(function(file) edit_json_chunk(file, "values/eager/0", function(values) {
    for (hash in names(values)) if (is.list(values[[hash]]) && length(values[[hash]]$pages) > 0L) values[[hash]]$pages[[1L]]$name <- "Page 9"
    values
  }), "its eager values holds a value that does not match its hash."),
  `unreferenced value` = list(function(file) edit_json_chunk(file, "values/lazy/0", function(values) {
    values[[cellucid:::.value_hash(list(unused = TRUE))]] <- list(unused = TRUE)
    values
  }), "it holds a value no step refers to."),
  `missing value` = list(function(file) edit_json_chunk(file, "values/lazy/0", function(values) values[-1L]),
                         "a step refers to a value the file does not hold."),
  `picture chunk` = list(manifest_edit(function(m) {
    m$chunks[[length(m$chunks) + 1L]] <- c(m$chunks[[length(m$chunks)]][c("id", "priority", "codec", "storedBytes", "uncompressedBytes")],
                                           list(kind = "picture", hash = strrep("a", 64L)))
    m
  }), "kind must be one of steps, parts, array, values, journal, chat."),
  `picture field` = list(step_edit(function(step) {
    step["picture"] <- list(NULL)
    step
  }), 'has an unknown field "picture"'),
  `no version` = list(step_edit(function(step) {
    step["v"] <- list(NULL)
    step
  }), "records the version of its action's parameters (v)"),
  `split turn` = list(function(file) edit_json_chunk(file, "steps/0/0", function(steps) {
    at <- which(vapply(steps, function(step) step$id, "") == "s3")
    steps[[at]]["group"] <- list(NULL)
    steps
  }), 'a group cuts the run of turn "py-1" between steps "s2" and "s3"; a group holds whole runs of a turn'),
  `squashed action` = list(function(file) edit_json_chunk(file, "steps/0/0", function(steps) {
    at <- which(vapply(steps, function(step) step$id, "") == "s14")
    steps[[at]]$type <- "color.by"
    steps
  }), "is squashed and records no action of its own"),
  `paused by nobody` = list(manifest_edit(function(m) {
    m$trackers[[2L]]$paused$by <- list()
    m
  }), "tracker.paused names who paused it (by)"),
  `stamps out of order` = list(manifest_edit(function(m) {
    m$actions[1:2] <- m$actions[2:1]
    m
  }), "session.actions must be sorted by type then version, each pair once"),
  `stamp unused` = list(manifest_edit(function(m) {
    m$actions[[1L]]$version <- 2L
    m
  }), "session.actions must stamp exactly each action and version its steps use"),
  `paused draft` = list(manifest_edit(function(m) {
    m$trackers[[2L]]$paused$draft <- list(camera = strrep("a", 64L))
    m
  }), "keeps a paused draft; a file squashes it into a step first"),
  `state of a history` = list(manifest_edit(function(m) {
    m$scope$content <- "state"
    m
  }), "session.scope of a state holds one tracker of one branch"),
  `head off its branch` = list(manifest_edit(function(m) {
    m$trackers[[1L]]$headId <- "s8"
    m
  }), "head must be on its active branch's path."),
  `duplicate branch names` = list(manifest_edit(function(m) {
    m$trackers[[1L]]$branches[[2L]]$name <- "MAIN"
    m
  }), 'branch "b0" and branch "b1" share the name "MAIN"'),
  `id-like name` = list(manifest_edit(function(m) {
    m$trackers[[1L]]$branches[[2L]]$name <- "B7"
    m
  }), 'A branch cannot be called "B7": that is the form of a branch id.'),
  tree = list(step_edit(function(step) {
    step$children <- c(step$children, list("s5"))
    step
  }), "lists a child that is not its own"),
  touches = list(step_edit(function(step) {
    step$touches <- c(step$touches, list("camera"))
    step
  }), "must hold a value before and after for exactly the domains it touches"),
  origin = list(step_edit(function(step) {
    step$by$origin <- "robot"
    step
  }), 'by.origin must be one of "user", "agent", "host", "session", "replay", "startup", "outside"'),
  `old step key` = list(step_edit(function(step) {
    step["note"] <- list(NULL)
    step
  }), 'has an unknown field "note"'),
  `branch membership` = list(step_edit(function(step) {
    step$branch <- "b3"
    step
  }), "steps must form one line starting right after its fork step"),
  active = list(manifest_edit(function(m) {
    m$activeTracker <- "t-0"
    m
  }), "its active tracker is not one it holds."),
  `notes left out` = list(manifest_edit(function(m) {
    m$scope$notes <- FALSE
    m
  }), "holds notes or bookmarks its scope leaves out (notes)."),
  `dates only` = list(manifest_edit(function(m) {
    m$scope$times <- "dates"
    m
  }), "holds a time of day at"),
  `source name` = list(manifest_edit(function(m) {
    m$scope$source <- "name"
    m
  }), "shared with the dataset's name only"),
  `scope trackers` = list(manifest_edit(function(m) {
    m$scope$trackers <- rev(m$scope$trackers)
    m
  }), "session.scope must name exactly the trackers it holds, in order"),
  order = list(manifest_edit(function(m) {
    last <- length(m$chunks)
    m$chunks <- c(m$chunks[last], m$chunks[-last])
    m
  }), 'its chunk "steps/0/0" is out of order (eager before lazy; steps, parts, array, values, journal, chat).'),
  requires = list(manifest_edit(function(m) {
    m$requires$obs <- list("z", "a")
    m
  }), "must be sorted and without duplicates."),
  `chrome missing` = list(manifest_edit(function(m) {
    m$chrome <- NULL
    m
  }), "missing chrome"),
  `chrome not a record` = list(manifest_edit(function(m) {
    m$chrome <- list()
    m
  }), "its page chrome must be exactly {pressed}."),
  `chrome key` = list(manifest_edit(function(m) {
    m$chrome$extra <- 1L
    m
  }), "its page chrome must be exactly {pressed}."),
  `chrome pressed` = list(manifest_edit(function(m) {
    m$chrome$pressed <- list(1L)
    m
  }), "its page chrome must map pressed-button groups to their pressed button."),
  `chrome empty group` = list(manifest_edit(function(m) {
    m$chrome$pressed <- c(m$chrome$pressed, stats::setNames(list("lasso"), ""))
    m
  }), 'its page chrome must name the pressed button of the group "" by its id.'),
  `chrome empty button` = list(manifest_edit(function(m) {
    m$chrome$pressed[["highlight-mode"]] <- ""
    m
  }), 'its page chrome must name the pressed button of the group "highlight-mode" by its id.'),
  `chrome not text` = list(manifest_edit(function(m) {
    m$chrome$pressed[["highlight-mode"]] <- 3L
    m
  }), 'its page chrome must name the pressed button of the group "highlight-mode" by its id.')
)

test_that("a damaged file is refused whole, naming the damage", {
  for (name in names(damages)) {
    damage <- damages[[name]]
    error <- tryCatch({
      cellucid_read_session(rewritten(damage[[1L]]))
      NULL
    }, error = function(error) conditionMessage(error))
    expect_false(is.null(error), info = name)
    if (is.null(error)) next
    if (!is.null(damage[[2L]])) expect_true(grepl(damage[[2L]], error, fixed = TRUE), info = paste(name, error))
  }
  # Trailing bytes, a truncated file, a manifest length past the file.
  bytes <- golden_bytes()
  path <- tempfile(fileext = ".cellucid-session")
  writeBin(c(bytes, as.raw(0)), path)
  expect_error(cellucid_read_session(path), "trailing bytes after declared chunks")
  writeBin(bytes[seq_len(length(bytes) - 10L)], path)
  expect_error(cellucid_read_session(path), "truncated|Unexpected end of file")
})

test_that("a file from a newer Cellucid is newer, never damaged", {
  newer <- rewritten(manifest_edit(function(m) {
    m$app$version <- "99.0.0"
    # A field this package's model does not know (a newer web app's).
    m$trackers[[1L]]$colour <- "teal"
    m
  }))
  package <- as.character(utils::packageVersion("cellucid"))
  expect_error(cellucid_read_session(newer), paste0(
    "This session was written by a newer Cellucid (99.0.0) than this package reads (", package,
    "); update the cellucid package to read it."
  ), fixed = TRUE)
  # A newer file this package reads is read.
  fine <- rewritten(manifest_edit(function(m) {
    m$app$version <- "99.0.0"
    m
  }))
  expect_identical(cellucid_read_session(fine)$app$version, "99.0.0")
  # The same unknown field from this version or an older one is damage.
  older <- rewritten(manifest_edit(function(m) {
    m$app$version <- "0.9.0"
    m$trackers[[1L]]$colour <- "teal"
    m
  }))
  expect_error(cellucid_read_session(older), "^This session file is damaged: ")
})

test_that("a value nested past what R follows is refused in a sentence, as cellucid-python refuses it", {
  # One eager value replaced by lists nested 5000 deep, written as text: JSON
  # jsonlite reads, deeper than R's expression stack follows (before: R's
  # own "evaluation nested too deeply" error).
  skip_if(getRversion() < "4.4.0", "R names a stack overflow's class from 4.4.0 on")
  deep <- rewritten(function(file) {
    at <- which(vapply(file$manifest$chunks, function(entry) entry$kind, "") == "values")[[1L]]
    meta <- file$manifest$chunks[[at]]
    values <- jsonlite::parse_json(rawToChar(session_payload(meta, file$chunks[[at]])), simplifyVector = FALSE)
    values[[1L]] <- "DEEP"
    text <- sub("\"DEEP\"", paste0(strrep("[", 5000L), strrep("]", 5000L)), cellucid:::.js_json(values), fixed = TRUE)
    stored <- session_store(meta[setdiff(names(meta), c("codec", "storedBytes", "uncompressedBytes"))], charToRaw(text))
    file$manifest$chunks[[at]] <- stored$meta
    file$chunks[[at]] <- stored$bytes
    file
  })
  expect_error(cellucid_read_session(deep), cellucid:::.session_too_deep_message, fixed = TRUE)
})

test_that("another format is refused with the web app's sentence", {
  sentences <- list(
    "CELLUCID_SESSION/3\n" = paste(
      "This session was saved by an earlier version of Cellucid (session format 3), which this version no",
      "longer reads: a session file now holds the history of every step (session format 4). It was not opened.",
      "Open it in the version of Cellucid that saved it, or re-create the session in this version."
    ),
    "CELLUCID_SESSION/5\n" = paste(
      "This session was saved by a newer version of Cellucid (session format 5), which this version cannot",
      "read. Open it in that version of Cellucid."
    ),
    "CELLUCID_SESSION/0\n" = "This file is not a Cellucid session file of format 4 (it does not start with CELLUCID_SESSION/4).",
    "PK\003\004 not a session at all" = "This file is not a Cellucid session file of format 4 (it does not start with CELLUCID_SESSION/4)."
  )
  for (prefix in names(sentences)) {
    path <- tempfile(fileext = ".cellucid-session")
    writeBin(c(charToRaw(prefix), as.raw(c(0x10, 0, 0, 0)), charToRaw("{}"), raw(64)), path)
    expect_error(cellucid_read_session(path), sentences[[prefix]], fixed = TRUE)
    expect_error(cellucid_apply_session(data.frame(), path), sentences[[prefix]], fixed = TRUE)
  }
  expect_error(cellucid_read_session(tempfile()), "No session file at")
})

# =============================================================================
# cellucid_pull() of a session file, and a session from a viewer
# =============================================================================

test_that("pulling a session file writes the state it shows", {
  skip_without_sce()
  out <- cellucid_pull(golden_sce(), golden_path(), return_summary = TRUE)
  expect_setequal(names(out$summary$columns), c("cellucid_highlight__highlight_1", "score", "cell_type (merged)"))
  # A pull writes results; renaming the object's own columns is apply's.
  expect_identical(levels(out$data$cell_type), c("B cell", "T cell", "NK cell"))
  expect_identical(names(out$summary$gene_sets), "T cell markers")
  record <- S4Vectors::metadata(out$data)$cellucid$pulled$cellucid_highlight__highlight_1
  expect_identical(c(record$kind, record$id, record$pulled_from), c("highlight", "highlight_1", "history_golden.cellucid-session"))
  again <- cellucid_pull(out$data, golden_path(), return_summary = TRUE)
  expect_setequal(names(again$summary$columns), names(out$summary$columns))
  other <- golden_sce()
  S4Vectors::metadata(other)$cellucid <- list(id = "another")
  expect_error(cellucid_pull(other, golden_path()), "saved on the dataset 'history_golden'")
})

test_that("a session requested from a viewer is collected, then applied", {
  skip_without_sce()
  viewer <- test_viewer(golden_sce(), dataset_id = "history_golden", dataset_name = "History golden")
  client <- ws_welcomed(viewer)
  request <- cellucid_request_session(viewer)
  command <- ws_receive(client)$json
  expect_identical(command$type, "requestSessionBundle")
  header <- charToRaw(paste0("{\"type\":\"session_bundle\",\"requestId\":\"", command$requestId, "\"}"))
  ws_send_binary(client, c(cellucid:::.encode_uint32(length(header)), header, golden_bytes()))
  plan <- cellucid_apply_session(golden_sce(), request, timeout = 10)
  expect_true("score" %in% plan$rows$name)
  expect_error(cellucid_apply_session(golden_sce(), viewer), "Ask the viewer for its session first")
  ws_close(client)
})
