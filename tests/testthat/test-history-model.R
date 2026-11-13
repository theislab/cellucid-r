# The history model R reads sessions with, generated from the web app's
# (cellucid-python's tests/test_history_model.py, case for case).
#
# inst/history_model.schema.json is the web app's historyModelDocument
# (app/history/model.js) byte for byte: its home is the web repository's
# tests/fixtures/history-model.schema.json, and the test compares the two
# files. R/history-model-generated.R is written from it by
# tools/generate-history-model.R and must be current; the reading model
# (R/history.R, the session, the apply plan) exposes every field of every
# record the schema names; the names rule and the constants are the schema's.

# What to run once the schema is the web app's again: the model's R file is
# written from it.
.generate_history_model <- ", then `Rscript tools/generate-history-model.R` here"

schema_path <- function() file.path(cellucid_repository_root(), "inst", "history_model.schema.json")
history_schema <- function() jsonlite::read_json(schema_path(), simplifyVector = FALSE)

test_that("the schema is the web app's, byte for byte", {
  expect_shared_file("inst/history_model.schema.json", "tests/fixtures/history-model.schema.json", then = .generate_history_model)
})

test_that("the schema is the model the web app's source makes now", {
  # Against app/history/model.js itself, not only the web app's committed
  # file: a model changed (a domain added) without writing the copies again
  # fails here, as it fails the web suite.
  schema <- readChar(schema_path(), file.size(schema_path()), useBytes = TRUE)
  expect_identical(cellucid_rendered_by_web("scripts/generate-history-model.mjs", "renderHistoryModel"), schema,
                   info = paste0("Run `", cellucid_sync_shared_files, "` in the web repository (cellucid)", .generate_history_model, "."))
})

test_that("the generated model is current", {
  # The generator is in the checkout only (.Rbuildignore).
  skip_if_not(cellucid_is_source_checkout(), "the generator is in the checkout only")
  generator <- new.env()
  sys.source(file.path(cellucid_repository_root(), "tools", "generate-history-model.R"), envir = generator)
  target <- file.path(cellucid_repository_root(), "R", "history-model-generated.R")
  expect_identical(generator$render_history_model(cellucid_repository_root()), readChar(target, file.size(target), useBytes = TRUE),
                   info = "R/history-model-generated.R is stale: run Rscript tools/generate-history-model.R")
})

test_that("every record of the schema has its fields", {
  schema <- history_schema()
  fields <- cellucid:::.history_fields
  expect_setequal(names(fields), names(schema$`$defs`))
  for (name in names(schema$`$defs`)) {
    definition <- schema$`$defs`[[name]]
    keys <- vapply(fields[[name]], function(field) field$key, "")
    optional <- vapply(fields[[name]], function(field) field$optional, TRUE)
    expect_setequal(keys, names(definition$properties))
    expect_identical(keys[!optional], as.character(unlist(definition$required)), info = name)
  }
})

test_that("the reading model exposes every field of the schema", {
  session <- cellucid_read_session(testthat::test_path("fixtures", "history_golden", "history_golden.cellucid-session"))
  tracker <- cellucid_tracker(session)
  records <- list(
    Tracker = tracker, Branch = tracker$branches[[2L]], Step = cellucid_step(tracker, "s5"), Turn = tracker$turns[[2L]],
    Note = tracker$notes[[1L]], JournalEntry = tracker$journal[[1L]], Start = cellucid_tracker(session, "Annotation (replay)")$start,
    Scope = session$scope, AppStamp = session$app, ActionStamp = session$actions[[1L]], Requires = session$requires, By = tracker$by
  )
  for (name in names(records)) {
    for (field in cellucid:::.history_fields[[name]]) {
      expect_true(field$attr %in% names(records[[name]]) || field$key %in% names(records[[name]]),
                  info = paste0(name, ".", field$key, " is not exposed as ", field$attr))
    }
  }
  # A step table's columns are cellucid-python's.
  expect_identical(names(cellucid_steps(tracker)), cellucid:::.history_step_columns)
})

test_that("names are compared as the web app compares them", {
  expect_identical(cellucid:::.history_normalize_name("  Without  alpha \t"), "Without alpha")
  expect_identical(cellucid:::.history_name_key("ÉQUIPE  One"), "équipe one")
  expect_identical(cellucid:::.history_name_problem("branch", "   "), "A branch needs a name.")
  expect_identical(cellucid:::.history_name_problem("branch", "a\u0007b"), "A branch name cannot hold control characters.")
  expect_identical(cellucid:::.history_name_problem("branch", "B7"),
                   'A branch cannot be called "B7": that is the form of a branch id. Choose another name.')
  expect_null(cellucid:::.history_name_problem("branch", "B2 cells"))
  expect_null(cellucid:::.history_name_problem("tracker", "T-cells"))
  expect_null(cellucid:::.history_name_problem("bookmark", "b7"))
  # The white space set is the schema's, code point for code point.
  expect_identical(sprintf("U+%04X", cellucid:::.history_name_white_space), unlist(history_schema()$`x-names`$whiteSpace))
})

test_that("the constants are the schema's", {
  schema <- history_schema()
  constants <- schema$`x-constants`
  expect_identical(c(cellucid:::.history_root_step, cellucid:::.history_first_branch, cellucid:::.history_first_branch_name),
                   c(constants$rootStep, constants$firstBranch, constants$firstBranchName))
  expect_identical(cellucid:::.history_lane_count, constants$laneCount)
  expect_identical(cellucid:::.history_outside_change_type, constants$outsideChangeType)
  expect_identical(as.list(cellucid:::.history_root_labels), constants$rootLabels)
  expect_identical(cellucid:::.history_domains, unlist(schema$`x-domains`))
  expect_identical(cellucid:::.history_journal_ops, unlist(schema$`x-journalOps`))
  expect_identical(as.list(cellucid:::.history_id_patterns), schema$`x-ids`)
})

test_that("a record that breaks the model is refused naming the path and the rule", {
  by <- list(origin = "agent", agent = NULL, model = NULL, host = NULL, person = NULL, turn = NULL)
  expect_error(cellucid:::.history_check_by(by, "by"), 'by names an agent exactly when its origin is "agent"', class = "cellucid_history_model_error")
  expect_error(cellucid:::.history_check_record("By", c(by, list(colour = 1)), "by"), 'by has an unknown field "colour"')
  user <- by
  user$origin <- "user"
  expect_error(cellucid:::.history_check_record("Branch", list(
    id = "b1", name = "X", base = "s1", tip = "s2", createdAt = 0, by = user, lane = 9L, locked = FALSE, archived = FALSE, notes = list()
  ), "branch"), "branch.lane must be a whole number from 0 to 7", fixed = TRUE)
})

# --- The rules of a paused recording, a squash and a group edit --------------------
# (the web app's tests/history-model-contract.test.mjs and cellucid-python's
# test_history_model.py, on the golden history's own records)

golden_raw_trackers <- function() {
  path <- test_path("fixtures", "history_golden", "history_golden.cellucid-session")
  file <- session_unframe(readBin(path, "raw", file.size(path)))
  ids <- vapply(file$manifest$chunks, function(entry) entry$id, "")
  json <- function(id) {
    at <- match(id, ids)
    jsonlite::parse_json(rawToChar(session_payload(file$manifest$chunks[[at]], file$chunks[[at]])), simplifyVector = FALSE)
  }
  lapply(file$manifest$trackers, function(record) {
    record$steps <- do.call(c, lapply(unlist(record$steps), json))
    record$journal <- if (is.null(record$journal)) list() else json(record$journal)
    record$chats <- list()
    record
  })
}

raw_step_at <- function(tracker, id) which(vapply(tracker$steps, function(step) step$id, "") == id)

test_that("a pause, a squash and a group edit keep their rules", {
  trackers <- golden_raw_trackers()
  where <- paste0("tracker \"", trackers[[1L]]$id, "\"")
  s14 <- raw_step_at(trackers[[1L]], "s14")
  s16 <- raw_step_at(trackers[[1L]], "s16")
  cases <- list(
    list(2L, function(t) { t$paused$by <- list(); t }, "tracker.paused names who paused it (by)"),
    list(2L, function(t) { t$paused$by <- c(t$paused$by, t$paused$by); t }, "tracker.paused.by names one author twice"),
    list(1L, function(t) { t$steps[[s14]]$squashed$carried$bookmarks <- list("Kept"); t },
         'tracker.steps["s14"].squashed carries notes or bookmarks only for squashed steps (cause steps)'),
    list(1L, function(t) { t$steps[[s16]]$squashed$carried$notes <- list("n6", "n6"); t }, 'tracker.steps["s16"].squashed.carried names a note twice'),
    list(1L, function(t) { t$steps[[s16]]$squashed$carried$bookmarks <- list("Before tidying", "BEFORE  tidying"); t },
         'tracker.steps["s16"].squashed.carried names a bookmark twice'),
    list(1L, function(t) { t$steps[[s16]]$notes <- list(); t }, paste0(where, ': step "s16" carries a note it does not hold')),
    list(1L, function(t) { t$steps[[s16]]$squashed$carried$bookmarks <- list("Markers saved"); t },
         paste0(where, ': step "s5" and step "s16" share the name "Markers saved"')),
    list(1L, function(t) { t$journal[[7L]]$after$steps[[2L]]$squashed$count <- 3L; t },
         'tracker.journal["j6"]: its count adds up the steps it replaced (a squashed one counts its own)'),
    list(1L, function(t) { t$journal[[7L]]$after$steps[[2L]]$squashed$carried$bookmarks <- list(); t },
         'tracker.journal["j6"]: the bookmarks of the steps it replaced are carried by the squashed step'),
    list(1L, function(t) { t$journal[[7L]]$after$steps[[2L]]$squashed$carried$notes <- list(); t },
         'tracker.journal["j6"]: the notes of the steps it replaced are carried by the squashed step'),
    list(1L, function(t) { t$journal[[7L]]$before$steps[[2L]]$branch <- "b1"; t }, 'tracker.journal["j6"] replaces consecutive steps of one branch'),
    list(1L, function(t) { t$journal[[7L]]$after$steps[[2L]]$group <- "g0"; t },
         'tracker.journal["j6"]: a squash inside a group keeps it, and one that covers whole groups leaves none'),
    list(1L, function(t) {
      t$journal[[7L]]$after$steps[[2L]]$squashed$cause <- "clear"
      t$journal[[7L]]$after$steps[[2L]]$squashed$carried <- list(notes = list(), bookmarks = list())
      t
    }, 'tracker.journal["j6"] makes exactly one squashed step (cause steps)'),
    list(1L, function(t) { t$journal[[3L]]$after$memberships[[1L]]$from <- "s4"; t },
         'tracker.journal["j2"].after.memberships[0] runs from a later step to an earlier one'),
    list(1L, function(t) { t$journal[[3L]]$after$memberships <- rep(t$journal[[3L]]$after$memberships, 2L); t },
         'tracker.journal["j2"].after.memberships names one group twice'),
    list(1L, function(t) { t$journal[[3L]]$after$steps <- list(t$steps[[raw_step_at(t, "s2")]]); t },
         'tracker.journal["j2"] records which steps a group held (memberships), never the steps')
  )
  for (case in cases) {
    tracker <- trackers[[case[[1L]]]]
    expect_silent(cellucid:::.history_check_tracker_rules(tracker, "tracker"))
    error <- tryCatch({ cellucid:::.history_check_tracker_rules(case[[2L]](tracker), "tracker"); NULL },
                      cellucid_history_model_error = function(error) conditionMessage(error))
    expect_identical(error, case[[3L]])
  }
})
