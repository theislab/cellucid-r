# The history's rules, as the web app keeps them (app/history/model.js), for
# the session reader: every record a session file holds is checked here before
# anything reads it -- its shape against the model's fields
# (R/history-model-generated.R, written from inst/history_model.schema.json),
# then the rules the types cannot say (who acted, a step's own rules by its
# kind, a squashed step's, a group's, a turn's, a journal entry's, and a
# tracker's whole tree: one root, branches that are lines from their fork step
# to their tip, groups that are runs of consecutive steps holding whole turns,
# unique names, the head on the active branch's path, every note and turn
# accounted for; and what a session's scope says it holds). The messages are
# the web app's and cellucid-python's, word for word.
#
# Names (branches, groups, trackers, bookmarks) follow the model's rule: the stored
# form has the model's white space trimmed and each inner run one space, and
# the key is that form lower-cased. Base R has no Unicode normalization, so
# names are compared without NFC here (a file's names are stored in NFC by the
# writer already; one typed in another normalization finds nothing and is
# refused naming what there is).

# A record that breaks the model: an error of class
# cellucid_history_model_error, its message the path and the rule (in the
# session's own encoding, as .abort() writes one).
.history_model_error <- function(message) {
  structure(class = c("cellucid_history_model_error", "error", "condition"), list(message = message, call = NULL))
}

.history_fail <- function(...) stop(.history_model_error(.condition_text(...)))

.history_number <- function(id) as.numeric(substring(id, 2L))

.history_is_id <- function(kind, value) {
  .is_json_string(value) && grepl(.history_id_patterns[[kind]], value, perl = TRUE)
}

.history_id_examples <- c(
  tracker = "t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10", step = "s12", branch = "b3", group = "g1", note = "n4", journal = "j2"
)

# --- Names ------------------------------------------------------------------------

.history_white_class <- paste0("[", intToUtf8(.history_name_white_space), "]")

.history_normalize_name <- function(text) {
  text <- enc2utf8(text)
  text <- gsub(paste0("^", .history_white_class, "+|", .history_white_class, "+$"), "", text, perl = TRUE)
  gsub(paste0(.history_white_class, "+"), " ", text, perl = TRUE)
}

.history_name_key <- function(text) tolower(.history_normalize_name(text))

.history_name_problem <- function(kind, text) {
  name <- .history_normalize_name(text)
  if (!nzchar(name)) return(paste0("A ", kind, " needs a name."))
  if (grepl("[\\x{1}-\\x{1f}\\x{7f}-\\x{9f}]", name, perl = TRUE)) {
    return(paste0("A ", kind, " name cannot hold control characters."))
  }
  if (kind %in% c("branch", "group", "tracker") && grepl(.history_id_patterns[[kind]], .history_name_key(name), perl = TRUE)) {
    return(paste0("A ", kind, " cannot be called \"", name, "\": that is the form of a ", kind, " id. Choose another name."))
  }
  NULL
}

.history_quoted_list <- function(items) {
  quoted <- paste0("\"", items, "\"")
  if (length(quoted) <= 1L) return(paste(quoted, collapse = ""))
  paste0(paste(quoted[-length(quoted)], collapse = ", "), " and ", quoted[[length(quoted)]])
}

.history_no_branch <- function(ref, names) {
  paste0("There is no branch \"", ref, "\"; the branches are ", .history_quoted_list(names), ".")
}

.history_no_tracker <- function(ref, names) {
  paste0("There is no tracker \"", ref, "\"; the trackers are ", .history_quoted_list(names), ".")
}

# The web app's noStep with this package's call in place of the agents'
# history_get (cellucid-python's says Tracker.steps).
.history_no_step <- function(step, tracker) {
  paste0("There is no step \"", step, "\" in the tracker \"", tracker, "\"; cellucid_steps(tracker) lists its steps.")
}

.history_no_group <- function(ref, names) {
  if (length(names) == 0L) return(paste0("There is no group \"", ref, "\"; the tracker has no group."))
  paste0("There is no group \"", ref, "\"; the groups are ", .history_quoted_list(names), ".")
}

# --- Shapes -----------------------------------------------------------------------

.history_expected <- function(spec) {
  base <- switch(spec$type,
    text = "a non-empty text",
    string = if (!is.null(spec$pattern)) paste0("text of the form ", spec$pattern) else "text",
    name = "a name in its stored form (NFC, trimmed, single-spaced)",
    time = "a time in milliseconds since 1970 or null",
    boolean = "true or false",
    hash = "a value hash (64 lower-case hex digits)",
    id = paste0("a ", spec$kind, " id like \"", .history_id_examples[[spec$kind]], "\""),
    integer = if (is.null(spec$maximum)) paste0("a whole number from ", spec$minimum) else
      paste0("a whole number from ", spec$minimum, " to ", spec$maximum),
    number = paste0("a number from ", spec$minimum),
    enum = paste0("one of ", paste0("\"", spec$values, "\"", collapse = ", ")),
    array = "a list",
    pair = "a pair",
    record = paste0("a ", spec$name, " record"),
    object = "a record",
    "a value"
  )
  if (isTRUE(spec$nullable) && !identical(spec$type, "time")) paste0(base, " or null") else base
}

.history_is_number <- function(value) {
  (is.numeric(value) && length(value) == 1L && !is.na(value) && is.finite(value))
}

.history_item_path <- function(path, item, index) {
  if (.is_json_object(item) && .is_json_string(item$id)) paste0(path, "[\"", item$id, "\"]") else paste0(path, "[", index - 1L, "]")
}

.history_check_value <- function(spec, value, path) {
  type <- spec$type
  if (is.null(value) && (isTRUE(spec$nullable) || type == "time" || type == "json")) return(invisible())
  wrong <- function() .history_fail(path, " must be ", .history_expected(spec))
  switch(type,
    text = if (!.is_json_string(value) || !nzchar(value)) wrong(),
    string = if (!.is_json_string(value) || (!is.null(spec$pattern) && !grepl(spec$pattern, value, perl = TRUE))) wrong(),
    name = if (!.is_json_string(value) || !nzchar(value) || !identical(.history_normalize_name(value), enc2utf8(value))) wrong(),
    time = if (!.history_is_number(value)) wrong(),
    boolean = if (!is.logical(value) || length(value) != 1L || is.na(value)) wrong(),
    hash = if (!.is_json_string(value) || !grepl(.history_value_hash_pattern, value, perl = TRUE)) wrong(),
    id = if (!.history_is_id(spec$kind, value)) wrong(),
    integer = {
      ok <- .history_is_number(value) && value == floor(value) && abs(value) <= 2^53 - 1 && value >= spec$minimum &&
        (is.null(spec$maximum) || value <= spec$maximum)
      if (!ok) wrong()
    },
    number = if (!.history_is_number(value) || value < spec$minimum) wrong(),
    enum = if (!.is_json_string(value) || !value %in% spec$values) wrong(),
    array = {
      if (!.is_json_array(value)) wrong()
      for (index in seq_along(value)) .history_check_value(spec$item, value[[index]], .history_item_path(path, value[[index]], index))
    },
    pair = {
      if (!.is_json_array(value) || length(value) != 2L) wrong()
      for (index in 1:2) .history_check_value(spec$items[[index]], value[[index]], paste0(path, "[", index - 1L, "]"))
    },
    record = .history_check_record(spec$name, value, path),
    object = if (!.is_json_object(value)) wrong(),
    json = invisible(),
    .history_fail(path, " has an unknown spec")
  )
  invisible()
}

# Each record's field keys, in order (looked up once, not per record checked).
.history_field_keys <- lapply(.history_fields, function(fields) vapply(fields, function(field) field$key, ""))

.history_check_record <- function(name, value, path) {
  fields <- .history_fields[[name]]
  if (!.is_json_object(value)) .history_fail(path, " must be a ", name, " record")
  keys <- .history_field_keys[[name]]
  unknown <- setdiff(names(value), keys)
  if (length(unknown) > 0L) .history_fail(path, " has an unknown field \"", unknown[[1L]], "\"")
  present <- keys %in% names(value)
  for (at in seq_along(fields)) {
    field <- fields[[at]]
    if (!present[[at]]) {
      if (field$optional) next
      .history_fail(path, " lacks \"", field$key, "\"")
    }
    .history_check_value(field$spec, value[[field$key]], paste0(path, ".", field$key))
  }
  invisible()
}

# --- Record rules -----------------------------------------------------------------

.history_check_by <- function(by, path) {
  origin <- by$origin
  if (identical(origin, "agent") != !is.null(by$agent)) .history_fail(path, " names an agent exactly when its origin is \"agent\"")
  if (!is.null(by$model) && !identical(origin, "agent")) .history_fail(path, " names a model only for an agent")
  if (identical(origin, "host") != !is.null(by$host)) .history_fail(path, " names a host exactly when its origin is \"host\"")
  if (!is.null(by$turn) && !origin %in% c("agent", "host")) .history_fail(path, " belongs to a turn only for an agent or a host")
}

.history_same_by <- function(a, b) {
  all(vapply(.history_fields$By, function(field) identical(a[[field$key]], b[[field$key]]), TRUE))
}

.history_check_notes <- function(notes, path) {
  for (note in notes) {
    label <- paste0(path, "[\"", note$id, "\"]")
    .history_check_by(note$by, paste0(label, ".by"))
    if (isTRUE(note$fromReply) && !identical(note$by$origin, "agent")) {
      .history_fail(label, " is from an agent's reply only when an agent is its author")
    }
  }
}

.history_same_keys <- function(record, keys) {
  setequal(names(record) %||% character(), unlist(keys) %||% character()) &&
    length(names(record) %||% character()) == length(unlist(keys) %||% character())
}

.history_check_step_rules <- function(step, path) {
  .history_check_by(step$by, paste0(path, ".by"))
  .history_check_notes(step$notes, paste0(path, ".notes"))
  number <- .history_number(step$id)
  children <- unlist(step$children) %||% character()
  if (anyDuplicated(children)) .history_fail(path, " lists a child twice")
  if (any(.history_number(children) <= number)) .history_fail(path, " has a child numbered before it")
  if (!is.null(step$at) && !is.null(step$until) && step$until < step$at) .history_fail(path, " ends (until) before it starts (at)")
  if (identical(step$status, "unapplicable") && is.null(step$reason)) .history_fail(path, " no longer applies and must say why (reason)")
  kind <- step$kind
  if (is.null(step$parent) != identical(kind, "root")) .history_fail(path, " is the root (kind \"root\") exactly when it has no parent")
  if (identical(kind, "squashed") != !is.null(step$squashed)) .history_fail(path, " says what it stands for (squashed) exactly when it is squashed")
  if (!identical(kind, "action") && !is.null(step$v)) .history_fail(path, " records a parameters version (v) only for an action")
  if (is.null(step$parent)) {
    if (!identical(step$id, .history_root_step)) .history_fail(path, " has no parent, and only the root (", .history_root_step, ") has none")
    if (!is.null(step$type) || !is.null(step$params) || !is.null(step$code) || length(step$touches) > 0L ||
        length(step$before) > 0L || length(step$reads) > 0L || !is.null(step$outcome) || !is.null(step$coalesceKey) ||
        !is.null(step$replayOf) || step$merged != 1 || !is.null(step$ms) || !identical(step$status, "applied") ||
        !is.null(step$reason) || !is.null(step$group)) {
      .history_fail(path, " is the root and records no action")
    }
    return(invisible())
  }
  if (.history_number(step$parent) >= number) .history_fail(path, " is numbered before its parent")
  if (identical(kind, "outside") != identical(step$type, .history_outside_change_type)) {
    .history_fail(path, " is of kind \"outside\" exactly when it is a change made outside an action (", .history_outside_change_type, ")")
  }
  if (identical(step$type, .history_outside_change_type) != identical(step$by$origin, "outside")) {
    .history_fail(path, " is by \"outside\" exactly when it is a change made outside an action (", .history_outside_change_type, ")")
  }
  if (identical(kind, "squashed")) {
    if (!is.null(step$type) || !is.null(step$params) || !is.null(step$code) || !is.null(step$coalesceKey) || !is.null(step$replayOf)) {
      .history_fail(path, " is squashed and records no action of its own (type, params, code, coalesceKey and replayOf null)")
    }
    .history_check_squash_rules(step$squashed, paste0(path, ".squashed"))
  } else {
    if (is.null(step$type)) .history_fail(path, " names no action")
    if (is.null(step$params)) .history_fail(path, " records an action and its parameters (params)")
  }
  if (identical(kind, "action") && is.null(step$v)) .history_fail(path, " records the version of its action's parameters (v)")
  touches <- unlist(step$touches) %||% character()
  if (anyDuplicated(touches)) .history_fail(path, " lists a touched domain twice")
  if (!identical(.history_domains[.history_domains %in% touches], touches)) .history_fail(path, " lists its touched domains out of domain order")
  if (!.history_same_keys(step$before, touches) || !.history_same_keys(step$after, touches)) {
    .history_fail(path, " must hold a value before and after for exactly the domains it touches")
  }
  invisible()
}

# What a squashed step stands for (checkSquashRules).
.history_check_squash_rules <- function(squash, path) {
  counted <- squash$cause %in% c("steps", "clear")
  if (counted != !is.null(squash$count)) {
    .history_fail(path, " counts the steps it replaced (count) exactly when they were recorded and shown (steps, clear)")
  }
  if (identical(squash$cause, "share") != (length(squash$by) == 0L)) .history_fail(path, " names who acted (by) exactly unless it is a shared state")
  .history_check_distinct_bys(squash$by, paste0(path, ".by"))
  notes <- unlist(squash$carried$notes) %||% character()
  bookmarks <- unlist(squash$carried$bookmarks) %||% character()
  if (!identical(squash$cause, "steps") && (length(notes) > 0L || length(bookmarks) > 0L)) {
    .history_fail(path, " carries notes or bookmarks only for squashed steps (cause steps)")
  }
  if (anyDuplicated(notes)) .history_fail(path, ".carried names a note twice")
  if (anyDuplicated(vapply(bookmarks, .history_name_key, "", USE.NAMES = FALSE))) .history_fail(path, ".carried names a bookmark twice")
}

.history_check_distinct_bys <- function(bys, path) {
  for (index in seq_along(bys)) .history_check_by(bys[[index]], paste0(path, "[", index - 1L, "]"))
  for (index in seq_along(bys)[-1L]) {
    if (any(vapply(bys[seq_len(index - 1L)], function(other) .history_same_by(other, bys[[index]]), TRUE))) {
      .history_fail(path, " names one author twice")
    }
  }
}

# A tracker's paused recording (checkPauseRules).
.history_check_pause_rules <- function(pause, path) {
  if (length(pause$by) == 0L) .history_fail(path, " names who paused it (by)")
  .history_check_distinct_bys(pause$by, paste0(path, ".by"))
}

# How many recorded steps a step stands for: a squashed step its count (one
# when none was recorded), any other one.
.history_step_count <- function(step) {
  if (identical(step$kind, "squashed") && !is.null(step$squashed$count)) step$squashed$count else 1
}

.history_carried <- function(step, what) {
  if (identical(step$kind, "squashed")) unlist(step$squashed$carried[[what]]) %||% character() else character()
}

# A squash's journal entry held to the squash rules (checkSquashEntry): the
# squashed step takes the id of the range's last step; the range is one line
# of one branch; its count adds up what it replaced; the replaced steps' notes
# and bookmarks are carried; it stays inside one group or leaves none.
.history_check_squash_entry <- function(entry, path) {
  before_ids <- vapply(entry$before$steps, function(step) step$id, "")
  made <- Filter(function(step) {
    if (!identical(step$kind, "squashed") || !identical(step$squashed$cause, "steps")) return(FALSE)
    at <- match(step$id, before_ids)
    if (is.na(at)) return(TRUE)
    was <- entry$before$steps[[at]]
    !identical(was$kind, "squashed") || !identical(was$squashed, step$squashed)
  }, entry$after$steps)
  if (length(made) != 1L) .history_fail(path, " makes exactly one squashed step (cause steps)")
  squashed <- made[[1L]]
  kept <- vapply(entry$after$steps, function(step) step$id, "")
  replaced <- Filter(function(step) !step$id %in% kept || identical(step$id, squashed$id), entry$before$steps)
  replaced_ids <- vapply(replaced, function(step) step$id, "")
  last <- replaced[match(squashed$id, replaced_ids)]
  if (anyNA(match(squashed$id, replaced_ids)) || any(.history_number(replaced_ids) > .history_number(squashed$id))) {
    .history_fail(path, ": the squashed step takes the id of the last step it replaced")
  }
  last <- last[[1L]]
  ordered <- replaced[order(.history_number(replaced_ids))]
  for (index in seq_along(ordered)) {
    if (!identical(ordered[[index]]$branch, squashed$branch) || (index > 1L && !identical(ordered[[index]]$parent, ordered[[index - 1L]]$id))) {
      .history_fail(path, " replaces consecutive steps of one branch")
    }
  }
  if (!identical(squashed$parent, ordered[[1L]]$parent) ||
      !identical(unlist(squashed$children) %||% character(), unlist(last$children) %||% character())) {
    .history_fail(path, ": the squashed step keeps the first replaced step's parent and the last one's children")
  }
  if (!identical(as.numeric(squashed$squashed$count), as.numeric(sum(vapply(replaced, .history_step_count, 0))))) {
    .history_fail(path, ": its count adds up the steps it replaced (a squashed one counts its own)")
  }
  notes <- unlist(lapply(replaced, function(step) c(vapply(step$notes, function(note) note$id, ""), .history_carried(step, "notes")))) %||% character()
  bookmarks <- unlist(lapply(replaced, function(step) c(step$bookmark, .history_carried(step, "bookmarks")))) %||% character()
  own <- vapply(squashed$notes, function(note) note$id, "")
  carried_notes <- unlist(squashed$squashed$carried$notes) %||% character()
  if (!all(notes %in% own) || !identical(sort(unique(notes), method = "radix"), sort(carried_notes, method = "radix"))) {
    .history_fail(path, ": the notes of the steps it replaced are carried by the squashed step")
  }
  if (!identical(sort(bookmarks, method = "radix"), sort(unlist(squashed$squashed$carried$bookmarks) %||% character(), method = "radix"))) {
    .history_fail(path, ": the bookmarks of the steps it replaced are carried by the squashed step")
  }
  groups <- unique(lapply(replaced, function(step) step$group))
  if (!is.null(squashed$group) && !(length(groups) == 1L && identical(groups[[1L]], squashed$group))) {
    .history_fail(path, ": a squash inside a group keeps it, and one that covers whole groups leaves none")
  }
}

# A group's own rules (checkGroupRules).
.history_check_group_rules <- function(group, path) {
  .history_check_by(group$by, paste0(path, ".by"))
  .history_check_notes(group$notes, paste0(path, ".notes"))
}

.history_check_turn_rules <- function(turn, path) {
  .history_check_by(turn$by, paste0(path, ".by"))
  if (!turn$by$origin %in% c("agent", "host")) .history_fail(path, " is an agent's or a notebook's")
  if (!identical(turn$by$turn, turn$id)) .history_fail(path, ".by.turn must be the turn's own id")
  if (!is.null(turn$startedAt) && !is.null(turn$endedAt) && turn$endedAt < turn$startedAt) .history_fail(path, " ends before it starts")
  .history_check_notes(turn$notes, paste0(path, ".notes"))
}

.history_check_start_rules <- function(start, path) {
  replay <- identical(start$kind, "replay")
  linked <- replay || identical(start$kind, "import")
  if (linked != !is.null(start$from)) .history_fail(path, " names where it comes from (from) exactly for a replay or an import")
  if (replay != !is.null(start$report)) .history_fail(path, " keeps a report exactly for a replay")
  if (!replay && !is.null(start$mapping)) .history_fail(path, " keeps a name mapping only for a replay")
  if (identical(start$kind, "import") && is.null(start$from$file)) .history_fail(path, " of an import names its file")
}

.history_check_patch_rules <- function(patch, path) {
  for (step in patch$steps) .history_check_step_rules(step, paste0(path, ".steps[\"", step$id, "\"]"))
  for (branch in patch$branches) {
    .history_check_by(branch$by, paste0(path, ".branches[\"", branch$id, "\"].by"))
    .history_check_notes(branch$notes, paste0(path, ".branches[\"", branch$id, "\"].notes"))
  }
  for (group in patch$groups) .history_check_group_rules(group, paste0(path, ".groups[\"", group$id, "\"]"))
  for (index in seq_along(patch$memberships)) {
    membership <- patch$memberships[[index]]
    if (.history_number(membership$from) > .history_number(membership$to)) {
      .history_fail(path, ".memberships[", index - 1L, "] runs from a later step to an earlier one")
    }
  }
  if (anyDuplicated(vapply(patch$memberships, function(membership) membership$group, ""))) .history_fail(path, ".memberships names one group twice")
  for (turn in patch$turns) .history_check_turn_rules(turn, paste0(path, ".turns[\"", turn$id, "\"]"))
  for (kind in c("steps", "branches", "groups", "turns")) {
    ids <- vapply(patch[[kind]], function(item) item$id, "")
    if (anyDuplicated(ids)) .history_fail(path, ".", kind, " holds one record twice")
  }
}

.history_check_journal_rules <- function(entry, path) {
  .history_check_by(entry$by, paste0(path, ".by"))
  if (entry$purged != is.null(entry$before) || entry$purged != is.null(entry$after)) {
    .history_fail(path, " keeps what it replaced and made (before, after) exactly until it is purged")
  }
  if (entry$purged) return(invisible())
  .history_check_patch_rules(entry$before, paste0(path, ".before"))
  .history_check_patch_rules(entry$after, paste0(path, ".after"))
  if (!.history_same_keys(entry$before$tracker, names(entry$after$tracker))) {
    .history_fail(path, " must replace the same tracker fields it puts back")
  }
  if (startsWith(entry$op, "group.") && (length(entry$before$steps) > 0L || length(entry$after$steps) > 0L)) {
    .history_fail(path, " records which steps a group held (memberships), never the steps")
  }
  if (identical(entry$op, "squash")) .history_check_squash_entry(entry, path)
}

.history_in_id_order <- function(items) {
  numbers <- vapply(items, function(item) .history_number(item$id), 0)
  length(numbers) < 2L || all(diff(numbers) > 0)
}

.history_check_unique_names <- function(items, kind, label, name_of, path) {
  seen <- new.env(parent = emptyenv())
  for (item in items) {
    name <- name_of(item)
    if (is.null(name)) next
    problem <- .history_name_problem(kind, name)
    if (!is.null(problem)) .history_fail(path, ": ", problem)
    key <- .history_name_key(name)
    if (!is.null(seen[[key]])) .history_fail(path, ": ", label(seen[[key]]), " and ", label(item), " share the name \"", name, "\"")
    assign(key, item, envir = seen)
  }
}

# One tracker's rules beyond its records' shapes (checked already).
.history_check_tracker_rules <- function(tracker, path) {
  where <- paste0(path, " \"", tracker$id, "\"")
  .history_check_by(tracker$by, paste0(path, ".by"))
  .history_check_start_rules(tracker$start, paste0(path, ".start"))
  if (!is.null(tracker$removed)) .history_check_by(tracker$removed$by, paste0(path, ".removed.by"))
  .history_check_notes(tracker$notes, paste0(path, ".notes"))

  steps <- new.env(parent = emptyenv())
  for (step in tracker$steps) {
    if (!is.null(steps[[step$id]])) .history_fail(where, " holds the step \"", step$id, "\" twice")
    assign(step$id, step, envir = steps)
    .history_check_step_rules(step, paste0(path, ".steps[\"", step$id, "\"]"))
    if (.history_number(step$id) >= tracker$nextStep) .history_fail(where, ": step \"", step$id, "\" is numbered past nextStep")
  }
  if (!.history_in_id_order(tracker$steps)) .history_fail(where, " lists its steps out of id order")
  if (!identical(tracker$rootId, .history_root_step)) .history_fail(where, "'s root must be ", .history_root_step)
  root <- steps[[.history_root_step]]
  if (is.null(root)) .history_fail(where, " has no root step ", .history_root_step)
  if (!.history_same_by(root$by, tracker$by)) .history_fail(where, "'s root is by whoever started it (tracker.by)")
  for (step in tracker$steps) {
    if (!is.null(step$parent)) {
      parent <- steps[[step$parent]]
      if (is.null(parent) || !step$id %in% unlist(parent$children)) .history_fail(where, ": step \"", step$id, "\" and its parent disagree")
    }
    for (child in unlist(step$children)) {
      found <- steps[[child]]
      if (is.null(found) || !identical(found$parent, step$id)) .history_fail(where, ": step \"", step$id, "\" lists a child that is not its own")
    }
  }

  branches <- new.env(parent = emptyenv())
  for (branch in tracker$branches) {
    if (!is.null(branches[[branch$id]])) .history_fail(where, " holds the branch \"", branch$id, "\" twice")
    assign(branch$id, branch, envir = branches)
    .history_check_by(branch$by, paste0(path, ".branches[\"", branch$id, "\"].by"))
    .history_check_notes(branch$notes, paste0(path, ".branches[\"", branch$id, "\"].notes"))
    if (.history_number(branch$id) >= tracker$nextBranch) .history_fail(where, ": branch \"", branch$id, "\" is numbered past nextBranch")
  }
  if (!.history_in_id_order(tracker$branches)) .history_fail(where, " lists its branches out of id order")
  .history_check_unique_names(tracker$branches, "branch", function(item) paste0("branch \"", item$id, "\""), function(item) item$name, where)
  # Each branch's steps, in id order (one pass: a history of many steps on one branch costs what its steps cost).
  step_ids <- vapply(tracker$steps, function(step) step$id, "")
  step_branches <- vapply(tracker$steps, function(step) step$branch, "")
  branch_ids <- vapply(tracker$branches, function(branch) branch$id, "")
  stray <- which(!step_branches %in% branch_ids)
  if (length(stray) > 0L) {
    .history_fail(where, ": step \"", step_ids[[stray[[1L]]]], "\" is on the branch \"", step_branches[[stray[[1L]]]], "\", which it does not hold")
  }
  own_steps <- split(step_ids, factor(step_branches, levels = branch_ids))
  root_branches <- Filter(function(branch) is.null(branch$base), tracker$branches)
  if (length(root_branches) != 1L) .history_fail(where, " must have exactly one branch without a fork step (the one holding the root)")
  if (!identical(root$branch, root_branches[[1L]]$id)) .history_fail(where, ": the root must be on the branch without a fork step")
  for (branch in tracker$branches) {
    label <- paste0(where, ": branch \"", branch$id, "\"")
    own <- own_steps[[branch$id]]
    if (length(own) == 0L) .history_fail(label, " has no step of its own")
    if (!is.null(branch$base)) {
      base <- steps[[branch$base]]
      if (is.null(base)) .history_fail(label, " forks from a step the tracker does not hold")
      if (identical(base$branch, branch$id)) .history_fail(label, " forks from one of its own steps")
    }
    firsts <- Filter(function(id) {
      step <- steps[[id]]
      is.null(step$parent) || !identical(steps[[step$parent]]$branch, branch$id)
    }, own)
    if (length(firsts) != 1L || !identical(steps[[firsts[[1L]]]]$parent, branch$base)) {
      .history_fail(label, "'s steps must form one line starting right after its fork step")
    }
    count <- 1L
    current <- steps[[firsts[[1L]]]]
    repeat {
      following <- Filter(function(child) identical(steps[[child]]$branch, branch$id), unlist(current$children))
      if (length(following) > 1L) .history_fail(label, " splits at step \"", current$id, "\"; a branch is one line")
      if (length(following) == 0L) break
      current <- steps[[following[[1L]]]]
      count <- count + 1L
    }
    if (count != length(own)) .history_fail(label, "'s steps must form one line starting right after its fork step")
    if (!identical(current$id, branch$tip)) .history_fail(label, "'s tip must be its last step (\"", current$id, "\")")
  }

  head <- steps[[tracker$headId]]
  if (is.null(head)) .history_fail(where, "'s head is not one of its steps")
  active <- branches[[tracker$activeBranch]]
  if (is.null(active)) .history_fail(where, "'s active branch is not one of its branches")
  on_path <- FALSE
  walker <- steps[[active$tip]]
  while (!is.null(walker)) {
    if (identical(walker$id, head$id)) {
      on_path <- TRUE
      break
    }
    walker <- if (is.null(walker$parent)) NULL else steps[[walker$parent]]
  }
  if (!on_path) .history_fail(where, "'s head must be on its active branch's path")
  if (isTRUE(active$archived) && !identical(head$branch, active$id)) .history_fail(where, "'s active branch is archived, so the head must be on it")
  main <- branches[[tracker$mainBranch]]
  if (is.null(main)) .history_fail(where, "'s main branch is not one of its branches")
  if (isTRUE(main$archived)) .history_fail(where, "'s main branch cannot be archived")
  if (!is.null(tracker$paused)) .history_check_pause_rules(tracker$paused, paste0(path, ".paused"))

  group_ids <- vapply(tracker$groups, function(group) group$id, "")
  if (anyDuplicated(group_ids)) .history_fail(where, " holds the group \"", group_ids[anyDuplicated(group_ids)], "\" twice")
  for (group in tracker$groups) {
    .history_check_group_rules(group, paste0(path, ".groups[\"", group$id, "\"]"))
    if (.history_number(group$id) >= tracker$nextGroup) .history_fail(where, ": group \"", group$id, "\" is numbered past nextGroup")
  }
  if (!.history_in_id_order(tracker$groups)) .history_fail(where, " lists its groups out of id order")
  .history_check_unique_names(tracker$groups, "group", function(item) paste0("group \"", item$id, "\""), function(item) item$name, where)
  step_groups <- vapply(tracker$steps, function(step) step$group %||% NA_character_, "")
  stray <- which(!is.na(step_groups) & !step_groups %in% group_ids)
  if (length(stray) > 0L) {
    .history_fail(where, ": step \"", step_ids[[stray[[1L]]]], "\" is in the group \"", step_groups[[stray[[1L]]]], "\", which it does not hold")
  }
  grouped <- which(!is.na(step_groups))
  members <- split(grouped, factor(step_groups[grouped], levels = group_ids))
  for (group_id in group_ids) {
    own <- tracker$steps[members[[group_id]]]
    if (length(own) == 0L) .history_fail(where, ": group \"", group_id, "\" has no step")
    # Steps are in id order, and along one branch a child's id is above its parent's.
    for (index in seq_along(own)[-1L]) {
      if (!identical(own[[index]]$branch, own[[1L]]$branch) || !identical(own[[index]]$parent, own[[index - 1L]]$id)) {
        .history_fail(where, ": group \"", group_id, "\"'s steps must be consecutive steps of one branch")
      }
    }
  }
  # A group never cuts a run of one turn's consecutive steps on its branch: a turn that acts
  # again after someone else's step starts a new run, which the group need not hold.
  for (step in tracker$steps) {
    if (is.null(step$by$turn) || is.null(step$parent)) next
    parent <- steps[[step$parent]]
    if (!is.null(parent) && identical(parent$branch, step$branch) && identical(parent$by$turn, step$by$turn) &&
        !identical(parent$group, step$group)) {
      .history_fail(where, ": a group cuts the run of turn \"", step$by$turn, "\" between steps \"", parent$id, "\" and \"", step$id,
                    "\"; a group holds whole runs of a turn")
    }
  }

  bookmarked <- unlist(lapply(tracker$steps, function(step) {
    lapply(c(step$bookmark, .history_carried(step, "bookmarks")), function(name) list(id = step$id, name = name))
  }), recursive = FALSE)
  .history_check_unique_names(bookmarked, "bookmark", function(item) paste0("step \"", item$id, "\""), function(item) item$name, where)
  for (step in tracker$steps) {
    if (!identical(step$kind, "squashed")) next
    own <- vapply(step$notes, function(note) note$id, "")
    if (!all(.history_carried(step, "notes") %in% own)) .history_fail(where, ": step \"", step$id, "\" carries a note it does not hold")
  }
  every_note <- c(
    tracker$notes,
    unlist(lapply(tracker$branches, function(branch) branch$notes), recursive = FALSE),
    unlist(lapply(tracker$groups, function(group) group$notes), recursive = FALSE),
    unlist(lapply(tracker$steps, function(step) step$notes), recursive = FALSE),
    unlist(lapply(tracker$turns, function(turn) turn$notes), recursive = FALSE)
  )
  note_ids <- vapply(every_note, function(note) note$id, "")
  if (anyDuplicated(note_ids)) .history_fail(where, " holds the note \"", note_ids[anyDuplicated(note_ids)], "\" twice")
  past <- note_ids[.history_number(note_ids) >= tracker$nextNote]
  if (length(past) > 0L) .history_fail(where, ": note \"", past[[1L]], "\" is numbered past nextNote")
  turn_ids <- vapply(tracker$turns, function(turn) turn$id, "")
  if (anyDuplicated(turn_ids)) .history_fail(where, " holds the turn \"", turn_ids[anyDuplicated(turn_ids)], "\" twice")
  for (turn in tracker$turns) {
    .history_check_turn_rules(turn, paste0(path, ".turns[\"", turn$id, "\"]"))
    if (is.null(branches[[turn$branch]])) .history_fail(where, ": turn \"", turn$id, "\" acts on a branch the tracker does not hold")
  }
  turn_of <- vapply(tracker$steps, function(step) step$by$turn %||% NA_character_, "")
  stray <- which(!is.na(turn_of) & !turn_of %in% turn_ids)
  if (length(stray) > 0L) {
    .history_fail(where, ": step \"", step_ids[[stray[[1L]]]], "\" belongs to the turn \"", turn_of[[stray[[1L]]]], "\", which it does not hold")
  }
  # A turn that only acted outside the state (a vote published) keeps that on record.
  for (turn in tracker$turns) {
    if (!turn$id %in% turn_of && length(turn$outside) == 0L) .history_fail(where, ": turn \"", turn$id, "\" has no step and did nothing outside the state")
  }

  journal_ids <- vapply(tracker$journal, function(entry) entry$id, "")
  if (anyDuplicated(journal_ids)) .history_fail(where, " holds the journal entry \"", journal_ids[anyDuplicated(journal_ids)], "\" twice")
  for (entry in tracker$journal) {
    if (.history_number(entry$id) >= tracker$nextJournal) .history_fail(where, ": journal entry \"", entry$id, "\" is numbered past nextJournal")
    .history_check_journal_rules(entry, paste0(path, ".journal[\"", entry$id, "\"]"))
  }
  if (!.history_in_id_order(tracker$journal)) .history_fail(where, " lists its journal out of id order")
  if (!is.null(tracker$removed) &&
      !any(vapply(tracker$journal, function(entry) identical(entry$op, "tracker.remove") && !isTRUE(entry$purged), TRUE))) {
    .history_fail(where, " is removed, so its journal must keep the removal (to restore it)")
  }
  invisible(tracker)
}

# A session file as read: its header records, every tracker, and its scope
# naming exactly the trackers it holds (validateSession).
.history_validate_session <- function(session) {
  .history_check_record("Session", session, "session")
  ids <- character()
  for (tracker in session$trackers) {
    .history_check_tracker_rules(tracker, "tracker")
    if (tracker$id %in% ids) .history_fail("session holds the tracker \"", tracker$id, "\" twice")
    ids <- c(ids, tracker$id)
  }
  .history_check_unique_names(
    Filter(function(tracker) is.null(tracker$removed), session$trackers), "tracker",
    function(item) paste0("tracker \"", item$id, "\""), function(item) item$name, "session"
  )
  active <- Filter(function(tracker) identical(tracker$id, session$activeTracker), session$trackers)
  if (length(active) == 0L) .history_fail("session's active tracker is not one of its trackers")
  if (!is.null(active[[1L]]$removed)) .history_fail("session's active tracker cannot be removed")
  if (!identical(unlist(session$scope$trackers) %||% character(), ids)) {
    .history_fail("session.scope must name exactly the trackers it holds, in order")
  }
  for (entry in session$scope$branches) {
    tracker <- Filter(function(candidate) identical(candidate$id, entry$tracker), session$trackers)
    if (length(tracker) == 0L) .history_fail("session.scope names branches of a tracker it does not hold (\"", entry$tracker, "\")")
    known <- vapply(tracker[[1L]]$branches, function(branch) branch$id, "")
    if (!all(unlist(entry$branches) %in% known)) .history_fail("session.scope names a branch the tracker \"", entry$tracker, "\" does not hold")
  }
  types <- vapply(session$actions, function(stamp) stamp$type, "")
  versions <- vapply(session$actions, function(stamp) as.numeric(stamp$version), 0)
  ordered <- order(types, versions, method = "radix")
  if (!identical(ordered, seq_along(types)) || anyDuplicated(data.frame(types, versions))) {
    .history_fail("session.actions must be sorted by type then version, each pair once")
  }
  used <- character()
  for (tracker in session$trackers) {
    journal_steps <- unlist(lapply(Filter(function(entry) !isTRUE(entry$purged), tracker$journal),
                                   function(entry) c(entry$before$steps, entry$after$steps)), recursive = FALSE)
    for (step in c(tracker$steps, journal_steps)) {
      if (identical(step$kind, "action")) used <- c(used, paste0(step$type, "@", format(step$v, scientific = FALSE)))
    }
    if (!is.null(tracker$paused) && !is.null(tracker$paused$draft)) {
      .history_fail("session: tracker \"", tracker$id, "\" keeps a paused draft; a file squashes it into a step first")
    }
  }
  stamped <- if (length(types) == 0L) character() else paste0(types, "@", format(versions, scientific = FALSE, trim = TRUE))
  if (!setequal(stamped, unique(used)) || length(stamped) != length(unique(used))) {
    .history_fail("session.actions must stamp exactly each action and version its steps use")
  }
  .history_check_scope_holds(session$scope, session$trackers)
  invisible(session)
}

# What the share chose is what the file holds (checkScopeHolds).
.history_check_scope_holds <- function(scope, trackers) {
  if (identical(scope$content, "steps") && !identical(scope$results, "recompute")) {
    .history_fail("session.scope of steps only holds no results (results is recompute)")
  }
  for (tracker in trackers) {
    where <- paste0("session: tracker \"", tracker$id, "\"")
    if (!isTRUE(scope$chats) && length(tracker$chats) > 0L) .history_fail(where, " holds conversations its scope leaves out (chats)")
    if (!isTRUE(scope$journal) && length(tracker$journal) > 0L) .history_fail(where, " holds Recently deleted, which its scope leaves out (journal)")
    if (!isTRUE(scope$notes)) {
      items <- c(list(tracker), tracker$branches, tracker$groups, tracker$steps, tracker$turns)
      noted <- any(vapply(items, function(item) length(item$notes) > 0L, TRUE)) ||
        any(vapply(tracker$steps, function(step) !is.null(step$bookmark), TRUE))
      if (noted) .history_fail(where, " holds notes or bookmarks its scope leaves out (notes)")
    }
  }
  if (!identical(scope$content, "state")) return(invisible())
  if (length(trackers) != 1L || length(trackers[[1L]]$branches) != 1L) .history_fail("session.scope of a state holds one tracker of one branch")
  tracker <- trackers[[1L]]
  steps <- tracker$steps
  if (length(steps) != 2L || !identical(steps[[2L]]$kind, "squashed") || !identical(steps[[2L]]$squashed$cause, "share") ||
      !identical(tracker$headId, steps[[2L]]$id)) {
    .history_fail("session.scope of a state holds the root and one squashed step (cause share), shown")
  }
  if (length(tracker$turns) > 0L || length(tracker$groups) > 0L || length(tracker$notes) > 0L || length(steps[[1L]]$notes) > 0L ||
      length(steps[[2L]]$notes) > 0L || length(tracker$journal) > 0L || length(tracker$chats) > 0L) {
    .history_fail("session.scope of a state holds no turn, group, note, journal or conversation")
  }
  # Privacy: a state shared without its history names nobody (no agent, model, notebook, person or turn).
  anonymous <- function(by) {
    by$origin %in% c("user", "startup") && is.null(by$agent) && is.null(by$model) && is.null(by$host) && is.null(by$person) &&
      is.null(by$turn)
  }
  if (!all(vapply(c(list(tracker$by, steps[[1L]]$by, steps[[2L]]$by), tracker$paused$by), anonymous, TRUE))) {
    .history_fail("session.scope of a state names nobody: every by is the person or the opening, with no agent, model, notebook, name or turn")
  }
  invisible()
}
