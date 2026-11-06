# A Cellucid history, read: trackers, branches, groups, steps, notes, turns and
# the journal -- the reading model cellucid-python's cellucid.history is, function
# for function. Every record is the history model's (inst/history_model.schema.json,
# the web app's app/history/model.js), its fields named in snake case
# (headId is head_id), built from the field table R/history-model-generated.R
# holds: cellucid_read_session() (R/session-file.R) builds them.
#
# - cellucid_trackers(session) and cellucid_tracker(session, tracker): the
#   trackers, one by name or id;
# - cellucid_branches(tracker), cellucid_branch(tracker, branch);
# - cellucid_groups(tracker), cellucid_group(tracker, group): named runs of
#   consecutive steps of one branch (each step names its group);
# - cellucid_steps(tracker, branch = NULL): every step as a data frame (the
#   same columns as Python's, in step order: ids in order, parents before
#   children), or a branch's path from the root to its tip; a step is of a
#   kind (root, action with its parameters' version v, outside, or squashed:
#   no action of its own, `squashed` says why and who acted);
# - cellucid_step(tracker, step), cellucid_step_state(), cellucid_step_diff();
# - cellucid_notes(), cellucid_turns(), cellucid_journal(), cellucid_stale();
# - cellucid_methods() and cellucid_pseudocode(): a branch in the web app's
#   words (app/history/exports/{methods,pseudocode}.js's rules: the abstract,
#   group headings and markers, bookmarks, squashed steps said plainly; the
#   pseudocode from each step's stored line).

# The columns of a tracker's step table and note table, in order (the same in Python).
.history_step_columns <- c(
  "id", "branch", "parent", "kind", "type", "v", "label", "code", "origin", "agent", "model", "host", "turn", "at",
  "status", "stale", "bookmark", "group", "n_notes"
)
.history_note_columns <- c("id", "on", "target", "text", "origin", "agent", "at", "edited_at", "from_reply")

# Each state domain's name, as a person reads it (app/actions/domains.js DOMAINS).
.history_domain_names <- c(
  cells = "Open cells", mapping = "Data mapping", `field-overlays` = "Fields and labels", `gene-sets` = "Gene sets",
  views = "Views", `active-field` = "Active field", filters = "Filters", colors = "Colors", visibility = "Hidden cells",
  highlights = "Highlights", camera = "Camera", `camera-path` = "Camera path", settings = "Settings",
  analysis = "Analysis", layout = "Layout"
)

# --- Records ----------------------------------------------------------------------

# A record from its JSON (already checked), its fields named in snake case,
# lists kept as lists; a map record keeps its keys.
.history_build <- function(name, value) {
  if (is.null(value)) return(NULL)
  fields <- .history_fields[[name]]
  convert <- function(spec, item) {
    if (is.null(item)) return(NULL)
    switch(spec$type,
      record = .history_build(spec$name, item),
      array = lapply(item, function(entry) convert(spec$item, entry)),
      pair = list(convert(spec$items[[1L]], item[[1L]]), convert(spec$items[[2L]], item[[2L]])),
      item
    )
  }
  if (name %in% .history_map_records) {
    specs <- stats::setNames(lapply(fields, function(field) field$spec), vapply(fields, function(field) field$key, ""))
    out <- lapply(names(value), function(key) convert(specs[[key]], value[[key]]))
    return(stats::setNames(out, names(value)))
  }
  out <- list()
  for (field in fields) {
    if (!field$key %in% names(value)) next
    item <- convert(field$spec, value[[field$key]])
    out[field$attr] <- list(item)
  }
  out
}

.require_tracker <- function(tracker) {
  if (!inherits(tracker, "cellucid_tracker")) .abort("tracker must be what cellucid_tracker() returned.")
  tracker
}

.history_reading <- function(tracker) attr(tracker, "cellucid_reading", exact = TRUE)

.history_steps_by_id <- function(tracker) attr(tracker, "cellucid_steps_by_id", exact = TRUE)

# A tracker from its record as the model's JSON (checked): its steps, branches,
# groups, turns, notes and journal, and what it reads its values from.
.history_tracker <- function(record, reading) {
  tracker <- .history_build("Tracker", record)
  by_id <- new.env(parent = emptyenv())
  for (step in tracker$steps) assign(step$id, step, envir = by_id)
  structure(tracker, class = c("cellucid_tracker", "list"), cellucid_reading = reading, cellucid_steps_by_id = by_id)
}

cellucid_trackers <- function(session) {
  .require_session(session)
  trackers <- unname(session$trackers)
  data.frame(
    id = vapply(trackers, function(item) item$id, ""),
    name = vapply(trackers, function(item) item$name, ""),
    active = vapply(trackers, function(item) identical(item$id, session$active_tracker), TRUE),
    removed = vapply(trackers, function(item) !is.null(item$removed), TRUE),
    start = vapply(trackers, function(item) item$start$kind, ""),
    head = vapply(trackers, function(item) item$head_id, ""),
    active_branch = vapply(trackers, function(item) item$active_branch, ""),
    main_branch = vapply(trackers, function(item) item$main_branch, ""),
    # Derived from the tracker's `paused` (NULL while it records): whether it records, and since when it does not.
    recording = vapply(trackers, function(item) is.null(item$paused), TRUE),
    paused_at = .history_time(vapply(trackers, function(item) .history_ms(item$paused$at), 0)),
    steps = vapply(trackers, function(item) length(item$steps), 0L),
    branches = vapply(trackers, function(item) length(item$branches), 0L),
    groups = vapply(trackers, function(item) length(item$groups), 0L),
    created_at = .history_time(vapply(trackers, function(item) .history_ms(item$created_at), 0)),
    stringsAsFactors = FALSE
  )
}

cellucid_tracker <- function(session, tracker = NULL) {
  .require_session(session)
  if (is.null(tracker)) return(session$trackers[[session$active_tracker]])
  if (!.is_json_string(tracker)) .abort("tracker must be a tracker's name or id, or NULL.")
  found <- session$trackers[[tracker]]
  if (!is.null(found)) return(found)
  key <- .history_name_key(tracker)
  live <- Filter(function(item) is.null(item$removed), session$trackers)
  named <- Filter(function(item) identical(.history_name_key(item$name), key), live)
  if (length(named) == 1L) return(named[[1L]])
  .abort(.history_no_tracker(tracker, vapply(unname(live), function(item) item$name, "")))
}

cellucid_branches <- function(tracker) {
  .require_tracker(tracker)
  branches <- tracker$branches
  data.frame(
    id = vapply(branches, function(branch) branch$id, ""),
    name = vapply(branches, function(branch) branch$name, ""),
    base = vapply(branches, function(branch) branch$base %||% NA_character_, ""),
    tip = vapply(branches, function(branch) branch$tip, ""),
    main = vapply(branches, function(branch) identical(branch$id, tracker$main_branch), TRUE),
    locked = vapply(branches, function(branch) isTRUE(branch$locked), TRUE),
    archived = vapply(branches, function(branch) isTRUE(branch$archived), TRUE),
    n_notes = vapply(branches, function(branch) length(branch$notes), 0L),
    created_at = .history_time(vapply(branches, function(branch) .history_ms(branch$created_at), 0)),
    stringsAsFactors = FALSE
  )
}

cellucid_branch <- function(tracker, branch = NULL) {
  .require_tracker(tracker)
  if (is.null(branch)) branch <- tracker$active_branch
  if (!.is_json_string(branch)) .abort("branch must be a branch's name or id, or NULL.")
  found <- Filter(function(item) identical(item$id, branch), tracker$branches)
  if (length(found) == 0L) {
    key <- .history_name_key(branch)
    found <- Filter(function(item) identical(.history_name_key(item$name), key), tracker$branches)
  }
  if (length(found) == 0L) .abort(.history_no_branch(branch, vapply(tracker$branches, function(item) item$name, "")))
  c(found[[1L]], list(main = identical(found[[1L]]$id, tracker$main_branch)))
}

# The ids of a group's steps, in step order.
.history_group_steps <- function(tracker, group_id) {
  ids <- vapply(tracker$steps, function(step) step$id, "")
  ids[vapply(tracker$steps, function(step) identical(step$group, group_id), TRUE)]
}

cellucid_groups <- function(tracker) {
  .require_tracker(tracker)
  groups <- tracker$groups
  members <- lapply(groups, function(group) .history_group_steps(tracker, group$id))
  data.frame(
    id = vapply(groups, function(group) group$id, ""),
    name = vapply(groups, function(group) group$name, ""),
    branch = vapply(members, function(ids) cellucid_step(tracker, ids[[1L]])$branch, ""),
    first = vapply(members, function(ids) ids[[1L]], ""),
    last = vapply(members, function(ids) ids[[length(ids)]], ""),
    steps = vapply(members, length, 0L),
    n_notes = vapply(groups, function(group) length(group$notes), 0L),
    created_at = .history_time(vapply(groups, function(group) .history_ms(group$created_at), 0)),
    stringsAsFactors = FALSE
  )
}

cellucid_group <- function(tracker, group) {
  .require_tracker(tracker)
  if (!.is_json_string(group)) .abort("group must be a group's name or id.")
  found <- Filter(function(item) identical(item$id, group), tracker$groups)
  if (length(found) == 0L) {
    key <- .history_name_key(group)
    found <- Filter(function(item) identical(.history_name_key(item$name), key), tracker$groups)
  }
  if (length(found) == 0L) .abort(.history_no_group(group, vapply(tracker$groups, function(item) item$name, "")))
  c(found[[1L]], list(steps = .history_group_steps(tracker, found[[1L]]$id)))
}

cellucid_step <- function(tracker, step) {
  .require_tracker(tracker)
  found <- if (.is_json_string(step)) .history_steps_by_id(tracker)[[step]] else NULL
  if (is.null(found)) .abort(.history_no_step(as.character(step)[1L], tracker$name))
  found
}

# --- Paths, states, staleness ---------------------------------------------------------

.history_path <- function(tracker, step) {
  by_id <- .history_steps_by_id(tracker)
  path <- list()
  current <- cellucid_step(tracker, step)
  while (!is.null(current)) {
    path[[length(path) + 1L]] <- current
    current <- if (is.null(current$parent)) NULL else by_id[[current$parent]]
  }
  rev(path)
}

.history_hashes <- function(tracker, step) {
  hashes <- list()
  for (item in .history_path(tracker, step)) for (domain in names(item$after)) hashes[[domain]] <- item$after[[domain]]
  hashes
}

.history_raw_state <- function(tracker, step) {
  reading <- .history_reading(tracker)
  if (!reading$values_included) {
    .abort("This session was shared as its steps only, without the values they leave behind: ",
           "replay it in Cellucid to see its states.")
  }
  lapply(.history_hashes(tracker, step), reading$value)
}

cellucid_step_state <- function(tracker, step) {
  .require_tracker(tracker)
  .history_raw_state(tracker, step)
}

# The results stale at the head (history-v2 3.4): a step on the head's path
# whose output is still what the head shows while an input it read changed
# since, the first later step that changed it, and the inputs.
.history_stale <- function(tracker) {
  path <- .history_path(tracker, tracker$head_id)
  head <- .history_hashes(tracker, tracker$head_id)
  out <- list()
  for (index in seq_along(path)) {
    step <- path[[index]]
    read <- step$reads[setdiff(names(step$reads), "history")]
    if (length(read) == 0L || !identical(step$status, "applied")) next
    current <- all(vapply(unlist(step$touches), function(domain) identical(head[[domain]], step$after[[domain]]), TRUE))
    if (!current) next
    changed <- names(read)[vapply(names(read), function(domain) !identical(head[[domain]], read[[domain]]), TRUE)]
    if (length(changed) == 0L) next
    cause <- NULL
    for (later in path[-seq_len(index)]) {
      if (any(unlist(later$touches) %in% changed)) {
        cause <- later$id
        break
      }
    }
    out[[step$id]] <- list(cause = cause, domains = .history_domains[.history_domains %in% changed])
  }
  out
}

cellucid_stale <- function(tracker) {
  .require_tracker(tracker)
  stale <- .history_stale(tracker)
  unname(lapply(names(stale), function(id) c(list(step = id), stale[[id]])))
}

# --- Tables -------------------------------------------------------------------------

.history_ms <- function(value) if (is.null(value)) NA_real_ else as.double(value)

.history_time <- function(milliseconds) as.POSIXct(milliseconds / 1000, origin = "1970-01-01", tz = "UTC")

.history_text <- function(value) if (is.null(value)) NA_character_ else value

.history_step_table <- function(tracker, steps) {
  stale <- names(.history_stale(tracker))
  data.frame(
    id = vapply(steps, function(step) step$id, ""),
    branch = vapply(steps, function(step) step$branch, ""),
    parent = vapply(steps, function(step) .history_text(step$parent), ""),
    kind = vapply(steps, function(step) step$kind, ""),
    type = vapply(steps, function(step) .history_text(step$type), ""),
    v = vapply(steps, function(step) if (is.null(step$v)) NA_integer_ else as.integer(step$v), 0L),
    label = vapply(steps, function(step) step$label, ""),
    code = vapply(steps, function(step) .history_text(step$code), ""),
    origin = vapply(steps, function(step) step$by$origin, ""),
    agent = vapply(steps, function(step) .history_text(step$by$agent), ""),
    model = vapply(steps, function(step) .history_text(step$by$model), ""),
    host = vapply(steps, function(step) .history_text(step$by$host), ""),
    turn = vapply(steps, function(step) .history_text(step$by$turn), ""),
    at = .history_time(vapply(steps, function(step) .history_ms(step$at), 0)),
    status = vapply(steps, function(step) step$status, ""),
    stale = vapply(steps, function(step) step$id %in% stale, TRUE),
    bookmark = vapply(steps, function(step) .history_text(step$bookmark), ""),
    group = vapply(steps, function(step) .history_text(step$group), ""),
    n_notes = vapply(steps, function(step) length(step$notes), 0L),
    stringsAsFactors = FALSE
  )
}

cellucid_steps <- function(tracker, branch = NULL) {
  .require_tracker(tracker)
  if (is.null(branch)) return(.history_step_table(tracker, tracker$steps))
  .history_step_table(tracker, .history_path(tracker, cellucid_branch(tracker, branch)$tip))
}

cellucid_notes <- function(tracker) {
  .require_tracker(tracker)
  entries <- c(
    lapply(tracker$notes, function(note) list(on = "tracker", target = tracker$id, note = note)),
    unlist(lapply(tracker$branches, function(branch) lapply(branch$notes, function(note) list(on = "branch", target = branch$id, note = note))), recursive = FALSE),
    unlist(lapply(tracker$groups, function(group) lapply(group$notes, function(note) list(on = "group", target = group$id, note = note))), recursive = FALSE),
    unlist(lapply(tracker$steps, function(step) lapply(step$notes, function(note) list(on = "step", target = step$id, note = note))), recursive = FALSE),
    unlist(lapply(tracker$turns, function(turn) lapply(turn$notes, function(note) list(on = "turn", target = turn$id, note = note))), recursive = FALSE)
  )
  entries <- entries[order(vapply(entries, function(entry) .history_number(entry$note$id), 0))]
  data.frame(
    id = vapply(entries, function(entry) entry$note$id, ""),
    on = vapply(entries, function(entry) entry$on, ""),
    target = vapply(entries, function(entry) entry$target, ""),
    text = vapply(entries, function(entry) entry$note$text, ""),
    origin = vapply(entries, function(entry) entry$note$by$origin, ""),
    agent = vapply(entries, function(entry) .history_text(entry$note$by$agent), ""),
    at = .history_time(vapply(entries, function(entry) .history_ms(entry$note$at), 0)),
    edited_at = .history_time(vapply(entries, function(entry) .history_ms(entry$note$edited_at), 0)),
    from_reply = vapply(entries, function(entry) isTRUE(entry$note$from_reply), TRUE),
    stringsAsFactors = FALSE
  )
}

cellucid_turns <- function(tracker) {
  .require_tracker(tracker)
  turns <- tracker$turns
  steps <- tracker$steps
  data.frame(
    id = vapply(turns, function(turn) turn$id, ""),
    label = vapply(turns, function(turn) turn$label, ""),
    origin = vapply(turns, function(turn) turn$by$origin, ""),
    agent = vapply(turns, function(turn) .history_text(turn$by$agent), ""),
    model = vapply(turns, function(turn) .history_text(turn$by$model), ""),
    host = vapply(turns, function(turn) .history_text(turn$by$host), ""),
    branch = vapply(turns, function(turn) turn$branch, ""),
    conversation = vapply(turns, function(turn) .history_text(turn$conversation), ""),
    message = vapply(turns, function(turn) .history_text(turn$message), ""),
    started_at = .history_time(vapply(turns, function(turn) .history_ms(turn$started_at), 0)),
    ended_at = .history_time(vapply(turns, function(turn) .history_ms(turn$ended_at), 0)),
    steps = vapply(turns, function(turn) sum(vapply(steps, function(step) identical(step$by$turn, turn$id), TRUE)), 0L),
    outside = vapply(turns, function(turn) length(turn$outside), 0L),
    summary = vapply(turns, function(turn) if (length(turn$notes) > 0L) turn$notes[[1L]]$text else NA_character_, ""),
    stringsAsFactors = FALSE
  )
}

cellucid_journal <- function(tracker) {
  .require_tracker(tracker)
  entries <- tracker$journal
  data.frame(
    id = vapply(entries, function(entry) entry$id, ""),
    op = vapply(entries, function(entry) entry$op, ""),
    at = .history_time(vapply(entries, function(entry) .history_ms(entry$at), 0)),
    origin = vapply(entries, function(entry) entry$by$origin, ""),
    label = vapply(entries, function(entry) entry$label, ""),
    purged = vapply(entries, function(entry) isTRUE(entry$purged), TRUE),
    stringsAsFactors = FALSE
  )
}

# --- What differs between two steps -------------------------------------------------

.history_format_count <- function(value) formatC(as.double(value), format = "d", big.mark = ",")

.history_count_of <- function(value, singular, plural = paste0(singular, "s")) {
  paste(.history_format_count(value), if (value == 1) singular else plural)
}

.history_json <- function(value) .js_json(.history_unarray(value))

# A value as the web app's JSON.stringify sees it in a comparison (typed arrays by their bytes).
.history_unarray <- function(value) {
  if (is.list(value)) {
    out <- lapply(value, .history_unarray)
    if (!is.null(names(value)) && length(out) == 0L) out <- structure(list(), names = character())
    return(out)
  }
  if (!is.null(attr(value, "cellucid_array"))) return(I(as.double(value)))
  value
}

.history_same_value <- function(a, b) identical(.history_unarray_exact(a), .history_unarray_exact(b))

.history_unarray_exact <- function(value) {
  if (is.list(value)) return(lapply(value, .history_unarray_exact))
  value
}

.history_label_text <- function(value) {
  if (is.logical(value)) return(if (isTRUE(value)) "true" else "false")
  if (is.numeric(value)) return(.js_numbers(value))
  as.character(value)
}

.history_describe_cells <- function(a, b) {
  opened <- function(open) {
    if (is.null(open)) return("every cell")
    paste0(.history_count_of(open$k, "cell"), " (", open$method, ", seed ", open$seed, ")")
  }
  genes <- function(count) if (is.null(count)) "the opening's sample" else .history_count_of(count, "cell")
  lines <- character()
  if (!identical(.history_json(a$open), .history_json(b$open))) lines <- c(lines, paste0("Cells open: ", opened(a$open), " \u2192 ", opened(b$open)))
  if (!identical(a$genes, b$genes)) lines <- c(lines, paste0("Genes read for: ", genes(a$genes), " \u2192 ", genes(b$genes)))
  lines
}

.history_describe_user_field <- function(before, after) {
  if (is.null(before)) return(paste0(if (identical(after$kind, "category")) "Annotation" else "Field", " \"", after$key, "\" made"))
  if (is.null(after)) return(paste0("Field \"", before$key, "\" removed"))
  lines <- character()
  if (!identical(before$key, after$key)) lines <- c(lines, paste0("Field \"", before$key, "\" renamed \"", after$key, "\""))
  if (!identical(before$deleted, after$deleted)) lines <- c(lines, paste0("Field \"", after$key, "\" ", if (isTRUE(after$deleted)) "deleted" else "restored"))
  if (!identical(before$purged, after$purged) && isTRUE(after$purged)) lines <- c(lines, paste0("Field \"", after$key, "\" deleted for good"))
  if (!.history_same_value(before$categories, after$categories)) {
    lines <- c(lines, paste0("Field \"", after$key, "\": categories changed"))
  } else if (!.history_same_value(before$codes, after$codes) || !.history_same_value(before$values, after$values)) {
    lines <- c(lines, paste0("Field \"", after$key, "\": values changed"))
  }
  lines
}

.history_describe_field_overlays <- function(a, b) {
  lines <- character()
  field_ids <- function(entries) stats::setNames(vapply(entries, function(entry) entry$name, ""), vapply(entries, function(entry) paste0(entry$source, ":", entry$key), ""))
  renamed_a <- field_ids(a$renamedFields)
  renamed_b <- field_ids(b$renamedFields)
  strip <- function(id) sub("^[^:]*:", "", id)
  for (id in names(renamed_b)) if (!identical(unname(renamed_a[id]), unname(renamed_b[id]))) lines <- c(lines, paste0("Field \"", strip(id), "\" shown as \"", renamed_b[[id]], "\""))
  for (id in names(renamed_a)) if (!id %in% names(renamed_b)) lines <- c(lines, paste0("Field \"", strip(id), "\" shown under its own name"))
  category_id <- function(entry) .js_json(list(entry$source, entry$field, entry$category))
  categories_a <- stats::setNames(a$renamedCategories, vapply(a$renamedCategories, category_id, ""))
  categories_b <- stats::setNames(b$renamedCategories, vapply(b$renamedCategories, category_id, ""))
  for (id in names(categories_b)) {
    entry <- categories_b[[id]]
    if (!identical(categories_a[[id]]$name, entry$name)) {
      lines <- c(lines, paste0(entry$field, ": \"", .history_label_text(entry$category), "\" shown as \"", .history_label_text(entry$name), "\""))
    }
  }
  for (id in names(categories_a)) {
    if (!id %in% names(categories_b)) lines <- c(lines, paste0(categories_a[[id]]$field, ": \"", .history_label_text(categories_a[[id]]$category), "\" shown under its own name"))
  }
  deleted_ids <- function(entries) vapply(entries, function(entry) paste0(entry$source, ":", entry$key), "")
  deleted_a <- stats::setNames(a$deletedFields, deleted_ids(a$deletedFields))
  deleted_b <- stats::setNames(b$deletedFields, deleted_ids(b$deletedFields))
  for (id in names(deleted_b)) if (!id %in% names(deleted_a)) lines <- c(lines, paste0("Field \"", deleted_b[[id]]$key, "\" deleted"))
  for (id in names(deleted_a)) if (!id %in% names(deleted_b)) lines <- c(lines, paste0("Field \"", deleted_a[[id]]$key, "\" restored"))
  user_a <- stats::setNames(a$userFields, vapply(a$userFields, function(item) item$id, ""))
  user_b <- stats::setNames(b$userFields, vapply(b$userFields, function(item) item$id, ""))
  for (id in unique(c(names(user_a), names(user_b)))) {
    before <- user_a[[id]]
    after <- user_b[[id]]
    if (!.history_same_value(before, after)) lines <- c(lines, .history_describe_user_field(before, after))
  }
  lines
}

.history_describe_gene_sets <- function(a, b) {
  before <- stats::setNames(a$sets, vapply(a$sets, function(entry) entry$name, ""))
  after <- stats::setNames(b$sets, vapply(b$sets, function(entry) entry$name, ""))
  lines <- character()
  for (name in names(after)) {
    was <- before[[name]]
    count <- .history_count_of(length(after[[name]]$genes), "gene")
    if (is.null(was)) {
      lines <- c(lines, paste0("Gene set ", name, " added (", count, ")"))
    } else if (!identical(.history_json(was), .history_json(after[[name]]))) {
      lines <- c(lines, paste0("Gene set ", name, " changed (", count, ")"))
    }
  }
  for (name in names(before)) if (!name %in% names(after)) lines <- c(lines, paste0("Gene set ", name, " removed"))
  lines
}

.history_describe_highlights <- function(a, b) {
  lines <- character()
  pages_a <- stats::setNames(a$pages, vapply(a$pages, function(page) page$id, ""))
  pages_b <- stats::setNames(b$pages, vapply(b$pages, function(page) page$id, ""))
  for (id in names(pages_b)) {
    page <- pages_b[[id]]
    was <- pages_a[[id]]
    if (is.null(was)) {
      lines <- c(lines, paste0("Highlight page \"", page$name, "\" added"))
      next
    }
    if (!identical(was$name, page$name)) lines <- c(lines, paste0("Highlight page \"", was$name, "\" renamed \"", page$name, "\""))
    if (!identical(was$color, page$color)) lines <- c(lines, paste0("Highlight page \"", page$name, "\" recoloured"))
    groups_a <- stats::setNames(was$groups, vapply(was$groups, function(group) group$id, ""))
    groups_b <- stats::setNames(page$groups, vapply(page$groups, function(group) group$id, ""))
    for (group_id in names(groups_b)) {
      group <- groups_b[[group_id]]
      before <- groups_a[[group_id]]
      count <- .history_count_of(length(group$cells), "cell")
      if (is.null(before)) {
        lines <- c(lines, paste0("\"", group$label, "\" highlighted (", count, ") on \"", page$name, "\""))
      } else if (!identical(before$enabled, group$enabled)) {
        lines <- c(lines, paste0("\"", group$label, "\" ", if (isTRUE(group$enabled)) "shown" else "hidden", " on \"", page$name, "\""))
      } else if (!identical(before$label, group$label)) {
        lines <- c(lines, paste0("\"", before$label, "\" renamed \"", group$label, "\""))
      } else if (length(before$cells) != length(group$cells)) {
        lines <- c(lines, paste0("\"", group$label, "\" now ", count))
      }
    }
    for (group_id in names(groups_a)) if (!group_id %in% names(groups_b)) lines <- c(lines, paste0("\"", groups_a[[group_id]]$label, "\" removed from \"", page$name, "\""))
  }
  for (id in names(pages_a)) if (!id %in% names(pages_b)) lines <- c(lines, paste0("Highlight page \"", pages_a[[id]]$name, "\" removed"))
  if (!identical(a$active, b$active) && .is_json_string(b$active) && b$active %in% names(pages_b)) {
    lines <- c(lines, paste0("Highlight page \"", pages_b[[b$active]]$name, "\" shown"))
  }
  lines
}

.history_describe <- function(domain, a, b) {
  lines <- switch(domain,
    cells = .history_describe_cells(a, b),
    `field-overlays` = .history_describe_field_overlays(a, b),
    `gene-sets` = .history_describe_gene_sets(a, b),
    highlights = .history_describe_highlights(a, b),
    character()
  )
  if (length(lines) == 0L) paste(.history_domain_names[[domain]], "changed") else lines
}

cellucid_step_diff <- function(tracker, a, b, other = NULL) {
  .require_tracker(tracker)
  other <- if (is.null(other)) tracker else .require_tracker(other)
  values_a <- .history_hashes(tracker, cellucid_step(tracker, a)$id)
  values_b <- .history_hashes(other, cellucid_step(other, b)$id)
  reading_a <- .history_reading(tracker)
  reading_b <- .history_reading(other)
  rows <- list()
  for (domain in .history_domains) {
    if (identical(values_a[[domain]], values_b[[domain]])) next
    lines <- if (reading_a$values_included && reading_b$values_included) {
      .history_describe(domain, reading_a$value(values_a[[domain]]), reading_b$value(values_b[[domain]]))
    } else {
      paste(.history_domain_names[[domain]], "changed")
    }
    for (line in lines) rows[[length(rows) + 1L]] <- c(domain = domain, name = .history_domain_names[[domain]], line = line)
  }
  out <- data.frame(
    domain = vapply(rows, function(row) row[["domain"]], ""),
    name = vapply(rows, function(row) row[["name"]], ""),
    line = vapply(rows, function(row) row[["line"]], ""),
    stringsAsFactors = FALSE
  )
  # The latest step both paths share, within one tracker; none across two.
  common <- NULL
  if (identical(tracker$id, other$id)) {
    on_a <- vapply(.history_path(tracker, a), function(step) step$id, "")
    on_b <- vapply(.history_path(other, b), function(step) step$id, "")
    common <- utils::tail(on_b[on_b %in% on_a], 1L)
  }
  attr(out, "common") <- common
  out
}

# --- Who took a step, the methods paragraph, the pseudocode ---------------------------

.history_method_domains <- c("mapping", "field-overlays", "gene-sets", "active-field", "filters", "visibility", "highlights", "analysis")

# Who acted (a By record: step$by, note$by): a key for its look and a short name (originOf).
.history_origin_of <- function(by) {
  switch(by$origin,
    outside = c("outside", "Not recorded"),
    user = c("person", "Person"),
    agent = {
      name <- by$agent %||% "Agent"
      c(if (tolower(name) %in% c("claude", "chatgpt")) tolower(name) else "agent", name)
    },
    host = {
      host <- tolower(by$host %||% "")
      if (host %in% c("python", "r")) c(host, if (identical(host, "python")) "Python" else "R") else c("host", "Notebook")
    },
    session = c("session", "Session"),
    replay = c("replay", "Replay"),
    c("start", "Start")
  )
}

.history_origin_phrase <- function(by) {
  origin <- .history_origin_of(by)
  switch(origin[[1L]],
    outside = NULL,
    person = "the user",
    claude = , chatgpt = , agent = paste0("an AI assistant (", origin[[2L]],
                                          if (!is.null(by$model) && nzchar(by$model)) paste0(", ", by$model) else "",
                                          ") at the user's request"),
    python = , r = paste0("a ", origin[[2L]], " notebook through the Cellucid viewer API"),
    host = "a notebook through the Cellucid viewer API",
    session = "loading a session file",
    replay = "replaying a recorded history",
    "opening the dataset"
  )
}

# Everyone who acted in a step: a squashed step's squashed$by when it names anyone, else its by.
.history_actors_of <- function(step) {
  if (identical(step$kind, "squashed") && length(step$squashed$by) > 0L) step$squashed$by else list(step$by)
}

# What a squashed step stands for, as a sentence opening (squashPhrase).
.history_squash_phrase <- function(squashed) {
  switch(squashed$cause,
    paused = "Changes made with recording paused",
    steps = paste(.history_count_of(squashed$count, "step"), "merged into one"),
    clear = "A fresh start from the state reached before the history was cleared",
    share = "A shared state, received without its history",
    .abort("A squashed step is made for paused, steps, clear or share, not ", squashed$cause)
  )
}

# Consecutive steps that share a group, in order (groupRuns): list(group = <Group or NULL>, steps = list()).
.history_group_runs <- function(steps, groups) {
  by_id <- stats::setNames(groups, vapply(groups, function(group) group$id, ""))
  runs <- list()
  for (step in steps) {
    last <- length(runs)
    if (last > 0L && identical(runs[[last]]$group$id, step$group)) {
      runs[[last]]$steps[[length(runs[[last]]$steps) + 1L]] <- step
      next
    }
    if (!is.null(step$group) && is.null(by_id[[step$group]])) {
      .abort("Step ", substring(step$id, 2L), " names the group \"", step$group, "\", which the tracker does not hold")
    }
    runs[[last + 1L]] <- list(group = if (is.null(step$group)) NULL else by_id[[step$group]], steps = list(step))
  }
  runs
}

# Each run's heading (runHeadings): its group's name, or "After <group>" / "Before <group>" for
# steps in no group; NA on a path with no group.
.history_run_headings <- function(runs) {
  vapply(seq_along(runs), function(index) {
    group <- runs[[index]]$group
    if (!is.null(group)) return(group$name)
    if (index > 1L) return(paste("After", runs[[index - 1L]]$group$name))
    if (length(runs) > 1L) paste("Before", runs[[2L]]$group$name) else NA_character_
  }, "")
}

.history_trim_end <- function(text) sub("[.;,[:space:]]+$", "", trimws(text))

.history_in_sentence <- function(text) {
  if (nchar(text) < 2L) return(text)
  first <- strsplit(text, "[[:space:]]")[[1L]][[1L]]
  if (nchar(first) > 1L && identical(first, toupper(first))) return(text)
  paste0(tolower(substring(text, 1L, 1L)), substring(text, 2L))
}

.history_number_ranges <- function(numbers) {
  numbers <- sort(unique(numbers))
  runs <- list()
  for (number in numbers) {
    last <- length(runs)
    if (last > 0L && number == runs[[last]][[2L]] + 1) runs[[last]][[2L]] <- number else runs[[last + 1L]] <- c(number, number)
  }
  parts <- vapply(runs, function(run) if (run[[1L]] == run[[2L]]) as.character(run[[1L]]) else paste0(run[[1L]], "\u2013", run[[2L]]), "")
  if (length(parts) <= 1L) return(paste(parts, collapse = ""))
  paste0(paste(parts[-length(parts)], collapse = ", "), " and ", parts[[length(parts)]])
}

.history_subsample_methods <- c(
  blocks = "in blocks of neighbouring cells", uniform = "uniformly at random", stratified = "stratified by a field",
  spatial = "spatially"
)

.history_dataset_sentence <- function(dataset) {
  name <- if (!is.null(dataset$name)) paste0("\u201c", .history_trim_end(dataset$name), "\u201d") else NULL
  size <- c(if (!is.null(dataset$cells)) .history_count_of(dataset$cells, "cell"), if (!is.null(dataset$genes)) .history_count_of(dataset$genes, "gene"))
  facts <- paste(c(if (length(size) > 0L) paste(size, collapse = " and ")), collapse = "; ")
  subject <- if (is.null(name)) "The dataset" else paste("The dataset", name)
  how <- "was opened in full"
  open <- dataset$open_cells
  if (!is.null(open) && open$open < open$of) {
    percent <- 100 * open$open / open$of
    share <- if (percent >= 1) paste0(floor(percent + 0.5), "%") else paste0(.js_numbers(signif(percent, 2)), "%")
    method <- .history_subsample_methods[open$method]
    seed <- if (!is.null(open$seed)) paste0(" with seed ", open$seed) else ""
    how <- paste0("was opened on ", .history_format_count(open$open), " of its ", .history_count_of(open$of, "cell"), " (", share, ")",
                  if (!is.na(method)) paste0(", chosen ", method, seed) else seed)
  }
  paste0(subject, if (nzchar(facts)) paste0(" (", facts, ")") else "", " ", how, ".")
}

.history_lower_domains <- function(domains) paste(tolower(.history_domain_names[domains]), collapse = ", ")

.history_step_item <- function(step) {
  touches <- unlist(step$touches)
  if (identical(step$kind, "outside")) {
    return(paste0(.history_lower_domains(intersect(touches, .history_method_domains)), " changed outside the recorded steps"))
  }
  if (identical(step$kind, "squashed")) {
    return(paste0(.history_in_sentence(.history_squash_phrase(step$squashed)),
                  if (length(touches) > 0L) paste0(" (", .history_lower_domains(touches), ")") else ""))
  }
  .history_in_sentence(.history_trim_end(step$label))
}

# "Step 2 was taken by ..., and steps 4-6 by ...", or NULL when the person took them all.
.history_origin_sentence <- function(steps, numbers) {
  phrases <- character()
  by_phrase <- list()
  for (index in seq_along(steps)) {
    for (by in .history_actors_of(steps[[index]])) {
      if (identical(.history_origin_of(by)[[1L]], "person")) next
      phrase <- .history_origin_phrase(by)
      if (is.null(phrase)) next
      if (!phrase %in% phrases) {
        phrases <- c(phrases, phrase)
        by_phrase[[phrase]] <- integer()
      }
      if (!numbers[[index]] %in% by_phrase[[phrase]]) by_phrase[[phrase]] <- c(by_phrase[[phrase]], numbers[[index]])
    }
  }
  if (length(phrases) == 0L) return(NULL)
  clauses <- vapply(seq_along(phrases), function(position) {
    phrase <- phrases[[position]]
    these <- by_phrase[[phrase]]
    noun <- if (length(these) == 1L) "step" else "steps"
    if (position > 1L) return(paste0(noun, " ", .history_number_ranges(these), " by ", phrase))
    paste0(toupper(substring(noun, 1L, 1L)), substring(noun, 2L), " ", .history_number_ranges(these), " ",
           if (length(these) == 1L) "was" else "were", " taken by ", phrase)
  }, "")
  if (length(clauses) == 1L) return(paste0(clauses, "."))
  paste0(paste(clauses[-length(clauses)], collapse = ", "), ", and ", clauses[[length(clauses)]], ".")
}

.history_list_sentences <- function(steps, numbers) {
  items <- vapply(seq_along(steps), function(index) paste0("(", numbers[[index]], ") ", .history_step_item(steps[[index]])), "")
  n <- length(steps)
  c(
    paste0("The following ", .history_count_of(n, "step"), " ", if (n == 1L) "was" else "were", " applied",
           if (n == 1L) "" else ", in this order", ": ", paste(items, collapse = "; "), "."),
    .history_origin_sentence(steps, numbers)
  )
}

# A branch's path of steps as methods text, as the web app writes it (app/history/exports/methods.js).
.history_methods_paragraph <- function(steps, dataset, version, notes, groups = list()) {
  taken <- Filter(function(step) !identical(step$kind, "root"), steps)
  is_flagged <- function(step) identical(step$kind, "action") && identical(step$status, "unapplicable")
  flagged <- Filter(is_flagged, taken)
  applied <- Filter(function(step) !is_flagged(step) && length(step$touches) > 0L, taken)
  listed <- Filter(function(step) any(unlist(step$touches) %in% .history_method_domains), applied)
  presentation <- length(applied) - length(listed)
  opening <- c(
    paste0("Single-cell data were explored interactively in Cellucid", if (!is.null(version)) paste0(" (version ", version, ")"), "."),
    .history_dataset_sentence(dataset)
  )
  if (length(listed) == 0L) opening <- c(opening, "No recorded step changed the data, its annotation or what was computed.")
  closing <- character()
  squashed <- which(vapply(listed, function(step) identical(step$kind, "squashed"), TRUE))
  if (length(squashed) == 1L) {
    closing <- c(closing, paste0("For step ", squashed, ", only the state it reached was kept, not the operations that led to it."))
  } else if (length(squashed) > 1L) {
    closing <- c(closing, paste0("For steps ", .history_number_ranges(squashed), ", only the states they reached were kept, not the operations that led to them."))
  }
  if (presentation > 0L) {
    closing <- c(closing, paste0(.history_count_of(presentation, "further step"), " only changed the camera, the panel layout, colors or display settings and ",
                                 if (presentation == 1L) "is" else "are", " not listed."))
  }
  if (length(flagged) > 0L) {
    closing <- c(closing, paste0(.history_count_of(length(flagged), "recorded step"), " no longer applied after an earlier step was removed and ",
                                 if (length(flagged) == 1L) "is" else "are", " not included."))
  }
  closing <- c(closing, "The complete history, with every parameter, can be reopened from the accompanying Cellucid session file.")

  number_of <- stats::setNames(seq_along(listed), vapply(listed, function(step) step$id, ""))
  runs <- .history_group_runs(listed, groups)
  headings <- .history_run_headings(runs)
  sentences_of <- function(run) .history_list_sentences(run$steps, unname(number_of[vapply(run$steps, function(step) step$id, "")]))
  paragraphs <- if (all(is.na(headings))) {
    paste(c(opening, unlist(lapply(runs, sentences_of)), closing), collapse = " ")
  } else {
    summary <- if (length(listed) == 1L) "The step applied is listed below." else paste0("The ", .history_count_of(length(listed), "step"), " applied are listed below, in order.")
    c(paste(c(opening, summary, closing), collapse = " "),
      vapply(seq_along(runs), function(index) paste0(headings[[index]], "\n", paste(sentences_of(runs[[index]]), collapse = " ")), ""))
  }
  paste(c(if (length(notes) > 0L) notes[[1L]]$text, paragraphs), collapse = "\n\n")
}

cellucid_methods <- function(tracker, branch = NULL) {
  .require_tracker(tracker)
  target <- cellucid_branch(tracker, branch %||% tracker$main_branch)
  dataset <- .history_reading(tracker)$dataset
  .history_methods_paragraph(.history_path(tracker, target$tip), dataset, dataset$version, tracker$notes, tracker$groups)
}

.history_comment <- function(text) gsub("[[:space:]]+", " ", trimws(text))

.history_pseudocode_line <- function(step) {
  names <- paste(.history_domain_names[unlist(step$touches)], collapse = ", ")
  if (identical(step$kind, "outside")) return(paste0("# Changed outside the history", if (nzchar(names)) paste0(": ", names) else ""))
  if (identical(step$kind, "squashed")) {
    return(paste0("# ", .history_squash_phrase(step$squashed), if (nzchar(names)) paste0(": ", names) else "", " (values only)"))
  }
  if (identical(step$status, "unapplicable")) return(paste0("# No longer applies (", .history_comment(step$reason %||% "not applied"), "): ", step$code))
  origin <- .history_origin_of(step$by)
  if (identical(origin[[1L]], "person")) step$code else paste0(step$code, "  # ", origin[[2L]])
}

# A branch's path of steps as pseudocode, by app/history/exports/pseudocode.js's rules.
.history_pseudocode_text <- function(steps, dataset_name, version, tracker, branch, groups = list()) {
  taken <- Filter(function(step) !identical(step$kind, "root"), steps)
  app <- if (!is.null(version)) paste("Cellucid", version) else "Cellucid"
  lines <- c(
    paste0("# ", app, if (!is.null(dataset_name)) paste0(" \u00b7 ", .history_comment(dataset_name)) else ""),
    paste0("# Tracker \"", .history_comment(tracker), "\", branch \"", .history_comment(branch), "\": ",
           .history_count_of(length(taken), "step"))
  )
  for (run in .history_group_runs(taken, groups)) {
    name <- if (is.null(run$group)) NULL else .history_comment(run$group$name)
    if (!is.null(name)) lines <- c(lines, paste0("# --- ", name, " ---"))
    for (step in run$steps) {
      lines <- c(lines, .history_pseudocode_line(step))
      if (!is.null(step$bookmark)) lines <- c(lines, paste0("# Bookmark: ", .history_comment(step$bookmark)))
      if (identical(step$kind, "squashed")) {
        for (bookmark in unlist(step$squashed$carried$bookmarks)) {
          lines <- c(lines, paste0("# Bookmark: ", .history_comment(bookmark), " (on a step merged into this one)"))
        }
      }
    }
    if (!is.null(name)) lines <- c(lines, paste0("# --- end of ", name, " ---"))
  }
  paste0(paste(lines, collapse = "\n"), "\n")
}

cellucid_pseudocode <- function(tracker, branch = NULL) {
  .require_tracker(tracker)
  target <- cellucid_branch(tracker, branch %||% tracker$main_branch)
  dataset <- .history_reading(tracker)$dataset
  .history_pseudocode_text(.history_path(tracker, target$tip), dataset$name, dataset$version, tracker$name, target$name, tracker$groups)
}

print.cellucid_tracker <- function(x, ...) {
  cat("<Cellucid tracker> ", x$name, " (", x$id, ")", if (!is.null(x$removed)) " [removed]", "\n",
      "  steps:    ", length(x$steps), ", head ", x$head_id, "\n",
      "  branches: ", paste(vapply(x$branches, function(branch) branch$name, ""), collapse = ", "),
      " (active ", x$active_branch, ", main ", x$main_branch, ")\n", sep = "")
  invisible(x)
}
