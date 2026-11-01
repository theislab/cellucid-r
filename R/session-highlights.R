# The obs columns that carry a session's highlight pages: one TRUE/FALSE
# column per highlight group, from the cell set its chunk holds
# (src/session_codecs.c). A category group may be stored relative to its
# category's cells: those are read from the view's obs, and a column whose
# category no longer holds the cells the session recorded is refused.
# cellucid-python's anndata_session/_highlights.py, rule for rule.

.session_check_keys <- function(value, expected, label) {
  if (!.is_json_object(value)) .abort(label, " must be an object")
  if (!setequal(names(value), expected) || anyDuplicated(names(value))) {
    .abort(label, " has invalid fields (", .session_key_report(names(value), expected), ")")
  }
  value
}

.session_text <- function(value, label) {
  if (!.is_json_string(value) || !nzchar(value)) .abort(label, " must be a non-empty string")
  value
}

.session_identity <- function(value, prefix, label) {
  .session_text(value, label)
  if (!grepl(paste0("^", prefix, "[1-9][0-9]*$"), value)) {
    .abort(label, " must use the current ", prefix, "<positive integer> identity")
  }
  value
}

.session_category_value <- function(value, label) {
  ok <- (is.character(value) || is.logical(value) || is.numeric(value)) && length(value) == 1L && !is.na(value)
  if (!ok || (is.numeric(value) && !is.finite(value))) {
    .abort(label, " must be a string, finite number, or boolean")
  }
  value
}

.highlight_group_keys <- c("id", "type", "label", "enabled", "cellCount")

.validate_highlight_group <- function(group, label) {
  type <- group[["type"]]
  if (!.is_json_string(type) || !type %in% c("annotation", "category", "combined", "knn", "lasso", "proximity", "range")) {
    .abort(label, ".type is unsupported")
  }
  expected <- switch(type,
    category = c(.highlight_group_keys, "fieldKey", "fieldSource", "categoryName"),
    range = c(.highlight_group_keys, "fieldKey", "fieldSource", "rangeMin", "rangeMax"),
    .highlight_group_keys
  )
  .session_check_keys(group, expected, label)
  .session_identity(group[["id"]], "highlight_", paste0(label, ".id"))
  .session_text(group[["label"]], paste0(label, ".label"))
  if (!is.logical(group[["enabled"]]) || length(group[["enabled"]]) != 1L || is.na(group[["enabled"]])) {
    .abort(label, ".enabled must be exactly True or False")
  }
  count <- group[["cellCount"]]
  if (!.is_json_integer(count)) .abort(label, ".cellCount must be an integer")
  if (count < 1) .abort(label, ".cellCount must be positive")
  if (type %in% c("category", "range")) {
    .session_text(group[["fieldKey"]], paste0(label, ".fieldKey"))
    if (!.is_json_string(group[["fieldSource"]]) || !group[["fieldSource"]] %in% c("obs", "var")) {
      .abort(label, ".fieldSource must be 'obs' or 'var'")
    }
  }
  if (identical(type, "category")) .session_category_value(group[["categoryName"]], paste0(label, ".categoryName"))
  if (identical(type, "range")) {
    low <- group[["rangeMin"]]
    high <- group[["rangeMax"]]
    for (bound in list(list(low, "rangeMin"), list(high, "rangeMax"))) {
      if (!is.numeric(bound[[1L]]) || length(bound[[1L]]) != 1L || !is.finite(bound[[1L]])) {
        .abort(label, ".", bound[[2L]], " must be a finite number")
      }
    }
    if (low > high) .abort(label, ".rangeMin must not exceed rangeMax")
  }
  invisible(group)
}

# The 0-based cells of a category group's category, ascending.
.category_reference <- function(obs, group, label) {
  if (!identical(group[["type"]], "category") || !identical(group[["fieldSource"]], "obs")) {
    .abort(label, " is stored relative to a category, and it is not a category group on an obs field")
  }
  key <- group[["fieldKey"]]
  if (sum(names(obs) == key) != 1L) {
    .abort(label, " is stored relative to the category '", group[["categoryName"]], "' of obs field '", key,
           "', which must exist exactly once in the obs")
  }
  values <- obs[[key]]
  category <- group[["categoryName"]]
  members <- if (is.logical(category)) {
    is.logical(values) & !is.na(values) & values == category
  } else if (is.numeric(category)) {
    is.numeric(values) & !is.na(values) & values == category
  } else {
    !is.na(values) & as.character(values) == category & (is.factor(values) | is.character(values))
  }
  as.integer(which(members) - 1L)
}

# Returns list(plans, uns, meta): one planned column per group (when
# `materialize`), the highlights record for the object's metadata, and the
# decoded highlights/meta chunk.
.plan_highlights <- function(bundle, chunks, view, prefix, existing, materialize) {
  ids <- names(chunks)
  memberships <- ids[startsWith(ids, .highlight_cells_prefix)]
  if (!"highlights/meta" %in% ids) {
    if (length(memberships) > 0L) .abort("Session contains highlight membership chunks without highlights/meta")
    return(list(plans = list(), uns = NULL, meta = NULL, existing = existing))
  }
  meta <- .session_check_keys(.session_chunk(bundle, "highlights/meta"), c("pages", "activePageId"), "highlights/meta")
  pages <- meta[["pages"]]
  if (!.is_json_array(pages)) .abort("highlights/meta.pages must be an array")
  if (length(pages) == 0L) .abort("highlights/meta.pages must be a non-empty array")
  active <- .session_identity(meta[["activePageId"]], "page_", "highlights/meta.activePageId")
  plans <- list()
  groups_uns <- list()
  page_ids <- character()
  group_ids <- character()
  for (page_index in seq_along(pages)) {
    page_label <- paste0("highlights/meta.pages[", page_index - 1L, "]")
    page <- .session_check_keys(pages[[page_index]], c("id", "name", "color", "highlightedGroups"), page_label)
    page_id <- .session_identity(page[["id"]], "page_", paste0(page_label, ".id"))
    if (page_id %in% page_ids) .abort("Duplicate highlight page id: ", page_id)
    page_ids <- c(page_ids, page_id)
    page_name <- .session_text(page[["name"]], paste0(page_label, ".name"))
    color <- .session_text(page[["color"]], paste0(page_label, ".color"))
    if (!grepl("^#[0-9A-Fa-f]{6}$", color)) .abort(page_label, ".color must be a six-digit hex color")
    groups <- page[["highlightedGroups"]]
    if (!.is_json_array(groups)) .abort(page_label, ".highlightedGroups must be an array")
    for (group_index in seq_along(groups)) {
      group_label <- paste0(page_label, ".highlightedGroups[", group_index - 1L, "]")
      group <- groups[[group_index]]
      if (!.is_json_object(group)) .abort(group_label, " must be an object")
      .validate_highlight_group(group, group_label)
      group_id <- group[["id"]]
      if (group_id %in% group_ids) .abort("Duplicate highlight group id: ", group_id)
      group_ids <- c(group_ids, group_id)
      chunk_id <- paste0(.highlight_cells_prefix, group_id)
      if (!chunk_id %in% ids) {
        .abort("Highlight group '", group_id, "' is missing required chunk '", chunk_id, "'")
      }
      if (!identical(chunks[[chunk_id]][["label"]], paste0("Highlight cells: ", group[["label"]]))) {
        .abort(chunk_id, " label must match highlight group '", group_id, "'")
      }
      payload <- .session_chunk(bundle, chunk_id)
      reference <- if (length(payload) > 0L && as.integer(payload[[1L]]) %in% 3:4) {
        .category_reference(view$obs, group, group_label)
      }
      cells <- .native_decode_cell_set(payload, view$n_cells, reference)
      if (length(cells) != group[["cellCount"]]) {
        .abort(chunk_id, " contains ", length(cells), " ", .plural(length(cells), "index", "indices"), ", but metadata declares ", group[["cellCount"]])
      }
      column <- paste0(prefix, group_id)
      if (materialize) {
        if (column %in% existing) .abort("Column already exists: ", column)
        existing <- c(existing, column)
        plans[[length(plans) + 1L]] <- list(
          name = column, values = .membership_column(cells, NULL, view$n_cells),
          metadata = list(kind = "highlight", group_id = group_id, page_id = page_id)
        )
      }
      groups_uns[[group_id]] <- list(
        obs_column = if (materialize) column else NULL, page_id = page_id, page_name = page_name, group = group
      )
    }
  }
  if (!active %in% page_ids) {
    .abort("highlights/meta.activePageId '", active, "' is not in the page inventory")
  }
  undeclared <- setdiff(memberships, paste0(.highlight_cells_prefix, group_ids))
  if (length(undeclared) > 0L) {
    .abort("Session contains undeclared highlight membership chunks: ", paste(sort(undeclared, method = "radix"), collapse = ", "))
  }
  if (!identical(memberships, paste0(.highlight_cells_prefix, group_ids))) {
    .abort("Highlight membership chunks must match the exact metadata group order")
  }
  list(plans = plans, uns = list(groups = groups_uns), meta = meta, existing = existing)
}
