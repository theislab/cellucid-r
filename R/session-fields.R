# The obs columns that carry a session's user-defined fields: a categorical
# field's codes and categories (src/session_codecs.c), or a continuous alias
# of an obs column or a gene, copied as the viewer copied it (NA stays NA).
# Each is checked against the field overlay the session declares first.
# cellucid-python's anndata_session/_fields.py, rule for rule.

.user_field_common_keys <- c(
  "id", "source", "kind", "key", "isDeleted", "isPurged", "sourceField", "operation", "createdAt"
)
.user_field_category_keys <- c(
  .user_field_common_keys, "categories", "codesLength", "codesType", "sourcePages", "overlapStrategy",
  "overlapLabel", "intersectionLabels", "uncoveredLabel"
)
.codes_maximum <- c(Uint8Array = 255, Uint16Array = 65535, Uint32Array = 4294967295)

# A field name as the viewer accepts one: no surrounding white space, no
# control or zero-width character; a ':' is allowed (ATAC peaks).
.session_field_key <- function(value, label) {
  key <- .session_text(value, label)
  if (!.is_trimmed_text(key)) .abort(label, " cannot have leading or trailing whitespace")
  points <- .code_points(key)
  for (rule in list(
    list(points <= 0x1F | (points >= 0x7F & points <= 0x9F), "control"),
    list(points %in% c(0x200B, 0x2060, 0xFEFF), "zero-width")
  )) {
    if (any(rule[[1L]])) {
      .abort(label, " cannot contain the ", rule[[2L]], " character U+", sprintf("%04X", points[rule[[1L]]][[1L]]))
    }
  }
  key
}

.session_flag <- function(value, label) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) .abort(label, " must be exactly True or False")
  value
}

.validate_source_field <- function(value, label) {
  if (!.is_json_object(value)) .abort(label, " must be an object")
  if (!"sourceKey" %in% names(value) || length(setdiff(names(value), c("sourceKey", "sourceIndex", "kind"))) > 0L) {
    .abort(label, " has invalid fields")
  }
  out <- list(sourceKey = .session_text(value[["sourceKey"]], paste0(label, ".sourceKey")))
  if ("sourceIndex" %in% names(value)) {
    index <- value[["sourceIndex"]]
    if (!.is_json_integer(index)) .abort(label, ".sourceIndex must be an integer")
    if (index < 0) .abort(label, ".sourceIndex must be non-negative")
    out$sourceIndex <- index
  }
  if ("kind" %in% names(value)) out$kind <- .session_text(value[["kind"]], paste0(label, ".kind"))
  out
}

.validate_operation <- function(value, label) {
  if (is.null(value)) return(invisible(NULL))
  if (!.is_json_object(value)) .abort(label, " must be an object or null")
  .session_text(value[["type"]], paste0(label, ".type"))
  invisible(NULL)
}

# A category label as the viewer shows it: JavaScript's String() of the value.
.category_label <- function(value) {
  if (is.logical(value)) return(if (value) "true" else "false")
  if (is.numeric(value)) return(if (value == floor(value) && abs(value) < 1e21) format(value, scientific = FALSE, trim = TRUE) else as.character(value))
  value
}

.validate_category_field <- function(field, label) {
  categories <- field[["categories"]]
  if (!.is_json_array(categories)) .abort(label, ".categories must be an array")
  if (length(categories) == 0L) .abort(label, ".categories must be a non-empty array")
  identities <- character()
  for (index in seq_along(categories)) {
    value <- .session_category_value(categories[[index]], paste0(label, ".categories[", index - 1L, "]"))
    identity <- paste0(if (is.logical(value)) "b:" else if (is.numeric(value)) "n:" else "s:", as.character(value))
    if (identity %in% identities) .abort(label, ".categories must contain unique exact values")
    identities <- c(identities, identity)
  }
  length <- field[["codesLength"]]
  if (!.is_json_integer(length)) .abort(label, ".codesLength must be an integer")
  if (length < 0) .abort(label, ".codesLength must be non-negative")
  type <- field[["codesType"]]
  if (!.is_json_string(type) || !type %in% names(.codes_maximum)) {
    .abort(label, ".codesType must be Uint8Array, Uint16Array or Uint32Array")
  }
  if (!.is_json_array(field[["sourcePages"]])) .abort(label, ".sourcePages must be an array")
  if (!.is_json_string(field[["overlapStrategy"]]) ||
      !field[["overlapStrategy"]] %in% c("first", "last", "overlap-label", "intersections")) {
    .abort(label, ".overlapStrategy is invalid")
  }
  for (key in c("overlapLabel", "uncoveredLabel")) {
    value <- field[[key]]
    if (!is.null(value) && !.is_json_string(value)) .abort(label, ".", key, " must be a string or null")
    if (identical(value, "")) .abort(label, ".", key, " must be non-empty or null")
  }
  if (!is.null(field[["intersectionLabels"]]) && !.is_json_object(field[["intersectionLabels"]])) {
    .abort(label, ".intersectionLabels must be an object or null")
  }
  if (!is.null(field[["sourceField"]])) .validate_source_field(field[["sourceField"]], paste0(label, ".sourceField"))
  .validate_operation(field[["operation"]], paste0(label, ".operation"))
  if (length > 0 && length(categories) > .codes_maximum[[type]]) {
    .abort(label, " has too many categories for ", type, " codes")
  }
  vapply(categories, .category_label, character(1L))
}

.validate_string_record <- function(value, label) {
  if (!.is_json_object(value)) .abort(label, " must be an object")
  for (key in names(value)) {
    if (!nzchar(key)) .abort(label, " key must be a non-empty string")
    .session_text(value[[key]], paste0(label, " value for '", key, "'"))
  }
}

.unique_strings <- function(value, label) {
  if (!.is_json_array(value)) .abort(label, " must be an array")
  out <- character()
  for (index in seq_along(value)) {
    item <- .session_text(value[[index]], paste0(label, "[", index - 1L, "]"))
    if (item %in% out) .abort(label, " duplicates '", item, "'")
    out <- c(out, item)
  }
  out
}

# One continuous alias's values: an obs column, or a gene named the way the
# session's mapping named it, from the matrix the view reads.
.continuous_alias_values <- function(view, mapping, source, source_field, label) {
  key <- source_field$sourceKey
  if (identical(source, "obs")) {
    index <- which(names(view$obs) == key)
    if (length(index) != 1L) .abort(label, " source obs field '", key, "' must exist exactly once")
    values <- view$obs[[index]]
  } else {
    if (is.null(mapping$expression)) {
      .abort(label, " cannot copy gene values: The data mapping reads no expression matrix.")
    }
    index <- which(.session_gene_names(view, mapping) == key)
    if (length(index) != 1L) .abort(label, " source gene '", key, "' must exist exactly once")
    values <- .session_gene_values(view, mapping, index)
  }
  declared <- source_field$sourceIndex
  if (!is.null(declared) && declared != index - 1L) {
    .abort(label, " sourceIndex ", declared, " does not match exact source '", key, "' at index ", index - 1L)
  }
  if (is.factor(values) || !(is.numeric(values) && !is.logical(values))) {
    .abort(label, " must contain real numeric values")
  }
  as.numeric(values)
}

# Returns list(plans, overlays, existing): one planned column per field the
# session applies, and the decoded core/field-overlays chunk.
.plan_user_fields <- function(bundle, chunks, view, mapping, prefix, existing, materialize, include_deleted) {
  ids <- names(chunks)
  code_ids <- ids[startsWith(ids, .user_codes_prefix)]
  if (!"core/field-overlays" %in% ids) {
    if (length(code_ids) > 0L) .abort("Session contains user-defined code chunks without core/field-overlays")
    return(list(plans = list(), overlays = NULL, existing = existing))
  }
  overlays <- .session_check_keys(.session_chunk(bundle, "core/field-overlays"), c("renames", "deletedFields", "userDefinedFields"), "core/field-overlays")
  renames <- .session_check_keys(overlays[["renames"]], c("fields", "categories"), "core/field-overlays.renames")
  .validate_string_record(renames[["fields"]], "core/field-overlays.renames.fields")
  .validate_string_record(renames[["categories"]], "core/field-overlays.renames.categories")
  deleted_fields <- .session_check_keys(overlays[["deletedFields"]], c("deleted", "purged"), "core/field-overlays.deletedFields")
  deleted <- .unique_strings(deleted_fields[["deleted"]], "core/field-overlays.deletedFields.deleted")
  purged <- .unique_strings(deleted_fields[["purged"]], "core/field-overlays.deletedFields.purged")
  if (length(setdiff(purged, deleted)) > 0L) .abort("core/field-overlays.deletedFields.purged must be a subset of deleted")
  fields <- overlays[["userDefinedFields"]]
  if (!.is_json_array(fields)) .abort("core/field-overlays.userDefinedFields must be an array")

  plans <- list()
  field_ids <- character()
  referenced <- character()
  categorical_ids <- character()
  for (index in seq_along(fields)) {
    label <- paste0("core/field-overlays.userDefinedFields[", index - 1L, "]")
    field <- fields[[index]]
    if (!.is_json_object(field)) .abort(label, " must be an object")
    kind <- field[["kind"]]
    .session_check_keys(field, if (identical(kind, "category")) .user_field_category_keys else .user_field_common_keys, label)
    if (!.is_json_string(kind) || !kind %in% c("category", "continuous")) .abort(label, ".kind must be 'category' or 'continuous'")
    field_id <- .session_text(field[["id"]], paste0(label, ".id"))
    if (!grepl(.wire_id_pattern, field_id)) {
      .abort(label, ".id must use only ASCII letters, digits, '.', '_', or '-' and be at most 180 characters")
    }
    if (field_id %in% field_ids) .abort("Duplicate user-defined field id: ", field_id)
    field_ids <- c(field_ids, field_id)
    source <- field[["source"]]
    if (!.is_json_string(source) || !source %in% c("obs", "var")) .abort(label, ".source must be 'obs' or 'var'")
    key <- .session_field_key(field[["key"]], paste0(label, ".key"))
    is_deleted <- .session_flag(field[["isDeleted"]], paste0(label, ".isDeleted"))
    is_purged <- .session_flag(field[["isPurged"]], paste0(label, ".isPurged"))
    if (is_purged && !is_deleted) .abort(label, " purged fields must also be deleted")
    created <- field[["createdAt"]]
    if (!is.null(created) && (!is.numeric(created) || length(created) != 1L || !is.finite(created) || created < 0)) {
      .abort(label, ".createdAt must be a finite non-negative number or null")
    }
    apply_it <- materialize && (include_deleted || !is_deleted) && !is_purged
    column <- paste0(prefix, key)
    if (identical(kind, "continuous")) {
      source_field <- .validate_source_field(field[["sourceField"]], paste0(label, ".sourceField"))
      .validate_operation(field[["operation"]], paste0(label, ".operation"))
      if (apply_it) {
        if (column %in% existing) .abort("Column already exists: ", column)
        existing <- c(existing, column)
        values <- .continuous_alias_values(view, mapping, source, source_field, paste0("Continuous field '", field_id, "'"))
        plans[[length(plans) + 1L]] <- list(
          name = column, values = values,
          metadata = list(kind = "continuous", field_id = field_id, source = source, source_field = source_field)
        )
      }
      next
    }
    labels <- .validate_category_field(field, label)
    categorical_ids <- c(categorical_ids, field_id)
    codes_type <- field[["codesType"]]
    chunk_id <- paste0(.user_codes_prefix, field_id)
    if (field[["codesLength"]] != view$n_cells) {
      .abort("User-defined categorical field '", field_id, "' has codesLength ", field[["codesLength"]],
             ", expected ", view$n_cells, " cell-aligned codes")
    }
    if (view$n_cells == 0L) {
      if (chunk_id %in% ids) .abort("Empty user-defined field '", field_id, "' must not contain a codes chunk")
      codes <- numeric()
    } else {
      if (!chunk_id %in% ids) {
        .abort("User-defined categorical field '", field_id, "' is missing required chunk '", chunk_id, "'")
      }
      if (!identical(chunks[[chunk_id]][["label"]], paste0("User-defined codes: ", key))) {
        .abort(chunk_id, " label must match user-defined field '", field_id, "'")
      }
      codes <- .native_decode_user_codes(.session_chunk(bundle, chunk_id), field[["codesLength"]], codes_type)
      referenced <- c(referenced, chunk_id)
    }
    if (is.null(field[["uncoveredLabel"]])) codes[codes == .codes_maximum[[codes_type]]] <- -1
    invalid <- which(codes < -1 | codes >= length(labels))
    if (length(invalid) > 0L) {
      .abort("User-defined categorical field '", field_id, "' contains code ", format(codes[[invalid[[1L]]]], scientific = FALSE),
             " at position ", invalid[[1L]] - 1L, ", but declares ", length(labels), " ",
             .plural(length(labels), "category", "categories"))
    }
    if (apply_it) {
      if (column %in% existing) .abort("Column already exists: ", column)
      if (anyDuplicated(labels)) {
        .abort("User-defined categorical field '", field_id, "' has categories that read as the same label: ",
               .format_value_list(unique(labels[duplicated(labels)])))
      }
      existing <- c(existing, column)
      codes[codes == -1] <- NA
      plans[[length(plans) + 1L]] <- list(
        name = column, values = factor(labels[codes + 1], levels = labels),
        metadata = list(kind = "category", field_id = field_id, source = source)
      )
    }
  }
  undeclared <- setdiff(code_ids, referenced)
  if (length(undeclared) > 0L) {
    .abort("Session contains undeclared user-defined code chunks: ", paste(sort(undeclared, method = "radix"), collapse = ", "))
  }
  for (priority in c("eager", "lazy")) {
    actual <- substring(code_ids[vapply(chunks[code_ids], function(meta) identical(meta[["priority"]], priority), logical(1L))],
                        nchar(.user_codes_prefix) + 1L)
    if (!identical(actual, categorical_ids[categorical_ids %in% actual])) {
      .abort("User-defined ", priority, " code chunks must match the exact field-overlay inventory order")
    }
  }
  list(plans = plans, overlays = overlays, existing = existing)
}
