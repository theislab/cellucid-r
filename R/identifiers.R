# What a name may be. Two rules apply, and which one applies is decided by
# whether the name reaches a path: dataset_id names the export directory and
# obeys the portable-filename rule, while a field identifier names one entry
# in a manifest and obeys the display-text and uniqueness rules instead.

# dataset_id is the one identifier that still names a path: the export
# directory a producer publishes and every data path addresses it by. Payload
# names on the obs, var, and vectors axes are integer indices, so no field
# identifier reaches this rule.
.safe_filename_component <- function(name, what = "Filename component") {
  # The portable-identifier rule below already rejects an empty string and any
  # whitespace, wherever it sits, so it is the rule that must speak: telling a
  # caller who passed '' that the value carries edge whitespace names a defect
  # the value does not have. Only the type check runs first, because a value
  # that is not one string cannot be quoted back into a message about it.
  name <- .validate_string(name, what)
  if (
    !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", name) ||
      grepl("\\.$", name) ||
      nchar(name, type = "bytes") > 180L
  ) {
    stop(
      what, " '", name, "' is not a portable identifier. Use ",
      "1-180 ASCII letters, numbers, '.', '_', or '-', beginning with a ",
      "letter or number and not ending with '.'.",
      call. = FALSE
    )
  }
  device_base <- toupper(strsplit(name, ".", fixed = TRUE)[[1]][[1]])
  if (
    device_base %in% c("CON", "PRN", "AUX", "NUL") ||
      grepl("^(COM|LPT)[1-9]$", device_base)
  ) {
    stop(
      what, " '", name, "' is reserved on Windows.",
      call. = FALSE
    )
  }
  name
}

# A payload filename is an integer index, so an identifier is never a path and
# carries no filename rule at all: no ASCII restriction, no case-collision
# rule, no reserved Windows name. What survives is what the identity is
# actually for. It names one field in the manifest, so it must be a non-empty
# string, distinct within its axis so a lookup resolves one field, and text the
# viewer can draw exactly as it is stored -- the same rule every category label
# obeys, because a gene name and a category label are shown in the same legend.
# The cellucid Python package enforces the identical rule.
#
# `what` names one field in the singular -- "Gene", "Observation field",
# "Vector field" -- because it is read as the subject of a sentence about a
# single identifier: "Gene identifier at position 1 ...". The checks that speak
# about the whole axis add their own plural, so no caller has to guess which
# number a message will need. cellucid-python passes the same singular nouns to
# _require_field_identities().
.require_field_identities <- function(values, what) {
  .validate_character_vector(values, what = paste0(what, " keys"))
  for (position in seq_along(values)) {
    .validate_display_text(
      values[[position]],
      paste0(what, " identifier at position ", position - 1L),
      allow_empty = FALSE
    )
  }
  .require_unique_identifiers(values, what = what)
}

# Identity, not path. A repeated key names no single field, so a lookup by key
# resolves to more than one row whatever the payloads are called. The reported
# key is the first repeat in supplied order, which is the one
# cellucid-python's _require_unique_identifiers() reports.
.require_unique_identifiers <- function(keys, what) {
  .validate_character_vector(keys, what = paste0(what, " keys"))
  repeated <- which(duplicated(keys))
  if (length(repeated) > 0L) {
    stop(
      what,
      " key '",
      keys[[repeated[[1L]]]],
      "' is duplicated.",
      call. = FALSE
    )
  }
  invisible(keys)
}
