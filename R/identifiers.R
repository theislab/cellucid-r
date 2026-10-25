# dataset_id is the one identifier with a filename rule: catalogs and hosts
# address a dataset by it, so it must be a portable file name component on
# every platform.
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
    .abort(
      what, " '", name, "' is not a portable identifier. Use ",
      "1-180 ASCII letters, numbers, '.', '_', or '-', beginning with a ",
      "letter or number and not ending with '.'."
    )
  }
  device_base <- toupper(strsplit(name, ".", fixed = TRUE)[[1]][[1]])
  if (
    device_base %in% c("CON", "PRN", "AUX", "NUL") ||
      grepl("^(COM|LPT)[1-9]$", device_base)
  ) {
    .abort(
      what, " '", name, "' is reserved on Windows."
    )
  }
  name
}
