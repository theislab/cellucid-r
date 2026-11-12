# What a dataset says about itself, stored in uns['cellucid']: its name,
# its portable identifier, its description, its citation and its source.

# The white space prose loses at either end: ASCII white space, the Unicode
# space separators, the line and paragraph separators and the byte order mark
# (the web app's PROSE_EDGE, data/anndata-cellucid-metadata.js;
# cellucid-python's _contracts.PROSE_EDGE).
.prose_edge_pattern <- paste0(
  "(*UTF)",
  "^[\\t\\n\\x{0b}\\f\\r\\x{20}\\x{a0}\\x{1680}\\x{2000}-\\x{200a}\\x{2028}\\x{2029}\\x{202f}\\x{205f}\\x{3000}\\x{feff}]+",
  "|[\\t\\n\\x{0b}\\f\\r\\x{20}\\x{a0}\\x{1680}\\x{2000}-\\x{200a}\\x{2028}\\x{2029}\\x{202f}\\x{205f}\\x{3000}\\x{feff}]+$"
)

# Descriptive text -- a dataset's name, description, citation text, source
# name -- trimmed (.prose_edge_pattern), then drawable as stored
# (.validate_display_text): a space typed at either end of a paragraph is
# removed, never refused, as the web app reads it. Identifiers keep the strict
# rule.
.validate_prose_text <- function(value, arg, allow_empty = FALSE) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    value <- gsub(.prose_edge_pattern, "", enc2utf8(value), perl = TRUE)
  }
  .validate_display_text(value, arg, allow_empty = allow_empty)
}

.validate_dataset_name <- function(value) {
  .validate_prose_text(value, "dataset_name")
}

.validate_dataset_id <- function(value) {
  .safe_filename_component(value, what = "dataset_id")
}

.validate_optional_text <- function(value, arg, allow_empty = FALSE) {
  if (is.null(value)) {
    return(NULL)
  }
  .validate_display_text(value, arg, allow_empty = allow_empty)
}

# The viewer links to it, so it must be an absolute web address.
.validate_url <- function(value, arg) {
  value <- .validate_optional_text(value, arg)
  if (!is.null(value) && !grepl("^https?://[^/?#[:space:]]+", value, ignore.case = TRUE)) {
    .abort(
      arg, " must be an absolute http:// or https:// address; got '", value, "'."
    )
  }
  value
}

# A link is stored only beside the text it belongs to: a source URL needs
# the source's name, a citation URL needs the citation's text.
.validate_linked_text <- function(text, url, text_arg, url_arg, text_key) {
  text <- if (!is.null(text)) .validate_prose_text(text, text_arg)
  url <- .validate_url(url, url_arg)
  if (is.null(text)) {
    if (!is.null(url)) {
      .abort(text_arg, " is required whenever ", url_arg, " is supplied.")
    }
    return(NULL)
  }
  entry <- list()
  entry[[text_key]] <- .uns_string(text)
  if (!is.null(url)) {
    entry$url <- .uns_string(url)
  }
  entry
}
