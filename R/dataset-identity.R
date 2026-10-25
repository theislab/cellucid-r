# What a dataset says about itself, stored in uns['cellucid']: its name,
# its portable identifier, its description, its citation and its source.

.validate_dataset_name <- function(value) {
  .validate_display_text(value, "dataset_name")
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
  text <- .validate_optional_text(text, text_arg)
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
