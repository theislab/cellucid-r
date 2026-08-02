# The identity a dataset publishes about itself: its name, its portable
# identifier, its description, when it was created, and where it came from.

.validate_dataset_name <- function(value) {
  .validate_display_text(value, "dataset_name")
}

.validate_dataset_id <- function(value) {
  .safe_filename_component(value, what = "dataset_id")
}

# dataset_description is the one identity string whose empty value carries
# meaning: a dataset with nothing to say about itself and a dataset described
# as "" are the same dataset, and both are published as "". NULL is what an R
# caller writes for an argument they are not supplying -- every other optional
# identity argument in this signature already defaults to it -- so NULL and ""
# are both accepted and both normalize to the "" that reaches the manifest.
# The cellucid Python package accepts None and "" for it in the same way.
.validate_dataset_description <- function(value) {
  if (is.null(value)) {
    return("")
  }
  .validate_display_text(value, "dataset_description", allow_empty = TRUE)
}

.resolve_created_at <- function(value) {
  if (is.null(value)) {
    return(format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC",
      usetz = FALSE
    ))
  }
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value)
  ) {
    stop(
      "created_at must be NULL or exactly one string in ",
      "'YYYY-MM-DDTHH:MM:SSZ' UTC format.",
      call. = FALSE
    )
  }
  if (!grepl(
    "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$",
    value,
    perl = TRUE,
    useBytes = TRUE
  )) {
    stop(
      "created_at must use exact 'YYYY-MM-DDTHH:MM:SSZ' UTC format.",
      call. = FALSE
    )
  }

  year <- as.integer(substr(value, 1L, 4L))
  month <- as.integer(substr(value, 6L, 7L))
  day <- as.integer(substr(value, 9L, 10L))
  hour <- as.integer(substr(value, 12L, 13L))
  minute <- as.integer(substr(value, 15L, 16L))
  second <- as.integer(substr(value, 18L, 19L))
  leap_year <- year %% 400L == 0L ||
    (year %% 4L == 0L && year %% 100L != 0L)
  days_in_month <- c(
    31L,
    if (leap_year) 29L else 28L,
    31L,
    30L,
    31L,
    30L,
    31L,
    31L,
    30L,
    31L,
    30L,
    31L
  )
  if (
    year < 1L ||
      month < 1L ||
      month > 12L ||
      day < 1L ||
      day > days_in_month[[month]] ||
      hour < 0L ||
      hour > 23L ||
      minute < 0L ||
      minute > 59L ||
      second < 0L ||
      second > 59L
  ) {
    stop(
      "created_at must be a valid UTC calendar timestamp in ",
      "'YYYY-MM-DDTHH:MM:SSZ' format.",
      call. = FALSE
    )
  }
  value
}

.validate_source_identity <- function(name, url, citation) {
  if (!is.null(name)) {
    name <- .validate_display_text(name, "source_name")
  }
  if (!is.null(url)) {
    url <- .validate_display_text(url, "source_url")
  }
  if (!is.null(citation)) {
    citation <- .validate_display_text(
      citation,
      "source_citation"
    )
  }
  if (is.null(name) && is.null(url) && is.null(citation)) {
    return(NULL)
  }
  if (is.null(name)) {
    stop(
      "source_name is required whenever source_url or source_citation ",
      "is supplied.",
      call. = FALSE
    )
  }
  source <- list(name = name)
  if (!is.null(url)) {
    source$url <- url
  }
  if (!is.null(citation)) {
    source$citation <- citation
  }
  source
}
