# Where the export is allowed to be written, and creating the directories
# that hold it.

.validate_output_path <- function(path) {
  path <- .validate_required_string(path, "out_dir")
  expanded <- path.expand(path)
  resolved <- normalizePath(
    expanded,
    winslash = "/",
    mustWork = FALSE
  )
  protected <- unique(c(
    "/",
    normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    normalizePath(path.expand("~"), winslash = "/", mustWork = TRUE)
  ))
  if (resolved %in% protected) {
    stop(
      "out_dir must name a dedicated dataset output directory.",
      call. = FALSE
    )
  }
  leaf <- basename(expanded)
  if (!nzchar(leaf) || leaf %in% c(".", "..")) {
    stop(
      "out_dir must name a dedicated dataset output directory.",
      call. = FALSE
    )
  }
  parent <- normalizePath(
    dirname(expanded),
    winslash = "/",
    mustWork = FALSE
  )
  final_path <- file.path(parent, leaf)
  if (final_path %in% protected) {
    stop(
      "out_dir must name a dedicated dataset output directory.",
      call. = FALSE
    )
  }
  final_path
}

.dir_create <- function(path) {
  if (!dir.exists(path)) {
    created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!created && !dir.exists(path)) {
      stop("Could not create directory: ", path, call. = FALSE)
    }
  }
  invisible(path)
}
