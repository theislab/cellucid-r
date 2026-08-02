# Guarded filesystem operations shared by the generation lock and the export
# transaction. Every one of them states what the path must already be before
# it acts, so a symlink, a hard link, or an unexpected file type stops the
# export instead of being followed.

.export_path_info <- function(path) {
  info <- .native_export_path_info(path)
  if (
    !is.double(info) ||
      length(info) != 2L ||
      anyNA(info) ||
      info[[1L]] < 0 ||
      info[[1L]] > 4 ||
      info[[1L]] != floor(info[[1L]]) ||
      info[[2L]] < 0 ||
      info[[2L]] != floor(info[[2L]])
  ) {
    stop("Native export path inspection returned invalid state.", call. = FALSE)
  }
  list(
    kind = as.integer(info[[1L]]),
    links = info[[2L]]
  )
}

.require_export_directory_or_absent <- function(path, label) {
  info <- .export_path_info(path)
  if (info$kind == 0L) {
    return(FALSE)
  }
  if (info$kind != 2L) {
    stop(
      label,
      " must be an ordinary non-symbolic directory or absent: ",
      path,
      call. = FALSE
    )
  }
  TRUE
}

.require_export_regular_file <- function(path, label) {
  info <- .export_path_info(path)
  if (info$kind == 0L) {
    stop(label, " is missing: ", path, call. = FALSE)
  }
  if (info$kind != 1L || info$links != 1) {
    stop(
      label,
      " must be one non-linked, non-symbolic regular file: ",
      path,
      call. = FALSE
    )
  }
  invisible(info)
}

.fsync_export_directory <- function(path) {
  status <- .native_export_sync_directory(path)
  if (
    !is.logical(status) ||
      length(status) != 1L ||
      is.na(status)
  ) {
    stop(
      "Native export directory synchronization returned invalid state.",
      call. = FALSE
    )
  }
  invisible(status)
}

.rename_export_path <- function(source, destination) {
  if (!isTRUE(file.rename(source, destination))) {
    stop(
      "Could not rename export transaction path ",
      source,
      " to ",
      destination,
      ".",
      call. = FALSE
    )
  }
  .fsync_export_directory(dirname(destination))
  invisible(destination)
}

.remove_export_tree <- function(path) {
  if (!.require_export_directory_or_absent(
    path,
    "Export transaction directory"
  )) {
    return(invisible(path))
  }
  status <- unlink(path, recursive = TRUE, force = TRUE)
  if (status != 0L || .export_path_info(path)$kind != 0L) {
    stop(
      "Could not remove export transaction directory: ",
      path,
      call. = FALSE
    )
  }
  .fsync_export_directory(dirname(path))
  invisible(path)
}

.remove_export_control_file <- function(path) {
  .require_export_regular_file(
    path,
    "Export transaction control file"
  )
  status <- unlink(path, recursive = FALSE, force = TRUE)
  if (status != 0L || .export_path_info(path)$kind != 0L) {
    stop(
      "Could not remove export transaction control file: ",
      path,
      call. = FALSE
    )
  }
  .fsync_export_directory(dirname(path))
  invisible(path)
}
