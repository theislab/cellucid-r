# One export generation at a time per output directory. The native advisory
# lock is the cross-process half; the in-process registry below is the half
# that catches a second generation started from the same session, and that
# discards handles a fork inherited.

.export_generation_lock_registry <- new.env(parent = emptyenv())
.export_generation_lock_registry$pid <- Sys.getpid()
.export_generation_lock_registry$paths <- new.env(parent = emptyenv())

.native_export_lock_cleanup_status <- function(lock_handle) {
  status <- .native_export_lock_release(lock_handle)
  if (
    !is.integer(status) ||
      length(status) != 3L ||
      anyNA(status)
  ) {
    stop("Native export lock cleanup returned an invalid status.", call. = FALSE)
  }
  status
}

.assert_native_export_lock_cleanup <- function(status, lock_path) {
  if (status[[1L]] != 0L || status[[2L]] != 0L) {
    stop(
      "Could not completely release the export generation lock ",
      lock_path,
      " (unlock error ",
      status[[1L]],
      ", close error ",
      status[[2L]],
      ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.refresh_export_generation_lock_registry <- function() {
  process_id <- Sys.getpid()
  if (!identical(.export_generation_lock_registry$pid, process_id)) {
    inherited_registry <- .export_generation_lock_registry$paths
    inherited_keys <- ls(envir = inherited_registry, all.names = TRUE)
    .export_generation_lock_registry$pid <- process_id
    .export_generation_lock_registry$paths <- new.env(parent = emptyenv())
    cleanup_failures <- character()
    for (lock_key in inherited_keys) {
      record <- get(lock_key, envir = inherited_registry, inherits = FALSE)
      if (!is.null(record$handle)) {
        status <- tryCatch(
          .native_export_lock_cleanup_status(record$handle),
          error = identity
        )
        if (inherits(status, "error")) {
          cleanup_failures <- c(
            cleanup_failures,
            paste0(record$path, ": ", conditionMessage(status))
          )
        } else if (status[[1L]] != 0L || status[[2L]] != 0L) {
          cleanup_failures <- c(
            cleanup_failures,
            paste0(
              record$path,
              ": unlock error ",
              status[[1L]],
              ", close error ",
              status[[2L]]
            )
          )
        }
      }
    }
    if (length(cleanup_failures) != 0L) {
      stop(
        "Could not discard every fork-inherited export lock handle: ",
        paste(cleanup_failures, collapse = "; "),
        call. = FALSE
      )
    }
  }
  process_id
}

.export_generation_lock_path <- function(final_path) {
  parent <- normalizePath(
    dirname(final_path),
    winslash = "/",
    mustWork = TRUE
  )
  file.path(
    parent,
    paste0(".", basename(final_path), ".cellucid.lock")
  )
}

.validate_export_generation_lock_path <- function(lock_path) {
  path_info <- .export_path_info(lock_path)
  if (path_info$kind == 3L) {
    stop(
      "Export lock path must not be a symbolic link or reparse point: ",
      lock_path,
      call. = FALSE
    )
  }
  if (.Platform$OS.type != "windows") {
    link_target <- suppressWarnings(Sys.readlink(lock_path))
    if (
      length(link_target) == 1L &&
        !is.na(link_target) &&
        nzchar(link_target)
    ) {
      stop(
        "Export lock path must not be a symbolic link: ",
        lock_path,
        call. = FALSE
      )
    }
  }
  if (
    path_info$kind != 0L &&
      (path_info$kind != 1L || !utils::file_test("-f", lock_path))
  ) {
    stop(
      "Export lock path must identify a regular file: ",
      lock_path,
      call. = FALSE
    )
  }
  invisible(lock_path)
}

.canonical_export_generation_lock_key <- function(lock_path) {
  key <- normalizePath(lock_path, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows") {
    key <- tolower(key)
  }
  key
}

.acquire_export_generation_lock <- function(final_path) {
  lock_path <- .export_generation_lock_path(final_path)
  .validate_export_generation_lock_path(lock_path)
  lock_key <- .canonical_export_generation_lock_key(lock_path)
  process_id <- .refresh_export_generation_lock_registry()
  registry <- .export_generation_lock_registry$paths

  if (exists(lock_key, envir = registry, inherits = FALSE)) {
    stop(
      "An export generation is already active for ",
      final_path,
      ".",
      call. = FALSE
    )
  }

  assign(
    lock_key,
    list(handle = NULL, path = lock_path, pid = process_id),
    envir = registry
  )
  reservation_is_owned <- TRUE
  on.exit({
    if (
      reservation_is_owned &&
        exists(lock_key, envir = registry, inherits = FALSE)
    ) {
      rm(list = lock_key, envir = registry)
    }
  }, add = TRUE)

  lock <- .native_export_lock_acquire(lock_path)
  if (is.null(lock)) {
    stop(
      "An export generation is already active for ",
      final_path,
      ".",
      call. = FALSE
    )
  }

  export_lock <- list(
    handle = lock,
    key = lock_key,
    path = lock_path,
    pid = process_id
  )
  assign(lock_key, export_lock, envir = registry)
  reservation_is_owned <- FALSE
  export_lock
}

.release_export_generation_lock <- function(export_lock) {
  process_id <- .refresh_export_generation_lock_registry()
  if (!identical(export_lock$pid, process_id)) {
    return(invisible(FALSE))
  }

  registry <- .export_generation_lock_registry$paths
  if (!exists(export_lock$key, envir = registry, inherits = FALSE)) {
    stop(
      "Export generation lock ownership is not active: ",
      export_lock$path,
      call. = FALSE
    )
  }
  active_lock <- get(export_lock$key, envir = registry, inherits = FALSE)
  if (!identical(active_lock$handle, export_lock$handle)) {
    stop(
      "Export generation lock ownership does not match: ",
      export_lock$path,
      call. = FALSE
    )
  }

  status <- .native_export_lock_cleanup_status(export_lock$handle)
  if (status[[3L]] == 1L) {
    rm(list = export_lock$key, envir = registry)
  }
  .assert_native_export_lock_cleanup(status, export_lock$path)
  invisible(TRUE)
}
