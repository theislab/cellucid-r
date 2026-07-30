.onUnload <- function(libpath) {
  cleanup_status <- .Call(C_cellucid_export_lock_drain)
  if (
    !is.integer(cleanup_status) ||
      !is.matrix(cleanup_status) ||
      ncol(cleanup_status) != 3L ||
      anyNA(cleanup_status) ||
      any(!(cleanup_status[, 3L] %in% c(0L, 1L)))
  ) {
    stop(
      "Native export lock drain returned an invalid status matrix.",
      call. = FALSE
    )
  }

  cleanup_failures <- which(
    cleanup_status[, 1L] != 0L |
      cleanup_status[, 2L] != 0L |
      cleanup_status[, 3L] != 1L
  )
  if (length(cleanup_failures) > 0L) {
    details <- vapply(
      cleanup_failures,
      function(index) {
        paste0(
          "#",
          index,
          " unlock/close/closed=",
          paste(cleanup_status[index, ], collapse = "/")
        )
      },
      character(1)
    )
    stop(
      "Cannot unload Cellucid after native export lock cleanup failures: ",
      paste(details, collapse = "; "),
      call. = FALSE
    )
  }

  .export_generation_lock_registry$pid <- Sys.getpid()
  .export_generation_lock_registry$paths <- new.env(parent = emptyenv())
  library.dynam.unload("cellucid", libpath)
}
