.create_test_directory_reparse <- function(source, destination) {
  if (.Platform$OS.type == "windows") {
    return(isTRUE(Sys.junction(source, destination)))
  }
  isTRUE(file.symlink(source, destination))
}

.remove_test_reparse <- function(path, directory = FALSE) {
  stopifnot(
    is.logical(directory),
    length(directory) == 1L,
    !is.na(directory)
  )
  if (.Platform$OS.type == "windows" && !directory) {
    # R 4.6.1's recursive unlink path hard-codes the directory-junction
    # reparse tag. file.symlink() uses the symbolic-link tag instead, while
    # file.remove() delegates file-link removal to the Windows CRT.
    return(if (isTRUE(file.remove(path))) 0L else 1L)
  }
  unlink(path, recursive = TRUE, force = TRUE)
}
