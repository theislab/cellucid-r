# Where the store is written, and how it reaches that path: the whole store is
# written beside its target under a temporary name, and moved into place only
# once it is complete, so an interrupted or failed call never leaves a
# half-written store at the path a viewer would open.

.dir_create <- function(path) {
  if (!dir.exists(path)) {
    created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!created && !dir.exists(path)) {
      .abort("Could not create directory: ", path)
    }
  }
  invisible(path)
}

# The suffix tells the web app, a host, and a person what the path holds, so
# out_path must end in the suffix of the chosen format: a folder ending in
# ".zarr", or one ZIP file ending in ".zarr.zip".
.resolve_output_path <- function(out_path, format, force) {
  out_path <- .validate_required_string(out_path, "out_path")
  zip <- identical(format, "zarr.zip")
  suffix <- .output_formats[[format]]
  name <- basename(out_path)
  lowered <- tolower(name)
  present <- if (endsWith(lowered, ".zarr.zip")) {
    ".zarr.zip"
  } else if (endsWith(lowered, ".zarr")) {
    ".zarr"
  } else if (endsWith(lowered, ".h5ad")) {
    ".h5ad"
  } else {
    NA_character_
  }
  if (!identical(present, suffix) || nchar(name) <= nchar(suffix)) {
    stem <- if (is.na(present)) name else substr(name, 1L, nchar(name) - nchar(present))
    if (!nzchar(stem)) stem <- "data"
    instead <- if (grepl("[/\\\\]", out_path)) {
      file.path(dirname(out_path), paste0(stem, suffix))
    } else {
      paste0(stem, suffix)
    }
    problem <- if (is.na(present)) {
      paste0("does not end in \"", suffix, "\"")
    } else if (identical(present, suffix)) {
      "is only the suffix"
    } else if (identical(present, ".h5ad")) {
      paste0(
        "ends in \".h5ad\", which this package does not write (Python's ",
        "cellucid.prepare(format = \"h5ad\") does)"
      )
    } else {
      paste0(
        "ends in \"", present, "\", which format = \"",
        names(.output_formats)[.output_formats == present], "\" writes"
      )
    }
    .abort(
      "format = \"", format, "\" writes ",
      if (zip) "one .zarr.zip file" else "a .zarr folder",
      ", and out_path '", out_path, "' ", problem, ": write '", instead,
      "' instead."
    )
  }
  parent <- normalizePath(dirname(path.expand(out_path)), winslash = "/", mustWork = FALSE)
  path <- file.path(parent, basename(out_path))
  if (file.exists(path)) {
    if (!force) {
      .abort(path, " already exists; pass force = TRUE to replace it.")
    }
    is_folder <- dir.exists(path)
    if (zip && is_folder) {
      .abort(
        path, " is a folder, and format = \"zarr.zip\" writes a file; remove it ",
        "yourself if it should go."
      )
    }
    is_store <- is_folder && (file.exists(file.path(path, ".zgroup")) || file.exists(file.path(path, "zarr.json")))
    if (!zip && !is_store) {
      .abort(
        path, " exists and is not a Zarr store, so force = TRUE does not ",
        "replace it; remove it yourself if it should go."
      )
    }
  }
  list(path = path, parent = parent, zip = zip)
}

# What a writer killed before it finished (its process gone) left beside the
# target, removed, as cellucid-python's writer removes it: its staging
# (.<name>.<pid>.partial), and the target it had moved aside
# (.<name>.<pid>.replaced) when the target is there again -- when it is not,
# the moved-aside copy is the only one, and it is left.
.remove_stale_staging <- function(target) {
  name <- basename(target$path)
  pattern <- paste0("^\\.", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", name), "\\.([0-9]+)\\.(partial|replaced)$")
  for (entry in list.files(target$parent, all.files = TRUE, no.. = TRUE)) {
    found <- regmatches(entry, regexec(pattern, entry))[[1L]]
    if (length(found) != 3L) next
    pid <- as.numeric(found[[2L]])
    if (pid == Sys.getpid() || .native_process_alive(pid)) next
    if (found[[3L]] == "replaced" && !file.exists(target$path)) next
    unlink(file.path(target$parent, entry), recursive = TRUE, force = TRUE)
  }
  invisible(NULL)
}

# Write the store to a temporary sibling of the target, then move it into
# place, replacing what was there only after the new store is complete.
.write_output <- function(target, write) {
  .dir_create(target$parent)
  .remove_stale_staging(target)
  tag <- paste0(".", basename(target$path), ".", Sys.getpid())
  staging <- file.path(target$parent, paste0(tag, ".partial"))
  replaced <- file.path(target$parent, paste0(tag, ".replaced"))
  unlink(c(staging, replaced), recursive = TRUE, force = TRUE)
  store <- .store_open(staging, zip = target$zip, zarr_format = target$zarr_format)
  on.exit({
    .store_abort(store)
    unlink(staging, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  write(store)
  .store_close(store)

  had_target <- file.exists(target$path)
  if (had_target && !file.rename(target$path, replaced)) {
    .abort("Could not move the existing ", target$path, " aside to replace it.")
  }
  if (!file.rename(staging, target$path)) {
    if (had_target) file.rename(replaced, target$path)
    .abort("Could not move the finished store to ", target$path, ".")
  }
  if (had_target) {
    unlink(replaced, recursive = TRUE, force = TRUE)
  }
  invisible(target$path)
}
