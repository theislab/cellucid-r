# What a name may be when it becomes an element of the store: an obs or var
# column, an obsm key. Each one is a path component in the Zarr store -- a
# folder on disk, or an entry name in a .zarr.zip that someone may unpack --
# so it obeys the rules a portable path component obeys, on every platform,
# whichever platform writes it: a store has to open where it is copied to.

# Characters Windows refuses in a file name, plus the path separators.
.element_key_forbidden <- "[/\\\\<>:\"|?*\u0001-\u001f]"
.zarr_metadata_names <- c(".zarray", ".zattrs", ".zgroup", ".zmetadata", "zarr.json")

.validate_element_key <- function(key, what) {
  if (!is.character(key) || length(key) != 1L || is.na(key) || !nzchar(key)) {
    stop(what, " names must be non-empty strings.", call. = FALSE)
  }
  problem <- NULL
  if (grepl(.element_key_forbidden, key, perl = TRUE)) {
    problem <- paste0(
      "contains a character that cannot appear in a file name on every ",
      "platform (/ \\ < > : \" | ? * or a control character)"
    )
  } else if (key %in% c(".", "..") || key %in% .zarr_metadata_names) {
    problem <- "is reserved by the Zarr store layout"
  } else if (grepl("[ .]$", key)) {
    problem <- "ends with a space or a dot, which Windows removes from file names"
  } else {
    base <- toupper(strsplit(key, ".", fixed = TRUE)[[1L]][[1L]])
    if (base %in% c("CON", "PRN", "AUX", "NUL") || grepl("^(COM|LPT)[1-9]$", base)) {
      problem <- "is a reserved device name on Windows"
    }
  }
  if (!is.null(problem)) {
    stop(
      what, " '", key, "' ", problem, ". It becomes a path in the Zarr ",
      "store, so rename it, for example with names(x)[names(x) == '", key,
      "'] <- 'new_name'.",
      call. = FALSE
    )
  }
  invisible(key)
}

# Names that live side by side in one group must differ by more than letter
# case: on macOS and Windows "Sex" and "sex" are the same folder, and one
# would overwrite the other.
.require_distinct_element_keys <- function(keys, what) {
  folded <- tolower(keys)
  repeated <- which(duplicated(folded))
  if (length(repeated) > 0L) {
    later <- repeated[[1L]]
    first <- match(folded[[later]], folded)
    if (identical(keys[[first]], keys[[later]])) {
      stop(what, " name '", keys[[later]], "' is used twice.", call. = FALSE)
    }
    stop(
      what, " names '", keys[[first]], "' and '", keys[[later]], "' differ ",
      "only in letter case. They become folders of one Zarr group, and on ",
      "macOS and Windows those are the same folder, so rename one of them.",
      call. = FALSE
    )
  }
  invisible(keys)
}

.validate_frame_columns <- function(frame, what) {
  columns <- names(frame)
  for (column in columns) {
    .validate_element_key(column, what)
  }
  if ("_index" %in% columns) {
    stop(
      what, " '_index' is reserved: AnnData stores the row names under ",
      "that name. Rename the column.",
      call. = FALSE
    )
  }
  .require_distinct_element_keys(columns, what)
  invisible(columns)
}

# --- The vector-field key rule ---------------------------------------------

# scVelo's naming convention, the rule the web app applies when a store opens
# (parseVectorFieldKey) and cellucid-python's parse_vector_field_key: an obsm
# key '<field>_<basis>' is the vector field '<field>' drawn over the embedding
# 'X_<basis>' ('velocity_umap' over 'X_umap', 'velocity_branches' over
# 'X_branches'). A key that starts with 'X_' names an embedding, never a
# vector field. Otherwise every split at one '_' into a non-empty field and a
# non-empty basis is tried from the left, and the first whose 'X_<basis>' is
# among `obsm_keys` wins: the longest basis. Matching is exact and
# case-sensitive. Returns list(field, basis), or NULL for an ordinary array.
# The shared data-mapping conformance fixture pins it (vectorFieldKeys).
.embedding_key_prefix <- "X_"

.parse_vector_field_key <- function(key, obsm_keys) {
  if (startsWith(key, .embedding_key_prefix)) {
    return(NULL)
  }
  characters <- strsplit(key, "", fixed = TRUE)[[1L]]
  n <- length(characters)
  for (split in which(characters == "_")) {
    if (split == 1L || split == n) {
      next
    }
    basis <- paste(characters[(split + 1L):n], collapse = "")
    if (paste0(.embedding_key_prefix, basis) %in% obsm_keys) {
      return(list(field = paste(characters[seq_len(split - 1L)], collapse = ""), basis = basis))
    }
  }
  NULL
}
