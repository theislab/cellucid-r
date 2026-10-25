# The chunks a current session holds, checked before any is applied: every
# generic chunk carries exactly the metadata of its profile, every chunk id
# is one this version knows (highlight cell sets and user-defined codes by
# their exact id patterns, so an id cannot smuggle a path), and the
# contributors appear in their registered order. cellucid-python's
# anndata_session/_schema.py and _chunks.py, rule for rule.

.session_static_profiles <- list(
  list(id = "core/field-overlays", contributorId = "field-overlays", priority = "eager", kind = "json", codec = "gzip", label = "Field overlays", datasetDependent = TRUE),
  list(id = "core/state", contributorId = "core-state", priority = "eager", kind = "json", codec = "gzip", label = "Core state", datasetDependent = TRUE),
  list(id = "ui/dockable-layout", contributorId = "dockable-layout", priority = "eager", kind = "json", codec = "gzip", label = "Floating panels", datasetDependent = FALSE),
  list(id = "ui/figure-export", contributorId = "figure-export", priority = "eager", kind = "json", codec = "gzip", label = "Figure export", datasetDependent = FALSE),
  list(id = "analysis/windows", contributorId = "analysis-windows", priority = "eager", kind = "binary", codec = "gzip", label = "Analysis windows", datasetDependent = TRUE),
  list(id = "analysis/docked", contributorId = "analysis-docked", priority = "eager", kind = "binary", codec = "gzip", label = "Docked analysis", datasetDependent = TRUE),
  list(id = "highlights/meta", contributorId = "highlights-meta", priority = "eager", kind = "json", codec = "gzip", label = "Highlight metadata", datasetDependent = TRUE),
  list(id = "cinematic/camera", contributorId = "cinematic-camera", priority = "eager", kind = "json", codec = "gzip", label = "Cinematic camera path", datasetDependent = TRUE)
)
.session_static_ids <- vapply(.session_static_profiles, function(profile) profile$id, character(1L))
.session_contributor_order <- c(
  "field-overlays", "user-defined-codes", "cinematic-camera", "core-state", "dockable-layout",
  "figure-export", "analysis-windows", "analysis-docked", "highlights-meta", "highlights-cells"
)
.highlight_cells_prefix <- "highlights/cells/"
.user_codes_prefix <- "user-defined/codes/"
.wire_id_pattern <- "^[A-Za-z0-9][A-Za-z0-9._-]{0,179}$"

.expected_chunk_contributor <- function(id) {
  if (id %in% .session_static_ids) {
    return(.session_static_profiles[[match(id, .session_static_ids)]]$contributorId)
  }
  if (startsWith(id, .highlight_cells_prefix)) {
    group <- substring(id, nchar(.highlight_cells_prefix) + 1L)
    return(if (grepl("^highlight_[1-9][0-9]*$", group)) "highlights-cells" else NULL)
  }
  if (startsWith(id, .user_codes_prefix)) {
    field <- substring(id, nchar(.user_codes_prefix) + 1L)
    return(if (grepl(.wire_id_pattern, field)) "user-defined-codes" else NULL)
  }
  NULL
}

# Check the session's chunk listing; returns each chunk's metadata by id.
.validate_session_chunks <- function(bundle) {
  last <- c(eager = 0L, lazy = 0L)
  saw_lazy <- FALSE
  for (meta in bundle$manifest$chunks) {
    id <- meta[["id"]]
    contributor <- .expected_chunk_contributor(id)
    if (is.null(contributor)) .abort("Unknown current session chunk '", id, "'")
    if (!identical(meta[["contributorId"]], contributor)) {
      .abort("Session chunk '", id, "' requires contributor '", contributor, "', received '",
             meta[["contributorId"]], "'")
    }
    priority <- meta[["priority"]]
    if (identical(priority, "lazy")) {
      saw_lazy <- TRUE
    } else if (saw_lazy) {
      .abort("Session eager chunks must precede lazy chunks")
    }
    rank <- match(contributor, .session_contributor_order)
    if (rank < last[[priority]]) {
      .abort("Session ", priority, " contributor groups must match their registered order")
    }
    last[[priority]] <- rank
    if (id %in% .session_static_ids) {
      profile <- .session_static_profiles[[match(id, .session_static_ids)]]
      for (key in names(profile)) {
        if (!identical(meta[[key]], profile[[key]])) {
          .abort("Session chunk '", id, "' must match its exact current profile")
        }
      }
    } else if (identical(contributor, "user-defined-codes")) {
      if (!identical(meta[["kind"]], "binary") || !identical(meta[["codec"]], "gzip") || !isTRUE(meta[["datasetDependent"]])) {
        .abort("Session chunk '", id, "' must be binary, gzip, and dataset-dependent")
      }
    } else if (!identical(priority, "lazy") || !identical(meta[["kind"]], "binary") ||
               !identical(meta[["codec"]], "gzip") || !isTRUE(meta[["datasetDependent"]])) {
      .abort("Session chunk '", id, "' must be lazy, binary, gzip, and dataset-dependent")
    }
  }
  missing <- setdiff(.session_static_ids, names(bundle$chunks))
  if (length(missing) > 0L) {
    .abort("Current session requires singleton chunks: ", paste(missing, collapse = ", "))
  }
  bundle$chunks
}
