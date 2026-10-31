# cellucid_apply_session(): a saved session's highlight groups and
# user-defined fields, written back onto the object it was saved on, as
# metadata columns (Seurat meta.data, SingleCellExperiment and
# SummarizedExperiment colData). cellucid-python's
# apply_cellucid_session_to_anndata() with the same rules: every check runs
# before the object is changed, so a session that fails leaves nothing
# half-written. The columns are written by the writer cellucid_pull() uses
# (R/round-trip.R). The help page is man/cellucid_apply_session.Rd.

cellucid_apply_session <- function(
    data,
    session,
    ...,
    expected_dataset_id = NULL,
    add_highlights = TRUE,
    highlights_prefix = "cellucid_highlight__",
    add_user_defined_fields = TRUE,
    user_defined_prefix = "",
    include_deleted_user_defined_fields = FALSE,
    store_metadata = TRUE,
    return_summary = FALSE,
    timeout = 60
) {
  if (missing(data) || missing(session)) {
    .abort("cellucid_apply_session() needs data (the object) and session (a .cellucid-session file or a viewer).")
  }
  flags <- list(
    add_highlights = add_highlights, add_user_defined_fields = add_user_defined_fields,
    include_deleted_user_defined_fields = include_deleted_user_defined_fields,
    store_metadata = store_metadata, return_summary = return_summary
  )
  for (name in names(flags)) .validate_flag(flags[[name]], name)
  for (name in c("highlights_prefix", "user_defined_prefix")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L || is.na(value)) .abort(name, " must be one string.")
  }
  extra <- list(...)
  if (length(extra) > 0L && (is.null(names(extra)) || any(!nzchar(names(extra))))) {
    .abort("Arguments after session must be named: the ones cellucid_view() was given.")
  }
  unknown <- setdiff(names(extra), .view_plan_arguments)
  if (length(unknown) > 0L) {
    .abort("cellucid_apply_session() takes the arguments cellucid_view() was given (",
           .format_value_list(.view_plan_arguments), "); not ", .format_value_list(unknown), ".")
  }

  if (inherits(session, "cellucid_viewer")) {
    .abort("Ask the viewer for its session first: request <- cellucid_request_session(viewer), then ",
           "cellucid_apply_session(object, request) (in a notebook, in a later cell).")
  }
  if (inherits(session, "cellucid_session_request")) {
    # Asked with cellucid_request_session(); in a notebook, applied in a
    # later cell.
    if (is.null(expected_dataset_id)) expected_dataset_id <- session$viewer$dataset_id
    path <- cellucid_session_result(session, timeout = timeout)
    on.exit(unlink(path), add = TRUE)
    session <- path
  }
  if (is.null(expected_dataset_id)) {
    .abort("expected_dataset_id is required: the dataset id the session must have been saved on ",
           "(the viewer's dataset_id).")
  }
  expected_dataset_id <- .validate_required_string(expected_dataset_id, "expected_dataset_id")
  planned <- .plan_session(data, session, extra, expected_dataset_id, add_highlights, highlights_prefix,
                           add_user_defined_fields, user_defined_prefix, include_deleted_user_defined_fields,
                           .object_column_names(data))
  plans <- planned$plans
  record <- if (store_metadata) {
    session_record <- list(
      # Without `source`: where the viewer opened the dataset (a link, a server
      # address, a local file name) is not part of the object, and would travel
      # with every file written from it.
      manifest = planned$bundle$manifest[setdiff(names(planned$bundle$manifest), "source")],
      dataset_fingerprint = planned$fingerprint,
      applied = list(expected_dataset_id = expected_dataset_id, contract = "exact"),
      materialized_fields = stats::setNames(lapply(plans, `[[`, "metadata"), vapply(plans, `[[`, "", "name")),
      chunks = c(
        if (!is.null(planned$highlights$meta)) list(`highlights/meta` = planned$highlights$meta),
        if (!is.null(planned$fields$overlays)) list(`core/field-overlays` = planned$fields$overlays)
      )
    )
    if (!is.null(planned$highlights$uns)) session_record$highlights <- planned$highlights$uns
    session_record
  }
  data <- .object_with_columns(data, plans)
  if (!is.null(record)) {
    held <- .object_record(data)
    if (!is.null(held$session) && !is.list(held$session)) .abort("The object's cellucid$session must be a list if present.")
    held$session <- record
    data <- .object_with_record(data, held)
  }
  added <- vapply(plans, `[[`, "", "name")
  if (return_summary) list(data = data, added_columns = added) else data
}

# Refuse an object that lacks a field the session uses, naming every one.
.require_session_fields <- function(view, requires, mapping) {
  missing_obs <- setdiff(requires$obs, names(view$obs))
  missing_var <- if (length(requires$var) > 0L) setdiff(requires$var, .session_gene_names(view, mapping)) else character()
  if (length(missing_obs) == 0L && length(missing_var) == 0L) return(invisible(NULL))
  quote <- function(keys) paste0("'", keys, "'", collapse = ", ")
  parts <- c(
    if (length(missing_obs) > 0L) paste0("the obs column", if (length(missing_obs) > 1L) "s", " ", quote(missing_obs)),
    if (length(missing_var) > 0L) paste0("the gene", if (length(missing_var) > 1L) "s", " ", quote(missing_var))
  )
  .abort("This session uses ", paste(parts, collapse = " and "), ", which this dataset does not have.")
}

# A session checked against the object and its columns planned, nothing
# written: the session file (a path), the arguments cellucid_view() was given
# (`extra`), and `existing`, the column names a planned column must not take.
# Returns the plans, the bundle, its fingerprint and the highlight and field
# plans (their metadata chunks). cellucid_pull() settles the names itself.
.plan_session <- function(data, session, extra, expected_dataset_id, add_highlights, highlights_prefix,
                          add_user_defined_fields, user_defined_prefix, include_deleted_user_defined_fields,
                          existing) {
  bundle <- .read_session_bundle(session)
  # An object read from a file in priority order (R/priority.R) holds its
  # cells in the stored order; a session names cells in their original order.
  # The record says so only when it describes these cells (.stored_order()):
  # one that does not leaves each row its own original position, as the web
  # app read the file the session was saved on.
  stored <- .stored_order_of(data)
  if (!is.null(stored$notice)) .warn(stored$notice)
  if (!is.null(stored$positions)) {
    .abort(
      "This object holds its cells in a Cellucid-ready file's stored order (uns['cellucid']['priority']), ",
      "and a session names cells in their original order. ",
      "Put them back first: cellucid_restore_original_order(object)."
    )
  }
  # Every per-cell chunk names its cells by original position, ascending (the
  # viewer's rule, whatever order the file or the opening held them in): every
  # cell, cell i is column i; a subsample, cell i is its kept position
  # kept[i]. A subsample's session is checked and planned against those
  # columns of the object, then each column is placed at them (the others NA).
  cells <- if ("cells" %in% names(bundle$manifest)) {
    .session_cells(bundle$manifest$cells, bundle$manifest$datasetFingerprint$cellCount)
  }
  saved_on <- data
  if (!is.null(cells)) {
    n <- length(.object_cells(data))
    if (n != cells$spec$nObs) {
      .abort("Dataset fingerprint mismatch: the session was saved on ", cells$spec$k, " of ", cells$spec$nObs, " ",
             .plural(cells$spec$nObs, "cell"), ", and the object has ", n, ".")
    }
    saved_on <- data[, cells$kept + 1]
  }
  view <- .session_view(saved_on, extra)

  fingerprint <- bundle$manifest$datasetFingerprint
  mismatches <- character()
  if (!identical(fingerprint$datasetId, expected_dataset_id)) {
    mismatches <- c(mismatches, paste0("datasetId '", fingerprint$datasetId %||% "null", "' != expected_dataset_id '", expected_dataset_id, "'"))
  }
  if (fingerprint$cellCount != view$n_cells) {
    mismatches <- c(mismatches, paste0("cellCount ", fingerprint$cellCount, " != the object's ", view$n_cells, " ", .plural(view$n_cells, "cell")))
  }
  if (length(mismatches) > 0L) .abort("Dataset fingerprint mismatch: ", paste(mismatches, collapse = "; "))
  mapping <- fingerprint$dataMapping
  tryCatch(
    .require_mapping_modality(mapping, view$modality_names),
    cellucid_data_mapping_error = function(error) .abort(conditionMessage(error))
  )
  .require_mapping_matches_view(view, mapping)
  genes <- .session_var_count(view, mapping)
  if (fingerprint$varCount != genes) {
    .abort("Dataset fingerprint mismatch: varCount ", fingerprint$varCount, " != ",
           if (is.null(mapping$expression)) "genes on a mapping without expression" else "the mapped genes", " ",
           genes)
  }
  requires <- lapply(bundle$manifest$requires, function(keys) as.character(unlist(keys)))
  .require_session_fields(view, requires, mapping)
  chunks <- .validate_session_chunks(bundle)

  highlights <- .plan_highlights(bundle, chunks, view, highlights_prefix, existing, add_highlights)
  fields <- .plan_user_fields(
    bundle, chunks, view, mapping, user_defined_prefix, highlights$existing,
    add_user_defined_fields, include_deleted_user_defined_fields
  )
  plans <- c(highlights$plans, fields$plans)
  if (!is.null(cells)) {
    n <- length(.object_cells(data))
    plans <- lapply(plans, function(plan) {
      plan$values <- .place_rows(plan$values, cells$kept, n)
      plan
    })
  }
  list(plans = plans, bundle = bundle, fingerprint = fingerprint, highlights = highlights, fields = fields)
}
