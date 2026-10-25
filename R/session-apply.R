# cellucid_apply_session(): a saved session's highlight groups and
# user-defined fields, written back onto the object it was saved on, as
# metadata columns (Seurat meta.data, SingleCellExperiment and
# SummarizedExperiment colData). cellucid-python's
# apply_cellucid_session_to_anndata() with the same rules: every check runs
# before the object is changed, so a session that fails leaves nothing
# half-written. The help page is man/cellucid_apply_session.Rd.

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
    if (is.null(expected_dataset_id)) expected_dataset_id <- session$dataset_id
    path <- cellucid_session(session, timeout = timeout)
    on.exit(unlink(path), add = TRUE)
    session <- path
  }
  if (is.null(expected_dataset_id)) {
    .abort("expected_dataset_id is required: the dataset id the session must have been saved on ",
           "(the viewer's dataset_id).")
  }
  expected_dataset_id <- .validate_required_string(expected_dataset_id, "expected_dataset_id")
  bundle <- .read_session_bundle(session)
  view <- .session_view(data, extra)

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

  existing <- .object_column_names(data)
  highlights <- .plan_highlights(bundle, chunks, view, highlights_prefix, existing, add_highlights)
  fields <- .plan_user_fields(
    bundle, chunks, view, mapping, user_defined_prefix, highlights$existing,
    add_user_defined_fields, include_deleted_user_defined_fields
  )
  plans <- c(highlights$plans, fields$plans)
  record <- if (store_metadata) {
    session_record <- list(
      manifest = bundle$manifest,
      dataset_fingerprint = fingerprint,
      applied = list(expected_dataset_id = expected_dataset_id, contract = "exact"),
      materialized_fields = stats::setNames(lapply(plans, `[[`, "metadata"), vapply(plans, `[[`, "", "name")),
      chunks = c(
        if (!is.null(highlights$meta)) list(`highlights/meta` = highlights$meta),
        if (!is.null(fields$overlays)) list(`core/field-overlays` = fields$overlays)
      )
    )
    if (!is.null(highlights$uns)) session_record$highlights <- highlights$uns
    session_record
  }
  data <- .object_with_columns(data, plans, record)
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

.object_column_names <- function(data) {
  if (methods::is(data, "Seurat")) names(methods::slot(data, "meta.data")) else names(SummarizedExperiment::colData(data))
}

# The object with each planned column added, and the session recorded under
# misc$cellucid$session (Seurat) or metadata(x)$cellucid$session.
.object_with_columns <- function(data, plans, record) {
  if (methods::is(data, "Seurat")) {
    meta <- methods::slot(data, "meta.data")
    for (plan in plans) meta[[plan$name]] <- plan$values
    methods::slot(data, "meta.data") <- meta
    if (!is.null(record)) {
      misc <- methods::slot(data, "misc")
      methods::slot(data, "misc") <- .with_session_record(misc, record, "misc")
    }
    return(data)
  }
  col_data <- SummarizedExperiment::colData(data)
  for (plan in plans) col_data[[plan$name]] <- plan$values
  SummarizedExperiment::colData(data) <- col_data
  if (!is.null(record)) {
    S4Vectors::metadata(data) <- .with_session_record(S4Vectors::metadata(data), record, "metadata(x)")
  }
  data
}

.with_session_record <- function(holder, record, where) {
  cellucid <- holder$cellucid
  if (is.null(cellucid)) cellucid <- list()
  if (!is.list(cellucid)) .abort(where, "$cellucid must be a list if present.")
  if (!is.null(cellucid$session) && !is.list(cellucid$session)) {
    .abort(where, "$cellucid$session must be a list if present.")
  }
  cellucid$session <- record
  holder$cellucid <- cellucid
  holder
}
