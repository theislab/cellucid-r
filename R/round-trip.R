# The round trip: what a person made in the Cellucid viewer, written onto the
# object in one call (cellucid_pull()), and columns and gene sets pushed into
# the viewer (cellucid_push(), cellucid_push_gene_set()). cellucid-python's
# anndata_session/_results.py and AnnDataViewer.pull()/push(), rule for rule:
# the same column names, dtypes and missing-value rules, the same results.
#
# The viewer sends every result it holds (viewer-protocol.json, `result`) once
# its dataset opens and whenever its channel opens again, then each change;
# everything here reads what it sent and never waits (a Jupyter notebook hears
# the viewer only between cells):
# highlight groups, user-defined categorical fields (labels, annotations,
# clusters), continuous fields and scores, gene sets, the visible cells and
# the cells it opened. Cells are 1-based positions here, as everywhere in
# this package; the wire's are 0-based.
#
# Column names, one scheme:
# - a highlight group: <highlight_prefix><group id>, TRUE/FALSE;
# - a categorical or continuous field, a score: <field_prefix><key>;
# - the visible cells (visibility = TRUE): cellucid_visible.
# When the viewer opened a subsample, a cell it did not open is NA. Columns an
# earlier pull wrote (recorded under cellucid$pulled) are replaced, and those
# whose result is gone removed (prune); any other column of the same name is
# a conflict (on_conflict). Help pages: man/cellucid_pull.Rd and
# man/cellucid_push.Rd.

.pull_on_conflict <- c("error", "replace", "suffix")
.visible_column <- "cellucid_visible"
.obs_name_column <- "obs_name"

# --- The columns, one rule per kind (shared with cellucid_apply_session()) ---

# TRUE for `cells` (0-based, NULL for every cell); NA outside `open_cells`.
.membership_column <- function(cells, open_cells, n) {
  values <- if (is.null(cells)) rep(TRUE, n) else logical(n)
  if (!is.null(cells) && length(cells) > 0L) values[cells + 1] <- TRUE
  if (!is.null(open_cells)) {
    known <- logical(n)
    known[open_cells + 1] <- TRUE
    if (is.null(cells)) values <- known
    values[!known] <- NA
  }
  values
}

# A factor from codes over `cells` (0-based, NULL for every cell);
# `missing_code` names no category; a cell outside `cells` has none.
.categorical_column <- function(codes, categories, cells, n, missing_code) {
  codes <- as.double(codes)
  codes[codes == missing_code] <- NA
  if (any(!is.na(codes) & (codes < 0 | codes >= length(categories)))) .abort("A code names no category.")
  placed <- if (is.null(cells)) codes else {
    out <- rep(NA_real_, n)
    out[cells + 1] <- codes
    out
  }
  factor(categories[placed + 1], levels = categories)
}

# Numbers over `cells` (0-based, NULL for every cell), NA elsewhere.
.continuous_column <- function(values, cells, n) {
  values <- as.double(values)
  values[is.nan(values)] <- NA
  if (is.null(cells)) return(values)
  out <- rep(NA_real_, n)
  out[cells + 1] <- values
  out
}

# A column of some cells (one value per 0-based position in `rows`) placed in
# a column of all `n`: the others NA.
.place_rows <- function(values, rows, n) {
  if (is.factor(values)) {
    out <- factor(rep(NA_character_, n), levels = levels(values))
    out[rows + 1] <- values
    return(out)
  }
  out <- if (is.logical(values)) rep(NA, n) else rep(NA_real_, n)
  out[rows + 1] <- values
  out
}

# --- The object: its columns and its record -----------------------------------

.object_kind <- function(data) {
  if (methods::is(data, "Seurat")) return("seurat")
  if (methods::is(data, "SummarizedExperiment")) return("se")
  if (is.data.frame(data)) return("data.frame")
  .abort("Pull into a Seurat, SingleCellExperiment or SummarizedExperiment object, or a data frame of its cells.")
}

.object_cells <- function(data) {
  switch(.object_kind(data),
    seurat = rownames(methods::slot(data, "meta.data")),
    se = colnames(data),
    data.frame = rownames(data)
  )
}

.object_column_names <- function(data) {
  switch(.object_kind(data),
    seurat = names(methods::slot(data, "meta.data")),
    se = names(SummarizedExperiment::colData(data)),
    data.frame = names(data)
  )
}

# What cellucid keeps on the object: misc$cellucid (Seurat),
# metadata(x)$cellucid, or the "cellucid" attribute of a data frame.
.object_record <- function(data) {
  record <- switch(.object_kind(data),
    seurat = methods::slot(data, "misc")$cellucid,
    se = S4Vectors::metadata(data)$cellucid,
    data.frame = attr(data, "cellucid", exact = TRUE)
  )
  if (is.null(record)) list() else if (is.list(record)) record else .abort("The object's cellucid record must be a list.")
}

.object_with_record <- function(data, record) {
  switch(.object_kind(data),
    seurat = {
      misc <- methods::slot(data, "misc")
      misc$cellucid <- record
      methods::slot(data, "misc") <- misc
      data
    },
    se = {
      S4Vectors::metadata(data)$cellucid <- record
      data
    },
    data.frame = {
      attr(data, "cellucid") <- record
      data
    }
  )
}

# The one writer: each planned column added (or replacing the column of its
# name), the columns in `remove` dropped, in one assignment. R copies on
# change: the object passed in is not changed.
.object_with_columns <- function(data, plans, remove = character()) {
  write <- function(table) {
    for (name in remove) table[[name]] <- NULL
    for (plan in plans) table[[plan$name]] <- plan$values
    table
  }
  switch(.object_kind(data),
    seurat = {
      methods::slot(data, "meta.data") <- write(methods::slot(data, "meta.data"))
      data
    },
    se = {
      SummarizedExperiment::colData(data) <- write(SummarizedExperiment::colData(data))
      data
    },
    data.frame = write(data)
  )
}

# --- The viewer's results -----------------------------------------------------

# Every result the viewer holds, as of the last event it sent: once its
# dataset opens (and whenever its channel opens again) the viewer sends every
# result and the inventory, then each change. Nothing is asked: in a Jupyter
# notebook the viewer is heard only between cells, so a lasso drawn and
# cellucid_pull() run in the next cell finds the lasso. The results are one
# epoch's (one page's account of one dataset it opened), and must describe
# the dataset this viewer serves.
cellucid_results <- function(viewer) {
  .require_viewer(viewer)
  inventory <- viewer$results_inventory
  if (is.null(inventory)) {
    .abort("The viewer has not sent what it holds yet: it does once its dataset has opened. ",
           "Run this again when the viewer shows the data.")
  }
  if (is.null(viewer$results_dataset)) {
    .abort("The viewer shows a dataset this session did not open: another one was opened in it. ",
           "Open '", viewer$dataset_id, "' in the viewer again (or show the viewer again) to pull from it.")
  }
  if (!identical(viewer$results_dataset, viewer$dataset_id)) {
    .abort("The viewer's results describe the dataset '", viewer$results_dataset, "', and this viewer serves '",
           viewer$dataset_id, "': open '", viewer$dataset_id, "' in the viewer again to pull from it.")
  }
  missing <- setdiff(inventory, names(viewer$results))
  if (length(missing) > 0L) {
    .abort("The viewer's results arrived incomplete (", paste(utils::head(missing, 3L), collapse = ", "),
           " missing); reload the viewer's page to have it send them again.")
  }
  unname(viewer$results)
}

# One delivered result as R holds it: cells 1-based (NULL for every cell),
# read from the viewer (`source`; a file's results carry the file's name).
.result_record <- function(payload, cell_names) {
  cells <- if (is.null(payload$cells)) NULL else payload$cells + 1
  record <- list(
    kind = payload$kind, id = payload$id, revision = payload$revision, n_obs = payload$nObs,
    meta = payload$meta, cells = cells, values = payload$values, value_type = payload$value_type,
    source = .viewer_source
  )
  if (!is.null(cell_names) && !is.null(cells)) record$cell_names <- cell_names[cells]
  record
}

.result_key <- function(result) paste(result$kind, result$id, sep = ":")

# --- Results to columns -------------------------------------------------------

# Where a pull read a result: the live viewer (a file's results carry the
# file's name). Every column a pull writes is recorded under cellucid$pulled
# with pulled_from, kind and id: the result it holds.
.viewer_source <- "viewer"

# The result a pulled record names ("<pulled_from>\r<kind>\r<id>"), or NA for
# a record that names none (never replaced or removed by a pull).
.pulled_key <- function(info) {
  parts <- list(info$pulled_from, info$kind, info$id)
  if (!all(vapply(parts, function(part) is.character(part) && length(part) == 1L && !is.na(part), TRUE))) {
    return(NA_character_)
  }
  paste(unlist(parts), collapse = "\r")
}

.pulled_from_of_key <- function(key) strsplit(key, "\r", fixed = TRUE)[[1L]][[1L]]

# A column name, unique among the pull's own. A result an earlier pull wrote
# takes the column it wrote then -- its own name, or that name with the suffix
# it took (unless on_conflict = "replace" asks for the name itself) -- so a
# pull run again writes the same column. Otherwise its name is used when no
# column has it, or when a pull of the same source wrote it; any other
# existing column is a conflict, settled by on_conflict.
.pull_name <- function(state, wanted, key) {
  earlier <- state$ours[key][[1L]]
  suffixed <- !is.null(earlier) && grepl(paste0("^", .regex_literal(wanted), "_[0-9]+$"), earlier)
  name <- if (!is.null(earlier) && !earlier %in% state$used &&
              (identical(earlier, wanted) || (!identical(state$on_conflict, "replace") && suffixed))) {
    earlier
  } else if (wanted %in% state$used) {
    .suffixed_name(wanted, c(state$existing, state$used))
  } else if (!wanted %in% state$existing || identical(state$written_by[wanted][[1L]], .pulled_from_of_key(key)) ||
             identical(state$on_conflict, "replace")) {
    wanted
  } else if (identical(state$on_conflict, "suffix")) {
    .suffixed_name(wanted, c(state$existing, state$used))
  } else {
    state$conflicts <- c(state$conflicts, wanted)
    wanted
  }
  state$used <- c(state$used, name)
  name
}

.regex_literal <- function(text) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", text)

.suffixed_name <- function(wanted, taken) {
  number <- 2L
  while (paste0(wanted, "_", number) %in% taken) number <- number + 1L
  paste0(wanted, "_", number)
}

# The columns of `results` before naming: each a key, the name it wants, its
# values and its record; and the results skipped (with why).
.result_candidates <- function(data, results, options) {
  n <- length(.object_cells(data))
  for (result in results) {
    if (result$n_obs != n) {
      .abort("The viewer shows a dataset of ", result$n_obs, " ", .plural(result$n_obs, "cell"), " and this object has ", n,
             ": pull into the object the viewer was opened on, with its cells in the same order.")
    }
  }
  open <- NULL
  for (result in results) if (identical(result$kind, "subset") && !is.null(result$cells)) open <- result$cells - 1
  candidates <- list()
  skipped <- list()
  ordered <- results[order(vapply(results, `[[`, "", "kind"), vapply(results, `[[`, 0, "revision"))]
  for (result in ordered) {
    meta <- result$meta
    cells <- if (is.null(result$cells)) NULL else result$cells - 1
    column <- result$column %||% NULL
    source <- result$source %||% .viewer_source
    provenance <- list(kind = result$kind, id = result$id, pulled_from = source)
    if (identical(result$kind, "highlight")) {
      if (!options$highlights) next
      wanted <- column %||% paste0(options$highlight_prefix, result$id)
      values <- .membership_column(cells, if (is.null(result$known)) open else result$known, n)
      info <- c(provenance, list(label = meta$label, tool = meta$tool, page = meta$pageName, enabled = meta$enabled))
    } else if (result$kind %in% c("categorical", "continuous")) {
      origin <- meta$origin
      if (identical(origin, "notebook") && !options$include_pushed) {
        skipped[[.result_key(result)]] <- "pushed from a notebook (include_pushed = TRUE brings it back)"
        next
      }
      if (identical(origin, "analysis") && !options$scores) next
      if (!identical(origin, "analysis") && !options$fields) next
      wanted <- column %||% paste0(options$field_prefix, meta$key)
      values <- if (identical(result$kind, "categorical")) {
        missing <- 2^(8 * .viewer_value_types[[result$value_type]]) - 1
        .categorical_column(result$values, as.character(unlist(meta$categories)), cells, n, missing)
      } else {
        .continuous_column(result$values, cells, n)
      }
      info <- c(provenance, list(key = meta$key, origin = origin))
    } else if (identical(result$kind, "visibility")) {
      if (!options$visibility) next
      wanted <- .visible_column
      values <- .membership_column(cells, open, n)
      info <- provenance
    } else {
      next
    }
    candidates[[length(candidates) + 1L]] <- list(key = .pulled_key(info), wanted = wanted, values = values, metadata = info)
  }
  list(candidates = candidates, skipped = skipped)
}

# The plans of `candidates`, named (.pull_name), and the pulled record they
# are named against.
.name_candidates <- function(data, candidates, on_conflict) {
  previous <- .object_record(data)$pulled %||% list()
  columns <- .object_column_names(data)
  state <- new.env(parent = emptyenv())
  state$existing <- columns
  state$ours <- list()
  state$written_by <- list()
  for (column in intersect(names(previous), columns)) {
    key <- .pulled_key(previous[[column]])
    if (is.na(key)) next
    state$ours[[key]] <- column
    state$written_by[[column]] <- .pulled_from_of_key(key)
  }
  state$used <- character()
  state$conflicts <- character()
  state$on_conflict <- on_conflict
  plans <- lapply(candidates, function(candidate) {
    list(name = .pull_name(state, candidate$wanted, candidate$key), values = candidate$values,
         metadata = candidate$metadata)
  })
  if (length(state$conflicts) > 0L) {
    .abort("These columns exist already and no pull from the same source wrote them: ",
           .format_value_list(state$conflicts),
           ". Pass on_conflict = \"replace\" to overwrite them, on_conflict = \"suffix\" to write beside them, ",
           "or another highlight_prefix/field_prefix.")
  }
  list(plans = plans, previous = previous)
}

# The columns of an earlier viewer pull to remove: those whose result the
# viewer no longer holds (`held`, keys), and those of a result this pull wrote
# under another name. A column a file pull wrote is never removed.
.pruned <- function(data, previous, plans, held) {
  written <- vapply(plans, `[[`, "", "name")
  written_keys <- vapply(plans, function(plan) .pulled_key(plan$metadata), "")
  remove <- character()
  for (column in intersect(names(previous), .object_column_names(data))) {
    key <- .pulled_key(previous[[column]])
    if (is.na(key) || !identical(.pulled_from_of_key(key), .viewer_source) || column %in% written) next
    if (!key %in% held || key %in% written_keys) remove <- c(remove, column)
  }
  sort(remove, method = "radix")
}

# `name`, or "name (2)", "name (3)"... -- the first not in `taken`: the one
# rule for a gene set name used twice, in the web app, cellucid-python and
# here alike.
.unique_gene_set_name <- function(name, taken) {
  candidate <- name
  number <- 2L
  while (candidate %in% taken) {
    candidate <- paste0(name, " (", number, ")")
    number <- number + 1L
  }
  candidate
}

# Each gene set of `results`: its key, a name unique among them and its entry.
.gene_sets_of <- function(results) {
  sets <- list()
  taken <- character()
  for (result in results) {
    if (!identical(result$kind, "geneSet")) next
    name <- .unique_gene_set_name(result$meta$name, taken)
    taken <- c(taken, name)
    source <- result$source %||% .viewer_source
    sets[[length(sets) + 1L]] <- list(
      key = paste(source, "geneSet", result$id, sep = "\r"), name = name,
      entry = c(list(genes = as.character(unlist(result$meta$genes))),
                if (!is.null(result$meta$symbols)) list(symbols = as.character(unlist(result$meta$symbols))),
                list(description = result$meta$description, origin = result$meta$origin, pulled_from = source,
                     id = result$id))
    )
  }
  sets
}

# cellucid$gene_sets with `brought`: a set brought again replaces its earlier
# entry; with `held` (a viewer pull that prunes), a set of an earlier viewer
# pull the viewer no longer holds goes. Names are unique.
.merge_gene_sets <- function(existing, brought, held) {
  brought_keys <- vapply(brought, `[[`, "", "key")
  merged <- list()
  removed <- character()
  for (name in names(existing)) {
    entry <- existing[[name]]
    key <- if (is.list(entry)) .pulled_key(c(entry, list(kind = "geneSet"))) else NA_character_
    if (!is.na(key) && key %in% brought_keys) next
    if (!is.na(key) && !is.null(held) && identical(.pulled_from_of_key(key), .viewer_source) && !key %in% held) {
      removed <- c(removed, name)
      next
    }
    merged[[name]] <- entry
  }
  names_brought <- list()
  for (set in brought) {
    final <- .unique_gene_set_name(set$name, names(merged))
    merged[[final]] <- set$entry
    names_brought[[final]] <- set$entry$genes
  }
  list(sets = merged, removed = sort(removed, method = "radix"), brought = names_brought)
}

.pull_options <- function(highlights, fields, scores, visibility, gene_sets, include_pushed, highlight_prefix,
                          field_prefix, on_conflict, prune) {
  options <- list(
    highlights = highlights, fields = fields, scores = scores, visibility = visibility, gene_sets = gene_sets,
    include_pushed = include_pushed, prune = prune
  )
  for (name in names(options)) .validate_flag(options[[name]], name)
  for (name in c("highlight_prefix", "field_prefix")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L || is.na(value)) .abort(name, " must be one string.")
  }
  if (!is.character(on_conflict) || length(on_conflict) != 1L || !on_conflict %in% .pull_on_conflict) {
    .abort("on_conflict must be one of ", .format_value_list(.pull_on_conflict), ".")
  }
  c(options, list(highlight_prefix = highlight_prefix, field_prefix = field_prefix, on_conflict = on_conflict))
}

# Write the named columns, remove `remove`, merge the gene sets (NULL: leave
# them) and record what was written: the object and the summary.
.write_pulled <- function(data, plans, previous, remove, sets, held_sets, skipped) {
  data <- .object_with_columns(data, plans, remove)
  record <- .object_record(data)
  pulled <- previous[setdiff(intersect(names(previous), .object_column_names(data)), remove)]
  for (plan in plans) pulled[[plan$name]] <- plan$metadata
  record$pulled <- pulled
  merged <- list(removed = character(), brought = list())
  if (!is.null(sets)) {
    merged <- .merge_gene_sets(record$gene_sets %||% list(), sets, held_sets)
    record$gene_sets <- merged$sets
  }
  data <- .object_with_record(data, record)
  summary <- list(
    columns = stats::setNames(lapply(plans, `[[`, "metadata"), vapply(plans, `[[`, "", "name")),
    removed = remove,
    gene_sets = merged$brought,
    removed_gene_sets = merged$removed,
    skipped = skipped
  )
  list(data = data, summary = summary)
}

cellucid_pull <- function(
    data,
    source,
    ...,
    expected_dataset_id = NULL,
    highlights = TRUE,
    fields = TRUE,
    scores = TRUE,
    visibility = FALSE,
    gene_sets = TRUE,
    include_pushed = FALSE,
    highlight_prefix = "cellucid_highlight__",
    field_prefix = "",
    on_conflict = "error",
    prune = TRUE,
    return_summary = FALSE
) {
  if (missing(data) || missing(source)) {
    .abort("cellucid_pull() needs data (the object) and source (a viewer, or files the web app exported).")
  }
  .validate_flag(return_summary, "return_summary")
  options <- .pull_options(highlights, fields, scores, visibility, gene_sets, include_pushed, highlight_prefix,
                           field_prefix, on_conflict, prune)
  if (inherits(source, "cellucid_viewer")) {
    if (length(list(...)) > 0L) .abort("Arguments after source are for a session file; a viewer takes none.")
    .require_viewer_order(data, source)
    results <- cellucid_results(source)
    found <- .result_candidates(data, results, options)
    named <- .name_candidates(data, found$candidates, on_conflict)
    held <- vapply(results, function(result) paste(.viewer_source, result$kind, result$id, sep = "\r"), "")
    remove <- if (prune) .pruned(data, named$previous, named$plans, held) else character()
    written <- .write_pulled(data, named$plans, named$previous, remove, if (gene_sets) .gene_sets_of(results),
                             if (prune) held, found$skipped)
    return(if (return_summary) written else written$data)
  }
  paths <- .export_paths(source)
  files <- .results_of_files(paths, data)
  found <- .result_candidates(data, files$results, options)
  candidates <- found$candidates
  sets <- .gene_sets_of(files$results)
  extra <- list(...)
  for (session in files$sessions) {
    dataset_id <- expected_dataset_id %||% .object_record(data)$id
    if (is.null(dataset_id)) {
      .abort(basename(session), ": pass expected_dataset_id (the dataset id the session was saved on).")
    }
    # Planned against no existing column: its names are then settled with
    # every other column of this pull (on_conflict, reuse).
    planned <- .plan_session(data, session, extra, .validate_required_string(dataset_id, "expected_dataset_id"),
                             highlights, highlight_prefix, fields, field_prefix, FALSE, character())
    for (plan in planned$plans) {
      info <- c(plan$metadata, list(id = plan$metadata$group_id %||% plan$metadata$field_id,
                                    pulled_from = basename(session)))
      candidates[[length(candidates) + 1L]] <- list(key = .pulled_key(info), wanted = plan$name,
                                                     values = plan$values, metadata = info)
    }
    taken <- vapply(sets, `[[`, "", "name")
    for (entry in .session_gene_sets(planned$bundle)) {
      name <- .unique_gene_set_name(entry$name, taken)
      taken <- c(taken, name)
      sets[[length(sets) + 1L]] <- list(
        key = paste(basename(session), "geneSet", entry$name, sep = "\r"), name = name,
        entry = list(genes = entry$genes, description = entry$description, origin = "notebook",
                     pulled_from = basename(session), id = entry$name)
      )
    }
  }
  named <- .name_candidates(data, candidates, on_conflict)
  written <- .write_pulled(data, named$plans, named$previous, character(), if (gene_sets) sets, NULL, found$skipped)
  if (return_summary) written else written$data
}

# The gene sets a notebook pushed, saved with a session: each its name,
# description and genes.
.session_gene_sets <- function(bundle) {
  if (!.gene_sets_profile$id %in% names(bundle$chunks)) return(list())
  payload <- .parse_json_bytes(.session_chunk(bundle, .gene_sets_profile$id), .gene_sets_profile$id)
  lapply(payload, function(entry) {
    .session_check_keys(entry, c("name", "description", "genes"), "gene-sets/notebook entry")
    list(name = .session_text(entry$name, "gene-sets/notebook name"), description = entry$description,
         genes = as.character(unlist(entry$genes)))
  })
}

# --- The web app's exports ----------------------------------------------------

.export_paths <- function(source) {
  if (!is.character(source) || length(source) == 0L || anyNA(source)) {
    .abort("source must be a viewer, or the paths of files the web app exported (.csv, .gmt, .cellucid-session).")
  }
  paths <- character()
  for (item in path.expand(source)) {
    if (dir.exists(item)) {
      found <- list.files(item, pattern = "\\.(csv|gmt|cellucid-session)$", full.names = TRUE, ignore.case = TRUE)
      paths <- c(paths, sort(found, method = "radix"))
    } else if (file.exists(item)) {
      paths <- c(paths, item)
    } else {
      .abort(item, " does not exist.")
    }
  }
  if (length(paths) == 0L) .abort(.format_value_list(source), " holds no .csv, .gmt or .cellucid-session file to pull.")
  paths
}

.results_of_files <- function(paths, data) {
  results <- list()
  sessions <- character()
  cells <- .object_cells(data)
  for (path in paths) {
    lower <- tolower(basename(path))
    if (endsWith(lower, ".csv")) {
      results <- c(results, .read_obs_csv(path, cells))
    } else if (endsWith(lower, ".gmt")) {
      sets <- cellucid_read_gmt(path)
      for (name in names(sets)) {
        results[[length(results) + 1L]] <- list(
          kind = "geneSet", id = name, revision = length(results) + 1, n_obs = length(cells),
          meta = list(name = name, description = "", genes = as.list(sets[[name]]), origin = "user"),
          cells = numeric(0), values = NULL, value_type = NULL, source = basename(path)
        )
      }
    } else if (endsWith(lower, ".cellucid-session")) {
      sessions <- c(sessions, path)
    } else {
      .abort(basename(path), ": cellucid_pull() reads the web app's exports -- .csv (obs columns), .gmt (gene sets) ",
             "and .cellucid-session files.")
    }
  }
  list(results = results, sessions = sessions)
}

# What a column of an obs-column CSV holds (its kind row).
.csv_kinds <- c("highlight", "categorical", "continuous")

# An obs-column CSV the web app exports: a header row (obs_name, then each
# column's name), a kind row (kind, then each column's kind: highlight,
# categorical or continuous), then a row per cell. Each column is written as
# named, as its kind says -- never guessed from its values; a cell the file
# does not list is unknown (NA).
.read_obs_csv <- function(path, cells) {
  table <- utils::read.csv(path, header = FALSE, colClasses = "character", check.names = FALSE,
                           na.strings = character(), encoding = "UTF-8")
  file <- basename(path)
  if (nrow(table) < 2L || !identical(table[1L, 1L], .obs_name_column) || !identical(table[2L, 1L], "kind")) {
    .abort(file, ": an obs-column CSV starts with a header row whose first cell is '", .obs_name_column,
           "' and a kind row whose first cell is 'kind' (then each column's kind: ", paste(.csv_kinds, collapse = ", "), ").")
  }
  header <- as.character(table[1L, ])
  kinds <- as.character(table[2L, ])
  if (anyDuplicated(header)) .abort(file, ": a column is named twice.")
  for (index in seq_along(header)[-1L]) {
    if (!kinds[[index]] %in% .csv_kinds) {
      .abort(file, ": the column '", header[[index]], "' is of kind '", kinds[[index]], "', not one of ",
             paste(.csv_kinds, collapse = ", "), ".")
    }
  }
  rows <- table[-(1:2), , drop = FALSE]
  names_given <- rows[[1L]]
  if (anyDuplicated(names_given)) .abort(file, ": the cell '", names_given[anyDuplicated(names_given)], "' is listed twice.")
  if (anyDuplicated(cells)) .abort("This object's cell names are not unique.")
  positions <- match(names_given, cells)
  if (anyNA(positions)) {
    .abort(file, ": ", sum(is.na(positions)), " of its cells are not in this object (first: '",
           names_given[is.na(positions)][[1L]], "').")
  }
  order_of <- order(positions)
  results <- list()
  for (index in seq_along(header)[-1L]) {
    column <- header[[index]]
    kind <- kinds[[index]]
    text <- trimws(rows[[index]])
    present <- nzchar(text)
    base <- list(id = column, revision = index - 1, n_obs = length(cells), column = column, value_type = NULL,
                 source = file)
    if (identical(kind, "highlight")) {
      lowered <- tolower(text)
      unreadable <- present & !lowered %in% c("true", "false")
      if (any(unreadable)) {
        .abort(file, ": the highlight column '", column, "' holds '", text[unreadable][[1L]],
               "'; a highlight column holds true, false or nothing.")
      }
      results[[length(results) + 1L]] <- c(base, list(
        kind = "highlight", cells = sort(positions[present & lowered == "true"]), values = NULL,
        known = sort(positions[present]) - 1,
        meta = list(label = column, tool = "annotation", pageId = "csv", pageName = file, color = "#000000", enabled = TRUE)
      ))
    } else if (identical(kind, "continuous")) {
      numbers <- suppressWarnings(as.numeric(ifelse(present, text, NA)))
      unreadable <- present & (is.na(numbers) | is.infinite(numbers))
      if (any(unreadable)) {
        .abort(file, ": the continuous column '", column, "' holds '", text[unreadable][[1L]],
               "'; a continuous column holds finite numbers or nothing.")
      }
      results[[length(results) + 1L]] <- c(base, list(
        kind = "continuous", cells = positions[order_of], values = numbers[order_of],
        meta = list(key = column, origin = "user")
      ))
    } else {
      categories <- unique(text[present])
      codes <- ifelse(present, match(text, categories) - 1, 4294967295)
      base$value_type <- "uint32"
      results[[length(results) + 1L]] <- c(base, list(
        kind = "categorical", cells = positions[order_of], values = codes[order_of],
        meta = list(key = column, categories = as.list(categories), origin = "user")
      ))
    }
  }
  results
}

# --- GMT ------------------------------------------------------------------------

# Text a GMT cell can hold: each run of tabs and line breaks one space,
# trimmed (the web app's formatGmt and cellucid-python's rule).
.clean_gmt_text <- function(value) trimws(gsub("[\t\r\n]+", " ", value))

cellucid_write_gmt <- function(gene_sets, file, descriptions = NULL) {
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
    .abort("gene_sets must be a named list: each set's name and its genes.")
  }
  file <- .validate_required_string(file, "file")
  lines <- character(length(gene_sets))
  taken <- character()
  for (index in seq_along(gene_sets)) {
    name <- .clean_gmt_text(names(gene_sets)[[index]])
    if (!nzchar(name)) .abort("A gene set needs a name.")
    name <- .unique_gene_set_name(name, taken)
    taken <- c(taken, name)
    genes <- unique(.clean_gmt_text(as.character(gene_sets[[index]])))
    description <- if (is.null(descriptions[[names(gene_sets)[[index]]]])) "" else .clean_gmt_text(descriptions[[names(gene_sets)[[index]]]])
    lines[[index]] <- paste(c(name, description, genes[nzchar(genes)]), collapse = "\t")
  }
  con <- file(path.expand(file), "wb")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE)
  invisible(normalizePath(path.expand(file), winslash = "/"))
}

cellucid_read_gmt <- function(file) {
  file <- .validate_required_string(file, "file")
  lines <- readLines(path.expand(file), encoding = "UTF-8", warn = FALSE)
  sets <- list()
  for (line in lines) {
    if (!nzchar(trimws(line))) next
    cells <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    name <- trimws(cells[[1L]])
    if (!nzchar(name)) name <- paste0("Set ", length(sets) + 1L)
    name <- .unique_gene_set_name(name, names(sets))
    genes <- trimws(cells[-c(1L, 2L)])
    sets[[name]] <- unique(genes[nzchar(genes)])
  }
  sets
}

# --- Selections, subsets and gene sets of a viewer ----------------------------

# The cells of a highlight group (by id or label), the visible cells
# ("visible"), or the last selection (NULL), as a subset of `data`.
# Refuse an object whose cells are not the viewer's in order. The viewer
# names a cell by its original position: an object it was given, by the
# object's own column (viewer$cell_names, in that order); a file in a verified
# priority order, by the original position its record gives the row
# (.stored_order()). An object holding the same cells in another order --
# read from such a file as stored, or reordered since -- would take each
# result in the wrong cells. When the viewer's names are not known here (a
# file it serves from disk), an object in a file's stored order is still told
# by its priority record.
.require_viewer_order <- function(data, viewer) {
  # In a verified stored order that moves a cell: a record that does not
  # describe the cells leaves each its own original position. The names must
  # be the viewer's in order, repeated names included; the record tells an
  # object in a file's stored order whatever its names.
  positions <- .stored_order_of(data)$positions
  stored <- !is.null(positions) && !identical(as.double(positions), as.double(seq_along(positions) - 1))
  names <- viewer$cell_names
  given <- .object_cells(data)
  if (!is.null(names) && !is.null(given)) {
    if (length(given) != length(names)) return(invisible())
    if (identical(as.character(given), as.character(names)) && !stored) return(invisible())
    if (!anyDuplicated(names) && !anyDuplicated(given)) {
      missing <- is.na(match(given, names))
      if (any(missing)) {
        .abort(sum(missing), " of this object's cells are not cells of the viewer's dataset (first: '",
               given[missing][[1L]], "'): pass the object the viewer shows.")
      }
    }
  } else if (!stored) {
    return(invisible())
  }
  .abort(
    "This object holds the viewer's cells in another order",
    if (stored) " (a Cellucid-ready file's stored order: uns['cellucid']['priority'])",
    ", and the viewer names cells in their original order. ",
    if (stored) "Put them back first: cellucid_restore_original_order(object)."
    else "Reorder it first, to the viewer's cells in order."
  )
}

cellucid_subset <- function(data, viewer, which = NULL) {
  .require_viewer(viewer)
  .require_viewer_order(data, viewer)
  n <- length(.object_cells(data))
  if (is.null(which)) {
    selection <- cellucid_selection(viewer)
    if (is.null(selection)) .abort("The viewer has sent no selection yet; select cells, or name a highlight group.")
    cells <- selection$cells
  } else {
    which <- .validate_required_string(which, "which")
    results <- cellucid_results(viewer)
    matches <- Filter(function(result) {
      if (identical(which, "visible")) identical(result$kind, "visibility")
      else identical(result$kind, "highlight") && (identical(result$id, which) || identical(result$meta$label, which))
    }, results)
    if (length(matches) != 1L) {
      groups <- vapply(Filter(function(result) identical(result$kind, "highlight"), results),
                       function(result) paste0(result$id, " (", result$meta$label, ")"), "")
      .abort("'", which, "' names ", length(matches), " highlight groups; name one by id. The viewer holds: ",
             if (length(groups) == 0L) "no highlight group" else paste(groups, collapse = ", "))
    }
    if (matches[[1L]]$n_obs != n) .abort("The viewer shows ", matches[[1L]]$n_obs, " ", .plural(matches[[1L]]$n_obs, "cell"), " and this object has ", n, ".")
    cells <- matches[[1L]]$cells %||% seq_len(n)
  }
  if (is.data.frame(data)) data[cells, , drop = FALSE] else data[, cells]
}

# Each gene set's genes by identifier -- the var index, rownames(object) --
# or, with symbols = TRUE, by the names the viewer shows for them: those of
# the var column its data mapping names the genes by, else the identifiers.
cellucid_gene_sets <- function(viewer, symbols = FALSE) {
  .require_viewer(viewer)
  .validate_flag(symbols, "symbols")
  sets <- .gene_sets_of(cellucid_results(viewer))
  stats::setNames(lapply(sets, function(set) if (symbols) set$entry$symbols %||% set$entry$genes else set$entry$genes),
                  vapply(sets, `[[`, "", "name"))
}

# --- Pushing into the viewer --------------------------------------------------

# The setField command of one column: a named vector or factor aligned by its
# names to the viewer's cells (some or all), or one value per cell in order.
.field_command <- function(values, key, n, cell_names) {
  cells <- NULL
  given <- names(values)
  if (!is.null(given) && !is.null(cell_names) && !identical(given, cell_names)) {
    positions <- match(given, cell_names)
    if (anyNA(positions)) {
      .abort(sum(is.na(positions)), " of the column's names are not cells of the viewer's dataset (first: '",
             given[is.na(positions)][[1L]], "').")
    }
    if (anyDuplicated(positions)) .abort("The column names a cell twice.")
    order_of <- order(positions)
    values <- values[order_of]
    cells <- positions[order_of] - 1
    if (length(cells) == n) cells <- NULL
  } else if (length(values) != n) {
    .abort("The column has ", length(values), " ", .plural(length(values), "value"), " for the viewer's ", n, " ", .plural(n, "cell"),
           "; give one per cell in order, or name them by cell.")
  }
  if (is.logical(values)) {
    # A logical column's labels are "false" and "true", in every language.
    labels <- ifelse(values, "true", "false")
    values <- factor(labels, levels = intersect(c("false", "true"), labels))
  }
  if (is.factor(values) || is.character(values)) {
    factor_values <- if (is.factor(values)) values else factor(values)
    categories <- levels(factor_values)
    if (length(categories) == 0L) .abort("The column has no value to push.")
    type <- if (length(categories) < 255) "uint8" else if (length(categories) < 65535) "uint16" else "uint32"
    missing <- 2^(8 * .viewer_value_types[[type]]) - 1
    codes <- as.integer(factor_values) - 1
    codes[is.na(codes)] <- missing
    return(list(type = "setField", requestId = "", key = key, kind = "categorical", categories = as.list(categories),
                nObs = n, cells = cells, values = codes, valueType = type))
  }
  if (!is.numeric(values)) .abort("Push labels (a factor, characters, TRUE/FALSE) or numbers.")
  if (any(is.infinite(values))) .abort("The column holds an infinity; the viewer takes finite numbers (NA for no value).")
  list(type = "setField", requestId = "", key = key, kind = "continuous", categories = NULL, nObs = n,
       cells = cells, values = as.double(values), valueType = "float32")
}

# The cells of the dataset this viewer serves: read by R, or reported by the
# viewer (its ready event); NULL before either.
.dataset_cell_count <- function(viewer) viewer$n_cells %||% viewer$ready$n_obs

# A pushed column is sent at once, or -- for a format only the browser reads
# -- once the viewer has said how many cells it has (its ready event): a
# cell run before that ("Run All") is fine. Nothing waits for the viewer.
cellucid_push <- function(viewer, values, key = NULL) {
  .require_viewer(viewer)
  if (is.null(key)) .abort("Name the field in the viewer: key = \"...\".")
  key <- .validate_required_string(trimws(key), "key")
  n <- .dataset_cell_count(viewer)
  if (is.null(n)) {
    viewer$pending_pushes <- c(viewer$pending_pushes, list(list(values = values, key = key)))
    return(invisible(NULL))
  }
  command <- .field_command(values, key, n, viewer$cell_names)
  command$requestId <- paste0("p", .random_hex(6L))
  .viewer_send(viewer, command)
}

# The columns pushed before the viewer's cell count was known: each sent, or
# reported as a command_error of setField when it does not fit the dataset.
.send_pending_pushes <- function(viewer) {
  pending <- viewer$pending_pushes
  n <- .dataset_cell_count(viewer)
  if (length(pending) == 0L || is.null(n)) return(invisible(NULL))
  viewer$pending_pushes <- NULL
  for (push in pending) {
    command <- tryCatch(.field_command(push$values, push$key, n, viewer$cell_names),
                        error = function(error) conditionMessage(error))
    if (is.character(command)) {
      .report_local_command_error(viewer, "setField", command)
      next
    }
    command$requestId <- paste0("p", .random_hex(6L))
    .viewer_send(viewer, command)
  }
  invisible(NULL)
}

cellucid_push_gene_set <- function(viewer, name, genes, description = "") {
  .require_viewer(viewer)
  name <- .validate_required_string(name, "name")
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  .viewer_send(viewer, list(
    type = "setGeneSet", requestId = paste0("p", .random_hex(6L)), name = name, description = description,
    genes = as.list(genes)
  ))
}
