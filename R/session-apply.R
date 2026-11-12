# cellucid_apply_session(): what applying one step of a saved history
# (session format 4, R/session-file.R) to the object it was saved on would
# change -- a plan, one row per change -- or the change made: its labels
# (renamed obs columns and categories), its user-defined fields, its gene sets
# and its highlight groups, written as cell metadata (Seurat meta.data,
# SingleCellExperiment and SummarizedExperiment colData, a data frame's
# columns), and a record of each apply in the object's cellucid record
# (Seurat misc$cellucid$history_applied, metadata(x)$cellucid$history_applied,
# a data frame's "cellucid" attribute). The same plan has two other writers,
# rule for rule, row for row: cellucid-python's
# apply_cellucid_session_to_anndata() (anndata_session/_apply.py, _labels.py,
# _fields.py, _highlights.py, _cells.py) and the web app's "Save AnnData with
# this state..." copy (app/history/exports/anndata-copy.js); the three are
# compared with one shared expected plan, the var columns each apply writes
# included (tests/testthat/fixtures/history_golden/history_golden.expected.json).
#
# - Rows, in domain order: field-overlays (renamed columns, renamed or merged
#   categories, deleted fields, user-defined fields), gene-sets, highlights
#   (one TRUE/FALSE column per group, cellucid_highlight__<page>__<group>
#   from the page's name and the group's label), and one presentation row per
#   other domain whose value at the step differs from the tracker's opening
#   (filters, colors, views, the camera, analyses: how the data was shown,
#   never applied).
# - Statuses: write; unchanged (the object holds it already: applying again
#   changes nothing and says so); replace (a column or gene set an earlier
#   apply wrote, holding other values, with overwrite = TRUE); not-applied,
#   with the reason (a column the person made is never replaced).
# - The record of an apply (the history model's Provenance) holds the step's
#   expression values (expression: list(normalize, log1p), NULL when the
#   session holds none): what the person looked at the genes in, recorded and
#   never applied -- the object's expression values stay as stored.
# - Gene sets go to the cellucid record's gene_sets by name (origin
#   "session") and, one TRUE/FALSE column each, to the feature metadata of
#   the gene axis the data mapping reads: the modality it reads (a Seurat
#   assay's meta.features or Assay5 meta.data; rowData of the main
#   experiment, or of the altExp of that name). The column is
#   cellucid_gene_set__<name> (the highlight columns' rule, "gene_set" when
#   nothing is left; _2, _3 on clashes among the step's sets, named in its
#   order), TRUE on each row whose gene, as the mapping names genes, the set
#   lists; genes the dataset lacks are left out, with a warning. Only a set
#   written or replaced gets its column: a new one after the frame's columns,
#   an earlier apply's (replace) in place; a column of that name beside a set
#   no apply wrote is the dataset's own, never replaced (a warning says so).
#   A mapping that reads no genes, or the raw matrix's, has no var to mark:
#   the sets stay in the record (a warning says so). A set the step no
#   longer holds has no row: its record and its column stay. The var columns
#   are planned with the rest and are not rows of the plan.
# - Colors go with the columns (R/round-trip.R .colors_state), where the
#   object keeps what an AnnData keeps in uns, under Scanpy's name
#   <column>_colors: a categorical field written or replaced with the colors
#   it is shown in (.apply_shown_colors): in the view in focus at the step,
#   the color the person or an agent gave a category (the colors value), else
#   the color it opens with -- its palette (a field made from highlight
#   pages: its pages' colors), else Cellucid's own category colors by
#   position (R/category-colors.R); any other column written with none; a
#   renamed column's carried along; merged levels keeping their first level's
#   color. A dataset column recolored in the view in focus (a factor; its
#   name after the apply's renames, its levels after its relabels) takes the
#   colors it is shown in (.apply_column_colors): the color the person or an
#   agent gave a category, else the color it opens with -- the colors held
#   for it (.apply_opening_palette), a merged level taking its first source
#   level's, else Cellucid's own by position -- one column-colors row each
#   (domain colors), in the colors value's order, before the colors row. The
#   colors row stays not-applied (every other color is how the data was
#   shown), its note saying so when a field's or a column's colors went with
#   its column (.apply_colors_note).
# - obs['cellucid_original_position'] (R/priority.R) is never written,
#   renamed or replaced: a plan that would is refused before anything is
#   written, dry runs included (an object in a verified stored order is
#   refused before that).
# - What is checked, only what the chosen step needs: the dataset (cell
#   count, the data mapping's keys and gene count, the dataset id -- the one
#   the object was shown under (dataset_id = , as cellucid_view() was given
#   it; a viewer's request names the id it served the object under itself)
#   or else the one the object declares in its cellucid record, when it
#   names one -- the cell order when the object has the embedding the
#   session names) and the fields the chosen step's state uses.
#
# The state values read here are the web app's (app/state/domains/README.md
# and each module's header), each checked against its closed key set. A value
# names cells by original position; a step's state knows the cells the viewer
# had open then: a cell outside them is unknown, NA, never "not selected".

# The closed key sets of the state values a session is applied from.
.session_domain_keys <- list(
  cells = c("open", "genes"),
  overlays = c("renamedFields", "renamedCategories", "deletedFields", "userFields"),
  renamed_field = c("source", "key", "name"),
  renamed_category = c("source", "field", "category", "name"),
  deleted_field = c("source", "key", "purged"),
  user_field = c("id", "source", "kind", "key", "deleted", "purged", "createdAt", "sourceField", "operation"),
  user_category = c("categories", "codes", "sourcePages", "overlapStrategy", "overlapLabel", "intersectionLabels",
                    "uncoveredLabel", "palette"),
  user_values = "values",
  column = c("cells", "values"),
  source_field = c("sourceKey", "sourceIndex", "kind"),
  highlights = c("active", "pages"),
  page = c("id", "name", "color", "groups"),
  group = c("id", "type", "label", "enabled", "cells"),
  field_ref = c("source", "key"),
  range = c("min", "max"),
  gene_sets = "sets",
  gene_set = c("name", "description", "genes"),
  expression = c("normalize", "log1p")
)
.session_overlap_strategies <- c("first", "last", "overlap-label", "intersections")
.session_group_types <- c("annotation", "category", "combined", "knn", "lasso", "proximity", "range")
# A categorical field's code types: each one's largest value marks a cell
# with no category.
.session_code_types <- c(Uint8Array = 255, Uint16Array = 65535, Uint32Array = 4294967295)


# The arguments cellucid_view() was given that choose what is shown.
.session_view_arguments <- function(extra) {
  if (length(extra) > 0L && (is.null(names(extra)) || any(!nzchar(names(extra))))) {
    .abort("The arguments after the session must be named: the ones cellucid_view() was given.")
  }
  unknown <- setdiff(names(extra), .view_plan_arguments)
  if (length(unknown) > 0L) {
    .abort("cellucid_apply_session() takes the arguments cellucid_view() was given (",
           .format_value_list(.view_plan_arguments), "); not ", .format_value_list(unknown), ".")
  }
  extra
}

# --- Exact checks on the state values read --------------------------------------

.session_refuse <- function(label, ...) .abort("This session's ", label, " ", ...)

.session_record <- function(value, label) {
  if (!.is_json_object(value)) .session_refuse(label, "must be a record.")
  value
}

.session_list <- function(value, label) {
  if (!.is_json_array(value)) .session_refuse(label, "must be a list.")
  value
}

.session_keys <- function(value, expected, label) {
  .session_record(value, label)
  missing <- sort(setdiff(expected, names(value)), method = "radix")
  unknown <- sort(setdiff(names(value), expected), method = "radix")
  if (length(missing) > 0L || length(unknown) > 0L || anyDuplicated(names(value))) {
    details <- c(if (length(missing) > 0L) paste0("missing ", paste(missing, collapse = ", ")),
                 if (length(unknown) > 0L) paste0("unknown ", paste(unknown, collapse = ", ")))
    .session_refuse(label, "has other fields than this version reads (", paste(details, collapse = "; "), ").")
  }
  value
}

.session_bool <- function(value, label) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) .session_refuse(label, "must be true or false.")
  value
}

.session_text <- function(value, label) {
  if (!.is_text(value)) .session_refuse(label, "must be a non-empty text.")
  value
}

.session_optional_text <- function(value, label) if (is.null(value)) NULL else .session_text(value, label)

.session_source <- function(value, label) {
  if (!.is_json_string(value) || !value %in% c("obs", "var")) .session_refuse(label, "must be \"obs\" or \"var\".")
  value
}

# A field name as the viewer accepts one (app/utils/state-validator.js
# validateFieldKey): no surrounding white space, no control or zero-width
# character; a ':' is allowed (ATAC peaks).
.session_field_key <- function(value, label) {
  key <- .session_text(value, label)
  if (!.is_trimmed_text(key)) .session_refuse(label, "cannot have leading or trailing whitespace.")
  points <- .code_points(key)
  for (rule in list(
    list(points <= 0x1F | (points >= 0x7F & points <= 0x9F), "control"),
    list(points %in% c(0x200B, 0x2060, 0xFEFF), "zero-width")
  )) {
    if (any(rule[[1L]])) {
      .session_refuse(label, "cannot contain the ", rule[[2L]], " character U+", sprintf("%04X", points[rule[[1L]]][[1L]]), ".")
    }
  }
  key
}

# A category label: a string, finite number, or boolean (validateCategoryLabel).
.session_category_primitive <- function(value, label) {
  ok <- length(value) == 1L && !is.list(value) && !is.na(value) &&
    (is.character(value) || is.logical(value) || (is.numeric(value) && is.finite(value)))
  if (!ok) .session_refuse(label, "must be a string, finite number, or boolean.")
  value
}

# A category label's identity as JavaScript tells labels apart: 1 and "1"
# differ, 1 and 1.0 do not.
.session_category_identity <- function(value) {
  if (is.logical(value)) return(paste0("boolean:", if (value) "true" else "false"))
  if (is.character(value)) return(paste0("string:", value))
  paste0("number:", .js_numbers(as.double(value)))
}

# A category label as Python's str() writes it in a label's name
# ("renamedCategories:obs:cell_type:NK cell"), so both packages name it alike.
.session_category_text <- function(value) {
  if (is.logical(value)) return(if (value) "True" else "False")
  if (is.character(value)) return(value)
  if (is.integer(value)) as.character(value) else .js_numbers(value)
}

# A category label as the viewer shows it: JavaScript's String() of the value.
.category_label <- function(value) {
  if (is.logical(value)) return(if (value) "true" else "false")
  if (is.numeric(value)) return(.js_numbers(value))
  value
}

.session_category_list <- function(value, label) {
  items <- .session_list(value, label)
  if (length(items) == 0L) .session_refuse(label, "must list at least one category.")
  identities <- character()
  for (index in seq_along(items)) {
    exact <- .session_category_primitive(items[[index]], paste0(label, "[", index - 1L, "]"))
    identity <- .session_category_identity(exact)
    if (identity %in% identities) .session_refuse(label, "lists the category ", .format_value(exact), " twice.")
    identities <- c(identities, identity)
  }
  items
}

# A typed array of one of `types` (JavaScript type names), as a plain vector.
.session_typed <- function(value, types, label) {
  type <- attr(value, "cellucid_array")
  if (is.null(type) || !type %in% types) .session_refuse(label, "must be a ", paste(sort(types), collapse = " or "), ".")
  list(type = type, data = as.vector(value))
}

# Cells by original position: an ascending Uint32Array, each below `n`.
.session_cell_list <- function(value, n, label) {
  cells <- .session_typed(value, "Uint32Array", label)$data
  if (length(cells) > 1L && is.unsorted(cells, strictly = TRUE)) .session_refuse(label, "must list cells strictly ascending.")
  if (length(cells) > 0L && cells[[length(cells)]] >= n) {
    .session_refuse(label, "names the cell ", format(cells[[length(cells)]], scientific = FALSE), ", and the dataset has ", n, ".")
  }
  cells
}

# --- The cells a step knows ---------------------------------------------------------

# The cells (ascending 0-based original positions) a step's state knows; NULL
# for every cell: the session's own subsample (the manifest's kept cells),
# intersected with the subsample the `cells` domain says was open at the
# step, chosen again by the subsampling specification (R/subsample.R).
.session_known_cells <- function(saved_on, cells_value, view, n, step) {
  label <- paste0("cells value at step ", step)
  value <- .session_keys(cells_value, .session_domain_keys$cells, label)
  saved <- saved_on$kept
  if (is.null(value$open)) return(saved)
  open <- value$open
  if (!.is_json_object(open) || "nObs" %in% names(open)) {
    .session_refuse(label, "names its open cells in a form this version does not read.")
  }
  spec <- .validate_subsample_spec(c(open, list(nObs = n)), paste("Session", label))
  same <- function(a, b) {
    a <- a[setdiff(names(a), "nObs")]
    b <- b[setdiff(names(b), "nObs")]
    identical(a[order(names(a))], b[order(names(b))])
  }
  if (!is.null(saved_on) && same(spec, saved_on$spec)) return(saved)
  chosen <- tryCatch(
    as.double(switch(spec$method,
      uniform = cellucid_subsample_cells(n, n_cells = spec$k, method = "uniform", seed = spec$seed),
      blocks = cellucid_subsample_cells(n, n_cells = spec$k, method = "blocks", seed = spec$seed, block_size = spec$blockSize),
      stratified = {
        if (!spec$field %in% names(view$obs)) .abort("obs has no column '", spec$field, "'")
        cellucid_subsample_cells(n, n_cells = spec$k, method = "stratified", seed = spec$seed, groupby = view$obs[[spec$field]],
                                 allocation = spec$allocation, min_per_group = spec$minPerCategory)
      },
      spatial = {
        basis <- view$obsm[[spec$embedding]]
        if (is.null(basis)) .abort("obsm has no key '", spec$embedding, "'")
        cellucid_subsample_cells(n, n_cells = spec$k, method = "spatial", seed = spec$seed, basis = as.matrix(basis),
                                 columns = spec$columns + 1, grid_size = spec$gridSize)
      }
    )) - 1,
    error = function(error) {
      .abort("The step ", step, " was taken on a ", spec$method, " subsample of ", spec$k, " of the ", n, " ",
             .plural(n, "cell"), ", which cannot be chosen again on this object (", conditionMessage(error), "). ",
             "Apply a step taken on other cells, or the object the session was saved on.")
    }
  )
  chosen <- sort(chosen, method = "radix")
  if (is.null(saved)) chosen else intersect(saved, chosen)
}

# The renames and deletions of a field-overlays value, checked.
.session_label_entries <- function(overlays) {
  keys <- .session_domain_keys
  fields <- .session_list(overlays$renamedFields, "field-overlays value renamedFields")
  seen <- character()
  renamed_fields <- list()
  for (at in seq_along(fields)) {
    label <- paste0("field-overlays value renamedFields[", at - 1L, "]")
    entry <- .session_keys(fields[[at]], keys$renamed_field, label)
    out <- list(source = .session_source(entry$source, paste0(label, ".source")),
                key = .session_field_key(entry$key, paste0(label, ".key")),
                name = .session_field_key(entry$name, paste0(label, ".name")))
    identity <- paste(out$source, out$key, sep = "\r")
    if (identity %in% seen) .abort("This session's field-overlays value renames the ", out$source, " field '", out$key, "' twice.")
    seen <- c(seen, identity)
    renamed_fields[[length(renamed_fields) + 1L]] <- out
  }
  categories <- .session_list(overlays$renamedCategories, "field-overlays value renamedCategories")
  seen <- character()
  renamed_categories <- list()
  for (at in seq_along(categories)) {
    label <- paste0("field-overlays value renamedCategories[", at - 1L, "]")
    entry <- .session_keys(categories[[at]], keys$renamed_category, label)
    out <- list(source = .session_source(entry$source, paste0(label, ".source")),
                field = .session_field_key(entry$field, paste0(label, ".field")),
                category = .session_category_primitive(entry$category, paste0(label, ".category")),
                name = .session_category_primitive(entry$name, paste0(label, ".name")))
    identity <- paste(out$source, out$field, .session_category_identity(out$category), sep = "\r")
    if (identity %in% seen) {
      .abort("This session's field-overlays value renames the category ", .format_value(out$category), " of '", out$field, "' twice.")
    }
    seen <- c(seen, identity)
    renamed_categories[[length(renamed_categories) + 1L]] <- out
  }
  deleted <- .session_list(overlays$deletedFields, "field-overlays value deletedFields")
  deleted_fields <- lapply(seq_along(deleted), function(at) {
    label <- paste0("field-overlays value deletedFields[", at - 1L, "]")
    entry <- .session_keys(deleted[[at]], keys$deleted_field, label)
    list(source = .session_source(entry$source, paste0(label, ".source")),
         key = .session_field_key(entry$key, paste0(label, ".key")),
         purged = .session_bool(entry$purged, paste0(label, ".purged")))
  })
  list(renamed_fields = renamed_fields, renamed_categories = renamed_categories, deleted_fields = deleted_fields)
}

# A categorical column with categories renamed (merged where two read alike,
# in the place of the first), whether each rename found its category, and
# for each level after the place (1-based) of the first level before it was
# made from (its colors follow it). A factor's levels are text, so only a
# text label names one of them ("1" is not 1).
.renamed_factor <- function(values, renames) {
  levels <- levels(values)
  labels <- levels
  found <- logical(length(renames))
  for (at in seq_along(renames)) {
    category <- renames[[at]]$category
    hits <- if (is.character(category)) which(levels == category) else integer()
    labels[hits] <- .category_label(renames[[at]]$name)
    found[[at]] <- length(hits) > 0L
  }
  # A level renamed to a label another holds merges with it.
  levels(values) <- labels
  list(values = values, found = found, sources = match(levels(values), labels))
}

# A text column with the values each rename names replaced.
.renamed_text <- function(values, renames) {
  # A logical column with a label renamed is text: the others as the viewer
  # names them ("false", "true").
  out <- if (is.logical(values)) ifelse(values, "true", "false") else values
  found <- logical(length(renames))
  for (at in seq_along(renames)) {
    category <- renames[[at]]$category
    named <- is.character(category) || (is.logical(category) && is.logical(values))
    hits <- if (named) !is.na(values) & values == category else rep(FALSE, length(values))
    found[[at]] <- any(hits)
    out[hits] <- .category_label(renames[[at]]$name)
  }
  list(values = out, found = found)
}

# The cell table of the object: Seurat meta.data, colData, a data frame.
.object_cell_table <- function(data) {
  switch(.object_kind(data),
    seurat = methods::slot(data, "meta.data"),
    se = SummarizedExperiment::colData(data),
    data.frame = data
  )
}

# The colors (R/round-trip.R .colors_state) the object's columns hold once
# its relabeled columns' values are replaced and its columns renamed:
# `sources` gives, for each relabeled factor, each level's first level
# before; a text column has none.
.relabel_colors <- function(data, replaced, renamed, sources) {
  colors <- .colors_state(data)
  table <- .object_cell_table(data)
  for (column in names(replaced)) {
    from <- sources[[column]]
    if (is.null(from)) {
      colors$put(column, NULL)
      next
    }
    before <- length(levels(table[[column]]))
    if (length(from) == before) next
    held <- .text_colors(colors$held(column))
    colors$put(column, if (is.null(held) || length(held) < before) NULL else held[from])
  }
  moved <- lapply(names(renamed), function(old) .text_colors(colors$held(old)))
  for (old in names(renamed)) colors$put(old, NULL)
  for (at in seq_along(renamed)) colors$put(renamed[[at]], moved[[at]])
  colors
}

# The object with its relabeled columns' values replaced and its columns
# renamed in place (all at once: two columns may swap names), their colors
# following them (.relabel_colors), and each dataset column the session
# recolored (`recolor`, by its name after the renames: .apply_column_colors)
# given the colors it is shown in.
.object_with_relabels <- function(data, replaced, renamed, sources, recolor) {
  if (length(replaced) == 0L && length(renamed) == 0L && length(recolor) == 0L) return(data)
  colors <- .relabel_colors(data, replaced, renamed, sources)
  for (name in names(recolor)) colors$put(name, recolor[[name]])
  write <- function(table) {
    for (column in names(replaced)) table[[column]] <- replaced[[column]]
    original <- names(table)
    columns <- original
    for (old in names(renamed)) columns[original == old] <- renamed[[old]]
    names(table) <- columns
    table
  }
  data <- switch(.object_kind(data),
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
  .with_colors(data, colors)
}

# --- User-defined fields --------------------------------------------------------------

# A per-cell column {cells, values}: its cells (NULL: every cell), its values
# and their type.
.session_column_part <- function(value, types, n, label) {
  column <- .session_keys(value, .session_domain_keys$column, label)
  cells <- if (is.null(column$cells)) NULL else .session_cell_list(column$cells, n, paste0(label, ".cells"))
  values <- .session_typed(column$values, types, paste0(label, ".values"))
  expected <- if (is.null(cells)) n else length(cells)
  held <- length(values$data)
  if (held != expected) {
    .session_refuse(label, "holds ", held, " ", .plural(held, "value"), " for ",
                    if (is.null(cells)) "every one of the dataset's" else "its", " ", expected, " ", .plural(expected, "cell"), ".")
  }
  list(cells = cells, values = values$data, type = values$type)
}

.session_source_field <- function(value, label) {
  record <- .session_record(value, label)
  if (!"sourceKey" %in% names(record) || length(setdiff(names(record), .session_domain_keys$source_field)) > 0L) {
    .session_refuse(label, "must name sourceKey, and may name sourceIndex and kind only.")
  }
  out <- list(sourceKey = .session_text(record$sourceKey, paste0(label, ".sourceKey")))
  if ("sourceIndex" %in% names(record)) {
    index <- record$sourceIndex
    if (is.logical(index) || !.is_safe_integer(index) || index < 0) .session_refuse(label, ".sourceIndex must be a non-negative integer.")
    out$sourceIndex <- index
  }
  if ("kind" %in% names(record)) out$kind <- .session_text(record$kind, paste0(label, ".kind"))
  out
}

# One continuous field copied from another: an obs column, or a gene named
# the way the session's mapping names it, read from the matrix the view reads
# (the viewer reads the source's values, as here), as numbers (NA kept).
.continuous_alias_values <- function(view, mapping, source, source_field, label) {
  key <- source_field$sourceKey
  if (identical(source, "obs")) {
    index <- which(names(view$obs) == key)
    if (length(index) != 1L) .abort(label, " copies the obs column '", key, "', which must exist exactly once")
    values <- view$obs[[index]]
  } else {
    if (is.null(mapping$expression)) {
      .abort(label, " cannot copy gene values: the data mapping reads no expression matrix.")
    }
    index <- which(.session_gene_names(view, mapping) == key)
    if (length(index) != 1L) .abort(label, " copies the gene '", key, "', which must exist exactly once")
    values <- .session_gene_values(view, mapping, index)
  }
  declared <- source_field$sourceIndex
  if (!is.null(declared) && declared != index - 1L) {
    .abort(label, " sourceIndex ", declared, " does not match exact source '", key, "' at index ", index - 1L)
  }
  if (is.factor(values) || !(is.numeric(values) && !is.logical(values))) .abort(label, " must contain real numeric values")
  as.double(values)
}

# A column with every cell outside `known` (ascending 0-based positions;
# NULL: every cell) missing.
.without_unknown <- function(values, known, n) {
  if (is.null(known)) return(values)
  unknown <- rep(TRUE, n)
  unknown[known + 1] <- FALSE
  values[unknown] <- NA
  values
}

# Values over `cells` (0-based, NULL for every cell) placed among all `n`
# cells, NA elsewhere.
.placed_values <- function(values, cells, n) {
  if (is.null(cells)) return(values)
  out <- rep(NA_real_, n)
  out[cells + 1] <- values
  out
}

.session_group_metadata <- function(group, label) {
  keys <- .session_domain_keys
  type <- group$type
  if (!.is_json_string(type) || !type %in% .session_group_types) {
    .session_refuse(label, ".type is not a highlight tool this version reads.")
  }
  expected <- c(keys$group, if (identical(type, "category")) c("field", "category"), if (identical(type, "range")) c("field", "range"))
  .session_keys(group, expected, label)
  metadata <- list(label = .session_text(group$label, paste0(label, ".label")), type = type,
                   enabled = .session_bool(group$enabled, paste0(label, ".enabled")))
  if (type %in% c("category", "range")) {
    field <- .session_keys(group$field, keys$field_ref, paste0(label, ".field"))
    metadata$field <- list(source = .session_source(field$source, paste0(label, ".field.source")),
                           key = .session_field_key(field$key, paste0(label, ".field.key")))
  }
  if (identical(type, "category")) metadata$category <- .session_category_primitive(group$category, paste0(label, ".category"))
  if (identical(type, "range")) {
    bounds <- .session_keys(group$range, keys$range, paste0(label, ".range"))
    for (bound in c("min", "max")) {
      if (is.logical(bounds[[bound]]) || !.is_finite_number(bounds[[bound]])) {
        .session_refuse(label, ".range.", bound, " must be a finite number.")
      }
    }
    if (bounds$min > bounds$max) .session_refuse(label, ".range has its min above its max.")
    metadata$range <- list(min = bounds$min, max = bounds$max)
  }
  metadata
}

# --- The words a plan uses (the same as cellucid-python's) --------------------------

.apply_presentation_note <- "Presentation: not applied to the object."
.apply_colors_note <- paste0(
  "The colors of the user fields written and of the dataset columns recolored go with their columns ",
  "(uns['<column>_colors']); every other color is how the data was shown: not applied to the object."
)
.apply_modality_colors_note <- function(column) paste0("\"", column, "\" is a column of a MuData modality; colors apply to the global obs.")
.apply_not_categorical_colors_note <- function(column) {
  paste0("The obs column \"", column, "\" is not categorical; Scanpy keeps colors for a categorical column only.")
}
.apply_missing_categories_note <- function(column, labels) {
  paste0("Not written: the ", .plural(length(labels), "color"), " of ", paste0("\"", labels, "\"", collapse = ", "),
         ", which the obs column \"", column, "\" does not have.")
}
.apply_gene_name_note <- "A gene's display name is recorded, never applied to the object."
.apply_gene_category_note <- "A gene has no categories."
.apply_deletion_note <- "A deletion is recorded, never applied to the object: the column stays."
.apply_deleted_note <- "Deleted in the session."
.apply_purged_note <- "Deleted for good in the session: it holds no values."
.apply_absent_note <- function(column) paste0("The object has no obs column \"", column, "\".")
.apply_modality_note <- function(column) paste0("\"", column, "\" is a column of a MuData modality; labels apply to the global obs.")
.apply_taken_note <- function(column) paste0("The object already has an obs column \"", column, "\" that Cellucid did not write; it is never replaced.")
.apply_numbers_note <- function(column) paste0("The obs column \"", column, "\" holds numbers, not categories.")
.apply_no_category_note <- function(column, category) paste0("The obs column \"", column, "\" has no category \"", .category_label(category), "\".")
.apply_owned_note <- function(name, step) paste0("An earlier apply (step ", step, ") wrote \"", name, "\" with other values; overwrite replaces it.")
.apply_taken_set_note <- function(name) paste0("The object already has a gene set \"", name, "\" that no session wrote; it is never replaced.")
.apply_steps_only <- paste0(
  "This session was shared as its steps only, without the values they leave behind, so it cannot be applied to an object; ",
  "replay it in Cellucid on its dataset instead."
)
.apply_data_domains <- c("cells", "field-overlays", "gene-sets", "highlights")

# The expression values a step's state holds (its "expression" value), as the
# provenance records them (the history model's ExpressionValues):
# list(normalize, log1p), or NULL when the session holds none. Recorded so a
# notebook knows what the person looked at the genes in; never applied (the
# object's values stay as stored).
.apply_expression_values <- function(value, step) {
  if (is.null(value)) return(NULL)
  label <- paste0("expression value at step ", step)
  .session_keys(value, .session_domain_keys$expression, label)
  list(normalize = .session_bool(value$normalize, paste0(label, "'s normalize")),
       log1p = .session_bool(value$log1p, paste0(label, "'s log1p")))
}

.apply_row <- function(kind, domain, name, status, column = NA_character_, from = character(), cells = NA_integer_, note = NA_character_) {
  list(kind = kind, domain = domain, name = name, column = column, from = from, cells = cells, status = status, note = note)
}

# --- Choosing the step ---------------------------------------------------------------

# The tracker, branch and step an apply takes: the step asked for (on the
# branch asked for, when both are), a branch's tip, or the tracker's head on
# its active branch.
.apply_choose <- function(session, tracker, branch, step) {
  found <- cellucid_tracker(session, tracker)
  if (!is.null(step)) {
    chosen <- cellucid_step(found, step)
    if (is.null(branch)) return(list(tracker = found, branch = chosen$branch, step = chosen))
    on <- cellucid_branch(found, branch)
    path <- vapply(.history_path(found, on$tip), function(item) item$id, "")
    if (!chosen$id %in% path) {
      .abort("Step ", chosen$id, " is not on the branch \"", on$name, "\" (its path runs from the root to step ", on$tip,
             "); it was made on \"", cellucid_branch(found, chosen$branch)$name, "\".")
    }
    return(list(tracker = found, branch = on$id, step = chosen))
  }
  if (!is.null(branch)) {
    on <- cellucid_branch(found, branch)
    return(list(tracker = found, branch = on$id, step = cellucid_step(found, on$tip)))
  }
  list(tracker = found, branch = found$active_branch, step = cellucid_step(found, found$head_id))
}

# --- The dataset the session was saved on --------------------------------------------

# The cell-order digest of some cells' stored embedding values
# (data/cell-order-digest.js): two MurmurHash3 lanes per cell over its
# original position and its float32 values on up to three axes, summed,
# finalized with the cell count; 16 lower-case hex digits. 32-bit unsigned
# words are held as doubles (exact below 2^53), their bit operations done on
# 16-bit halves (R's integers are signed and hold no -2^31).
.apply_cell_order_digest <- function(values, positions = NULL) {
  values <- as.matrix(values)
  cells <- nrow(values)
  width <- ncol(values)
  if (width < 1L || width > 3L) .abort("A cell-order digest is taken over 1, 2 or 3 columns.")
  two32 <- 4294967296
  u32 <- function(x) x %% two32
  xor32 <- function(a, b) {
    bitwXor(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 + bitwXor(as.integer(a %% 65536), as.integer(b %% 65536))
  }
  mul32 <- function(a, k) u32(a * (k %% 65536) + ((a * (k %/% 65536)) %% 65536) * 65536)
  rotl <- function(a, r) u32(a * 2^r) + a %/% 2^(32 - r)
  shr <- function(a, r) a %/% 2^r
  scramble <- function(word) mul32(rotl(mul32(word, 0xCC9E2D51), 15L), 0x1B873593)
  step <- function(h, k) u32(mul32(rotl(xor32(h, k), 13L), 5) + 0xE6546B64)
  fmix <- function(h) {
    h <- xor32(h, shr(h, 16L))
    h <- mul32(h, 0x85EBCA6B)
    h <- xor32(h, shr(h, 13L))
    h <- mul32(h, 0xC2B2AE35)
    xor32(h, shr(h, 16L))
  }
  bits <- function(column) {
    raw <- writeBin(as.double(column), raw(), size = 4L, endian = "little")
    word <- readBin(raw, "integer", n = length(column), size = 4L, endian = "little")
    word <- ifelse(is.na(word), two32 / 2, ifelse(word < 0, word + two32, as.double(word)))
    nan <- (word %/% 8388608) %% 256 == 255 & word %% 8388608 != 0
    word[nan] <- 0x7FC00000
    word
  }
  zero <- scramble(rep(0, cells))
  words <- lapply(1:3, function(axis) if (axis <= width) scramble(bits(values[, axis])) else zero)
  position <- if (is.null(positions)) seq_len(cells) - 1 else as.double(positions)
  scrambled <- scramble(u32(position))
  lane <- function(seed) {
    h <- rep(seed, cells)
    h <- step(h, scrambled)
    for (word in words) h <- step(h, word)
    mixed <- fmix(xor32(h, rep(16, cells)))
    total <- 0
    for (start in seq(1L, max(cells, 1L), by = 1048576L)) {
      end <- min(cells, start + 1048575L)
      if (end >= start && cells > 0L) total <- u32(total + sum(mixed[start:end]))
    }
    total
  }
  low <- lane(0x9E3779B1)
  high <- lane(0x85EBCA6B)
  first <- fmix(xor32(high, fmix(xor32(cells, 0x297A2D39))))
  second <- fmix(xor32(low, fmix(xor32(cells, 0x2C1B3C6D))))
  hex <- function(value) sprintf("%04x%04x", as.integer(value %/% 65536), as.integer(value %% 65536))
  paste0(hex(first), hex(second))
}

# The cell-order digest of the embedding the session names, over the cells it
# was saved on, when the object has that embedding.
.apply_check_cell_order <- function(fingerprint, view, kept) {
  order <- fingerprint$cellOrder
  key <- order$embedding$key
  columns <- as.integer(unlist(order$embedding$columns))
  values <- view$obsm[[key]]
  if (is.null(values)) return(invisible())
  values <- as.matrix(values)
  if (max(columns) >= ncol(values)) return(invisible())
  rows <- if (is.null(kept)) values[, columns + 1L, drop = FALSE] else values[kept + 1, columns + 1L, drop = FALSE]
  digest <- .apply_cell_order_digest(rows, kept)
  if (!identical(digest, order$digest)) {
    .abort("This object's cells are not the cells the session was saved on, in the same order: the cell-order digest of obsm['",
           key, "'] (columns ", paste(columns, collapse = ", "), ") is ", digest, ", and the session's is ", order$digest,
           ". Put the cells in their original order, or apply the session to the object it was saved on.")
  }
}

# The dataset fields a step's state uses (collectRequiredFields): every obs
# column and gene a view colors by or filters on, or a highlight group is
# bound to; fields the session made itself are not listed.
.apply_required_fields <- function(state) {
  obs <- character()
  var <- character()
  add <- function(ref) {
    if (!.is_json_object(ref) || !.is_json_string(ref$key)) return()
    if (identical(ref$source, "obs")) obs <<- c(obs, ref$key)
    if (identical(ref$source, "var")) var <<- c(var, ref$key)
  }
  for (domain in c("filters", "colors")) for (view in state[[domain]]$views) for (entry in view) add(entry$field)
  for (ref in state[["active-field"]]$views) add(ref)
  for (page in state$highlights$pages) for (group in page$groups) add(group$field)
  order <- function(keys) {
    keys <- unique(keys)
    keys[.utf16_order(keys)]
  }
  list(obs = order(obs), var = order(var))
}

.apply_require_fields <- function(view, requires, mapping, step) {
  missing_obs <- setdiff(requires$obs, names(view$obs))
  missing_var <- if (length(requires$var) > 0L) setdiff(requires$var, .session_gene_names(view, mapping)) else character()
  if (length(missing_obs) == 0L && length(missing_var) == 0L) return(invisible(NULL))
  quote <- function(keys) paste0("'", keys, "'", collapse = ", ")
  parts <- c(
    if (length(missing_obs) > 0L) paste0("the obs column", if (length(missing_obs) > 1L) "s", " ", quote(missing_obs)),
    if (length(missing_var) > 0L) paste0("the gene", if (length(missing_var) > 1L) "s", " ", quote(missing_var))
  )
  .abort("The state at step ", step, " uses ", paste(parts, collapse = " and "), ", which this dataset does not have.")
}

# --- What earlier applies wrote --------------------------------------------------------

cellucid_history_applied <- function(data) {
  record <- .object_record(data)
  entries <- record$history_applied
  if (is.null(entries)) return(list())
  if (!is.list(entries) || !is.null(names(entries))) .abort("The object's cellucid$history_applied must be a list of apply records.")
  entries
}

.apply_same_column <- function(existing, planned) {
  strip <- function(values) {
    attributes(values) <- attributes(values)[intersect(names(attributes(values)), c("levels", "class"))]
    values
  }
  identical(strip(existing), strip(planned))
}

# --- Labels ------------------------------------------------------------------------------

.apply_labels_of <- function(values) {
  if (is.factor(values)) return(levels(values))
  if (is.character(values)) return(unique(values[!is.na(values)]))
  # A logical column's categories are false and true, as the viewer reads it.
  if (is.logical(values)) return(c(FALSE, TRUE)[c(FALSE, TRUE) %in% values])
  NULL
}

# The renames and deletions of a field-overlays value as rows against the
# object's cells; `earlier` lists the renames earlier applies recorded.
.apply_plan_labels <- function(overlays, data, earlier, modality_columns) {
  entries <- .session_label_entries(overlays)
  table <- .object_cell_table(data)
  present <- names(table)
  rows <- list()
  rename <- stats::setNames(character(), character())
  renamed <- list()
  earlier_columns <- Filter(function(entry) identical(entry$kind, "column"), earlier)
  earlier_categories <- Filter(function(entry) identical(entry$kind, "category"), earlier)
  absent <- function(column) {
    if (column %in% modality_columns && !column %in% present) return(.apply_modality_note(column))
    if (!column %in% present) return(.apply_absent_note(column))
    NULL
  }
  applicable <- Filter(function(entry) identical(entry$source, "obs") && !identical(entry$name, entry$key) && entry$key %in% present,
                       entries$renamed_fields)
  kept <- setdiff(present, vapply(applicable, function(entry) entry$key, ""))
  targets <- character()
  for (entry in entries$renamed_fields) {
    key <- entry$key
    name <- entry$name
    if (identical(entry$source, "var")) {
      rows[[length(rows) + 1L]] <- .apply_row("column-renamed", "field-overlays", name, "not-applied", from = key, note = .apply_gene_name_note)
      next
    }
    if (identical(name, key)) next
    reason <- absent(key)
    if (is.null(reason)) {
      if (name %in% kept || name %in% targets) {
        rows[[length(rows) + 1L]] <- .apply_row("column-renamed", "field-overlays", name, "not-applied", from = key, note = .apply_taken_note(name))
        next
      }
      targets <- c(targets, name)
      rename[[key]] <- name
      renamed[[length(renamed) + 1L]] <- list(kind = "column", column = NULL, from = key, to = name)
      rows[[length(rows) + 1L]] <- .apply_row("column-renamed", "field-overlays", name, "write", from = key)
    } else if (name %in% present && any(vapply(earlier_columns, function(item) identical(item$from, key) && identical(item$to, name), TRUE))) {
      rows[[length(rows) + 1L]] <- .apply_row("column-renamed", "field-overlays", name, "unchanged", from = key)
    } else {
      rows[[length(rows) + 1L]] <- .apply_row("column-renamed", "field-overlays", name, "not-applied", from = key, note = reason)
    }
  }
  earlier_names <- stats::setNames(vapply(earlier_columns, function(item) item$to, ""), vapply(earlier_columns, function(item) item$from, ""))
  targets_of <- list()
  for (entry in entries$renamed_categories) {
    if (identical(entry$source, "obs")) targets_of[[entry$field]] <- c(targets_of[[entry$field]], .session_category_identity(entry$name))
  }
  by_column <- list()
  for (entry in entries$renamed_categories) {
    field <- entry$field
    category <- entry$category
    name <- .category_label(entry$name)
    from <- .category_label(category)
    if (identical(entry$source, "var")) {
      rows[[length(rows) + 1L]] <- .apply_row("categories-renamed", "field-overlays", name, "not-applied", column = field, from = from,
                                              note = .apply_gene_category_note)
      next
    }
    before <- if (field %in% present || is.na(earlier_names[field])) field else earlier_names[[field]]
    after <- if (before %in% names(rename)) rename[[before]] else before
    reason <- absent(before)
    if (!is.null(reason)) {
      rows[[length(rows) + 1L]] <- .apply_row("categories-renamed", "field-overlays", name, "not-applied", column = field, from = from, note = reason)
      next
    }
    values <- table[[before]]
    labels <- .apply_labels_of(values)
    if (is.null(labels)) {
      rows[[length(rows) + 1L]] <- .apply_row("categories-renamed", "field-overlays", name, "not-applied", column = after, from = from,
                                              note = .apply_numbers_note(after))
      next
    }
    identities <- vapply(labels, .session_category_identity, "")
    has_from <- .session_category_identity(category) %in% identities
    has_to <- .session_category_identity(entry$name) %in% identities
    shared <- sum(targets_of[[field]] == .session_category_identity(entry$name)) > 1L
    merged <- shared || (has_from && has_to && !identical(.session_category_identity(category), .session_category_identity(entry$name)))
    kind <- if (merged) "categories-merged" else "categories-renamed"
    if (has_from) {
      by_column[[before]] <- c(by_column[[before]], list(entry))
      renamed[[length(renamed) + 1L]] <- list(kind = "category", column = after, from = from, to = name)
      shown <- if (is.logical(values)) ifelse(values, "true", "false") else as.character(values)
      count <- sum(!is.na(values) & shown == from)
      rows[[length(rows) + 1L]] <- .apply_row(kind, "field-overlays", name, "write", column = after, from = from, cells = as.integer(count))
    } else if (has_to && any(vapply(earlier_categories, function(item) {
      identical(item$column, after) && identical(item$from, from) && identical(item$to, name)
    }, TRUE))) {
      rows[[length(rows) + 1L]] <- .apply_row(kind, "field-overlays", name, "unchanged", column = after, from = from)
    } else {
      rows[[length(rows) + 1L]] <- .apply_row(kind, "field-overlays", name, "not-applied", column = after, from = from,
                                              note = .apply_no_category_note(after, category))
    }
  }
  replaced <- list()
  sources <- list()
  for (column in names(by_column)) {
    values <- table[[column]]
    if (is.factor(values)) {
      renamed_factor <- .renamed_factor(values, by_column[[column]])
      replaced[[column]] <- renamed_factor$values
      sources[[column]] <- renamed_factor$sources
    } else {
      replaced[[column]] <- .renamed_text(values, by_column[[column]])$values
    }
  }
  for (entry in entries$deleted_fields) {
    rows[[length(rows) + 1L]] <- .apply_row("field-deleted", "field-overlays", entry$key, "not-applied", note = .apply_deletion_note)
  }
  list(rows = rows, rename = rename, replaced = replaced, sources = sources, renamed = renamed, earlier_names = earlier_names)
}

# --- User-defined fields and highlights as candidate columns -----------------------------

# The colors a step shows (cellucid-python's anndata_session/_colors.py, by
# the same rule). The colors value (app/state/domains/colors.js) holds, per
# view (live first, then kept views), the fields whose colors differ from
# their opening colors; a categorical field's entry lists each category
# recolored ({category, color} by label). The views value names the view in
# focus (focus, its position; 0 is the live view). `colors` and `views` are
# the step's values (NULL when the session holds none: nothing recolored, the
# live view in focus).
#
# Each field the view in focus recolors categories of, in the value's order:
# list(source, key (NULL for a user field), user (its id, NULL for a dataset
# field), recolored: each category's identity (.session_category_identity) to
# its "#rrggbb", labels: each identity to the category as the viewer shows
# it).
.apply_focused_entries <- function(colors, views) {
  if (is.null(colors)) return(list())
  record <- .session_record(colors, "colors value")
  per_view <- .session_list(record[["views"]], "colors value views")
  focus <- 0
  if (!is.null(views)) {
    shown <- .session_record(views, "views value")
    focus <- shown[["focus"]]
    if (is.logical(focus) || !is.numeric(focus) || length(focus) != 1L || is.na(focus) || focus < 0 || focus != floor(focus)) {
      .session_refuse("views value focus", "must be a non-negative whole number.")
    }
  }
  if (focus >= length(per_view)) {
    .session_refuse("colors value", "holds ", length(per_view), " ", .plural(length(per_view), "view"), ", and the view in focus is view ", focus, ".")
  }
  entries <- .session_list(per_view[[focus + 1]], paste0("colors value views[", focus, "]"))
  out <- list()
  for (position in seq_along(entries)) {
    label <- paste0("colors value views[", focus, "][", position - 1L, "]")
    entry <- .session_record(entries[[position]], label)
    field <- .session_record(entry[["field"]], paste0(label, ".field"))
    source <- .session_source(field[["source"]], paste0(label, ".field.source"))
    user <- field[["user"]]
    key <- field[["key"]]
    if (!is.null(user)) {
      if (!.is_text(user)) .session_refuse(paste0(label, ".field.user"), "must be a non-empty text.")
    } else if (!.is_text(key)) {
      .session_refuse(paste0(label, ".field.key"), "must be a non-empty text.")
    }
    if (!"colors" %in% names(entry)) next
    recolored <- list()
    labels <- character()
    pairs <- .session_list(entry[["colors"]], paste0(label, ".colors"))
    for (index in seq_along(pairs)) {
      at <- paste0(label, ".colors[", index - 1L, "]")
      pair <- .session_record(pairs[[index]], at)
      if (!setequal(names(pair), c("category", "color")) || length(pair) != 2L) .session_refuse(at, "must name category and color only.")
      category <- .session_category_primitive(pair[["category"]], paste0(at, ".category"))
      color <- pair[["color"]]
      if (!.is_json_string(color) || !grepl("^#[0-9a-f]{6}$", color)) .session_refuse(paste0(at, ".color"), "must be one #rrggbb color.")
      identity <- .session_category_identity(category)
      recolored[[identity]] <- color
      labels[[identity]] <- .category_label(category)
    }
    out[[length(out) + 1L]] <- list(source = source, key = if (is.null(user)) key else NULL, user = user,
                                    recolored = recolored, labels = labels)
  }
  out
}

# The categories of user fields recolored in the view in focus, by
# "<source>\u0001<field id>": each category's identity to its "#rrggbb".
.apply_focused_recolors <- function(colors, views) {
  out <- list()
  for (entry in .apply_focused_entries(colors, views)) {
    if (!is.null(entry$user)) out[[paste0(entry$source, "\u0001", entry$user)]] <- entry$recolored
  }
  out
}

# The categories of the dataset's obs columns recolored in the view in focus,
# by the column's key as the dataset names it (a column renamed in the
# session keeps its key), in the value's order: each an entry of
# .apply_focused_entries. A gene's colors (a color map, a range) are never
# categories: not here.
.apply_focused_column_recolors <- function(colors, views) {
  out <- list()
  for (entry in .apply_focused_entries(colors, views)) {
    if (is.null(entry$user) && identical(entry$source, "obs") && length(entry$recolored) > 0L) out[[entry$key]] <- entry
  }
  out
}

# The colors a categorical user field's categories are shown in, one
# "#rrggbb" a category in its order: the view in focus's recolors, else its
# palette, else Cellucid's own colors by position.
.apply_shown_colors <- function(categories, palette, recolored) {
  opening <- if (is.null(palette)) .default_category_colors(length(categories)) else palette
  vapply(seq_along(categories), function(index) {
    given <- recolored[[.session_category_identity(categories[[index]])]]
    if (is.null(given)) opening[[index]] else given
  }, "")
}

# A file's color as the viewer draws it: "#rrggbb" in lowercase (#rgb
# widened, the alpha of #rrggbbaa dropped).
.apply_drawn_color <- function(text) {
  digits <- substring(text, 2L)
  digits <- ifelse(nchar(digits) == 3L, gsub("(.)", "\\1\\1", digits), substring(digits, 1L, 6L))
  tolower(paste0("#", digits))
}

# The colors the `count` categories (factor levels) of a dataset column open
# with in the viewer from the colors held for it (`held`, its
# <column>_colors): the first `count`, as drawn, when there are that many and
# each is a hex color (R/dataset-plan.R .hex_file_color; palettes.js
# scanpyCategoryPalette); NULL otherwise, and the column opens in Cellucid's
# own colors by position.
.apply_opening_palette <- function(held, count) {
  colors <- .text_colors(held)
  if (is.null(colors) || count == 0L || length(colors) < count) return(NULL)
  first <- colors[seq_len(count)]
  if (!all(grepl(.hex_file_color, first))) return(NULL)
  .apply_drawn_color(first)
}

# Each dataset column the view in focus recolored categories of
# (.apply_focused_column_recolors, by the key the dataset names it), written
# with the colors it is shown in: each level, in the column's order once the
# apply's labels are applied, takes the color the person or an agent gave it,
# else the color it opens with -- the colors held for the column
# (.apply_opening_palette; a level merged by the apply takes its first source
# level's), else Cellucid's own colors by position. A row (column-colors,
# domain colors, in the colors value's order) is write when those differ
# from the colors the column holds after the apply's other changes,
# unchanged when they are the same, not-applied (with the reason) for a
# column the object lacks, a MuData modality's column, a column that is not a
# factor (Scanpy keeps colors for a categorical column only), or one with
# none of the categories recolored. `recolor`: the colors to write, by the
# column's name after the apply's renames.
.apply_column_colors <- function(recolored, data, labels, modality_columns) {
  table <- .object_cell_table(data)
  present <- names(table)
  would_hold <- .relabel_colors(data, labels$replaced, labels$rename, labels$sources)
  rows <- list()
  recolor <- list()
  row <- function(name, status, note = NA_character_) .apply_row("column-colors", "colors", name, status, note = note)
  for (key in names(recolored)) {
    given <- recolored[[key]]
    before <- if (key %in% present || is.na(labels$earlier_names[key])) key else labels$earlier_names[[key]]
    after <- if (before %in% names(labels$rename)) labels$rename[[before]] else before
    if (!before %in% present) {
      note <- if (key %in% modality_columns) .apply_modality_colors_note(key) else .apply_absent_note(key)
      rows[[length(rows) + 1L]] <- row(after, "not-applied", note)
      next
    }
    values <- table[[before]]
    if (!is.factor(values)) {
      rows[[length(rows) + 1L]] <- row(after, "not-applied", .apply_not_categorical_colors_note(after))
      next
    }
    count <- nlevels(values)
    relabeled <- labels$replaced[[before]]
    categories <- if (is.null(relabeled)) levels(values) else levels(relabeled)
    sources <- if (is.null(labels$sources[[before]])) seq_len(count) else labels$sources[[before]]
    opening <- .apply_opening_palette(.object_uns(data, .colors_key(before)), count)
    if (is.null(opening)) opening <- .default_category_colors(count)
    identities <- vapply(categories, .session_category_identity, "", USE.NAMES = FALSE)
    shown <- vapply(seq_along(categories), function(at) {
      color <- given$recolored[[identities[[at]]]]
      if (is.null(color)) opening[[sources[[at]]]] else color
    }, "")
    asked <- names(given$recolored)
    missing <- unname(given$labels[asked[!asked %in% identities]])
    note <- if (length(missing) > 0L) .apply_missing_categories_note(after, missing) else NA_character_
    if (length(missing) == length(asked)) {
      rows[[length(rows) + 1L]] <- row(after, "not-applied", note)
      next
    }
    if (identical(.text_colors(would_hold$held(after)), shown)) {
      rows[[length(rows) + 1L]] <- row(after, "unchanged", note)
      next
    }
    recolor[[after]] <- shown
    rows[[length(rows) + 1L]] <- row(after, "write", note)
  }
  list(rows = rows, recolor = recolor)
}

# One candidate column per user-defined field the state holds, every field
# checked before any is built: its values (NULL when the session does not
# apply it: a field deleted or purged in it), the cells it labels.
.apply_user_fields <- function(context, fields) {
  keys <- .session_domain_keys
  n <- context$n
  ids <- character()
  out <- list()
  fields <- .session_list(fields, "field-overlays value userFields")
  for (position in seq_along(fields)) {
    label <- paste0("field-overlays value userFields[", position - 1L, "]")
    record <- .session_record(fields[[position]], label)
    kind <- record$kind
    if (!.is_json_string(kind) || !kind %in% c("category", "continuous")) .session_refuse(label, ".kind must be \"category\" or \"continuous\".")
    purged <- .session_bool(record$purged, paste0(label, ".purged"))
    copy_of <- identical(kind, "continuous") && !is.null(record$sourceField)
    expected <- if (purged || copy_of) keys$user_field else if (identical(kind, "category")) c(keys$user_field, keys$user_category) else c(keys$user_field, keys$user_values)
    field <- .session_keys(record, expected, label)
    id <- .session_text(field$id, paste0(label, ".id"))
    if (id %in% ids) .abort("This session's field-overlays value holds the user field '", id, "' twice.")
    ids <- c(ids, id)
    source <- .session_source(field$source, paste0(label, ".source"))
    key <- .session_field_key(field$key, paste0(label, ".key"))
    deleted <- .session_bool(field$deleted, paste0(label, ".deleted"))
    if (purged && !deleted) .session_refuse(label, "is purged and not deleted.")
    created <- field$createdAt
    if (!is.null(created) && (is.logical(created) || !.is_finite_number(created) || created < 0)) {
      .session_refuse(label, ".createdAt must be a finite non-negative number or null.")
    }
    if (!is.null(field$operation)) .session_text(.session_record(field$operation, paste0(label, ".operation"))$type, paste0(label, ".operation.type"))
    source_field <- if (is.null(field$sourceField)) NULL else .session_source_field(field$sourceField, paste0(label, ".sourceField"))
    metadata <- list(kind = kind, field_id = id, key = key, source = source, deleted = deleted, source_field = source_field)
    what <- paste0("The user field '", key, "' (", id, ")")
    unapplied <- function(note) list(kind = "user-field", name = key, values = NULL, cells = NA_integer_, note = note, metadata = metadata)
    if (purged) {
      out[[length(out) + 1L]] <- unapplied(.apply_purged_note)
      next
    }
    colors <- NULL
    if (identical(kind, "category")) {
      categories <- .session_category_list(field$categories, paste0(label, ".categories"))
      .session_list(field$sourcePages, paste0(label, ".sourcePages"))
      if (!.is_json_string(field$overlapStrategy) || !field$overlapStrategy %in% .session_overlap_strategies) {
        .session_refuse(label, ".overlapStrategy is not one this version reads.")
      }
      .session_optional_text(field$overlapLabel, paste0(label, ".overlapLabel"))
      if (!is.null(field$intersectionLabels)) .session_record(field$intersectionLabels, paste0(label, ".intersectionLabels"))
      .session_optional_text(field$uncoveredLabel, paste0(label, ".uncoveredLabel"))
      # The colors its categories open with in the viewer (a field made from
      # highlight pages: its pages' colors), or null; written with the view
      # in focus's recolors (.apply_shown_colors).
      palette <- field$palette
      if (!is.null(palette) && (!is.list(palette) || length(palette) != length(categories) ||
        !all(vapply(palette, function(color) .is_json_string(color) && grepl("^#[0-9a-f]{6}$", color), logical(1L))))) {
        .session_refuse(label, ".palette must be null or one #rrggbb color a category.")
      }
      colors <- .apply_shown_colors(categories, if (is.null(palette)) NULL else vapply(palette, identity, ""),
                                    context$recolors[[paste0(source, "\u0001", id)]])
      part <- .session_column_part(field$codes, names(.session_code_types), n, paste0(label, ".codes"))
      missing <- .session_code_types[[part$type]]
      if (length(categories) > missing) .session_refuse(label, "has more categories than its code type (", part$type, ") can name.")
      invalid <- which(part$values >= length(categories) & part$values != missing)
      if (length(invalid) > 0L) {
        count <- length(categories)
        .session_refuse(label, ".codes holds the code ", format(part$values[[invalid[[1L]]]], scientific = FALSE),
                        ", and the field has ", count, " ", .plural(count, "category", "categories"), ".")
      }
      labels <- vapply(categories, .category_label, "")
      if (anyDuplicated(labels)) .abort(what, " has categories that read as the same label in R: ", .format_value_list(unique(labels[duplicated(labels)])), ".")
      if (deleted) {
        out[[length(out) + 1L]] <- unapplied(.apply_deleted_note)
        next
      }
      codes <- as.double(part$values)
      codes[codes == missing] <- NA
      values <- factor(labels[.placed_values(codes, part$cells, n) + 1], levels = labels)
    } else if (copy_of) {
      if (deleted) {
        out[[length(out) + 1L]] <- unapplied(.apply_deleted_note)
        next
      }
      values <- .continuous_alias_values(context$view, context$mapping, source, source_field, what)
    } else {
      part <- .session_column_part(field$values, "Float32Array", n, paste0(label, ".values"))
      if (any(is.infinite(part$values))) .session_refuse(label, ".values holds an infinity; a field holds finite values or NaN.")
      if (deleted) {
        out[[length(out) + 1L]] <- unapplied(.apply_deleted_note)
        next
      }
      numbers <- part$values
      numbers[is.nan(numbers)] <- NA
      values <- .placed_values(numbers, part$cells, n)
    }
    values <- .without_unknown(values, context$known, n)
    out[[length(out) + 1L]] <- list(kind = "user-field", name = key, values = values, cells = as.integer(sum(!is.na(values))),
                                    note = NA_character_, metadata = metadata, colors = colors)
  }
  out
}

# A name as part of a column name: each run of characters other than ASCII
# letters and digits one underscore, none at either end; `id` when nothing
# is left.
.apply_column_part <- function(text, id) {
  part <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", text, perl = TRUE))
  if (nzchar(part)) part else id
}

.apply_unique_column <- function(name, taken) {
  candidate <- name
  number <- 2L
  while (candidate %in% taken) {
    candidate <- paste0(name, "_", number)
    number <- number + 1L
  }
  candidate
}

# One candidate TRUE/FALSE column per highlight group of every page, checked
# whole before any is built: cellucid_highlight__<page>__<group>, TRUE for its
# cells, NA for every cell the step does not know.
.apply_highlights <- function(context, value) {
  keys <- .session_domain_keys
  highlights <- .session_keys(value, keys$highlights, "highlights value")
  pages <- .session_list(highlights$pages, "highlights value pages")
  page_ids <- character()
  group_ids <- character()
  taken <- character()
  out <- list()
  for (page_index in seq_along(pages)) {
    page_label <- paste0("highlights value pages[", page_index - 1L, "]")
    page <- .session_keys(pages[[page_index]], keys$page, page_label)
    page_id <- .session_text(page$id, paste0(page_label, ".id"))
    if (page_id %in% page_ids) .abort("This session's highlights value holds the page '", page_id, "' twice.")
    page_ids <- c(page_ids, page_id)
    page_name <- .session_text(page$name, paste0(page_label, ".name"))
    color <- .session_text(page$color, paste0(page_label, ".color"))
    if (!grepl("^#[0-9A-Fa-f]{6}$", color)) .session_refuse(page_label, ".color must be a six-digit hex color.")
    groups <- .session_list(page$groups, paste0(page_label, ".groups"))
    for (group_index in seq_along(groups)) {
      group_label <- paste0(page_label, ".groups[", group_index - 1L, "]")
      group <- .session_record(groups[[group_index]], group_label)
      metadata <- .session_group_metadata(group, group_label)
      group_id <- .session_text(group$id, paste0(group_label, ".id"))
      if (group_id %in% group_ids) .abort("This session's highlights value holds the group '", group_id, "' twice.")
      group_ids <- c(group_ids, group_id)
      cells <- .session_cell_list(group$cells, context$n, paste0(group_label, ".cells"))
      name <- .apply_unique_column(paste0("cellucid_highlight__", .apply_column_part(page_name, page_id), "__",
                                          .apply_column_part(metadata$label, group_id)), taken)
      taken <- c(taken, name)
      values <- .membership_column(cells, context$known, context$n)
      out[[length(out) + 1L]] <- list(
        kind = "highlight", name = name, values = values, cells = as.integer(sum(values, na.rm = TRUE)), note = NA_character_,
        metadata = c(list(kind = "highlight", group_id = group_id, page_id = page_id, page_name = page_name), metadata)
      )
    }
  }
  active <- highlights$active
  if (!is.null(active) && !(.is_json_string(active) && active %in% page_ids)) {
    .abort("This session's highlights value shows the page ", .format_value(active), ", which it does not hold.")
  }
  out
}

# The gene sets of a gene-sets value: list(name, description, genes) each.
.apply_read_gene_sets <- function(value) {
  record <- .session_keys(value, .session_domain_keys$gene_sets, "gene-sets value")
  sets <- .session_list(record$sets, "gene-sets value sets")
  lapply(seq_along(sets), function(index) {
    label <- paste0("gene-sets value sets[", index - 1L, "]")
    entry <- .session_keys(sets[[index]], .session_domain_keys$gene_set, label)
    name <- .session_text(entry$name, paste0(label, ".name"))
    if (!.is_json_string(entry$description)) .session_refuse(label, ".description must be a text.")
    genes <- .session_list(entry$genes, paste0(label, ".genes"))
    genes <- vapply(seq_along(genes), function(at) .session_text(genes[[at]], paste0(label, ".genes[", at - 1L, "]")), "")
    list(name = name, description = entry$description, genes = genes)
  })
}

# The gene sets the state holds: a row each, and the entries to write
# (cellucid$gene_sets, as cellucid_pull() keeps them; origin "session").
.apply_gene_sets <- function(sets, data, session, overwrite) {
  current <- .object_record(data)$gene_sets %||% list()
  rows <- list()
  writes <- list()
  for (entry in sets) {
    name <- entry$name
    genes <- entry$genes
    written <- list(genes = genes, description = entry$description, origin = "session", pulled_from = basename(session$path), id = name)
    held <- current[[name]]
    if (is.null(held)) {
      writes[[name]] <- written
      rows[[length(rows) + 1L]] <- .apply_row("gene-set", "gene-sets", name, "write")
      next
    }
    same <- is.list(held) && identical(as.character(unlist(held$genes)), genes) && identical(held$description, entry$description)
    owned <- is.list(held) && identical(held$origin, "session")
    if (same) {
      rows[[length(rows) + 1L]] <- .apply_row("gene-set", "gene-sets", name, "unchanged")
    } else if (!owned) {
      rows[[length(rows) + 1L]] <- .apply_row("gene-set", "gene-sets", name, "not-applied", note = .apply_taken_set_note(name))
    } else if (overwrite) {
      writes[[name]] <- written
      rows[[length(rows) + 1L]] <- .apply_row("gene-set", "gene-sets", name, "replace", note = "Written by an earlier apply; replaced.")
    } else {
      rows[[length(rows) + 1L]] <- .apply_row("gene-set", "gene-sets", name, "not-applied",
                                              note = paste0("An earlier apply wrote the gene set \"", name, "\" with other genes; overwrite replaces it."))
    }
  }
  list(rows = rows, writes = writes)
}

# --- Gene sets as feature-metadata columns ----------------------------------------------

.apply_gene_set_prefix <- "cellucid_gene_set__"
.apply_gene_sets_record_only <- paste0(
  "The gene sets are kept in the cellucid record's gene_sets only: the data mapping reads its genes from another ",
  "matrix than var, so they have no var columns."
)
.apply_own_var_column_note <- function(name, column) {
  paste0("The gene set \"", name, "\" has no var column: the var column \"", column, "\" is the dataset's own.")
}
.apply_missing_genes_note <- function(name, missing) {
  count <- length(missing)
  paste0("The gene set \"", name, "\" names ", count, " ", .plural(count, "gene"), " the dataset does not have (",
         paste(utils::head(missing, 3L), collapse = ", "), if (count > 3L) ", \u2026", "), left out of its var column.")
}

# Each gene set's column, by gene-set name: cellucid_gene_set__<name> (the
# highlight columns' rule; "gene_set" when nothing of the name is left), _2,
# _3 on clashes, given to every set of the step in its order, so a set's
# column is the same whichever of them an apply writes (the web app's
# geneSetColumns).
.apply_gene_set_column_names <- function(sets) {
  taken <- character()
  out <- character()
  for (entry in sets) {
    column <- .apply_unique_column(paste0(.apply_gene_set_prefix, .apply_column_part(entry$name, "gene_set")), taken)
    taken <- c(taken, column)
    out[[entry$name]] <- column
  }
  out
}

# The feature metadata of the gene axis the session's data mapping reads, as
# the part of the object that holds it: `where` is "none" (a data frame of
# cells has no feature axis), "elsewhere" (the mapping reads no genes, the raw
# matrix's, or a matrix cellucid_view() was given), "seurat" (the assay of
# that name: meta.features, an Assay5's meta.data), "main" (rowData of the
# experiment) or "altexp" (rowData of altExp(x, name)). The modality is the
# view's (.session_modality): a Seurat assay, or the main experiment (its
# mainExpName, else "main") or an altExp of a SingleCellExperiment, by name;
# one feature space is the only one. `genes` names each of the frame's rows
# as the mapping names genes (NA for a row the view does not show);
# `columns` are the frame's own columns.
.apply_feature_frame <- function(data, view, mapping, extra) {
  kind <- .object_kind(data)
  if (identical(kind, "data.frame")) return(list(where = "none"))
  expression <- mapping$expression
  if (is.null(expression) || identical(expression, "raw") || !is.null(extra$gene_expression)) return(list(where = "elsewhere"))
  modality <- .session_modality(view, mapping)
  name <- modality$name %||% modality$source
  genes <- .session_all_gene_names(view, mapping)
  if (identical(kind, "seurat")) {
    assay <- data[[name]]
    # The view's rows are the assay's matrix rows; the frame's, every feature of the assay.
    rows <- rownames(assay)
    return(list(where = "seurat", name = name, columns = names(assay[[]]),
                genes = genes[match(rows, rownames(.modality_matrix(modality)))]))
  }
  alt <- methods::is(data, "SingleCellExperiment") && name %in% SingleCellExperiment::altExpNames(data)
  experiment <- if (alt) SingleCellExperiment::altExp(data, name) else data
  # An experiment's view holds its rows in their order, one each.
  list(where = if (alt) "altexp" else "main", name = name, columns = names(SummarizedExperiment::rowData(experiment)),
       genes = genes)
}

# The column of each gene set the plan writes or replaces, and why any has
# none (warned when the plan is written).
.apply_gene_set_columns <- function(sets, gene_sets, frame) {
  out <- list(columns = list(), notes = character())
  if (length(gene_sets$writes) == 0L || identical(frame$where, "none")) return(out)
  if (identical(frame$where, "elsewhere")) {
    out$notes <- .apply_gene_sets_record_only
    return(out)
  }
  column_of <- .apply_gene_set_column_names(sets)
  status_of <- vapply(gene_sets$rows, function(row) row$status, "")
  names(status_of) <- vapply(gene_sets$rows, function(row) row$name, "")
  for (name in names(gene_sets$writes)) {
    column <- column_of[[name]]
    if (column %in% frame$columns && !identical(status_of[[name]], "replace")) {
      out$notes <- c(out$notes, .apply_own_var_column_note(name, column))
      next
    }
    genes <- gene_sets$writes[[name]]$genes
    missing <- genes[!genes %in% frame$genes]
    if (length(missing) > 0L) out$notes <- c(out$notes, .apply_missing_genes_note(name, missing))
    out$columns[[length(out$columns) + 1L]] <- list(name = column, values = !is.na(frame$genes) & frame$genes %in% genes)
  }
  out
}

# The object with each gene-set column written to the feature metadata of
# `frame`: a column of its name replaced in place, a new one after the
# frame's columns.
.object_with_feature_columns <- function(data, frame, columns) {
  if (length(columns) == 0L) return(data)
  write <- function(table) {
    for (column in columns) table[[column$name]] <- column$values
    table
  }
  switch(frame$where,
    seurat = {
      assay <- data[[frame$name]]
      rows <- rownames(assay)
      # SeuratObject's own accessor, for an Assay and an Assay5 alike: its
      # values are placed by feature name.
      for (column in columns) assay[[column$name]] <- stats::setNames(column$values, rows)
      data[[frame$name]] <- assay
      data
    },
    main = {
      SummarizedExperiment::rowData(data) <- write(SummarizedExperiment::rowData(data))
      data
    },
    altexp = {
      experiment <- SingleCellExperiment::altExp(data, frame$name)
      SummarizedExperiment::rowData(experiment) <- write(SummarizedExperiment::rowData(experiment))
      SingleCellExperiment::altExp(data, frame$name) <- experiment
      data
    }
  )
}

.apply_original_position_message <- paste0(
  "obs['", .original_position_column, "'] holds each row's original position in a file stored in priority order, ",
  "so an apply never writes, renames or replaces it."
)

# --- The plan --------------------------------------------------------------------------

.apply_rows_frame <- function(rows) {
  data.frame(
    kind = vapply(rows, function(row) row$kind, ""),
    domain = vapply(rows, function(row) row$domain, ""),
    name = vapply(rows, function(row) row$name, ""),
    column = vapply(rows, function(row) row$column, ""),
    from = I(lapply(rows, function(row) as.character(row$from))),
    cells = vapply(rows, function(row) as.integer(row$cells), 0L),
    status = vapply(rows, function(row) row$status, ""),
    note = vapply(rows, function(row) row$note, ""),
    stringsAsFactors = FALSE
  )
}

.apply_plan <- function(data, session, extra, tracker, branch, step, overwrite, keep_session, shown_as) {
  .require_session(session)
  if (identical(session$scope$content, "steps")) .abort(.apply_steps_only)
  chosen <- .apply_choose(session, tracker, branch, step)
  # An object read from a file in priority order (R/priority.R) holds its
  # cells in the stored order; a session names cells in their original order.
  stored <- .stored_order_of(data)
  if (!is.null(stored$notice)) .warn(stored$notice)
  if (!is.null(stored$positions)) {
    .abort(
      "This object holds its cells in a Cellucid-ready file's stored order (uns['cellucid']['priority']), ",
      "and a session names cells in their original order. ",
      "Put them back first: cellucid_restore_original_order(object)."
    )
  }
  n <- length(.object_cells(data))
  saved_on <- session$cells
  open_cells <- n
  kept <- NULL
  if (!is.null(saved_on)) {
    if (n != saved_on$spec$nObs) {
      .abort("This session was saved on ", saved_on$spec$k, " of ", saved_on$spec$nObs, " ", .plural(saved_on$spec$nObs, "cell"),
             ", and this object has ", n, ". Apply it to the dataset it was saved on.")
    }
    kept <- as.double(saved_on$kept)
    open_cells <- length(kept)
  }
  fingerprint <- session$dataset_fingerprint
  # The id the object was shown under (a viewer given dataset_id = serves it
  # under that id, and saves its session under it), else the one it declares.
  own_id <- .object_record(data)$id
  if (!is.null(fingerprint$datasetId) && !is.null(shown_as) && !identical(shown_as, fingerprint$datasetId)) {
    .abort("This session was saved on the dataset '", fingerprint$datasetId, "', and this object was shown as '", shown_as,
           "'. Apply it to the dataset it was saved on.")
  }
  if (is.null(shown_as) && .is_text(own_id) && !is.null(fingerprint$datasetId) && !identical(own_id, fingerprint$datasetId)) {
    .abort("This session was saved on the dataset '", fingerprint$datasetId, "', and this object is '", own_id,
           "' (its cellucid record's id). Apply it to the dataset it was saved on; when it is that dataset, shown under ",
           "another id (dataset_id = when it was shown), name that id: dataset_id = '", fingerprint$datasetId, "'.")
  }
  if (fingerprint$cellCount != open_cells) {
    .abort("This session was saved on ", fingerprint$cellCount, " ", .plural(fingerprint$cellCount, "cell"), ", and this object has ",
           open_cells, ". Apply it to the dataset it was saved on.")
  }
  view <- .session_view(data, extra)
  mapping <- fingerprint$dataMapping
  tryCatch(.require_mapping_modality(mapping, view$modality_names), cellucid_data_mapping_error = function(error) .abort(conditionMessage(error)))
  .require_mapping_matches_view(view, mapping)
  genes <- .session_var_count(view, mapping)
  if (fingerprint$varCount != genes) {
    .abort("This session was saved on a dataset of ", fingerprint$varCount, " ", .plural(fingerprint$varCount, "gene"),
           ", and this object has ", genes, " ", if (is.null(mapping$expression)) "genes on a mapping without expression" else "genes",
           ". Apply it to the dataset it was saved on.")
  }
  .apply_check_cell_order(fingerprint, view, kept)
  tracker_record <- chosen$tracker
  state <- .history_raw_state(tracker_record, chosen$step$id)
  hashes <- .history_hashes(tracker_record, chosen$step$id)
  root <- .history_hashes(tracker_record, tracker_record$root_id)
  for (domain in .apply_data_domains) {
    if (is.null(state[[domain]])) .abort("This session's state at step ", chosen$step$id, " holds no ", domain, " value.")
  }
  requires <- .apply_required_fields(state)
  .apply_require_fields(view, requires, mapping, chosen$step$id)

  applied <- cellucid_history_applied(data)
  owners <- list()
  for (record in applied) for (column in unlist(record$columns)) owners[[column]] <- record$step
  earlier <- unlist(lapply(applied, function(record) record$renamed), recursive = FALSE)

  known <- .session_known_cells(saved_on, state$cells, view, n, chosen$step$id)
  overlays <- .session_keys(state[["field-overlays"]], .session_domain_keys$overlays, "field-overlays value")
  modality_columns <- setdiff(names(view$obs), .object_column_names(data))
  labels <- .apply_plan_labels(overlays, data, earlier, modality_columns)
  present_after <- c(setdiff(.object_column_names(data), names(labels$rename)), unname(labels$rename))
  table <- .object_cell_table(data)
  columns <- list()
  settle <- function(candidate, domain) {
    name <- candidate$name
    if (is.null(candidate$values)) return(.apply_row(candidate$kind, domain, name, "not-applied", cells = candidate$cells, note = candidate$note))
    if (!name %in% present_after) {
      columns[[length(columns) + 1L]] <<- list(name = name, values = candidate$values, metadata = candidate$metadata,
                                               colors = candidate$colors)
      return(.apply_row(candidate$kind, domain, name, "write", cells = candidate$cells))
    }
    existing_name <- names(labels$rename)[labels$rename == name]
    existing <- table[[if (length(existing_name) > 0L) existing_name[[1L]] else name]]
    if (is.null(owners[[name]])) return(.apply_row(candidate$kind, domain, name, "not-applied", cells = candidate$cells, note = .apply_taken_note(name)))
    if (.apply_same_column(existing, candidate$values)) return(.apply_row(candidate$kind, domain, name, "unchanged", cells = candidate$cells))
    if (overwrite) {
      columns[[length(columns) + 1L]] <<- list(name = name, values = candidate$values, metadata = candidate$metadata,
                                               colors = candidate$colors)
      return(.apply_row(candidate$kind, domain, name, "replace", cells = candidate$cells,
                        note = paste0("Written by an earlier apply (step ", owners[[name]], "); replaced.")))
    }
    .apply_row(candidate$kind, domain, name, "not-applied", cells = candidate$cells, note = .apply_owned_note(name, owners[[name]]))
  }
  recolors <- .apply_focused_recolors(state[["colors"]], state[["views"]])
  context <- list(view = view, mapping = mapping, n = n, known = known, recolors = recolors)
  fields <- .apply_user_fields(context, overlays$userFields)
  field_rows <- lapply(fields, settle, domain = "field-overlays")
  # Whether a field written or replaced took colors the step's colors value gave it.
  recolored_written <- any(vapply(seq_along(fields), function(at) {
    field <- fields[[at]]$metadata
    field_rows[[at]]$status %in% c("write", "replace") && length(recolors[[paste0(field$source, "\u0001", field$field_id)]]) > 0L
  }, logical(1L)))
  sets <- .apply_read_gene_sets(state[["gene-sets"]])
  gene_sets <- .apply_gene_sets(sets, data, session, overwrite)
  feature_frame <- if (length(gene_sets$writes) > 0L) .apply_feature_frame(data, view, mapping, extra) else list(where = "none")
  var_columns <- .apply_gene_set_columns(sets, gene_sets, feature_frame)
  highlight_rows <- lapply(.apply_highlights(context, state$highlights), settle, domain = "highlights")
  # The priority-order column is each row's place in a stored file: nothing
  # an apply plans may write, rename or replace it.
  touched <- c(vapply(columns, function(column) column$name, ""), names(labels$replaced), names(labels$rename), unname(labels$rename))
  if (.original_position_column %in% touched) .abort(.apply_original_position_message)
  column_colors <- .apply_column_colors(.apply_focused_column_recolors(state[["colors"]], state[["views"]]), data, labels,
                                        modality_columns)
  # Whether any colors went with a column: a field's, or a recolored dataset column's.
  colors_went <- recolored_written || any(vapply(column_colors$rows, function(row) row$status != "not-applied", logical(1L)))
  by_domain <- list(`field-overlays` = c(labels$rows, field_rows), `gene-sets` = gene_sets$rows, highlights = highlight_rows)
  rows <- list()
  for (domain in .history_domains) {
    presentation <- NULL
    if (!domain %in% .apply_data_domains && !identical(hashes[[domain]], root[[domain]])) {
      note <- if (identical(domain, "colors") && colors_went) .apply_colors_note else .apply_presentation_note
      presentation <- list(.apply_row("presentation", domain, .history_domain_names[[domain]], "not-applied", note = note))
    }
    # A dataset column's colors, then the colors row.
    if (identical(domain, "colors")) rows <- c(rows, column_colors$rows)
    rows <- c(rows, if (!is.null(by_domain[[domain]])) by_domain[[domain]] else presentation)
  }
  structure(
    list(
      session_sha256 = session$sha256,
      tracker = tracker_record$id,
      branch = chosen$branch,
      step = chosen$step$id,
      requires = requires,
      rows = .apply_rows_frame(rows)
    ),
    class = "cellucid_apply_plan",
    cellucid_context = list(
      data = data, session = session, extra = extra, tracker = tracker, branch = branch, step = step, overwrite = overwrite,
      keep_session = keep_session, shown_as = shown_as, expression = .apply_expression_values(state[["expression"]], chosen$step$id),
      writes = list(columns = columns, replaced = labels$replaced, sources = labels$sources, rename = labels$rename, gene_sets = gene_sets$writes,
                    feature_frame = feature_frame, var_columns = var_columns$columns, var_notes = var_columns$notes,
                    renamed = labels$renamed, recolor = column_colors$recolor)
    )
  )
}

.apply_iso_now <- function() {
  now <- Sys.time()
  paste0(format(now, "%Y-%m-%dT%H:%M:%S", tz = "UTC"), ".", sprintf("%03d", as.integer((as.numeric(now) %% 1) * 1000)), "Z")
}

# Write a plan made just now: every column at once, then the object's record.
.apply_write <- function(plan) {
  context <- attr(plan, "cellucid_context", exact = TRUE)
  data <- context$data
  rows <- plan$rows
  writing <- rows$status %in% c("write", "replace")
  if (!any(writing)) {
    attr(data, "cellucid_apply_report") <- list(plan = plan, applied = FALSE, provenance = NULL)
    return(data)
  }
  writes <- context$writes
  provenance <- list(
    session_sha256 = plan$session_sha256,
    tracker = plan$tracker,
    branch = plan$branch,
    step = plan$step,
    applied_at = .apply_iso_now(),
    cellucid_version = .session_package_version(),
    package = "r",
    columns = as.list(rows$name[rows$kind %in% c("user-field", "highlight") & writing]),
    renamed = writes$renamed,
    skipped = lapply(which(rows$status == "not-applied"), function(at) list(kind = rows$kind[[at]], name = rows$name[[at]], note = rows$note[[at]])),
    expression = context$expression,
    session = if (isTRUE(context$keep_session)) {
      path <- context$session$path
      jsonlite::base64_enc(readBin(path, "raw", n = file.size(path)))
    } else NULL
  )
  json <- provenance
  json$renamed <- lapply(json$renamed, function(entry) {
    entry["column"] <- list(entry$column)
    entry
  })
  json["expression"] <- list(provenance$expression)
  json["session"] <- list(provenance$session)
  tryCatch(.history_check_record("Provenance", json, "provenance"), cellucid_history_model_error = function(error) .abort(conditionMessage(error)))
  record <- .object_record(data)
  if (length(writes$gene_sets) > 0L) {
    sets <- record$gene_sets %||% list()
    for (name in names(writes$gene_sets)) sets[[name]] <- writes$gene_sets[[name]]
    record$gene_sets <- sets
  }
  record$history_applied <- c(cellucid_history_applied(data), list(provenance))
  data <- .object_with_relabels(data, writes$replaced, writes$rename, writes$sources, writes$recolor)
  data <- .object_with_columns(data, writes$columns)
  data <- .object_with_feature_columns(data, writes$feature_frame, writes$var_columns)
  data <- .object_with_record(data, record)
  for (note in writes$var_notes) .warn(note)
  attr(data, "cellucid_apply_report") <- list(plan = plan, applied = TRUE, provenance = provenance)
  data
}

# The highlight groups and user-defined fields a session shows (its active
# tracker's head) as columns of `data`, and its gene sets: what
# cellucid_pull() writes from a session file, checked as an apply checks it.
.session_columns <- function(data, session, extra) {
  plan <- .apply_plan(data, session, extra, NULL, NULL, NULL, FALSE, FALSE, .apply_shown_as(extra))
  tracker <- cellucid_tracker(session, plan$tracker)
  state <- .history_raw_state(tracker, plan$step)
  n <- length(.object_cells(data))
  view <- .session_view(data, extra)
  known <- .session_known_cells(session$cells, state$cells, view, n, plan$step)
  overlays <- .session_keys(state[["field-overlays"]], .session_domain_keys$overlays, "field-overlays value")
  context <- list(view = view, mapping = session$dataset_fingerprint$dataMapping, n = n, known = known,
                  recolors = .apply_focused_recolors(state[["colors"]], state[["views"]]))
  columns <- c(.apply_highlights(context, state$highlights), .apply_user_fields(context, overlays$userFields))
  list(columns = Filter(function(column) !is.null(column$values), columns), gene_sets = .apply_read_gene_sets(state[["gene-sets"]]))
}

# The id the object was shown under, as cellucid_view() was given it
# (dataset_id = among the arguments after the session), or NULL.
.apply_shown_as <- function(extra) {
  if (is.null(extra$dataset_id)) NULL else .validate_dataset_id(extra$dataset_id)
}

.apply_session_of <- function(session, timeout) {
  if (inherits(session, "cellucid_viewer")) {
    .abort("Ask the viewer for its session first: request <- cellucid_request_session(viewer), then ",
           "cellucid_apply_session(object, request) (in a notebook, in a later cell).")
  }
  if (inherits(session, "cellucid_session_request")) return(cellucid_session_result(session, timeout = timeout))
  if (inherits(session, "cellucid_session")) return(session)
  cellucid_read_session(session)
}

cellucid_apply_session <- function(
    data,
    session,
    tracker = NULL,
    branch = NULL,
    step = NULL,
    ...,
    dry_run = TRUE,
    overwrite = FALSE,
    keep_session = FALSE,
    timeout = 60
) {
  if (missing(data) || missing(session)) {
    .abort("cellucid_apply_session() needs data (the object) and session (a .cellucid-session file, a session ",
           "cellucid_read_session() read, or a request).")
  }
  for (name in c("dry_run", "overwrite", "keep_session")) .validate_flag(get(name), name)
  for (name in c("tracker", "branch", "step")) {
    value <- get(name)
    if (!is.null(value) && !.is_json_string(value)) .abort(name, " must be one string (a name or an id) or NULL.")
  }
  extra <- .session_view_arguments(list(...))
  shown_as <- .apply_shown_as(extra)
  # A viewer's request: the viewer saved its session under the id it served
  # its data under, which dataset_id = may have made another than the one the
  # data declares; an object declaring that same id is the one it served.
  if (is.null(shown_as) && inherits(session, "cellucid_session_request")) {
    own_id <- .object_record(data)$id
    if (.is_text(own_id) && identical(own_id, session$viewer$declared_id)) shown_as <- session$viewer$dataset_id
  }
  plan <- .apply_plan(data, .apply_session_of(session, timeout), extra, tracker, branch, step, overwrite, keep_session, shown_as)
  if (dry_run) return(plan)
  .apply_write(plan)
}

cellucid_apply_plan <- function(plan) {
  if (!inherits(plan, "cellucid_apply_plan")) .abort("plan must be what cellucid_apply_session(..., dry_run = TRUE) returned.")
  context <- attr(plan, "cellucid_context", exact = TRUE)
  again <- .apply_plan(context$data, context$session, context$extra, context$tracker, context$branch, context$step,
                       context$overwrite, context$keep_session, context$shown_as)
  if (!identical(again$rows, plan$rows)) {
    .abort("The object changed since this plan was made; make the plan again (dry_run = TRUE) and apply that.")
  }
  .apply_write(again)
}

print.cellucid_apply_plan <- function(x, ...) {
  counts <- table(factor(x$rows$status, levels = .history_enums$APPLY_ROW_STATUSES))
  summary <- paste(paste(counts[counts > 0L], names(counts)[counts > 0L]), collapse = ", ")
  cat("<Cellucid apply plan> step ", x$step, " of branch ", x$branch, " of tracker ", x$tracker, ": ",
      if (nzchar(summary)) summary else "no change", "\n", sep = "")
  if (nrow(x$rows) > 0L) {
    shown <- x$rows
    shown$from <- vapply(shown$from, function(names) paste(names, collapse = ", "), "")
    print(shown[, c("kind", "name", "column", "from", "cells", "status", "note")], row.names = FALSE)
  }
  cat("Apply it with cellucid_apply_plan(plan).\n")
  invisible(x)
}
