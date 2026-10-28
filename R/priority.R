# Priority order: the row order of a Cellucid-ready file (docs: "Cellucid-ready
# AnnData specification", priority order; cellucid-python's
# prepare_data/_priority.py writes the same).
#
# Every cell has a priority, its `uniform` rank in the subsampling
# specification (R/subsample.R): rank(seed, p), a bijection of the uint32
# positions for one seed, so no two cells share one. A file in priority order
# stores its cells by ascending rank: row r holds the cell of the r-th lowest
# rank. The k cells a `uniform` subsample of size k keeps with the file's seed
# are then its first k rows in every array with one row per cell -- obs,
# obsm, a dense X, both axes of obsp -- and, in a gene-major X whose rows
# ascend within each gene, the first values of every gene: a subsample reads
# one range of each array and never a byte of a cell it left out.
#
# What the file records:
#   uns['cellucid']['priority'] = {seed, specVersion, nObs, levels}: the seed,
#     the subsampling specification's version the ranks follow, the cell count
#     the file was written with, and the prefix sizes (ascending,
#     .prefix_levels(nObs)) varm['cellucid_prefix_counts'] counts at;
#   obs['cellucid_original_position']: each row's position in the order the
#     cells were given (uint32), a permutation of 0 .. n - 1 (in a MuData
#     modality, its own rows' original rows);
#   varm['cellucid_prefix_counts'] beside a sparse gene-major X: n_vars x
#     length(levels) uint32, column j each gene's values in rows below
#     levels[j] (a MuData modality: below the global level, counted in the
#     modality's rows of those cells).
#
# The record is the claim, verified before anything relies on it
# (.stored_order(), the rule cellucid-python's stored_order() and the web
# app's readPriorityLayout() share): its nObs is the object's cell count, its
# levels those of that count, and the position column, when there, the
# seed's order. An object subset, combined or reordered after it was read
# keeps a record that no longer describes it: it is then taken as cells
# Cellucid did not write -- each row its own original position -- with a
# warning, never refused and never misread. A position column without a
# record names nothing.

.original_position_column <- "cellucid_original_position"
.prefix_counts_key <- "cellucid_prefix_counts"

# The version of the subsampling specification the ranks follow (as
# cellucid.subsample.SUBSAMPLE_SPEC_VERSION).
.subsample_spec_version <- 1L

# The prefix sizes counted, ascending: n stepped down by 3/4 (each level
# floor(previous * 3 / 4), exact in doubles below 2^51) to 1. n itself is not
# a level, and a level the rounding repeats is kept once. The factor is
# measured (cellucid-python's LEVEL_RATIO): a prefix read up to the next
# level fetches at most 1.01x its exact bytes, for 0.09-0.27% of X.
.prefix_levels <- function(n) {
  levels <- numeric(0)
  level <- as.double(n)
  repeat {
    level <- floor(level * 3 / 4)
    if (level < 1 || (length(levels) > 0L && level == levels[[length(levels)]])) break
    levels <- c(levels, level)
  }
  rev(levels)
}

# A uint32 column (written as "<u4", .ad_column()).
.uint32_column <- function(values) structure(as.double(values), class = "cellucid_uint32")

# The priority order of `n` cells for `seed`, or NULL to store them as given.
# The cells are in their original order (an object in a verified stored order
# was put back first, .in_original_order()), so each row is its own original
# position. `order[r]` is the given row (1-based) row r holds; `original[r]`
# its original position (0-based).
.cell_priority <- function(n, seed) {
  if (is.null(seed)) return(NULL)
  original <- seq_len(n) - 1
  # Ranks never tie for one seed: any sort gives the one order.
  order <- order(.subsample_rank(seed, original), method = "radix")
  new_row <- integer(n)
  new_row[order] <- seq_len(n)
  list(seed = seed, order = order, original = original[order], new_row = new_row, levels = .prefix_levels(n))
}

# uns['cellucid']['priority'].
.priority_record <- function(priority) {
  list(
    seed = .uns_integer(priority$seed),
    specVersion = .uns_integer(.subsample_spec_version),
    nObs = .uns_integer(length(priority$order)),
    levels = .uns_integers(priority$levels)
  )
}

# obs in file order, the original positions last.
.obs_in_priority <- function(obs, priority) {
  obs <- obs[priority$order, , drop = FALSE]
  obs[[.original_position_column]] <- NULL
  obs[[.original_position_column]] <- .uint32_column(priority$original)
  rownames(obs) <- NULL
  obs
}

# One feature space (an AnnData's, or a MuData modality's) with its rows in
# the order of their cells' priority: every cell (`positions` NULL) in file
# order, or a modality's cells ordered by their new global rows. Its prefix
# levels are its rows of cells below each global level, and `original` each
# row's original row in the space.
.space_in_priority <- function(space, priority) {
  n <- length(priority$order)
  if (is.null(space$positions)) {
    rows <- priority$order
    placed <- NULL
  } else {
    global <- priority$new_row[space$positions]
    rows <- order(global, method = "radix")
    placed <- global[rows]
  }
  # The expression matrix is read in this order, never copied in it
  # (expression-input.R): a served store holds the object and the order.
  if (!is.null(space$expression)) space$expression <- .expression_in_cell_order(space$expression, rows)
  if (!is.null(space$obs_index)) space$obs_index <- space$obs_index[rows]
  space$positions <- if (!is.null(placed) && !(length(placed) == n && identical(placed, seq_len(n)))) placed
  space$original <- rows - 1
  space$prefix_levels <- if (is.null(placed)) priority$levels else findInterval(priority$levels, placed)
  space
}

# The neighbour graph with both axes in priority order, without a reordered
# copy: the graph and the order are held, and its row-compressed arrays are
# made chunk by chunk as they are written or read (.ad_graph()).
.graph_in_priority <- function(graph, priority) {
  structure(list(matrix = graph, order = priority$order, new_row = priority$new_row), class = "cellucid_graph_in_order")
}

# varm['cellucid_prefix_counts'] of a feature space's sparse X: a
# n_genes x length(levels) matrix of counts, NULL for a dense X.
.expression_prefix_counts <- function(expression, levels) {
  if (!expression$sparse) return(NULL)
  n_genes <- .expression_dims(expression)[[2L]]
  matrix <- expression$matrix
  cells <- if (is.null(expression$cells)) integer(0) else expression$cells - 1L
  if (expression$feature_major) {
    .native_prefix_counts(matrix, levels, n_genes, TRUE, .expression_gene_map(expression)$map, cells = cells)
  } else {
    new_row <- if (is.null(expression$cells)) integer(0) else .expression_new_rows(expression)
    .native_prefix_counts(matrix, levels, n_genes, FALSE, new_row = new_row)
  }
}

# One whole number of a priority record (a number, a one-element vector, or a
# leaf as cellucid writes it), or NULL.
.record_whole_number <- function(value) {
  if (inherits(value, "cellucid_uns_leaf")) value <- unclass(value)$value
  value <- unclass(value)
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value != floor(value)) return(NULL)
  as.double(value)
}

# The one rule for a priority record (cellucid-python's stored_order(); the
# web app's readPriorityLayout()): are `n` rows in the order `record`
# (uns['cellucid']['priority'], NULL when absent) claims? It verifies when its
# seed is a uint32, its specVersion the one this package follows, its nObs is
# `n`, its levels (when present) are .prefix_levels(n), and `column`
# (obs$cellucid_original_position, NULL when absent) is exactly the seed's
# order. Then `positions` is that order (0-based original position of each
# row), computed from the seed. Anything else is data Cellucid did not write:
# `positions` NULL (each row its own original position), `notice` says why.
.stored_order <- function(record, n, column) {
  unused <- function(reason) {
    list(positions = NULL, notice = paste0(
      "uns['cellucid']['priority'] does not describe these cells: ", reason, ". They are taken in the ",
      "order they are in, as cells Cellucid did not write (each row its own original position); the web ",
      "app reads such a file the same way."
    ))
  }
  if (is.null(record)) {
    if (is.null(column)) return(list(positions = NULL, notice = NULL))
    return(list(positions = NULL, notice = paste0(
      .original_position_column, " has no uns['cellucid']['priority'] beside it (misc$cellucid$priority of a ",
      "Seurat object, metadata(x)$cellucid$priority of a SingleCellExperiment, the \"cellucid\" attribute of ",
      "a data frame, where anndataR and zellkonverter keep uns), so it names no order: the cells are taken in ",
      "the order they are in, as the web app reads such a file."
    )))
  }
  if (!is.list(record)) return(unused("it is not a list"))
  seed <- .record_whole_number(record$seed)
  if (is.null(seed) || seed < 0 || seed > 4294967295) return(unused("its seed is not a whole number from 0 to 4294967295"))
  version <- .record_whole_number(record$specVersion)
  if (!identical(version, as.double(.subsample_spec_version))) {
    return(unused(paste0(
      "it follows version ", if (is.null(version)) "unknown" else version,
      " of the subsampling specification, and this package follows ", .subsample_spec_version
    )))
  }
  recorded <- .record_whole_number(record$nObs)
  if (is.null(recorded)) return(unused("it records no cell count (nObs)"))
  if (recorded != n) {
    return(unused(paste0(
      "it was written for ", recorded, " ", .plural(recorded, "cell"), ", and there ", .plural(n, "is", "are"), " ", n,
      " (cells were added or removed since)"
    )))
  }
  if (!is.null(record$levels)) {
    levels <- record$levels
    if (inherits(levels, "cellucid_uns_leaf")) levels <- unclass(levels)$value
    levels <- unclass(levels)
    if (!is.numeric(levels) || !identical(as.double(levels), as.double(.prefix_levels(n)))) {
      return(unused(paste0("its levels are not those of ", n, " ", .plural(n, "cell"))))
    }
  }
  positions <- order(.subsample_rank(seed, seq_len(n) - 1), method = "radix") - 1
  if (!is.null(column)) {
    stored <- suppressWarnings(as.double(unclass(column)))
    if (!is.numeric(unclass(column)) || length(stored) != n || anyNA(stored) || !identical(stored, as.double(positions))) {
      return(unused(paste0(.original_position_column, " is not its seed's order (the cells were reordered since)")))
    }
  }
  list(positions = positions, notice = NULL)
}

# .stored_order() of a data frame (an AnnData's obs), a Seurat object or a
# SingleCellExperiment / SummarizedExperiment: its record where cellucid keeps
# one on the object (.object_record()), its column from its cell metadata.
.stored_order_of <- function(x) {
  frame <- switch(.object_kind(x),
    seurat = methods::slot(x, "meta.data"),
    se = SummarizedExperiment::colData(x),
    data.frame = x
  )
  column <- if (.original_position_column %in% names(frame)) frame[[.original_position_column]]
  .stored_order(.object_record(x)$priority, nrow(frame), column)
}

# A Seurat object with its cells in `order` (indices of its cells): every
# assay's layers and cell map (Assay5) or matrices (Assay), every reduction's
# cell embeddings, every graph (both axes), every neighbour index (rows, and
# the cells they point at), the cell metadata and the identities. Subsetting
# (x[, cells]) keeps a Seurat object's own order, so it cannot do this.
# Images keep their cells by name, as an image's cells need not be the
# object's in order.
.seurat_in_cell_order <- function(x, order) {
  cells <- colnames(x)
  ordered <- cells[order]
  for (name in names(x@assays)) {
    assay <- x@assays[[name]]
    if (methods::is(assay, "Assay5")) {
      map <- unclass(assay@cells)
      for (layer in names(assay@layers)) {
        held <- rownames(map)[map[, layer]]
        assay@layers[[layer]] <- assay@layers[[layer]][, match(ordered[ordered %in% held], held), drop = FALSE]
      }
      assay@cells <- methods::new("LogMap", map[ordered, , drop = FALSE])
    } else {
      for (slot in c("counts", "data", "scale.data")) {
        values <- methods::slot(assay, slot)
        if (ncol(values) > 0L) methods::slot(assay, slot) <- values[, ordered[ordered %in% colnames(values)], drop = FALSE]
      }
    }
    x@assays[[name]] <- assay
  }
  for (name in names(x@reductions)) {
    reduction <- x@reductions[[name]]
    embeddings <- reduction@cell.embeddings
    reduction@cell.embeddings <- embeddings[ordered[ordered %in% rownames(embeddings)], , drop = FALSE]
    x@reductions[[name]] <- reduction
  }
  for (name in names(x@graphs)) {
    graph <- x@graphs[[name]]
    reordered <- methods::as(graph, "dgCMatrix")[ordered, ordered, drop = FALSE]
    x@graphs[[name]] <- methods::new(class(graph), reordered, assay.used = graph@assay.used)
  }
  for (name in names(x@neighbors)) {
    neighbor <- x@neighbors[[name]]
    rows <- match(ordered, neighbor@cell.names)
    new_row <- integer(length(rows))
    new_row[rows] <- seq_along(rows)
    index <- neighbor@nn.idx[rows, , drop = FALSE]
    index[] <- new_row[index]
    neighbor@nn.idx <- index
    neighbor@nn.dist <- neighbor@nn.dist[rows, , drop = FALSE]
    neighbor@cell.names <- ordered
    x@neighbors[[name]] <- neighbor
  }
  x@meta.data <- x@meta.data[ordered, , drop = FALSE]
  x@active.ident <- x@active.ident[ordered]
  methods::validObject(x)
  x
}

# The cells of a Cellucid-ready file back in the order they were given
# (man/cellucid_restore_original_order.Rd), as cellucid-python's
# restore_original_order() puts an AnnData back: by the rule .stored_order()
# states; an object whose record does not describe it is returned as it is,
# with a warning.
cellucid_restore_original_order <- function(x) {
  if (!(is.data.frame(x) || inherits(x, "Seurat") || methods::is(x, "SummarizedExperiment"))) {
    .abort(
      "x must be a data.frame (an AnnData's obs), a Seurat object or a SingleCellExperiment / ",
      "SummarizedExperiment read from a Cellucid-ready file; got ", class(x)[[1L]], "."
    )
  }
  order <- .stored_order_of(x)
  if (!is.null(order$notice)) .warn(order$notice)
  if (is.null(order$positions)) return(x)
  .restored(x, order$positions)
}

# `x` in its original order for cellucid_prepare(): an object in a verified
# stored order is put back first, so a file prepared from it -- in priority
# order or not -- names every cell as the first file did.
.in_original_order <- function(x) {
  if (!(is.data.frame(x) || inherits(x, "Seurat") || methods::is(x, "SummarizedExperiment"))) return(x)
  order <- .stored_order_of(x)
  if (!is.null(order$notice)) .warn(order$notice)
  if (is.null(order$positions)) x else .restored(x, order$positions)
}

# `x`, its rows in a verified priority order (`positions` each row's 0-based
# original position), in the original order, without the position column and
# the priority record (which describe the stored order only).
.restored <- function(x, positions) {
  rows <- order(positions, method = "radix")
  column <- .original_position_column
  record <- .object_record(x)
  record$priority <- NULL
  if (length(record) == 0L) record <- NULL
  if (is.data.frame(x)) {
    restored <- x[rows, , drop = FALSE]
    restored[[column]] <- NULL
  } else if (inherits(x, "Seurat")) {
    restored <- .seurat_in_cell_order(x, rows)
    restored[[column]] <- NULL
  } else {
    restored <- x[, rows]
    SummarizedExperiment::colData(restored)[[column]] <- NULL
  }
  .object_with_record(restored, record)
}

# A cell metadata frame without a position column cellucid_prepare() was
# given: an object in a verified stored order was put back first
# (.in_original_order()), so a column still there names no order, and written
# again it would name cells the file does not hold where it says.
.without_positions <- function(obs) {
  if (is.null(obs) || !.original_position_column %in% names(obs)) return(obs)
  .warn(
    "obs$", .original_position_column, " is not written: without a priority record that describes these ",
    "cells (uns['cellucid']['priority']) it names no order, and cellucid_prepare() records the cells' order itself."
  )
  obs[[.original_position_column]] <- NULL
  obs
}
