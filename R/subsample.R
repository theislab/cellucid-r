# Subsampling: the cells a Cellucid subsample keeps, the counterpart of
# cellucid-python's subsample_cells() and of the web app's subsampled opening
# (Open data, Cells to open). The help page is hand-written in
# man/cellucid_subsample_cells.Rd, as every page of this package is.
#
# The arithmetic is the subsampling specification's (docs: "Subsampling
# specification"), done in src/subsample.c: uint32 hashes, exact 64-bit
# allocations, and the spatial grid in IEEE 754 binary64 operations never
# fused. tests/testthat/fixtures/subsample_conformance.json, shared byte for
# byte with the web app and cellucid-python, pins the result.

.subsample_methods <- c("uniform", "blocks", "stratified", "spatial")

# The cell count of what cellucid_subsample_cells() is given, whose columns
# (rows of a data frame) are each its own original position. An object in a
# Cellucid-ready file's verified stored order (.stored_order()) is refused, as
# cellucid_pull(), cellucid_subset() and a session refuse it: the cells are
# chosen by original position, and a column of such an object is not one, so
# x[, kept] would take other cells than the viewer's. A record that no longer
# describes the cells leaves each its own original position, with a warning.
.subsample_n_obs <- function(x) {
  if (is.numeric(x) && length(x) == 1L && !is.object(x)) {
    return(.subsample_uint32(x, "x"))
  }
  if (!(inherits(x, "Seurat") || methods::is(x, "SummarizedExperiment") || is.data.frame(x))) {
    .abort("x must be a Seurat object, a SingleCellExperiment or SummarizedExperiment, a data.frame of cells, or the cell count.")
  }
  stored <- .stored_order_of(x)
  if (!is.null(stored$notice)) .warn(stored$notice)
  if (!is.null(stored$positions)) {
    .abort(
      "This object holds its cells in a Cellucid-ready file's stored order (uns['cellucid']['priority']), ",
      "and a subsample names cells by their original position, so x[, kept] would take other cells than the ",
      "viewer opens. Put them back first: x <- cellucid_restore_original_order(x)."
    )
  }
  as.double(if (is.data.frame(x)) nrow(x) else ncol(x))
}

.subsample_uint32 <- function(value, label, minimum = 0) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value != floor(value) ||
      value < minimum || value > 4294967295) {
    .abort(label, " must be one whole number in [", minimum, ", 2^32 - 1].")
  }
  as.double(value)
}

# The cells' metadata column `name` (Seurat meta.data, colData, a data.frame).
.subsample_column <- function(x, name) {
  frame <- if (inherits(x, "Seurat")) {
    x@meta.data
  } else if (methods::is(x, "SummarizedExperiment")) {
    as.data.frame(SummarizedExperiment::colData(x))
  } else if (is.data.frame(x)) {
    x
  } else {
    .abort("groupby names a column, so x must be an object with cell metadata.")
  }
  if (!name %in% names(frame)) .abort("The cell metadata has no column '", name, "'.")
  frame[[name]]
}

# The embedding `name` (a Seurat reduction, a reducedDim).
.subsample_embedding <- function(x, name) {
  if (inherits(x, "Seurat")) {
    if (!name %in% names(x@reductions)) .abort("The Seurat object has no reduction '", name, "'.")
    return(x@reductions[[name]]@cell.embeddings)
  }
  if (methods::is(x, "SingleCellExperiment")) {
    if (!name %in% SingleCellExperiment::reducedDimNames(x)) .abort("The object has no reducedDim '", name, "'.")
    return(SingleCellExperiment::reducedDim(x, name))
  }
  .abort("basis names an embedding, so x must be a Seurat object or a SingleCellExperiment.")
}

# The subsample size a fraction of n_obs cells stands for, exactly: the
# fraction as the rational its shortest round-trip decimal names (as the web
# app's String(x) and Python's repr, a tie to the even digit; found in exact
# integer arithmetic in C, src/subsample.c, since neither R's own text
# conversion nor the C library's printf gives it on every build),
# round-half-up(n_obs * fraction) in integer arithmetic (decimal digits in
# doubles, every step exact), at least 1 and at most n_obs.
.subsample_count_from_fraction <- function(fraction, n_obs) {
  if (!is.numeric(fraction) || length(fraction) != 1L || !is.finite(fraction) || fraction <= 0 || fraction > 1) {
    .abort("fraction must be one number in (0, 1].")
  }
  text <- .native_shortest_decimal(fraction)
  parts <- regmatches(text, regexec("^(\\d+)(?:\\.(\\d*))?(?:[eE]([+-]?\\d+))?$", text))[[1]]
  decimals <- parts[3]
  mantissa <- sub("^0+(?=.)", "", paste0(parts[2], decimals), perl = TRUE)
  scale <- nchar(decimals) - (if (nzchar(parts[4])) as.integer(parts[4]) else 0L)
  if (scale <= 0) return(as.double(n_obs))
  # N = 2 * n_obs * mantissa + 10^scale, as decimal digits (least significant first).
  digits <- rev(as.numeric(strsplit(mantissa, "")[[1]])) * (2 * n_obs)
  width <- max(length(digits), scale + 1L) + 12L
  number <- c(digits, rep(0, width - length(digits)))
  number[scale + 1L] <- number[scale + 1L] + 1
  carry <- 0
  for (i in seq_along(number)) {
    total <- number[i] + carry
    number[i] <- total %% 10
    carry <- (total - number[i]) / 10
  }
  # floor(N / (2 * 10^scale)): halve, then drop the last `scale` digits.
  remainder <- 0
  for (i in rev(seq_along(number))) {
    current <- remainder * 10 + number[i]
    number[i] <- current %/% 2
    remainder <- current %% 2
  }
  kept <- number[-seq_len(scale)]
  k <- sum(kept * 10^(seq_along(kept) - 1L))
  as.double(min(n_obs, max(1, k)))
}

cellucid_subsample_cells <- function(x, n_cells = NULL, fraction = NULL, method = "uniform", seed = 0,
                                     groupby = NULL, allocation = "proportional", min_per_group = 0,
                                     block_size = NULL, basis = NULL, columns = c(1, 2), grid_size = 32) {
  n_obs <- .subsample_n_obs(x)
  seed <- .subsample_uint32(seed, "seed")
  if (!is.character(method) || length(method) != 1L || !method %in% .subsample_methods) {
    .abort("method must be one of ", paste0("'", .subsample_methods, "'", collapse = ", "), ".")
  }
  if (is.null(n_cells) == is.null(fraction)) .abort("Give exactly one of n_cells and fraction.")
  k <- if (is.null(fraction)) .subsample_uint32(n_cells, "n_cells", minimum = 1) else .subsample_count_from_fraction(fraction, n_obs)
  positions <- if (k >= n_obs) {
    seq_len(n_obs) - 1
  } else if (method == "uniform") {
    .native_subsample_uniform(n_obs, k, seed)
  } else if (method == "blocks") {
    if (is.null(block_size)) .abort("method 'blocks' needs block_size (the web app shows the one it used).")
    .native_subsample_blocks(n_obs, k, seed, .subsample_uint32(block_size, "block_size", minimum = 1))
  } else if (method == "stratified") {
    if (is.null(groupby)) .abort("method 'stratified' needs groupby.")
    if (!is.character(allocation) || length(allocation) != 1L || !allocation %in% c("proportional", "balanced")) {
      .abort("allocation must be 'proportional' or 'balanced'.")
    }
    values <- if (is.character(groupby) && length(groupby) == 1L && n_obs != 1) .subsample_column(x, groupby) else groupby
    if (length(values) != n_obs) .abort("groupby has ", length(values), " ", .plural(length(values), "value"), " for ", n_obs, " ", .plural(n_obs, "cell"), ".")
    # Categories in code-point order (a factor keeps its own), as Python's
    # pandas Categorical and the web app order them: never the locale's.
    groups <- if (is.factor(values)) values else factor(values, levels = sort(unique(values[!is.na(values)]), method = "radix"))
    n_categories <- nlevels(groups)
    strata <- as.double(as.integer(groups) - 1L)
    strata[is.na(strata)] <- n_categories
    .native_subsample_split(strata, n_categories + 1, k, seed, allocation,
                            .subsample_uint32(min_per_group, "min_per_group"))
  } else {
    if (is.null(basis)) .abort("method 'spatial' needs basis.")
    coordinates <- if (is.character(basis) && length(basis) == 1L) .subsample_embedding(x, basis) else basis
    coordinates <- as.matrix(coordinates)
    if (nrow(coordinates) != n_obs) .abort("basis has ", nrow(coordinates), " ", .plural(nrow(coordinates), "row"), " for ", n_obs, " ", .plural(n_obs, "cell"), ".")
    if (length(columns) < 1L || length(columns) > 3L) .abort("columns names one to three embedding columns.")
    if (!is.numeric(columns) || anyNA(columns) || any(columns != floor(columns)) || any(columns < 1) || any(columns > ncol(coordinates))) {
      .abort("columns must name embedding columns from 1 to ", ncol(coordinates), "; got c(", paste(columns, collapse = ", "), ").")
    }
    coordinates <- coordinates[, columns, drop = FALSE]
    storage.mode(coordinates) <- "double"
    if (as.double(grid_size)^ncol(coordinates) >= 2^53) {
      .abort("grid_size^", ncol(coordinates), " boxes are more than are numbered exactly (2^53).")
    }
    boxes <- .native_subsample_boxes(as.vector(t(coordinates)), ncol(coordinates),
                                     .subsample_uint32(grid_size, "grid_size", minimum = 1))
    occupied <- sort(unique(boxes))
    .native_subsample_split(as.double(match(boxes, occupied) - 1L), length(occupied), k, seed, "balanced", 0)
  }
  # R's own 1-based positions; doubles only past what an R integer counts.
  if (n_obs <= .Machine$integer.max) as.integer(positions + 1) else positions + 1
}

# The specification's cell-set digest of 0-based kept positions (NULL: every
# cell): what a session saved on a subsample records (manifest cells.digest).
.subsample_digest <- function(n_obs, kept) .native_subsample_digest(n_obs, kept)

# A cell set (0-based positions) as a session stores it, and back.
.subsample_encode_cell_set <- function(kept) .native_subsample_encode(kept)
.subsample_decode_cell_set <- function(bytes, n_obs) .native_subsample_decode(bytes, n_obs)

# The rank of 0-based positions under a seed (the conformance test reads it).
.subsample_rank <- function(seed, positions) .native_subsample_rank(seed, positions)

# The `cells` record of a session saved on a subsample, checked and decoded:
# list(spec, kept = 0-based positions ascending, digest). Every per-cell
# chunk of a session names its cells by original position, ascending: cell i
# of a subsample's session is obs row kept[i], 0-based.
.session_cells <- function(cells, cell_count) {
  if (!is.list(cells) || !setequal(names(cells), c("spec", "kept", "digest")) || anyDuplicated(names(cells))) {
    .session_manifest_error("cells must have exactly spec, kept and digest")
  }
  spec <- cells[["spec"]]
  if (!is.list(spec) || !is.character(spec[["method"]]) || !spec[["method"]] %in% .subsample_methods) {
    .session_manifest_error("cells.spec must name a subsample method")
  }
  # Exactly the keys of its method, as the web app's requireSubsampleSpec.
  method_keys <- list(uniform = character(), blocks = "blockSize", stratified = c("field", "allocation", "minPerCategory"),
                      spatial = c("embedding", "columns", "gridSize"))[[spec[["method"]]]]
  if (!setequal(names(spec), c("method", "nObs", "k", "seed", method_keys)) || anyDuplicated(names(spec))) {
    .session_manifest_error("cells.spec of ", spec[["method"]], " must have exactly method, nObs, k, seed",
                            paste0(", ", method_keys, collapse = ""))
  }
  n_obs <- spec[["nObs"]]
  k <- spec[["k"]]
  if (!is.numeric(n_obs) || !is.numeric(k) || length(n_obs) != 1L || length(k) != 1L ||
      !(k >= 1 && k < n_obs && n_obs <= 4294967295)) {
    .session_manifest_error("cells.spec must keep 1 <= k < nObs cells")
  }
  if (!identical(as.numeric(k), as.numeric(cell_count))) {
    .session_manifest_error("cells must keep exactly the fingerprint cell count")
  }
  if (!is.character(cells[["kept"]]) || !grepl("^[A-Za-z0-9+/]*={0,2}$", cells[["kept"]])) {
    .session_manifest_error("cells.kept must be base64")
  }
  kept <- .subsample_decode_cell_set(jsonlite::base64_dec(cells[["kept"]]), n_obs)
  if (length(kept) != k) .session_manifest_error("cells.kept must list exactly k cells")
  if (!identical(.subsample_digest(n_obs, kept), cells[["digest"]])) {
    .session_manifest_error("cells.digest does not match its cells")
  }
  if (identical(spec[["method"]], "blocks") &&
      !(is.numeric(spec[["blockSize"]]) && length(spec[["blockSize"]]) == 1L && spec[["blockSize"]] >= 1)) {
    .session_manifest_error("cells.spec must give a blocks subsample its blockSize")
  }
  list(spec = spec, kept = kept, digest = cells[["digest"]])
}

# The `genes` record of a session saved while genes were read for a sample
# of the open cells: genes were read for the first `cells` open cells in
# their order (the order of the `blocks` subsample of every cell `order`
# names -- every cell's, or a share's spread through a file in its own
# order, its cells filtered -- or, with `order` null, a subsample's priority
# order), whose original positions
# have the cell-set digest `digest`. Checked as the web app's
# decodeGeneSampleRecord and cellucid-python's _decode_genes check it.
.session_genes <- function(genes, cell_count) {
  if (!is.list(genes) || !setequal(names(genes), c("cells", "order", "digest")) || anyDuplicated(names(genes))) {
    .session_manifest_error("genes must have exactly cells, order and digest")
  }
  cells <- genes[["cells"]]
  if (!is.numeric(cells) || length(cells) != 1L || cells != floor(cells) || !(cells >= 1 && cells < cell_count)) {
    .session_manifest_error("genes.cells must be a count of open cells below the fingerprint cell count")
  }
  order <- genes[["order"]]
  if (!is.null(order)) {
    if (!is.list(order) || !setequal(names(order), c("method", "seed", "blockSize")) || anyDuplicated(names(order)) ||
        !identical(order[["method"]], "blocks")) {
      .session_manifest_error("genes.order must be {method: 'blocks', seed, blockSize}")
    }
    seed <- order[["seed"]]
    if (!is.numeric(seed) || length(seed) != 1L || seed != floor(seed) || seed < 0 || seed > 4294967295) {
      .session_manifest_error("genes.order.seed must be a uint32")
    }
    block_size <- order[["blockSize"]]
    if (!is.numeric(block_size) || length(block_size) != 1L || block_size != floor(block_size) || block_size < 1) {
      .session_manifest_error("genes.order.blockSize must be a positive count")
    }
  }
  digest <- genes[["digest"]]
  if (!is.character(digest) || length(digest) != 1L || !grepl("^[0-9a-f]{16}$", digest)) {
    .session_manifest_error("genes.digest must be 16 lowercase hex digits")
  }
  list(cells = cells, order = order, digest = digest)
}
