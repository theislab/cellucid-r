# Expression matrices held outside R's memory: BPCells (an IterableMatrix --
# a matrix directory, a 10x or AnnData HDF5 file, or operations on them) and
# DelayedArray backends (HDF5Array, TENxMatrix, ...). Such a matrix is read
# a block at a time, never whole, never as a dgCMatrix of the whole matrix,
# and never in R's memory beyond one block (src/backed.c):
#
#   1. counting, when it is planned (.expression_input()): each gene's stored
#      values, their extent (the dtype, a value that cannot be stored refused
#      at once), one pass;
#   2. placing, when X is written or served: each value at its gene's
#      cursor, its cell's output row beside it, then each gene sorted by row
#      (both on every core, the genes shared among the threads), its prefix
#      counts taken from its sorted rows. The
#      genes a pass places hold at most a share of the machine's memory
#      (.backed_group_values()), and the pass reads only their rows; a larger
#      matrix takes one pass per share. (Prefix counts asked for without a
#      placing take one pass of their own.)
#
# A written store (cellucid_prepare()) places one block of genes, cuts its
# chunks and lets it go (X's values are written, then its cells, so a matrix
# of more than one share is placed twice); a served store lays the whole
# matrix out once into memory-mapped scratch files, as it lays out a matrix
# held in memory (expression-input.R), but on R's thread, before the viewer
# opens: R's thread is the only one that may read these backends. Positions
# and counts are doubles, so a matrix may hold more than 2^31 - 1 values.
# The bytes are those a dgCMatrix of the same values gives
# (test-backed-matrix.R).

# Whether `x` is a matrix held outside R's memory this package reads: the
# package that reads it ("BPCells" or "DelayedArray"), loaded, or NULL.
.backed_package <- function(x) {
  package <- attr(class(x), "package")
  if (isS4(x) && !is.null(package) && !package %in% c("methods", "Matrix", ".GlobalEnv")) {
    .require_object_package(package, class(x)[[1L]])
  }
  if (inherits(x, "IterableMatrix")) return("BPCells")
  if (methods::is(x, "DelayedArray") && length(dim(x)) == 2L) return("DelayedArray")
  NULL
}

# Values of one block, read at a time: about 2^24 positions of the matrix
# (whole cells, or whole genes), a scheduling choice that keeps a block's
# dgCMatrix small (a dense backend's block is made dense first), not a limit
# on the matrix.
.backed_block_values <- 16777216

# Whether the backend keeps genes (rows) together rather than cells: a
# BPCells matrix in row storage order (transpose_storage_order(), or the
# transpose of a cells x genes matrix), a DelayedArray whose chunks hold
# whole rows. Its blocks are then whole genes, so no block reads the whole
# matrix to find its cells.
.backed_by_rows <- function(matrix) {
  if (inherits(matrix, "IterableMatrix")) {
    return(identical(getExportedValue("BPCells", "storage_order")(matrix), "row"))
  }
  chunks <- DelayedArray::chunkdim(matrix)
  !is.null(chunks) && chunks[[2L]] >= ncol(matrix) && chunks[[1L]] < nrow(matrix)
}

# The blocks of a genes x cells matrix held outside R's memory, in the
# order the backend stores them, over its rows `rows` (1-based, ascending;
# NULL: every row): list(n, read(k)), read(k) giving the k-th block as the
# compressed arrays the native passes read -- list(index, pointers, values,
# by_rows, rows, col0): its major lines its columns (cells; the first one
# col0, 0-based) or, `by_rows`, its rows; `rows` the matrix rows its rows
# are. A pass that needs some genes reads only their rows.
#
# BPCells writes each block into memory as its own compressed arrays
# (write_matrix_memory(), values as doubles), read as they are: its
# dgCMatrix conversion builds the matrix through Matrix::sparseMatrix(),
# which took two thirds of a pass (44M values: 2.2 s against 0.7 s).
.backed_blocks <- function(matrix, rows = NULL) {
  dims <- dim(matrix)
  by_rows <- .backed_by_rows(matrix)
  all_rows <- is.null(rows)
  if (all_rows) rows <- seq_len(dims[[1L]])
  along <- if (by_rows) length(rows) else dims[[2L]]
  across <- if (by_rows) dims[[2L]] else length(rows)
  per_block <- max(1, floor(.backed_block_values / max(1, across)))
  starts <- if (along == 0) numeric(0) else seq(1, along, by = per_block)
  bpcells <- inherits(matrix, "IterableMatrix")
  read <- function(k) {
    first <- starts[[k]]
    index <- first:min(along, first + per_block - 1)
    block_rows <- if (by_rows) rows[index] else rows
    block <- if (by_rows) {
      matrix[block_rows, , drop = FALSE]
    } else if (all_rows) {
      matrix[, index, drop = FALSE]
    } else {
      matrix[block_rows, index, drop = FALSE]
    }
    col0 <- if (by_rows) 0 else first - 1
    if (bpcells) {
      held <- getExportedValue("BPCells", "write_matrix_memory")(
        getExportedValue("BPCells", "convert_matrix_type")(block, "double"), compress = FALSE
      )
      return(list(
        index = methods::slot(held, "index"), pointers = as.double(methods::slot(held, "idxptr")),
        values = methods::slot(held, "val"), by_rows = isTRUE(methods::slot(held, "transpose")),
        rows = block_rows, col0 = col0
      ))
    }
    block <- .without_omp_pause_warning(methods::as(block, "dgCMatrix"))
    list(index = block@i, pointers = as.double(block@p), values = block@x, by_rows = FALSE, rows = block_rows, col0 = col0)
  }
  list(n = length(starts), read = read)
}

# SparseArray's conversion of a block to a dgCMatrix asks OpenMP to let its
# threads go after each call, and an OpenMP runtime that cannot (macOS's
# libomp) warns "omp_pause_resource() failed to relinquish resources" every
# time: one warning per block, nothing wrong with the values. Only that
# warning is muffled, and only around that conversion.
.without_omp_pause_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(condition) {
    if (grepl("omp_pause_resource() failed to relinquish resources", conditionMessage(condition), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# The counting pass (.expression_input()): each row's (gene's) stored
# non-zero values (`count`), their smallest and largest (`low`, `high`) and
# whether all are whole numbers (`whole`), and the extent of them all as
# .native_extent() gives it. A value that is not finite or beyond the
# float32 range is refused, as one in a matrix held in memory is.
.backed_counts <- function(matrix, name) {
  .inform("Reading ", name, " (", class(matrix)[[1L]], ") block by block, never whole.")
  tally <- .native_backed_tally(nrow(matrix))
  blocks <- .backed_blocks(matrix)
  for (k in seq_len(blocks$n)) .native_backed_count(tally, blocks$read(k))
  counted <- .native_backed_counted(tally)
  summary <- counted[[5L]]
  non_finite <- c(nan = summary[[2L]], missing = summary[[3L]], positive = summary[[4L]], negative = summary[[5L]])
  if (sum(non_finite) > 0) {
    .abort(name, " must contain only finite values. ", .describe_non_finite_counts(non_finite, summary[[1L]], name))
  }
  stats <- list(count = counted[[1L]], low = counted[[2L]], high = counted[[3L]], whole = counted[[4L]])
  stats$extent <- .backed_extent(stats, seq_len(nrow(matrix)), zeros = summary[[6L]] > 0)
  if (any(stats$count > 0) && max(abs(stats$extent[1:2])) > .float32_max_magnitude) {
    .abort(
      name, " holds values beyond the float32 range (magnitude above ",
      "3.4e38), which would be stored as infinities."
    )
  }
  stats
}

# The extent (.native_extent()'s five numbers) of the stored values of rows
# `rows`.
.backed_extent <- function(stats, rows, zeros = FALSE) {
  held <- rows[stats$count[rows] > 0]
  if (length(held) == 0L) return(c(Inf, -Inf, 0, as.numeric(zeros), 1))
  c(min(stats$low[held]), max(stats$high[held]), 0, as.numeric(zeros), as.numeric(all(stats$whole[held])))
}

# Each column's (cell's) 0-based output row, -1 for a column left out: the
# columns chosen (.expression_select_cells()) taken in the cells' order
# (.expression_in_cell_order()).
.backed_new_rows <- function(expression) {
  n_columns <- ncol(expression$matrix)
  columns <- expression$columns %||% seq_len(n_columns)
  source <- if (is.null(expression$cells)) columns else columns[expression$cells]
  new_row <- rep.int(-1L, n_columns)
  new_row[source] <- seq_along(source) - 1L
  new_row
}

# Values one placing pass holds: a quarter of the machine's memory
# (.machine_memory_bytes()), so the rest stays for R, the backend's reads, the
# compressor and the system; a scheduling choice from the machine's real
# memory, not a limit on the matrix, which is placed in as many passes as it
# takes.
.backed_group_values <- function(bytes_per_value) {
  max(1, floor(.machine_memory_bytes() / 4 / bytes_per_value))
}

# The output genes (0-based first and last) of each placing pass: genes
# whose first position falls in the same share of `group_values` positions.
.backed_gene_groups <- function(pointers, group_values) {
  n_genes <- length(pointers) - 1L
  share <- floor(pointers[seq_len(n_genes)] / group_values)
  firsts <- which(!duplicated(share)) - 1
  Map(c, firsts, c(firsts[-1L] - 1, n_genes - 1))
}

# varm['cellucid_prefix_counts'] of a backed X: those its placing passes
# counted when X was written or laid out (every gene sorted by row; X is
# put before varm), else one pass of its own (src/backed.c).
.backed_prefix_counts <- function(expression, levels) {
  placed <- expression$placed$prefix_counts
  if (!is.null(placed) && identical(expression$placed$levels, as.double(levels)) && !anyNA(placed)) return(placed)
  genes <- .expression_gene_map(expression)
  tally <- .native_backed_prefix_tally(length(genes$origin), levels)
  new_row <- .backed_new_rows(expression)
  rows <- sort(genes$origin + 1L)
  blocks <- .backed_blocks(expression$matrix, rows)
  for (k in seq_len(blocks$n)) {
    read <- blocks$read(k)
    .native_backed_prefix_count(tally, read, genes$map[read$rows], new_row)
  }
  .native_backed_prefix_counted(tally)
}

# The gene-major csc_matrix arrays of a backed X (`shape` cells x genes):
# the pointers from the counts, and the indices and values placed. A written
# store places a block of positions (a whole number of its chunks) when its
# first chunk is written, and lets the previous one go; a served store lays
# the whole matrix out into scratch files now, one pass per share of memory
# (the module comment). Each pass reads only its genes' rows. The prefix
# counts of each gene sorted are kept for .backed_prefix_counts().
.backed_arrays <- function(store, path, expression, shape, dtype) {
  matrix <- expression$matrix
  genes <- .expression_gene_map(expression)
  pointers <- c(0, cumsum(expression$stats$count[genes$origin + 1L]))
  n_values <- pointers[[length(pointers)]]
  index_size <- .zarr_item_size(.sparse_index_dtype("csc_matrix", shape))
  sizes <- c(index_size, .zarr_item_size(dtype))
  if (n_values == 0) return(list(pointers = pointers, indices = integer(0), values = numeric(0)))
  n_cells <- shape[[1L]]
  new_row <- .backed_new_rows(expression)
  levels <- expression$prefix_levels
  placed <- expression$placed
  if (!is.null(placed) && !is.null(levels)) {
    placed$levels <- as.double(levels)
    placed$prefix_counts <- matrix(NA_real_, length(genes$origin), length(levels))
  }
  # The cores the placing and sorting share: those a written store
  # compresses on (.deflate_threads(); two under R CMD check).
  threads <- .deflate_threads()
  place <- function(target, first_gene, last_gene) {
    blocks <- .backed_blocks(matrix, sort(genes$origin[(first_gene:last_gene) + 1L] + 1L))
    for (k in seq_len(blocks$n)) {
      read <- blocks$read(k)
      .native_backed_place(target, first_gene, last_gene, read, genes$map[read$rows], new_row, threads)
    }
    counts <- .native_backed_finish(target, first_gene, last_gene, levels, threads)
    if (!is.null(placed$prefix_counts)) placed$prefix_counts[(first_gene:last_gene) + 1L, ] <- counts
    invisible(NULL)
  }
  group_values <- .backed_group_values(sum(sizes))
  part <- function(which, target_of) structure(list(
    length = n_values,
    chunk = function(first, count, chunk_length, ...) {
      target_of(first - 1, first - 1 + count, which, function(target) {
        .native_backed_slice(target, which, first - 1, count, chunk_length)
      })
    }
  ), class = "cellucid_chunked")
  if (isTRUE(store$memory)) {
    paths <- tempfile(c("cellucid-layout-indices-", "cellucid-layout-values-"))
    target <- .native_backed_target(pointers, 0, n_values, n_cells, index_size, dtype, paths)
    store$backed_layouts <- c(store$backed_layouts, list(target))
    for (group in .backed_gene_groups(pointers, group_values)) place(target, group[[1L]], group[[2L]])
    laid_out <- function(first, end, which, cut) cut(target)
    return(list(pointers = pointers, indices = part(1L, laid_out), values = part(2L, laid_out)))
  }
  chunk <- .sliced_chunk_elements(store)
  block_values <- max(1, floor(group_values / chunk)) * chunk
  # The block held, let go once the last chunk of both parts is cut (R's
  # collector does not see memory held outside it).
  held <- new.env(parent = emptyenv())
  held$ended <- logical(2)
  release <- function() {
    if (!is.null(held$target)) .native_backed_release(held$target)
    held$target <- NULL
    held$index <- NULL
  }
  block_of <- function(first, end, which, cut) {
    index <- first %/% block_values
    if (!identical(held$index, index)) {
      release()
      start <- index * block_values
      target <- .native_backed_target(pointers, start, min(n_values, start + block_values), n_cells, index_size, dtype)
      genes_held <- .native_backed_target_genes(target)
      place(target, genes_held[[1L]], genes_held[[2L]])
      held$target <- target
      held$index <- index
    }
    bytes <- cut(held$target)
    if (end >= n_values) held$ended[[which]] <- TRUE
    if (all(held$ended)) release()
    bytes
  }
  list(pointers = pointers, indices = part(1L, block_of), values = part(2L, block_of))
}

# Let a served store's backed layouts go (the viewer stopped).
.backed_release_layouts <- function(store) {
  for (target in store$backed_layouts) .native_backed_release(target)
  store$backed_layouts <- NULL
  invisible(store)
}
