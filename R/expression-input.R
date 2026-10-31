# One expression matrix as the plan holds it: the matrix as the caller or the
# object gave it, and its orientation. A caller passes cells x genes (rows are
# cells); an R object holds genes x cells (a Seurat assay, a
# SummarizedExperiment assay), which is marked .feature_major() and never
# transposed whole: the store's metadata (shape, number of values, dtype) is
# the same either way, and the gene-major arrays of a sparse matrix -- its
# transpose's row indices and values -- are made block by block from the
# matrix as it is (.gene_major_arrays()). A dense matrix is read the same way,
# each chunk's block taken from it and transposed alone (.transposed() in
# zarr-store.R).

.feature_major <- function(matrix) {
  structure(list(matrix = matrix), class = "cellucid_feature_major")
}

# Checked and normalized: a dgCMatrix with no stored zeros, or a dense double
# matrix, every value finite and within the float32 range. A matrix held
# outside R's memory (BPCells, a DelayedArray backend; backed-matrix.R) stays
# where it is, genes x cells (a cells x genes one is read through its
# transpose, which both backends make without reading it), its values counted
# in one pass; it is written and served sparse, as its dgCMatrix would be.
.expression_input <- function(x, name) {
  feature_major <- inherits(x, "cellucid_feature_major")
  matrix <- if (feature_major) x$matrix else x
  if (!is.null(.backed_package(matrix))) {
    # Matrix's generic, on which both backends define their lazy transpose.
    if (!feature_major) matrix <- Matrix::t(matrix)
    return(list(
      matrix = matrix, stats = .backed_counts(matrix, name),
      feature_major = TRUE, sparse = TRUE, backed = TRUE
    ))
  }
  checked <- .as_expression_matrix(matrix, name)
  list(
    matrix = checked$matrix,
    extent = checked$extent,
    feature_major = feature_major,
    sparse = inherits(checked$matrix, "sparseMatrix")
  )
}

# c(cells, genes).
.expression_dims <- function(expression) {
  shape <- dim(expression$matrix)
  if (!is.null(expression$columns)) shape[[2L]] <- length(expression$columns)
  if (expression$feature_major) shape <- rev(shape)
  if (!is.null(expression$genes)) shape[[2L]] <- length(expression$genes)
  shape
}

.expression_cell_names <- function(expression) {
  names <- if (expression$feature_major) colnames(expression$matrix) else rownames(expression$matrix)
  if (is.null(expression$columns) || is.null(names)) names else names[expression$columns]
}

.expression_gene_names <- function(expression) {
  names <- if (expression$feature_major) rownames(expression$matrix) else colnames(expression$matrix)
  if (is.null(expression$genes) || is.null(names)) names else names[expression$genes]
}

# The genes kept, by position (gene_identifiers). A genes x cells sparse
# matrix (in memory or held outside it) keeps them as positions read by its
# gene-major arrays (.gene_major_arrays(), .backed_arrays()), never as a copy
# of its rows: at a million cells the row subset copied every value of the
# matrix. Every other matrix is narrowed to its columns (cells x genes) or
# rows (dense genes x cells).
.expression_select_genes <- function(expression, positions) {
  if (expression$feature_major && expression$sparse) {
    current <- if (is.null(expression$genes)) seq_len(nrow(expression$matrix)) else expression$genes
    expression$genes <- current[positions]
  } else if (expression$feature_major) {
    expression$matrix <- expression$matrix[positions, , drop = FALSE]
  } else {
    expression$matrix <- expression$matrix[, positions, drop = FALSE]
  }
  expression$extent <- NULL
  expression
}

# A genes x cells matrix's rows as the gene-major arrays read them: `map`,
# each row's 0-based output gene or -1 when it is left out, and `origin`,
# each output gene's 0-based row.
.expression_gene_map <- function(expression) {
  n_rows <- nrow(expression$matrix)
  kept <- if (is.null(expression$genes)) seq_len(n_rows) else expression$genes
  map <- rep.int(-1L, n_rows)
  map[kept] <- seq_along(kept) - 1L
  list(map = map, origin = as.integer(kept - 1L))
}

# The matrix with its cells (rows of a cells x genes matrix, columns of a
# feature-major one) in the order `rows` gives: a reordering keeps every
# value, so the measured extent still holds. A matrix held outside R's memory
# keeps the columns chosen (`columns`), read through them, never a reordered
# copy of the backend (whose blocks are read in the order it stores them).
.expression_select_cells <- function(expression, rows) {
  if (isTRUE(expression$backed)) {
    expression$columns <- (expression$columns %||% seq_len(ncol(expression$matrix)))[rows]
    return(expression)
  }
  expression$matrix <- if (expression$feature_major) {
    expression$matrix[, rows, drop = FALSE]
  } else {
    expression$matrix[rows, , drop = FALSE]
  }
  expression
}

# The matrix with its cells in the order `rows` gives (1-based), as a store
# in priority order presents it (R/priority.R), without a reordered copy: a
# sparse genes x cells matrix is read in that order (its gene-major blocks
# visit the cells so, src/chunks.c), a sparse cells x genes one renumbered
# (each gene's values sorted by new row when its block is made), a dense one
# read row block by row block through the order. Only the order is held.
.expression_in_cell_order <- function(expression, rows) {
  expression$cells <- as.integer(rows)
  expression
}

# The 0-based new row of each cell of a cells x genes matrix in cell order.
.expression_new_rows <- function(expression) {
  new_row <- integer(length(expression$cells))
  new_row[expression$cells] <- seq_along(expression$cells) - 1L
  new_row
}

# The stored values' dtype, from the extent already measured when no genes
# were left out since.
.expression_value_dtype <- function(expression) {
  if (isTRUE(expression$backed)) {
    # The kept genes' counted extent (.backed_counts()).
    rows <- if (is.null(expression$genes)) seq_len(nrow(expression$matrix)) else expression$genes
    extent <- .backed_extent(expression$stats, rows)
    if (!is.finite(extent[[1L]]) || max(abs(extent[1:2])) <= .float32_exact_integer_limit) return("<f4")
    return(if (extent[[5L]] == 1) "<f8" else "<f4")
  }
  if (!is.null(expression$genes)) {
    # The kept genes' values, read in place (.native_kept_extent()).
    map <- .expression_gene_map(expression)$map
    extent <- .native_kept_extent(expression$matrix, map)
    if (!is.finite(extent[[1L]]) || max(abs(extent[1:2])) <= .float32_exact_integer_limit) {
      return("<f4")
    }
    return(if (.native_kept_extent(expression$matrix, map, whole = TRUE)[[5L]] == 1) "<f8" else "<f4")
  }
  .expression_dtype(.expression_stored_values(expression), expression$extent)
}

# The values the dtype is decided on: a sparse matrix's stored values, a
# dense matrix's every value.
.expression_stored_values <- function(expression) {
  if (expression$sparse) expression$matrix@x else expression$matrix
}

# Write one expression matrix gene-major at `path`: a csc_matrix of cells x
# genes from a sparse matrix, a dense array in tall gene-major chunks from a
# dense one.
.expression_put <- function(store, path, expression, dtype) {
  matrix <- expression$matrix
  shape <- .expression_dims(expression)
  if (!expression$sparse) {
    # A genes x cells matrix is read block by block through its transpose,
    # never transposed whole; cells in an order are read through it.
    values <- if (expression$feature_major) .transposed(matrix) else matrix
    if (!is.null(expression$cells)) values <- .rows_in_order(values, expression$cells)
    .ad_dense(store, path, values, gene_major = TRUE, dtype = dtype)
    return(invisible(NULL))
  }
  if (isTRUE(expression$backed)) {
    arrays <- .backed_arrays(store, path, expression, shape, dtype)
    pointers <- arrays$pointers
    indices <- arrays$indices
    values <- arrays$values
  } else if (expression$feature_major) {
    arrays <- .gene_major_arrays(store, path, matrix, shape, dtype, .expression_gene_map(expression),
                                 cells = expression$cells)
    pointers <- arrays$pointers
    indices <- arrays$indices
    values <- arrays$values
  } else if (!is.null(expression$cells)) {
    sizes <- c(.zarr_item_size(.sparse_index_dtype("csc_matrix", shape)), .zarr_item_size(dtype))
    arrays <- .renumbered_arrays(
      store, path, matrix, .expression_new_rows(expression), sizes, dtype,
      block_values = .renumbered_block_values(store, shape[[1L]], length(matrix@x), sum(sizes))
    )
    pointers <- matrix@p
    indices <- arrays$indices
    values <- arrays$values
  } else {
    pointers <- matrix@p
    indices <- matrix@i
    values <- matrix@x
  }
  .ad_sparse(
    store, path, "csc_matrix", shape,
    pointers = pointers, indices = indices, values = values,
    dtype = dtype, sliced = TRUE
  )
}

# The gene-major csc_matrix of a genes x cells dgCMatrix (`shape` is cells x
# genes) without its transpose: the pointers counted from the row indices
# (one pass, src/extent.c), and the indices (each value's cell) and data
# made by one pass over the cells (src/chunks.c): each cell's values are
# ascending by gene, so its values in a block's genes are one run of its
# column, found by a binary search, plus a write per value it holds.
#
# cellucid_prepare() writes every chunk in order, so its store makes one block
# of every value, once, on every core, and holds it while X is written (8
# bytes a value with float32 data, where Matrix::t() made a dgCMatrix of 12).
# A served store
# lays the matrix out once, in the background (.laid_out_arrays()), as
# cellucid-python's served store does: a gene read made a block of the genes
# around it, kept within cache_bytes, and genes read across a matrix larger
# than the cache made a block again for nearly every gene (1M cells x 100
# values: 80 ms a gene, against 1 ms laid out).
.gene_major_arrays <- function(store, path, matrix, shape, dtype, genes, cells = NULL) {
  visit <- if (is.null(cells)) integer(0) else as.integer(cells - 1L)
  n_cells <- shape[[1L]]
  n_genes <- shape[[2L]]
  pointers <- .deferred(function() .native_row_pointers(matrix@i, genes$map, n_genes), length = n_genes + 1)
  # Every value when every gene is kept; else the kept genes' (one pass).
  n_values <- if (n_genes == nrow(matrix)) length(matrix@x) else .realize(pointers)[[n_genes + 1L]]
  index_dtype <- .sparse_index_dtype("csc_matrix", shape)
  sizes <- c(.zarr_item_size(index_dtype), .zarr_item_size(dtype))
  if (store$memory) {
    arrays <- .laid_out_arrays(store, path, n_values, sizes, function() list(
      0L, matrix@i, matrix@p, matrix@x, genes$map, genes$origin, .realize(pointers), visit, as.integer(sizes[[1L]]), dtype
    ))
    return(c(list(pointers = pointers), arrays))
  }
  # Every value in one block, made once, as X is written chunk by chunk.
  block_values <- max(1, n_values)
  held <- new.env(parent = emptyenv())
  block <- function(index) {
    if (is.null(held$block)) {
      held$block <- .native_gene_major_block(
        matrix, genes$map, genes$origin, .realize(pointers), 0, n_values, sizes[[1L]], dtype, visit,
        threads = .deflate_threads()
      )
    }
    held$block
  }
  list(pointers = pointers, indices = .block_array(n_values, block_values, sizes, block, 1L),
       values = .block_array(n_values, block_values, sizes, block, 2L))
}

# The gene-major arrays (indices, values) of an expression matrix of a
# served store, laid out once in the background (src/chunks.c,
# cellucid_layouts_start()): the matrix is queued now (`job()` its pass's
# arguments) and its pass starts when the store is closed (.store_close()),
# one matrix after another in the order they were planned, X first, each on
# every core, into memory-mapped scratch files (the operating system keeps
# their pages or writes them to disk; nothing of it in R's memory). A chunk,
# or a range of one, is a copy of the bytes it covers once the matrix is laid
# out. Until then the chunk is made whole from the layout's progress
# (cellucid_layout_block()): the layouts stand while it is made, the values
# already laid out are copied and only those of the cells (or genes) not laid
# out yet are made -- found by one scan of the matrix's values in its own
# order, on every core -- the bytes the layout will hold; the viewer waits
# for a layout only when it ends sooner. A matrix of no values has no chunk
# and is not queued.
.laid_out_arrays <- function(store, path, n_values, sizes, job) {
  index <- NA_integer_
  if (n_values > 0) {
    index <- length(store$layout_jobs)
    store$layout_jobs[[index + 1L]] <- job
    store$layout_labels <- c(store$layout_labels, path)
  }
  # The last chunk made before the layout ended, both its parts, whole: a
  # gene's values are asked for right after its indices, the next genes of a
  # chunk after it, and one pass makes every value a chunk holds (as
  # cellucid-python's served store keeps the chunk it made).
  made <- new.env(parent = emptyenv())
  part <- function(which) structure(list(
    length = n_values,
    chunk = function(first, count, chunk_length, chunk_first = first, chunk_last = first + count - 1) {
      size <- sizes[[which]]
      if (.layout_done(store, index)) {
        made$block <- NULL
        return(.native_layout_slice(store$layouts, index, which, (first - 1) * size, count * size, chunk_length * size))
      }
      start <- chunk_first - 1
      if (!identical(made$start, start)) {
        made$block <- .native_layout_block(store$layouts, index, start, chunk_last)
        made$start <- start
      }
      .native_raw_slice(made$block[[which]], (first - 1 - start) * size, count * size, chunk_length * size)
    }
  ), class = "cellucid_chunked")
  list(indices = part(1L), values = part(2L))
}

# Start laying out the matrices a memory store queued (.laid_out_arrays()):
# every job's arguments made now (the pointers counted, on R's thread), its
# scratch files in R's session folder (each removed once mapped), each
# matrix on every core (.deflate_threads(); two under R CMD check).
# `park` (tests): every thread stands that many units into its share until
# .native_layouts_unpark().
.layouts_start <- function(store, park = -1) {
  jobs <- lapply(store$layout_jobs, function(job) job())
  paths <- tempfile(rep(c("cellucid-layout-indices-", "cellucid-layout-values-"), length(jobs)))
  store$layouts <- .native_layouts_start(jobs, paths, .deflate_threads(), park)
  store$layout_jobs <- NULL
  invisible(store)
}

# Whether matrix `index` (0-based) of the store's layouts is laid out; stops
# with the reason when its layout failed (the viewer was stopped, or memory
# or disk ran out).
.layout_done <- function(store, index) {
  state <- .native_layout_state(store$layouts, index)
  if (state[[1L]] == 3) .abort(store$layout_labels[[index + 1L]], " could not be laid out gene-major for the viewer (the viewer was stopped, or out of memory or disk).")
  state[[1L]] == 2
}

# Wait until every matrix of the store's layouts is laid out (or `timeout`
# seconds), as cellucid-python's PriorityStore.laid_out() does; stops with a
# failed one's error.
.store_laid_out <- function(store, timeout = Inf) {
  if (is.null(store$layouts)) return(invisible(store))
  deadline <- proc.time()[["elapsed"]] + timeout
  for (index in seq_along(store$layout_labels) - 1L) {
    repeat {
      state <- .native_layout_state(store$layouts, index)
      if (state[[1L]] >= 2) break
      if (proc.time()[["elapsed"]] > deadline) .abort(store$layout_labels[[index + 1L]], " was not laid out within ", timeout, " s.")
      Sys.sleep(0.005)
    }
    if (state[[1L]] == 3) .abort(store$layout_labels[[index + 1L]], " could not be laid out gene-major for the viewer.")
  }
  invisible(store)
}

# One part (1: indices, 2: values) of gene-major arrays made by blocks of
# `block_values` positions: each chunk cut from its block.
.block_array <- function(n_values, block_values, sizes, block, part) {
  structure(list(
    length = n_values,
    chunk = function(first, count, chunk_length, ...) {
      position <- first - 1
      index <- position %/% block_values
      size <- sizes[[part]]
      .native_raw_slice(
        block(index)[[part]], (position - index * block_values) * size, count * size, chunk_length * size
      )
    }
  ), class = "cellucid_chunked")
}

# The gene-major arrays of a cells x genes dgCMatrix whose rows (cells) are
# renumbered (`new_row`, 0-based) -- an expression matrix in priority order,
# its genes in their order: each block's genes read whole, renumbered and
# sorted by new row (src/chunks.c), never a reordered copy of the matrix. A
# served store lays it out once in the background (.laid_out_arrays()).
# `sizes` are the bytes of an index and of a value (`dtype`).
.renumbered_arrays <- function(store, path, matrix, new_row, sizes, dtype, block_values) {
  n_values <- length(matrix@x)
  if (store$memory) {
    return(.laid_out_arrays(store, path, n_values, sizes, function() list(
      1L, matrix@i, matrix@p, matrix@x, as.integer(new_row), integer(0), matrix@p, as.integer(sizes[[1L]]), dtype
    )))
  }
  held <- new.env(parent = emptyenv())
  block <- function(index) {
    if (!identical(held$index, index)) {
      held$block <- NULL
      first <- index * block_values
      held$block <- .native_renumbered_block(matrix, new_row, first, min(n_values, first + block_values), sizes[[1L]], dtype)
      held$index <- index
    }
    held$block
  }
  list(indices = .block_array(n_values, block_values, sizes, block, 1L),
       values = .block_array(n_values, block_values, sizes, block, 2L))
}

# How many values one renumbered block of an expression matrix a written
# store writes holds: a whole number of the sliced chunks (a chunk never
# spans two blocks), at most eight values a cell. A block sorts every gene
# reaching into it whole, and a gene holds at most one value a cell, so at
# most one gene in eight is sorted twice. A written store holds one block at
# a time -- X's values are written block by block, then its cells, each block
# made once for each -- so writing in priority order holds one block (8
# bytes a value with float32 data), not a copy of X.
.renumbered_block_values <- function(store, n_cells, n_values, bytes_per_value) {
  chunk <- .sliced_chunk_elements(store)
  max(1, min(ceiling(n_values / chunk), ceiling(8 * n_cells / chunk))) * chunk
}

