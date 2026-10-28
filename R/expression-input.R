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
# matrix, every value finite and within the float32 range.
.expression_input <- function(x, name) {
  feature_major <- inherits(x, "cellucid_feature_major")
  matrix <- if (feature_major) x$matrix else x
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
  if (expression$feature_major) shape <- rev(shape)
  if (!is.null(expression$genes)) shape[[2L]] <- length(expression$genes)
  shape
}

.expression_cell_names <- function(expression) {
  if (expression$feature_major) colnames(expression$matrix) else rownames(expression$matrix)
}

.expression_gene_names <- function(expression) {
  names <- if (expression$feature_major) rownames(expression$matrix) else colnames(expression$matrix)
  if (is.null(expression$genes) || is.null(names)) names else names[expression$genes]
}

# The genes kept, by position (gene_identifiers). A genes x cells sparse
# matrix keeps them as positions read by its gene-major arrays
# (.gene_major_arrays()), never as a copy of its rows: at a million cells the
# row subset copied every value of the matrix. Every other matrix is
# narrowed to its columns (cells x genes) or rows (dense genes x cells).
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
# value, so the measured extent still holds.
.expression_select_cells <- function(expression, rows) {
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
  if (expression$feature_major) {
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
# made one block of positions at a time (src/chunks.c), each chunk then cut
# from its block. A block costs one binary search per cell -- each cell's
# values are ascending by gene, so its values in the block's genes are one
# run of its column -- plus a write per value it holds; it never reads the
# rest of the matrix's values.
#
# cellucid_prepare() writes every chunk in order, so its store makes one block
# of every value, once, and holds it while X is written (8 bytes a value with
# float32 data, where Matrix::t() made a dgCMatrix of 12). A served store
# makes a block when the viewer first reads one of its genes: the first gene
# costs one block, not the whole transpose, and a block is kept in the
# store's cache, so the blocks held never exceed cache_bytes.
.gene_major_arrays <- function(store, path, matrix, shape, dtype, genes, cells = NULL) {
  visit <- if (is.null(cells)) integer(0) else as.integer(cells - 1L)
  n_cells <- shape[[1L]]
  n_genes <- shape[[2L]]
  pointers <- .deferred(function() .native_row_pointers(matrix@i, genes$map, n_genes), length = n_genes + 1)
  # Every value when every gene is kept; else the kept genes' (one pass).
  n_values <- if (n_genes == nrow(matrix)) length(matrix@x) else .realize(pointers)[[n_genes + 1L]]
  index_dtype <- .sparse_index_dtype("csc_matrix", shape)
  sizes <- c(.zarr_item_size(index_dtype), .zarr_item_size(dtype))
  block_values <- .gene_major_block_values(store, n_cells, n_values, sum(sizes))
  held <- new.env(parent = emptyenv())
  block <- function(index) {
    make <- function() {
      first <- index * block_values
      .native_gene_major_block(
        matrix, genes$map, genes$origin, .realize(pointers), first, min(n_values, first + block_values),
        sizes[[1L]], dtype, visit
      )
    }
    if (store$memory) {
      # A control character: never one of the store's own keys.
      return(.cached(store$cache, paste0(path, "\001", index), make))
    }
    if (!identical(held$index, index)) {
      held$block <- NULL
      held$block <- make()
      held$index <- index
    }
    held$block
  }
  list(pointers = pointers, indices = .block_array(n_values, block_values, sizes, block, 1L),
       values = .block_array(n_values, block_values, sizes, block, 2L))
}

# One part (1: indices, 2: values) of gene-major arrays made by blocks of
# `block_values` positions: each chunk cut from its block.
.block_array <- function(n_values, block_values, sizes, block, part) {
  structure(list(
    length = n_values,
    chunk = function(first, count, chunk_length) {
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
# sorted by new row (src/chunks.c), never a reordered copy of the matrix; the
# blocks kept as .gene_major_arrays() keeps them. `sizes` are the bytes of an
# index and of a value (`dtype`).
.renumbered_arrays <- function(store, path, matrix, new_row, sizes, dtype, block_values) {
  n_values <- length(matrix@x)
  held <- new.env(parent = emptyenv())
  block <- function(index) {
    make <- function() {
      first <- index * block_values
      .native_renumbered_block(matrix, new_row, first, min(n_values, first + block_values), sizes[[1L]], dtype)
    }
    if (store$memory) return(.cached(store$cache, paste0(path, "\001", index), make))
    if (!identical(held$index, index)) {
      held$block <- NULL
      held$block <- make()
      held$index <- index
    }
    held$block
  }
  list(indices = .block_array(n_values, block_values, sizes, block, 1L),
       values = .block_array(n_values, block_values, sizes, block, 2L))
}

# How many values one renumbered block of an expression matrix holds: a whole
# number of the sliced chunks (a chunk never spans two blocks), at most eight
# values a cell. A block sorts every gene reaching into it whole, and a gene
# holds at most one value a cell, so at most one gene in eight is sorted
# twice. A written store holds one block at a time -- X's values are written
# block by block, then its cells, each block made once for each -- so writing
# in priority order holds one block (8 bytes a value with float32 data), not
# a copy of X; a served store keeps its blocks in its cache, each at most
# half of it.
.renumbered_block_values <- function(store, n_cells, n_values, bytes_per_value) {
  chunk <- .chunk_elements_sliced
  blocks <- min(ceiling(n_values / chunk), ceiling(8 * n_cells / chunk))
  if (store$memory) blocks <- min(blocks, floor(store$cache$budget / 2 / bytes_per_value / chunk))
  max(1, blocks) * chunk
}

# How many values one gene-major block holds: a whole number of the sliced
# chunks (a chunk never spans two blocks). A written store holds every value
# in one block. A served store's block is the smaller of two sizes. One is
# eight values a cell: a block's search costs about 50 ns a cell, and writing
# its values about 9 ns a value (measured on one million cells of 150
# values), so a block of eight values a cell spends about 40% of its time
# searching, and the first gene of a million cells costs about 0.12 s rather
# than the 1.3 s of a whole transpose. The other is half of cache_bytes, so
# two blocks are always kept.
.gene_major_block_values <- function(store, n_cells, n_values, bytes_per_value) {
  chunk <- .chunk_elements_sliced
  whole <- max(1, ceiling(n_values / chunk))
  if (!store$memory) {
    return(whole * chunk)
  }
  per_cells <- ceiling(8 * n_cells / chunk)
  in_cache <- floor(store$cache$budget / 2 / bytes_per_value / chunk)
  max(1, min(whole, per_cells, in_cache)) * chunk
}
