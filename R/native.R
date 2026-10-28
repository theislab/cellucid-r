# The complete surface of src/ (registered in src/init.c): every .Call() of
# the package, so the R side of the native contract is readable in one screen.

# src/crc32.c: the CRC-32 of a raw vector, continuing from `previous` (the CRC
# of the bytes before it, 0 for none). The ZIP writer in zip.R needs it, and
# so does the gzip framing of a Zarr format 3 chunk (zarr-store.R).
.native_crc32 <- function(bytes, previous = 0) {
  .Call(C_cellucid_crc32, bytes, as.double(previous))
}

# src/random.c: `n` bytes from the operating system's cryptographic random
# generator; R's own generator (and the caller's .Random.seed) is not used.
.native_random_bytes <- function(n) {
  .Call(C_cellucid_random_bytes, as.double(n))
}

# src/sha256.c: the SHA-256 of a raw vector, as lower-case hex; the web
# build's inventory declares one for every file (web-cache.R).
.native_sha256 <- function(bytes) {
  .Call(C_cellucid_sha256, bytes)
}

# src/extent.c: one pass over a numeric vector: c(smallest finite value,
# largest finite value, count of values that are not finite, 1 when one is
# zero, 1 when every finite value is whole -- checked only when `whole`).
.native_extent <- function(values, whole = FALSE) {
  .Call(C_cellucid_extent, values, whole)
}

# src/extent.c: the pointers of a dgCMatrix's transpose over the rows `map`
# keeps (map[row + 1] is its 0-based output row, or -1), from its row indices
# alone: the cumulative count of stored values per output row.
.native_row_pointers <- function(i, map, n_rows) {
  .Call(C_cellucid_row_pointers, i, map, as.double(n_rows))
}

# src/extent.c: .native_extent() over the values of the rows `map` keeps.
.native_kept_extent <- function(matrix, map, whole = FALSE) {
  .Call(C_cellucid_kept_extent, matrix@x, matrix@i, map, whole)
}

# src/latent_outliers.c: one categorical obs column's outlier quantiles in the
# latent space (R/outlier-quantiles.R).
.native_latent_outlier_quantiles <- function(latent, codes, n_categories, min_points) {
  .Call(C_cellucid_latent_outlier_quantiles, latent, codes, n_categories, min_points)
}

# src/session_codecs.c: one gzip (`gzip = TRUE`) or zlib stream inflated to
# exactly `expected` bytes (negative: any size), refusing a truncated,
# longer, shorter or concatenated stream. Base R's memDecompress() never
# returns on two concatenated gzip members, so every inflate goes here.
.native_inflate <- function(bytes, gzip, expected = -1) {
  .Call(C_cellucid_inflate, bytes, gzip, as.double(expected))
}

# src/session_codecs.c: the web app's order-independent 32-bit digest of a
# set of 0-based cell indices (codecs/cell-set.js digestCellSet).
.native_digest_cell_set <- function(indices) {
  .Call(C_cellucid_digest_cell_set, indices)
}

# src/session_codecs.c: one highlight membership of a session as 0-based
# cell indices, ascending; `reference` the 0-based cells of the group's
# category when the set is stored relative to it, else NULL.
.native_decode_cell_set <- function(payload, cell_count, reference = NULL) {
  .Call(C_cellucid_decode_cell_set, payload, as.double(cell_count), reference)
}

# src/session_codecs.c: one user-defined categorical field's codes of a
# session, as doubles (a Uint32Array code can exceed R's integers).
.native_decode_user_codes <- function(payload, expected_length, codes_type) {
  .Call(C_cellucid_decode_user_codes, payload, as.double(expected_length), codes_type)
}

# src/chunks.c: one chunk's bytes -- values first .. first + count - 1 of
# `values` (1-based) in `dtype`, little-endian, then zeros to `chunk_length`
# values (the fill value of every fixed-width dtype the package writes).
.native_encode_values <- function(values, first, count, chunk_length, dtype) {
  .Call(C_cellucid_encode_values, values, as.double(first), as.double(count), as.double(chunk_length), dtype)
}

# src/chunks.c: positions first .. last - 1 (0-based) of the gene-major
# arrays of a genes x cells dgCMatrix over the genes it keeps (`map`: each
# row's 0-based output gene or -1; `origin`: each output gene's 0-based row),
# as list(cell index bytes, value bytes), made without transposing the matrix
# or copying the rows it keeps (R/expression-input.R).
# `cells`: the 0-based columns in the order the output rows take them
# (integer(0): the matrix's own).
.native_gene_major_block <- function(matrix, map, origin, pointers, first, last, index_size, dtype,
                                     cells = integer(0)) {
  .Call(C_cellucid_gene_major_block, matrix@i, matrix@p, matrix@x, map, origin, pointers, as.integer(cells),
        as.double(first), as.double(last), as.integer(index_size), dtype)
}

# src/chunks.c: positions first .. last - 1 (0-based) of the column-compressed
# arrays of a dgCMatrix whose rows are renumbered (`new_row`, 0-based) and
# whose columns are taken in `column_order` (0-based; integer(0) for their
# own), each column's values ascending by new row; `pointers` are the
# output's column pointers (the matrix's own @p when its columns keep their
# order). A cells x genes expression matrix in priority order, or both axes
# of a symmetric graph. `part` 0 returns list(index bytes, value bytes); 1 or
# 2 only the indices or only the values, zero-padded to `padded` items (one
# chunk).
.native_renumbered_block <- function(matrix, new_row, first, last, index_size, dtype,
                                     column_order = integer(0), pointers = matrix@p, part = 0L, padded = 0) {
  .Call(C_cellucid_renumbered_block, matrix@i, matrix@p, matrix@x, as.integer(new_row),
        as.integer(column_order), pointers, as.double(first), as.double(last), as.integer(index_size), dtype,
        as.integer(part), as.double(padded))
}

# src/chunks.c: bytes offset .. offset + count - 1 (0-based) of `bytes`,
# then zeros to `total` bytes.
.native_raw_slice <- function(bytes, offset, count, total) {
  .Call(C_cellucid_raw_slice, bytes, as.double(offset), as.double(count), as.double(total))
}

# src/chunks.c: one binary command frame -- the header's length, the header
# bytes, then `cells` as little-endian uint32 -- or NULL when the cells are
# not whole numbers from 0 to 2^32 - 1, strictly ascending, at least one
# (none allowed with `allow_empty`).
.native_cells_frame <- function(header, cells, allow_empty = FALSE) {
  .Call(C_cellucid_cells_frame, header, cells, allow_empty)
}

# src/subsample.c: the subsampling specification (R/subsample.R). Positions
# are 0-based doubles both ways (a dataset holds up to 2^32 - 1 cells).
.native_subsample_rank <- function(seed, positions) {
  .Call(C_cellucid_subsample_rank, as.double(seed), as.double(positions))
}
.native_subsample_uniform <- function(n_obs, k, seed) {
  .Call(C_cellucid_subsample_uniform, as.double(n_obs), as.double(k), as.double(seed))
}
.native_subsample_blocks <- function(n_obs, k, seed, block_size) {
  .Call(C_cellucid_subsample_blocks, as.double(n_obs), as.double(k), as.double(seed), as.double(block_size))
}
.native_subsample_split <- function(strata, n_strata, k, seed, allocation, minimum) {
  .Call(C_cellucid_subsample_split, as.double(strata), as.double(n_strata), as.double(k), as.double(seed),
        as.character(allocation), as.double(minimum))
}
.native_subsample_boxes <- function(coordinates, dimensions, grid_size) {
  .Call(C_cellucid_subsample_boxes, as.double(coordinates), as.integer(dimensions), as.double(grid_size))
}
.native_subsample_digest <- function(n_obs, kept) {
  .Call(C_cellucid_subsample_digest, as.double(n_obs), if (is.null(kept)) NULL else as.double(kept))
}
.native_subsample_encode <- function(kept) {
  .Call(C_cellucid_subsample_encode, as.double(kept))
}
.native_subsample_decode <- function(bytes, n_obs) {
  .Call(C_cellucid_subsample_decode, bytes, as.double(n_obs))
}

# src/priority.c: each gene's stored values in the cells before each prefix
# level (R/priority.R). `matrix` is a dgCMatrix, genes x cells
# (`cells_by_column`, `map` each row's 0-based output gene, -1 for none) or
# cells x genes. `cells` takes a genes x cells matrix's columns in that order
# (0-based), `new_row` renumbers a cells x genes matrix's rows (0-based).
.native_prefix_counts <- function(matrix, levels, n_genes, cells_by_column, map = integer(0),
                                  cells = integer(0), new_row = integer(0)) {
  .Call(C_cellucid_prefix_counts, matrix@p, matrix@i, as.double(levels), as.double(n_genes),
        cells_by_column, as.integer(map), as.integer(cells), as.integer(new_row))
}
