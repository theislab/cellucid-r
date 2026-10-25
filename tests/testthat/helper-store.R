# Reading back what cellucid_prepare() wrote, with nothing but base R and
# jsonlite: the JSON metadata of one element, and the values of one array,
# decompressed chunk by chunk and decoded from its declared dtype. The reader
# is deliberately independent of the writer's code, so a test compares the
# store with the caller's input and not the writer with itself.

store_json <- function(store, key) {
  jsonlite::read_json(file.path(store, key), simplifyVector = FALSE)
}

decode_chunk <- function(bytes, dtype, count) {
  switch(dtype,
    "<f4" = readBin(bytes, "double", n = count, size = 4L, endian = "little"),
    "<f8" = readBin(bytes, "double", n = count, size = 8L, endian = "little"),
    "<i4" = readBin(bytes, "integer", n = count, size = 4L, endian = "little"),
    "<u4" = {
      words <- as.double(readBin(bytes, "integer", n = count, size = 4L, endian = "little"))
      ifelse(words < 0, words + 4294967296, words)
    },
    "<i2" = readBin(bytes, "integer", n = count, size = 2L, endian = "little"),
    "|i1" = readBin(bytes, "integer", n = count, size = 1L, endian = "little"),
    "|b1" = as.logical(as.integer(bytes[seq_len(count)])),
    "<i8" = {
      words <- readBin(bytes, "integer", n = 2L * count, size = 4L, endian = "little")
      low <- words[c(TRUE, FALSE)]
      low <- ifelse(low < 0, low + 4294967296, low)
      low + words[c(FALSE, TRUE)] * 4294967296
    },
    "|O" = {
      n <- readBin(bytes[1:4], "integer", size = 4L, endian = "little")
      out <- character(n)
      position <- 5L
      for (index in seq_len(n)) {
        length <- readBin(bytes[position:(position + 3L)], "integer", size = 4L, endian = "little")
        position <- position + 4L
        out[[index]] <- if (length == 0L) "" else rawToChar(bytes[position:(position + length - 1L)])
        Encoding(out[[index]]) <- "UTF-8"
        position <- position + length
      }
      out
    },
    stop("unknown dtype ", dtype)
  )
}

# What one array's metadata says, in either Zarr format: its dtype (as the
# Zarr format 2 name), shape, chunk shape, whether its chunks are compressed
# (a served store's are not), the byte size its values are shuffled by (0
# for none: format 2's shuffle filter, format 3's blosc codec), whether its
# chunks are Blosc frames, and where its chunks are.
array_meta <- function(store, key) {
  if (file.exists(file.path(store, key, ".zarray"))) {
    meta <- store_json(store, file.path(key, ".zarray"))
    shuffle <- Filter(function(filter) identical(filter$id, "shuffle"), meta$filters)
    return(list(
      dtype = meta$dtype,
      shape = unlist(meta$shape),
      chunks = unlist(meta$chunks),
      compressed = !is.null(meta$compressor),
      shuffle = if (length(shuffle)) shuffle[[1L]]$elementsize else 0L,
      blosc = FALSE,
      key = function(position) if (length(position)) paste(position, collapse = ".") else "0"
    ))
  }
  meta <- store_json(store, file.path(key, "zarr.json"))
  types <- c(
    float32 = "<f4", float64 = "<f8", int8 = "|i1", int16 = "<i2", int32 = "<i4", uint32 = "<u4",
    int64 = "<i8", bool = "|b1", string = "|O"
  )
  blosc <- Filter(function(codec) identical(codec$name, "blosc"), meta$codecs)
  list(
    dtype = types[[meta$data_type]],
    shape = unlist(meta$shape),
    chunks = unlist(meta$chunk_grid$configuration$chunk_shape),
    compressed = any(vapply(meta$codecs, function(codec) codec$name %in% c("blosc", "gzip"), logical(1L))),
    shuffle = if (length(blosc)) blosc[[1L]]$configuration$typesize else 0L,
    blosc = length(blosc) > 0L,
    key = function(position) paste(c("c", position), collapse = "/")
  )
}

# The bytes of a Blosc 1 frame of byte-shuffled, zlib-compressed blocks
# (format 3's blosc codec): the 16-byte header (version, codec version,
# flags, value size, then the byte size, block size and frame size as
# little-endian uint32), the start of every block, and in each block one
# stream per byte of a value (one for a last, shorter block, or with the
# no-split flag 0x10), each a 4-byte length and its zlib data, or the stream
# as it is when the length is the stream's size. A whole copy (flag 0x02)
# holds the bytes after the header, unshuffled.
blosc_decode <- function(frame) {
  word <- function(at) readBin(frame[at + 0:3], "integer", size = 4L, endian = "little")
  flags <- as.integer(frame[[3L]])
  item_size <- as.integer(frame[[4L]])
  n_bytes <- word(5L)
  block <- word(9L)
  stopifnot(as.integer(frame[[1L]]) == 2L, word(13L) == length(frame), bitwAnd(flags, 0x01L) == 1L)
  if (bitwAnd(flags, 0x02L) != 0L) {
    stopifnot(length(frame) == n_bytes + 16L)
    return(frame[16L + seq_len(n_bytes)])
  }
  split <- bitwAnd(flags, 0x10L) == 0L && item_size <= 16L && block %/% item_size >= 128L
  n_blocks <- ceiling(n_bytes / block)
  at <- 17L + 4L * n_blocks
  out <- raw(0)
  for (index in seq_len(n_blocks)) {
    stopifnot(word(17L + 4L * (index - 1L)) == at - 1L)
    size <- min(block, n_bytes - (index - 1L) * block)
    n_streams <- if (split && size == block) item_size else 1L
    shuffled <- raw(0)
    for (stream in seq_len(n_streams)) {
      length <- word(at)
      data <- frame[at + 3L + seq_len(length)]
      at <- at + 4L + length
      plane <- if (length == size / n_streams) data else memDecompress(data, type = "gzip")
      stopifnot(length(plane) == size / n_streams)
      shuffled <- c(shuffled, plane)
    }
    out <- c(out, as.vector(t(matrix(shuffled, ncol = item_size))))
  }
  stopifnot(at == length(frame) + 1L)
  out
}

read_chunk <- function(store, key, meta, position, count) {
  path <- file.path(store, key, meta$key(position))
  bytes <- readBin(path, "raw", file.size(path))
  if (!meta$compressed) {
    stopifnot(meta$shuffle == 0L)
  } else if (meta$blosc) {
    bytes <- blosc_decode(bytes)
  } else {
    bytes <- memDecompress(bytes, type = "gzip")
    if (meta$shuffle > 1L) {
      bytes <- as.vector(t(matrix(bytes, ncol = meta$shuffle)))
    }
  }
  decode_chunk(bytes, meta$dtype, count)
}

# Every value of a 0-, 1- or 2-dimensional array, as a vector or a matrix.
store_array <- function(store, key) {
  meta <- array_meta(store, key)
  shape <- meta$shape
  chunks <- meta$chunks
  dtype <- meta$dtype
  if (length(shape) == 0L) {
    return(read_chunk(store, key, meta, integer(0), 1L))
  }
  if (length(shape) == 1L) {
    if (shape == 0) return(decode_chunk(raw(0), dtype, 0L))
    values <- unlist(lapply(seq_len(ceiling(shape / chunks)) - 1L, function(index) {
      read_chunk(store, key, meta, index, chunks)
    }))
    return(values[seq_len(shape)])
  }
  out <- matrix(NA_real_, nrow = shape[[1L]], ncol = shape[[2L]])
  for (row_chunk in seq_len(ceiling(shape[[1L]] / chunks[[1L]])) - 1L) {
    for (col_chunk in seq_len(ceiling(shape[[2L]] / chunks[[2L]])) - 1L) {
      values <- read_chunk(store, key, meta, c(row_chunk, col_chunk), prod(chunks))
      block <- matrix(values, nrow = chunks[[1L]], ncol = chunks[[2L]], byrow = TRUE)
      rows <- row_chunk * chunks[[1L]] + seq_len(chunks[[1L]])
      cols <- col_chunk * chunks[[2L]] + seq_len(chunks[[2L]])
      keep_rows <- rows <= shape[[1L]]
      keep_cols <- cols <= shape[[2L]]
      out[rows[keep_rows], cols[keep_cols]] <- block[keep_rows, keep_cols, drop = FALSE]
    }
  }
  out
}

# Values as float32 holds them.
as_float32 <- function(values) {
  readBin(writeBin(as.double(values), raw(), size = 4L), "double", n = length(values), size = 4L)
}

# A small dataset with every kind of input the writer handles.
prepare_fixture <- function(out_path) {
  expression <- Matrix::sparseMatrix(
    i = c(1, 2, 4, 6, 1, 3, 5, 2, 6),
    j = c(1, 1, 1, 1, 2, 2, 3, 4, 4),
    x = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
    dims = c(6L, 4L),
    dimnames = list(paste0("cell_", 1:6), c("ENSG1", "ENSG2", "ENSG3", "ENSG4"))
  )
  knn <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 5, 6), j = c(2, 1, 4, 3, 6, 5),
    x = c(0.25, 0.25, 0.5, 0.5, 1, 1), dims = c(6L, 6L)
  )
  list(
    out_path = out_path,
    format = "zarr",
    embeddings = list(
      X_umap = matrix(c(0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0) / 4, ncol = 2),
      X_pca = matrix(seq_len(18) / 8, ncol = 3)
    ),
    dataset_id = "fixture",
    dataset_name = "Fixture",
    obs = data.frame(
      cluster = factor(c("b", "a", "b", "c", "a", "b"), levels = c("c", "b", "a")),
      batch = c("x", "y", "x", "y", "x", "y"),
      barcode = c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC"),
      score = c(0.5, NA, 1.25, -2, 3, 0),
      flag = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
      maybe = c(TRUE, NA, FALSE, TRUE, NA, FALSE),
      n_genes = c(10L, 20L, NA, 40L, 50L, 60L),
      stringsAsFactors = FALSE
    ),
    var = data.frame(
      symbol = c("CD3E", "CD8A", "MS4A1", "CD3E"),
      row.names = c("ENSG1", "ENSG2", "ENSG3", "ENSG4")
    ),
    gene_expression = expression,
    connectivities = knn
  )
}

as_float32_matrix <- function(matrix) {
  out <- matrix(as_float32(matrix), nrow(matrix))
  dimnames(out) <- dimnames(matrix)
  out
}
