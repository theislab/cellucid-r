# A Zarr store, written by this package itself, in Zarr format 2 or 3: the
# JSON metadata documents (format 2: .zgroup, .zattrs, .zarray; format 3: one
# zarr.json per group and array), one file per chunk, and the consolidated
# metadata (format 2: .zmetadata; format 3: the consolidated_metadata member
# of the root zarr.json) that lets a web host open the whole store with one
# request. The store goes either to a folder or, entry by entry, into one
# stored (uncompressed) ZIP archive whose top level is the store, which is
# what a .zarr.zip is.
#
# A written store (cellucid_prepare()) compresses every chunk. Every chunk of
# a multi-byte array is byte-shuffled before it is compressed -- the first
# byte of every value, then every second byte, and so on -- so the similar
# high bytes of neighbouring values compress together; the expression values
# (X's data, a dense X) are not, because deflate finds their repeated values
# only unshuffled (shuffled, real expression data grew 14-45%). Every chunk is
# compressed with DEFLATE level 6 (zlib's default) into a zlib stream (RFC
# 1950: a 0x78 header, deflate data, an Adler-32 trailer) by the zlib this
# package links, a batch of chunks at a time on every core (src/deflate.c;
# each chunk alone, so the bytes never depend on the number of cores): what
# base R's memCompress(type = "gzip") and numcodecs' Zlib codec write (the
# same bytes wherever they link a zlib that deflates alike; zlib 1.3.2
# deflates some inputs differently from 1.2.12, so the writers are compared by
# metadata and values, and a chunk by its bytes once inflated). Format 2 declares numcodecs' Shuffle filter {"id":
# "shuffle", "elementsize": n} (a format 2 filter is any numcodecs codec) and
# the compressor {"id": "zlib", "level": 6}. Format 3 has no shuffle codec of
# its own (numcodecs.shuffle is outside its specification), so a shuffled
# array is compressed with the specification's blosc codec -- byte shuffle and
# zlib in a Blosc frame (.blosc_frame()) -- and every other array with its
# gzip codec, the gzip framing (RFC 1952) around the same deflate data: a
# ten-byte header, then the CRC-32 and length of the uncompressed bytes.
# Browsers inflate all of these natively (DecompressionStream). These are the
# Python writer's choices (prepare_data/_store.py), so both writers produce the
# same store.
#
# A store held in memory (cellucid_view(), cellucid_serve()) declares no
# compressor and no shuffle filter, and its chunks are the values' bytes as
# they are. It is read over the loopback interface, where compressing costs
# far more than the bytes it saves: deflate runs on R's one main thread at
# about 20 MB/s, the same thread that answers every request, while the browser
# reads loopback at gigabytes a second. The shuffle is dropped with the
# compressor, because it only ever helped the compressor and would cost a
# permutation of every chunk on that thread.
#
# The store is a key plan: every element puts its keys, the metadata
# documents as bytes and every chunk as a function that makes its bytes. A
# folder or a ZIP backend calls each chunk's function at once, compresses its
# bytes and writes them (cellucid_prepare()); the memory backend keeps the
# functions and calls one when its key is asked for, so a served store holds
# exactly the values cellucid_prepare() writes, in the same arrays and chunks
# but for an expression matrix's entries, which it chunks for loopback
# (.chunk_elements_served) and lays out once in the background
# (expression-input.R, .laid_out_arrays()) (tests/testthat/test-served-store.R
# reads every array of both and compares them, dtype, shape and values); a
# Range of one of its chunks is made alone (.store_body()). A chunk is made from its R vector in one
# pass (src/chunks.c), straight into the vector the server answers with; R
# gives no raw vector that shares an integer or double vector's memory, so
# that one copy is the least a chunk can cost. Only chunks made by a
# transform (text, a dense matrix's block) are kept, in the store's
# least-recently-used cache.

.deflate_level <- 6L

# How many elements one chunk holds. A chunk is the unit of every read, so its
# size is a scheduling choice: large enough that a whole-column read is a
# handful of requests, small enough that reading one gene out of a sparse
# matrix does not fetch much more than that gene. Arrays with one row per cell
# (obs columns and index, embeddings, masks, cell maps, outlier quantiles) use
# 2^16 elements per chunk, so a subsample -- the first rows of a store in
# priority order -- reads little beyond its own cells; the other arrays read
# whole (var, X's pointers, the neighbour graph) 2^18 (1 MiB of float32); the
# data and indices of the expression matrix, read one gene's range (or its
# first values) at a time, 2^14. These are the Python writer's constants
# (CELL_CHUNK_ELEMENTS, WHOLE_CHUNK_ELEMENTS, SLICED_CHUNK_ELEMENTS in
# prepare_data/_store.py), with the measurements behind them.
.chunk_elements_whole <- 262144
.chunk_elements_cell <- 65536
.chunk_elements_sliced <- 16384

# A served store (held in memory) chunks the data and indices of its
# expression matrices in 2^18 values (1 MiB of float32): a request on this
# machine costs a fixed time whatever its length, at both ends -- measured on
# an Apple M1 Pro over loopback, six connections: these handlers on R's thread
# answer about 1,600 requests a second (0.6 ms each), cellucid-python's
# http.server about 4,700, and a browser's fetch from a worker moves 46-79
# MB/s in requests of 64 KiB (Chromium, Firefox, WebKit), 300-510 MB/s in
# requests of 1 MiB and no more in larger ones -- so at 2^14 values a chunk a
# pass over a matrix (an analysis, many genes) was bound by requests. 2^18
# values is the least at which every engine's loopback rate has reached its
# plateau. The viewer reads a served chunk in part (it is stored uncompressed: a
# Range of one is only those bytes, made alone, .store_body()), so a gene, or
# a subsample's first entries of one, reads no more than with smaller chunks.
# A scheduling choice, measured; cellucid-python's served stores use the same
# (SERVED_ENTRY_CHUNK_ELEMENTS in _served_store.py). A written store keeps
# .chunk_elements_sliced: its chunks are compressed, and a chunk is the least
# a reader of it decodes.
.chunk_elements_served <- 262144

# The values one chunk of an expression matrix's data and indices holds in
# `store`: a served store's (held in memory), else a written one's.
.sliced_chunk_elements <- function(store) {
  if (isTRUE(store$memory)) .chunk_elements_served else .chunk_elements_sliced
}

# --- JSON -------------------------------------------------------------------

# Vectors that must stay JSON arrays even with one element are wrapped in
# I(); NULL is written as null.
.json_text <- function(value) {
  as.character(jsonlite::toJSON(
    value,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = NA
  ))
}

.json_object <- function() structure(list(), names = character())

# --- The store --------------------------------------------------------------

# `path` is the staging location: a folder, or the file the ZIP is written to.
.store_open <- function(path, zip, zarr_format) {
  store <- new.env(parent = emptyenv())
  store$memory <- FALSE
  store$compressed <- TRUE
  store$zip <- zip
  store$path <- path
  store$zarr_format <- zarr_format
  # Format 2: every metadata document's key and text, for .zmetadata.
  # Format 3: every node's path and its consolidated document's text, and the
  # root group's attributes, written last in the root zarr.json.
  store$metadata_keys <- character()
  store$metadata_texts <- character()
  store$root_attrs <- .json_object()
  # Chunks made and not yet compressed (.zarr_put_chunk()), and their bytes.
  store$pending <- list()
  store$pending_bytes <- 0
  if (zip) {
    store$connection <- file(path, open = "wb")
    store$offset <- 0
    store$central <- list()
    store$deferred <- list()
  } else {
    .dir_create(path)
  }
  store
}

# A store held in memory: every key and its bytes, or the function that
# produces them (a chunk), in the order they were put. Nothing is written, and
# nothing is compressed. `cache_bytes` bounds what the store keeps of the
# chunks it makes by a transform (.chunk_cache()).
.store_open_memory <- function(zarr_format, cache_bytes = 0) {
  store <- new.env(parent = emptyenv())
  store$memory <- TRUE
  store$compressed <- FALSE
  store$cache <- .chunk_cache(cache_bytes)
  store$zip <- FALSE
  store$zarr_format <- zarr_format
  store$metadata_keys <- character()
  store$metadata_texts <- character()
  store$root_attrs <- .json_object()
  store$entries <- new.env(parent = emptyenv(), hash = TRUE)
  store$keys <- character()
  # The expression matrices laid out in the background (.laid_out_arrays()).
  store$layout_jobs <- list()
  store$layout_labels <- character()
  store
}

# One key of a memory store as a server answers it -- its length, and its
# bytes [start, end) made when asked (`read`) -- or NULL for a key it does
# not hold. A chunk that is a slice of its vector (.ranged()) is made on each
# request, and only the values a range of it covers: making it is one copy at
# memory speed, and keeping it would hold a second copy of the data; a Range
# of a served chunk (a gene, a subsample's first entries of one) makes those
# bytes alone. A chunk made by a transform is made whole and kept in the
# store's cache. The entries are keyed by .symbol_keys(): a key names an obs
# column, which may hold any character.
.store_body <- function(store, key) {
  value <- get0(.symbol_keys(key), envir = store$entries, inherits = FALSE)
  if (is.null(value)) return(NULL)
  item <- attr(value, "item")
  if (is.function(value) && !is.null(item)) {
    return(list(size = attr(value, "elements") * item, read = function(start, end) {
      from <- start %/% item
      to <- ceiling(end / item)
      bytes <- value(from, to)
      if (start == from * item && end == to * item) bytes else .native_raw_slice(bytes, start - from * item, end - start, end - start)
    }))
  }
  .raw_body(if (!is.function(value)) value else if (isTRUE(attr(value, "transform"))) .cached(store$cache, key, value) else value())
}

# Bytes held whole, as a body (.store_body()).
.raw_body <- function(bytes) {
  list(size = length(bytes), read = function(start, end) {
    if (start == 0 && end == length(bytes)) bytes else .native_raw_slice(bytes, start, end - start, end - start)
  })
}

# One key of a memory store whole (.store_body()), or NULL.
.store_get <- function(store, key) {
  body <- .store_body(store, key)
  if (is.null(body)) NULL else body$read(0, body$size)
}

# A chunk of a vector whose values [from, to) (offsets in the chunk; the
# whole chunk by default) `make(from, to)` makes alone, padded to `to - from`:
# a slice of a numeric vector or of a gene-major array, `elements` values of
# `item` bytes a chunk (.store_body()).
.ranged <- function(make, elements, item) {
  structure(make, elements = elements, item = item)
}

# Bytes made by a transform, kept least recently used first out within
# `budget` bytes: a chunk of text, a dense matrix's block. Bytes larger than
# the budget are made each time and not kept.
.chunk_cache <- function(budget) {
  cache <- new.env(parent = emptyenv())
  cache$budget <- budget
  cache$used <- 0
  cache$entries <- new.env(parent = emptyenv(), hash = TRUE)
  cache$sizes <- numeric(0)
  cache$touched <- numeric(0)
  cache$clock <- 0
  cache
}

# What `make()` returns for `key` (a raw vector), from the cache when it is there.
.cached <- function(cache, key, make) {
  key <- .symbol_keys(key)
  cache$clock <- cache$clock + 1
  hit <- get0(key, envir = cache$entries, inherits = FALSE)
  if (!is.null(hit)) {
    cache$touched[[key]] <- cache$clock
    return(hit)
  }
  value <- make()
  size <- length(value)
  if (size > cache$budget) return(value)
  while (cache$used + size > cache$budget && length(cache$sizes) > 0L) {
    oldest <- names(cache$touched)[[which.min(cache$touched)]]
    cache$used <- cache$used - cache$sizes[[oldest]]
    rm(list = oldest, envir = cache$entries)
    cache$sizes <- cache$sizes[names(cache$sizes) != oldest]
    cache$touched <- cache$touched[names(cache$touched) != oldest]
  }
  assign(key, value, envir = cache$entries)
  cache$sizes[[key]] <- size
  cache$touched[[key]] <- cache$clock
  cache$used <- cache$used + size
  value
}

# In a .zarr.zip the metadata documents are held back and written after every
# chunk (.store_close), the consolidated metadata last: a reader opening the
# archive reads its end, and with it the documents that describe the store.
# The Python writer orders a .zarr.zip's entries the same way. `bytes` is a raw
# vector, or for a memory store a function returning one (a chunk), which it
# keeps.
.store_put <- function(store, key, bytes) {
  # Chunks still to be compressed go first, so a ZIP holds every entry in the
  # order it was put.
  if (length(store$pending) > 0L) .store_flush(store)
  if (store$memory) {
    slot <- .symbol_keys(key)
    if (!exists(slot, envir = store$entries, inherits = FALSE)) {
      store$keys <- c(store$keys, key)
    }
    assign(slot, bytes, envir = store$entries)
    return(invisible(NULL))
  }
  if (store$zip && basename(key) %in% .zarr_metadata_names) {
    store$deferred[[length(store$deferred) + 1L]] <- list(key = key, bytes = bytes)
  } else if (store$zip) {
    .zip_add_entry(store, key, bytes)
  } else {
    target <- file.path(store$path, key)
    .dir_create(dirname(target))
    writeBin(bytes, target)
  }
  invisible(NULL)
}

.store_remember <- function(store, key, text) {
  store$metadata_keys <- c(store$metadata_keys, key)
  store$metadata_texts <- c(store$metadata_texts, text)
  invisible(NULL)
}

# A format 2 document, written and remembered for .zmetadata.
.store_put_metadata <- function(store, key, value) {
  text <- .json_text(value)
  .store_put(store, key, .utf8_bytes(text))
  .store_remember(store, key, text)
}

# The consolidated metadata, written last: in a .zarr.zip it is then the
# archive's last entry, read together with the central directory when the
# archive opens.
.store_close <- function(store) {
  .store_flush(store)
  if (store$zip) {
    # The documents of the root group last, after every nested one.
    deferred <- store$deferred
    at_root <- vapply(deferred, function(entry) !grepl("/", entry$key, fixed = TRUE), logical(1))
    for (entry in c(deferred[!at_root], deferred[at_root])) {
      .zip_add_entry(store, entry$key, entry$bytes)
    }
    store$deferred <- list()
  }
  members <- paste0(
    vapply(store$metadata_keys, .json_text, character(1), USE.NAMES = FALSE),
    ":",
    store$metadata_texts,
    collapse = ","
  )
  if (store$zarr_format == 2L) {
    text <- paste0("{\"metadata\":{", members, "},\"zarr_consolidated_format\":1}")
    .store_put_last(store, ".zmetadata", .utf8_bytes(text))
  } else {
    # The root group's own document, with every other node's document inline,
    # as zarr.consolidate_metadata() writes it.
    text <- paste0(
      "{\"attributes\":", .json_text(store$root_attrs),
      ",\"zarr_format\":3,\"consolidated_metadata\":{\"kind\":\"inline\",",
      "\"must_understand\":false,\"metadata\":{", members, "}},",
      "\"node_type\":\"group\"}"
    )
    .store_put_last(store, "zarr.json", .utf8_bytes(text))
  }
  if (store$zip) {
    .zip_finish(store)
    close(store$connection)
    store$connection <- NULL
  }
  if (isTRUE(store$memory) && length(store$layout_jobs) > 0L) .layouts_start(store)
  invisible(NULL)
}

# The consolidated metadata: the store's last file, or the archive's last entry.
.store_put_last <- function(store, key, bytes) {
  if (store$memory) {
    .store_put(store, key, bytes)
  } else if (store$zip) {
    .zip_add_entry(store, key, bytes)
  } else {
    writeBin(bytes, file.path(store$path, key))
  }
  invisible(NULL)
}

.store_abort <- function(store) {
  store$pending <- list()
  if (isTRUE(store$zip) && !is.null(store$connection)) {
    close(store$connection)
    store$connection <- NULL
  }
  invisible(NULL)
}

.zarr_key <- function(...) {
  parts <- c(...)
  paste(parts[nzchar(parts)], collapse = "/")
}

# --- Groups and arrays ------------------------------------------------------

.zarr_group <- function(store, path, attrs = NULL) {
  if (store$zarr_format == 2L) {
    .store_put_metadata(store, .zarr_key(path, ".zgroup"), list(zarr_format = 2L))
    if (!is.null(attrs)) {
      .store_put_metadata(store, .zarr_key(path, ".zattrs"), attrs)
    }
    return(invisible(NULL))
  }
  if (is.null(attrs)) attrs <- .json_object()
  if (!nzchar(path)) {
    store$root_attrs <- attrs
    return(invisible(NULL))
  }
  .store_put(store, .zarr_key(path, "zarr.json"), .utf8_bytes(.json_text(list(
    attributes = attrs, zarr_format = 3L, node_type = "group"
  ))))
  # In the consolidated metadata a group carries its own (empty) inline
  # consolidated metadata, as zarr-python writes it.
  .store_remember(store, path, .json_text(list(
    attributes = attrs,
    zarr_format = 3L,
    consolidated_metadata = list(kind = "inline", must_understand = FALSE, metadata = .json_object()),
    node_type = "group"
  )))
}

.zarr_fill_value <- function(dtype) {
  switch(dtype,
    "|O" = "",
    "|b1" = FALSE,
    "<f4" = 0,
    "<f8" = 0,
    0L
  )
}

# The bytes in one value of a fixed-width dtype; 0 for variable-length text.
.zarr_item_size <- function(dtype) {
  switch(dtype,
    "|O" = 0L,
    "|b1" = 1L,
    "|i1" = 1L,
    "<i2" = 2L,
    "<i4" = 4L,
    "<u4" = 4L,
    "<f4" = 4L,
    "<i8" = 8L,
    "<f8" = 8L
  )
}

.zarr_v3_data_type <- function(dtype) {
  switch(dtype,
    "|O" = "string",
    "|b1" = "bool",
    "|i1" = "int8",
    "<i2" = "int16",
    "<i4" = "int32",
    "<u4" = "uint32",
    "<i8" = "int64",
    "<f4" = "float32",
    "<f8" = "float64"
  )
}

# The bytes a shuffle filter rearranges an array's values by: their size,
# or 0 for no shuffle (one-byte values, text, the expression values, and
# every array of an uncompressed store).
.zarr_shuffle_size <- function(store, dtype, shuffle) {
  item_size <- .zarr_item_size(dtype)
  if (store$compressed && shuffle && item_size > 1L) item_size else 0L
}

.zarr_array_metadata <- function(store, path, dtype, shape, chunks, attrs, shuffle = TRUE) {
  # Every array has attributes, an empty object when it carries none (the
  # arrays inside a sparse matrix), as zarr-python writes them.
  if (is.null(attrs)) attrs <- .json_object()
  item_size <- .zarr_shuffle_size(store, dtype, shuffle)
  if (store$zarr_format == 2L) {
    filters <- if (identical(dtype, "|O")) {
      list(list(id = "vlen-utf8"))
    } else if (item_size > 0L) {
      list(list(elementsize = item_size, id = "shuffle"))
    } else {
      NULL
    }
    .store_put_metadata(store, .zarr_key(path, ".zarray"), list(
      chunks = I(as.numeric(chunks)),
      compressor = if (store$compressed) list(id = "zlib", level = .deflate_level),
      dimension_separator = ".",
      dtype = dtype,
      fill_value = .zarr_fill_value(dtype),
      filters = filters,
      order = "C",
      shape = I(as.numeric(shape)),
      zarr_format = 2L
    ))
    .store_put_metadata(store, .zarr_key(path, ".zattrs"), attrs)
    return(invisible(NULL))
  }
  serializer <- if (identical(dtype, "|O")) {
    list(name = "vlen-utf8", configuration = .json_object())
  } else if (.zarr_item_size(dtype) == 1L) {
    list(name = "bytes")
  } else {
    list(name = "bytes", configuration = list(endian = "little"))
  }
  compressor <- if (!store$compressed) {
    NULL
  } else if (item_size > 0L) {
    list(name = "blosc", configuration = list(
      typesize = item_size, cname = "zlib", clevel = .deflate_level, shuffle = "shuffle", blocksize = 0L
    ))
  } else {
    list(name = "gzip", configuration = list(level = .deflate_level))
  }
  codecs <- c(list(serializer), if (!is.null(compressor)) list(compressor))
  text <- .json_text(list(
    shape = I(as.numeric(shape)),
    data_type = .zarr_v3_data_type(dtype),
    chunk_grid = list(name = "regular", configuration = list(chunk_shape = I(as.numeric(chunks)))),
    chunk_key_encoding = list(name = "default", configuration = list(separator = "/")),
    fill_value = .zarr_fill_value(dtype),
    codecs = codecs,
    attributes = attrs,
    zarr_format = 3L,
    node_type = "array",
    storage_transformers = I(list())
  ))
  .store_put(store, .zarr_key(path, "zarr.json"), .utf8_bytes(text))
  .store_remember(store, path, text)
}

# One chunk's key: format 2 joins the chunk's grid position with "." ("0",
# "3.0"), format 3 puts it under "c/" ("c/0", "c/3/0"; "c" for a scalar).
.zarr_chunk_key <- function(store, position) {
  if (store$zarr_format == 2L) {
    if (length(position) == 0L) "0" else paste(position, collapse = ".")
  } else {
    paste(c("c", position), collapse = "/")
  }
}

# A byte shuffle (numcodecs' Shuffle, Blosc's within a block): the first byte
# of every value, then the second, ...
.shuffle_bytes <- function(bytes, item_size) {
  if (item_size <= 1L || length(bytes) == 0L) return(bytes)
  as.vector(t(matrix(bytes, nrow = item_size)))
}

# A written store's chunks are compressed a batch at a time, every chunk of
# the batch on its own core (src/deflate.c): about this many uncompressed
# bytes, a scheduling choice that keeps a batch small beside the data, not a
# limit on it. Each chunk's bytes are what compressing it alone gives,
# whatever the number of cores.
.deflate_batch_bytes <- 64 * 2^20

# The threads a written store compresses its chunks on: every core, as
# cellucid-python's writer uses every core; at most two under R CMD check,
# as CRAN's policy asks of a check (_R_CHECK_LIMIT_CORES_).
.deflate_threads <- function() {
  threads <- suppressWarnings(parallel::detectCores())
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) || threads < 1) threads <- 1L
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(limit) && !identical(tolower(limit), "false")) threads <- min(threads, 2L)
  as.integer(threads)
}

# DEFLATE in the framing the store's format declares: zlib (format 2) or gzip
# (format 3; header with no name and modification time 0, then the CRC-32 and
# size of the uncompressed bytes, both little-endian). `bytes` is the chunk,
# `zlib` its zlib stream (src/deflate.c, what memCompress(type = "gzip")
# makes).
.framed <- function(store, bytes, zlib) {
  if (store$zarr_format == 2L) return(zlib)
  deflate <- zlib[3:(length(zlib) - 4L)]
  c(
    as.raw(c(0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff)),
    deflate,
    .zip_le(.native_crc32(bytes), 4),
    .zip_le(length(bytes) %% 4294967296, 4)
  )
}

# Whether c-blosc splits a block of `block` bytes into one stream per byte of
# a value, as it does for zlib: values of at most 16 bytes, at least 128 per
# block.
.blosc_splits <- function(item_size, block) {
  item_size <= 16 && block %/% item_size >= 128
}

# c-blosc's automatic block size (compute_blocksize()) for zlib at level 6:
# the whole buffer below 32 KiB, otherwise 512 KiB (32 KiB, doubled for
# zlib, times 8 at level 6), which a split block replaces with 256 KiB values
# held between 64 KiB and 1 MiB; never more than the buffer, and a multiple of
# the value size.
.blosc_block_size <- function(n_bytes, item_size) {
  block <- n_bytes
  if (n_bytes >= 32768) block <- 524288
  if (.blosc_splits(item_size, block)) {
    block <- min(max(min(block, 262144) * item_size, 65536), 1048576)
  }
  block <- min(block, n_bytes)
  if (block > item_size) block <- block %/% item_size * item_size
  block
}

# The Blosc 1 frame of Zarr format 3's blosc codec as the writers declare it
# (cname zlib, clevel 6, byte shuffle, automatic block size), the frame
# c-blosc writes compressing in one thread, as numcodecs does under
# zarr-python -- the same header, blocks and streams, each stream deflated by
# the zlib this package links: a 16-byte header (format 2, zlib format 1,
# flags, value size, then the byte size, block size and frame size as
# little-endian uint32), the start of every block, and each block shuffled and
# cut into one stream per byte of a value (a single stream for a last, shorter
# block), each stream a 4-byte length and its zlib data, or the stream as it
# is when zlib gains nothing. A buffer under 128 bytes, or one whose streams
# would outgrow it, is copied whole after the header instead. The flags: byte
# shuffle (0x01), a whole copy (0x02), no split (0x10), zlib (3 << 5).
#
# In two steps, so the streams of many frames are deflated together: the
# planes to deflate, then `finish(packed)` of their zlib streams, in order.
.blosc_planes <- function(bytes, item_size) {
  n_bytes <- length(bytes)
  block <- .blosc_block_size(n_bytes, item_size)
  split <- .blosc_splits(item_size, block)
  flags <- 0x01 + (if (split) 0 else 0x10) + 0x60
  header <- function(flags, total) {
    c(
      as.raw(c(2L, 1L, flags, item_size)),
      .zip_le(n_bytes, 4), .zip_le(block, 4), .zip_le(total, 4)
    )
  }
  whole_copy <- function() c(header(flags + 0x02, n_bytes + 16), bytes)
  if (n_bytes < 128) return(list(planes = list(), finish = function(packed) whole_copy()))
  n_blocks <- ceiling(n_bytes / block)
  planes <- list()
  for (index in seq_len(n_blocks)) {
    first <- (index - 1) * block
    size <- min(block, n_bytes - first)
    shuffled <- .shuffle_bytes(bytes[first + seq_len(size)], item_size)
    n_streams <- if (split && size == block) item_size else 1L
    each <- size / n_streams
    for (stream in seq_len(n_streams)) planes[[length(planes) + 1L]] <- shuffled[(stream - 1) * each + seq_len(each)]
  }
  finish <- function(packed) {
    limit <- n_bytes + 16
    written <- 16 + 4 * n_blocks
    starts <- numeric(n_blocks)
    body <- vector("list", n_blocks)
    at <- 0L
    for (index in seq_len(n_blocks)) {
      starts[[index]] <- written
      first <- (index - 1) * block
      size <- min(block, n_bytes - first)
      n_streams <- if (split && size == block) item_size else 1L
      each <- size / n_streams
      streams <- vector("list", n_streams)
      for (stream in seq_len(n_streams)) {
        at <- at + 1L
        written <- written + 4
        room <- min(each, limit - written)
        if (room <= 0) return(whole_copy())
        stored <- packed[[at]]
        if (length(stored) > room || length(stored) == each) {
          if (written + each > limit) return(whole_copy())
          stored <- planes[[at]]
        }
        streams[[stream]] <- c(.zip_le(length(stored), 4), stored)
        written <- written + length(stored)
      }
      body[[index]] <- unlist(streams, use.names = FALSE)
    }
    c(
      header(flags, written),
      unlist(lapply(starts, .zip_le, bytes = 4), use.names = FALSE),
      unlist(body, use.names = FALSE)
    )
  }
  list(planes = planes, finish = finish)
}

# One chunk: `values` is a function returning the chunk's little-endian
# bytes. A written store calls it now and holds its bytes until a batch of
# them is compressed (.store_flush()); a memory store keeps the function,
# marked `transform` when its bytes are made by more than a slice of a vector
# (they are then kept in the store's cache, .store_get()).
.zarr_put_chunk <- function(store, path, position, values, dtype, shuffle = TRUE, transform = FALSE) {
  key <- .zarr_key(path, .zarr_chunk_key(store, position))
  if (!store$compressed) {
    if (transform) attr(values, "transform") <- TRUE
    return(.store_put(store, key, values))
  }
  item_size <- .zarr_shuffle_size(store, dtype, shuffle)
  bytes <- values()
  pending <- if (store$zarr_format == 3L && item_size > 0L) {
    .blosc_planes(bytes, item_size)
  } else {
    shuffled <- .shuffle_bytes(bytes, item_size)
    list(planes = list(shuffled), finish = function(packed) .framed(store, shuffled, packed[[1L]]))
  }
  pending$key <- key
  store$pending[[length(store$pending) + 1L]] <- pending
  store$pending_bytes <- store$pending_bytes + length(bytes)
  if (store$pending_bytes >= .deflate_batch_bytes) .store_flush(store)
  invisible(NULL)
}

# Compress the chunks held (.zarr_put_chunk()), every plane on its own core,
# and put each in the order it was made.
.store_flush <- function(store) {
  pending <- store$pending
  if (length(pending) == 0L) return(invisible(NULL))
  store$pending <- list()
  store$pending_bytes <- 0
  planes <- unlist(lapply(pending, `[[`, "planes"), recursive = FALSE)
  packed <- .native_deflate_chunks(planes, .deflate_level, .deflate_threads())
  at <- 0L
  for (entry in pending) {
    count <- length(entry$planes)
    .store_put(store, entry$key, entry$finish(packed[at + seq_len(count)]))
    at <- at + count
  }
  invisible(NULL)
}

# --- Values computed when a chunk needs them ---------------------------------

# A vector whose values are computed only when a chunk is made, once, and
# then kept: the gene pointers of an R object's genes x cells matrix
# (expression-input.R). Its length is known without computing it, so the
# store's metadata is planned first.
.deferred <- function(compute, length) {
  memo <- new.env(parent = emptyenv())
  structure(list(compute = compute, length = length, memo = memo), class = "cellucid_deferred")
}

.is_deferred <- function(values) inherits(values, "cellucid_deferred")

# A sparse matrix's gene-major array (its data or its indices) made block by
# block from a genes x cells matrix (expression-input.R): `chunk(first,
# count, chunk_length, chunk_first =, chunk_last =)` returns positions
# first .. first + count - 1 padded to `chunk_length` values, of the chunk
# holding positions chunk_first .. chunk_last (a source that makes a chunk
# whole keeps it for the chunk's other ranges).
.is_chunked <- function(values) inherits(values, "cellucid_chunked")

.value_length <- function(values) {
  if (.is_deferred(values) || .is_chunked(values)) values$length else length(values)
}

.value_dim <- function(values) {
  if (inherits(values, "cellucid_rows_in_order")) return(.value_dim(values$matrix))
  if (.is_transposed(values)) rev(dim(values$matrix)) else dim(values)
}

# A matrix read as its transpose without making it: a dense genes x cells
# matrix presented as cells x genes, each chunk's block taken from the
# matrix as it is and transposed alone, so no transposed copy of the whole
# matrix is ever held.
.transposed <- function(matrix) structure(list(matrix = matrix), class = "cellucid_transposed")

.is_transposed <- function(values) inherits(values, "cellucid_transposed")

# A matrix (or a .transposed() one) read with its rows in the order `order`
# gives (1-based): row r of it is row order[r] of the matrix, each block taken
# through the order, never a reordered copy.
.rows_in_order <- function(matrix, order) structure(list(matrix = matrix, order = order), class = "cellucid_rows_in_order")

# Rows `rows` and columns `cols` of a matrix, a .transposed() one, or one
# read through .rows_in_order().
.matrix_block <- function(matrix, rows, cols) {
  if (inherits(matrix, "cellucid_rows_in_order")) {
    return(.matrix_block(matrix$matrix, matrix$order[rows], cols))
  }
  if (.is_transposed(matrix)) {
    return(t(matrix$matrix[cols, rows, drop = FALSE]))
  }
  matrix[rows, cols, drop = FALSE]
}

.realize <- function(values) {
  if (!.is_deferred(values)) {
    return(values)
  }
  memo <- values$memo
  if (is.null(memo$value)) memo$value <- values$compute()
  memo$value
}

# --- Element encodings ------------------------------------------------------

# Little-endian bytes of `values` (numbers or logicals) in one fixed-width
# dtype, one pass in C (src/chunks.c). Integer dtypes take whole numbers;
# "<i8" is exact for every value of magnitude below 2^53. Every NaN and NA
# of a floating dtype is the one canonical quiet NaN numpy writes
# (0x7FC00000, 0x7FF8000000000000), so the bytes do not depend on the
# platform (R's NA is a NaN with a payload, and its NaN has the sign bit set
# on x86-64).
.encode_values <- function(values, dtype) {
  .native_encode_values(values, 1, length(values), length(values), dtype)
}

# numcodecs' VLenUTF8 (Zarr format 3: the vlen-utf8 codec, the same bytes):
# the item count, then each item's byte length and bytes, every count a
# little-endian uint32.
.encode_vlen_utf8 <- function(strings) {
  strings <- enc2utf8(as.character(strings))
  lengths <- nchar(strings, type = "bytes")
  total <- 4 + 4 * length(strings) + sum(lengths)
  out <- raw(total)
  out[1:4] <- writeBin(length(strings), raw(), size = 4L, endian = "little")
  if (length(strings) == 0L) {
    return(out)
  }
  starts <- 5 + cumsum(c(0, utils::head(4 + lengths, -1L)))
  header_positions <- as.vector(outer(0:3, starts, "+"))
  out[header_positions] <- writeBin(as.integer(lengths), raw(), size = 4L, endian = "little")
  if (sum(lengths) > 0) {
    is_text <- rep(TRUE, total)
    is_text[1:4] <- FALSE
    is_text[header_positions] <- FALSE
    out[is_text] <- charToRaw(paste(strings, collapse = ""))
  }
  out
}

.chunk_length <- function(n, target) {
  max(1, min(n, target))
}

# A one-dimensional array. Every chunk is stored at the full chunk length,
# the last one padded with the fill value, as both Zarr formats require.
# `values` is a vector, a .deferred() one, or a gene-major array
# (.is_chunked()). A chunk of numbers is one slice of the vector, made in one
# pass; a chunk of text is encoded (a transform, kept in a memory store's
# cache).
.zarr_vector <- function(store, path, values, dtype, attrs = NULL,
                         target = .chunk_elements_whole, shuffle = TRUE) {
  n <- .value_length(values)
  chunk <- .chunk_length(n, target)
  .zarr_array_metadata(store, path, dtype, n, chunk, attrs, shuffle = shuffle)
  if (n == 0) {
    return(invisible(NULL))
  }
  text <- identical(dtype, "|O")
  item <- if (text) NULL else .zarr_item_size(dtype)
  chunk_values <- function(first, last) {
    force(first)
    force(last)
    if (text) {
      return(function() {
        part <- .realize(values)[first:last]
        .encode_vlen_utf8(c(part, rep("", chunk - length(part))))
      })
    }
    # Values [from, to) of the chunk: those of the vector, then the padding.
    .ranged(function(from = 0, to = chunk) {
      count <- max(0, min(last, first + to - 1) - (first + from) + 1)
      if (count == 0) return(raw((to - from) * item))
      # The chunk's own positions too: a source that makes a chunk whole keeps it for its other ranges.
      if (.is_chunked(values)) return(values$chunk(first + from, count, to - from, chunk_first = first, chunk_last = last))
      .native_encode_values(.realize(values), first + from, count, to - from, dtype)
    }, chunk, item)
  }
  n_chunks <- ceiling(n / chunk)
  for (index in seq_len(n_chunks)) {
    first <- (index - 1) * chunk + 1
    last <- min(n, index * chunk)
    .zarr_put_chunk(store, path, index - 1, chunk_values(first, last), dtype, shuffle = shuffle,
                    transform = text)
  }
  invisible(NULL)
}

# A zero-dimensional array: one value in one chunk.
.zarr_scalar <- function(store, path, value, dtype, attrs) {
  .zarr_array_metadata(store, path, dtype, integer(0), integer(0), attrs)
  bytes <- if (identical(dtype, "|O")) {
    .encode_vlen_utf8(value)
  } else {
    .encode_values(value, dtype)
  }
  .zarr_put_chunk(store, path, integer(0), function() bytes, dtype)
}

# A two-dimensional array of one floating dtype ("<f4" or "<f8") from a
# numeric matrix (or a .transposed() one), in chunks of `chunk_rows` x
# `chunk_cols`, each chunk in C (row-major) order and padded to the full chunk
# shape at the edges. Making a chunk is a transform (R's matrices are
# column-major), so a memory store keeps the chunks it made in its cache.
.zarr_matrix <- function(store, path, matrix, dtype, chunk_rows, chunk_cols, attrs = NULL,
                         shuffle = TRUE) {
  shape <- .value_dim(matrix)
  n_rows <- shape[[1L]]
  n_cols <- shape[[2L]]
  .zarr_array_metadata(store, path, dtype, c(n_rows, n_cols), c(chunk_rows, chunk_cols), attrs, shuffle = shuffle)
  if (n_rows == 0 || n_cols == 0) {
    return(invisible(NULL))
  }
  chunk_values <- function(rows, cols) {
    force(rows)
    force(cols)
    function() {
      block <- .matrix_block(matrix, rows, cols)
      if (length(rows) < chunk_rows || length(cols) < chunk_cols) {
        padded <- matrix(0, nrow = chunk_rows, ncol = chunk_cols)
        padded[seq_along(rows), seq_along(cols)] <- block
        block <- padded
      }
      .encode_values(as.vector(t(block)), dtype)
    }
  }
  for (row_index in seq_len(ceiling(n_rows / chunk_rows))) {
    rows <- ((row_index - 1) * chunk_rows + 1):min(n_rows, row_index * chunk_rows)
    for (col_index in seq_len(ceiling(n_cols / chunk_cols))) {
      cols <- ((col_index - 1) * chunk_cols + 1):min(n_cols, col_index * chunk_cols)
      .zarr_put_chunk(
        store,
        path,
        c(row_index - 1, col_index - 1),
        chunk_values(rows, cols),
        dtype,
        shuffle = shuffle,
        transform = TRUE
      )
    }
  }
  invisible(NULL)
}
