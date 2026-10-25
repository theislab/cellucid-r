# The .zarr.zip container: one ZIP archive whose top level is the Zarr
# store, every entry stored (method 0, not compressed again), so a reader can
# fetch any chunk with one range request at the offset the central directory
# gives. The bytes are fully determined by the entries: names are UTF-8
# (general-purpose flag bit 11), every timestamp is 1980-01-01 00:00, and the
# permissions are those of a 0644 regular file. ZIP64 records are written only
# where a size, an offset, or the entry count exceeds the classic fields.

.zip_u32_max <- 4294967295
.zip_u16_max <- 65535
.zip_dos_date_1980_01_01 <- 33
.zip_utf8_flag <- 2048
.zip_external_attributes <- 33188 * 65536 # a 0644 regular file (0o100644 << 16)

# `value` as `bytes` little-endian bytes; exact for every integer below 2^53.
.zip_le <- function(value, bytes) {
  value <- as.numeric(value)
  out <- raw(bytes)
  for (index in seq_len(bytes)) {
    out[[index]] <- as.raw(value %% 256)
    value <- value %/% 256
  }
  out
}

.zip_local_header <- function(name, crc, size) {
  zip64 <- size >= .zip_u32_max
  extra <- if (zip64) {
    c(.zip_le(1, 2), .zip_le(16, 2), .zip_le(size, 8), .zip_le(size, 8))
  } else {
    raw(0)
  }
  classic_size <- if (zip64) .zip_u32_max else size
  c(
    .zip_le(67324752, 4), # 0x04034b50
    .zip_le(if (zip64) 45 else 20, 2),
    .zip_le(.zip_utf8_flag, 2),
    .zip_le(0, 2),
    .zip_le(0, 2),
    .zip_le(.zip_dos_date_1980_01_01, 2),
    .zip_le(crc, 4),
    .zip_le(classic_size, 4),
    .zip_le(classic_size, 4),
    .zip_le(length(name), 2),
    .zip_le(length(extra), 2),
    name,
    extra
  )
}

.zip_central_record <- function(name, crc, size, offset) {
  fields <- numeric(0)
  if (size >= .zip_u32_max) fields <- c(fields, size, size)
  if (offset >= .zip_u32_max) fields <- c(fields, offset)
  extra <- if (length(fields) > 0L) {
    c(
      .zip_le(1, 2),
      .zip_le(8 * length(fields), 2),
      unlist(lapply(fields, .zip_le, bytes = 8), use.names = FALSE)
    )
  } else {
    raw(0)
  }
  version <- if (length(fields) > 0L) 45 else 20
  c(
    .zip_le(33639248, 4), # 0x02014b50
    .zip_le(3 * 256 + version, 2),
    .zip_le(version, 2),
    .zip_le(.zip_utf8_flag, 2),
    .zip_le(0, 2),
    .zip_le(0, 2),
    .zip_le(.zip_dos_date_1980_01_01, 2),
    .zip_le(crc, 4),
    .zip_le(min(size, .zip_u32_max), 4),
    .zip_le(min(size, .zip_u32_max), 4),
    .zip_le(length(name), 2),
    .zip_le(length(extra), 2),
    .zip_le(0, 2),
    .zip_le(0, 2),
    .zip_le(0, 2),
    .zip_le(.zip_external_attributes, 4),
    .zip_le(min(offset, .zip_u32_max), 4),
    name,
    extra
  )
}

.zip_end_records <- function(count, cd_size, cd_offset) {
  zip64 <- count >= .zip_u16_max || cd_size >= .zip_u32_max ||
    cd_offset >= .zip_u32_max
  out <- raw(0)
  if (zip64) {
    out <- c(
      .zip_le(101075792, 4), # 0x06064b50
      .zip_le(44, 8),
      .zip_le(3 * 256 + 45, 2),
      .zip_le(45, 2),
      .zip_le(0, 4),
      .zip_le(0, 4),
      .zip_le(count, 8),
      .zip_le(count, 8),
      .zip_le(cd_size, 8),
      .zip_le(cd_offset, 8),
      .zip_le(117853008, 4), # 0x07064b50
      .zip_le(0, 4),
      .zip_le(cd_offset + cd_size, 8),
      .zip_le(1, 4)
    )
  }
  classic_count <- min(count, .zip_u16_max)
  c(
    out,
    .zip_le(101010256, 4), # 0x06054b50
    .zip_le(0, 2),
    .zip_le(0, 2),
    .zip_le(classic_count, 2),
    .zip_le(classic_count, 2),
    .zip_le(min(cd_size, .zip_u32_max), 4),
    .zip_le(min(cd_offset, .zip_u32_max), 4),
    .zip_le(0, 2)
  )
}

.zip_add_entry <- function(store, key, bytes) {
  name <- .utf8_bytes(key)
  size <- length(bytes)
  crc <- .native_crc32(bytes)
  header <- .zip_local_header(name, crc, size)
  writeBin(header, store$connection)
  writeBin(bytes, store$connection)
  store$central[[length(store$central) + 1L]] <-
    .zip_central_record(name, crc, size, store$offset)
  store$offset <- store$offset + length(header) + size
  invisible(NULL)
}

.zip_finish <- function(store) {
  central <- unlist(store$central, use.names = FALSE)
  if (is.null(central)) central <- raw(0)
  writeBin(central, store$connection)
  writeBin(
    .zip_end_records(length(store$central), length(central), store$offset),
    store$connection
  )
  invisible(NULL)
}
