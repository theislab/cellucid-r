# One input must produce one set of bytes on every machine, so a compressed
# payload carries a fixed ten-byte gzip header rather than the one the
# platform's zlib would stamp on it.

# RFC 1952 section 2.3.1: `XFL` records which end of the speed/ratio range the
# compressor was asked for. zlib's own gzip writer and `gzip.GzipFile` in
# cellucid-python both write 2 for the slowest level and 4 for the fastest;
# R's `gzfile()` hardcodes 0 for every level, which is a third spelling of the
# same fact.
.gzip_extra_flags <- function(compression) {
  if (compression >= 9L) {
    return(as.raw(0x02))
  }
  if (compression <= 1L) {
    return(as.raw(0x04))
  }
  as.raw(0x00)
}

# The ten bytes `gzip.GzipFile(filename = "", mtime = 0)` produces in
# cellucid-python: the magic, deflate, no optional header fields, a fixed
# Unix-epoch timestamp, the level's extra flags, and `OS = 0xff`, the
# "unknown" code. `0xff` is the only `OS` value that does not vary with the
# machine, so it is what makes one input produce one set of bytes everywhere.
.canonical_gzip_header <- function(compression) {
  c(
    as.raw(c(0x1f, 0x8b, 0x08, 0x00)),
    as.raw(c(0x00, 0x00, 0x00, 0x00)),
    .gzip_extra_flags(compression),
    as.raw(0xff)
  )
}

# `gzfile()` gives no control over the header it writes, so the header is
# replaced in place afterwards. It is a fixed ten bytes with no optional
# fields, and only the `OS` and `XFL` bytes carry a value this writer needs to
# change; the deflate stream and the CRC32/ISIZE trailer zlib produced are
# left exactly as written. Anything other than that fixed form means the
# member was not produced by the writer this rewrite was verified against, and
# the export stops rather than stamping a header onto bytes it cannot read.
.write_canonical_gzip_header <- function(path, compression) {
  con <- file(path, open = "r+b")
  on.exit(close(con), add = TRUE)
  header <- readBin(con, what = "raw", n = 10L)
  if (
    length(header) != 10L ||
      !identical(header[1:3], as.raw(c(0x1f, 0x8b, 0x08))) ||
      !identical(header[[4L]], as.raw(0x00))
  ) {
    stop(
      "Compressed payload ",
      path,
      " does not begin with the fixed ten-byte gzip header, so its header ",
      "cannot be made platform independent.",
      call. = FALSE
    )
  }
  seek(con, where = 0L, origin = "start")
  writeBin(.canonical_gzip_header(compression), con)
  invisible(path)
}
