# Every scientific payload byte the viewer reads is written here, in the
# little-endian widths the browser contract names, through one function that
# also makes a compressed member's header platform independent.
#
# Byte order is guaranteed twice, in two different ways, and neither is the
# host's. The float widths pass `endian = "little"` to `writeBin()`, which is
# what stops a big-endian machine from publishing a payload the viewer would
# read backwards. The unsigned widths never let `writeBin()` see a number at
# all: `.pack_uint16()` and `.pack_uint32()` lay the bytes out themselves,
# least significant first, and hand `writeBin()` a `raw` vector, which has no
# byte order to get wrong. `endian =` on a `raw` write does nothing, so it is
# not written there -- an argument that looks like the guarantee but is inert
# would hide where the guarantee actually is.

# The two ends of the viewer's float32 domain, named once.
#
# Overflow is the loud end: a value past the largest float32 becomes Inf, and
# every finiteness check in the package already catches it. Underflow is the
# quiet one: a nonzero value below the smallest subnormal becomes 0, which is
# finite, so a finiteness-only check publishes a payload whose value has
# silently been replaced by zero. Both ends are one rule, so both live in one
# predicate -- a second copy of it is how the quantized path came to accept
# what the unquantized path refused.

.float32_max_magnitude <- (2 - 2^-23) * 2^127
.float32_min_subnormal <- 2^-149

.outside_finite_float32 <- function(values) {
  abs(values) > .float32_max_magnitude |
    (values != 0 & abs(values) < .float32_min_subnormal)
}

# Gate an input the exporter transforms before writing.
#
# `.write_float32_vector()` measures the bytes it is about to write, which is
# the right place for everything published unchanged. It is the wrong place for
# an input that is normalized first: embedding coordinates are centred and
# scaled into roughly [-1, 1], so a source coordinate of 1e300 or 1e-320
# arrives at the writer already inside the float32 range and is published as a
# number the caller never supplied. cellucid-python refuses both at the source
# in `_require_finite_embedding_source()`, and one export format cannot depend
# on which writer wrote it.
.require_finite_float32_source <- function(values, label) {
  if (any(.outside_finite_float32(values), na.rm = TRUE)) {
    stop(
      label,
      " contains values outside the finite float32 range: every nonzero ",
      "value must be representable between 2^-149 and (2 - 2^-23) * 2^127 ",
      "in absolute magnitude.",
      call. = FALSE
    )
  }
  invisible(values)
}

.write_float32_matrix_row_major <- function(path, mat, compression = NULL) {
  .write_float32_vector(path, as.vector(t(mat)), compression = compression)
}

.write_float32_vector <- function(
    path,
    values,
    compression = NULL,
    allow_generated_nan = FALSE
) {
  if (
    !is.logical(allow_generated_nan) ||
      length(allow_generated_nan) != 1L ||
      is.na(allow_generated_nan)
  ) {
    stop("allow_generated_nan must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.numeric(values) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop(
      "float32 values must be a native numeric vector.",
      call. = FALSE
    )
  }
  has_forbidden_missing <- any(is.na(values) & !is.nan(values))
  if (
    has_forbidden_missing ||
      any(is.infinite(values)) ||
      (!allow_generated_nan && any(is.nan(values)))
  ) {
    stop(
      "float32 values must contain only finite float32 values; generated ",
      "NaN is accepted only when allow_generated_nan is TRUE.",
      call. = FALSE
    )
  }

  if (any(.outside_finite_float32(values), na.rm = TRUE)) {
    stop(
      "Every nonzero finite float32 value must be representable between ",
      "2^-149 and (2 - 2^-23) * 2^127 in absolute magnitude.",
      call. = FALSE
    )
  }

  .write_binary_payload(path, compression, function(con) {
    writeBin(as.double(values), con, size = 4L, endian = "little")
  })
}

.write_float64_vector <- function(path, values, compression = NULL) {
  if (
    !is.double(values) ||
      is.object(values) ||
      !is.null(dim(values)) ||
      anyNA(values) ||
      any(!is.finite(values))
  ) {
    stop(
      "float64 values must be a finite native double vector.",
      call. = FALSE
    )
  }
  .write_binary_payload(path, compression, function(con) {
    writeBin(values, con, size = 8L, endian = "little")
  })
}

.write_uint8 <- function(path, values, compression = NULL) {
  if (
    !(is.integer(values) || is.double(values)) ||
      is.object(values) ||
      !is.null(dim(values)) ||
      anyNA(values) ||
      any(!is.finite(values)) ||
      any(values != floor(values)) ||
      any(values < 0 | values > 255)
  ) {
    stop(
      "uint8 values must be a numeric vector of finite integers in the ",
      "inclusive range 0 to 255.",
      call. = FALSE
    )
  }
  # One byte per value, so there is no byte order to choose.
  raw_vals <- as.raw(as.integer(values))
  .write_binary_payload(path, compression, function(con) {
    writeBin(raw_vals, con)
  })
}

# The uint16 byte order, in one place. Two bytes are laid out least
# significant first by arithmetic on the value, so the result does not depend
# on the byte order of the machine that ran it: 258 is `0x02 0x01` on every
# host. `.write_uint16()` then only has to hand these bytes over unaltered.
.pack_uint16 <- function(values) {
  if (
    !(is.integer(values) || is.double(values)) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop("uint16 values must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(values) || any(!is.finite(values))) {
    stop("uint16 values must all be finite.", call. = FALSE)
  }
  if (any(values != floor(values))) {
    stop("uint16 values must all be integers.", call. = FALSE)
  }
  if (any(values < 0 | values > 65535)) {
    stop(
      "uint16 values must all be in the inclusive range 0 to 65,535.",
      call. = FALSE
    )
  }
  integer_values <- as.integer(values)
  lo <- as.raw(bitwAnd(integer_values, 0xFFL))
  hi <- as.raw(bitwAnd(bitwShiftR(integer_values, 8L), 0xFFL))
  as.raw(rbind(lo, hi))
}

.write_uint16 <- function(path, values, compression = NULL) {
  raw_vals <- .pack_uint16(values)
  .write_binary_payload(path, compression, function(con) {
    writeBin(raw_vals, con)
  })
}

# The uint32 byte order, in one place, on the same terms as `.pack_uint16()`:
# four bytes, least significant first, from arithmetic rather than from the
# host's representation.
.pack_uint32 <- function(values) {
  if (
    !(is.integer(values) || is.double(values)) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop("uint32 values must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(values) || any(!is.finite(values))) {
    stop("uint32 values must all be finite.", call. = FALSE)
  }
  if (any(values != floor(values))) {
    stop("uint32 values must all be integers.", call. = FALSE)
  }
  if (any(values < 0 | values > 4294967295)) {
    stop(
      "uint32 values must all be in the inclusive range 0 to 4,294,967,295.",
      call. = FALSE
    )
  }
  b0 <- as.raw(floor(values) %% 256)
  b1 <- as.raw(floor(values / 256) %% 256)
  b2 <- as.raw(floor(values / 65536) %% 256)
  b3 <- as.raw(floor(values / 16777216) %% 256)
  as.raw(rbind(b0, b1, b2, b3))
}

.write_uint32 <- function(path, values, compression = NULL) {
  raw_vals <- .pack_uint32(values)
  .write_binary_payload(path, compression, function(con) {
    writeBin(raw_vals, con)
  })
}

# Compression is visible in the export as one suffix on every payload path, so
# the file this writer opens and the path patterns the manifests publish put
# the same characters in the same place.
.gz_suffix <- function(compression) {
  if (!is.null(compression)) ".gz" else ""
}

.final_binary_path <- function(path, compression) {
  suffix <- .gz_suffix(compression)
  if (nzchar(suffix)) paste0(path, suffix) else path
}

.binary_connection <- function(path, compression = NULL) {
  final <- .final_binary_path(path, compression)
  if (!is.null(compression)) {
    return(gzfile(final, open = "wb", compression = compression))
  }
  file(final, open = "wb")
}

# Every scientific payload is written here, so a compressed one cannot reach
# the export without its header being made platform independent. The payload
# is still streamed through the connection: nothing holds the compressed bytes
# in memory, and `gzfile()` still chooses the deflate encoding for the level
# the caller asked for.
.write_binary_payload <- function(path, compression, write_payload) {
  final <- .final_binary_path(path, compression)
  con <- .binary_connection(path, compression = compression)
  open_connection <- TRUE
  on.exit(if (open_connection) close(con), add = TRUE)
  write_payload(con)
  # The gzip member is only complete, and its header only on disk, once the
  # connection is closed, so the header is rewritten after the close rather
  # than from the deferred one.
  close(con)
  open_connection <- FALSE
  if (!is.null(compression)) {
    .write_canonical_gzip_header(final, compression)
  }
  invisible(final)
}
