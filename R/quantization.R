# Continuous values narrowed to an integer code range. The reader
# reconstructs a value from the code and the bounds published beside it, so
# the bounds are measured in the viewer's float32 domain, and the constant
# field is a named case rather than a division by a zero-width range.

.float32_quantization_chunk_size <- 1048576L

.roundtrip_finite_float32_chunk <- function(values, field_name) {
  encoded <- writeBin(
    as.double(values),
    raw(),
    size = 4L,
    endian = "little"
  )
  rounded <- readBin(
    encoded,
    what = double(),
    n = length(values),
    size = 4L,
    endian = "little"
  )
  if (
    length(rounded) != length(values) ||
      anyNA(rounded) ||
      any(!is.finite(rounded))
  ) {
    stop(
      "Cannot quantize field '",
      field_name,
      "': values must remain finite in the viewer's float32 domain.",
      call. = FALSE
    )
  }
  # Underflow is the other end of the same range and this path was missing it:
  # a nonzero value below the smallest float32 subnormal round-trips to 0,
  # which is finite, so the check above accepted it and the field was published
  # with that value silently replaced by zero. The rule itself lives in
  # `.outside_finite_float32()` beside the writer that already enforced it --
  # a second copy of the rule is exactly how this path drifted -- and only the
  # message, which names the field, belongs here.
  if (any(.outside_finite_float32(values), na.rm = TRUE)) {
    stop(
      "Cannot quantize field '",
      field_name,
      "': values must remain finite in the viewer's float32 domain.",
      call. = FALSE
    )
  }
  rounded
}

# Name compact_v1's constant-field case from one payload's bounds.
#
# A gene detected in no published cell, or an obs column a subset flattened, is
# ordinary scientific data, so the format encodes it rather than rejecting it.
# The case is declared by equal bounds -- the entry stays exactly
# [index, key, minValue, maxValue], so neither its shape nor its length moves --
# and every code is written as 0. Writer and reader then take the same named
# branch instead of the general arithmetic: the writer never divides by
# (maxValue - minValue), and the reader fills minValue straight through. The
# value returns bit-exact, not within a quantization step.
#
# Sole derivation point on this writer, and the exact peer of
# _is_constant_continuous_range() in the Python exporter.
.is_constant_continuous_range <- function(min_val, max_val) {
  min_val == max_val
}

.quantize_continuous <- function(
    values,
    bits,
    field_name,
    .chunk_size = .float32_quantization_chunk_size
) {
  bits <- .validate_quantization_bits(bits, arg = "bits")
  field_name <- .validate_required_string(field_name, "field_name")
  .chunk_size <- .validate_positive_integer(
    .chunk_size,
    ".chunk_size"
  )
  if (
    !is.numeric(values) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop(
      "Cannot quantize field '",
      field_name,
      "': values must be a native numeric vector.",
      call. = FALSE
    )
  }
  if (length(values) == 0L) {
    stop(
      "Cannot quantize field '",
      field_name,
      "': it contains no finite values.",
      call. = FALSE
    )
  }

  min_val <- Inf
  max_val <- -Inf
  starts <- seq.int(1L, length(values), by = .chunk_size)
  for (start in starts) {
    end <- min(
      length(values),
      start + as.double(.chunk_size) - 1
    )
    source_chunk <- values[start:end]
    if (anyNA(source_chunk) || any(!is.finite(source_chunk))) {
      stop(
        "Cannot quantize field '",
        field_name,
        "': values must contain only finite values.",
        call. = FALSE
      )
    }
    rounded <- .roundtrip_finite_float32_chunk(
      source_chunk,
      field_name
    )
    min_val <- min(min_val, min(rounded))
    max_val <- max(max_val, max(rounded))
  }

  if (bits == 8L) {
    max_quant <- 254L
  } else {
    max_quant <- 65534L
  }

  if (.is_constant_continuous_range(min_val, max_val)) {
    # The constant case, taken before any scale exists. The general path
    # divides by (max_val - min_val), which is exactly zero here, so a constant
    # field must never reach it. This also covers native-double variation that
    # collapses to one value in the viewer's float32 domain: one float32 value
    # is one constant. Mirrors _quantize_continuous() in the Python exporter.
    return(list(
      quantized = integer(length(values)),
      min_val = min_val,
      max_val = max_val,
      scale = 0
    ))
  }

  scale <- max_quant / (max_val - min_val)
  quantized <- integer(length(values))
  for (start in starts) {
    end <- min(
      length(values),
      start + as.double(.chunk_size) - 1
    )
    rounded <- .roundtrip_finite_float32_chunk(
      values[start:end],
      field_name
    )
    normalized <- (rounded - min_val) * scale
    normalized[normalized < 0] <- 0
    normalized[normalized > max_quant] <- max_quant
    quantized[start:end] <- as.integer(
      floor(normalized + 1e-8)
    )
  }

  list(
    quantized = quantized,
    min_val = min_val,
    max_val = max_val,
    scale = scale
  )
}

.quantize_nullable_outlier_quantiles <- function(values, bits, field_name) {
  bits <- .validate_quantization_bits(bits, arg = "bits")
  field_name <- .validate_required_string(field_name, "field_name")
  if (
    !is.numeric(values) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop(
      "Cannot quantize outlier field '",
      field_name,
      "': values must be a native numeric vector.",
      call. = FALSE
    )
  }
  if (any(is.na(values) & !is.nan(values))) {
    stop(
      "Cannot quantize outlier field '",
      field_name,
      "': only generated NaN values may represent missing quantiles.",
      call. = FALSE
    )
  }
  if (any(is.infinite(values))) {
    stop(
      "Cannot quantize outlier field '",
      field_name,
      "': values must contain only finite values or generated NaN.",
      call. = FALSE
    )
  }

  missing_mask <- is.nan(values)
  finite_quantized <- .quantize_continuous(
    values[!missing_mask],
    bits = bits,
    field_name = field_name
  )
  missing_value <- if (bits == 8L) 255L else 65535L
  quantized <- rep(missing_value, length(values))
  quantized[!missing_mask] <- finite_quantized$quantized

  list(
    quantized = quantized,
    min_val = finite_quantized$min_val,
    max_val = finite_quantized$max_val,
    scale = finite_quantized$scale
  )
}
