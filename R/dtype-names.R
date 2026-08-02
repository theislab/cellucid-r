# One storage width, spelled three ways: the number of bits the caller asks
# for, the dtype the manifest publishes beside a payload, and the extension the
# payload file itself carries. All three spellings are part of the export
# format the viewer reads, and the schema that declares a path pattern and the
# check that re-expands it against the directory must agree on them, so they
# are derived here and nowhere else.

.quantization_dtype <- function(bits) {
  if (bits == 8L) "uint8" else "uint16"
}

.quantization_ext <- function(bits) {
  if (bits == 8L) "u8" else "u16"
}

.codes_extension_by_dtype <- c(uint8 = "u8", uint16 = "u16")
