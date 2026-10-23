# The complete surface of src/ (registered in src/init.c): every .Call() of
# the package, so the R side of the native contract is readable in one screen.

# src/crc32.c: the CRC-32 of a raw vector, continuing from `previous` (the CRC
# of the bytes before it, 0 for none). The ZIP writer in zip.R needs it, and
# so does the gzip framing of a Zarr format 3 chunk (zarr-store.R).
.native_crc32 <- function(bytes, previous = 0) {
  .Call(C_cellucid_crc32, bytes, as.double(previous))
}

# src/latent_outliers.c: one categorical obs column's outlier quantiles in the
# latent space (R/outlier-quantiles.R).
.native_latent_outlier_quantiles <- function(latent, codes, n_categories, min_points) {
  .Call(C_cellucid_latent_outlier_quantiles, latent, codes, n_categories, min_points)
}
