# Normalizing one embedding into the viewer's coordinate range. Each exported
# embedding is centred and scaled on its own, with one factor shared by its
# axes so its aspect ratio survives; the factor is kept because the vector
# fields drawn on top of it have to use the same one.

.normalize_embedding <- function(arr) {
  if (
    !is.matrix(arr) ||
      !is.numeric(arr) ||
      anyNA(arr) ||
      any(!is.finite(arr))
  ) {
    stop(
      "Embedding coordinates must be a finite numeric matrix.",
      call. = FALSE
    )
  }
  axis_mins <- apply(arr, 2, min)
  axis_maxs <- apply(arr, 2, max)
  axis_ranges <- axis_maxs - axis_mins
  max_range <- max(axis_ranges)
  if (!is.finite(max_range) || max_range <= 0) {
    stop(
      "Embedding coordinates must span a positive finite coordinate range.",
      call. = FALSE
    )
  }
  center <- (axis_mins + axis_maxs) / 2
  scale_factor <- 2 / max_range

  centered <- sweep(arr, 2, center, FUN = "-")
  coords <- centered * scale_factor

  list(
    coords = coords,
    info = list(
      original_range = as.numeric(max_range),
      center = as.numeric(center),
      scale_factor = as.numeric(scale_factor)
    )
  )
}
