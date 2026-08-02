# The coordinates the viewer draws. They arrive per dimension, are checked
# against each other for one shared cell axis, and are centered and scaled
# into the viewer's coordinate range; the scale is kept because the vector
# fields drawn on top of them have to use the same one.

.collect_embeddings <- function(embeddings_in) {
  embeddings <- list()
  for (dim_chr in names(embeddings_in)) {
    x <- embeddings_in[[dim_chr]]
    if (is.null(x)) {
      next
    }
    dim_int <- as.integer(dim_chr)
    if (
      is.numeric(x) &&
        !is.object(x) &&
        is.null(dim(x)) &&
        dim_int == 1L
    ) {
      x <- matrix(x, ncol = 1L)
    } else if (.is_numeric_data_frame(x)) {
      x <- as.matrix(x)
    } else if (inherits(x, "Matrix")) {
      x <- as.matrix(x)
    } else if (
      !is.matrix(x) ||
        !is.numeric(x) ||
        is.object(x)
    ) {
      stop(
        "X_umap_", dim_chr,
        "d must be a finite real numeric matrix.",
        call. = FALSE
      )
    } else {
      x <- unclass(x)
    }
    storage.mode(x) <- "double"
    embeddings[[dim_chr]] <- x
  }
  if (length(embeddings) == 0) {
    stop("At least one embedding must be provided: X_umap_1d, X_umap_2d, or X_umap_3d.")
  }
  embeddings
}

.validate_embeddings <- function(embeddings) {
  n_cells <- NULL
  for (dim_chr in names(embeddings)) {
    dim_int <- as.integer(dim_chr)
    arr <- embeddings[[dim_chr]]
    if (length(dim(arr)) != 2L) {
      stop("X_umap_", dim_int, "d must be a 2D array.")
    }
    if (ncol(arr) != dim_int) {
      stop("X_umap_", dim_int, "d must have exactly ", dim_int, " columns, got ", ncol(arr), ".")
    }
    if (is.null(n_cells)) {
      n_cells <- nrow(arr)
    } else if (nrow(arr) != n_cells) {
      stop("All embeddings must have the same number of cells (rows).")
    }
  }
  if (
    is.null(n_cells) ||
      !is.numeric(n_cells) ||
      length(n_cells) != 1L ||
      n_cells < 1L ||
      n_cells > .Machine$integer.max
  ) {
    stop(
      "Embeddings must contain between 1 and 2,147,483,647 cells.",
      call. = FALSE
    )
  }
  as.integer(n_cells)
}

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
