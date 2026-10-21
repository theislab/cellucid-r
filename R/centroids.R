# Where the viewer puts a category's label, and how far each cell sits from
# its own category's center. Both are measured from unrounded coordinates.

.compute_centroids_for_field <- function(coords, codes, categories, outlier_quantile = 0.95, min_points = 10) {
  if (nrow(coords) != length(codes)) {
    stop("coords and codes must have the same length.")
  }

  outlier_quantile <- .validate_centroid_outlier_quantile(
    outlier_quantile
  )
  if (is.null(outlier_quantile)) {
    stop(
      "outlier_quantile must be a finite number when computing centroids.",
      call. = FALSE
    )
  }

  centroids <- list()
  for (code in seq_along(categories) - 1L) {
    label <- categories[[code + 1L]]
    idx <- which(codes == code)
    n <- length(idx)
    if (n < min_points) {
      next
    }

    pts <- coords[idx, , drop = FALSE]
    center <- colMeans(pts)

    if (n > min_points) {
      diffs <- sweep(pts, 2, center, FUN = "-")
      dists <- sqrt(rowSums(diffs * diffs))
      thr <- as.numeric(stats::quantile(dists, probs = outlier_quantile, names = FALSE, type = 7))
      inlier_mask <- dists <= thr
      n_in <- sum(inlier_mask)
      if (n_in >= min_points) {
        pts_in <- pts[inlier_mask, , drop = FALSE]
        center <- colMeans(pts_in)
        used_count <- n_in
      } else {
        used_count <- n
      }
    } else {
      used_count <- n
    }

    centroids[[length(centroids) + 1L]] <- list(
      category = label,
      position = I(as.numeric(center)),
      n_points = as.integer(used_count)
    )
  }

  centroids
}

# One centroid list per exported embedding, keyed by embedding id in item
# order -- the centroidsByEmbedding object of a categorical obs entry.
.compute_centroids_for_all_embeddings <- function(items, codes, categories, outlier_quantile = 0.95, min_points = 10) {
  result <- list()
  for (item in items) {
    result[[item$id]] <- .compute_centroids_for_field(
      coords = item$coords,
      codes = codes,
      categories = categories,
      outlier_quantile = outlier_quantile,
      min_points = min_points
    )
  }
  result
}

# Per-cell outlier quantiles: how far a cell lies from its category's centre
# in the latent space, as the share of its category at that distance or
# closer. The rule is cellucid-python's and the web app's, bit for bit (see
# src/latent_outliers.c): float32 values widened to double, centroid sums in
# cell order, squared distances in column order, the share as float32. Cells
# with a missing code, or in a category under `min_points`, are NaN.
.compute_latent_space_quantiles <- function(latent, codes, categories, min_points = 10) {
  latent <- as.matrix(latent)
  storage.mode(latent) <- "double"
  codes <- as.integer(codes)
  .native_latent_outlier_quantiles(
    latent,
    codes,
    as.integer(length(categories)),
    as.integer(min_points)
  )
}
