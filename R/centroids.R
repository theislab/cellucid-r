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

.compute_centroids_for_all_dimensions <- function(embeddings, codes, categories, outlier_quantile = 0.95, min_points = 10) {
  result <- list()
  for (dim in names(embeddings)) {
    coords <- embeddings[[dim]]
    result[[dim]] <- .compute_centroids_for_field(
      coords = coords,
      codes = codes,
      categories = categories,
      outlier_quantile = outlier_quantile,
      min_points = min_points
    )
  }
  result
}

.compute_latent_space_quantiles <- function(latent, codes, categories, min_points = 10) {
  n_cells <- nrow(latent)
  quantiles <- rep(NaN, n_cells)

  for (code in seq_along(categories) - 1L) {
    idx <- which(codes == code)
    n <- length(idx)
    if (n < min_points) {
      next
    }
    pts <- latent[idx, , drop = FALSE]
    centroid <- colMeans(pts)
    diffs <- sweep(pts, 2, centroid, FUN = "-")
    dists <- sqrt(rowSums(diffs * diffs))
    ranks <- rank(dists, ties.method = "max")
    quantiles[idx] <- as.numeric(ranks) / n
  }

  quantiles
}
