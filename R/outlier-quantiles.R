# Latent-space outlier quantiles, which cellucid_prepare() stores with the
# latent space it is given (latent_key) in uns['cellucid']['outlierQuantiles'].
#
# For each categorical obs column, each cell's quantile says how far it lies
# from the centre of its category in the latent space (an embedding such as
# scVI or PCA, never the 2D/3D layout): the share of its category at that
# distance or closer. The web app's outlier filter hides the cells above a
# threshold, and reads these instead of the whole latent space. The rule is
# src/latent_outliers.c, bit for bit the one cellucid-python's prepare()
# stores (prepare_data/_outliers.py) and the web app computes for a file that
# stores none (assets/js/data/latent-outliers.js), all three checked against
# the web repository's independent numpy reference
# (tests/testthat/fixtures/latent_outlier_conformance.json).
#
# Which columns: every obs column the viewer shows as categorical -- a
# factor, a character column (stored as a categorical or as text), a logical
# column with or without NA -- partitioned by its values (NA belongs to no
# category). A column none of whose categories reaches `min_points` cells
# gets no entry; a numeric column never does.

# The fewest cells a category needs for its cells to get outlier quantiles,
# unless centroid_min_points says otherwise; the web app's default too.
.default_outlier_min_points <- 10L

# One obs column as 0-based category codes (NA for a missing value) and the
# category count; NULL when the viewer does not show it as categorical. Which
# category gets which code does not change a quantile, only which cells share
# a category does.
.outlier_partition <- function(values) {
  if (is.factor(values)) {
    return(list(codes = as.integer(values) - 1L, n_categories = nlevels(values)))
  }
  if (is.object(values) || !is.null(dim(values))) {
    return(NULL)
  }
  if (is.logical(values)) {
    return(list(codes = as.integer(values), n_categories = 2L))
  }
  if (is.character(values)) {
    categories <- unique(values[!is.na(values)])
    return(list(codes = match(values, categories) - 1L, n_categories = length(categories)))
  }
  NULL
}

# uns['cellucid']['outlierQuantiles']: the latent key, the fewest cells a
# category needs, and one float32 quantile per cell (NaN where none) for every
# categorical obs column with a category of at least `min_points` cells, in
# obs column order.
.stored_outlier_quantiles <- function(obs, latent_key, latent, min_points) {
  fields <- list()
  for (column in names(obs)) {
    partition <- .outlier_partition(obs[[column]])
    if (is.null(partition)) next
    codes <- partition$codes
    counts <- tabulate(codes[!is.na(codes)] + 1L, nbins = partition$n_categories)
    if (!any(counts >= min_points)) next
    fields[[column]] <- .uns_floats32(.native_latent_outlier_quantiles(
      latent, codes, partition$n_categories, min_points
    ))
  }
  list(
    latentKey = .uns_string(latent_key),
    minPoints = .uns_integer(min_points),
    fields = fields
  )
}
