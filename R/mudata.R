# The Cellucid-ready MuData: several feature spaces over the same cells (a
# Seurat object's assays, a SingleCellExperiment's main and alternative
# experiments), written as MuData's own Zarr layout
# (https://github.com/scverse/mudata, docs/io/spec.md), which mudata.read_zarr()
# reads and the web app opens (data/anndata/views/mudata.js):
#
#   /            encoding-type "MuData", encoding-version "0.1.0", axis 0
#   /obs         the cells' annotations (every cell, in the dataset's order)
#   /var         every modality's features, in modality order (no columns)
#   /obsm        the embeddings and vector fields, over every cell, and
#                mudata's masks obsm/<m>: bool (n_cells, 1), which cells <m>
#                measured
#   /obsp        the neighbour graph, over every cell
#   /varm, /varp empty
#   /uns         title and uns['cellucid'] (identity, citation, launch
#                defaults, outlier quantiles): the dataset's, so at the root
#   /obsmap/<m>  uint32, one per cell: its 1-based row in modality <m>, 0 when
#                <m> did not measure it
#   /varmap/<m>  uint32, one per root var row: its 1-based row in <m>, else 0
#   /mod         attribute mod-order: the modality names, in order
#   /mod/<m>     a Cellucid-ready AnnData minus uns['cellucid']: X gene-major,
#                obs (the modality's cells, no columns), var, empty layers,
#                obsm, obsp, varm, varp and uns
#
# R writes MuData as Zarr only (no .h5mu: no HDF5 library); cellucid-python's
# prepare(data = mdata) writes the same store, which
# tests/testthat/test-cross-writer.R checks element by element.

.write_mudata <- function(store, shared, spaces, names) {
  n_cells <- length(shared$obs_index)
  .zarr_group(store, "", .encoding("MuData", "0.1.0", axis = 0L))
  .ad_dataframe(store, "obs", shared$obs, shared$obs_index, "obs")
  feature_names <- unlist(lapply(spaces, function(space) space$index), use.names = FALSE)
  if (is.null(feature_names)) feature_names <- character(0)
  empty <- data.frame(row = seq_along(feature_names))[, 0L, drop = FALSE]
  .ad_dataframe(store, "var", empty, feature_names, "var")
  .ad_dict(store, "obsm")
  for (name in names(shared$matrices)) {
    .ad_dense(store, .zarr_key("obsm", name), shared$matrices[[name]])
  }
  for (name in names(shared$vectors$arrays)) {
    .ad_dense(store, .zarr_key("obsm", name), shared$vectors$arrays[[name]])
  }
  # mudata's per-modality masks: which cells each modality measured, one
  # boolean column (the viewer shows them as <name>:measured).
  for (index in seq_along(spaces)) {
    measured <- if (is.null(spaces[[index]]$positions)) {
      rep(TRUE, n_cells)
    } else {
      seq_len(n_cells) %in% spaces[[index]]$positions
    }
    .ad_mask(store, .zarr_key("obsm", names[[index]]), measured)
  }
  .ad_dict(store, "obsp")
  if (!is.null(shared$graph)) {
    graph <- shared$graph
    .ad_sparse(
      store, "obsp/connectivities", "csr_matrix", dim(graph),
      pointers = graph@p, indices = graph@i, values = graph@x
    )
  }
  .ad_dict(store, "varm")
  .ad_dict(store, "varp")
  .ad_uns_value(store, "uns", shared$uns)

  .ad_dict(store, "obsmap")
  for (index in seq_along(spaces)) {
    space <- spaces[[index]]
    map <- integer(n_cells)
    rows <- if (is.null(space$positions)) seq_len(n_cells) else space$positions
    map[rows] <- seq_along(rows)
    .ad_numeric_array(store, .zarr_key("obsmap", names[[index]]), map, "<u4")
  }
  .ad_dict(store, "varmap")
  offset <- 0L
  for (index in seq_along(spaces)) {
    n_features <- length(spaces[[index]]$index)
    map <- integer(length(feature_names))
    map[offset + seq_len(n_features)] <- seq_len(n_features)
    offset <- offset + n_features
    .ad_numeric_array(store, .zarr_key("varmap", names[[index]]), map, "<u4")
  }

  .zarr_group(store, "mod", list(`mod-order` = I(names)))
  no_columns <- function(n) data.frame(row = seq_len(n))[, 0L, drop = FALSE]
  for (index in seq_along(spaces)) {
    space <- spaces[[index]]
    prefix <- .zarr_key("mod", names[[index]])
    .zarr_group(store, prefix, .encoding("anndata", "0.1.0"))
    .write_anndata_group(
      store, prefix, space,
      obs = no_columns(length(space$obs_index)), obs_index = space$obs_index,
      matrices = list(), vectors = list(arrays = list()),
      graph = NULL, uns = structure(list(), names = character())
    )
  }
  invisible(NULL)
}
