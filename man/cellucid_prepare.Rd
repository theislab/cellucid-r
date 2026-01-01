% Generated manually for Bioconductor submission.
\name{cellucid_prepare}
\alias{cellucid_prepare}
\alias{prepare}
\title{Export single-cell data to the Cellucid viewer format}
\usage{
cellucid_prepare(
  latent_space = NULL,
  obs = NULL,
  var = NULL,
  gene_expression = NULL,
  var_gene_id_column = "index",
  gene_identifiers = NULL,
  connectivities = NULL,
  out_dir = file.path(getwd(), "exports"),
  obs_keys = NULL,
  centroid_outlier_quantile = 0.95,
  centroid_min_points = 10,
  obs_manifest_filename = "obs_manifest.json",
  obs_binary_dirname = "obs",
  var_manifest_filename = "var_manifest.json",
  var_binary_dirname = "var",
  connectivity_manifest_filename = "connectivity_manifest.json",
  connectivity_binary_dirname = "connectivity",
  force = FALSE,
  var_quantization = NULL,
  obs_continuous_quantization = NULL,
  obs_categorical_dtype = "auto",
  compression = NULL,
  dataset_name = NULL,
  dataset_description = NULL,
  dataset_id = NULL,
  source_name = NULL,
  source_url = NULL,
  source_citation = NULL,
  X_umap_1d = NULL,
  X_umap_2d = NULL,
  X_umap_3d = NULL,
  X_umap_4d = NULL,
  vector_fields = NULL
)

prepare(
  latent_space = NULL,
  obs = NULL,
  var = NULL,
  gene_expression = NULL,
  var_gene_id_column = "index",
  gene_identifiers = NULL,
  connectivities = NULL,
  out_dir = file.path(getwd(), "exports"),
  obs_keys = NULL,
  centroid_outlier_quantile = 0.95,
  centroid_min_points = 10,
  obs_manifest_filename = "obs_manifest.json",
  obs_binary_dirname = "obs",
  var_manifest_filename = "var_manifest.json",
  var_binary_dirname = "var",
  connectivity_manifest_filename = "connectivity_manifest.json",
  connectivity_binary_dirname = "connectivity",
  force = FALSE,
  var_quantization = NULL,
  obs_continuous_quantization = NULL,
  obs_categorical_dtype = "auto",
  compression = NULL,
  dataset_name = NULL,
  dataset_description = NULL,
  dataset_id = NULL,
  source_name = NULL,
  source_url = NULL,
  source_citation = NULL,
  X_umap_1d = NULL,
  X_umap_2d = NULL,
  X_umap_3d = NULL,
  X_umap_4d = NULL,
  vector_fields = NULL
)
}
\arguments{
\item{latent_space}{Numeric matrix-like of shape \code{(n_cells, n_dims)}. Required for per-cell outlier quantile calculations for categorical fields.}
\item{obs}{A \code{data.frame} of cell metadata with \code{n_cells} rows.}
\item{var}{Optional \code{data.frame} of feature metadata with \code{n_genes} rows. Required if \code{gene_expression} is provided.}
\item{gene_expression}{Optional numeric matrix-like of shape \code{(n_cells, n_genes)}.}
\item{var_gene_id_column}{Column name in \code{var} containing gene identifiers, or \code{"index"} (default) to use \code{rownames(var)}.}
\item{gene_identifiers}{Optional character vector specifying which genes to export. If \code{NULL}, all genes are exported.}
\item{connectivities}{Optional adjacency/connectivity matrix of shape \code{(n_cells, n_cells)}. Sparse matrices from the \pkg{Matrix} package are supported.}
\item{out_dir}{Output directory path. Defaults to \code{file.path(getwd(), "exports")}.}
\item{obs_keys}{Optional character vector of column names in \code{obs} to export. If \code{NULL}, all columns are exported.}
\item{centroid_outlier_quantile}{Quantile of distances to keep as inliers when computing centroids for categorical fields. Set to \code{NULL} to disable centroid computation.}
\item{centroid_min_points}{Minimum number of points in a category required to compute centroids and outlier quantiles.}
\item{obs_manifest_filename}{Filename for the obs manifest JSON.}
\item{obs_binary_dirname}{Directory name (under \code{out_dir}) for obs binaries.}
\item{var_manifest_filename}{Filename for the var manifest JSON.}
\item{var_binary_dirname}{Directory name (under \code{out_dir}) for var binaries.}
\item{connectivity_manifest_filename}{Filename for the connectivity manifest JSON.}
\item{connectivity_binary_dirname}{Directory name (under \code{out_dir}) for connectivity binaries.}
\item{force}{If \code{TRUE}, overwrite existing files; if \code{FALSE}, skip writing files whose manifest already exists (and skip embedding/vector files that already exist).}
\item{var_quantization}{Integer bits for gene expression quantization (\code{8} or \code{16}), or \code{NULL} for full float32 export.}
\item{obs_continuous_quantization}{Integer bits for continuous obs field quantization (\code{8} or \code{16}), or \code{NULL} for full float32 export.}
\item{obs_categorical_dtype}{One of \code{"auto"}, \code{"uint8"}, or \code{"uint16"} for categorical obs codes.}
\item{compression}{Optional gzip compression level (\code{1}-\code{9}). Use \code{NULL} or \code{0} to disable compression.}
\item{dataset_name}{Optional dataset name for \code{dataset_identity.json}.}
\item{dataset_description}{Optional dataset description.}
\item{dataset_id}{Optional dataset identifier. If \code{NULL}, a filesystem-safe identifier is derived from \code{dataset_name} or the \code{out_dir} name.}
\item{source_name}{Optional data source name for metadata.}
\item{source_url}{Optional data source URL for metadata.}
\item{source_citation}{Optional citation text for metadata.}
\item{X_umap_1d}{Optional 1D embedding matrix of shape \code{(n_cells, 1)}.}
\item{X_umap_2d}{Optional 2D embedding matrix of shape \code{(n_cells, 2)}.}
\item{X_umap_3d}{Optional 3D embedding matrix of shape \code{(n_cells, 3)}.}
\item{X_umap_4d}{Reserved for future development. Must be \code{NULL}.}
\item{vector_fields}{Optional named list of per-cell displacement vectors. Each element must be a numeric vector (1D) or matrix with 1/2/3 columns and \code{n_cells} rows.}
}
\value{
Invisibly returns \code{NULL}. Called for its side effects (file export).
}
\description{
\code{cellucid_prepare()} mirrors the \code{prepare()} function from the \pkg{cellucid} Python package and writes a directory of binary and JSON manifest files consumed by the Cellucid WebGL viewer.
}
\details{
The export directory includes:
\itemize{
  \item \code{points_1d.bin}, \code{points_2d.bin}, \code{points_3d.bin} (embedding coordinates; float32 row-major; optionally gzipped)
  \item \code{obs_manifest.json} and \code{obs/} binaries for continuous values, categorical codes, and per-cell outlier quantiles
  \item \code{var_manifest.json} and \code{var/} binaries for gene expression (one file per gene)
  \item \code{connectivity_manifest.json} and \code{connectivity/} binaries for unique undirected edges
  \item \code{dataset_identity.json} describing the export and available embeddings
}

Note on gene expression orientation: \code{gene_expression} must be shaped
\code{(n_cells, n_genes)} (cells as rows, genes as columns). Many R single-cell
containers store expression as \code{(n_genes, n_cells)}; transpose before
passing to \code{cellucid_prepare()}. See the vignette for Seurat and
SingleCellExperiment extraction examples.

Note on filename safety: \code{cellucid_prepare()} writes per-field and per-gene
files using filesystem-safe identifiers derived from obs keys and gene IDs. To
avoid silent overwrites, gene IDs and obs keys must be unique after
sanitization; duplicates or sanitization collisions are rejected.
}
\examples{
latent <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
obs <- data.frame(
  cluster = factor(c("A", "A", "B")),
  score = c(0.1, 0.2, 0.3)
)
umap2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
expr <- matrix(c(0, 1, 2, 3, 4, 5), nrow = 3, ncol = 2)
var <- data.frame(symbol = c("G1", "G2"))
rownames(var) <- var$symbol

out <- file.path(tempdir(), "cellucid_example")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

cellucid_prepare(
  latent_space = latent,
  obs = obs,
  var = var,
  gene_expression = expr,
  X_umap_2d = umap2,
  out_dir = out,
  centroid_min_points = 1,
  force = TRUE
)
}
