% Generated manually for CRAN submission.
\name{cellucid_prepare}
\alias{cellucid_prepare}
\title{Export single-cell data to the Cellucid viewer format}
\usage{
cellucid_prepare(
  latent_space = NULL,
  obs = NULL,
  var = NULL,
  gene_expression = NULL,
  var_gene_id_column = NULL,
  gene_identifiers = NULL,
  connectivities = NULL,
  out_dir = file.path(getwd(), "exports"),
  obs_keys = NULL,
  centroid_outlier_quantile = 0.95,
  centroid_min_points = 10,
  force = FALSE,
  var_quantization = NULL,
  obs_continuous_quantization = NULL,
  obs_categorical_dtype,
  compression = NULL,
  dataset_name = NULL,
  dataset_description = "",
  dataset_id = NULL,
  source_name = NULL,
  source_url = NULL,
  source_citation = NULL,
  X_umap_1d = NULL,
  X_umap_2d = NULL,
  X_umap_3d = NULL,
  vector_fields = NULL,
  vector_field_default = NULL
)
}
\arguments{
\item{latent_space}{Numeric matrix-like of shape \code{(n_cells, n_dims)}. Required for per-cell outlier quantile calculations for categorical fields.}
\item{obs}{A \code{data.frame} of cell metadata with \code{n_cells} rows.}
\item{var}{Optional \code{data.frame} of feature metadata with \code{n_genes} rows. Required if \code{gene_expression} is provided.}
\item{gene_expression}{Optional numeric matrix-like of shape \code{(n_cells, n_genes)}. Every value must be finite and, for unquantized output, representable as a nonzero finite float32 value unless it is exactly zero. Rows are always interpreted as cells; Cellucid never transposes the matrix or infers its semantic orientation from \code{var}.}
\item{var_gene_id_column}{Column name in \code{var} containing gene identifiers. Use \code{NULL} (default) to use \code{rownames(var)}. Every string, including \code{"index"}, names that exact column.}
\item{gene_identifiers}{Optional character vector specifying which genes to export. If \code{NULL}, all genes are exported.}
\item{connectivities}{Optional adjacency/connectivity matrix of shape \code{(n_cells, n_cells)}. Sparse matrices from the \pkg{Matrix} package are supported. Weights must be finite and non-negative, the matrix must be exactly symmetric in topology and weight, sparse coordinates must be unique and omit stored zero values, and its diagonal must be exactly zero.}
\item{out_dir}{Output directory path. Defaults to \code{file.path(getwd(), "exports")}.}
\item{obs_keys}{Optional character vector of column names in \code{obs} to export. If \code{NULL}, all columns are exported. Numeric observation fields must contain only finite values and, for unquantized output, every nonzero value must remain nonzero and finite as float32.}
\item{centroid_outlier_quantile}{Quantile of distances to keep as inliers when computing centroids for categorical fields. Set to \code{NULL} to disable centroid computation.}
\item{centroid_min_points}{Minimum number of points in a category required to compute centroids and outlier quantiles.}
\item{force}{If \code{TRUE}, atomically replace an existing output generation. If \code{FALSE}, require that \code{out_dir} does not already exist.}
\item{var_quantization}{Integer bits for gene expression quantization (\code{8} or \code{16}), or \code{NULL} for full float32 export.}
\item{obs_continuous_quantization}{Integer bits for continuous obs field quantization (\code{8} or \code{16}), or \code{NULL} for full float32 export.}
\item{obs_categorical_dtype}{Required exact categorical code storage: \code{"uint8"} or \code{"uint16"}.}
\item{compression}{Optional gzip compression level (\code{1}-\code{9}). Use \code{NULL} to disable compression.}
\item{dataset_name}{Required non-empty dataset name without leading or trailing
whitespace for \code{dataset_identity.json}.}
\item{dataset_description}{Dataset description. Defaults to an empty string.}
\item{dataset_id}{Required portable dataset identifier.}
\item{source_name}{Optional data source name. If supplied, \code{source_url} and \code{source_citation} are also required.}
\item{source_url}{Optional data source URL. If supplied, \code{source_name} and \code{source_citation} are also required.}
\item{source_citation}{Optional citation text. If supplied, \code{source_name} and \code{source_url} are also required.}
\item{X_umap_1d}{Optional 1D embedding matrix of shape \code{(n_cells, 1)}.}
\item{X_umap_2d}{Optional 2D embedding matrix of shape \code{(n_cells, 2)}.}
\item{X_umap_3d}{Optional 3D embedding matrix of shape \code{(n_cells, 3)}.}
\item{vector_fields}{Optional named list of per-cell displacement vectors. Each element must be a numeric vector (1D) or matrix with 1/2/3 columns and \code{n_cells} rows.}
\item{vector_field_default}{Exact vector-field identifier to select by default. Required when \code{vector_fields} contains more than one field.}
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
  \item \code{connectivity_manifest.json} and aligned \code{connectivity/edges.src.bin}, \code{edges.dst.bin}, and little-endian float64 \code{edges.weights.f64.bin} payloads for unique undirected edges
  \item \code{dataset_identity.json} describing the export and available embeddings
}

Note on gene expression orientation: \code{gene_expression} must be shaped
\code{(n_cells, n_genes)} (cells as rows, genes as columns). Many R single-cell
containers store expression as \code{(n_genes, n_cells)}; transpose before
passing to \code{cellucid_prepare()}. See the vignette for Seurat and
SingleCellExperiment extraction examples. Shape checks cannot detect an
incorrectly oriented square matrix, so callers must verify orientation.

Continuous observation fields, gene expression, and vector fields must contain
only finite values. Every nonzero value written as float32, after required
embedding/vector scaling, must have an absolute magnitude from \code{2^-149}
through \code{(2 - 2^-23) * 2^127}; values that would become zero or infinity
in the file are rejected before publication. Categorical outlier quantiles are
generated internally; cells in categories smaller than
\code{centroid_min_points} carry \code{NaN}, encoded with the reserved integer
marker only when those generated quantiles are quantized.

Note on filename safety: gene IDs, obs keys, dataset IDs, and vector-field IDs
must already be 1--180-byte portable ASCII filename components. They begin with
a letter or digit, contain only letters, digits, \code{.}, \code{_}, or
\code{-}, do not end with \code{.}, are not Windows device names, and are
unique under case-insensitive filesystem comparison.
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
  dataset_id = "cellucid-example",
  dataset_name = "Cellucid example",
  latent_space = latent,
  obs = obs,
  var = var,
  gene_expression = expr,
  X_umap_2d = umap2,
  out_dir = out,
  centroid_min_points = 1,
  force = TRUE,
  obs_categorical_dtype = "uint16"
)
}
