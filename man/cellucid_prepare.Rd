% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds the
% \usage block below into a function and compares it to the real formals of
% cellucid_prepare(), so the two cannot drift apart.
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
  dataset_name,
  dataset_description = NULL,
  dataset_id,
  source_name = NULL,
  source_url = NULL,
  source_citation = NULL,
  X_umap_1d = NULL,
  X_umap_2d = NULL,
  X_umap_3d = NULL,
  vector_fields = NULL,
  vector_field_default = NULL,
  created_at = NULL
)
}
\arguments{
\item{latent_space}{Numeric matrix-like of shape \code{(n_cells, n_dims)}. Required for per-cell outlier quantile calculations for categorical fields.}
\item{obs}{A \code{data.frame} of cell metadata with \code{n_cells} rows.
Every string category label must be text the viewer can show as written; see
the note on displayed text below.}
\item{var}{Optional \code{data.frame} of feature metadata with \code{n_genes} rows. Required if \code{gene_expression} is provided.}
\item{gene_expression}{Optional numeric matrix-like of shape \code{(n_cells, n_genes)}. Every value must be finite and, for unquantized output, representable as a nonzero finite float32 value unless it is exactly zero. Rows are always interpreted as cells; Cellucid never transposes the matrix or infers its semantic orientation from \code{var}.}
\item{var_gene_id_column}{Column name in \code{var} containing the gene names to
record. Use \code{NULL} (default) to use \code{rownames(var)}; that requires
row names a caller actually set, because a \code{data.frame} that never had any
reports the automatic sequence \code{"1"} to \code{"n_genes"} and those are row
numbers, not gene identifiers. Every string,
including \code{"index"}, names that exact column. Whatever this selects is
recorded faithfully: Cellucid performs no symbol lookup and ships no mapping, so
the caller decides what a gene is called. Every identifier in the selected axis
must be a non-empty string and must be distinct, because \code{gene_identifiers}
addresses \code{var} rows by identifier and a repeated one names no single row.}
\item{gene_identifiers}{Optional character vector specifying which genes to
export. If \code{NULL}, all genes are exported. Identifiers must be distinct and
must all be present in \code{var}. Payload files are named by integer index, so
an identifier is never a path; the exported ones must additionally be text the
viewer can draw exactly as stored, exactly as \code{obs_keys} narrows which
observation keys are held to that rule. A \code{var} row this argument leaves
out reaches no manifest and is not checked against it.}
\item{connectivities}{Optional adjacency/connectivity matrix of shape \code{(n_cells, n_cells)}. Sparse matrices from the \pkg{Matrix} package are supported. Weights must be finite and non-negative, the matrix must be exactly symmetric in topology and weight, sparse coordinates must be unique and omit stored zero values, and its diagonal must be exactly zero.}
\item{out_dir}{Output directory path. Defaults to \code{file.path(getwd(), "exports")}.}
\item{obs_keys}{Optional character vector of column names in \code{obs} to export. If \code{NULL}, all columns are exported. Exported keys must be distinct and must be text the viewer can draw exactly as stored; see the note on displayed text below. Numeric observation fields must contain only finite values and, for unquantized output, every nonzero value must remain nonzero and finite as float32.}
\item{centroid_outlier_quantile}{Quantile of distances to keep as inliers when computing centroids for categorical fields. Set to \code{NULL} to disable centroid computation.}
\item{centroid_min_points}{Minimum number of points in a category required to compute centroids and outlier quantiles.}
\item{force}{If \code{TRUE}, atomically replace an existing output
generation. If \code{FALSE}, require that \code{out_dir} does not already
exist. Independent Python or R exporters targeting the same directory are
always rejected immediately while another live process owns the generation.
The hidden sibling coordination file is intentionally persistent and is not
part of the published dataset.}
\item{var_quantization}{Integer bits for gene expression quantization (\code{8} or
\code{16}), or \code{NULL} for full float32 export. Quantization is derived from the
exact float32 values materialized by the viewer; a native-double range
that collapses to one float32 value is rejected. Individual source values
may round to zero when the float32 range remains non-collapsed.}
\item{obs_continuous_quantization}{Integer bits for continuous obs field
quantization (\code{8} or \code{16}), or \code{NULL} for full float32 export.
Quantization is derived from the exact float32 values materialized by the
viewer; a native-double range that collapses to one float32 value is
rejected. Individual source values may round to zero when the float32
range remains non-collapsed.}
\item{obs_categorical_dtype}{Required exact categorical code storage:
\code{"uint8"} or \code{"uint16"}. It has no default, so omitting it is
reported as a missing argument rather than as a value of the wrong type.}
\item{compression}{Optional gzip compression level (\code{1}-\code{9}). Use
\code{NULL} to disable compression. Every compressed payload carries the same
fixed gzip header on every platform: no filename, a Unix-epoch timestamp, and
the \code{OS} code \code{0xff}, which records no operating system. It is the
header the \code{cellucid} Python package writes. The deflate stream itself is
whatever the zlib underneath R produces at the requested level.}
\item{dataset_name}{Required non-empty dataset name for
\code{dataset_identity.json}, shown to the reader verbatim; see the note on
displayed text below. It has no default, so omitting it is reported as a
missing argument rather than as a value of the wrong type.}
\item{dataset_description}{Optional dataset description shown to the reader
verbatim. \code{NULL} (the default) and \code{""} both mean the dataset carries
no description and both publish \code{""}; any other value must satisfy the same
rule as a category label, see the note on displayed text below.}
\item{dataset_id}{Required portable dataset identifier: 1--180 ASCII letters,
numbers, \code{.}, \code{_}, or \code{-}, beginning with a letter or number,
not ending with \code{.}, and not a Windows device name. It has no default, so
omitting it is reported as a missing argument rather than as a value of the
wrong type.}
\item{source_name}{Optional non-empty data source name. Required whenever \code{source_url} or \code{source_citation} is supplied.}
\item{source_url}{Optional non-empty data source URL. May be supplied independently when \code{source_name} is present.}
\item{source_citation}{Optional non-empty citation text. May be supplied independently when \code{source_name} is present.}
\item{X_umap_1d}{Optional 1D embedding matrix of shape \code{(n_cells, 1)}. A
plain length-\code{n_cells} numeric vector is also accepted as the single column
it declares.}
\item{X_umap_2d}{Optional 2D embedding matrix of shape \code{(n_cells, 2)}.}
\item{X_umap_3d}{Optional 3D embedding matrix of shape \code{(n_cells, 3)}.}
\item{vector_fields}{Optional named list of per-cell displacement vectors. Each
name must exactly match \code{<field>_umap_<1|2|3>d}; the key declares the
dimension and it is never inferred from the array. Each element must be a
numeric matrix with \code{n_cells} rows and exactly the declared number of
columns; a \code{_1d} key also accepts a plain length-\code{n_cells} numeric
vector as the single column it declares.}
\item{vector_field_default}{Exact vector-field identifier to select by default. Required when \code{vector_fields} contains more than one field.}
\item{created_at}{Optional exact UTC-seconds timestamp in
\code{YYYY-MM-DDTHH:MM:SSZ} format. Defaults to the current UTC time. Supply
an exact value for reproducible \code{dataset_identity.json} metadata.}
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
  \item \code{var_manifest.json} and \code{var/} binaries for gene expression (one file per exported gene)
  \item \code{connectivity_manifest.json} and aligned \code{connectivity/edges.src.bin}, \code{edges.dst.bin}, and little-endian float64 \code{edges.weights.f64.bin} payloads for unique undirected edges
  \item \code{vectors/} binaries for each exported vector field and dimension
  \item \code{dataset_identity.json} describing the export and available embeddings
}

Note on payload names: no exported filename carries dataset content. A payload
on the \code{obs}, \code{var}, or \code{vectors} axis is named by the integer
position of its field in the export -- \code{var/0.values.f32},
\code{obs/1.codes.u8}, \code{vectors/0_2d.bin} -- and \code{connectivity/} uses
fixed neutral names. Every dataset therefore has the same directory listing.
Each manifest entry declares its own index as its first element, so
\code{fields} entries read \code{[index, name]} or
\code{[index, name, minValue, maxValue]}, \code{_continuousFields} entries read
\code{[index, key]} or \code{[index, key, minValue, maxValue]}, and
\code{_categoricalFields} entries read \code{[index, key, categories, dtype,
sentinel, centroids]} with \code{outlierMinValue} and \code{outlierMaxValue}
appended when the generated quantiles are quantized. The schema path patterns
\code{pathPattern}, \code{codesPathPattern}, and \code{outlierPathPattern}
substitute \code{\{index\}}. Within one axis directory the indices are exactly
\code{0} through \code{N-1} with each used once; \code{obs/} is written by both
manifest arrays, so those two arrays share one index space. The
\code{cellucid} Python package writes the identical layout.

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
in the file are rejected before publication. Quantized continuous values and
their manifest bounds are derived only after round-tripping through the
viewer's float32 domain. A field whose exported values are all one float32
value -- including source variation that collapses to one under float32
rounding -- is published as a constant field: equal manifest bounds and every
code \code{0}, which the viewer decodes back to that exact value. An
individual nonzero source value may round to zero when the resulting float32
range remains non-collapsed. Categorical outlier quantiles are
generated internally; cells in categories smaller than
\code{centroid_min_points} carry \code{NaN}, encoded with the reserved integer
marker only when those generated quantiles are quantized.

Note on required arguments: \code{obs_categorical_dtype}, \code{dataset_name},
and \code{dataset_id} carry no default, because no value of any of them is a
reasonable guess about someone else's dataset. Leaving one out is reported as
what it is -- \code{cellucid_prepare() is missing 1 required argument:} --
followed by the argument's name and what a valid value of it is, so a missing
argument is never described as a value of the wrong type. Every argument that
is left out is named in that one message. The \code{cellucid} Python package
requires the same three and reports the same set in one \code{TypeError}.

Note on displayed text: a category label, an exported gene name, an exported
\code{obs} key, a vector-field id, \code{dataset_name},
\code{dataset_description}, \code{source_name}, \code{source_url}, and
\code{source_citation} are drawn verbatim in the viewer's legend, its field
selector, and every exported figure, so a character with no glyph changes the
stored value without changing the picture. Each must be free of control
characters (\code{U+0001}--\code{U+001F}, \code{U+007F}--\code{U+009F}), free of
the zero-width characters \code{U+200B}, \code{U+2060}, and \code{U+FEFF}, and
without leading or trailing whitespace of any kind, \code{U+00A0} NO-BREAK
SPACE included. An empty category label is rejected, as are two labels in one
field that a whitespace-collapsing renderer draws identically, such as
\code{"T cell"} and \code{"T  cell"}. Nothing is trimmed: trimming would
rewrite an annotation the caller never asked to change, and would merge
\code{"Liver "} into a separate \code{"Liver"} category, moving cells between
them. The rejection names every offending label in the field at once. The
\code{cellucid} Python package enforces the identical rule.

Note on field identity: a payload filename is an integer index, so a gene name,
an \code{obs} key, and a vector-field id are never paths and carry no filename
rule -- no ASCII restriction, no case-insensitive collision rule, and no
Windows device-name rule. What each one must be is what it is for: a non-empty
string, distinct within its axis so a lookup resolves one field, and text the
viewer can draw exactly as it is stored, because a gene name and a category
label appear in the same legend. The rule covers the identifiers that are
exported: \code{obs_keys} and \code{gene_identifiers} narrow the export and
narrow this rule with it, and a \code{var} row left out reaches no manifest.
Gene IDs carry one further rule that spans the whole \code{var} and is
therefore not narrowed -- every ID must be a non-empty string and must be
distinct, because \code{gene_identifiers} addresses \code{var} rows by
identifier. \code{dataset_id} is the one identifier that still names a
directory, so it must be a 1--180-byte portable ASCII filename component: it
begins with a letter or digit, contains only letters, digits, \code{.},
\code{_}, or \code{-}, does not end with \code{.}, and is not a Windows device
name.
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
