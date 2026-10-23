% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds the
% \usage block below into a function and compares it to the real formals of
% cellucid_prepare(), so the two cannot drift apart.
\name{cellucid_prepare}
\alias{cellucid_prepare}
\title{Write single-cell data as a Cellucid-ready AnnData store}
\usage{
cellucid_prepare(
  out_path,
  format,
  embeddings,
  dataset_id,
  dataset_name,
  obs = NULL,
  var = NULL,
  gene_expression = NULL,
  connectivities = NULL,
  vector_fields = NULL,
  gene_identifiers = NULL,
  dataset_description = NULL,
  citation_text = NULL,
  citation_url = NULL,
  source_name = NULL,
  source_url = NULL,
  embedding_defaults = NULL,
  latent_key = NULL,
  obs_keys = NULL,
  var_gene_id_column = NULL,
  vector_field_default = NULL,
  normalize_embeddings = NULL,
  centroid_min_points = NULL,
  zarr_format = NULL,
  force = FALSE
)
}
\arguments{
\item{out_path}{Required path of the store to write, ending in the suffix of
\code{format}: a folder whose name ends in \code{".zarr"} for
\code{format = "zarr"}, a file whose name ends in \code{".zarr.zip"} for
\code{format = "zarr.zip"}. A path that does not match \code{format} is refused
with the path to write instead. Missing parent folders are created.}
\item{format}{Required kind of output, with no default, because each suits a
different place and only you know where the store goes: \code{"zarr"} writes a
\code{.zarr} folder of small files, the way to publish on GitHub (no file over
100 MB, no Git LFS); \code{"zarr.zip"} writes the same store as one
\code{.zarr.zip} file, for everywhere else: single-file hosts (Zenodo,
Figshare, Hugging Face, S3, a lab server), local files, \code{cellucid serve}
and Jupyter. See the section \emph{Which format?}. Matched as
\code{match.arg()} matches. An \code{.h5ad} file is written by Python's
\code{cellucid.prepare(format = "h5ad")}; this package writes no HDF5, so it
depends on no HDF5 library.}
\item{embeddings}{Required named list mapping each \code{obsm} key to a numeric
matrix with one row per cell, such as \code{list(X_umap = umap, X_pca = pca)}.
A numeric \code{data.frame}, a \pkg{Matrix} matrix, or a plain numeric vector
(one column) is accepted too. Each one is stored whole, every column, as a
float32 array at its own coordinates; the viewer chooses the columns it plots
and scales them to its view. Every value must be finite and within the float32
range. Each key becomes a path in the store, see the note on names below.}
\item{dataset_id}{Required portable identifier of the dataset, stored as
\code{uns['cellucid']['id']}: 1--180 ASCII letters, numbers, \code{.},
\code{_}, or \code{-}, beginning with a letter or number, not ending with
\code{.}, and not a Windows device name. Sessions and community annotations are
keyed on it.}
\item{dataset_name}{Required non-empty name of the dataset, shown to the reader
verbatim and stored as \code{uns['cellucid']['name']} and \code{uns['title']};
see the note on displayed text below.}
\item{obs}{Optional \code{data.frame} of cell metadata, one row per cell. Every
column is stored, see the note on column encodings below.}
\item{var}{Optional \code{data.frame} of gene metadata, one row per column of
\code{gene_expression}. Every column is stored. Supplying it without
\code{gene_expression} is an error.}
\item{gene_expression}{Optional cells x genes matrix: a \pkg{Matrix} sparse
matrix (stored as a gene-major \code{csc_matrix}) or a dense numeric matrix
(stored as a dense array in tall, gene-major chunks). Values must be finite and
within the float32 range, and are stored as float32 -- or the whole matrix as
float64 when every value is a whole number and one exceeds \code{2^24} in
magnitude: integer counts float32 cannot hold exactly. That is decided on
the genes stored, after \code{gene_identifiers} selects them.
A stored zero of a sparse matrix is no value and is dropped; the matrix stays
sparse. Rows are always cells;
Cellucid never transposes the matrix, so transpose a genes x cells matrix with
\code{t()} first.}
\item{connectivities}{Optional cells x cells neighbour graph, a \pkg{Matrix}
sparse matrix or a numeric matrix, stored as
\code{obsp['connectivities']} (a float32 \code{csr_matrix}). It must be exactly
symmetric, with a zero diagonal and finite, non-negative weights; a stored zero
is not an edge and is dropped.}
\item{vector_fields}{Optional per-cell vectors drawn over an embedding, keyed by
field and then by the embedding each array belongs to:
\code{list(velocity = list(X_umap = v2, X_umap_3d = v3, X_branches = vb))}. Each
array has one row per cell and one column per column of its embedding, and is
stored under scVelo's name for it, \code{obsm['<field>_<basis>']} for the
embedding \code{X_<basis>} (\code{velocity_umap}, \code{velocity_umap_3d},
\code{velocity_branches}), which the viewer lists as the field
\code{'<field>'} and draws over that embedding. The embedding's key must
therefore start with \code{X_}, and each stored key must read back as its own
field and embedding: a field named like an embedding (\code{X_...}), or a key
another embedding's longer basis claims (field \code{T_fwd} on \code{X_umap}
beside an embedding \code{X_fwd_umap}), is refused.}
\item{gene_identifiers}{Optional character vector of the genes to store, in the
order given. \code{NULL} stores every gene. The identifiers are matched against
the gene names the viewer shows (\code{var_gene_id_column}, or the \code{var}
row names) after repeated names were made unique, and must be distinct.}
\item{dataset_description}{Optional description of the dataset, shown
verbatim; stored as \code{uns['cellucid']['description']}.}
\item{citation_text}{Optional text that says how to cite the dataset, stored as
\code{uns['cellucid']['citation']['text']}. Required whenever
\code{citation_url} is supplied.}
\item{citation_url}{Optional absolute \code{http://} or \code{https://} address
of the publication, stored as \code{uns['cellucid']['citation']['url']}.}
\item{source_name}{Optional name of where the data comes from, stored as
\code{uns['cellucid']['source']['name']}. Required whenever \code{source_url} is
supplied.}
\item{source_url}{Optional absolute \code{http://} or \code{https://} address of
the source, stored as \code{uns['cellucid']['source']['url']}.}
\item{embedding_defaults}{Optional launch default: the embedding each view
dimension shows first, by \code{"1d"}, \code{"2d"}, or \code{"3d"}. A value is an
\code{obsm} key, plotted by its first columns, such as
\code{c("2d" = "X_umap")}, or \code{list(key = , columns = )} naming the 1-based
columns to plot, such as \code{list("3d" = list(key = "X_pca", columns = 1:3))}.
Stored as \code{uns['cellucid']['defaults']['embedding<d>d']}, columns
0-based. Without it the viewer chooses.}
\item{latent_key}{Optional \code{embeddings} key of the latent space the
viewer's outlier filter measures in, such as \code{"X_pca"} or
\code{"X_scVI"}: a latent space such as scVI or PCA, never a 2D or 3D layout.
The outlier quantiles of every categorical \code{obs} column (each cell's
distance to its category's centre, as the share of its category at that
distance or closer) are computed here and stored in
\code{uns['cellucid']['outlierQuantiles']}, so the viewer filters outliers at
once instead of reading the whole latent space to compute them. See the
section \emph{Outlier quantiles}. Also stored as the launch default
\code{defaults['latentKey']}.}
\item{obs_keys}{Optional launch default: the \code{obs} columns the viewer
shows, in this order. Every column is stored either way. Stored as
\code{defaults['obsKeys']}.}
\item{var_gene_id_column}{Optional launch default: the \code{var} column whose
values name the genes in the viewer, such as gene symbols beside Ensembl row
names. \code{NULL} names genes by the \code{var} row names. A column of text or
a factor names genes as it is, an integer column by its digits; a double or
logical column cannot name genes and is refused. Stored as
\code{defaults['geneIdColumn']}.}
\item{vector_field_default}{Optional launch default: the field of
\code{vector_fields} selected first, by the name used there (\code{"velocity"}),
which is the id the viewer lists. One of its embeddings must be shown when the
file opens: one named in \code{embedding_defaults} or found by the viewer's
default rule, or any other 1 to 3 column embedding that is not
\code{latent_key}. Stored as \code{defaults['vectorFieldDefault']}.}
\item{normalize_embeddings}{Optional launch default: \code{FALSE} shows
embeddings at their own coordinates instead of scaled to the view. Stored as
\code{defaults['normalizeEmbeddings']}.}
\item{centroid_min_points}{Optional positive integer: the fewest cells a
category needs to get outlier quantiles (10 when \code{NULL}). Stored as the
launch default \code{defaults['centroidMinPoints']}, and with
\code{latent_key} as \code{uns['cellucid']['outlierQuantiles']['minPoints']}.}
\item{zarr_format}{Optional Zarr format of the store: \code{2} (the default,
when \code{NULL}), read by every \pkg{anndata} release since 0.8 and by
zarr-python 2 and 3, or \code{3}, read by \pkg{anndata} 0.11 and later
(zarr-python 3). The web app reads both equally fast. Python's
\code{cellucid.prepare(zarr_format = )} takes the same values.}
\item{force}{\code{TRUE} replaces an existing store at \code{out_path}, once the
new one is complete. An existing folder that is not a Zarr store is never
replaced.}
}
\value{
Invisibly, the absolute path of the written store.
}
\description{
\code{cellucid_prepare()} writes embeddings, cell and gene metadata, gene
expression, a neighbour graph, and vector fields as a \emph{Cellucid-ready
AnnData}: a standard AnnData store, written as AnnData Zarr (format 2, or 3
with \code{zarr_format = 3}), that Python's \pkg{anndata} reads without a
warning and the Cellucid web app opens directly, with the dataset's name,
citation, and launch defaults in \code{uns['cellucid']}.
}
\details{
The store is an AnnData Zarr group (Zarr format 2 unless \code{zarr_format = 3}):
\itemize{
  \item \code{X}: gene expression, gene-major float32 (float64 for integer
  counts above \code{2^24}), a \code{csc_matrix} with int32 \code{indices}
  and \code{indptr} and no stored zeros, or a dense array;
  \item \code{obs} and \code{var}: data frames, the index in \code{_index};
  \item \code{obsm}: every embedding and vector field, float32;
  \item \code{obsp/connectivities}: the neighbour graph, float32 \code{csr_matrix};
  \item \code{uns/title} and \code{uns/cellucid}: \code{schema} (1),
  \code{id}, \code{name}, and, when supplied, \code{description},
  \code{citation} (\code{text}, \code{url}), \code{source} (\code{name},
  \code{url}), \code{defaults} holding only the launch defaults that were
  set, and, with \code{latent_key}, \code{outlierQuantiles};
  \item the consolidated metadata (\code{.zmetadata}; format 3: the root
  \code{zarr.json}), so a store on a web host opens with one request.
}
Every element carries the \code{encoding-type} and \code{encoding-version}
attributes \pkg{anndata} 0.8 and later write. Every chunk of values wider
than one byte is byte-shuffled (numcodecs' \code{Shuffle} filter; format 3:
the \code{numcodecs.shuffle} codec), so the similar high bytes of neighbouring
values compress together -- except the expression values themselves, which
compress best as they are -- and every chunk is compressed with DEFLATE level 6
(the Zarr compressor \code{\{"id": "zlib", "level": 6\}}; format 3: the
\code{gzip} codec), which browsers decompress natively. A chunk is the unit of every read: arrays read whole (data
frame columns, embeddings as a chunk of rows with every column, the neighbour
graph) hold up to 262,144 (\code{2^18}) values per chunk, and the data and
indices of the expression matrix, which the viewer reads one gene at a time,
65,536 (\code{2^16}). Python's \code{cellucid.prepare()} writes the same store,
element for element, from the same input.

A \code{.zarr.zip} holds the same files, byte for byte, as entries of one ZIP
archive at its top level, each stored without further compression, so the web
app reads any chunk with one range request. The chunks come first and the
metadata documents last, the consolidated metadata the very last entry: the
web app reads the end of an archive when it opens it, and with it everything
that describes the store.

The store is written beside \code{out_path} under a temporary name and moved
into place only once it is complete: a failed or interrupted call leaves
nothing at \code{out_path}, and with \code{force = TRUE} the store it replaces
stays untouched until then.

Note on column encodings: the columns of \code{obs} and \code{var} are stored
with \pkg{anndata}'s own rules. A factor is a categorical with its levels in
their order (\code{ordered} for an ordered factor). A character column is a
categorical, its levels in code-point order, when a value repeats or is
missing, and a string array otherwise. Category codes take the integer type
pandas gives them: int8 for fewer than 127 categories, int16 for fewer than
32,767, else int32. Every other column keeps its own type, as in Python: a
logical column is a boolean array, or a nullable boolean when it holds
\code{NA}; an integer column is int32 (R's integer), or a nullable integer
when it holds \code{NA}; a double column is float64 whether or not its values
are whole numbers, with \code{NA} stored as NaN, and must hold no infinite
value. Any other column, a
\code{Date} for example, is refused with the conversion that makes it storable.

Note on names: cell names are the row names \code{obs},
\code{gene_expression}, and the embeddings agree on; where several are set they
must be the same names in the same order, and with none the cells are named
\code{"0"} to \code{"n - 1"}. Gene names are the \code{var} row names, else the
column names of \code{gene_expression}. Repeated names are made unique the way
\pkg{anndata}'s \code{obs_names_make_unique()} and
\code{var_names_make_unique()} do it -- the second \code{"MATR3"} becomes
\code{"MATR3-1"} -- and what changed is reported as a message. Every
\code{obs} and \code{var} column name and every \code{obsm} key becomes a path
in the store, so it must be a portable file name on every platform: no slash,
backslash, colon, asterisk, question mark, quotation mark, \code{<}, \code{>},
\code{|}, or control character, no trailing space or dot, not a
Windows device name, not \code{_index} for a column, and distinct from its
neighbours by more than letter case (\code{"Sex"} and \code{"sex"} are one
folder on macOS and Windows).

Note on displayed text: \code{dataset_name}, \code{dataset_description},
\code{citation_text}, and \code{source_name} are drawn verbatim, so each must be
free of control characters (\code{U+0001}--\code{U+001F},
\code{U+007F}--\code{U+009F}), of the zero-width characters \code{U+200B},
\code{U+2060}, and \code{U+FEFF}, and of leading or trailing whitespace of any
kind. Nothing is trimmed: trimming would rewrite text the caller never asked to
change.

Note on required arguments: \code{out_path}, \code{format},
\code{embeddings}, \code{dataset_id}, and \code{dataset_name} carry no default, because no value of
any of them is a reasonable guess about someone else's dataset. Leaving one out
is reported as what it is -- \code{cellucid_prepare() is missing 1 required
argument:} -- followed by the argument's name and what a valid value of it is.
}
\section{Which format?}{
Where the store goes decides \code{format}:
\itemize{
  \item a GitHub repository or its GitHub Pages site: \code{"zarr"}. GitHub
  refuses files over 100 MB, and a folder of small files publishes any size
  without Git LFS; GitHub answers over HTTP/2, 100 requests at a time, so a
  folder's one request per chunk costs nothing there.
  \item everything else -- a single-file host (Zenodo, Figshare, Hugging
  Face, S3, a lab web server), your own disk (the web app's file picker),
  \code{cellucid serve}, Jupyter: \code{"zarr.zip"}. One file is read in
  ranges that the web app merges into few requests; most single-file hosts
  answer over HTTP/1.1, six requests at a time, where a folder's one request
  per chunk makes many-gene and whole-matrix reads several times slower.
}
Measured on a 561,947-cell atlas, every read the web app makes, in total:
over GitHub's HTTP/2 (80 ms round trip, 50 MB/s) 8.0 s from the folder and
8.2 s from the \code{.zarr.zip}; over HTTP/1.1 (80 ms, 50 MB/s) 46.8 s and
9.0 s; from a local file in a browser 2.6 and 2.7 s.
In a browser, a local Zarr folder reads as fast as a single file. Python's
\code{cellucid.prepare(format = "h5ad")} writes an \code{.h5ad}, which goes
wherever a \code{.zarr.zip} goes; this package writes Zarr only, so it needs
no HDF5 library.
}
\section{Outlier quantiles}{
With \code{latent_key}, \code{uns['cellucid']['outlierQuantiles']} holds
\code{latentKey}, \code{minPoints} (\code{centroid_min_points}, or 10) and
\code{fields}: one float32 array of one value per cell for every \code{obs}
column the viewer shows as categorical -- a factor, a character column, a
logical column with or without \code{NA} -- in \code{obs} column order.
Its categories are its values; \code{NA} belongs to none. For each category of
at least \code{minPoints} cells, over the latent values as stored (float32):
the centroid is the column sums, accumulated cell by cell in cell order,
divided by the cell count; a cell's distance is the square root of its squared
differences from the centroid, accumulated column by column; its quantile is
the number of the category's cells at that distance or closer divided by the
cell count, both as float32. Every other cell is \code{NaN}. A column with no
category of \code{minPoints} cells gets no entry, and a numeric column never
does. The viewer's outlier filter hides the cells above a threshold of these
quantiles. Python's \code{cellucid.prepare()} stores the same bytes, and the
web app computes the same bytes for a file that stores none.
}
\examples{
umap <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
pca <- matrix(c(0, 3, 1, 2, 5, 1, 4, 0, 2, 1, 1, 1), ncol = 3, byrow = TRUE)
obs <- data.frame(
  cluster = factor(c("A", "A", "B", "B")),
  score = c(0.1, 0.2, 0.3, 0.4),
  row.names = paste0("cell", 1:4)
)
expr <- matrix(
  c(0, 1, 2, 3, 4, 5, 0, 0),
  nrow = 4,
  dimnames = list(rownames(obs), c("ENSG01", "ENSG02"))
)
var <- data.frame(symbol = c("CD3E", "MS4A1"), row.names = colnames(expr))

store <- cellucid_prepare(
  out_path = file.path(tempdir(), "example.zarr"),
  format = "zarr",
  embeddings = list(X_umap = umap, X_pca = pca),
  dataset_id = "cellucid-example",
  dataset_name = "Cellucid example",
  obs = obs,
  var = var,
  gene_expression = expr,
  embedding_defaults = c("2d" = "X_umap"),
  latent_key = "X_pca",
  var_gene_id_column = "symbol",
  force = TRUE
)
list.files(store, all.files = TRUE, no.. = TRUE)

# The same store as one file
cellucid_prepare(
  out_path = file.path(tempdir(), "example.zarr.zip"),
  format = "zarr.zip",
  embeddings = list(X_umap = umap),
  dataset_id = "cellucid-example",
  dataset_name = "Cellucid example",
  force = TRUE
)
}
