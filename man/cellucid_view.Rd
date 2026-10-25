% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_view}
\alias{cellucid_view}
\alias{cellucid_serve}
\alias{cellucid_stop}
\alias{cellucid_stop_all}
\alias{cellucid_viewers}
\title{View an R object or a file in the Cellucid viewer}
\usage{
cellucid_view(
  data,
  ...,
  pane = c("auto", "viewer", "browser", "none"),
  port = 8765L,
  host = "127.0.0.1",
  allowed_hosts = NULL,
  serve_web_ui = TRUE,
  web_cache_dir = NULL,
  web_source_url = "https://www.cellucid.com",
  cache_bytes = 268435456,
  quiet = FALSE
)

cellucid_serve(
  data,
  ...,
  port = 8765L,
  host = "127.0.0.1",
  allowed_hosts = NULL,
  serve_web_ui = TRUE,
  web_cache_dir = NULL,
  web_source_url = "https://www.cellucid.com",
  cache_bytes = 268435456,
  quiet = FALSE
)

cellucid_stop(viewer)

cellucid_stop_all()

cellucid_viewers()
}
\arguments{
\item{data}{A Seurat, SingleCellExperiment or SummarizedExperiment object
(read as \code{cellucid_prepare(data = )} reads it), or the path of a file or
folder the viewer opens: a Cellucid-ready or any AnnData \code{.h5ad},
\code{.zarr} folder or \code{.zarr.zip}, a MuData Zarr, \code{.h5mu},
\code{.loom}, a 10x \code{.h5}, an \code{.rds} or \code{.rda} (every format of
the web app's registry, \code{system.file("input_formats.json", package =
"cellucid")}). An h5Seurat file is not opened: it is refused before anything
is served, with the route that converts it (read it with
\code{SeuratDisk::LoadH5Seurat()} and pass the object, or save it with
\code{saveRDS()}).}
\item{...}{For an object, the arguments of \code{\link{cellucid_prepare}()}
that choose what is shown: \code{dataset_id}, \code{dataset_name},
\code{embeddings}, \code{obs}, \code{modalities}, \code{modality},
\code{latent_key}, \code{embedding_defaults}, \code{obs_keys} and the
others, except \code{out_path}, \code{format}, \code{zarr_format} and
\code{force}. For a path: \code{dataset_id}, \code{dataset_name} and the launch
defaults (\code{embedding_defaults}, \code{latent_key}, \code{obs_keys},
\code{var_gene_id_column}, \code{vector_field_default},
\code{modality}, \code{normalize_embeddings},
\code{centroid_min_points}), which the browser checks against the file. Without
\code{dataset_id} and \code{dataset_name} (and no record of them in the object,
\code{misc$cellucid} or \code{metadata(x)$cellucid}), the dataset is named after
the expression passed as \code{data}, or the file. Of a file, only a Zarr
store's own identity (\code{uns['cellucid']}) is read; an \code{.h5ad} or
\code{.h5mu} is named after the file (this package reads no HDF5), where
cellucid-python's \code{serve_anndata()} reads its \code{uns['cellucid']} and
refuses a file that names none. Sessions are kept under the id, so give
\code{dataset_id} to keep them across both.}
\item{pane}{Where the viewer opens: \code{"auto"} the RStudio or Positron
Viewer pane when there is one, else the browser; \code{"viewer"} the pane only;
\code{"browser"} the default browser; \code{"none"} nowhere (the address is
printed).}
\item{port}{The port to serve on. The default 8765 is the one local address
Google sign-in (community annotation) accepts; when it is taken, the next free
port is used and the message says so. A port you choose is used or refused.
\code{0} picks a free port.}
\item{host}{The address to listen on. The default listens on this machine
only; \code{"0.0.0.0"} on every interface (others on the network can then open
the viewer).}
\item{allowed_hosts}{Extra host names the server answers to, needed only behind
a proxy that forwards to it under another name (RStudio Server, Posit
Workbench, a reverse proxy): bare names such as \code{"rstudio.example.org"},
matched on any port.}
\item{serve_web_ui}{\code{FALSE} serves the data and the channel only, for a
viewer page served elsewhere.}
\item{web_cache_dir}{The folder holding the web app's verified build;
\code{tools::R_user_dir("cellucid", "cache")/web} when \code{NULL}.}
\item{web_source_url}{The site the web app is downloaded from, at first use and
whenever it publishes a new build. \code{NULL} serves the cached build as it
is, with no network access.}
\item{cache_bytes}{The memory the served store may keep of what it makes by
a transform, in bytes (least recently used first out): chunks of text, of a
dense matrix, and the blocks of an object's genes x cells sparse matrix read
gene by gene. A chunk that is a slice of a vector (most of them) is made again
each time, in one pass at memory speed, and never kept. An object's sparse
matrix is never transposed whole: when the viewer reads a gene, the block of
genes around it (about eight values a cell, at most half of
\code{cache_bytes}) is made from the matrix as it is -- one binary search a
cell plus the block's own values, a fraction of a second for a million cells -- and kept
within this budget; a block dropped from it is made again when it is read
again. A dense matrix is read through its transpose, block by block, and
never copied whole.}
\item{quiet}{\code{TRUE} prints nothing: neither the address nor the notes
reading the data writes (columns not stored, gene names made unique).}
\item{viewer}{A viewer \code{cellucid_view()} or \code{cellucid_serve()}
returned.}
}
\value{
\code{cellucid_view()} and \code{cellucid_serve()} return the viewer,
invisibly: pass it to \code{\link{cellucid_highlight}()},
\code{\link{cellucid_on}()} and the other functions that talk to it.
\code{cellucid_stop()} returns the viewer invisibly; \code{cellucid_stop_all()}
the number of viewers stopped; \code{cellucid_viewers()} a list of the
running viewers.
}
\description{
\code{cellucid_view()} opens a Seurat, SingleCellExperiment or
SummarizedExperiment object, or a file, in the Cellucid viewer: in the RStudio
or Positron Viewer pane, or a browser. \code{cellucid_serve()} serves it for a
browser elsewhere (a computing cluster reached through an SSH tunnel) and
prints the address. Neither needs Python.
}
\details{
An object is read once, when the viewer starts, into the store
\code{cellucid_prepare(format = "zarr")} would write -- the same plan, the
same arrays and values -- held in memory: its metadata is built at once and
every chunk is made only when the viewer asks for it, so a large object opens
in about the time its cells' annotations and first embedding take. The served
store is not compressed (it declares no compressor), since the browser reads
it over the loopback interface, where deflating each chunk on R's one thread
cost far more than the bytes it saved; files \code{cellucid_prepare()} writes
are compressed. The browser reads it as it reads any Cellucid-ready AnnData,
with range requests. An object
changed after the viewer started is not shown changed (R copies on change):
view it again. A path is served from disk as it is, with exact
\code{Range} answers; the browser reads the file itself.

The server runs in the background (\pkg{httpuv}): the R console stays usable
while the viewer is open, and the viewer's requests are answered whenever R is
idle. It answers only requests that name it (the \code{Host} check refuses
DNS-rebinding pages with 421) -- unless it listens on a non-loopback IP address
(\code{host = "0.0.0.0"}, say) with no \code{allowed_hosts}, an explicit request
for network exposure whose names cannot be known -- and opens the viewer channel -- a WebSocket that
carries commands to the viewer and its events back -- only for a page of its
own, a loopback or an \code{allowed_hosts} origin. The channel is pinged every
20 seconds, so a proxy between the page and R (RStudio Server, Posit
Workbench) does not close it for being idle, and a socket that does not say
who it is within 60 seconds is closed. The web app itself is the
published build of \url{https://www.cellucid.com}, downloaded once into
\code{web_cache_dir} and checked file by file against its SHA-256 inventory
before any byte is served; \code{\link{cellucid_clear_web_cache}()} removes it.

Several viewers run at once, each on its own port. \code{cellucid_stop()}
closes the viewer's channel (the page shows it was stopped) and frees its
port; \code{cellucid_stop_all()} stops every viewer of the session.

A Viewer pane that renders without WebGL stays blank: view it with
\code{pane = "browser"}.
}
\seealso{
\code{\link{cellucid_highlight}()} and \code{\link{cellucid_on}()} to talk
to the viewer, \code{\link{cellucid_prepare}()} to write the same store to
share.
}
\examples{
if (interactive() && requireNamespace("SeuratObject", quietly = TRUE)) {
  data("pbmc_small", package = "SeuratObject")
  viewer <- cellucid_view(SeuratObject::UpdateSeuratObject(pbmc_small))
  cellucid_highlight(viewer, colnames(pbmc_small)[1:10])
  cellucid_stop(viewer)
}

# Serve the data without the web app, on a free port, and stop it.
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  data("pbmc_small", package = "SeuratObject")
  viewer <- cellucid_serve(
    SeuratObject::UpdateSeuratObject(pbmc_small),
    port = 0L, serve_web_ui = FALSE, quiet = TRUE
  )
  viewer
  cellucid_stop(viewer)
}
}
