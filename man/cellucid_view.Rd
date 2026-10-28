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
  pane = c("auto", "viewer", "browser", "notebook", "none"),
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
\item{pane}{Where the viewer opens: \code{"auto"} the cell's output in a
Jupyter notebook (an IRkernel session), else the RStudio or Positron Viewer
pane when there is one, else the browser; \code{"notebook"} the cell's output
only; \code{"viewer"} the pane only; \code{"browser"} the default browser;
\code{"none"} nowhere (the address is printed). See "In a Jupyter notebook"
below.}
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
is, with no network access. In a Jupyter notebook the browser loads the app
from this site itself (it must then be an http(s) address).}
\item{cache_bytes}{The memory the served store may keep of what it makes by
a transform, in bytes (least recently used first out): chunks of text, of a
dense matrix, and the blocks of an object's sparse matrix read gene by gene.
A chunk that is a slice of a vector (most of them) is made again each time,
in one pass at memory speed, and never kept. An object's sparse matrix is
never transposed or reordered whole: when the viewer reads a gene, the block
of genes around it (about eight values a cell, at most half of
\code{cache_bytes}) is made from the matrix as it is -- for a genes x cells
matrix one binary search a cell plus the block's own values, for a cells x
genes one its genes' values sorted by priority, a fraction of a second for a
million cells -- and kept within this budget; a block dropped from it is made
again when it is read again. A dense matrix is read through its transpose,
block by block, and never copied whole. The neighbour graph is not copied
either: each chunk of it (\eqn{2^{18}} entries) is made from the graph as it
is, both axes in priority order, when the viewer reads it, and not kept (at
one million cells and 20 million edges the served store holds no second
graph).}
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
SummarizedExperiment object, or a file, in the Cellucid viewer: in a Jupyter
notebook cell's output, the RStudio or Positron Viewer pane, or a browser. \code{cellucid_serve()} serves it for a
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
with range requests. Its cells are in priority order (seed 0, as
\code{cellucid_prepare()} writes them by default), so a subsample the viewer
opens is the first rows of every array and the first entries of every gene,
and only its cells are sent; the original positions, which selections,
highlights and labels name, are the ones \code{cellucid_prepare()} gives:
an object read from a file in priority order whose record
(\code{misc$cellucid$priority}, \code{metadata(x)$cellucid$priority}) still
describes its cells is put back in its original order first; any other is
taken in its own order, and a \code{cellucid_original_position} column
without such a record names nothing (a warning says so). An object
changed after the viewer started is not shown changed (R copies on change):
view it again. A path is served from disk as it is, with exact
\code{Range} answers; the browser reads the file itself. R has no HDF5
reader and no reader of Zarr chunks to serve a file in priority order, so
for a file whose cells are not in priority order already (a Zarr store
\code{cellucid_prepare()} wrote is) the viewer says so as it starts: a
subsample then reads the storage blocks that hold its cells. To send only
its cells, read the file into R and view the object, write it once with
\code{cellucid_prepare()}, or serve it from Python (\code{cellucid serve}),
which serves \code{.h5ad}, \code{.h5mu}, Loom, Cell Ranger and AnnData and
MuData Zarr files in priority order. An \code{.rds}, \code{.rda} or
MatrixMarket folder is read whole once by the browser. A file served as it
is is listed to the viewer with its size and a strong \code{ETag} -- its
size, modification time, status-change time and identity (inode; on
Windows NTFS's change time and file index), read from the file opened when
the viewer asks -- the validator its bytes are answered with: the
browser keeps what it read of the file (its byte cache, gene indexes, an
\code{.rds} converted once) under that version across visits, and a file
rewritten while it is served, in place or replaced, is a new version, never
read into what was kept of the old one.

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
\section{In a Jupyter notebook}{
In a Jupyter notebook with the R kernel (IRkernel: JupyterLab, Jupyter
Notebook 7, VS Code), \code{cellucid_view()} shows the viewer in the cell's
output: a notebook widget -- the one cellucid-python's \code{show_anndata()}
shows -- that frames the Cellucid app loaded from \code{web_source_url} and
reaches this R session through the kernel's own connection (a Jupyter comm).
No server is started and no port is opened, so it works the same on a laptop
and on a remote JupyterHub. It needs the \pkg{IRkernel} and \pkg{IRdisplay}
packages and the \code{anywidget} Python package where the notebook server
runs (\code{pip install anywidget}), which shows the widget.

The kernel hands the widget the viewer's messages between cells: while a cell
runs, the viewer waits for data it has not read yet, and a call that waits
for the viewer (\code{cellucid_wait_ready()}, or
\code{cellucid_session_result()} before the viewer answered) stops at once
with what to do instead; commands never wait.
Hooks registered with \code{cellucid_on()} run, and \code{cellucid_selection()}
is current, whenever the kernel is idle. \code{port}, \code{host},
\code{allowed_hosts}, \code{serve_web_ui} and \code{web_cache_dir} do not
apply.

In Google Colab's R runtime (detected from its \code{COLAB_RELEASE_TAG}
environment variable) the widget's display also names Colab's widget manager,
as \code{google.colab.output.enable_custom_widget_manager()} does for Python:
Colab shows third-party widgets only then. The widget's module is then loaded
from jsDelivr.

A notebook front end that does not show widgets shows the widget's text
instead. A server is no way around it inside IRkernel, which never runs R's
event loop between cells: \code{cellucid_serve()} in a notebook answers
nothing until a cell runs. Serve from an R session outside the notebook (R in
a terminal, on the same machine) and open the viewer in a browser tab, or
write the data with \code{cellucid_prepare()} and open it at
\url{https://www.cellucid.com} or with the Python package's
\code{show_anndata()}.
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
