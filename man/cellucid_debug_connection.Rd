% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_debug_connection}
\alias{cellucid_debug_connection}
\title{Check what works between R and an open Cellucid viewer}
\usage{
cellucid_debug_connection(viewer, timeout = 5)
}
\arguments{
\item{viewer}{A viewer \code{\link{cellucid_view}()} or
\code{\link{cellucid_serve}()} returned.}
\item{timeout}{Seconds to wait for each answer.}
}
\value{
A list: \code{viewer_id}, \code{viewer_url}, \code{server_url},
\code{server_running}, \code{channel_connected}, \code{web_ui} (the web
build's folder, and whether it is served), \code{state} (the viewer's
\code{ready} event and its last \code{command_error}); the server's own
\code{server_health}, \code{server_info} and \code{server_datasets} documents
(or \code{server_health_error} and so on); \code{anndata_store_probe}, the
store's \code{HEAD} answer (status, \code{Content-Length},
\code{Accept-Ranges}, \code{ETag}) and its first byte asked with
\code{Range}; and \code{frontend_roundtrip} and
\code{frontend_debug_snapshot}, each \code{list(ok, answer)} or
\code{list(ok = FALSE, error)}: a \code{ping} answered by a \code{pong}, and
the viewer's own view of its address, origin and connection.
}
\description{
The report cellucid-python's \code{AnnDataViewer.debug_connection()} gives,
for a viewer of this R session: whether the server answers its routes,
whether the store is read with byte ranges, and whether the viewer's
channel (a WebSocket) carries a command and its answer. It is the first
thing to look at when a viewer stays blank or does not react, for example
behind a proxy that does not pass WebSockets.
}
\details{
The server is asked over HTTP from this R session, through R's event loop,
which is where it answers; nothing leaves the machine. Without an open
channel the two \code{frontend_} entries say so instead of waiting.
}
\examples{
if (interactive() && requireNamespace("SeuratObject", quietly = TRUE)) {
  data("pbmc_small", package = "SeuratObject")
  viewer <- cellucid_view(SeuratObject::UpdateSeuratObject(pbmc_small))
  str(cellucid_debug_connection(viewer))
  cellucid_stop(viewer)
}
}
