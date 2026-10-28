% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_on}
\alias{cellucid_on}
\alias{cellucid_off}
\alias{cellucid_selection}
\alias{cellucid_wait_ready}
\alias{cellucid_request_session}
\alias{cellucid_session_result}
\title{React to what happens in an open Cellucid viewer}
\usage{
cellucid_on(viewer, event, handler)

cellucid_off(viewer, event = NULL)

cellucid_selection(viewer)

cellucid_wait_ready(viewer, timeout = 30)

cellucid_request_session(viewer)

cellucid_session_result(request, file = NULL, timeout = 60)
}
\arguments{
\item{viewer}{A viewer \code{\link{cellucid_view}()} or
\code{\link{cellucid_serve}()} returned.}
\item{event}{One of \code{"selection"}, \code{"hover"}, \code{"click"},
\code{"ready"}, \code{"data_mapping"}, \code{"command_error"}, \code{"pong"},
\code{"debug_snapshot"}, \code{"session_bundle"}, \code{"result"} (a result
made or changed, as \code{\link{cellucid_results}()} lists them),
\code{"result_removed"}, \code{"results"} (every result the viewer holds) and
\code{"color_by"} (the field it colours by: \code{field}, a gene by its
identifier, the var index; \code{source}, \code{"obs"} or \code{"var"};
and \code{label}, the name the viewer shows for it). A command the viewer
refused, or one that waited for the viewer and could not be built once it
could, arrives as \code{"command_error"}; with no handler for it, R warns.
For \code{cellucid_off()},
\code{NULL} removes every handler.}
\item{handler}{A function of one argument, the event: a list of its fields.
Cells are 1-based positions in every event, as the commands take them:
\code{selection} gives \code{cells} and (for an object) \code{cell_names};
\code{hover} and \code{click} give \code{cell} (\code{NULL} for no cell under
the pointer) and \code{cell_name}. cellucid-python's events are 0-based, as
Python indexes. The \code{pong} and \code{debug_snapshot} events answer
\code{\link{cellucid_debug_connection}()}; the channel's own keep-alive pings
are answered without an event.}
\item{timeout}{Seconds to wait for the viewer.}
\item{request}{What \code{cellucid_request_session()} returned.}
\item{file}{Where to save the session; a temporary file when \code{NULL}.}
}
\value{
\code{cellucid_on()} and \code{cellucid_off()} return the viewer invisibly.
\code{cellucid_selection()} returns the last selection, \code{list(cells,
cell_names, source)}: the selected cells' 1-based positions, their names (for an
object), and the tool that made it; \code{NULL} before the first.
\code{cellucid_wait_ready()} returns the viewer's \code{ready} event,
\code{list(n_cells, n_obs, dimensions)}: the cells open (fewer than
\code{n_obs} when the viewer opened a subsample), the cells in the dataset,
and the dimensions shown. \code{cellucid_request_session()} returns the
request (printed as waiting, answered or failed);
\code{cellucid_session_result()} the path of the \code{.cellucid-session}
file. A session is always asked for, then collected -- the one way, in a
notebook and outside one (cellucid-python's
\code{viewer.request_session_bundle()} and \code{request.result()}). In a
Jupyter notebook the viewer answers between cells, so collect it in a later
cell (\code{cellucid_session_result()} says so at once before then), or
apply it there with \code{\link{cellucid_apply_session}(object, request)};
where R can wait (a browser or RStudio viewer) \code{cellucid_session_result()}
waits for it.
}
\description{
The viewer sends its events back over its channel: a selection (every
selected cell, however many), the cell under the pointer, a click, the data
mapping it applied, a command it could not carry out. \code{cellucid_on()} runs
a handler for each event of a kind, whenever R is idle; a handler that fails
gives a warning and the others still run. \code{cellucid_request_session()}
asks the viewer for its whole session (what you highlighted, coloured, annotated and
arranged), the file the viewer's Save session writes.
}
\details{
Handlers run on R's main thread while R is idle, so an RStudio or Positron
console stays usable. \code{cellucid_wait_ready()} and
\code{cellucid_session_result()} wait by running R's event loop (outside a
notebook), and stop with what to
do when the viewer does not answer in time.
}
\examples{
if (interactive() && requireNamespace("SeuratObject", quietly = TRUE)) {
  data("pbmc_small", package = "SeuratObject")
  object <- SeuratObject::UpdateSeuratObject(pbmc_small)
  viewer <- cellucid_view(object)
  cellucid_on(viewer, "selection", function(event) {
    message(length(event$cells), " cells selected with the ", event$source)
  })
  # Select cells in the viewer, then:
  selected <- cellucid_selection(viewer)
  if (!is.null(selected)) subset(object, cells = selected$cell_names)
  cellucid_stop(viewer)
}
}
