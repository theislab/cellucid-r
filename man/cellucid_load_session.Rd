% Hand-written and authoritative: this package runs no documentation generator
% for this page. Edit this file directly. tests/testthat/test-current-contract.R
% rebuilds every call of the \usage block below into a function and compares
% it to the real formals, so the two cannot drift apart.
\name{cellucid_load_session}
\alias{cellucid_load_session}
\alias{cellucid_load_session_result}
\title{Push a saved Cellucid session into the viewer}
\usage{
cellucid_load_session(viewer, path)

cellucid_load_session_result(request, timeout = 60)
}
\arguments{
\item{viewer}{What \code{\link{cellucid_view}()} or
\code{\link{cellucid_serve}()} returned.}
\item{path}{A \code{.cellucid-session} file saved on the dataset the viewer
shows.}
\item{request}{What \code{cellucid_load_session()} returned.}
\item{timeout}{Seconds to wait for the viewer's answer.}
}
\value{
\code{cellucid_load_session()}: a \code{cellucid_load_session_request},
invisibly, at once.

\code{cellucid_load_session_result()}: what the viewer shows once it loaded
the file, \code{list(tracker = list(id, name), head)}. It stops with the
viewer's reason when the viewer did not load it (another dataset's session,
a damaged file). In a Jupyter notebook the viewer answers between cells:
collect the result in a later cell.
}
\description{
The session's trackers are added to the viewer's history beside the person's
own, and the viewer shows the session's active tracker, as opening the file
in the viewer does. The file goes to the viewer as one binary
\code{load_session} command, answered by \code{session_loaded}.
cellucid-python's \code{viewer.load_session()} is the same call.
}
\seealso{
\code{\link{cellucid_turn}()}, \code{\link{cellucid_request_session}()},
\code{\link{cellucid_read_session}()}.
}
\examples{
\dontrun{
v <- cellucid_view(sce)
request <- cellucid_load_session(v, "pbmc.cellucid-session")
cellucid_load_session_result(request)
}
}
