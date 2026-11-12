% Hand-written and authoritative: this package runs no documentation generator
% for this page. Edit this file directly. tests/testthat/test-current-contract.R
% rebuilds every call of the \usage block below into a function and compares
% it to the real formals, so the two cannot drift apart.
\name{cellucid_read_session}
\alias{cellucid_read_session}
\alias{cellucid_conversations}
\title{Read a saved Cellucid session: its history, checked whole}
\usage{
cellucid_read_session(path)

cellucid_conversations(session, tracker = NULL)
}
\arguments{
\item{path}{A \code{.cellucid-session} file: what the viewer's Save session
or Share history writes, and what \code{\link{cellucid_session_result}()}
collects.}
\item{session}{What \code{cellucid_read_session()} or
\code{\link{cellucid_session_result}()} returned.}
\item{tracker}{A tracker of the session (by its name or id, or what
\code{\link{cellucid_tracker}()} returned); \code{NULL} for the active one.}
}
\value{
\code{cellucid_read_session()}: a \code{cellucid_session}, a list of
\code{path}, \code{sha256} (of the file), \code{manifest} (the file's
manifest as \pkg{jsonlite} reads JSON), \code{created_at}, \code{app} (the
web app that wrote it: \code{version}, \code{catalog} -- the digest of its
action catalog), \code{actions} (each action type the steps use, with the
digest of its parameters' schema), \code{scope} (what the share or save kept:
\code{trackers}, \code{branches}, \code{content} -- \code{"history"},
\code{"state"} for the current state only (one tracker: the dataset as
opened and one squashed step, read and applied as any session), or
\code{"steps"} for the steps only, no values --, \code{results}
(\code{"include"} or \code{"recompute"}), \code{notes}, \code{chats},
\code{journal}, \code{times} and \code{source}),
\code{dataset_fingerprint}, \code{source} (where the viewer opened the
dataset; \code{NULL} when the share left it out), \code{requires} (the obs
columns and genes the history uses), \code{chrome} (page chrome kept outside
the history, which the web app puts back after a load:
\code{list(pressed = list(<group id> = <button id>))}, the pressed button of
each pressed-button group, such as the selection tool; empty from a file
written without a page; read, never applied), \code{cells} and \code{genes} (the
subsamples it was saved on, or \code{NULL}), \code{active_tracker} (an id)
and \code{trackers} (each a \code{cellucid_tracker}, named by id: see
\code{\link{cellucid_trackers}()}). Printing it shows the dataset, the
trackers and the step shown.

\code{cellucid_conversations()}: the tracker's agent conversations the file
kept, each \code{list(conversation, transcript)}.
}
\description{
A session file (format 4) holds the history of a dataset in the viewer:
every step taken -- by a person, an agent, Python or R -- in one or more
trackers (independent histories), each a tree of steps whose branches are
the alternatives tried, with its notes, bookmarks, turns, journal and every
state the steps lead to. \code{cellucid_read_session()} reads and checks it
whole; \code{\link{cellucid_trackers}()} and its companions read the history;
\code{\link{cellucid_apply_session}()} writes the state at a step onto the
object it was saved on.
}
\details{
The file is untrusted input and is checked whole before anything is
returned, as the viewer checks it: every chunk inflated to exactly its
declared size, every stored array against its digest, every value against
its content hash once its parts are joined (a part, a sub-value several values
share and the file stores once, is joined once and its name never trusted; a
part that names itself through others, and a value whose text would pass
2^29 - 24 characters, are refused before anything is expanded), every record
against the history model
(\code{inst/history_model.schema.json}, the web app's, pinned by digest in
the web app, cellucid-python and this package), every tracker's tree whole
(one root, parents and children that agree, branches that are lines, groups
that are runs of consecutive steps holding whole turns, the head on its
active branch), what the share left out absent, the values the head shows
in the file's first (eager) chunks with every part and array they name, and
no value, part or array that nothing needs, none that something needs missing, and no trailing bytes. A damaged file is refused with "This session file is
damaged: " and what is wrong. A file written by a newer Cellucid than this
package, holding records this package does not know, is refused as newer
(update the package), never as damaged. A session file of an earlier format
(3, which held one state) or of a newer one is refused in the viewer's words.
}
\seealso{
\code{\link{cellucid_trackers}()}, \code{\link{cellucid_apply_session}()},
\code{\link{cellucid_request_session}()}; cellucid-python's
\code{read_session()}.
}
\examples{
# The shared golden history: three trackers (a replay, one removed), five
# branches (one forking after a squash in mid-branch), a notebook's and an
# agent's turn, notes, a bookmark.
path <- system.file("extdata", "history_golden.cellucid-session", package = "cellucid")
session <- cellucid_read_session(path)
session
session$scope$content
cellucid_conversations(session)[[1L]]$conversation
}
