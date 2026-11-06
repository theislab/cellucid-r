% Hand-written and authoritative: this package runs no documentation generator
% for this page. Edit this file directly. tests/testthat/test-current-contract.R
% rebuilds every call of the \usage block below into a function and compares
% it to the real formals, so the two cannot drift apart.
\name{cellucid_turn}
\alias{cellucid_turn}
\alias{cellucid_last_history_change}
\title{The viewer's history, live: notebook turns and what changed}
\usage{
cellucid_turn(viewer, label, code)

cellucid_last_history_change(viewer)
}
\arguments{
\item{viewer}{What \code{\link{cellucid_view}()} or
\code{\link{cellucid_serve}()} returned.}
\item{label}{What the turn does, in a few words (\code{"Normalise and
cluster"}): the State panel shows it as the turn's name. White space is
collapsed.}
\item{code}{The code to run inside the turn: every action it sends is one
step of the turn.}
}
\value{
\code{cellucid_turn()}: the turn's id (\code{"r-"} and 12 hexadecimal
digits), invisibly, for \code{\link{cellucid_history_turn_undo}()} or a note
on the turn (\code{\link{cellucid_history_note_add}(viewer, target =
list(turn = id), text = ...)}).

\code{cellucid_last_history_change()}: the latest change of the viewer's
history it told, \code{list(kind, tracker, head, step, branch)}, or
\code{NULL} before the first: \code{kind} is \code{"step"} (a step
recorded: \code{step} is its id), \code{"moved"} (the step shown moved:
undo, redo, a jump, a branch switched to), \code{"edited"} (a note, a
bookmark, a branch renamed, locked, archived or deleted), \code{"tracker"}
(a tracker started, used, renamed or removed) or \code{"loaded"} (a session
loaded); \code{tracker} the tracker in use, \code{head} the step shown and
\code{branch} the active branch.
}
\description{
A notebook turn makes every action sent while its code runs one turn of the
viewer's history: the State panel shows them as one group, and they are
undone together. Turns do not nest.

The viewer tells every change of its history as a \code{history_changed}
event: \code{cellucid_last_history_change()} is the latest, and
\code{\link{cellucid_on}(viewer, "history_changed", handler)} runs code on
each. The history is read as a session
(\code{\link{cellucid_request_session}()}, then
\code{\link{cellucid_trackers}()}) and changed through its actions, the
generated \code{cellucid_history_*()} functions
(\code{\link{cellucid_history_undo}()}, \code{\link{cellucid_history_jump}()},
\code{\link{cellucid_history_note_add}()}, \ldots), as the State panel and
agents change it. cellucid-python's \code{viewer.turn()} and
\code{viewer.history} are the same calls.
}
\seealso{
\code{\link{cellucid_load_session}()}, \code{\link{cellucid_on}()},
\code{\link{cellucid_request_session}()}.
}
\examples{
\dontrun{
v <- cellucid_view(sce)
turn <- cellucid_turn(v, "Mark the T cells", {
  cellucid_color_by(v, field = list(source = "obs", key = "cell_type"))
  cellucid_highlight_cells(v, which(sce$cell_type == "T cell"))
})
cellucid_history_note_add(v, target = list(turn = turn), text = "T cells by cell_type.")
cellucid_on(v, "history_changed", function(change) message(change$kind, " ", change$head))
cellucid_last_history_change(v)
cellucid_history_turn_undo(v, turn = turn)
}
}
