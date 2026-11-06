% Hand-written and authoritative: this package runs no documentation generator
% for this page. Edit this file directly. tests/testthat/test-current-contract.R
% rebuilds every call of the \usage block below into a function and compares
% it to the real formals, so the two cannot drift apart.
\name{cellucid_trackers}
\alias{cellucid_trackers}
\alias{cellucid_tracker}
\alias{cellucid_branches}
\alias{cellucid_branch}
\alias{cellucid_groups}
\alias{cellucid_group}
\alias{cellucid_steps}
\alias{cellucid_step}
\alias{cellucid_step_state}
\alias{cellucid_step_diff}
\alias{cellucid_stale}
\alias{cellucid_notes}
\alias{cellucid_turns}
\alias{cellucid_journal}
\alias{cellucid_methods}
\alias{cellucid_pseudocode}
\title{Read a Cellucid history: trackers, branches, groups, steps, notes, turns}
\usage{
cellucid_trackers(session)

cellucid_tracker(session, tracker = NULL)

cellucid_branches(tracker)

cellucid_branch(tracker, branch = NULL)

cellucid_groups(tracker)

cellucid_group(tracker, group)

cellucid_steps(tracker, branch = NULL)

cellucid_step(tracker, step)

cellucid_step_state(tracker, step)

cellucid_step_diff(tracker, a, b, other = NULL)

cellucid_stale(tracker)

cellucid_notes(tracker)

cellucid_turns(tracker)

cellucid_journal(tracker)

cellucid_methods(tracker, branch = NULL)

cellucid_pseudocode(tracker, branch = NULL)
}
\arguments{
\item{session}{What \code{\link{cellucid_read_session}()} or
\code{\link{cellucid_session_result}()} returned.}
\item{tracker}{For \code{cellucid_tracker()}: a tracker's name (as the web app
compares names: white space collapsed, case ignored; among the trackers not
removed) or id, \code{NULL} for the active one. Elsewhere: what
\code{cellucid_tracker()} returned.}
\item{branch}{A branch's name or id; \code{NULL} for the active branch
(\code{cellucid_branch()}), every step (\code{cellucid_steps()}) or the main
branch (\code{cellucid_methods()}, \code{cellucid_pseudocode()}).}
\item{group}{A group's name (compared as the web app compares names) or id
(\code{"g0"}).}
\item{step}{A step's id (\code{"s4"}).}
\item{a}{A step's id: what differs from it ...}
\item{b}{... to this step.}
\item{other}{The tracker \code{b} is in, when it is another one (a replay);
\code{NULL} for \code{tracker}.}
}
\value{
\code{cellucid_trackers()}: a data frame, one row per tracker: \code{id},
\code{name}, \code{active}, \code{removed}, \code{start} (\code{"opened"},
\code{"shown"}, \code{"replay"}, \ldots), \code{head}, \code{active_branch},
\code{main_branch}, \code{recording} (whether changes become steps) and
\code{paused_at} (since when they do not; \code{NA} while recording), both
derived from the tracker's \code{paused},
\code{steps}, \code{branches}, \code{groups} and \code{created_at}.

\code{cellucid_tracker()}: a \code{cellucid_tracker}, the history model's
Tracker with every field in snake case (\code{head_id}, \code{active_branch},
\code{main_branch}, \code{paused} (\code{NULL} while it records, else
\code{at}, since when, \code{by}, who paused it and everyone who acted while
it was, and \code{draft}, \code{NULL} in a file), \code{start},
\code{removed}, \code{branches}, \code{groups}, \code{turns}, \code{notes},
\code{journal}, \code{steps}: each a list of the model's records). A record not found is refused naming the ones there are.

\code{cellucid_branches()}: a data frame, one row per branch: \code{id},
\code{name}, \code{base} (the step it forks from), \code{tip}, \code{main},
\code{locked}, \code{archived}, \code{n_notes} and \code{created_at};
\code{cellucid_branch()}: one branch's record.

\code{cellucid_groups()}: a data frame, one row per group (a named run of
consecutive steps of one branch, a folder over them: nothing is merged):
\code{id}, \code{name}, \code{branch}, \code{first}, \code{last},
\code{steps} (how many), \code{n_notes} and \code{created_at};
\code{cellucid_group()}: one group's record, with \code{steps}, its steps'
ids. Each step names its group (\code{group}).

\code{cellucid_steps()}: a data frame, one row per step (with \code{branch},
its path from the root to its tip), the root first and every step before its
children: \code{id}, \code{branch} (the branch it was made on),
\code{parent}, \code{kind} (\code{"root"}, \code{"action"},
\code{"outside"} -- a change made outside an action -- or \code{"squashed"}),
\code{type} (the action), \code{v} (the version of the action's parameters
it was recorded with), \code{label} (what it did, in a sentence),
\code{code} (its line of pseudocode, as recorded), \code{origin} (\code{"user"},
\code{"agent"}, \code{"host"}, \code{"session"}, \code{"replay"},
\code{"startup"} or \code{"outside"}), \code{agent}, \code{model},
\code{host}, \code{turn}, \code{at} (a \code{POSIXct}), \code{status}
(\code{"applied"} or \code{"unapplicable"}), \code{stale} (its result is no
longer what its inputs give), \code{bookmark}, \code{group} (the group it is
in) and \code{n_notes} -- the columns of cellucid-python's
\code{Tracker.steps}.

A squashed step has no action of its own (\code{type}, \code{v} and
\code{code} are \code{NA}): it is the net change of several steps merged
into one, of the history cleared keeping the state, or of what changed while
recording was paused; its record's \code{squashed} says why
(\code{cause}: \code{"steps"}, \code{"clear"}, \code{"paused"} or
\code{"share"}), how many steps it replaced (\code{count}) and who acted in
it (\code{by}). Its state is exact, as any step's. A step's \code{code} is
its pseudocode as recorded: parameters of an earlier version (\code{v} below
the action's version in \code{inst/action_catalog.json}) are shown and
applied as they were; the web app upgrades them only to replay them.

\code{cellucid_step()}: one step's record: \code{by} (who took it),
\code{params}, \code{outcome}, \code{touches}, \code{before}, \code{after}
and \code{reads} (the domains' value hashes), \code{kind}, \code{v},
\code{squashed}, \code{notes}, \code{bookmark}, \code{group}, \ldots

\code{cellucid_step_state()}: the viewer's state after the step, a named list
by state domain (\code{"cells"}, \code{"mapping"}, \code{"field-overlays"},
\code{"gene-sets"}, \code{"highlights"}, \code{"camera"}, \ldots): each value
as the viewer captures it, lists as \pkg{jsonlite} reads JSON, every typed
array an R vector whose attribute \code{"cellucid_array"} names its
JavaScript type, the numbers JSON cannot hold \code{NaN}, \code{Inf},
\code{-Inf} and \code{-0}. Refused for a share of the steps only. No picture
of a step is ever stored: the web app draws one from the live view when it
is asked for.

\code{cellucid_step_diff()}: a data frame of \code{domain}, \code{name} and
\code{line} (what differs, in the web app's words), one row per change; its
attribute \code{"common"} is the latest step both paths share (within one
tracker; \code{NULL} across two).

\code{cellucid_stale()}: the results stale at the head, each
\code{list(step, cause, domains)}: a step whose output is still shown while
an input it read changed since, the first later step that changed it, and
those inputs.

\code{cellucid_notes()}: a data frame, one row per note (on the tracker, a
branch, a group, a step or a turn): \code{id}, \code{on}, \code{target}, \code{text}
(Markdown), \code{origin}, \code{agent}, \code{at}, \code{edited_at},
\code{from_reply}.

\code{cellucid_turns()}: a data frame, one row per turn (a notebook's or an
agent's group of steps): \code{id}, \code{label}, \code{origin},
\code{agent}, \code{model}, \code{host}, \code{branch},
\code{conversation}, \code{message}, \code{started_at}, \code{ended_at},
\code{steps} (how many), \code{outside} (how many acts outside the state
it did, such as a vote published, which Undo turn does not take back) and
\code{summary} (the turn's first note).

\code{cellucid_journal()}: a data frame, one row per edit of the history
itself (a note, a bookmark, a group, a squash, a branch locked, a step
deleted): \code{id},
\code{op}, \code{at}, \code{origin}, \code{label}, \code{purged}.

\code{cellucid_methods()}: the branch as methods text, in the web app's
words (its rules: \code{app/history/exports/methods.js}): the tracker's
first note (its abstract) opens it; a squashed step is said as what it stands
for ("changes made with recording paused (filters, colors)", "3 steps merged
into one (highlights)", a fresh start, a shared state) with a sentence that
only its resulting state was kept; when a listed step is in a group, each run
of steps follows under its heading (the group's name, or "Before"/"After" the
group). \code{cellucid_pseudocode()}: the branch as pseudocode, one stored
line per step, \code{# --- QC ---} and \code{# --- end of QC ---} around a
group's steps, \code{# Bookmark: <name>} after a bookmarked step (a squashed step's
bookmarks kept from the steps it replaced marked "(on a step merged into
this one)"), and a
flagged step, a change made outside an action and a squashed step (values
only) as comments.
}
\description{
The reading model of a Cellucid history, function for function
cellucid-python's \code{cellucid.history}: every record is the history
model's (\code{inst/history_model.schema.json}), its fields named in snake
case. Both packages are held to the same shared golden history and the same
expected reading.
}
\seealso{
\code{\link{cellucid_read_session}()}, \code{\link{cellucid_apply_session}()};
cellucid-python's \code{Tracker}, \code{Branch}, \code{Group} and \code{Step}.
}
\examples{
path <- system.file("extdata", "history_golden.cellucid-session", package = "cellucid")
session <- cellucid_read_session(path)
cellucid_trackers(session)[, c("name", "active", "removed", "head")]
tracker <- cellucid_tracker(session, "Annotation")
cellucid_branches(tracker)[, c("id", "name", "base", "tip", "locked", "archived")]
cellucid_groups(tracker)[, c("id", "name", "first", "last", "steps")]
cellucid_group(tracker, "QC")$steps
cellucid_steps(tracker, "Without score")[, c("id", "label", "origin")]
cellucid_turns(tracker)[, c("id", "label", "origin", "steps", "summary")]
cellucid_step_diff(tracker, "s5", "s13")
cellucid_stale(tracker)
cat(cellucid_methods(tracker))
cat(cellucid_pseudocode(tracker), sep = "\n")
state <- cellucid_step_state(tracker, "s13")
names(state)
}
