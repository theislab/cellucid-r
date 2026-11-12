% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_mcp}
\alias{cellucid_mcp}
\title{Let an agent app on this computer work in a Cellucid viewer}
\usage{
cellucid_mcp(
  open = NULL,
  viewer = NULL,
  port = NULL,
  browser = TRUE,
  web_source_url = "https://www.cellucid.com",
  web_cache_dir = NULL
)
}
\arguments{
\item{open}{A file or folder to serve and open in a viewer of this process
(any path \code{\link{cellucid_view}()} takes); its page opens in the
browser, and the viewer stops when the server ends. \code{NULL} drives a
viewer opened elsewhere on this computer.}
\item{viewer}{Which published viewer to drive: its viewer id, dataset id or
dataset name. \code{NULL} takes the most recently opened one whose page is
open, else the most recently opened one. Not with \code{open}.}
\item{port}{The port of the viewer \code{open} starts; \code{NULL} takes 8765,
or the next free one.}
\item{browser}{Whether \code{open} opens the viewer's page in the browser
(\code{FALSE}: its address is written to standard error).}
\item{web_source_url}{Where the viewer \code{open} starts loads its web app
from, as for \code{\link{cellucid_view}()}.}
\item{web_cache_dir}{Where that viewer keeps its copy of the web app;
\code{NULL} is the user's cache folder.}
}
\value{
\code{NULL}, invisibly, once the app closes the server's standard input.
}
\description{
A local MCP (Model Context Protocol) server on this process's standard input
and output, for an agent app on this computer -- Claude Code, Claude
Desktop, Cursor -- to work in a Cellucid viewer: the app starts it, and it
runs until the app ends it. It is the R counterpart of cellucid-python's
\code{cellucid mcp}, with the same tools, sentences and status, and either
drives a viewer of either package. Nothing passes through a relay or the
deployed site: the server talks to the viewer on this computer.

The tools are the web app's own (the web connector's, generated into the
package as \file{mcp_tools.json}; while a viewer page is open, the page's
own list, the same tools), and \code{cellucid_status}, which says which
viewers are found, which one is driven and what to tell the person when none
is open. The list is the same whatever the person's settings and never
changes while the server runs: the page tells the agent what the person's
settings allow after its first tool result, and again after a result once
they changed. The page runs every call as the app's agent (\dQuote{Claude Code
(local)}) under the person's permission settings: a change runs, or waits
for the person's answer in the viewer's confirm card -- asked in the app
too, as an MCP elicitation, when the app can ask -- or is refused with one
sentence. Every step is in the person's history, and the person can Stop
the agent's turn or undo it whole.
}
\details{
Register it with the app, for example
\preformatted{claude mcp add cellucid -- Rscript -e "cellucid::cellucid_mcp()"}
or, in an app's \code{"mcpServers"} settings,
\preformatted{"cellucid": {"command": "Rscript", "args": ["-e", "cellucid::cellucid_mcp()"]}}
and to have it open a file itself,
\code{Rscript -e "cellucid::cellucid_mcp(open = 'data.h5ad')"}.

\strong{Which viewer.} Every viewer with a server --
\code{\link{cellucid_view}()} outside a notebook, \code{\link{cellucid_serve}()}
-- is published to agent apps on this computer while it runs, as
cellucid-python publishes its viewers: one record in the agents folder
(\env{CELLUCID_AGENTS_DIR} when set, an absolute path; else
\file{\%LOCALAPPDATA\%/cellucid/agents} on Windows,
\file{~/Library/Caches/cellucid/agents} on macOS,
\file{$XDG_CACHE_HOME/cellucid/agents} or \file{~/.cache/cellucid/agents}
elsewhere), readable by this user only (0600, the folder made 0700), naming
this computer, where the viewer is served (never its page address, which
carries the viewer's secret: an agent is not given it), and the viewer's
endpoint on its own server with the secret token it takes. The endpoint
answers requests from this computer only, never a browser's. The record
goes when the viewer stops or R exits. On Linux and macOS a folder another
user owns, or others may write in, is refused: no viewer is published there
or read from it, and the agent is told why. A home folder shared between
computers is fine: a record of another computer is skipped and left alone.
\env{CELLUCID_AGENTS}\code{=off} (or \code{0}, \code{false}, \code{no}; any
case) publishes no viewer. A notebook widget's viewer (\code{cellucid_view()}
in an IRkernel notebook) has no server and is not published, nor is a viewer
bound to one non-loopback interface (\code{host}).

\code{cellucid_mcp()} drives the viewer \code{viewer} names, else the newest
whose page is open: when the page of the one it drives is gone and another
viewer's is open, it moves to that one; when the one it drives stops, it
finds another. A call whose page goes (reloaded or closed) before it ends
fails with a sentence that says whether it was done is not known, and the
app is told to list the tools again when a page takes the viewer back.

\strong{What the agent is sent.} What the agent reads (names, counts,
summaries, values, a picture of the view) goes to the agent app and from
there to the model provider the app uses; nothing passes through a Cellucid
server. A dataset's \dQuote{Agents see names and summaries only} switch, set
in the viewer (Permissions in the chat settings or the Agents panel), applies
to this agent as to every other: the page leaves values of single cells,
extremes (a minimum, a maximum, a median) and pictures out of its results,
and refuses summaries of groups an agent made by naming cells. It is not
anonymization: a summary of a very small group can still reveal its cells'
values. Text written by others (labels, notes, dataset text) reaches the
agent marked as data, never as instructions.

\strong{The protocol.} JSON-RPC 2.0, one message per line: the handshake of
protocol versions 2024-11-05 to 2025-11-25 (\code{initialize}), \code{ping},
\code{tools/list}, \code{tools/call} and \code{notifications/cancelled} (the
call is stopped in the viewer), and the server's \code{elicitation/create}
(a form, to an app whose elicitation capability is empty or names forms;
withdrawn with \code{notifications/cancelled} when the person answers in the
viewer first) and \code{notifications/progress} (when a call names a
progress token, counting 1, 2, ...); the server's tools capability says
\code{listChanged: false}, as its list never changes; and the
per-request protocol 2026-07-28 (\code{server/discover}), whose clients are
asked in the viewer's card only. Nothing but protocol messages is written to
standard output: what R prints while the server runs goes to standard error.

\strong{One thread.} Standard input is read without blocking while R's event
loop runs (\code{later::run_now()}), so the viewer \code{open} serves keeps
answering its page, and a call waiting for the person never stops the server
from reading the app's next message. A request to another process's viewer
runs R's event loop while it waits (for at most 10 seconds). On Windows the
app's pipe is read without blocking too; a console is refused (the server
is started by the app, not typed into).
}
\seealso{
\code{\link{cellucid_view}()} opens a viewer this server can drive.
}
\examples{
\dontrun{
# Run by the agent app, which talks to it over standard input and output:
#   Rscript -e "cellucid::cellucid_mcp()"
cellucid_mcp()
}
}
