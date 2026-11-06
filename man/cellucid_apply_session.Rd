% Hand-written and authoritative: this package runs no documentation generator
% for this page. Edit this file directly. tests/testthat/test-current-contract.R
% rebuilds every call of the \usage block below into a function and compares
% it to the real formals, so the two cannot drift apart.
\name{cellucid_apply_session}
\alias{cellucid_apply_session}
\alias{cellucid_apply_plan}
\alias{cellucid_history_applied}
\title{Apply one step of a saved Cellucid history to an object: a plan first}
\usage{
cellucid_apply_session(
  data,
  session,
  tracker = NULL,
  branch = NULL,
  step = NULL,
  ...,
  dry_run = TRUE,
  overwrite = FALSE,
  keep_session = FALSE,
  timeout = 60
)

cellucid_apply_plan(plan)

cellucid_history_applied(data)
}
\arguments{
\item{data}{The Seurat, SingleCellExperiment or SummarizedExperiment object
the session was saved on, its cells in their original order: a session names
every cell by its original position. An object read from a Cellucid-ready
store in its stored priority order is refused;
\code{\link{cellucid_restore_original_order}()} puts it back.}
\item{session}{A session \code{\link{cellucid_read_session}()} read, a
\code{.cellucid-session} file (what the viewer's Save session writes), or
what \code{\link{cellucid_request_session}()} returned: the session is then
collected first (in a Jupyter notebook in a later cell, since the viewer
answers between cells; elsewhere this waits for it).}
\item{tracker}{The tracker, by its name or id; \code{NULL} for the one the
viewer showed (the active tracker).}
\item{branch}{A branch of the tracker, by its name or id: its last step, or,
with \code{step}, the branch that step must be on.}
\item{step}{A step of the tracker, by its id (\code{"s4"};
\code{\link{cellucid_steps}()} lists them). With neither \code{branch} nor
\code{step}: the step the viewer showed (the tracker's head).}
\item{...}{The arguments \code{\link{cellucid_view}()} was given that choose
what is shown (\code{embeddings}, \code{obs}, \code{modalities}, \ldots), so
the object is checked against what the viewer showed.}
\item{dry_run}{\code{TRUE} (the default) returns the plan and changes
nothing; \code{FALSE} applies it at once.}
\item{overwrite}{Replace a column or gene set an earlier apply wrote with
other values. A column or gene set the person made is never replaced,
whatever this says.}
\item{keep_session}{Keep the whole session file, base64, in the apply's
record, so the history can be rebuilt from the object alone.}
\item{timeout}{Seconds to wait for a viewer's session, when \code{session}
is a request.}
\item{plan}{What \code{cellucid_apply_session(dry_run = TRUE)} returned.}
}
\value{
\code{cellucid_apply_session(dry_run = TRUE)}: a \code{cellucid_apply_plan},
a list of \code{session_sha256}, \code{tracker}, \code{branch} and
\code{step} (the ids chosen), \code{requires} (the obs columns and genes the
step's state uses: \code{obs}, \code{var}) and \code{rows}, a data frame of
one row per change, in domain order: \code{kind} (\code{"column-renamed"},
\code{"categories-renamed"}, \code{"categories-merged"},
\code{"field-deleted"}, \code{"user-field"}, \code{"gene-set"},
\code{"highlight"}, \code{"presentation"}), \code{domain}, \code{name} (the
column, category or gene set it writes), \code{column} (the column a category
is in), \code{from} (a list column: the old name or category), \code{cells}
(the cells it labels), \code{status} (\code{"write"}; \code{"unchanged"}: the
object holds it already; \code{"replace"}: an earlier apply's, with
\code{overwrite = TRUE}; \code{"not-applied"}) and \code{note} (why, in a
sentence). Printing it shows the rows.

\code{cellucid_apply_plan()}, and \code{cellucid_apply_session(dry_run =
FALSE)}: the object with the plan written -- its labels renamed in place,
each user-defined field and highlight group a column of its cell metadata
(Seurat \code{meta.data}; \code{colData}), its gene sets in
\code{cellucid$gene_sets} (origin \code{"session"}) and as columns of the
feature metadata the data mapping reads (see Details) -- and a record of the
apply appended to \code{cellucid$history_applied} (Seurat
\code{misc$cellucid}; \code{metadata(x)$cellucid}). Its attribute
\code{"cellucid_apply_report"} is \code{list(plan, applied, provenance)}:
\code{applied} is \code{FALSE} when nothing was written (applying the same
step again), and \code{provenance} the record. R copies on change: the object
passed in is not changed. A plan whose object changed since it was made is
refused: make it again.

\code{cellucid_history_applied()}: every apply recorded on the object, oldest
first, each a list: \code{session_sha256}, \code{tracker}, \code{branch},
\code{step}, \code{applied_at}, \code{cellucid_version}, \code{package}
(\code{"r"}), \code{columns} (the columns it wrote), \code{renamed} (each
rename it applied: \code{kind}, \code{column}, \code{from}, \code{to}),
\code{skipped} (each change not applied: \code{kind}, \code{name},
\code{note}) and \code{session} (the file, base64, with
\code{keep_session = TRUE}). cellucid-python's
\code{uns['cellucid']['history_applied']} holds the same records. A record
holds no \code{source} (where the viewer opened the dataset, which may be a
private link); with \code{keep_session = TRUE} it keeps the whole file,
its \code{source} included unless the file was shared without one.
}
\description{
What was made in the viewer, as it stood at one step of its history --
renamed fields and categories, user-defined fields (labels, annotations,
merged categories, pushed columns, copies of a column or a gene), gene sets
and highlight groups -- written onto the object it was saved on, after a
plan says what would change. The rules are cellucid-python's
\code{apply_cellucid_session_to_anndata()} and the web app's
\dQuote{Save AnnData with this state\ldots} copy
(\code{app/history/exports/anndata-copy.js}), row for row: the three are held
to one shared expected plan, the feature-metadata columns each writes
included.
}
\details{
Only what the chosen step needs is checked, before anything is written: the
dataset (its cell count, its id when the object's cellucid record names one,
the data mapping the session was saved with and its gene count, and the
cells' order when the object has the embedding the session names: a
cell-order digest of its stored values), the obs columns and genes the
step's state uses, and every value of the state against the keys this
version reads (a value it does not understand is refused, never applied in
part).

Applying the same step again changes nothing and says so (every row
\code{"unchanged"} or \code{"not-applied"}). A column the person made is
never replaced: its row says so and the rest applies. Presentation (filters,
colours, views, the camera, analyses) is reported, one row per domain that
differs from the tracker's opening, never applied. A session shared as its
steps only (without the values they leave behind) cannot be applied.

A highlight group's column is
\code{cellucid_highlight__<page>__<group>}, from the page's name and the
group's label (each run of characters other than letters and digits one
underscore); a user-defined field's is its name. A step's state knows the
cells the viewer had open then: a cell outside them is \code{NA} in a
highlight column and in a field's, never "not selected".

A gene set written or replaced is also a \code{TRUE}/\code{FALSE} column of
the feature metadata of the genes the session's data mapping reads: the
modality it reads (a Seurat assay's \code{meta.features}, or an
\code{Assay5}'s \code{meta.data}; \code{rowData} of the main experiment, or
of \code{altExp(x, name)} when the mapping reads that alternative
experiment). The column is \code{cellucid_gene_set__<name>} (the highlight
columns' rule; \code{gene_set} when nothing of the name is left;
\code{_2}, \code{_3} on clashes among the step's gene sets, named in their
order), \code{TRUE} on each feature whose gene name, as the mapping names
genes, the set lists; genes the object does not have are left out, with a
warning. A new column is added after the existing ones, and an earlier
apply's (\code{"replace"}) is replaced in place; a column of that name
beside a gene set no apply wrote is the object's own and is never replaced
(a warning says so). A mapping that reads no expression matrix has no
features to mark: the gene sets stay in \code{cellucid$gene_sets} only (a
warning says so). A gene set the step no longer holds has no row: its record
and its column stay. These columns are not rows of the plan, and
\code{columns} in the apply's record lists cell-metadata columns only.

\code{obs['cellucid_original_position']} (a Cellucid-ready store's record of
each cell's original position in priority order) is never written, renamed
or replaced: a plan that would is refused before anything is written, a dry
run included.
}
\seealso{
\code{\link{cellucid_read_session}()} and \code{\link{cellucid_trackers}()}
for what a session holds; \code{\link{cellucid_request_session}()} to ask a
viewer for its session; \code{\link{cellucid_pull}()} for the state a session
shows as columns; cellucid-python's \code{apply_cellucid_session_to_anndata()}.
}
\examples{
path <- system.file("extdata", "history_golden.cellucid-session", package = "cellucid")
session <- cellucid_read_session(path)
cellucid_steps(cellucid_tracker(session))[, c("id", "branch", "label")]
\dontrun{
# The object the session was saved on:
plan <- cellucid_apply_session(sce, session)
plan
sce <- cellucid_apply_plan(plan)
# A bookmarked step, or a branch's last step:
sce <- cellucid_apply_session(sce, session, step = "s5", dry_run = FALSE)
sce <- cellucid_apply_session(sce, session, branch = "Without score", dry_run = FALSE)
cellucid_history_applied(sce)

# From the viewer:
request <- cellucid_request_session(v)
sce <- cellucid_apply_session(sce, request, dry_run = FALSE)
}
}
