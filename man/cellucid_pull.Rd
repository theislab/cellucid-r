% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_pull}
\alias{cellucid_pull}
\alias{cellucid_results}
\alias{cellucid_subset}
\alias{cellucid_gene_sets}
\title{Bring what was made in Cellucid back onto the object, in one call}
\usage{
cellucid_pull(
  data,
  source,
  ...,
  expected_dataset_id = NULL,
  highlights = TRUE,
  fields = TRUE,
  scores = TRUE,
  visibility = FALSE,
  gene_sets = TRUE,
  include_pushed = FALSE,
  highlight_prefix = "cellucid_highlight__",
  field_prefix = "",
  on_conflict = "error",
  prune = TRUE,
  return_summary = FALSE
)

cellucid_results(viewer)

cellucid_subset(data, viewer, which = NULL)

cellucid_gene_sets(viewer, symbols = FALSE)
}
\arguments{
\item{data}{The object the viewer shows, cells in the same order: a Seurat
object (columns go to \code{meta.data}), a SingleCellExperiment or
SummarizedExperiment (\code{colData}), or a data frame with one row per cell
named by cell. From a viewer, an object holding the viewer's cells in another
order is refused (nothing is written), as is one in a Cellucid-ready file's
stored order (its \code{cellucid$priority} record describes its cells;
\code{\link{cellucid_restore_original_order}()} puts it back).}
\item{source}{An open viewer (\code{\link{cellucid_view}()}): what it holds,
as it last sent it -- every result as it is made, and all of them whenever
its dataset opens -- so nothing is waited for, and in a Jupyter notebook a
lasso drawn and \code{cellucid_pull()} run in the next cell finds the lasso;
or what the web app exports with no R session at all: obs-column
\code{.csv} files (a header row, \code{obs_name} then each column's name; a
kind row, \code{kind} then each column's kind, \code{highlight},
\code{categorical} or \code{continuous}; a row per cell; each column written
as named and as its kind says, never guessed from its values, and a cell the
file does not list \code{NA}), \code{.gmt} gene sets and
\code{.cellucid-session} files; one file, several, or a folder of them.}
\item{...}{For a session file: the arguments \code{\link{cellucid_view}()} was
given that choose what is shown (as \code{\link{cellucid_apply_session}()}
takes them).}
\item{expected_dataset_id}{For a session file: the dataset id it must have
been saved on; by default the object's own (\code{cellucid$id}).}
\item{highlights}{Write the highlight groups.}
\item{fields}{Write the user-defined fields (labels, annotations, clusters).}
\item{scores}{Write the scores an analysis made (a gene signature's).}
\item{visibility}{Write the cells the filters leave visible, as
\code{cellucid_visible} (off by default).}
\item{gene_sets}{Write the gene sets.}
\item{include_pushed}{Also write the fields this R session pushed
(\code{\link{cellucid_push}()}), which the object holds already.}
\item{highlight_prefix}{A highlight group's column is
\code{<highlight_prefix><group id>}.}
\item{field_prefix}{A field's or score's column is \code{<field_prefix><key>}.}
\item{on_conflict}{What a column the object has already does, when no pull
from the same source (the viewer, or the same file) wrote it: \code{"error"}
(nothing is written), \code{"replace"}, or \code{"suffix"} (written as
\code{<name>_2} beside it). A result an earlier pull wrote is written to the
same column again, suffix and all. Session files follow the same rule.}
\item{prune}{For a viewer: remove the columns, and the gene sets, an earlier
viewer pull wrote whose result the viewer no longer holds (a highlight
deleted since). A result left out of this pull (\code{highlights = FALSE}, a
pushed field) keeps its column, and files never remove anything.}
\item{return_summary}{\code{TRUE} returns \code{list(data, summary)}: the
columns written (\code{columns}) and removed (\code{removed}), the gene sets
brought (\code{gene_sets}) and removed (\code{removed_gene_sets}), and the
results left out (\code{skipped}).}
\item{viewer}{A viewer \code{\link{cellucid_view}()} or
\code{\link{cellucid_serve}()} returned.}
\item{which}{A highlight group, by id (\code{"highlight_3"}) or by its label;
\code{"visible"} for the visible cells; \code{NULL} for the last selection.}
\item{symbols}{\code{FALSE} names each gene set's genes by identifier -- the
var index, \code{rownames(object)}; \code{TRUE} by the names the viewer shows
for them: those of the var column its data mapping names the genes by (such
as \code{feature_name}), else the identifiers themselves.}
}
\value{
\code{cellucid_pull()}: the object with the columns added (R copies on
change: the object passed in is not changed); highlight groups are
\code{TRUE}/\code{FALSE}, categorical fields factors with the viewer's
categories, continuous fields and scores numbers (\code{NA} where a cell has
none). When the viewer opened a subsample, a cell it did not open is
\code{NA}: unknown, not "not selected". The gene sets go to
\code{misc$cellucid$gene_sets} (Seurat), \code{metadata(x)$cellucid$gene_sets},
or the \code{cellucid} attribute of a data frame -- each its \code{genes} by
identifier (the var index), and from a viewer whose data mapping names the
genes by a var column the \code{symbols} shown --
beside \code{pulled}, the record of each column written: where it came from
(\code{pulled_from}: \code{"viewer"} or the file's name) and the result it
holds (\code{kind}, \code{id}). A gene set name used twice is
\code{"name (2)"}, as the web app and cellucid-python name it.

\code{cellucid_results()}: every result the viewer holds, as of the last
event it sent -- one epoch's: a reloaded page, or the viewer shown again,
replaces what an older page sent -- (it stops with what to do before the
viewer sent them all once, and while the viewer shows a dataset this viewer
did not serve, opened in it instead), each a list of
\code{kind} (\code{"highlight"}, \code{"categorical"}, \code{"continuous"},
\code{"geneSet"}, \code{"visibility"}, \code{"subset"}), \code{id},
\code{revision}, \code{n_obs}, \code{meta}, \code{cells} (1-based, \code{NULL}
for every cell), \code{values} and, for an object, \code{cell_names}.

\code{cellucid_subset()}: those cells of \code{data}.
\code{cellucid_gene_sets()}: a named list of character vectors, one per gene
set (the R session's, the shown marker genes', the shown gene signature's),
by identifier or by symbol (\code{symbols}).
}
\description{
Everything a person made in the viewer -- every highlight group, every
user-defined categorical field (labels, annotations, clusters), every
continuous field and score, every gene set -- written onto the object as
metadata columns in one call, from a live viewer or from the web app's
exports. cellucid-python's \code{viewer.pull(adata)} and
\code{cellucid.pull(adata, source)} write the same columns.
}
\details{
\code{\link{cellucid_apply_session}()} writes a saved session's highlight
groups and fields with the same column rules; \code{cellucid_pull()} of a
\code{.cellucid-session} file plans its columns the same way and names them
with every other column of the pull.
}
\examples{
if (interactive() && requireNamespace("SeuratObject", quietly = TRUE)) {
  data("pbmc_small", package = "SeuratObject")
  object <- SeuratObject::UpdateSeuratObject(pbmc_small)
  viewer <- cellucid_view(object)
  # Lasso some cells, label a few clusters, score a gene signature, then:
  object <- cellucid_pull(object, viewer)
  head(object[[]])
  # From the web app's exports, without a viewer:
  # object <- cellucid_pull(object, "~/Downloads/cellucid-exports")
  cellucid_stop(viewer)
}
}
