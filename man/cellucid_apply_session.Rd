% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_apply_session}
\alias{cellucid_apply_session}
\title{Write a saved Cellucid session's highlights and fields onto an object}
\usage{
cellucid_apply_session(
  data,
  session,
  ...,
  expected_dataset_id = NULL,
  add_highlights = TRUE,
  highlights_prefix = "cellucid_highlight__",
  add_user_defined_fields = TRUE,
  user_defined_prefix = "",
  include_deleted_user_defined_fields = FALSE,
  store_metadata = TRUE,
  return_summary = FALSE,
  timeout = 60
)
}
\arguments{
\item{data}{The Seurat, SingleCellExperiment or SummarizedExperiment object
the session was saved on, its cells in their original order: a session names
every cell by its original position (its column in the order the cells were
given, whatever order a store keeps them in). An object read from a
Cellucid-ready store in its stored priority order (its
\code{cellucid$priority} record describes its cells) is refused;
\code{\link{cellucid_restore_original_order}()} puts it back. One whose
record no longer describes its cells (subset, combined or reordered since it
was read) is taken row by row as it is, with a warning, as the web app reads
such a file.}
\item{session}{A \code{.cellucid-session} file (what the viewer saves, and
\code{\link{cellucid_session_result}()} returns), or what
\code{\link{cellucid_request_session}()} returned: the session is asked for
first and applied here (in a Jupyter notebook in a later cell, since the
viewer answers between cells; elsewhere this waits for it).}
\item{...}{The arguments \code{\link{cellucid_view}()} was given that choose
what is shown (\code{embeddings}, \code{obs}, \code{modalities}, \ldots), so
the object is checked against what the viewer showed.}
\item{expected_dataset_id}{The dataset id the session must have been saved
on: the viewer's \code{dataset_id}. Required for a file; for a viewer, the
viewer's own.}
\item{add_highlights}{Add one \code{TRUE}/\code{FALSE} column per highlight
group.}
\item{highlights_prefix}{Prefix of each highlight column, followed by the
group id (\code{highlight_1}, \ldots).}
\item{add_user_defined_fields}{Add the session's user-defined fields: a
categorical field as a factor, a copied continuous field as numbers.}
\item{user_defined_prefix}{Prefix of each user-defined field's column,
followed by its name.}
\item{include_deleted_user_defined_fields}{Also add user-defined fields the
session deleted (never purged ones).}
\item{store_metadata}{Record the session -- its manifest, dataset
fingerprint, the columns added, its highlight pages and field overlays -- in
\code{misc$cellucid$session} (Seurat) or \code{metadata(x)$cellucid$session}
(SingleCellExperiment, SummarizedExperiment).}
\item{return_summary}{\code{TRUE} returns \code{list(data, added_columns)}
instead of the object alone.}
\item{timeout}{Seconds to wait for a viewer's session, when \code{session}
is a request.}
}
\value{
The object with the added columns in its cell metadata (Seurat
\code{meta.data}; \code{colData}), or \code{list(data, added_columns)} with
\code{return_summary = TRUE}. R copies on change: the object passed in is
not changed.
}
\description{
What was made in the viewer -- highlight groups (selections, lasso, category
picks) and user-defined fields (categorical fields built from highlight
pages, copies of a column or a gene) -- added to the object's cell metadata.
The rules are cellucid-python's \code{apply_cellucid_session_to_anndata()},
check for check: the same session applies to the same dataset the same way
in both languages.
}
\details{
Every check runs before anything is added, so a session that does not match
leaves nothing half-written: the dataset id, the cell count, the data
mapping the session was saved with (every embedding, graph, gene-name
column and expression matrix it names must exist on the object as the
viewer showed it), the obs columns and genes the session uses, every chunk
of the file, and every cell index. A highlight group stored relative to a
category is refused when that category no longer holds the cells the
session recorded. A column the object already has is never overwritten.

The object is read as \code{\link{cellucid_view}()} shows it: Seurat
\code{meta.data} (and \code{active_ident}), SingleCellExperiment
\code{colData}, reductions as \code{X_<name>}, the neighbour graph, and for
several assays or alternative experiments the modality the session's data
mapping names. Pass \code{\link{cellucid_view}()}'s own arguments in
\code{...} when it was given any that change what is shown.
}
\seealso{
\code{\link{cellucid_request_session}()} (in \code{\link{cellucid_on}}) to
ask a viewer for its session; cellucid-python's \code{apply_cellucid_session_to_anndata()}.
}
\examples{
\dontrun{
v <- cellucid_view(seu)
# ... highlight cells and build fields in the viewer ...
seu <- cellucid_apply_session(seu, v)

# Or from a file saved in the viewer:
seu <- cellucid_apply_session(seu, "pbmc.cellucid-session", expected_dataset_id = "seu")
head(seu[[]])
}
}
