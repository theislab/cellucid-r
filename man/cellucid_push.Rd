% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_push}
\alias{cellucid_push}
\alias{cellucid_write_gmt}
\alias{cellucid_read_gmt}
\title{Show a column in Cellucid, and read or write GMT files}
\usage{
cellucid_push(viewer, values, key = NULL)

cellucid_write_gmt(gene_sets, file, descriptions = NULL)

cellucid_read_gmt(file)
}
\arguments{
\item{viewer}{A viewer \code{\link{cellucid_view}()} or
\code{\link{cellucid_serve}()} returned.}
\item{values}{The column: one value per cell in order, or a vector named by
cell (some or all cells; the others get no value). A factor, characters or
\code{TRUE}/\code{FALSE} become a categorical field (\code{NA}: no category;
a logical column's labels are \code{"false"} and \code{"true"}, as
cellucid-python pushes a bool column); numbers a continuous field (\code{NA}:
no value).}
\item{key}{The field's name in the viewer. Pushing the same key again replaces
it; a field of that name a notebook did not push is never replaced: the
viewer refuses it.}
\item{gene_sets}{A named list: each set's name and its genes.}
\item{file}{The GMT file to write or read.}
\item{descriptions}{A named list of descriptions by set name; empty when not
given.}
}
\value{
\code{cellucid_push()}: the request of the viewer's \code{field.push} action
(\code{\link{cellucid_act}()}), at once and invisibly: the column is sent now, or when the
viewer's channel opens, so a notebook run top to bottom is fine.
\code{\link{cellucid_action_result}(request)} returns the field's name, kind
and id once the viewer added it (one step of its history, labelled as R's),
or stops with the viewer's reason when it refused the column (one that does
not fit the dataset); \code{\link{cellucid_results}()} lists it once the
viewer added it (in a notebook, from the next cell).
\code{cellucid_write_gmt()}: the file's path, invisibly.
\code{cellucid_read_gmt()}: a named list of character vectors.
}
\description{
\code{cellucid_push()} shows a column -- clusters, labels, a score -- in the
viewer as a field to colour, filter and analyse by: the viewer's
\code{field.push} action (\code{\link{cellucid_field_push}()}) with its
parameters made from the vector, its categories and codes (as uint8, uint16
or uint32) or its numbers (as float32), its cells by original position. It
does not wait: a refusal reaches \code{cellucid_action_result()}, or
\code{cellucid_on(viewer, "command_error", ...)}, or a warning. When the
viewer opened a subsample, the cells it holds take their values. A gene set
is the viewer's \code{gene-set.save} action:
\code{\link{cellucid_gene_set_save}()}. A GMT file (MSigDB) holds
one set per line, \code{name<TAB>description<TAB>gene...}: what the web app
exports and imports, and what cellucid-python's \code{write_gmt()} and
\code{read_gmt()} write and read, by one rule: each run of tabs and line
breaks in a name, description or gene is written as one space; a line
without a name reads as \code{"Set <n>"}; a name used again is
\code{"name (2)"}, \code{"name (3)"}...; each gene is kept once.
}
\examples{
sets <- list(`T cells` = c("CD3E", "CD3D"), `B cells` = "MS4A1")
path <- tempfile(fileext = ".gmt")
cellucid_write_gmt(sets, path)
identical(cellucid_read_gmt(path), sets)
if (interactive() && requireNamespace("SeuratObject", quietly = TRUE)) {
  data("pbmc_small", package = "SeuratObject")
  object <- SeuratObject::UpdateSeuratObject(pbmc_small)
  viewer <- cellucid_view(object)
  cellucid_push(viewer, object$RNA_snn_res.1, key = "clusters")
  cellucid_gene_set_save(viewer, "T cells", c("CD3E", "CD3D"))
  cellucid_stop(viewer)
}
}
