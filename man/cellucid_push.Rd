% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_push}
\alias{cellucid_push}
\alias{cellucid_push_gene_set}
\alias{cellucid_write_gmt}
\alias{cellucid_read_gmt}
\title{Show a column or a gene set in Cellucid, and read or write GMT files}
\usage{
cellucid_push(viewer, values, key = NULL)

cellucid_push_gene_set(viewer, name, genes, description = "")

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
it; a field of that name this R session did not push is never replaced: the
viewer refuses it with a \code{command_error} event.}
\item{name}{The gene set's name; a set of the same name pushed before is
replaced.}
\item{genes}{The genes, in order, by identifier (\code{rownames(object)}) or
by the symbol the viewer shows (identifiers first); each is kept once.}
\item{description}{The gene set's description.}
\item{gene_sets}{A named list: each set's name and its genes.}
\item{file}{The GMT file to write or read.}
\item{descriptions}{A named list of descriptions by set name; empty when not
given.}
}
\value{
\code{cellucid_push()} and \code{cellucid_push_gene_set()}: the viewer,
invisibly; the column or set is sent at once -- or, for a format only the
browser reads, the column once the viewer has said how many cells it has (its
\code{ready}), so a notebook run top to bottom is fine -- and
\code{\link{cellucid_results}()} lists it once the viewer added it (in a
notebook, from the next cell). A column that does not fit the dataset once
its cell count is known arrives as a \code{command_error}. \code{cellucid_write_gmt()}: the file's path, invisibly.
\code{cellucid_read_gmt()}: a named list of character vectors.
}
\description{
\code{cellucid_push()} shows a column -- clusters, labels, a score -- in the
viewer as a field to colour, filter and analyse by. \code{cellucid_push_gene_set()}
gives the viewer a gene set, offered by Gene signature and Dot plot ("From
notebook") and saved with the session. Neither waits: a refusal reaches
\code{cellucid_on(viewer, "command_error", ...)}. When the viewer opened a
subsample, the cells it holds take their values. A GMT file (MSigDB) holds
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
  cellucid_push_gene_set(viewer, "T cells", c("CD3E", "CD3D"))
  cellucid_stop(viewer)
}
}
