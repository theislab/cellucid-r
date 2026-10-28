% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds the
% \usage block below into a function and compares it to the real formals of
% cellucid_restore_original_order(), so the two cannot drift apart.
\name{cellucid_restore_original_order}
\alias{cellucid_restore_original_order}
\title{Put the cells of a Cellucid-ready file back in the order they were given}
\usage{
cellucid_restore_original_order(x)
}
\arguments{
\item{x}{The cells of a Cellucid-ready file as another package reads them: a
\code{data.frame} (its \code{obs}), a \pkg{Seurat} object, or a
\pkg{SingleCellExperiment} / \pkg{SummarizedExperiment} (for example from
\pkg{anndataR} or \pkg{zellkonverter}).}
}
\value{
\code{x} with its cells in the order they were given to
\code{\link{cellucid_prepare}()} (or Python's \code{cellucid.prepare()}),
without the \code{cellucid_original_position} column and the
\code{cellucid$priority} record; \code{x} itself when its rows are each its
own original position already.
}
\description{
\code{\link{cellucid_prepare}()} stores the cells of a file in priority order
(\code{priority_seed}), so a subsample is the first rows of every array,
records the order in \code{uns['cellucid']['priority']} (the seed and the
cell count it was written for) and keeps each row's original position in the
cell metadata column \code{cellucid_original_position}. A reader of the file
sees the cells in the stored order; this puts them back. Python's
\code{cellucid.restore_original_order()} does the same for an AnnData or
MuData.
}
\details{
The record is where the importers that keep \code{uns} put it:
\code{metadata(x)$cellucid$priority} of a SingleCellExperiment,
\code{misc$cellucid$priority} of a Seurat object, the \code{"cellucid"}
attribute of a data frame. It is verified before it is used, by the rule the
web app and Python share: its seed and specification version are valid, its
cell count is the object's, its levels are those of that count, and the
position column (when there) is the seed's order. An object subset, combined
or reordered after it was read keeps a record that no longer describes it,
and a position column with no record names no order: either is returned as
it is, with a warning, each row its own original position -- as the web app
reads such a file.

A Seurat object is reordered in every part: each assay's layers (or
matrices), its reductions' cell embeddings, its graphs (both axes), its
neighbour indices (rows, and the cells they point at), its cell metadata and
its identities. Images keep their cells by name.
}
\examples{
# Three cells stored in priority order under seed 7 (rows hold the cells 0, 2
# and 1 by original position), as a Cellucid-ready file records them.
positions <- c(0, 2, 1)
obs <- data.frame(cell_type = c("B", "T", "NK"), row.names = c("c1", "c2", "c3"))
obs <- obs[positions + 1, , drop = FALSE]
obs$cellucid_original_position <- positions
attr(obs, "cellucid") <- list(
  priority = list(seed = 7, specVersion = 1, nObs = 3, levels = c(1, 2))
)
cellucid_restore_original_order(obs)
}
