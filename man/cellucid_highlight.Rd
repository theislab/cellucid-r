% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_highlight}
\alias{cellucid_highlight}
\alias{cellucid_clear_highlights}
\alias{cellucid_color_by}
\alias{cellucid_set_visibility}
\alias{cellucid_reset_camera}
\alias{cellucid_freeze}
\alias{cellucid_set_data_mapping}
\title{Send commands to an open Cellucid viewer}
\usage{
cellucid_highlight(viewer, cells, color = "#e41a1c", label = NULL)

cellucid_clear_highlights(viewer)

cellucid_color_by(viewer, field)

cellucid_set_visibility(viewer, cells = NULL, visible = TRUE)

cellucid_reset_camera(viewer)

cellucid_freeze(viewer)

cellucid_set_data_mapping(viewer, ...)
}
\arguments{
\item{viewer}{A viewer \code{\link{cellucid_view}()} or
\code{\link{cellucid_serve}()} returned.}
\item{cells}{The cells: their names (an object's cell names), their 1-based
positions, or \code{TRUE}/\code{FALSE} for every cell; at least one. For
\code{cellucid_set_visibility()}, \code{NULL} means every cell.}
\item{color}{A six-digit hex colour, such as \code{"#e41a1c"}.}
\item{label}{The highlight group's name in the viewer; \code{NULL} names it
by its cell count. \code{\link{cellucid_pull}()} brings the group back as
the column \code{cellucid_highlight__<group id>}.}
\item{field}{The \code{obs} column (or gene) to colour the cells by.}
\item{visible}{\code{TRUE} shows the cells, \code{FALSE} hides them.}
\item{...}{The data-mapping fields to change, by name, such as
\code{modality = "ADT"} (a multimodal object) or \code{knn = "connectivities"};
a list value is sent as a JSON array or object. An object's store holds one
matrix per modality (its \code{data}, else \code{counts}) and no layers, so
\code{expression} stays \code{"X"} for an object; a file may hold layers
(\code{expression = "layer:counts"}).}
}
\value{
The viewer, invisibly; \code{cellucid_set_data_mapping()} \code{NULL},
invisibly.
}
\description{
Commands go to the viewer over its channel, in order. A command sent before the
viewer has opened its page waits and goes out as soon as it connects.
\code{cellucid_highlight()} and \code{cellucid_set_visibility()} send their
cells as one binary frame, four bytes a cell, so a highlight of millions of
cells is a few megabytes written in one pass.
\code{cellucid_set_data_mapping()} waits for nothing (a Jupyter notebook hears
the viewer only between cells): the merged mapping is checked and sent at once
when the viewer has reported its mapping (an invalid value stops here, the
viewer untouched), or as soon as it reports one, after its dataset opens (so a
notebook run top to bottom is fine). The viewer's answer becomes
\code{viewer$data_mapping}; a mapping it refuses arrives as a
\code{command_error} (a \code{\link{cellucid_on}()} hook, or a warning when
none handles it).
\code{cellucid_clear_highlights()} removes every highlight group on every
highlight page; the pages themselves stay.
}
\seealso{
\code{\link{cellucid_on}()} for what the viewer sends back.
}
\examples{
if (interactive() && requireNamespace("SeuratObject", quietly = TRUE)) {
  data("pbmc_small", package = "SeuratObject")
  viewer <- cellucid_view(SeuratObject::UpdateSeuratObject(pbmc_small))
  cellucid_color_by(viewer, "groups")
  cellucid_highlight(viewer, 1:20, color = "#377eb8")
  cellucid_set_visibility(viewer, pbmc_small$groups == "g2", visible = FALSE)
  cellucid_stop(viewer)
}
}
