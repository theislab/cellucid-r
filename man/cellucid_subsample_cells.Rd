% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds the
% \usage block below into a function and compares it to the real formals of
% cellucid_subsample_cells(), so the two cannot drift apart.
\name{cellucid_subsample_cells}
\alias{cellucid_subsample_cells}
\title{The cells a Cellucid subsample keeps}
\usage{
cellucid_subsample_cells(x, n_cells = NULL, fraction = NULL, method = "uniform", seed = 0,
  groupby = NULL, allocation = "proportional", min_per_group = 0,
  block_size = NULL, basis = NULL, columns = c(1, 2), grid_size = 32)
}
\arguments{
\item{x}{A \pkg{Seurat} object, a \pkg{SingleCellExperiment} or
\pkg{SummarizedExperiment}, a \code{data.frame} of cells, or the cell count
itself (enough for \code{"uniform"} and \code{"blocks"}).}
\item{n_cells}{How many cells to keep; or give \code{fraction}. A count of
every cell or more keeps every cell.}
\item{fraction}{The share of cells to keep, in (0, 1]; the count is
\code{n * fraction} rounded half up, at least 1, computed exactly from the
fraction's shortest decimal text, never as a floating-point product, as the
web app and Python compute it. Any count from 1 is a
selection here; the web app opens a subsample of at least 2 cells, since one
cell has no coordinate variation to draw.}
\item{method}{\code{"uniform"} (every cell alike), \code{"blocks"} (runs of
\code{block_size} neighbouring cells in file order, kept together),
\code{"stratified"} (split over the categories of \code{groupby}) or
\code{"spatial"} (split over the occupied boxes of a \code{grid_size} grid on
\code{basis}, balanced, so sparse regions keep their cells).}
\item{seed}{A whole number in [0, 2^32 - 1]. The same seed keeps the same
cells; with \code{"uniform"} and \code{"blocks"} a smaller subsample is inside
a larger one.}
\item{groupby}{\code{"stratified"}: a cell metadata column name
(\code{meta.data}, \code{colData}), or one value per cell (a factor keeps its
level order). Cells without a value are one more stratum.}
\item{allocation}{\code{"stratified"}: \code{"proportional"} (largest
remainder) or \code{"balanced"} (equal shares, water-filled).}
\item{min_per_group}{\code{"stratified"}: each category first gets this many
cells (or all it has).}
\item{block_size}{\code{"blocks"}: cells per block. The web app derives it
from the file's chunk layout and shows it; give that value to reproduce its
cells.}
\item{basis}{\code{"spatial"}: a reduction name (a Seurat reduction, a
\code{reducedDim}), or a matrix of coordinates, one row per cell.}
\item{columns}{\code{"spatial"}: the embedding columns used, one to three.}
\item{grid_size}{\code{"spatial"}: boxes per axis.}
}
\value{
The kept cells' positions (1-based, ascending): an integer vector, or a
double one for an object of more cells than an R integer counts.
}
\description{
\code{cellucid_subsample_cells()} gives the cells a subsample keeps: the cells
the Cellucid web app opens for the same request (Open data, Cells to open),
and the cells \code{cellucid.subsample_cells()} gives in Python. Use it to
take the cells the viewer showed, or to choose cells the same way before
opening.
}
\details{
All three follow one subsampling specification: a cell's rank is a 32-bit
hash of the seed and its position, and the lowest ranks are kept (within each
stratum, for \code{"stratified"} and \code{"spatial"}). The computation is
native code (exact 32- and 64-bit integer arithmetic, and the spatial grid in
binary64 operations that are never fused), held with the web app and the
Python package to one shared conformance fixture
(\code{tests/testthat/fixtures/subsample_conformance.json}).
}
\examples{
kept <- cellucid_subsample_cells(10000, fraction = 0.1, seed = 0)
length(kept)
}
