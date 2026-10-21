% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds the
% \usage block below into a function and compares it to the real formals of
% cellucid_embedding(), so the two cannot drift apart.
\name{cellucid_embedding}
\alias{cellucid_embedding}
\title{One embedding to export, with the columns it is read from}
\usage{
cellucid_embedding(
  values,
  columns = NULL,
  label = NULL,
  source_key = NULL
)
}
\arguments{
\item{values}{Numeric matrix (or numeric \code{data.frame}, or \pkg{Matrix}) with
one row per cell, or a plain numeric vector for a single column.}
\item{columns}{Optional numeric vector of the 1-based columns to plot, one per
axis, so its length is the embedding's dimension (1, 2, or 3). The columns must
be distinct whole numbers within the array. Required when \code{values} has more
than three columns; otherwise every column is plotted. The export records the
columns 0-based, as \code{cellucid-python} does.}
\item{label}{Optional text the viewer shows for this embedding, such as
\code{"t-SNE"}. It is drawn verbatim, so it follows the same display-text rule
as a category label. Defaults to the embedding's id.}
\item{source_key}{Optional non-empty string naming where the array came from,
such as \code{"X_pca"}. It is recorded as provenance in
\code{dataset_identity.json} and changes nothing else.}
}
\value{
An object of class \code{"cellucid_embedding"} to pass as one element of the
\code{embeddings} list of \code{\link{cellucid_prepare}()}.
}
\description{
\code{cellucid_embedding()} describes one embedding for
\code{\link{cellucid_prepare}()} when a plain matrix is not enough: a wide array
read at named columns, a display label, or a recorded source key. It mirrors
\code{cellucid.Embedding} in the \pkg{cellucid} Python package.
}
\examples{
pca <- matrix(c(0, 3, 1, 2, 5, 1, 4, 0, 2, 1, 1, 1), ncol = 4, byrow = TRUE)
cellucid_embedding(pca, columns = c(1, 2), label = "PCA", source_key = "X_pca")
}
