% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds the
% \usage block below into a function and compares it to the real formals of
% cellucid_transition_drift(), so the two cannot drift apart.
\name{cellucid_transition_drift}
\alias{cellucid_transition_drift}
\title{Compute transition drift in embedding space}
\usage{
cellucid_transition_drift(transition_matrix, embedding, normalize_rows)
}
\arguments{
\item{transition_matrix}{Required cells x cells cell-cell transition matrix
(CellRank's, for example): a \pkg{Matrix} sparse matrix with numeric values,
a dense \pkg{Matrix} matrix, or a base numeric matrix. Row \code{i} holds
the weights of the cells cell \code{i} moves to. Every value must be finite
and non-negative; a logical or pattern matrix is refused.}
\item{embedding}{Required coordinates the drift is computed in, one row per
cell, as \code{cellucid_prepare()} takes an embedding: a numeric matrix, a
numeric \code{data.frame}, a \pkg{Matrix} matrix, or a plain numeric vector
(one column). Every column is used, so the drift is as wide as the embedding.
Every value must be finite and representable in float32, the type the drift
is stored in.}
\item{normalize_rows}{Required \code{TRUE} or \code{FALSE}, with no default,
because only you know whether your matrix is already row-stochastic:
\code{TRUE} divides every row of \code{transition_matrix} by its sum (every
row sum must then be nonzero); \code{FALSE} uses the matrix as given.}
}
\value{
A double matrix of the same dimensions (and dimnames) as \code{embedding}:
the drift at the precision it was computed in. It is rounded to float32 once,
when \code{cellucid_prepare()} writes it.
}
\description{
\code{cellucid_transition_drift()} computes the drift of every cell in the
coordinates of an embedding: the transition-weighted mean position of its
successors minus its own position, \code{T \%*\% E - E}, with the product
divided by the row sums of \code{T} first when \code{normalize_rows = TRUE}.
Pass the result to \code{cellucid_prepare()} as a vector field on that
embedding, and the viewer draws it as an animated flow.
}
\details{
The drift is stored by \code{cellucid_prepare()}, as every vector field is:
\code{vector_fields = list(T_fwd = list(X_umap = drift))} writes
\code{obsm['T_fwd_umap']}, the field \code{T_fwd} over \code{X_umap}; the
embedding must be among \code{embeddings} under that key. Compute one
drift per embedding to draw the field on several
(\code{list(T_fwd = list(X_umap = d2, X_umap_3d = d3))}).

This is the computation of the Python package's
\code{cellucid.compute_transition_drift()}, operation for operation. Python's
\code{cellucid.add_transition_drift_to_obsm()}, which writes the result into
the \code{obsm} of an AnnData held in memory, has no R counterpart: this
package holds no AnnData in memory, and \code{cellucid_prepare()}'s
\code{vector_fields} names the \code{obsm} key and checks it by the same
vector-field key rule.

Both packages are held to one shared conformance fixture
(\code{tests/testthat/fixtures/transition_drift_conformance.json}, computed
with exact rational arithmetic): bit for bit where the product and the row
sums are exact in double precision, and otherwise within the forward-error
bound of the computation, since the order in which a BLAS, a sparse kernel or
\code{rowSums()} adds terms is not fixed.

The function refuses, with the reason: a matrix that is not cells x cells, a
negative, missing or infinite value, a zero row sum with
\code{normalize_rows = TRUE}, and an embedding value or a computed drift
outside the float32 range, including a nonzero magnitude so small that
float32 would store it as zero.
}
\examples{
umap <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
# Cell 1 moves to cells 2 and 3, cell 2 to cell 4, cells 3 and 4 stay.
transition <- matrix(
  c(0, 1, 1, 0,
    0, 0, 0, 1,
    0, 0, 1, 0,
    0, 0, 0, 1),
  nrow = 4, byrow = TRUE
)
drift <- cellucid_transition_drift(transition, umap, normalize_rows = TRUE)
drift

cellucid_prepare(
  out_path = file.path(tempdir(), "drift.zarr"),
  format = "zarr",
  embeddings = list(X_umap = umap),
  dataset_id = "cellucid-drift-example",
  dataset_name = "Cellucid drift example",
  vector_fields = list(T_fwd = list(X_umap = drift)),
  force = TRUE
)
}
