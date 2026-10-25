% Hand-written and authoritative: this package runs no documentation generator.
% Edit this file directly. tests/testthat/test-current-contract.R rebuilds
% every call of the \usage block below into a function and compares it to the
% real formals, so the two cannot drift apart.
\name{cellucid_clear_web_cache}
\alias{cellucid_clear_web_cache}
\title{Remove the cached Cellucid web app}
\usage{
cellucid_clear_web_cache(cache_dir = NULL)
}
\arguments{
\item{cache_dir}{The cache folder; \code{tools::R_user_dir("cellucid",
"cache")/web} when \code{NULL}. It is removed whole, so it must be a folder
of its own: the home folder, its parent, the working directory and the file
system root are refused.}
}
\value{
The folder removed, invisibly.
}
\description{
\code{\link{cellucid_view}()} downloads the Cellucid web app at first use into
a cache folder and checks every file against the published inventory's SHA-256
digests. Clearing the cache makes the next viewer download it again.
}
\examples{
cellucid_clear_web_cache(file.path(tempdir(), "cellucid-web-example"))
}
