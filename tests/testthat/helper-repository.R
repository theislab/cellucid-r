# Where the repository is, and what must be in it.
#
# Two layouts run this suite, and they are told apart by which one matches
# rather than by probing for one file that might simply have been deleted:
#
#   testthat::test_local()    tests/testthat/../.. is the git checkout.
#   R CMD check <tarball>     tests/testthat/../.. is <package>.Rcheck, and
#                             the unpacked sources are one level further down,
#                             in 00_pkg_src/cellucid.
#
# Resolving the root wrongly is not a visible failure: it makes every
# path-based assertion in test-current-contract.R disappear instead. Under
# R CMD check the root used to resolve to the .Rcheck directory, where none of
# DESCRIPTION, NEWS.md, README.md, NAMESPACE, man/ or R/ exists, and forty-two
# assertions about the release identity stopped running while the suite stayed
# green. An unrecognised layout is therefore an error here, never a quiet skip.
#
# The checkout holds every file. The tarball holds only what .Rbuildignore did
# not strip, so the pages that carry release and continuous-integration policy
# are genuinely absent there. Which of the two a path belongs to is declared in
# cellucid_required_paths() and asserted in both directions, so a deleted file
# fails by name instead of turning its assertions into no-ops.

cellucid_test_layout <- function() {
  base <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  if (file.exists(file.path(base, ".Rbuildignore"))) {
    return(list(root = base, checkout = TRUE))
  }
  unpacked <- file.path(base, "00_pkg_src", "cellucid")
  if (file.exists(file.path(unpacked, "DESCRIPTION"))) {
    return(list(
      root = normalizePath(unpacked, winslash = "/", mustWork = TRUE),
      checkout = FALSE
    ))
  }
  stop(
    "The cellucid test suite does not recognise this layout: it expects ",
    "either a git checkout, marked by an .Rbuildignore beside DESCRIPTION, ",
    "or an R CMD check directory holding 00_pkg_src/cellucid/DESCRIPTION, ",
    "and found neither under ",
    base,
    ".",
    call. = FALSE
  )
}

cellucid_repository_root <- function() {
  cellucid_test_layout()$root
}

cellucid_is_source_checkout <- function() {
  cellucid_test_layout()$checkout
}

# Every path test-current-contract.R reads, and which layout holds it. "both"
# means the file survives .Rbuildignore and ships inside the tarball too;
# "checkout" means .Rbuildignore strips it and only the git checkout has it.
cellucid_required_paths <- function() {
  specification <- matrix(
    c(
      "DESCRIPTION", "file", "both",
      "NAMESPACE", "file", "both",
      "NEWS.md", "file", "both",
      "README.md", "file", "both",
      "R", "directory", "both",
      "R/cellucid_prepare.R", "file", "both",
      "man", "directory", "both",
      "man/cellucid_prepare.Rd", "file", "both",
      "inst/CITATION", "file", "both",
      "inst/action_catalog.json", "file", "both",
      "vignettes", "directory", "both",
      "vignettes/installation.Rmd", "file", "both",
      ".Rbuildignore", "file", "checkout",
      "CITATION.cff", "file", "checkout",
      "publishing.md", "file", "checkout",
      "_pkgdown.yml", "file", "checkout",
      "tools/generate-actions.R", "file", "checkout",
      ".github/workflows/R-CMD-check.yaml", "file", "checkout",
      ".github/workflows/pkgdown.yaml", "file", "checkout",
      ".github/workflows/release.yaml", "file", "checkout"
    ),
    ncol = 3L,
    byrow = TRUE,
    dimnames = list(NULL, c("path", "kind", "layout"))
  )
  as.data.frame(specification, stringsAsFactors = FALSE)
}

# The web app's checkout: CELLUCID_WEB_ROOT when set (a worktree of it), else
# the sibling `cellucid` beside this repository. The sibling-copy tests compare
# the shared files (catalogue, history model, golden history, viewer protocol)
# with it, as the web suite reads this package through CELLUCID_R_ROOT.
cellucid_web_root <- function() {
  root <- Sys.getenv("CELLUCID_WEB_ROOT", "")
  if (nzchar(root)) {
    if (!dir.exists(file.path(root, "assets", "js"))) stop("CELLUCID_WEB_ROOT (", root, ") is not a checkout of the web app (it has no assets/js).")
    return(root)
  }
  file.path(cellucid_repository_root(), "..", "cellucid")
}

# A file of a sibling checkout ("cellucid" for the web app, as
# cellucid_web_root() finds it, or "cellucid-python"): the test is skipped,
# saying so, only when that checkout is absent (this repository checked out
# alone, as in its own CI); a checkout present that lacks the file fails it (a
# shared file renamed or deleted on one side is a divergence, never a skip).
cellucid_sibling_file <- function(repository, ...) {
  root <- if (identical(repository, "cellucid")) cellucid_web_root() else file.path(cellucid_repository_root(), "..", repository)
  present <- if (identical(repository, "cellucid")) dir.exists(file.path(root, "assets", "js")) else dir.exists(root)
  if (!present) testthat::skip(paste0(repository, " is not checked out beside this repository (CELLUCID_WEB_ROOT names the web app's)."))
  path <- file.path(root, ...)
  if (!file.exists(path)) testthat::fail(paste0("The ", repository, " checkout (", root, ") has no ", file.path(...), "."))
  path
}
