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
# R CMD check a root taken to be the .Rcheck directory holds none of
# DESCRIPTION, NEWS.md, README.md, NAMESPACE, man/ or R/, and every assertion
# about the release identity would stop running while the suite stays green.
# An unrecognized layout is therefore an error here, never a quiet skip.
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
    "The cellucid test suite does not recognize this layout: it expects ",
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

# The sibling repositories' checkouts: the web app ("cellucid") at
# CELLUCID_WEB_ROOT (a worktree of it) and cellucid-python at
# CELLUCID_PYTHON_ROOT when set, else the checkout beside this repository.
# The tests compare this package's copy of every shared file with them
# (expect_shared_file() below), as the web suite reads this package through
# CELLUCID_R_ROOT. A named folder that is there must be a checkout of that
# repository. A sibling that is not there (an installed package, a CRAN check,
# a variable naming a folder that does not exist) returns NULL, unless
# CELLUCID_REQUIRE_SIBLINGS is 1 (this repository's CI fetches both siblings
# at main and sets it, also under R CMD check, where the sources sit in
# <package>.Rcheck and nothing is beside them): then it is an error.
cellucid_siblings_required <- function() identical(Sys.getenv("CELLUCID_REQUIRE_SIBLINGS", ""), "1")

cellucid_sibling_variable <- function(repository) {
  switch(repository, cellucid = "CELLUCID_WEB_ROOT", `cellucid-python` = "CELLUCID_PYTHON_ROOT",
         stop("No sibling repository is called ", repository, ".", call. = FALSE))
}

# Where `repository` was looked for and not found, as one sentence.
cellucid_sibling_not_there <- function(repository) {
  variable <- cellucid_sibling_variable(repository)
  named <- Sys.getenv(variable, "")
  if (nzchar(named)) return(paste0(repository, " is not there: ", variable, " names ", named, ", which does not exist."))
  paste0(repository, " is not checked out beside this repository (set ", variable, " to name its checkout).")
}

cellucid_sibling_root <- function(repository) {
  variable <- cellucid_sibling_variable(repository)
  marker <- if (identical(repository, "cellucid")) file.path("assets", "js") else "pyproject.toml"
  named <- Sys.getenv(variable, "")
  root <- if (nzchar(named)) named else file.path(cellucid_repository_root(), "..", repository)
  if (file.exists(file.path(root, marker))) return(root)
  if (nzchar(named) && file.exists(root)) {
    stop(variable, " (", named, ") is not a checkout of ", repository, " (it has no ", marker, ").", call. = FALSE)
  }
  if (cellucid_siblings_required()) {
    stop("CELLUCID_REQUIRE_SIBLINGS is 1, and ", cellucid_sibling_not_there(repository), call. = FALSE)
  }
  NULL
}

# A file of a sibling checkout ("cellucid" for the web app, or
# "cellucid-python"): the test is skipped, saying so, only when that checkout
# is absent and not required; a checkout present that lacks the file fails it
# (a shared file renamed or deleted on one side is a divergence, never a skip).
cellucid_sibling_file <- function(repository, ...) {
  root <- cellucid_sibling_root(repository)
  if (is.null(root)) testthat::skip(cellucid_sibling_not_there(repository))
  path <- file.path(root, ...)
  if (!file.exists(path)) testthat::fail(paste0("The ", repository, " checkout (", root, ") has no ", file.path(...), "."))
  path
}

# The one command that makes every copy of every shared file agree, run in
# the web repository: it writes the generated files again and copies every
# home into the packages.
cellucid_sync_shared_files <- "node scripts/sync-shared-files.mjs"

# This package's `ours` is byte for byte `repository`'s `theirs` (both named
# from their repository's root): the web app's file for a file whose home is
# the web repository, cellucid-python's for one only the two packages share.
# No digest is written anywhere. Skipped, saying so, without that checkout,
# unless the siblings are required (cellucid_sibling_file()): so it is a test
# of its own, never the first expectation of a test that goes on. `then` names
# what to run after the copies agree again, for a file generated code is
# written from.
expect_shared_file <- function(ours, theirs, repository = "cellucid", then = "") {
  other <- cellucid_sibling_file(repository, theirs)
  bytes <- function(path) readBin(path, "raw", file.size(path))
  testthat::expect(
    identical(bytes(file.path(cellucid_repository_root(), ours)), bytes(other)),
    paste0(
      ours, " is not byte for byte ", repository, "'s ", theirs, ". Run `", cellucid_sync_shared_files,
      "` in the web repository (cellucid): it makes every copy of every shared file agree", then, "."
    )
  )
  invisible(NULL)
}

# The `node` that runs the web app's code for `purpose`: skipped without Node,
# unless the siblings are required, where its absence is an error.
cellucid_node <- function(purpose) {
  node <- Sys.which("node")
  if (nzchar(node)) return(node)
  if (cellucid_siblings_required()) {
    stop("CELLUCID_REQUIRE_SIBLINGS is 1 but Node is not installed: ", purpose, " cannot run.", call. = FALSE)
  }
  testthat::skip(paste0("Node is not installed: ", purpose, " cannot run."))
}

# What the web app's generator `script` renders now from its source
# (renderHistoryModel from app/history/model.js, renderActionCatalog and
# renderMcpTools from the action definitions): a shared file this package holds
# is compared with the source itself, not only with the web app's committed
# copy, so a source changed without regenerating the copies fails here too
# (cellucid-python's tests/_web_root.py rendered_by_web()). Skipped without a
# web checkout or Node, unless the siblings are required.
cellucid_rendered_by_web <- function(script, render) {
  path <- cellucid_sibling_file("cellucid", script)
  node <- cellucid_node("the web app's generator")
  # The first argument is a word, never the script: the script runs its own
  # main() when argv[1] names it.
  code <- paste0(
    "import { pathToFileURL } from 'node:url';",
    "const [script, name] = process.argv.slice(2);",
    "const module = await import(pathToFileURL(script).href);",
    "process.stdout.write(module[name]());"
  )
  output <- tempfile("cellucid_rendered_")
  on.exit(unlink(output), add = TRUE)
  status <- system2(node, c("--input-type=module", "-e", shQuote(code), "render", shQuote(normalizePath(path)), render),
                    stdout = output, stderr = output)
  text <- readChar(output, file.size(output), useBytes = TRUE)
  if (!identical(status, 0L)) testthat::fail(paste0(script, " ", render, "() failed: ", text))
  text
}
