# DESCRIPTION is the one in-repo source of the package version. Every other
# in-repo site (NEWS.md, CITATION.cff, README.md, vignettes/installation.Rmd)
# is checked against it in test-current-contract.R, so a partial bump fails.
#
# The root this reads from is resolved in helper-repository.R, which finds the
# real sources in both the git checkout and an R CMD check directory. Reading
# the version from anywhere else would let DESCRIPTION and the pages compared
# against it drift apart unnoticed.

cellucid_description_version <- function() {
  description_path <- file.path(cellucid_repository_root(), "DESCRIPTION")
  unname(read.dcf(description_path)[1L, "Version"])
}

cellucid_semantic_versions <- function(path) {
  text <- paste(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  # The generated citation region is not a version statement: once the
  # preprint is published its DOI (10.1101/YYYY.MM.DD.NNNNNN) holds digit
  # triples that are not package versions. test-citation-contract.R owns it.
  text <- gsub(
    "<!-- BEGIN cellucid-citation[\\s\\S]*?<!-- END cellucid-citation -->",
    "",
    text,
    perl = TRUE
  )
  unique(unlist(regmatches(
    text,
    gregexpr("[0-9]+\\.[0-9]+\\.[0-9]+", text)
  )))
}
