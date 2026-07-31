# DESCRIPTION is the one in-repo source of the package version. Every other
# in-repo site (NEWS.md, CITATION.cff, README.md, vignettes/installation.Rmd)
# is checked against it in test-current-contract.R, so a partial bump fails.

cellucid_repository_root <- function() {
  normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = FALSE
  )
}

cellucid_description_version <- function() {
  description_path <- file.path(cellucid_repository_root(), "DESCRIPTION")
  if (file.exists(description_path)) {
    unname(read.dcf(description_path)[1L, "Version"])
  } else {
    as.character(utils::packageVersion("cellucid"))
  }
}

cellucid_semantic_versions <- function(path) {
  text <- paste(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  unique(unlist(regmatches(
    text,
    gregexpr("[0-9]+\\.[0-9]+\\.[0-9]+", text)
  )))
}
