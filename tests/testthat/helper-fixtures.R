# The fixtures this suite shares byte for byte with cellucid-python, and the
# SHA-256 each repository pins for them. A change to either file changes the
# pin in both repositories together, so neither copy can change alone.

cellucid_fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

cellucid_read_fixture <- function(name) {
  jsonlite::fromJSON(cellucid_fixture_path(name), simplifyVector = FALSE)
}

# tools::sha256sum() exists from R 4.5; on an older R the digest is computed by
# the system's shasum or sha256sum, and NA means no SHA-256 tool is available.
cellucid_file_sha256 <- function(path) {
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE)) {
    return(unname(tools::sha256sum(path)))
  }
  for (tool in c("shasum", "sha256sum")) {
    binary <- Sys.which(tool)
    if (nzchar(binary)) {
      arguments <- if (identical(tool, "shasum")) {
        c("-a", "256", shQuote(path))
      } else {
        shQuote(path)
      }
      output <- system2(binary, arguments, stdout = TRUE)
      return(strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]])
    }
  }
  NA_character_
}

expect_pinned_fixture <- function(name, sha256) {
  digest <- cellucid_file_sha256(cellucid_fixture_path(name))
  skip_if(is.na(digest), "no SHA-256 tool on this system")
  expect_identical(digest, sha256, info = name)
}
