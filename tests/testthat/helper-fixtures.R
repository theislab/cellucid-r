# The fixtures this suite shares byte for byte with the web app and
# cellucid-python (compared with their homes by expect_shared_file(),
# helper-repository.R).

cellucid_fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

cellucid_read_fixture <- function(name) {
  jsonlite::fromJSON(cellucid_fixture_path(name), simplifyVector = FALSE)
}

# The SHA-256 of a file, by the package's own implementation (src/sha256.c,
# checked against the FIPS 180-4 test vectors in test-web-cache.R), on every R
# the package supports; tools::sha256sum() needs R 4.5.
cellucid_file_sha256 <- function(path) {
  cellucid:::.native_sha256(readBin(path, "raw", file.size(path)))
}
