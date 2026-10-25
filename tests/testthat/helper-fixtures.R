# The fixtures this suite shares byte for byte with cellucid-python, and the
# SHA-256 each repository pins for them. A change to either file changes the
# pin in both repositories together, so neither copy can change alone.

cellucid_fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

cellucid_read_fixture <- function(name) {
  jsonlite::fromJSON(cellucid_fixture_path(name), simplifyVector = FALSE)
}

# The SHA-256 of a file, by the package's own implementation (src/sha256.c,
# checked against the FIPS 180-4 test vectors in test-web-cache.R), so the pins
# hold on every R the package supports; tools::sha256sum() needs R 4.5.
cellucid_file_sha256 <- function(path) {
  cellucid:::.native_sha256(readBin(path, "raw", file.size(path)))
}

expect_pinned_fixture <- function(name, sha256) {
  expect_identical(cellucid_file_sha256(cellucid_fixture_path(name)), sha256, info = name)
}
