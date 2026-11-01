# The web app the viewer serves: the published build, downloaded once into a
# cache and checked file by file with cellucid-python's rules
# (web_cache.py): the inventory's exact schema, the build id stamped in
# index.html, every file's length and SHA-256. Downloads come from a local
# server of the fixture build (python/web_source.py), never the internet.

fixture_web <- function() testthat::test_path("fixtures", "web")
inventory_bytes <- function() {
  path <- file.path(fixture_web(), "cellucid-web-assets.json")
  readBin(path, "raw", file.size(path))
}
inventory_with <- function(change) {
  inventory <- jsonlite::fromJSON(rawToChar(inventory_bytes()), simplifyVector = FALSE)
  charToRaw(as.character(jsonlite::toJSON(change(inventory), auto_unbox = TRUE)))
}

test_that("SHA-256 is computed as FIPS 180-4 defines it", {
  sha <- cellucid:::.native_sha256
  expect_identical(sha(raw(0)), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")
  expect_identical(sha(charToRaw("abc")), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
  expect_identical(
    sha(charToRaw("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")),
    "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
  )
  expect_identical(sha(as.raw(rep(0x61, 1e6))), "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0")
  # The fixture inventory records index.html's SHA-256 as `shasum -a 256` gives it.
  recorded <- jsonlite::fromJSON(rawToChar(inventory_bytes()), simplifyVector = FALSE)$assets[[2L]]
  expect_identical(recorded$path, "index.html")
  expect_identical(sha(readBin(file.path(fixture_web(), "index.html"), "raw", 1e5)), recorded$sha256)
})

test_that("the inventory is parsed with its one exact schema", {
  parse <- cellucid:::.parse_web_inventory
  inventory <- parse(inventory_bytes())
  expect_identical(inventory$build_id, "fixture-build-1")
  expect_identical(names(inventory$assets), c("assets/app.js", "index.html"))
  expect_error(parse(inventory_with(function(x) { x$version <- 2L; x })), "version must be integer 1")
  expect_error(parse(inventory_with(function(x) { x$extra <- 1L; x })), "exactly version, build_id, and assets")
  expect_error(parse(inventory_with(function(x) { x$assets <- rev(x$assets); x })), "must be sorted")
  expect_error(parse(inventory_with(function(x) { x$assets <- x$assets[1]; x })), "must contain index.html")
  expect_error(parse(inventory_with(function(x) { x$assets[[1]]$path <- "../app.js"; x })), "empty or dot segment")
  expect_error(parse(inventory_with(function(x) { x$assets[[1]]$path <- "assets/CON.js"; x })), "not portable")
  expect_error(parse(inventory_with(function(x) { x$assets[[1]]$sha256 <- "abc"; x })), "invalid sha256")
  expect_error(parse(inventory_with(function(x) { x$assets[[1]]$content_type <- "Text/JS"; x })), "MIME type")
  expect_error(parse(charToRaw("{")), "valid JSON")
})

test_that("the WebAssembly decoders are listed with the type compileStreaming needs", {
  # The web app's inventory lists assets/wasm/codecs.wasm as application/wasm;
  # the viewer serves each file with the type its inventory records.
  inventory <- cellucid:::.parse_web_inventory(inventory_with(function(x) {
    wasm <- list(path = "assets/wasm/codecs.wasm", sha256 = strrep("ab", 32L), bytes = 8L, content_type = "application/wasm")
    x$assets <- c(x$assets[1L], list(wasm), x$assets[-1L])
    x
  }))
  expect_identical(names(inventory$assets), c("assets/app.js", "assets/wasm/codecs.wasm", "index.html"))
  expect_identical(inventory$assets[["assets/wasm/codecs.wasm"]]$content_type, "application/wasm")
})

test_that("index.html carries exactly one build id", {
  id <- cellucid:::.web_build_id_of
  expect_identical(id(charToRaw("<meta name=\"cellucid-web-build-id\" content=\"b1\">")), "b1")
  expect_error(id(charToRaw("<html></html>")), "exactly one cellucid-web-build-id")
  expect_error(id(charToRaw(strrep("<meta name='cellucid-web-build-id' content='b'>", 2))), "exactly one")
  expect_error(id(charToRaw("<meta name=\"cellucid-web-build-id\" content=\" b\">")), "invalid")
})

test_that("a cache is served only when every file checks out", {
  cache <- viewer_fixture_web_cache()
  inventory <- cellucid:::.verify_web_cache(cache)
  expect_identical(inventory$build_id, "fixture-build-1")
  writeLines("extra", file.path(cache, "stray.txt"))
  expect_error(cellucid:::.verify_web_cache(cache), "unexpected c\\(\"stray.txt\"\\)")
  unlink(file.path(cache, "stray.txt"))
  writeLines("changed", file.path(cache, "assets", "app.js"))
  expect_error(cellucid:::.verify_web_cache(cache), "byte length is|sha256 is")
})

test_that("the source must be https, or http on this machine only", {
  expect_error(cellucid:::.require_source_url("http://example.org"), "must use https")
  expect_error(cellucid:::.require_source_url("https://example.org/"), "must not end")
  expect_error(cellucid:::.require_source_url("https://user@example.org"), "no\\s+credentials")
  expect_identical(cellucid:::.require_source_url("http://127.0.0.1:4173"), "http://127.0.0.1:4173")
})

test_that("a build is downloaded, checked, kept, and used again offline", {
  skip_on_cran()
  python <- Sys.which("python3")
  skip_if(!nzchar(python), "no python3 to serve the fixture build")
  source <- tempfile("cellucid_source_")
  dir.create(source)
  file.copy(list.files(fixture_web(), full.names = TRUE), source, recursive = TRUE)
  port_file <- tempfile()
  system2(python, c(shQuote(testthat::test_path("python", "web_source.py")), shQuote(source), shQuote(port_file)), wait = FALSE)
  on.exit(writeLines("", file.path(source, ".stop")), add = TRUE)
  deadline <- Sys.time() + 20
  while (!file.exists(port_file) && Sys.time() < deadline) Sys.sleep(0.05)
  skip_if(!file.exists(port_file), "the fixture web source did not start")
  url <- paste0("http://127.0.0.1:", readLines(port_file, warn = FALSE))

  cache <- file.path(tempfile("cellucid_cache_"), "web")
  expect_message(
    build <- cellucid:::.establish_web_build(cache, url),
    "Downloading the Cellucid web app \\(2 files"
  )
  expect_identical(build$inventory$build_id, "fixture-build-1")
  expect_identical(cellucid:::.read_web_asset(build, "index.html")$content_type, "text/html; charset=utf-8")
  # The same published build: nothing is downloaded again.
  expect_silent(cellucid:::.establish_web_build(cache, url))

  # A new build replaces the cache whole.
  index <- sub("fixture-build-1", "fixture-build-2", readLines(file.path(source, "index.html")), fixed = TRUE)
  writeLines(index, file.path(source, "index.html"))
  bytes <- readBin(file.path(source, "index.html"), "raw", 1e5)
  inventory <- jsonlite::fromJSON(file.path(source, "cellucid-web-assets.json"), simplifyVector = FALSE)
  inventory$build_id <- "fixture-build-2"
  inventory$assets[[2]]$sha256 <- cellucid:::.native_sha256(bytes)
  inventory$assets[[2]]$bytes <- length(bytes)
  writeLines(jsonlite::toJSON(inventory, auto_unbox = TRUE), file.path(source, "cellucid-web-assets.json"))
  build <- suppressMessages(cellucid:::.establish_web_build(cache, url))
  expect_identical(build$inventory$build_id, "fixture-build-2")
  expect_identical(cellucid:::.verify_web_cache(cache)$build_id, "fixture-build-2")

  # A file the source serves wrong is refused, and the cache is kept.
  writeLines("tampered", file.path(source, "assets", "app.js"))
  inventory$build_id <- "fixture-build-3"
  writeLines(jsonlite::toJSON(inventory, auto_unbox = TRUE), file.path(source, "cellucid-web-assets.json"))
  expect_error(suppressMessages(cellucid:::.establish_web_build(cache, url)), "app.js' byte length is")
  expect_identical(cellucid:::.verify_web_cache(cache)$build_id, "fixture-build-2")
  expect_identical(list.files(dirname(cache), all.files = TRUE, no.. = TRUE), "web")

  # Offline: the verified cache is served; with no cache, the error says what to do.
  writeLines("", file.path(source, ".stop"))
  Sys.sleep(0.5)
  expect_message(build <- cellucid:::.establish_web_build(cache, url), "serving the cached web app")
  expect_identical(build$inventory$build_id, "fixture-build-2")
  expect_error(
    cellucid:::.establish_web_build(file.path(tempfile(), "web"), url),
    "could not be downloaded.*Connect to the internet once"
  )
})

test_that("clearing the cache removes its folder, and refuses the home folder", {
  cache <- viewer_fixture_web_cache()
  resolved <- normalizePath(cache, winslash = "/")
  expect_identical(cellucid_clear_web_cache(cache), resolved)
  expect_false(dir.exists(cache))
  expect_error(cellucid_clear_web_cache("~"), "dedicated web cache folder")
})
