# The Cellucid web app the viewer serves: one published build of
# www.cellucid.com, downloaded at first use into a local cache and checked
# file by file, with the rules of cellucid-python's web cache (web_cache.py),
# so both packages accept exactly the same builds:
#
# - The build publishes an inventory, cellucid-web-assets.json: exactly
#   {version: 1, build_id, assets}, every asset exactly {path, sha256, bytes,
#   content_type}, paths root-relative, portable, sorted and unique, and
#   index.html among them.
# - index.html carries the build's id in exactly one
#   <meta name="cellucid-web-build-id" content="..."> element, equal to the
#   inventory's build_id.
# - Every downloaded file must have the declared length and SHA-256, and a
#   response the declared media type; a cache holds exactly the declared
#   files, the inventory and the source record (.cellucid-web-source.json),
#   and no symbolic link.
# - A new build is staged in a folder beside the cache (the files the cache
#   already holds with the declared bytes copied, every other one downloaded,
#   six at once), checked whole, and only then swapped in, so a failed or
#   interrupted download never leaves a half-written cache.
#
# The cache is the published build of `source_url`: when the published build
# id equals the cached one and the cache checks out, nothing is downloaded;
# when the site cannot be reached, a cache that checks out is used as it is.
# The build is about 40 MB, far over what a CRAN package may ship, so it is
# never bundled. Its default home is tools::R_user_dir("cellucid", "cache").

.cellucid_web_url <- "https://www.cellucid.com"

# What a start says when the source cannot be reached, word for word as
# cellucid-python says it (web_cache.py WEB_APP_*; its
# tests/test_r_python_parity_contract.py pins the two): the cached build served
# offline, or why nothing can be served and how to get the app.
.web_app_offline_notice <- "Could not reach {source}; serving the cached web app (build {build_id})."
.web_app_unavailable <- paste0(
  "The Cellucid web app could not be downloaded from {source} ({reason}), and {cache}. ",
  "Connect to the internet once, or set web_cache_dir to a folder that holds a verified build."
)
.web_app_not_cached <- "no copy of it is cached in {cache_dir}"
.web_app_cache_unverified <- "the copy cached in {cache_dir} does not check out ({problem})"

# A message template with each {name} replaced by its value, as Python's str.format.
.fill_message <- function(template, ...) {
  values <- list(...)
  for (name in names(values)) {
    template <- gsub(paste0("{", name, "}"), as.character(values[[name]]), template, fixed = TRUE)
  }
  template
}
.web_inventory_file <- "cellucid-web-assets.json"
.web_source_file <- ".cellucid-web-source.json"
.web_build_id_max <- 200L
.windows_reserved_names <- c("CON", "PRN", "AUX", "NUL", paste0("COM", 1:9), paste0("LPT", 1:9))

.default_web_cache_dir <- function() {
  file.path(tools::R_user_dir("cellucid", "cache"), "web")
}

.require_source_url <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value) ||
      grepl("[[:space:]]", value)) {
    .abort("web_source_url must be one URL without spaces.")
  }
  if (endsWith(value, "/")) .abort("web_source_url must not end with '/'.")
  parts <- regmatches(value, regexec("^(https?)://([^/?#@]+)(/[^?#]*)?$", value))[[1L]]
  if (length(parts) == 0L) {
    .abort(
      "web_source_url must be an absolute http:// or https:// address with no ",
      "credentials, query or fragment."
    )
  }
  host <- sub(":[0-9]+$", "", parts[[3L]])
  host <- gsub("^\\[|\\]$", "", host)
  # Over plain http an on-path attacker could rewrite the assets and the
  # inventory's digests together; only a loopback host, where there is no
  # path to stand on, may use it (a locally built web app).
  if (identical(parts[[2L]], "http") && !.is_loopback_host(host)) {
    .abort(
      "web_source_url must use https; plain http is accepted only for a loopback ",
      "host (localhost, 127.0.0.0/8, ::1), not '", host, "'."
    )
  }
  value
}

.is_loopback_host <- function(host) {
  host <- tolower(host)
  host %in% c("localhost", "::1", "0:0:0:0:0:0:0:1") ||
    grepl("^127\\.[0-9]{1,3}\\.[0-9]{1,3}\\.[0-9]{1,3}$", host) ||
    grepl("^::ffff:127\\.", host)
}

.require_asset_path <- function(value) {
  bad <- function(why) .abort("Web asset path ", why, ": '", value, "'.")
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    .abort("Web asset path must be a non-empty string.")
  }
  if (startsWith(value, "/") || grepl("[\\\\?#]", value)) bad("is not root-relative POSIX")
  if (grepl("[[:space:]]", value)) bad("contains whitespace")
  if (grepl("[\001-\037\177]", value)) bad("contains a control character")
  if (grepl("[<>:\"|*]", value)) bad("is not portable across platforms")
  parts <- strsplit(value, "/", fixed = TRUE)[[1L]]
  if (endsWith(value, "/") || any(parts %in% c("", ".", ".."))) bad("contains an empty or dot segment")
  stems <- toupper(vapply(strsplit(parts, ".", fixed = TRUE), `[`, character(1L), 1L))
  if (any(endsWith(parts, ".")) || any(stems %in% .windows_reserved_names)) {
    bad("is not portable across platforms")
  }
  if (identical(value, .web_inventory_file)) {
    .abort("The web asset inventory cannot declare itself as an asset.")
  }
  value
}

.media_type_pattern <- paste0(
  "^[a-z0-9][a-z0-9!#$&^_.+-]*/[a-z0-9][a-z0-9!#$&^_.+-]*",
  "(; [a-z0-9][a-z0-9!#$&^_.+-]*=[a-z0-9][a-z0-9!#$&^_.+-]*)*$"
)

.require_content_type <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl(.media_type_pattern, value)) {
    .abort(label, " must be an exact normalized MIME type.")
  }
  value
}

# The names one kind of file goes by on hosts: JavaScript is text/javascript
# (RFC 9239), and its historic names are still what many hosts send (GitHub
# Pages and nginx application/javascript; Node's and Python's own servers
# text/javascript). The response check only tells a file of another kind (an
# HTML error page) apart.
.media_type_aliases <- c(
  "application/javascript" = "text/javascript",
  "application/x-javascript" = "text/javascript",
  "application/ecmascript" = "text/javascript",
  "text/ecmascript" = "text/javascript"
)

.media_type <- function(content_type) {
  bare <- tolower(trimws(strsplit(content_type, ";", fixed = TRUE)[[1L]][[1L]]))
  if (bare %in% names(.media_type_aliases)) unname(.media_type_aliases[[bare]]) else bare
}

# The one accepted inventory schema.
.parse_web_inventory <- function(bytes) {
  text <- .utf8_text(bytes)
  if (is.na(text)) .abort("Web asset inventory must be valid UTF-8.")
  raw <- tryCatch(
    jsonlite::parse_json(text, simplifyVector = FALSE),
    error = function(error) .abort("Web asset inventory must be valid JSON.")
  )
  if (!is.list(raw) || is.null(names(raw)) || !setequal(names(raw), c("version", "build_id", "assets")) ||
      length(raw) != 3L) {
    .abort("Web asset inventory must contain exactly version, build_id, and assets.")
  }
  if (!is.integer(raw$version) || !identical(raw$version, 1L)) {
    .abort("Web asset inventory version must be integer 1.")
  }
  build_id <- raw$build_id
  if (!is.character(build_id) || length(build_id) != 1L || !nzchar(build_id) ||
      !identical(build_id, trimws(build_id)) || nchar(build_id) > .web_build_id_max) {
    .abort("Web asset inventory build_id must be a non-empty exact string.")
  }
  if (!is.list(raw$assets) || length(raw$assets) == 0L || !is.null(names(raw$assets))) {
    .abort("Web asset inventory assets must be a non-empty array.")
  }
  assets <- lapply(seq_along(raw$assets), function(index) {
    asset <- raw$assets[[index]]
    if (!is.list(asset) || is.null(names(asset)) || length(asset) != 4L ||
        !setequal(names(asset), c("path", "sha256", "bytes", "content_type"))) {
      .abort(
        "Web asset inventory entry ", index - 1L, " must contain exactly path, ",
        "sha256, bytes, and content_type."
      )
    }
    path <- .require_asset_path(asset$path)
    if (!is.character(asset$sha256) || length(asset$sha256) != 1L || !grepl("^[0-9a-f]{64}$", asset$sha256)) {
      .abort("Web asset '", path, "' has an invalid sha256.")
    }
    size <- asset$bytes
    if (!is.numeric(size) || length(size) != 1L || is.na(size) || size < 0 || size != floor(size)) {
      .abort("Web asset '", path, "' has an invalid byte length.")
    }
    list(
      path = path,
      sha256 = asset$sha256,
      bytes = as.numeric(size),
      content_type = .require_content_type(asset$content_type, paste0("Web asset '", path, "' content_type"))
    )
  })
  paths <- vapply(assets, function(asset) asset$path, character(1L))
  if (!identical(paths, sort(paths, method = "radix"))) {
    .abort("Web asset inventory paths must be sorted.")
  }
  if (anyDuplicated(paths)) .abort("Web asset inventory paths must be unique.")
  if (!"index.html" %in% paths) .abort("Web asset inventory must contain index.html.")
  names(assets) <- paths
  list(version = 1L, build_id = build_id, assets = assets)
}

# The single exact build stamp of index.html.
.web_build_id_of <- function(bytes) {
  text <- .utf8_text(bytes)
  if (is.na(text)) .abort("index.html must be valid UTF-8.")
  tags <- regmatches(text, gregexpr("<meta\\b[^>]*>", text, ignore.case = TRUE, perl = TRUE))[[1L]]
  attribute <- function(tag, name) {
    found <- regmatches(tag, regexec(
      paste0("\\b", name, "\\s*=\\s*(\"([^\"]*)\"|'([^']*)')"), tag,
      ignore.case = TRUE, perl = TRUE
    ))[[1L]]
    if (length(found) == 0L) NULL else paste0(found[[3L]], found[[4L]])
  }
  ids <- character()
  for (tag in tags) {
    if (!identical(attribute(tag, "name"), "cellucid-web-build-id")) next
    content <- attribute(tag, "content")
    if (is.null(content)) {
      .abort("cellucid-web-build-id meta content must have a value.")
    }
    ids <- c(ids, content)
  }
  if (length(ids) != 1L) {
    .abort("index.html must contain exactly one cellucid-web-build-id meta element.")
  }
  id <- ids[[1L]]
  if (!nzchar(id) || !identical(id, trimws(id)) || nchar(id) > .web_build_id_max) {
    .abort("index.html contains an invalid cellucid-web-build-id.")
  }
  id
}

.read_bytes <- function(path) readBin(path, "raw", n = file.size(path))

.verify_asset_bytes <- function(asset, bytes) {
  if (length(bytes) != asset$bytes) {
    .abort(
      "Web asset '", asset$path, "' byte length is ", length(bytes), ", expected ",
      format(asset$bytes, scientific = FALSE), "."
    )
  }
  digest <- .native_sha256(bytes)
  if (!identical(digest, asset$sha256)) {
    .abort("Web asset '", asset$path, "' sha256 is ", digest, ", expected ", asset$sha256, ".")
  }
  invisible(bytes)
}

.source_record <- function(source_url) {
  paste0("{\"version\":1,\"source_url\":", .json_text(source_url), "}")
}

.read_cache_source <- function(cache_dir) {
  path <- file.path(cache_dir, .web_source_file)
  if (!file.exists(path)) .abort("Web cache source control is missing: ", path)
  raw <- tryCatch(
    .read_json_file(path),
    error = function(error) .abort("Web cache source control is invalid: ", path)
  )
  if (!is.list(raw) || !setequal(names(raw), c("version", "source_url")) || length(raw) != 2L ||
      !identical(raw$version, 1L)) {
    .abort("Web cache source control has a noncanonical schema.")
  }
  .require_source_url(raw$source_url)
}

# The inventory of a cache folder, its index checked against it; every
# declared file is checked when `files` is TRUE.
.verify_web_cache <- function(cache_dir, expected_source_url = NULL, files = TRUE) {
  if (!dir.exists(cache_dir)) .abort("Web cache directory does not exist: ", cache_dir)
  source <- .read_cache_source(cache_dir)
  if (!is.null(expected_source_url) && !identical(source, expected_source_url)) {
    .abort("Web cache source is '", source, "', expected '", expected_source_url, "'.")
  }
  inventory_path <- file.path(cache_dir, .web_inventory_file)
  if (!file.exists(inventory_path)) .abort("Web asset inventory is missing: ", inventory_path)
  inventory <- .parse_web_inventory(.read_bytes(inventory_path))
  index_path <- file.path(cache_dir, "index.html")
  if (!file.exists(index_path)) .abort("Web cache index is missing: ", index_path)
  index_id <- .web_build_id_of(.read_bytes(index_path))
  if (!identical(index_id, inventory$build_id)) {
    .abort(
      "Web cache index build id is '", index_id, "', inventory declares '",
      inventory$build_id, "'."
    )
  }
  if (!files) return(inventory)
  present <- list.files(cache_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  links <- present[nzchar(Sys.readlink(file.path(cache_dir, present)))]
  if (length(links) > 0L) {
    .abort("Web cache contains a symbolic link: ", file.path(cache_dir, links[[1L]]))
  }
  expected <- c(.web_inventory_file, .web_source_file, names(inventory$assets))
  missing <- setdiff(expected, present)
  unexpected <- setdiff(present, expected)
  if (length(missing) > 0L || length(unexpected) > 0L) {
    .abort(
      "Web cache file set differs from its inventory; missing ",
      .format_value_list(utils::head(missing, 5L)), ", unexpected ",
      .format_value_list(utils::head(unexpected, 5L)), "."
    )
  }
  for (asset in inventory$assets) {
    .verify_asset_bytes(asset, .read_bytes(file.path(cache_dir, asset$path)))
  }
  inventory
}

# The handle of one request to the source. The bytes are asked for as they are
# stored (Accept-Encoding: identity), as Python's web cache asks: a host that
# compresses what it sends (GitHub Pages gzips text) declares the length of
# the compressed body, which curl would otherwise decode to more bytes.
.web_handle <- function(timeout) {
  curl::new_handle(
    useragent = "cellucid-r/web-cache", timeout = timeout, followlocation = FALSE,
    accept_encoding = "identity"
  )
}

.web_url <- function(source_url, path) paste0(source_url, "/", utils::URLencode(path))

# One response of the source, read whole: status 200, a media type, and
# exactly the Content-Length it declares.
.check_web_response <- function(url, response) {
  if (!identical(as.integer(response$status_code), 200L)) {
    .abort("Expected HTTP 200 for ", url, ", received ", response$status_code, ".")
  }
  headers <- curl::parse_headers_list(response$headers)
  content_type <- .require_content_type(headers[["content-type"]], paste0("HTTP Content-Type for ", url))
  length_header <- headers[["content-length"]]
  if (is.null(length_header) || !grepl("^[0-9]+$", length_header)) {
    .abort("HTTP Content-Length is missing or invalid for ", url, ".")
  }
  if (length(response$content) != as.numeric(length_header)) {
    .abort(
      "HTTP body length for ", url, " is ", length(response$content), ", expected ",
      length_header, "."
    )
  }
  list(data = response$content, content_type = content_type)
}

.fetch_web <- function(source_url, path, timeout) {
  url <- .web_url(source_url, path)
  .check_web_response(url, curl::curl_fetch_memory(url, handle = .web_handle(timeout)))
}

# Downloads in flight at once, as web_cache.py: the six connections per host a
# browser opens to a site, so the source sees no more than a browser loading
# the app asks of it. One file at a time took over six minutes for the 1,209
# files of a build from www.cellucid.com (2026-11-11).
.concurrent_downloads <- 6L

# Several files of the source at once, each read whole and checked as
# .fetch_web checks one; named by path, in the order asked. The first failure
# stops the build.
.fetch_web_all <- function(source_url, paths, timeout) {
  pool <- curl::new_pool(total_con = .concurrent_downloads, host_con = .concurrent_downloads)
  responses <- list()
  failures <- list()
  for (path in paths) {
    local({
      asked <- path
      curl::curl_fetch_multi(
        .web_url(source_url, asked),
        done = function(response) responses[[asked]] <<- response,
        fail = function(message) failures[[asked]] <<- message,
        pool = pool,
        handle = .web_handle(timeout)
      )
    })
  }
  curl::multi_run(pool = pool)
  if (length(failures) > 0L) {
    .abort("Could not download ", .web_url(source_url, names(failures)[[1L]]), ": ", failures[[1L]])
  }
  stats::setNames(lapply(paths, function(path) .check_web_response(.web_url(source_url, path), responses[[path]])), paths)
}

# The bytes the cache from `source_url` holds for `asset` when they are
# exactly the ones the new inventory declares, else NULL: an unchanged file is
# not downloaded again.
.active_copy <- function(cache_dir, source_url, asset) {
  tryCatch(
    {
      if (!identical(.read_cache_source(cache_dir), source_url)) return(NULL)
      bytes <- .read_bytes(file.path(cache_dir, asset$path))
      .verify_asset_bytes(asset, bytes)
      bytes
    },
    error = function(error) NULL,
    warning = function(warning) NULL
  )
}

# Stage one build beside the cache -- each file the cache holds with the
# declared bytes copied from it, every other one downloaded (several at once)
# and checked -- check it whole, and swap it in.
.download_web_build <- function(cache_dir, source_url, inventory, inventory_bytes, timeout, quiet) {
  parent <- dirname(cache_dir)
  .dir_create(parent)
  stage <- tempfile(paste0(".", basename(cache_dir), ".staging-"), tmpdir = parent)
  .dir_create(stage)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  place <- function(path, bytes) {
    target <- file.path(stage, path)
    .dir_create(dirname(target))
    writeBin(bytes, target)
  }
  missing <- character()
  for (asset in inventory$assets) {
    kept <- .active_copy(cache_dir, source_url, asset)
    if (is.null(kept)) missing <- c(missing, asset$path) else place(asset$path, kept)
  }
  if (!quiet && length(missing) > 0L) {
    total <- sum(vapply(inventory$assets[missing], function(asset) asset$bytes, numeric(1L)))
    .inform(
      "Downloading the Cellucid web app (", length(missing), " ", .plural(length(missing), "file"), ", ",
      format(round(total / 1e6, 1), nsmall = 1), " MB) from ", source_url, " into ", cache_dir, "."
    )
  }
  responses <- .fetch_web_all(source_url, missing, timeout)
  for (path in missing) {
    asset <- inventory$assets[[path]]
    response <- responses[[path]]
    if (!identical(.media_type(response$content_type), .media_type(asset$content_type))) {
      .abort(
        "Web asset '", asset$path, "' response media type is '", .media_type(response$content_type),
        "', expected '", .media_type(asset$content_type), "'."
      )
    }
    .verify_asset_bytes(asset, response$data)
    place(asset$path, response$data)
  }
  writeBin(inventory_bytes, file.path(stage, .web_inventory_file))
  writeBin(charToRaw(.source_record(source_url)), file.path(stage, .web_source_file))
  .verify_web_cache(stage, expected_source_url = source_url)
  backup <- NULL
  if (dir.exists(cache_dir)) {
    backup <- tempfile(paste0(".", basename(cache_dir), ".previous-"), tmpdir = parent)
    if (!file.rename(cache_dir, backup)) .abort("Could not move the previous web cache aside: ", cache_dir)
  }
  if (!file.rename(stage, cache_dir)) {
    if (!is.null(backup)) file.rename(backup, cache_dir)
    .abort("Could not move the downloaded web build to ", cache_dir, ".")
  }
  if (!is.null(backup)) unlink(backup, recursive = TRUE, force = TRUE)
  invisible(cache_dir)
}

# The verified web build to serve: its folder and inventory. `source_url`
# NULL serves the cache as it is, with no network access.
.establish_web_build <- function(cache_dir, source_url, timeout = 15, quiet = FALSE) {
  if (is.null(source_url)) {
    return(.web_build(cache_dir, .verify_web_cache(cache_dir)))
  }
  source_url <- .require_source_url(source_url)
  published <- tryCatch(
    .fetch_web(source_url, .web_inventory_file, timeout),
    error = function(error) error
  )
  if (inherits(published, "error")) {
    # The source cannot be read (offline, blocked, down): a cache whose every
    # file matches its inventory's SHA-256 is that build exactly, so it is
    # served; else the error says what is missing and how to get it.
    cached <- tryCatch(.verify_web_cache(cache_dir, expected_source_url = source_url), error = function(error) error)
    if (inherits(cached, "error")) {
      cache <- if (dir.exists(cache_dir)) {
        .fill_message(.web_app_cache_unverified, cache_dir = cache_dir, problem = conditionMessage(cached))
      } else {
        .fill_message(.web_app_not_cached, cache_dir = cache_dir)
      }
      .abort(.fill_message(
        .web_app_unavailable,
        source = source_url, reason = conditionMessage(published), cache = cache
      ))
    }
    # Said even when quiet, as cellucid-python says it: the app is the cached
    # build, not what the source publishes now.
    .inform(.fill_message(.web_app_offline_notice, source = source_url, build_id = cached$build_id))
    return(.web_build(cache_dir, cached))
  }
  if (!identical(published$content_type, "application/json; charset=utf-8")) {
    .abort("Web asset inventory response Content-Type must be 'application/json; charset=utf-8'.")
  }
  inventory <- .parse_web_inventory(published$data)
  cached <- tryCatch(.verify_web_cache(cache_dir, expected_source_url = source_url), error = function(error) NULL)
  if (!is.null(cached) && identical(cached$build_id, inventory$build_id)) {
    return(.web_build(cache_dir, cached))
  }
  .download_web_build(cache_dir, source_url, inventory, published$data, timeout, quiet)
  .web_build(cache_dir, .parse_web_inventory(published$data))
}

# One declared asset, read and checked; NULL when the build does not declare it.
# The established build: its folder, its inventory, and the assets already
# checked against their digests (by size and modification time).
.web_build <- function(dir, inventory) {
  list(dir = dir, inventory = inventory, verified = new.env(parent = emptyenv()))
}

# One declared asset, checked against its digest the first time it is served
# (and again whenever its file changes), then served from disk: its bytes
# are neither read nor hashed again for every request. A file changes when
# its version does (.file_version(): size, modification and status-change
# time, identity), so one rewritten in place with its modification time put
# back, or replaced by a copy (cp -p), is checked again.
.read_web_asset <- function(build, path) {
  asset <- build$inventory$assets[[path]]
  if (is.null(asset)) return(NULL)
  file <- file.path(build$dir, path)
  version <- .file_version(file)
  if (is.null(version)) .abort("The web build lost ", path, "; clear it with cellucid_clear_web_cache().")
  stamp <- version$etag
  verified <- build$verified
  if (!is.environment(verified) || !identical(get0(.symbol_keys(path), envir = verified, inherits = FALSE), stamp)) {
    .verify_asset_bytes(asset, .read_bytes(file))
    if (is.environment(verified)) assign(.symbol_keys(path), stamp, envir = verified)
  }
  list(file = file, content_type = asset$content_type)
}

# Hand-written help: man/cellucid_clear_web_cache.Rd.
cellucid_clear_web_cache <- function(cache_dir = NULL) {
  cache_dir <- if (is.null(cache_dir)) .default_web_cache_dir() else .validate_required_string(cache_dir, "cache_dir")
  cache_dir <- normalizePath(path.expand(cache_dir), winslash = "/", mustWork = FALSE)
  protected <- normalizePath(c(path.expand("~"), dirname(path.expand("~")), getwd(), "/"), winslash = "/", mustWork = FALSE)
  if (cache_dir %in% protected || identical(dirname(cache_dir), cache_dir)) {
    .abort(
      "cache_dir must name a dedicated web cache folder, not ", cache_dir, ": clearing ",
      "removes the whole folder."
    )
  }
  if (dir.exists(cache_dir)) {
    if (nzchar(Sys.readlink(cache_dir))) .abort("The web cache must not be a symbolic link: ", cache_dir)
    unlink(cache_dir, recursive = TRUE, force = TRUE)
  }
  invisible(cache_dir)
}
