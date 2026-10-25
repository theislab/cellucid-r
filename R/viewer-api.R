# The viewer from R: cellucid_view() opens an object or a file in the RStudio
# or Positron Viewer pane or a browser, cellucid_serve() serves it for a
# browser elsewhere (a cluster through an SSH tunnel), and the functions
# below talk to an open viewer over its channel (viewer-protocol.R). The help
# pages are hand-written: man/cellucid_view.Rd and its siblings.

# The arguments cellucid_view() passes on to the dataset plan for an object.
.view_plan_arguments <- c(
  "embeddings", "dataset_id", "dataset_name", "obs", "var", "gene_expression",
  "connectivities", "vector_fields", "modalities", "gene_identifiers", "dataset_description",
  "citation_text", "citation_url", "source_name", "source_url", "embedding_defaults",
  "latent_key", "obs_keys", "var_gene_id_column", "vector_field_default", "modality",
  "normalize_embeddings", "centroid_min_points"
)

# For a path, the launch defaults the server states (the browser checks them
# against the file).
.view_path_arguments <- c(
  "dataset_id", "dataset_name", "embedding_defaults", "latent_key", "obs_keys",
  "var_gene_id_column", "vector_field_default", "modality", "normalize_embeddings",
  "centroid_min_points"
)

.default_port <- 8765L
.default_host <- "127.0.0.1"

cellucid_view <- function(
    data,
    ...,
    pane = c("auto", "viewer", "browser", "none"),
    port = 8765L,
    host = "127.0.0.1",
    allowed_hosts = NULL,
    serve_web_ui = TRUE,
    web_cache_dir = NULL,
    web_source_url = "https://www.cellucid.com",
    cache_bytes = 268435456,
    quiet = FALSE
) {
  if (missing(data)) {
    .abort("cellucid_view() needs data: a Seurat, SingleCellExperiment or SummarizedExperiment object, or a path.")
  }
  pane <- match.arg(pane)
  label <- paste(deparse(substitute(data), width.cutoff = 60L), collapse = " ")
  viewer <- .start_viewer(
    data, list(...), label = label, port = port, port_given = !missing(port), host = host,
    allowed_hosts = allowed_hosts, serve_web_ui = serve_web_ui, web_cache_dir = web_cache_dir, web_source_url = web_source_url,
    cache_bytes = cache_bytes, quiet = quiet
  )
  .open_pane(viewer, pane, quiet)
  invisible(viewer)
}

cellucid_serve <- function(
    data,
    ...,
    port = 8765L,
    host = "127.0.0.1",
    allowed_hosts = NULL,
    serve_web_ui = TRUE,
    web_cache_dir = NULL,
    web_source_url = "https://www.cellucid.com",
    cache_bytes = 268435456,
    quiet = FALSE
) {
  if (missing(data)) {
    .abort("cellucid_serve() needs data: a Seurat, SingleCellExperiment or SummarizedExperiment object, or a path.")
  }
  label <- paste(deparse(substitute(data), width.cutoff = 60L), collapse = " ")
  viewer <- .start_viewer(
    data, list(...), label = label, port = port, port_given = !missing(port), host = host,
    allowed_hosts = allowed_hosts, serve_web_ui = serve_web_ui, web_cache_dir = web_cache_dir, web_source_url = web_source_url,
    cache_bytes = cache_bytes, quiet = quiet
  )
  if (!quiet) {
    .inform(
      "Serving ", viewer$dataset_name, " at ", viewer$url, "\n",
      "Open ", viewer$viewer_url, " in a browser",
      if (!.is_loopback_host(host)) "" else paste0(
        " (from another machine: ssh -L ", viewer$port, ":127.0.0.1:", viewer$port,
        " <this host>, then open the address there)"
      ),
      ".\nStop it with cellucid_stop()."
    )
  }
  invisible(viewer)
}

# Hex text of `n_bytes` bytes from the system's cryptographic generator: the
# viewer id and token, and request ids. Never R's generator, which set.seed()
# makes predictable and whose use would move the caller's .Random.seed.
.random_hex <- function(n_bytes) {
  paste(as.character(.native_random_bytes(n_bytes)), collapse = "")
}

# A dataset id made from a label: portable, and named so the caller sees where
# it came from.
.label_id <- function(label) {
  id <- gsub("[^A-Za-z0-9._-]+", "-", label)
  id <- gsub("^[^A-Za-z0-9]+|[.-]+$", "", id)
  id <- substr(id, 1L, 180L)
  if (!nzchar(id) || toupper(id) %in% .windows_reserved_names) "r-object" else id
}

.start_viewer <- function(data, extra, label, port, port_given, host, allowed_hosts,
                          serve_web_ui, web_cache_dir, web_source_url, cache_bytes, quiet) {
  if (!requireNamespace("httpuv", quietly = TRUE)) {
    .abort("The viewer needs the httpuv package: install.packages(\"httpuv\").")
  }
  port <- .validate_port(port)
  host <- .validate_required_string(host, "host")
  allowed_hosts <- .validate_allowed_hosts(allowed_hosts)
  quiet <- .validate_flag(quiet, "quiet")
  if (!is.numeric(cache_bytes) || length(cache_bytes) != 1L || is.na(cache_bytes) || cache_bytes < 0) {
    .abort("cache_bytes must be one non-negative number of bytes.")
  }
  if (length(extra) > 0L && (is.null(names(extra)) || any(!nzchar(names(extra))))) {
    .abort("Arguments after data must be named, for example latent_key = \"X_pca\".")
  }

  viewer <- new.env(parent = emptyenv())
  class(viewer) <- "cellucid_viewer"
  viewer$host <- host
  viewer$allowed_hosts <- allowed_hosts
  viewer$viewer_id <- .random_hex(8L)
  viewer$viewer_token <- .random_hex(16L)
  viewer$etag <- paste0("\"", .random_hex(16L), "\"")
  viewer$queue <- list()
  viewer$hooks <- list()
  viewer$last <- list()
  viewer$session_requests <- list()
  viewer$running <- FALSE

  if (is.character(data) && length(data) == 1L && !is.na(data)) {
    allowed <- .view_path_arguments
    unknown <- setdiff(names(extra), allowed)
    if (length(unknown) > 0L) {
      .abort(
        "For a path, cellucid_view() takes ", .format_value_list(allowed), "; not ",
        .format_value_list(unknown), ". Pass an R object to choose its contents."
      )
    }
    .quietly(quiet, .serve_path(viewer, data, extra, quiet))
  } else {
    unknown <- setdiff(names(extra), .view_plan_arguments)
    if (length(unknown) > 0L) {
      .abort(
        "cellucid_view() has no argument ", .format_value_list(unknown), ". It takes the ",
        "arguments of cellucid_prepare() that choose what is shown: ",
        .format_value_list(.view_plan_arguments), "."
      )
    }
    .quietly(quiet, .serve_object(viewer, data, extra, label, quiet, cache_bytes))
  }

  viewer$web <- NULL
  if (.validate_flag(serve_web_ui, "serve_web_ui")) {
    if (is.null(web_cache_dir)) web_cache_dir <- .default_web_cache_dir()
    web_cache_dir <- normalizePath(path.expand(web_cache_dir), winslash = "/", mustWork = FALSE)
    viewer$web <- .establish_web_build(web_cache_dir, web_source_url, quiet = quiet)
  }

  .bind_viewer(viewer, port, port_given, quiet)
  viewer$health <- list(
    status = "ok", type = "anndata", version = as.character(utils::packageVersion("cellucid")),
    format = viewer$served$kind, n_cells = viewer$n_cells, n_genes = viewer$n_genes
  )
  viewer$info <- list(
    version = viewer$health$version, type = "anndata", format = viewer$served$kind,
    host = host, port = viewer$port, n_cells = viewer$n_cells, n_genes = viewer$n_genes,
    # The program serving the viewer, which the viewer names ("R").
    runtime = "r"
  )
  file <- list(url = paste0("_cellucid/anndata/", utils::URLencode(viewer$served$route_name, reserved = TRUE),
                            if (identical(viewer$served$layout, "folder")) "/" else ""),
               kind = viewer$served$kind, layout = viewer$served$layout)
  if (!identical(viewer$served$layout, "folder")) file$size <- viewer$served$size
  viewer$datasets <- list(datasets = list(list(
    id = viewer$dataset_id, name = viewer$dataset_name,
    anndata = list(file = file, defaults = .as_json_object(viewer$defaults))
  )))
  viewer
}

# `expr`, with its notes (columns not stored, gene names made unique, a
# matrix read block by block) silenced when quiet = TRUE.
.quietly <- function(quiet, expr) {
  if (isTRUE(quiet)) suppressMessages(expr) else expr
}

.validate_port <- function(port) {
  if (!is.numeric(port) || length(port) != 1L || is.na(port) || port != floor(port) || port < 0 || port > 65535) {
    .abort("port must be one whole number from 0 to 65535 (0 picks a free port).")
  }
  as.integer(port)
}

.validate_allowed_hosts <- function(hosts) {
  if (is.null(hosts)) return(character())
  .validate_character_vector(hosts, "allowed_hosts")
  hosts <- tolower(hosts)
  if (any(grepl("[:/*[:space:]]", hosts) & !grepl("^[0-9a-f:]+$", hosts))) {
    .abort(
      "allowed_hosts holds bare host names or IP addresses, with no scheme, port or wildcard, ",
      "for example \"rstudio.example.org\"."
    )
  }
  hosts
}

# An R object, planned once and held as a store in memory, uncompressed.
.serve_object <- function(viewer, data, extra, label, quiet, cache_bytes) {
  components <- .dataset_components(data, extra)
  if (is.null(components$dataset_name)) components$dataset_name <- label
  if (is.null(components$dataset_id)) {
    components$dataset_id <- .label_id(components$dataset_name)
    if (!quiet) {
      .inform(
        "Dataset id '", components$dataset_id, "' (sessions are kept under it); name it ",
        "with dataset_id = , or in misc$cellucid / metadata(x)$cellucid."
      )
    }
  }
  started <- proc.time()[["elapsed"]]
  plan <- .plan_dataset(components)
  store <- .store_open_memory(2L, cache_bytes)
  plan$write(store)
  .store_close(store)
  viewer$plan_seconds <- proc.time()[["elapsed"]] - started
  viewer$store <- store
  viewer$cell_names <- plan$obs_index
  viewer$dataset_id <- plan$dataset_id
  viewer$dataset_name <- plan$dataset_name
  viewer$n_cells <- plan$n_cells
  viewer$n_genes <- plan$n_genes
  viewer$modalities <- plan$modalities
  viewer$defaults <- list()
  viewer$served <- list(kind = plan$kind, layout = "folder", route_name = paste0(plan$dataset_id, ".zarr"), zarr_format = 2L)
  invisible(viewer)
}

# A file or a Zarr folder, served from disk as it is.
.serve_path <- function(viewer, path, extra, quiet) {
  classified <- .classify_path(path)
  info <- list(n_cells = NULL, n_genes = NULL)
  generated <- NULL
  if (identical(classified$layout, "folder")) {
    info <- .zarr_store_counts(classified$path, classified$zarr_format, classified$kind)
    generated <- .consolidated_zarr_metadata(classified$path, classified$zarr_format)
  }
  declared <- list(id = NULL, name = NULL)
  # A Zarr store's own identity (uns['cellucid']) names it unless the caller
  # did. An .h5ad or .h5mu file's is not read here (this package has no HDF5
  # reader): it is named after the file unless dataset_id/dataset_name are
  # given, where cellucid-python reads it from the file.
  if (classified$kind %in% c("zarr", "mudata-zarr") && (is.null(extra$dataset_id) || is.null(extra$dataset_name))) {
    declared <- .zarr_store_identity(classified$path, classified$layout, classified$zarr_format)
  }
  dataset_name <- extra$dataset_name %||% declared$name %||% basename(classified$path)
  dataset_id <- extra$dataset_id %||% declared$id %||%
    .label_id(sub("(\\.zarr\\.zip|\\.[A-Za-z0-9]+)$", "", basename(classified$path)))
  viewer$dataset_id <- .validate_dataset_id(dataset_id)
  viewer$dataset_name <- .validate_dataset_name(dataset_name)
  viewer$n_cells <- info$n_cells
  viewer$n_genes <- info$n_genes
  viewer$store <- NULL
  viewer$cell_names <- NULL
  viewer$defaults <- .served_defaults(extra)
  viewer$served <- list(
    kind = classified$kind, layout = classified$layout, path = classified$path,
    resolved = classified$path,
    route_name = paste0(viewer$dataset_id, if (identical(classified$layout, "folder")) ".zarr" else classified$ending),
    size = if (identical(classified$layout, "folder")) NULL else file.size(classified$path),
    generated = generated,
    zarr_format = classified$zarr_format
  )
  if (!quiet && is.null(extra$dataset_id) && is.null(declared$id)) {
    .inform("Dataset id '", viewer$dataset_id, "'; name it with dataset_id = to keep sessions across files.")
  }
  invisible(viewer)
}

# The launch defaults a server states for a file: the arguments given, in
# the served form (camelCase, 0-based columns).
.served_defaults <- function(extra) {
  out <- list()
  slots <- extra$embedding_defaults
  for (slot in names(slots)) {
    spec <- slots[[slot]]
    out[[paste0("embedding", slot)]] <- if (is.list(spec)) {
      list(key = spec$key, columns = I(as.integer(spec$columns) - 1L))
    } else {
      .validate_string(spec, "embedding_defaults")
    }
  }
  if (!is.null(extra$modality)) out$modality <- .validate_string(extra$modality, "modality")
  if (!is.null(extra$var_gene_id_column)) out$geneIdColumn <- .validate_string(extra$var_gene_id_column, "var_gene_id_column")
  if (!is.null(extra$latent_key)) out$latentKey <- .validate_string(extra$latent_key, "latent_key")
  if (!is.null(extra$vector_field_default)) out$vectorFieldDefault <- .validate_string(extra$vector_field_default, "vector_field_default")
  if (!is.null(extra$obs_keys)) {
    .validate_character_vector(extra$obs_keys, "obs_keys")
    out$obsKeys <- I(extra$obs_keys)
  }
  if (!is.null(extra$normalize_embeddings)) out$normalizeEmbeddings <- .validate_flag(extra$normalize_embeddings, "normalize_embeddings")
  if (!is.null(extra$centroid_min_points)) out$centroidMinPoints <- .validate_positive_integer(extra$centroid_min_points, "centroid_min_points")
  out
}

# Bind the server: `port`, or when the default 8765 is taken the next free
# port (where Google sign-in, which accepts only 8765, is unavailable). A port
# the caller chose is used or refused.
.bind_viewer <- function(viewer, port, port_given, quiet) {
  app <- list(
    call = function(req) .on_call(viewer, req),
    onHeaders = function(req) .on_headers(viewer, req),
    onWSOpen = function(ws) .on_ws_open(viewer, ws)
  )
  candidates <- if (port_given || port == 0L) port else seq.int(port, 65535L)
  server <- NULL
  for (candidate in candidates) {
    server <- tryCatch(
      httpuv::startServer(viewer$host, if (candidate == 0L) httpuv::randomPort(host = viewer$host) else candidate, app, quiet = TRUE),
      error = function(error) error
    )
    if (!inherits(server, "error")) break
  }
  if (inherits(server, "error")) {
    .abort(
      "Could not serve on ", viewer$host, ":", port, " (", conditionMessage(server), "). ",
      if (port_given) "Choose another port, or stop the viewer using it with cellucid_stop_all()." else ""
    )
  }
  viewer$server <- server
  viewer$port <- server$getPort()
  viewer$running <- TRUE
  shown_host <- if (.host_is_wildcard(viewer$host)) "127.0.0.1" else viewer$host
  authority <- if (grepl(":", shown_host, fixed = TRUE)) paste0("[", shown_host, "]") else shown_host
  viewer$url <- paste0("http://", authority, ":", viewer$port)
  viewer$viewer_url <- paste0(
    viewer$url, "/?viewer=true&viewerId=", viewer$viewer_id,
    "&viewerToken=", viewer$viewer_token, "&anndata=true"
  )
  assign(viewer$viewer_id, viewer, envir = .viewer_registry)
  if (!quiet && !port_given && viewer$port != .default_port) {
    .inform(
      "Port ", .default_port, " is in use; serving on ", viewer$port, ". Google sign-in ",
      "(community annotation) works only on port 8765."
    )
  }
  invisible(viewer)
}

.host_is_wildcard <- function(host) host %in% c("0.0.0.0", "::", "")

# Show the viewer: the IDE's Viewer pane (RStudio, Positron), else the
# browser; none for cellucid_serve() and pane = "none".
.open_pane <- function(viewer, pane, quiet) {
  if (identical(pane, "none")) return(invisible(NULL))
  viewer_option <- getOption("viewer")
  in_ide <- is.function(viewer_option)
  if (identical(pane, "viewer") && !in_ide) {
    .abort("pane = \"viewer\" needs RStudio or Positron; use pane = \"browser\".")
  }
  if (pane %in% c("auto", "viewer") && in_ide) {
    viewer_option(viewer$viewer_url)
    if (!quiet) {
      .inform(
        "Opened ", viewer$dataset_name, " in the Viewer pane. If it stays blank (no WebGL ",
        "in the pane), view it with pane = \"browser\"."
      )
    }
    return(invisible(NULL))
  }
  if (!interactive()) {
    if (!quiet) .inform("Open ", viewer$viewer_url, " in a browser.")
    return(invisible(NULL))
  }
  utils::browseURL(viewer$viewer_url)
  invisible(NULL)
}

print.cellucid_viewer <- function(x, ...) {
  state <- if (isTRUE(x$running)) "running" else "stopped"
  cat(
    "<Cellucid viewer> ", x$dataset_name, " (", x$dataset_id, "), ", state, "\n",
    "  address:  ", x$viewer_url %||% "", "\n",
    "  data:     ", x$served$kind, ", ", format(x$n_cells %||% NA), " ", .plural(x$n_cells %||% NA, "cell"),
    if (!is.null(x$modalities)) paste0(", modalities ", paste(x$modalities, collapse = ", ")), "\n",
    "  channel:  ", if (is.null(x$socket)) "waiting for the viewer" else "connected", "\n",
    sep = ""
  )
  invisible(x)
}

.require_viewer <- function(viewer) {
  if (!inherits(viewer, "cellucid_viewer")) {
    .abort("Pass the viewer cellucid_view() or cellucid_serve() returned.")
  }
  viewer
}

cellucid_stop <- function(viewer) {
  .require_viewer(viewer)
  if (isTRUE(viewer$running)) {
    viewer$running <- FALSE
    if (!is.null(viewer$socket)) {
      try(viewer$socket$close(.viewer_close_codes[["goingAway"]], "The R session stopped the viewer"), silent = TRUE)
      viewer$socket <- NULL
    }
    failure <- tryCatch({
      viewer$server$stop()
      NULL
    }, error = conditionMessage)
    for (request in names(viewer$session_requests)) viewer$session_requests[[request]] <- NULL
    viewer$store <- NULL
    if (!is.null(failure)) {
      # The server may still hold its port: the viewer stays listed, so
      # cellucid_stop_all() can try again.
      .abort(
        "The viewer's server on port ", viewer$port, " did not stop (", failure, "). ",
        "Try cellucid_stop_all() again, or restart R to free the port."
      )
    }
  }
  if (exists(viewer$viewer_id, envir = .viewer_registry, inherits = FALSE)) {
    rm(list = viewer$viewer_id, envir = .viewer_registry)
  }
  invisible(viewer)
}

cellucid_stop_all <- function() {
  viewers <- mget(ls(.viewer_registry), envir = .viewer_registry)
  for (viewer in viewers) cellucid_stop(viewer)
  invisible(length(viewers))
}

cellucid_viewers <- function() {
  unname(mget(ls(.viewer_registry), envir = .viewer_registry))
}

# --- Hooks and events -----------------------------------------------------------

.hook_events <- c(
  "selection", "hover", "click", "ready", "data_mapping", "command_error", "pong",
  "debug_snapshot", "session_bundle"
)

cellucid_on <- function(viewer, event, handler) {
  .require_viewer(viewer)
  event <- match.arg(event, .hook_events)
  if (!is.function(handler)) .abort("handler must be a function of one argument, the event.")
  viewer$hooks[[event]] <- c(viewer$hooks[[event]], handler)
  invisible(viewer)
}

cellucid_off <- function(viewer, event = NULL) {
  .require_viewer(viewer)
  if (is.null(event)) {
    viewer$hooks <- list()
  } else {
    viewer$hooks[[match.arg(event, .hook_events)]] <- NULL
  }
  invisible(viewer)
}

cellucid_selection <- function(viewer) {
  .require_viewer(viewer)
  viewer$selection
}

# Run R's event loop until `done()` holds or `timeout` seconds pass: the one
# place this package waits on the viewer, and only when the caller asked to.
.wait_for <- function(viewer, done, timeout, what) {
  deadline <- proc.time()[["elapsed"]] + timeout
  while (!done()) {
    if (!isTRUE(viewer$running)) .abort("The viewer was stopped while waiting for ", what, ".")
    if (proc.time()[["elapsed"]] > deadline) {
      .abort(
        "The viewer did not answer within ", timeout, " seconds (waiting for ", what, "). ",
        "Is it open? ", if (is.null(viewer$socket)) paste0("No viewer is connected: open ", viewer$viewer_url, ".") else ""
      )
    }
    later::run_now(0.05)
  }
  invisible(TRUE)
}

cellucid_wait_ready <- function(viewer, timeout = 30) {
  .require_viewer(viewer)
  .wait_for(viewer, function() !is.null(viewer$ready), timeout, "the viewer to be ready")
  viewer$ready
}

# A report of what works between this R session and the viewer, as
# cellucid-python's AnnDataViewer.debug_connection() gives one: the server's
# own routes asked over HTTP (through R's event loop, which answers them),
# the store's first byte, and a ping and a debug snapshot over the channel.
# Hand-written help: man/cellucid_debug_connection.Rd.
cellucid_debug_connection <- function(viewer, timeout = 5) {
  .require_viewer(viewer)
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) || timeout <= 0) {
    .abort("timeout must be one positive number of seconds.")
  }
  report <- list(
    viewer_id = viewer$viewer_id,
    viewer_url = viewer$viewer_url,
    server_url = viewer$url,
    server_running = isTRUE(viewer$running),
    channel_connected = !is.null(viewer$socket),
    web_ui = list(cache_dir = viewer$web$dir, served = !is.null(viewer$web)),
    state = list(ready = viewer$ready, last_command_error = viewer$last$command_error)
  )
  if (!isTRUE(viewer$running)) {
    return(report)
  }
  json <- function(route) {
    response <- .debug_fetch(paste0(viewer$url, "/", route), timeout)
    .parse_json_bytes(response$content, route)
  }
  for (route in c("health", "info", "datasets")) {
    key <- paste0("server_", route)
    value <- tryCatch(json(paste0("_cellucid/", route)), error = function(error) structure(conditionMessage(error), class = "probe_error"))
    if (inherits(value, "probe_error")) report[[paste0(key, "_error")]] <- unclass(value) else report[[key]] <- value
  }
  # The store, probed as the viewer opens it: HEAD for its size and identity,
  # then its first byte with Range (a Zarr folder at its root metadata).
  file <- viewer$datasets$datasets[[1L]]$anndata$file
  target <- paste0(viewer$url, "/", file$url)
  if (endsWith(target, "/")) {
    target <- paste0(target, if (identical(viewer$served$zarr_format, 3L)) "zarr.json" else ".zmetadata")
  }
  report$anndata_store_probe <- tryCatch({
    head <- .debug_fetch(target, timeout, method = "HEAD")
    headers <- curl::parse_headers_list(head$headers)
    ranged <- .debug_fetch(target, timeout, headers = list(Range = "bytes=0-0"))
    list(
      probed = target,
      head = list(status = head$status_code, content_length = headers[["content-length"]],
                  accept_ranges = headers[["accept-ranges"]], etag = headers[["etag"]]),
      range = list(status = ranged$status_code,
                   content_range = curl::parse_headers_list(ranged$headers)[["content-range"]])
    )
  }, error = function(error) list(probed = target, error = conditionMessage(error)))
  if (is.null(viewer$socket)) {
    reason <- paste0(
      "The viewer's channel is not connected: open ", viewer$viewer_url,
      ", and check that any proxy in between passes WebSockets."
    )
    report$frontend_roundtrip <- list(ok = FALSE, error = reason)
    report$frontend_debug_snapshot <- list(ok = FALSE, error = reason)
    return(report)
  }
  ask <- function(type) {
    request <- paste0("debug-", .random_hex(6L))
    tryCatch({
      .viewer_send(viewer, list(type = type, requestId = request))
      answer <- if (identical(type, "ping")) "pong" else "debug_snapshot"
      .wait_for(viewer, function() identical(viewer$last[[answer]]$requestId, request), timeout, paste0("the viewer's ", answer))
      list(ok = TRUE, answer = viewer$last[[answer]])
    }, error = function(error) list(ok = FALSE, error = conditionMessage(error)))
  }
  report$frontend_roundtrip <- ask("ping")
  report$frontend_debug_snapshot <- ask("debug_snapshot")
  report
}

# One HTTP request to this R session's own server: through curl's multi
# interface and R's event loop, which is where the server answers it.
.debug_fetch <- function(url, timeout, method = "GET", headers = list()) {
  handle <- curl::new_handle(customrequest = method, nobody = identical(method, "HEAD"))
  if (length(headers) > 0L) curl::handle_setheaders(handle, .list = headers)
  result <- NULL
  failure <- NULL
  pool <- curl::new_pool()
  curl::curl_fetch_multi(url, done = function(response) result <<- response,
                         fail = function(message) failure <<- message, pool = pool, handle = handle)
  deadline <- proc.time()[["elapsed"]] + timeout
  while (is.null(result) && is.null(failure)) {
    if (proc.time()[["elapsed"]] > deadline) .abort("no answer from ", url, " within ", timeout, " seconds")
    curl::multi_run(timeout = 0, poll = TRUE, pool = pool)
    later::run_now(0.01)
  }
  if (!is.null(failure)) .abort(failure)
  result
}

# --- Commands -------------------------------------------------------------------

# Cells as the viewer's 0-based indices: 1-based positions, cell names, or a
# logical vector over the cells.
.cell_indices <- function(viewer, cells, what) {
  n <- viewer$n_cells
  if (is.logical(cells)) {
    if (!is.null(n) && length(cells) != n) .abort(what, " as TRUE/FALSE needs one value per cell (", n, ").")
    if (anyNA(cells)) {
      .abort(what, " as TRUE/FALSE holds NA at cell ", which(is.na(cells))[[1L]], "; say TRUE or FALSE for every cell.")
    }
    cells <- which(cells)
  } else if (is.character(cells)) {
    if (is.null(viewer$cell_names)) .abort(what, " by name needs a viewer of an R object; give positions.")
    positions <- match(cells, viewer$cell_names)
    if (anyNA(positions)) {
      .abort(what, " names cells the dataset does not have: ", .format_value_list(utils::head(cells[is.na(positions)], 5L)), ".")
    }
    cells <- positions
  } else if (!is.numeric(cells)) {
    .refuse_cells(what, n)
  }
  if (length(cells) > 0L) {
    # One pass, no copy (src/extent.c): no NA, every position a whole number
    # from 1 to the number of cells.
    extent <- .native_extent(cells, whole = is.double(cells))
    if (extent[[3L]] > 0 || extent[[1L]] < 1 || (!is.null(n) && extent[[2L]] > n) ||
        extent[[2L]] > 4294967296 || (is.double(cells) && extent[[5L]] != 1)) {
      .refuse_cells(what, n)
    }
  }
  # Ascending and each once, as the wire carries them: positions given so
  # (which(), a range) are used as they are, others sorted once; the frame is
  # then written from them in one pass (.binary_command_frame()).
  if (is.unsorted(cells, strictly = TRUE)) {
    cells <- sort(cells, method = "radix")
    if (length(cells) > 1L) cells <- cells[c(TRUE, cells[-1L] != cells[-length(cells)])]
  }
  cells - 1L
}

.refuse_cells <- function(what, n) {
  .abort(what, " must be cell names, 1-based positions", if (!is.null(n)) paste0(" from 1 to ", n), ", or TRUE/FALSE per cell.")
}

cellucid_highlight <- function(viewer, cells, color = "#e41a1c") {
  .require_viewer(viewer)
  indices <- .cell_indices(viewer, cells, "cells")
  if (length(indices) == 0L) .abort("cells selects no cell to highlight.")
  .viewer_send(viewer, list(type = "highlight", cells = indices, color = color))
}

cellucid_clear_highlights <- function(viewer) {
  .viewer_send(.require_viewer(viewer), list(type = "clearHighlights"))
}

cellucid_color_by <- function(viewer, field) {
  .viewer_send(.require_viewer(viewer), list(type = "setColorBy", field = field))
}

cellucid_set_visibility <- function(viewer, cells = NULL, visible = TRUE) {
  .require_viewer(viewer)
  indices <- if (is.null(cells)) NULL else .cell_indices(viewer, cells, "cells")
  .viewer_send(viewer, list(type = "setVisibility", cells = indices, visible = visible))
}

cellucid_reset_camera <- function(viewer) {
  .viewer_send(.require_viewer(viewer), list(type = "resetCamera"))
}

cellucid_freeze <- function(viewer) {
  .viewer_send(.require_viewer(viewer), list(type = "freeze"))
}

# The viewer's data mapping with the named fields replaced, sent and
# answered: the mapping the viewer applied (its data_mapping event).
cellucid_set_data_mapping <- function(viewer, ..., timeout = 30) {
  .require_viewer(viewer)
  fields <- list(...)
  if (length(fields) == 0L || is.null(names(fields)) || any(!nzchar(names(fields)))) {
    .abort("Name the fields to change, for example expression = \"layer:counts\".")
  }
  if (is.null(viewer$data_mapping)) {
    .wait_for(viewer, function() !is.null(viewer$data_mapping), timeout, "the viewer's data mapping")
  }
  mapping <- viewer$data_mapping$mapping
  for (name in names(fields)) mapping[name] <- list(fields[[name]])
  request <- paste0("r", .random_hex(6L))
  # An answer to an earlier call (one that timed out, say) is not this one's.
  viewer$last$data_mapping <- NULL
  viewer$last$command_error <- NULL
  .viewer_send(viewer, list(type = "setDataMapping", requestId = request, mapping = mapping))
  answered <- function() {
    identical(viewer$last$data_mapping$requestId, request) ||
      identical(viewer$last$command_error$command, "setDataMapping")
  }
  .wait_for(viewer, answered, timeout, "the data mapping to apply")
  if (!identical(viewer$last$data_mapping$requestId, request)) {
    error <- viewer$last$command_error
    viewer$last$command_error <- NULL
    .abort("The viewer refused the data mapping: ", error$reason)
  }
  invisible(viewer$data_mapping$mapping)
}

# The session the viewer holds now, as a .cellucid-session file.
cellucid_session <- function(viewer, file = NULL, timeout = 60) {
  .require_viewer(viewer)
  request <- paste0("s", .random_hex(6L))
  viewer$session_requests[[request]] <- list(asked = Sys.time())
  on.exit(viewer$session_requests[[request]] <- NULL, add = TRUE)
  viewer$last$session_bundle <- NULL
  viewer$last$command_error <- NULL
  .viewer_send(viewer, list(type = "requestSessionBundle", requestId = request))
  answered <- function() {
    identical(viewer$last$session_bundle$requestId, request) ||
      identical(viewer$last$command_error$command, "requestSessionBundle")
  }
  .wait_for(viewer, answered, timeout, "the session")
  if (!identical(viewer$last$session_bundle$requestId, request)) {
    .abort(viewer$last$command_error$reason, " Try cellucid_session() again.")
  }
  path <- viewer$last$session_bundle$path
  if (!is.null(file)) {
    file <- .validate_required_string(file, "file")
    if (!file.copy(path, file, overwrite = TRUE)) .abort("Could not write the session to ", file, ".")
    unlink(path)
    path <- normalizePath(file, winslash = "/")
  }
  path
}
