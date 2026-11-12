# The viewer from R: cellucid_view() opens an object or a file in the RStudio
# or Positron Viewer pane, a browser, or -- in a Jupyter notebook -- the
# cell's output (viewer-widget.R), cellucid_serve() serves it for a browser
# elsewhere (a cluster through an SSH tunnel), and the functions below talk to
# an open viewer over its channel (viewer-protocol.R). The help pages are
# hand-written: man/cellucid_view.Rd and its siblings.

# The arguments cellucid_view() passes on to the dataset plan for an object.
.view_plan_arguments <- c(
  "embeddings", "dataset_id", "dataset_name", "obs", "obs_colors", "var", "gene_expression",
  "connectivities", "vector_fields", "modalities", "gene_identifiers", "total_counts", "total_counts_over", "dataset_description",
  "citation_text", "citation_url", "source_name", "source_url", "embedding_defaults",
  "latent_key", "obs_keys", "var_gene_id_column", "vector_field_default", "modality",
  "normalize_embeddings", "centroid_min_points"
)

# For a path, the launch defaults the server states (the browser checks them
# against the file). knn_key and expression_source choose among what a file
# holds (its obsp graphs; X, raw.X or a layer), as cellucid-python's
# serve_anndata() and cellucid serve --knn/--expression do; an object's
# graph and matrix are chosen by passing them (connectivities,
# gene_expression).
.view_path_arguments <- c(
  "dataset_id", "dataset_name", "embedding_defaults", "latent_key", "obs_keys",
  "var_gene_id_column", "vector_field_default", "modality", "normalize_embeddings",
  "centroid_min_points", "knn_key", "expression_source"
)

.default_port <- 8765L
.default_host <- "127.0.0.1"

cellucid_view <- function(
    data,
    ...,
    pane = c("auto", "viewer", "browser", "notebook", "none"),
    port = 8765L,
    host = "127.0.0.1",
    allowed_hosts = NULL,
    serve_web_ui = TRUE,
    web_cache_dir = NULL,
    web_source_url = "https://www.cellucid.com",
    cache_bytes = NULL,
    quiet = FALSE
) {
  if (missing(data)) {
    .abort("cellucid_view() needs data: a Seurat, SingleCellExperiment or SummarizedExperiment object, or a path.")
  }
  pane <- match.arg(pane)
  label <- paste(deparse(substitute(data), width.cutoff = 60L), collapse = " ")
  # In a Jupyter notebook the viewer is the cell's output, reached through
  # the kernel's own connection: no server, no port (viewer-widget.R).
  if (identical(pane, "notebook") && !.in_jupyter_kernel()) {
    .abort("pane = \"notebook\" shows the viewer in a Jupyter notebook's output (an IRkernel session); use pane = \"auto\" here.")
  }
  notebook <- .in_jupyter_kernel()
  if (notebook && pane %in% c("viewer", "browser", "none")) {
    .abort("In a Jupyter notebook the viewer is the cell's output: use pane = \"auto\" or \"notebook\".")
  }
  if (notebook) web_source_url <- .validate_widget_source_url(web_source_url)
  viewer <- .start_viewer(
    data, list(...), label = label, port = port, port_given = !missing(port), host = host,
    allowed_hosts = allowed_hosts, serve_web_ui = serve_web_ui && !notebook, web_cache_dir = web_cache_dir,
    web_source_url = web_source_url, cache_bytes = cache_bytes, quiet = quiet, bind = !notebook
  )
  if (notebook) {
    tryCatch(.widget_show(viewer, web_source_url), error = function(error) {
      cellucid_stop(viewer)
      stop(error)
    })
  } else {
    .open_pane(viewer, pane, quiet)
  }
  .release_replaced_viewers(viewer, quiet)
  invisible(viewer)
}

# The notebook cell running now: its front end's id (JupyterLab, Notebook 7
# and VS Code send `cellId` with each run) and that run's execution count,
# from IRkernel's executor on the call stack; NULL outside IRkernel, or where
# the front end sends no cell id.
.notebook_run <- function() {
  if (!.in_jupyter_kernel()) return(NULL)
  for (frame in rev(sys.frames())) {
    executor <- tryCatch(get(".self", envir = frame), error = function(error) NULL)
    if (!inherits(executor, "Executor")) next
    cell <- executor$current_request$metadata$cellId
    if (!is.character(cell) || length(cell) != 1L || is.na(cell) || !nzchar(cell)) return(NULL)
    return(list(cell = cell, count = executor$execution_count))
  }
  NULL
}

# The viewers each notebook cell showed in its last run, by cell id: that
# run's execution count and the viewers' ids (looked up in the registry, so
# a stopped viewer is not kept alive here).
.shown_by_cell <- new.env(parent = emptyenv())

# `viewer` was just shown by cellucid_view(). A cell run again replaces its
# output: a viewer the cell showed in an earlier run, and no other cell
# shows, is stopped (its object's matrices and chunk cache released), as
# cellucid-python stops the widget a re-run cell replaced. Where no cell id
# is known (a terminal, a console, a front end that sends none), showing a
# dataset again stops the earlier viewer of the same dataset id, and says so.
.release_replaced_viewers <- function(viewer, quiet) {
  run <- .notebook_run()
  if (is.null(run)) {
    for (other in cellucid_viewers()) {
      if (identical(other$viewer_id, viewer$viewer_id) || !isTRUE(other$running)) next
      if (!identical(other$dataset_id, viewer$dataset_id)) next
      cellucid_stop(other)
      if (!isTRUE(quiet)) {
        .inform("Stopped the earlier viewer of '", viewer$dataset_id, "': cellucid_view() showed it again.")
      }
    }
    return(invisible(NULL))
  }
  previous <- if (exists(run$cell, envir = .shown_by_cell, inherits = FALSE)) get(run$cell, envir = .shown_by_cell) else NULL
  if (!is.null(previous) && identical(previous$count, run$count)) {
    previous$ids <- c(previous$ids, viewer$viewer_id)
    assign(run$cell, previous, envir = .shown_by_cell)
    return(invisible(NULL))
  }
  assign(run$cell, list(count = run$count, ids = viewer$viewer_id), envir = .shown_by_cell)
  if (is.null(previous)) return(invisible(NULL))
  elsewhere <- unlist(lapply(setdiff(ls(.shown_by_cell), run$cell), function(cell) get(cell, envir = .shown_by_cell)$ids))
  for (id in setdiff(previous$ids, c(viewer$viewer_id, elsewhere))) {
    if (exists(id, envir = .viewer_registry, inherits = FALSE)) cellucid_stop(get(id, envir = .viewer_registry))
  }
  invisible(NULL)
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
    cache_bytes = NULL,
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

# What one served store keeps of what it makes, by default: a 64th of the
# machine's memory (.machine_memory_bytes()): at most what the web app gives
# its decoded chunks (a 32nd of what its working set leaves of half the
# device, data/memory-budget.js; a 64th of the device when nothing is drawn).
# cellucid-python's served stores follow the same rule
# (_served_store.default_cache_bytes()).
.default_cache_bytes <- function() {
  floor(.machine_memory_bytes() / 64)
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
                          serve_web_ui, web_cache_dir, web_source_url, cache_bytes, quiet,
                          bind = TRUE, runtime = "r") {
  if (bind && .in_jupyter_kernel()) {
    # IRkernel waits for the notebook's messages without running R's event
    # loop, which runs an httpuv server's handlers: a server started here
    # would answer nothing until a cell runs.
    .abort(
      "In a Jupyter notebook the R kernel cannot run a server between cells, so a browser could ",
      "not read from it. Show the viewer in the cell's output with cellucid_view(data), or serve ",
      "it from R outside the notebook (R in a terminal)."
    )
  }
  if (!requireNamespace("httpuv", quietly = TRUE)) {
    .abort("The viewer needs the httpuv package: install.packages(\"httpuv\").")
  }
  port <- .validate_port(port)
  host <- .validate_required_string(host, "host")
  allowed_hosts <- .validate_allowed_hosts(allowed_hosts)
  quiet <- .validate_flag(quiet, "quiet")
  if (is.null(cache_bytes)) cache_bytes <- .default_cache_bytes()
  if (!is.numeric(cache_bytes) || length(cache_bytes) != 1L || is.na(cache_bytes) || cache_bytes < 0) {
    .abort("cache_bytes must be NULL or one non-negative number of bytes.")
  }
  if (length(extra) > 0L && (is.null(names(extra)) || any(!nzchar(names(extra))))) {
    .abort("Arguments after data must be named, for example latent_key = \"X_pca\".")
  }

  viewer <- new.env(parent = emptyenv())
  class(viewer) <- "cellucid_viewer"
  # What opened the viewer, as _cellucid/info names it and every sentence of
  # the viewer names it (the shared protocol's runtimes): "r" (an R session:
  # cellucid_view(), cellucid_serve()), "r-mcp" (cellucid_mcp(open = )).
  viewer$runtime <- runtime
  viewer$host <- host
  viewer$allowed_hosts <- allowed_hosts
  viewer$viewer_id <- .random_hex(8L)
  viewer$viewer_token <- .random_hex(16L)
  viewer$etag <- paste0("\"", .random_hex(16L), "\"")
  viewer$queue <- list()
  viewer$hooks <- list()
  # The calls waiting for an event (cellucid_wait_event()).
  viewer$event_waiters <- list()
  viewer$last <- list()
  viewer$session_requests <- list()
  viewer$action_requests <- list()
  viewer$session_loads <- list()
  viewer$turn <- NULL
  viewer$running <- FALSE
  # Its end of local agent connections (cellucid_mcp(); R/agent-link.R).
  .agent_port(viewer)

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

  if (bind) {
    .bind_viewer(viewer, port, port_given, quiet)
  } else {
    # A notebook widget's viewer: answered through the widget, on no port.
    viewer$port <- NULL
    viewer$running <- TRUE
    assign(viewer$viewer_id, viewer, envir = .viewer_registry)
  }
  viewer$health <- list(
    status = "ok", type = "anndata", version = as.character(utils::packageVersion("cellucid")),
    format = viewer$served$kind, n_cells = viewer$n_cells, n_genes = viewer$n_genes
  )
  viewer$info <- list(
    version = viewer$health$version, type = "anndata", format = viewer$served$kind,
    host = if (bind) host, port = viewer$port, n_cells = viewer$n_cells, n_genes = viewer$n_genes,
    # What opened the viewer, which the viewer names in what it says.
    runtime = viewer$runtime
  )
  # A viewer with a server is published to agent apps on this computer
  # (cellucid_mcp(), cellucid-python's cellucid mcp); a notebook widget's,
  # which has none, is not (R/agent-link.R).
  if (bind) .agent_publish(viewer)
  viewer
}

# The dataset list the viewer reads (_cellucid/datasets): the dataset, where
# its store is and the launch defaults. A file served as it is is described
# as it is when the list is asked for -- its size and the validator its bytes
# are answered with (.file_version()) -- so a file rewritten while it is
# served is a new version to the viewer, which files what it keeps of the
# file (its byte cache, a gene index, a converted store) under that version.
# A store folder states neither: its metadata's validator names its version.
# NULL when the file served is gone (the route answers 404).
.datasets_document <- function(viewer) {
  served <- viewer$served
  file <- list(url = paste0("_cellucid/anndata/", utils::URLencode(served$route_name, reserved = TRUE),
                            if (identical(served$layout, "folder")) "/" else ""),
               kind = served$kind, layout = served$layout)
  if (!identical(served$layout, "folder")) {
    version <- .file_version(served$path)
    if (is.null(version)) return(NULL)
    file$size <- version$size
    file$etag <- version$etag
  }
  list(datasets = list(list(
    id = viewer$dataset_id, name = viewer$dataset_name,
    anndata = list(file = file, defaults = .as_json_object(viewer$defaults))
  )))
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

# The seed of a served object's priority order: the subsampling
# specification's default, so the viewer's default subsample of a served
# object is a prefix of every array (cellucid-python's SERVED_PRIORITY_SEED).
.served_priority_seed <- 0

# The id an object declares in its cellucid record (Seurat misc$cellucid,
# metadata(x)$cellucid), or NULL.
.served_declared_id <- function(data) {
  if (!(methods::is(data, "Seurat") || methods::is(data, "SummarizedExperiment"))) return(NULL)
  id <- .object_record(data)[["id"]]
  if (.is_text(id)) id else NULL
}

# An R object, planned once and held as a store in memory, uncompressed, its
# cells in priority order (R/priority.R): a subsample the viewer opens is the
# first rows of every array, so the bytes sent follow the cells opened. A
# cell's original position is the one cellucid_prepare() gives it, by the one
# rule (.stored_order()): an object in a verified stored order (read from a
# file in priority order, its record kept) is put back in its original order
# first; any other is taken in its own order, a position column without a
# record that describes the cells naming nothing (not served, with a
# warning). Pull, push and highlights name cells by that position.
.serve_object <- function(viewer, data, extra, label, quiet, cache_bytes) {
  if (!is.null(data)) data <- .in_original_order(data)
  components <- .dataset_components(data, extra)
  components[["obs"]] <- .without_positions(components[["obs"]])
  components$priority_seed <- .served_priority_seed
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
  # By original position: the object's own order, which the viewer names.
  viewer$cell_names <- plan$cell_names
  viewer$dataset_id <- plan$dataset_id
  # The id the object declares in its cellucid record, which dataset_id =
  # may have replaced: an object declaring it is the one served under
  # viewer$dataset_id (cellucid_apply_session() of this viewer's request).
  viewer$declared_id <- .served_declared_id(data)
  viewer$dataset_name <- plan$dataset_name
  viewer$n_cells <- plan$n_cells
  viewer$n_genes <- plan$n_genes
  viewer$modalities <- plan$modalities
  viewer$defaults <- list()
  viewer$served <- list(kind = plan$kind, layout = "folder", route_name = paste0(plan$dataset_id, ".zarr"), zarr_format = 2L)
  invisible(viewer)
}

# A file or a Zarr folder, served from disk as it is (.path_order_note()
# says what a subsample of it reads).
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
  # did, and is kept either way (viewer$declared_id: the object read from it
  # declares it, and takes back a session saved under the id it was served
  # under). A field the caller names is not needed from the store: one that
  # cannot be read there is left out rather than refusing the view. An .h5ad
  # or .h5mu file's is not read here (this package has no HDF5 reader): it is
  # named after the file unless dataset_id/dataset_name are given, where
  # cellucid-python reads it from the file.
  if (classified$kind %in% c("zarr", "mudata-zarr")) {
    declared <- .zarr_store_identity(
      classified$path, classified$layout, classified$zarr_format,
      named = c(if (!is.null(extra$dataset_id)) "id", if (!is.null(extra$dataset_name)) "name")
    )
  }
  dataset_name <- extra$dataset_name %||% declared$name %||% basename(classified$path)
  dataset_id <- extra$dataset_id %||% declared$id %||%
    .label_id(sub("(\\.zarr\\.zip|\\.[A-Za-z0-9]+)$", "", basename(classified$path)))
  viewer$dataset_id <- .validate_dataset_id(dataset_id)
  viewer$declared_id <- if (.is_text(declared$id)) declared$id else NULL
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
    generated = generated,
    zarr_format = classified$zarr_format
  )
  if (!quiet && is.null(extra$dataset_id) && is.null(declared$id)) {
    .inform("Dataset id '", viewer$dataset_id, "'; name it with dataset_id = to keep sessions across files.")
  }
  note <- .path_order_note(classified)
  if (!quiet && !is.null(note)) .inform(note)
  invisible(viewer)
}

# What a subsample of a file served as it is reads, said when it is served.
# A file in priority order (cellucid_prepare()'s, prepare()'s) is read as the
# first rows of every array; R reads the order of a Zarr store from its
# metadata (uns['cellucid']['priority']), and has no reader for an HDF5 file
# (.h5ad, .h5mu, Loom, Cell Ranger) to say or to put it in priority order --
# cellucid-python serves each of those, and an AnnData or MuData Zarr store,
# in priority order. A single-stream file (an R object, a MatrixMarket
# folder) is read whole by the browser once, into a store in priority order.
# NULL when there is nothing to say.
.path_order_note <- function(classified) {
  kind <- classified$kind
  if (kind %in% c("rds", "rda", "mtx")) return(NULL)
  priority <- c("uns/cellucid/priority/.zgroup", "uns/cellucid/priority/zarr.json")
  if (identical(classified$layout, "folder") && kind %in% c("zarr", "mudata-zarr")) {
    if (any(file.exists(file.path(classified$path, priority)))) return(NULL)
  }
  if (identical(classified$layout, "zip")) {
    members <- utils::unzip(classified$path, list = TRUE)$Name
    if (any(vapply(priority, function(marker) any(endsWith(members, marker)), logical(1L)))) return(NULL)
  }
  hdf5 <- !kind %in% c("zarr", "mudata-zarr")
  paste0(
    basename(classified$path), " is served as it is. Its cells are in their own order",
    if (kind %in% c("h5ad", "h5mu")) " (unless cellucid_prepare() or prepare() wrote it)" else "",
    ": the browser reads it in place, and a subsample reads the ",
    "storage blocks that hold its cells, with the cells around them (R has no ",
    if (hdf5) "HDF5 reader to serve it in priority order). " else "reader of Zarr chunks to serve it in priority order). ",
    "For a subsample that reads only its own cells, view the object instead (read it into R, then ",
    "cellucid_view(object): served in priority order), write it once with cellucid_prepare(), or serve ",
    "the file from Python (cellucid serve <path>), which serves .h5ad, .h5mu, Loom, Cell Ranger .h5 ",
    "and outs/ folders, and AnnData and MuData Zarr stores, in priority order."
  )
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
  if (!is.null(extra$knn_key)) out$knn <- .validate_key(extra$knn_key, "knn_key")
  if (!is.null(extra$expression_source)) out$expression <- .served_expression(extra$expression_source)
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

# A non-empty key (an obsp graph's name).
.validate_key <- function(value, arg) {
  if (!nzchar(.validate_string(value, arg))) .abort(arg, " must name a key; got \"\".")
  value
}

# expression_source as the viewer names the matrix, by its place in the store:
# "X", "raw" (raw.X) or "layer:<name>" (a layer, its name holding no "/"), as
# cellucid-python's expression_source is written (_launch_defaults.py) and the
# web app reads defaults.expression (requireAnnDataLaunchDefaults).
.served_expression <- function(value) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    if (value %in% c("X", "raw")) return(if (identical(value, "raw")) "raw/X" else "X")
    layer <- sub("^layer:", "", value)
    if (startsWith(value, "layer:") && nzchar(layer) && !grepl("/", layer, fixed = TRUE)) return(paste0("layers/", layer))
  }
  .abort(
    "expression_source must be \"X\", \"raw\" or \"layer:<name>\" (a layer's name, without \"/\"); got ",
    if (is.character(value)) .format_value_list(value) else paste0("a ", class(value)[[1L]]), "."
  )
}

# Bind the server: `port`, or when the default 8765 is taken the next free
# port (where Google sign-in, which accepts only 8765, is unavailable). A port
# the caller chose is used or refused.
.bind_viewer <- function(viewer, port, port_given, quiet) {
  app <- list(
    call = function(req) .on_call(viewer, req),
    onHeaders = function(req) .on_headers(viewer, req),
    onWSOpen = function(ws) .on_server_ws_open(viewer, ws)
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
    # httpuv says "Failed to create server" whatever failed: a host this
    # machine cannot serve from is told apart from a port taken by binding
    # any port on it. Said as cellucid-python's bind_refusal() says it, with
    # R's spelling of each setting.
    host_serves <- tryCatch({
      httpuv::randomPort(host = viewer$host)
      TRUE
    }, error = function(error) FALSE)
    if (!host_serves) {
      .abort(
        "Host '", viewer$host, "' is not an address this machine can serve from. Use host = \"127.0.0.1\" to ",
        "serve only this computer, which is the default, or host = \"0.0.0.0\" to accept connections from ",
        "other machines."
      )
    }
    if (port_given) {
      .abort(
        "Port ", port, " is already in use on ", viewer$host, ", so Cellucid could not start its server there. ",
        "Another program is listening on that port, often an earlier Cellucid that was never stopped. Serve ",
        "on a different port with port = 9000, or let the operating system pick a free one with port = 0 and ",
        "use the Viewer URL it prints."
      )
    }
    .abort("Could not serve on ", viewer$host, ":", port, " (", conditionMessage(server), ").")
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
      "Port ", .default_port, " is in use; serving on ", viewer$port, ". Google Drive ",
      "files open only on port 8765."
    )
  }
  invisible(viewer)
}

.host_is_wildcard <- function(host) host %in% c("0.0.0.0", "::", "")

# Show the viewer: the IDE's Viewer pane (RStudio, Positron), else the
# browser. pane = "none", and a session with no browser to open
# (Rscript), say the address to open instead.
.open_pane <- function(viewer, pane, quiet) {
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
  if (identical(pane, "none") || !interactive()) {
    if (!quiet) .inform("Open ", viewer$viewer_url, " in a browser.")
    return(invisible(NULL))
  }
  utils::browseURL(viewer$viewer_url)
  invisible(NULL)
}

# A viewer page's address without what grants access to it: the viewerToken
# query parameter is the viewer channel's secret. The address with it is said
# once, to the person, where they are told to open it (cellucid_serve(),
# cellucid_view(pane = "none") or outside an interactive session,
# cellucid_mcp(open = ) on
# stderr). Everything that may be logged or pasted (an error, a debug report,
# a printed viewer) says this form and names viewer$viewer_url instead
# (cellucid-python's public_viewer_url(), the same two patterns). `text` is an
# address or a line naming addresses.
.public_viewer_url <- function(text) {
  text <- gsub("([?&])viewerToken=[^&#\\s]*&", "\\1", text, perl = TRUE)
  gsub("[?&]viewerToken=[^&#\\s]*(?=[#\\s]|$)", "", text, perl = TRUE)
}

print.cellucid_viewer <- function(x, ...) {
  state <- if (isTRUE(x$running)) "running" else "stopped"
  cat(
    "<Cellucid viewer> ", x$dataset_name, " (", x$dataset_id, "), ", state, "\n",
    # Without its token: a printed viewer lands in logs and knitted
    # documents. viewer$viewer_url is the address that opens it.
    "  address:  ", if (!is.null(x$widget)) "the notebook cell's output" else if (is.null(x$viewer_url)) "" else
      paste0(.public_viewer_url(x$viewer_url), " (opened by viewer$viewer_url)"), "\n",
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
    # Agent apps waiting on this viewer hear it stopped, and it is no longer
    # published to them: told before its server stops (R/agent-link.R).
    .agent_port_stopped(viewer, "the viewer was stopped")
    .agent_tell_withdrawn(.agent_withdraw(viewer), .viewer_close_seconds)
    if (!is.null(viewer$socket)) {
      # An open page is told to freeze first, as cellucid-python's stop() tells
      # it: it keeps its last view instead of losing its data source.
      try(.viewer_send(viewer, list(type = "freeze")), silent = TRUE)
    }
    viewer$running <- FALSE
    if (!is.null(viewer$socket)) {
      socket <- viewer$socket
      viewer$socket <- NULL
      # The server stops only once the page has answered the close (or did
      # not, .viewer_close_seconds): httpuv writes the freeze and the close a
      # turn of its loop apart, and stopping it first drops the close. A
      # notebook widget's channel has no server and closes at once.
      answered <- is.null(viewer$server)
      if (!answered) try(socket$onClose(function() answered <<- TRUE), silent = TRUE)
      try(socket$close(.viewer_close_codes[["goingAway"]], "The R session stopped the viewer"), silent = TRUE)
      deadline <- proc.time()[["elapsed"]] + .viewer_close_seconds
      while (!answered && proc.time()[["elapsed"]] < deadline) later::run_now(0.01)
    }
    failure <- tryCatch({
      # A notebook widget's viewer has no server.
      if (!is.null(viewer$server)) viewer$server$stop()
      NULL
    }, error = conditionMessage)
    .stop_session_requests(viewer)
    .stop_action_requests(viewer)
    .stop_session_loads(viewer)
    # Its matrices' layouts stop now, their scratch files let go.
    if (!is.null(viewer$store$layouts)) .native_layouts_stop(viewer$store$layouts)
    if (!is.null(viewer$store)) .backed_release_layouts(viewer$store)
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

# The events a handler takes: every event the viewer sends (the protocol's,
# R/viewer-protocol.R) but an agent app's own, which reach its port
# (R/agent-link.R), and "message", which catches them all. cellucid-python's
# hooks take the same names (cellucid.jupyter.HOOK_EVENTS); both are pinned to
# the shared protocol fixture. A function, not a constant: R/viewer-protocol.R
# is collated after this file.
.hook_events <- function() {
  events <- c(names(.viewer_event_fields), names(.viewer_binary_event_fields))
  c(sort(setdiff(events, .viewer_agent_events), method = "radix"), "message")
}

# One event a handler takes, said exactly: no partial match (match.arg()'s
# "sel" for "selection") and a misspelt event refused, never kept unheard.
.require_hook_event <- function(event, catch_all = TRUE) {
  events <- .hook_events()
  # "message" is no event the viewer sends: nothing waits for it.
  if (!catch_all) events <- setdiff(events, "message")
  if (!is.character(event) || length(event) != 1L || is.na(event) || !event %in% events) {
    .abort("event must be one of ", paste0("\"", events, "\"", collapse = ", "), "; got ",
           if (is.character(event)) .format_value_list(event) else paste0("a ", class(event)[[1L]]), ".")
  }
  event
}

cellucid_on <- function(viewer, event, handler) {
  .require_viewer(viewer)
  event <- .require_hook_event(event)
  if (!is.function(handler)) .abort("handler must be a function of one argument, the event.")
  viewer$hooks[[event]] <- c(viewer$hooks[[event]], handler)
  invisible(viewer)
}

cellucid_off <- function(viewer, event = NULL) {
  .require_viewer(viewer)
  if (is.null(event)) {
    viewer$hooks <- list()
  } else {
    viewer$hooks[[.require_hook_event(event)]] <- NULL
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
  if (!is.null(viewer$widget) && !done()) {
    # IRkernel hands the widget the viewer's messages between cells, on this
    # same thread: nothing the viewer says can arrive while this cell runs.
    .abort(
      "Waiting for ", what, " cannot work in the cell that waits: in a Jupyter notebook the ",
      "viewer is heard only between cells. Run the call again in a later cell, or register ",
      "a hook with cellucid_on(), which runs when the viewer answers."
    )
  }
  deadline <- proc.time()[["elapsed"]] + timeout
  while (!done()) {
    if (!isTRUE(viewer$running)) .abort("The viewer was stopped while waiting for ", what, ".")
    if (proc.time()[["elapsed"]] > deadline) {
      .abort(
        "The viewer did not answer within ", timeout, " seconds (waiting for ", what, "). ",
        "Is it open? ", if (is.null(viewer$socket)) "No viewer is connected: open the viewer's address (viewer$viewer_url)." else ""
      )
    }
    later::run_now(0.05)
  }
  invisible(TRUE)
}

# The next event of a kind, as cellucid-python's viewer.wait_for_event()
# gives it: the first to arrive after the call that `predicate` (when given)
# accepts; one received before is not returned (cellucid_on() and the
# viewer's last events have those).
cellucid_wait_event <- function(viewer, event, timeout = 30, predicate = NULL) {
  .require_viewer(viewer)
  event <- .require_hook_event(event, catch_all = FALSE)
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) || timeout <= 0) {
    .abort("timeout must be one positive number of seconds.")
  }
  if (!is.null(predicate) && !is.function(predicate)) {
    .abort("predicate must be NULL or a function of one argument, the event, returning TRUE or FALSE.")
  }
  waiter <- new.env(parent = emptyenv())
  waiter$event <- event
  waiter$predicate <- predicate
  waiter$payload <- NULL
  waiter$error <- NULL
  viewer$event_waiters <- c(viewer$event_waiters, waiter)
  on.exit(viewer$event_waiters <- Filter(function(other) !identical(other, waiter), viewer$event_waiters), add = TRUE)
  .wait_for(viewer, function() !is.null(waiter$payload) || !is.null(waiter$error), timeout, paste0("a \"", event, "\" event"))
  if (!is.null(waiter$error)) .abort("predicate failed on a \"", event, "\" event: ", waiter$error)
  waiter$payload
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
  # The report is what a person pastes when asking for help: the page's
  # address without what grants access to it (the viewer token), as the
  # page's own debug snapshot gives it; viewer$viewer_url opens it.
  public_url <- if (!is.null(viewer$viewer_url)) .public_viewer_url(viewer$viewer_url)
  report <- list(
    viewer_id = viewer$viewer_id,
    viewer_url = public_url,
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
  listed <- .datasets_document(viewer)
  if (is.null(listed)) {
    report$anndata_store_probe <- list(error = paste0(basename(viewer$served$path), " is no longer readable"))
    return(report)
  }
  file <- listed$datasets[[1L]]$anndata$file
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
      "The viewer's channel is not connected: open the viewer's address (viewer$viewer_url), ",
      "and check that any proxy in between passes WebSockets."
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

# The one way to get the viewer's session, in a notebook and outside one
# (cellucid-python's viewer.request_session()): ask without waiting, then
# collect it with cellucid_session_result() -- in a Jupyter notebook in a
# later cell, since the viewer answers between cells; where R can wait (a
# browser or RStudio viewer) the same call waits for it. The viewer answers
# with a .cellucid-session file of format 4: its history, every tracker and
# every step (R/session-file.R).
cellucid_request_session <- function(viewer) {
  .require_viewer(viewer)
  id <- paste0("s", .random_hex(6L))
  record <- new.env(parent = emptyenv())
  record$asked <- Sys.time()
  record$waiting <- FALSE
  viewer$session_requests[[id]] <- record
  .viewer_send(viewer, list(type = "requestSessionBundle", requestId = id))
  structure(list(viewer = viewer, id = id), class = "cellucid_session_request")
}

cellucid_session_result <- function(request, file = NULL, timeout = 60) {
  if (!inherits(request, "cellucid_session_request")) .abort("request must be what cellucid_request_session() returned.")
  if (!is.null(file)) file <- .validate_required_string(file, "file")
  viewer <- request$viewer
  record <- viewer$session_requests[[request$id]]
  if (is.null(record)) .abort("This session request was collected or withdrawn already; ask again with cellucid_request_session().")
  answered <- function() !is.null(record$answer) || !is.null(record$failed) || isTRUE(record$stopped)
  if (!answered()) {
    if (!is.null(viewer$widget)) {
      .abort("The viewer has not sent the session yet: in a Jupyter notebook it answers between cells. ",
             "Run cellucid_session_result() again in a later cell.")
    }
    record$waiting <- TRUE
    on.exit(record$waiting <- FALSE, add = TRUE)
    .wait_for(viewer, answered, timeout, "the session")
  }
  viewer$session_requests[[request$id]] <- NULL
  if (!is.null(record$failed)) .abort(record$failed, " Ask again with cellucid_request_session().")
  if (isTRUE(record$stopped)) .abort("The viewer was stopped before it sent the session.")
  path <- record$answer$path
  if (!is.null(file)) {
    if (!file.copy(path, path.expand(file), overwrite = TRUE)) .abort("Could not write the session to ", file, ".")
    unlink(path)
    path <- normalizePath(path.expand(file), winslash = "/")
  }
  cellucid_read_session(path)
}

# Every session request still open when the viewer stops: none will be
# answered; a file sent and not collected is removed.
.stop_session_requests <- function(viewer) {
  for (id in names(viewer$session_requests)) {
    record <- viewer$session_requests[[id]]
    if (!is.null(record$answer$path)) unlink(record$answer$path)
    record$answer <- NULL
    record$stopped <- TRUE
  }
  viewer$session_requests <- list()
  invisible(NULL)
}

print.cellucid_session_request <- function(x, ...) {
  record <- x$viewer$session_requests[[x$id]]
  state <- if (is.null(record)) "collected" else if (!is.null(record$answer)) "answered" else if (!is.null(record$failed)) "failed" else "waiting"
  cat("<cellucid session request ", x$id, ": ", state, ">\n", sep = "")
  invisible(x)
}

# --- The history, live ------------------------------------------------------------

# The viewer tells every change of its history (a step recorded, the step
# shown moved, an edit of the history's structure, a tracker started, used,
# renamed or removed, a session loaded) as a history_changed event:
# cellucid_last_history_change() is the latest, cellucid_on(viewer,
# "history_changed", handler) runs code on each. The history is read as a
# session (cellucid_request_session(), then the reading model of
# R/history.R) and changed through the history's actions, the generated
# cellucid_history_*() functions (R/actions-generated.R), as the State panel
# and agents change it. cellucid-python's viewer.history (LiveHistory),
# viewer.turn() and viewer.load_session() are the same calls.
# Hand-written help: man/cellucid_turn.Rd, man/cellucid_load_session.Rd.

cellucid_last_history_change <- function(viewer) {
  .require_viewer(viewer)
  viewer$last$history_changed
}

# A notebook turn: every action sent while `code` runs is one turn of the
# viewer's history (the State panel shows them as one group and
# cellucid_history_turn_undo() takes them back together). Returns the turn's
# id ("r-" and 12 hexadecimal digits), invisibly.
cellucid_turn <- function(viewer, label, code) {
  .require_viewer(viewer)
  if (missing(code)) .abort("cellucid_turn() runs code inside the turn: cellucid_turn(viewer, \"QC\", { ... }).")
  if (!is.character(label) || length(label) != 1L || is.na(label)) .abort("A turn's label is one string: what it does.")
  label <- gsub("[[:space:]]+", " ", trimws(enc2utf8(label)))
  if (!nzchar(label)) .abort("A turn needs a label: what it does, such as \"Normalize and cluster\".")
  if (!is.null(viewer$turn)) .abort("A turn (\"", viewer$turn$label, "\") is open already; turns do not nest.")
  id <- paste0("r-", .random_hex(6L))
  viewer$turn <- list(id = id, label = label)
  on.exit(viewer$turn <- NULL, add = TRUE)
  force(code)
  invisible(id)
}

# Push a .cellucid-session file into the viewer: its trackers are added
# beside the person's own, and the viewer shows the session's active tracker.
# Returns a request; cellucid_load_session_result() collects what the viewer
# shows then (in a Jupyter notebook, in a later cell).
cellucid_load_session <- function(viewer, path) {
  .require_viewer(viewer)
  path <- .validate_required_string(path, "path")
  if (!file.exists(path) || dir.exists(path)) .abort("No session file at ", path, ".")
  bytes <- readBin(path, "raw", n = file.size(path))
  name <- basename(path)
  id <- paste0("l", .random_hex(6L))
  frame <- .encode_viewer_load_session(id, name, bytes)
  record <- new.env(parent = emptyenv())
  record$name <- name
  record$waiting <- FALSE
  viewer$session_loads[[id]] <- record
  .viewer_send_frame(viewer, frame)
  invisible(structure(list(viewer = viewer, id = id, name = name), class = "cellucid_load_session_request"))
}

cellucid_load_session_result <- function(request, timeout = 60) {
  if (!inherits(request, "cellucid_load_session_request")) .abort("request must be what cellucid_load_session() returned.")
  viewer <- request$viewer
  record <- viewer$session_loads[[request$id]]
  if (is.null(record)) .abort("This session load was collected already.")
  answered <- function() !is.null(record$answer) || !is.null(record$failed) || isTRUE(record$stopped)
  if (!answered()) {
    if (!is.null(viewer$widget)) {
      .abort("The viewer has not loaded ", request$name, " yet: in a Jupyter notebook it answers between cells. ",
             "Run cellucid_load_session_result() again in a later cell.")
    }
    record$waiting <- TRUE
    on.exit(record$waiting <- FALSE, add = TRUE)
    .wait_for(viewer, answered, timeout, paste0("the viewer to load ", request$name))
  }
  viewer$session_loads[[request$id]] <- NULL
  if (!is.null(record$failed)) .abort("The viewer did not load ", request$name, ": ", record$failed)
  if (isTRUE(record$stopped)) .abort("The viewer was stopped before it loaded ", request$name, ".")
  record$answer[c("tracker", "head")]
}

.stop_session_loads <- function(viewer) {
  for (id in names(viewer$session_loads)) viewer$session_loads[[id]]$stopped <- TRUE
  viewer$session_loads <- list()
  invisible(NULL)
}

print.cellucid_load_session_request <- function(x, ...) {
  record <- x$viewer$session_loads[[x$id]]
  state <- if (is.null(record)) "collected" else if (!is.null(record$answer)) "answered" else if (!is.null(record$failed)) "failed" else "waiting"
  cat("<cellucid session load ", x$name, ": ", state, ">\n", sep = "")
  invisible(x)
}
