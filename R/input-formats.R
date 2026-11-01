# The input format registry: every format the Cellucid viewer opens, by id.
# inst/input_formats.json is byte for byte the web app's
# tests/fixtures/input-formats.json (written from assets/js/data/input-formats.js)
# and cellucid-python's packaged input_formats.json; the three repositories
# pin the same SHA-256 (tests/testthat/test-input-formats.R).
#
# cellucid_view(path) and cellucid_serve(path) classify a file or folder with
# it, as the web app and Python's cellucid serve do: a folder by its markers
# (a Zarr store's root metadata; MuData by its root encoding), a file by the
# signature in its first bytes (HDF5, ZIP), and by its name where the
# signature leaves the format open. The server reports the format id as the
# store's `kind`, and the browser reads the bytes and checks the content.

.input_formats_schema <- 1L

.input_formats <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      path <- system.file("input_formats.json", package = "cellucid", mustWork = TRUE)
      registry <- .read_json_file(path)
      if (!identical(registry$schema, .input_formats_schema)) {
        .abort("input_formats.json has schema ", registry$schema, "; this package reads schema 1.")
      }
      cached <<- registry
    }
    cached
  }
})

.format_entry <- function(id) {
  Filter(function(entry) identical(entry$id, id), .input_formats()$formats)[[1L]]
}

# The format a gzip-compressed copy (`x.h5ad.gz`) is a copy of, or NULL: a
# copy of a file read in parts, as the web app's registry tells it
# (`compressedFormatOfName`); a format read as one stream names its
# compressed endings itself (`.rds.gz`).
.compressed_copy_of <- function(name) {
  if (!endsWith(tolower(name), ".gz")) return(NULL)
  named <- .format_of_name(substr(name, 1L, nchar(name) - 3L), c("file", "zip"))
  if (is.null(named)) return(NULL)
  entry <- .format_entry(named$id)
  if (identical(entry$access, "ranges")) entry else NULL
}

# The format and layout of the longest name ending `name` ends with.
.format_of_name <- function(name, layouts) {
  lower <- tolower(name)
  best <- NULL
  for (entry in .input_formats()$formats) {
    for (layout in layouts) {
      for (ending in unlist(entry$names[[layout]])) {
        if (endsWith(lower, ending) && (is.null(best) || nchar(ending) > nchar(best$ending))) {
          best <- list(id = entry$id, layout = layout, ending = ending)
        }
      }
    }
  }
  best
}

.signature_of_head <- function(head) {
  for (signature in .input_formats()$signatures) {
    magic <- as.raw(strtoi(substring(signature$hex, seq(1, nchar(signature$hex), 2), seq(2, nchar(signature$hex), 2)), 16L))
    for (offset in unlist(signature$offsets)) {
      at <- offset + seq_along(magic)
      if (length(head) >= max(at) && identical(head[at], magic)) {
        return(list(layout = signature$layout, container = signature$container))
      }
    }
  }
  NULL
}

.file_endings <- function() {
  endings <- unlist(lapply(.input_formats()$formats, function(entry) {
    unlist(lapply(setdiff(names(entry$names), "folder"), function(layout) unlist(entry$names[[layout]])[1L]))
  }))
  paste(endings[!is.na(endings)], collapse = ", ")
}

# One path as the server serves it: list(path, kind, layout, zarr_format,
# ending). Files of every registered format are served as bytes; folders are
# served key by key when they are Zarr stores (AnnData or MuData).
.classify_path <- function(path) {
  resolved <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
  if (!file.exists(resolved)) .abort("Nothing exists at ", resolved, ".")
  if (dir.exists(resolved)) {
    zarr_format <- if (file.exists(file.path(resolved, "zarr.json"))) {
      3L
    } else if (file.exists(file.path(resolved, ".zgroup"))) {
      2L
    } else {
      NULL
    }
    if (is.null(zarr_format)) {
      .abort(
        resolved, " is not a Zarr store: it has no zarr.json (Zarr v3) or .zgroup (Zarr v2) ",
        "at its root. A 10x outs/ or MatrixMarket folder opens in the web app directly ",
        "(Open > Folder...); in R, read it with Seurat::Read10X() and pass the object."
      )
    }
    attrs <- .zarr_root_attributes(resolved, zarr_format)
    kind <- if (tolower(as.character(attrs[["encoding-type"]] %||% "")) == "mudata") "mudata-zarr" else "zarr"
    return(list(path = resolved, kind = kind, layout = "folder", zarr_format = zarr_format, ending = ".zarr"))
  }
  named <- .format_of_name(basename(resolved), c("file", "zip"))
  head <- local({
    connection <- file(resolved, "rb")
    on.exit(close(connection))
    readBin(connection, "raw", n = 2056L)
  })
  sniffed <- .signature_of_head(head)
  compressed <- if (is.null(named) && is.null(sniffed)) .compressed_copy_of(basename(resolved)) else NULL
  if (!is.null(compressed)) {
    # The browser reads a served file in parts; a gzip stream has none.
    .abort(
      basename(resolved), " is a gzip-compressed copy of ", compressed$label, ": it cannot be ",
      "read in parts. Decompress it (gunzip ", basename(resolved), ") and serve the ",
      "decompressed file. The web app opens the compressed copy itself, from a ",
      "link or from this computer, reading it whole once."
    )
  }
  if (is.null(named) && is.null(sniffed)) {
    .abort(
      "A file Cellucid opens is one of ", .file_endings(), "; ", basename(resolved),
      " is none of them."
    )
  }
  if (!is.null(sniffed) && (is.null(named) || !identical(.format_entry(named$id)$container, sniffed$container))) {
    candidates <- Filter(function(entry) {
      identical(entry$container, sniffed$container) && sniffed$layout %in% names(entry$names)
    }, .input_formats()$formats)
    named <- list(id = candidates[[1L]]$id, layout = sniffed$layout, ending = NULL)
  }
  # A format the registry recognises but Cellucid never opens (h5Seurat) is
  # refused here, in the registry's words (the web app and cellucid serve use
  # the same ones), rather than served for the browser to refuse.
  refusal <- .format_entry(named$id)$refusal
  if (!is.null(refusal)) .abort(basename(resolved), ": ", refusal)
  ending <- named$ending
  if (is.null(ending)) ending <- if (identical(named$layout, "zip")) ".zip" else paste0(".", tools::file_ext(resolved))
  list(path = resolved, kind = named$id, layout = named$layout, zarr_format = NULL, ending = ending)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# A JSON document's bytes, parsed. JSON is UTF-8, and its strings are read as
# UTF-8 in every session (.utf8_text()).
.parse_json_bytes <- function(bytes, what) {
  text <- .utf8_text(bytes)
  if (is.na(text)) .abort(what, " is not UTF-8 text.")
  jsonlite::parse_json(text, simplifyVector = FALSE)
}

.read_json_file <- function(path) {
  .parse_json_bytes(.read_bytes(path), path)
}

.zarr_root_attributes <- function(folder, zarr_format) {
  if (zarr_format == 3L) {
    return(.read_json_file(file.path(folder, "zarr.json"))$attributes %||% list())
  }
  path <- file.path(folder, ".zattrs")
  if (file.exists(path)) .read_json_file(path) else list()
}

# The consolidated metadata a Zarr store folder lacks, built from its
# documents (cellucid-python's consolidated_zarr_metadata): the key to answer
# and its bytes, or NULL when the store has it. A browser cannot list a
# folder, so this one document is how it learns every array and group.
.consolidated_zarr_metadata <- function(folder, zarr_format) {
  documents <- function(names) {
    found <- list()
    walk <- function(relative) {
      directory <- if (nzchar(relative)) file.path(folder, relative) else folder
      files <- list.files(directory, all.files = TRUE, no.. = TRUE)
      is_array <- ".zarray" %in% files || ("zarr.json" %in% files && identical(
        tryCatch(.read_json_file(file.path(directory, "zarr.json"))$node_type, error = function(e) NULL), "array"
      ))
      for (name in sort(files[files %in% names], method = "radix")) {
        key <- if (nzchar(relative)) paste0(relative, "/", name) else name
        found[[key]] <<- .read_json_file(file.path(directory, name))
      }
      if (is_array) return(invisible(NULL))
      for (child in .sort_code_points(files[dir.exists(file.path(directory, files))])) {
        walk(if (nzchar(relative)) paste0(relative, "/", child) else child)
      }
    }
    walk("")
    found
  }
  if (zarr_format == 3L) {
    root <- .read_json_file(file.path(folder, "zarr.json"))
    if (is.list(root$consolidated_metadata)) return(NULL)
    nodes <- documents("zarr.json")
    nodes <- nodes[names(nodes) != "zarr.json"]
    names(nodes) <- sub("/zarr\\.json$", "", names(nodes))
    root$consolidated_metadata <- list(kind = "inline", must_understand = FALSE, metadata = .as_json_object(nodes))
    return(list(key = "zarr.json", bytes = charToRaw(enc2utf8(.json_text(root)))))
  }
  if (file.exists(file.path(folder, ".zmetadata"))) return(NULL)
  metadata <- documents(c(".zgroup", ".zattrs", ".zarray"))
  text <- .json_text(list(zarr_consolidated_format = 1L, metadata = .as_json_object(metadata)))
  list(key = ".zmetadata", bytes = charToRaw(enc2utf8(text)))
}

.as_json_object <- function(value) if (length(value) == 0L) .json_object() else value

# Cells and genes of a Zarr store folder, from its metadata alone: NULL when
# a count is not stated there.
.zarr_store_counts <- function(folder, zarr_format, kind) {
  length_of <- function(key) {
    path <- if (zarr_format == 3L) file.path(folder, key, "zarr.json") else file.path(folder, key, ".zarray")
    if (!file.exists(path)) return(NULL)
    shape <- tryCatch(.read_json_file(path)$shape, error = function(e) NULL)
    if (length(shape) < 1L) NULL else as.numeric(shape[[1L]])
  }
  n_cells <- length_of("obs/_index")
  n_genes <- if (identical(kind, "mudata-zarr")) {
    attrs <- if (zarr_format == 3L) {
      tryCatch(.read_json_file(file.path(folder, "mod", "zarr.json"))$attributes, error = function(e) NULL)
    } else {
      path <- file.path(folder, "mod", ".zattrs")
      if (file.exists(path)) .read_json_file(path) else NULL
    }
    names <- unlist(attrs[["mod-order"]])
    if (length(names) == 0L) NULL else length_of(file.path("mod", .default_modality(names), "var/_index"))
  } else {
    length_of("var/_index")
  }
  list(n_cells = n_cells, n_genes = n_genes)
}

# The identity a Zarr store declares (uns['cellucid'] id and name), read from
# its metadata and two small chunks: list(id, name), each NULL when the store
# does not state it or stores it with a codec this reader does not undo (a
# store written by cellucid_prepare() or cellucid.prepare() always can be).
.zarr_store_identity <- function(path, layout, zarr_format = NULL) {
  read_entry <- if (identical(layout, "zip")) {
    listing <- tryCatch(utils::unzip(path, list = TRUE)$Name, error = function(error) {
      .abort(basename(path), " cannot be read as a ZIP file (", conditionMessage(error), ").")
    })
    folder <- tempfile("cellucid_identity_")
    on.exit(unlink(folder, recursive = TRUE), add = TRUE)
    if (is.null(zarr_format)) zarr_format <- if ("zarr.json" %in% listing) 3L else 2L
    function(key) {
      if (!key %in% listing) return(NULL)
      utils::unzip(path, files = key, exdir = folder)
      .read_bytes(file.path(folder, key))
    }
  } else {
    function(key) {
      file <- file.path(path, key)
      if (file.exists(file)) .read_bytes(file) else NULL
    }
  }
  scalar <- function(key) {
    # An absent key declares nothing; a present one that cannot be read is a
    # damaged store, said so rather than passed over for the file's name.
    if (is.null(read_entry(paste0(key, if (zarr_format == 3L) "/zarr.json" else "/.zarray")))) return(NULL)
    tryCatch({
      if (zarr_format == 3L) {
        meta <- .parse_json_bytes(read_entry(paste0(key, "/zarr.json")), paste0(key, "/zarr.json"))
        codecs <- vapply(meta$codecs, function(codec) codec$name, character(1L))
        if (!identical(codecs, c("vlen-utf8", "gzip"))) return(NULL)
        bytes <- read_entry(paste0(key, "/c"))
        bytes <- .native_inflate(bytes, gzip = TRUE)
      } else {
        meta <- .parse_json_bytes(read_entry(paste0(key, "/.zarray")), paste0(key, "/.zarray"))
        filters <- vapply(meta$filters %||% list(), function(filter) filter$id, character(1L))
        if (!identical(filters, "vlen-utf8") || !meta$compressor$id %in% c("zlib", "gzip")) return(NULL)
        bytes <- .native_inflate(read_entry(paste0(key, "/0")), gzip = identical(meta$compressor$id, "gzip"))
      }
      .decode_vlen_utf8(bytes)[[1L]]
    }, error = function(error) {
      .abort(
        "The store ", basename(path), " holds ", key, " and it cannot be read (",
        conditionMessage(error), "). Write the store again, or name the dataset with ",
        "dataset_id = and dataset_name =."
      )
    })
  }
  list(id = scalar("uns/cellucid/id"), name = scalar("uns/cellucid/name"))
}

.decode_vlen_utf8 <- function(bytes) {
  n <- readBin(bytes[1:4], "integer", size = 4L, endian = "little")
  out <- character(n)
  at <- 5L
  for (index in seq_len(n)) {
    size <- readBin(bytes[at:(at + 3L)], "integer", size = 4L, endian = "little")
    at <- at + 4L
    out[[index]] <- if (size == 0L) "" else rawToChar(bytes[at:(at + size - 1L)])
    at <- at + size
  }
  Encoding(out) <- "UTF-8"
  out
}
