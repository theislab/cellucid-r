graph_fixture <- function() {
  graph <- matrix(0, nrow = 5, ncol = 5)
  for (edge in list(c(4L, 5L), c(1L, 4L), c(2L, 4L), c(1L, 3L))) {
    graph[edge[[1L]], edge[[2L]]] <- 1
    graph[edge[[2L]], edge[[1L]]] <- 1
  }
  graph
}

connectivity_prepare_args <- function(out_dir, connectivities) {
  list(
    latent_space = matrix(seq_len(10), nrow = 5, ncol = 2),
    obs = data.frame(score = seq_len(5)),
    connectivities = connectivities,
    X_umap_2d = matrix(seq_len(10), nrow = 5, ncol = 2),
    out_dir = out_dir,
    centroid_min_points = 1,
    obs_categorical_dtype = "uint16",
    dataset_id = "exact-connectivity",
    dataset_name = "Exact connectivity",
    force = TRUE
  )
}

read_uint16_file <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  readBin(
    connection,
    what = "integer",
    size = 2L,
    signed = FALSE,
    endian = "little",
    n = 100
  )
}

read_float64_file <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  readBin(
    connection,
    what = "double",
    size = 8L,
    endian = "little",
    n = 100
  )
}

test_that("base and Matrix graphs produce one canonical non-mutating sequence", {
  graph <- graph_fixture()
  coordinates <- which(graph == 1, arr.ind = TRUE)
  pattern_graph <- Matrix::sparseMatrix(
    i = coordinates[, 1L],
    j = coordinates[, 2L],
    dims = dim(graph)
  )

  for (candidate in list(
    graph,
    graph == 1,
    Matrix::Matrix(graph, sparse = TRUE),
    pattern_graph
  )) {
    before <- candidate
    edges <- cellucid:::.validate_connectivity_edges(candidate, n_cells = 5L)

    expect_equal(edges$sources, c(0, 0, 1, 3))
    expect_equal(edges$destinations, c(2, 3, 3, 4))
    expect_equal(edges$weights, rep(1, 4))
    expect_equal(edges$n_edges, 4L)
    expect_equal(edges$max_neighbors, 3L)
    expect_identical(edges$index_dtype, "uint16")
    expect_identical(edges$index_bytes, 2L)
    expect_identical(candidate, before)
  }
})

test_that("invalid graphs are rejected without scientific reinterpretation", {
  invalid <- list(
    symmetric = list(
      matrix(c(0, 0, 1, 0), nrow = 2),
      "symmetric"
    ),
    asymmetric_weight = list(
      matrix(c(0, 1, 2, 0), nrow = 2),
      "symmetric"
    ),
    negative = list(
      matrix(c(0, -1, -1, 0), nrow = 2),
      "non-negative"
    ),
    nonfinite = list(
      matrix(c(0, NaN, NaN, 0), nrow = 2),
      "finite"
    ),
    diagonal = list(
      matrix(c(1, 0, 0, 0), nrow = 2),
      "diagonal"
    ),
    shape = list(
      matrix(0, nrow = 2, ncol = 3),
      "shape"
    )
  )

  for (case in invalid) {
    for (candidate in list(case[[1L]], Matrix::Matrix(case[[1L]], sparse = TRUE))) {
      expect_error(
        cellucid:::.validate_connectivity_edges(candidate, n_cells = 2L),
        case[[2L]]
      )
    }
  }

  implicit_diagonal <- Matrix::Diagonal(2L, x = 1)
  expect_error(
    cellucid:::.validate_connectivity_edges(
      implicit_diagonal,
      n_cells = 2L
    ),
    "diagonal"
  )

  explicit_zero <- methods::new(
    "dgTMatrix",
    i = as.integer(c(0, 1)),
    j = as.integer(c(1, 0)),
    x = c(0, 0),
    Dim = as.integer(c(2, 2))
  )
  expect_error(
    cellucid:::.validate_connectivity_edges(
      explicit_zero,
      n_cells = 2L
    ),
    "omit zero-weight coordinates"
  )
})

test_that("duplicate sparse coordinates are rejected before Matrix coalescing", {
  duplicate_graph <- methods::new(
    "dgTMatrix",
    i = as.integer(c(0, 0, 1, 1)),
    j = as.integer(c(1, 1, 0, 0)),
    x = c(1, 0, 1, 0),
    Dim = as.integer(c(2, 2))
  )

  expect_error(
    cellucid:::.validate_connectivity_edges(
      duplicate_graph,
      n_cells = 2L
    ),
    "duplicate sparse coordinates at \\(0, 1\\)"
  )
})

test_that("unsigned packers reject non-integer and out-of-range values", {
  # The byte layout these packers produce is asserted in test-byte-order.R,
  # against values that are not byte palindromes. It was asserted here against
  # 0, 65535 and 4294967295, whose bytes are all identical to each other, so
  # every permutation of the byte lanes satisfied it -- including a reversed
  # packer. What is left here is what this test is named for.
  for (invalid in list(-1, 65536)) {
    expect_error(cellucid:::.pack_uint16(invalid), "0.*65,535")
  }
  for (invalid in list(-1, 4294967296)) {
    expect_error(cellucid:::.pack_uint32(invalid), "0.*4,294,967,295")
  }
  for (invalid in list(0.5, 1.5)) {
    expect_error(cellucid:::.pack_uint16(invalid), "integers")
    expect_error(cellucid:::.pack_uint32(invalid), "integers")
  }
  for (invalid in list(TRUE, "1", factor("1"))) {
    expect_error(cellucid:::.pack_uint16(invalid), "numeric vector")
    expect_error(cellucid:::.pack_uint32(invalid), "numeric vector")
  }
  for (invalid in list(NA_real_, NaN, Inf, -Inf)) {
    expect_error(cellucid:::.pack_uint16(invalid), "finite")
    expect_error(cellucid:::.pack_uint32(invalid), "finite")
  }
})

test_that("float64 weight writer is exact for plain and gzip payloads", {
  values <- c(0.5, 1 + (2^-30), 2)
  for (compression in list(NULL, 6L)) {
    path <- tempfile("cellucid-float64-weights-")
    on.exit(unlink(c(path, paste0(path, ".gz")), force = TRUE), add = TRUE)
    written <- cellucid:::.write_float64_vector(
      path,
      values,
      compression = compression
    )
    connection <- if (is.null(compression)) {
      file(written, open = "rb")
    } else {
      gzfile(written, open = "rb")
    }
    observed <- readBin(
      connection,
      what = "double",
      size = 8L,
      endian = "little",
      n = length(values)
    )
    close(connection)
    expect_identical(observed, values)
  }

  for (invalid in list(1L, TRUE, "1", structure(1, class = "weight"), NaN, Inf)) {
    expect_error(
      cellucid:::.write_float64_vector(tempfile(), invalid),
      "finite native double vector"
    )
  }
})

test_that("empty graphs and browser index-width boundaries are exact", {
  empty <- cellucid:::.validate_connectivity_edges(
    Matrix::Matrix(0, nrow = 4, ncol = 4, sparse = TRUE),
    n_cells = 4L
  )
  expect_length(empty$sources, 0L)
  expect_length(empty$destinations, 0L)
  expect_length(empty$weights, 0L)
  expect_equal(empty$n_edges, 0L)
  expect_equal(empty$max_neighbors, 0L)

  for (case in list(
    list(n_cells = 65535L, dtype = "uint16", bytes = 2L),
    list(n_cells = 65536L, dtype = "uint16", bytes = 2L),
    list(n_cells = 65537L, dtype = "uint32", bytes = 4L)
  )) {
    n_cells <- case$n_cells
    graph <- Matrix::sparseMatrix(
      i = c(1L, n_cells),
      j = c(n_cells, 1L),
      x = c(1, 1),
      dims = c(n_cells, n_cells)
    )
    edges <- cellucid:::.validate_connectivity_edges(
      graph,
      n_cells = n_cells
    )
    expect_equal(edges$sources, 0)
    expect_equal(edges$destinations, n_cells - 1)
    expect_equal(edges$weights, 1)
    expect_identical(edges$index_dtype, case$dtype)
    expect_identical(edges$index_bytes, case$bytes)
    expect_equal(edges$max_neighbors, 1L)
  }

  expect_identical(
    cellucid:::.connectivity_index_dtype(4294967296),
    list(index_dtype = "uint32", index_bytes = 4L)
  )
  # Two separate rules refuse these, and each input is bound to the one that
  # owns it. A single "positive|uint32" alternation was satisfied by either
  # message for every input, so it could not tell the rules apart -- a count
  # rejected for exceeding the index domain and a count rejected for not being
  # a positive integer read the same to it.
  invalid <- list(
    list(value = 0L, message = "^n_cells must be one positive integer\\.$"),
    list(value = -1L, message = "^n_cells must be one positive integer\\.$"),
    list(value = TRUE, message = "^n_cells must be one positive integer\\.$"),
    list(value = 1.5, message = "^n_cells must be one positive integer\\.$"),
    list(
      value = 4294967297,
      message = "^Connectivity cannot exceed 4,294,967,296 cells"
    )
  )
  for (case in invalid) {
    expect_error(
      cellucid:::.connectivity_index_dtype(case$value),
      case$message,
      info = format(case$value)
    )
  }
})

test_that("prepared writer emits exact base and Matrix graph payloads", {
  graph <- graph_fixture()

  for (candidate in list(graph, Matrix::Matrix(graph, sparse = TRUE))) {
    out <- tempfile("cellucid-connectivity-")
    dir.create(out, recursive = TRUE)
    on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)

    do.call(cellucid_prepare, connectivity_prepare_args(out, candidate))

    manifest <- jsonlite::read_json(
      file.path(out, "connectivity_manifest.json"),
      simplifyVector = TRUE
    )
    identity <- jsonlite::read_json(
      file.path(out, "dataset_identity.json"),
      simplifyVector = TRUE
    )
    expect_equal(
      read_uint16_file(file.path(out, "connectivity", "edges.src.bin")),
      c(0, 0, 1, 3)
    )
    expect_equal(
      read_uint16_file(file.path(out, "connectivity", "edges.dst.bin")),
      c(2, 3, 3, 4)
    )
    expect_equal(
      read_float64_file(
        file.path(out, "connectivity", "edges.weights.f64.bin")
      ),
      rep(1, 4)
    )
    expect_identical(manifest$format, "edge_pairs")
    expect_equal(manifest$n_cells, 5L)
    expect_equal(manifest$n_edges, 4L)
    expect_equal(manifest$max_neighbors, 3L)
    expect_equal(manifest$index_bytes, 2L)
    expect_identical(manifest$index_dtype, "uint16")
    expect_identical(
      manifest$weightsPath,
      "connectivity/edges.weights.f64.bin"
    )
    expect_identical(manifest$weight_dtype, "float64")
    expect_equal(manifest$weight_bytes, 8L)
    expect_true(identity$stats$has_connectivity)
    expect_equal(identity$stats$n_edges, 4L)
  }
})

test_that("prepared writer preserves exact symmetric positive weights", {
  graph <- matrix(0, nrow = 5, ncol = 5)
  expected <- c(0.5, 1 + (2^-30), 2)
  for (edge in list(
    list(source = 1L, destination = 2L, weight = expected[[1L]]),
    list(source = 1L, destination = 4L, weight = expected[[2L]]),
    list(source = 3L, destination = 5L, weight = expected[[3L]])
  )) {
    graph[edge$source, edge$destination] <- edge$weight
    graph[edge$destination, edge$source] <- edge$weight
  }

  for (candidate in list(graph, Matrix::Matrix(graph, sparse = TRUE))) {
    out <- tempfile("cellucid-weighted-connectivity-")
    dir.create(out, recursive = TRUE)
    on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)

    edges <- cellucid:::.validate_connectivity_edges(candidate, n_cells = 5L)
    expect_equal(edges$sources, c(0, 0, 2))
    expect_equal(edges$destinations, c(1, 3, 4))
    expect_identical(edges$weights, expected)

    do.call(cellucid_prepare, connectivity_prepare_args(out, candidate))
    expect_identical(
      read_float64_file(
        file.path(out, "connectivity", "edges.weights.f64.bin")
      ),
      expected
    )
  }
})

test_that("compressed prepared graph declares and writes aligned gzip payloads", {
  graph <- graph_fixture()
  graph[1L, 3L] <- 1 + (2^-30)
  graph[3L, 1L] <- 1 + (2^-30)
  out <- tempfile("cellucid-compressed-connectivity-")
  dir.create(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)

  arguments <- connectivity_prepare_args(out, graph)
  arguments$compression <- 6L
  do.call(cellucid_prepare, arguments)

  manifest <- jsonlite::read_json(
    file.path(out, "connectivity_manifest.json"),
    simplifyVector = TRUE
  )
  expect_identical(
    manifest$sourcesPath,
    "connectivity/edges.src.bin.gz"
  )
  expect_identical(
    manifest$destinationsPath,
    "connectivity/edges.dst.bin.gz"
  )
  expect_identical(
    manifest$weightsPath,
    "connectivity/edges.weights.f64.bin.gz"
  )
  expect_equal(manifest$compression, 6L)

  weight_connection <- gzfile(
    file.path(out, manifest$weightsPath),
    open = "rb"
  )
  observed <- readBin(
    weight_connection,
    what = "double",
    size = 8L,
    endian = "little",
    n = manifest$n_edges
  )
  close(weight_connection)
  expect_identical(observed, c(1 + (2^-30), 1, 1, 1))
  expect_false(file.exists(
    file.path(out, "connectivity", "edges.weights.f64.bin")
  ))
})

test_that("empty prepared graph remains an advertised exact capability", {
  out <- tempfile("cellucid-empty-connectivity-")
  dir.create(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)

  do.call(
    cellucid_prepare,
    connectivity_prepare_args(out, matrix(0, nrow = 5, ncol = 5))
  )

  manifest <- jsonlite::read_json(
    file.path(out, "connectivity_manifest.json"),
    simplifyVector = TRUE
  )
  identity <- jsonlite::read_json(
    file.path(out, "dataset_identity.json"),
    simplifyVector = TRUE
  )
  expect_equal(file.info(
    file.path(out, "connectivity", "edges.src.bin")
  )$size, 0)
  expect_equal(file.info(
    file.path(out, "connectivity", "edges.dst.bin")
  )$size, 0)
  expect_equal(file.info(
    file.path(out, "connectivity", "edges.weights.f64.bin")
  )$size, 0)
  expect_equal(manifest$n_edges, 0L)
  expect_equal(manifest$max_neighbors, 0L)
  expect_true(identity$stats$has_connectivity)
  expect_equal(identity$stats$n_edges, 0L)
})

test_that("rejected graph cannot mutate graph or identity publication", {
  out <- tempfile("cellucid-rejected-connectivity-")
  graph_dir <- file.path(out, "connectivity")
  dir.create(graph_dir, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)

  sentinels <- c(
    "connectivity_manifest.json" = "manifest-before",
    "connectivity/edges.src.bin" = "sources-before",
    "connectivity/edges.dst.bin" = "destinations-before",
    "connectivity/edges.weights.f64.bin" = "weights-before",
    "dataset_identity.json" = "identity-before"
  )
  for (relative_path in names(sentinels)) {
    writeBin(
      charToRaw(sentinels[[relative_path]]),
      file.path(out, relative_path)
    )
  }

  invalid <- graph_fixture()
  invalid[1L, 3L] <- -0.25
  invalid[3L, 1L] <- -0.25
  before <- invalid
  expect_error(
    do.call(cellucid_prepare, connectivity_prepare_args(out, invalid)),
    "non-negative"
  )
  expect_identical(invalid, before)
  for (relative_path in names(sentinels)) {
    expect_identical(
      readChar(
        file.path(out, relative_path),
        nchars = nchar(sentinels[[relative_path]]),
        useBytes = TRUE
      ),
      sentinels[[relative_path]]
    )
  }
})
