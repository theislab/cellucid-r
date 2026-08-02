# The neighbourhood graph, published as edge pairs. The matrix must be
# exactly symmetric with a zero diagonal, so each undirected edge is taken
# once from the strict upper triangle and the index width follows the cell
# count.

.connectivity_triplet_values <- function(triplet) {
  if ("x" %in% names(triplet)) {
    return(triplet$x)
  }
  rep.int(1, nrow(triplet))
}

.connectivity_triplet <- function(matrix) {
  triplet_matrix <- if (inherits(matrix, "TsparseMatrix")) {
    matrix
  } else {
    methods::as(matrix, "TsparseMatrix")
  }
  rows <- triplet_matrix@i
  columns <- triplet_matrix@j
  if (length(rows) > 1L) {
    coordinate_order <- order(rows, columns, method = "radix")
    ordered_rows <- rows[coordinate_order]
    ordered_columns <- columns[coordinate_order]
    duplicate_offsets <- which(
      ordered_rows[-1L] == ordered_rows[-length(ordered_rows)] &
        ordered_columns[-1L] == ordered_columns[-length(ordered_columns)]
    )
    if (length(duplicate_offsets) > 0L) {
      duplicate_offset <- duplicate_offsets[[1L]]
      stop(
        "Connectivity contains duplicate sparse coordinates at (",
        ordered_rows[[duplicate_offset]], ", ",
        ordered_columns[[duplicate_offset]], ").",
        call. = FALSE
      )
    }
  }
  Matrix::summary(triplet_matrix)
}

.validate_connectivity_edges <- function(connectivities, n_cells) {
  info <- .connectivity_index_dtype(n_cells)
  expected_shape <- c(n_cells, n_cells)

  if (inherits(connectivities, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Matrix package is required to validate sparse connectivity matrices.")
    }
    if (!(
      inherits(connectivities, "dMatrix") ||
        inherits(connectivities, "lMatrix") ||
        inherits(connectivities, "nMatrix")
    )) {
      stop(
        "Connectivity values must be real numeric or logical values.",
        call. = FALSE
      )
    }
    if (!identical(as.numeric(dim(connectivities)), as.numeric(expected_shape))) {
      stop(
        "Connectivity matrix shape must exactly match the cell axis: expected (",
        n_cells, ", ", n_cells, "), got (",
        nrow(connectivities), ", ", ncol(connectivities), ").",
        call. = FALSE
      )
    }

    triplet <- .connectivity_triplet(connectivities)
    values <- .connectivity_triplet_values(triplet)
    if (anyNA(values) || any(!is.finite(values))) {
      stop("Connectivity values must all be finite.", call. = FALSE)
    }
    if (any(values < 0)) {
      stop(
        "Connectivity weights must all be non-negative.",
        call. = FALSE
      )
    }
    if (any(values == 0)) {
      stop(
        "Sparse connectivity storage must omit zero-weight coordinates.",
        call. = FALSE
      )
    }
    diagonal_index <- seq_len(n_cells)
    diagonal_values <- connectivities[cbind(diagonal_index, diagonal_index)]
    if (any(diagonal_values != 0)) {
      stop(
        "Connectivity diagonal values must all be exactly zero.",
        call. = FALSE
      )
    }
    if (length(Matrix::which(connectivities != Matrix::t(connectivities))) != 0L) {
      stop("Connectivity matrix must be exactly symmetric.", call. = FALSE)
    }

    upper_triplet <- .connectivity_triplet(
      Matrix::triu(connectivities, k = 1L)
    )
    weights <- as.double(.connectivity_triplet_values(upper_triplet))
    sources <- as.numeric(upper_triplet$i - 1L)
    destinations <- as.numeric(upper_triplet$j - 1L)
  } else {
    if (
      !is.matrix(connectivities) ||
        !(is.numeric(connectivities) || is.logical(connectivities)) ||
        is.object(connectivities) ||
        is.complex(connectivities)
    ) {
      stop(
        "connectivities must be a real numeric or logical matrix or a ",
        "Matrix::Matrix sparse matrix.",
        call. = FALSE
      )
    }
    if (!identical(as.numeric(dim(connectivities)), as.numeric(expected_shape))) {
      stop(
        "Connectivity matrix shape must exactly match the cell axis: expected (",
        n_cells, ", ", n_cells, "), got (",
        nrow(connectivities), ", ", ncol(connectivities), ").",
        call. = FALSE
      )
    }
    if (anyNA(connectivities) || any(!is.finite(connectivities))) {
      stop("Connectivity values must all be finite.", call. = FALSE)
    }
    if (any(connectivities < 0)) {
      stop(
        "Connectivity weights must all be non-negative.",
        call. = FALSE
      )
    }
    if (any(diag(connectivities) != 0)) {
      stop(
        "Connectivity diagonal values must all be exactly zero.",
        call. = FALSE
      )
    }
    if (!all(connectivities == t(connectivities))) {
      stop("Connectivity matrix must be exactly symmetric.", call. = FALSE)
    }

    coordinates <- which(
      upper.tri(connectivities) & connectivities > 0,
      arr.ind = TRUE
    )
    sources <- as.numeric(coordinates[, 1L] - 1L)
    destinations <- as.numeric(coordinates[, 2L] - 1L)
    weights <- as.double(connectivities[coordinates])
  }

  if (length(sources) > 0L) {
    order_index <- order(sources, destinations)
    sources <- sources[order_index]
    destinations <- destinations[order_index]
    weights <- weights[order_index]
  }
  degrees <- tabulate(
    c(sources + 1, destinations + 1),
    nbins = n_cells
  )
  max_neighbors <- if (length(degrees) > 0L) max(degrees) else 0L

  list(
    sources = sources,
    destinations = destinations,
    weights = weights,
    n_edges = length(sources),
    max_neighbors = max_neighbors,
    index_dtype = info$index_dtype,
    index_bytes = info$index_bytes
  )
}

.write_connectivity_edges <- function(connectivity_edges, out_dir, compression) {
  sources_fname <- "edges.src.bin"
  dests_fname <- "edges.dst.bin"
  weights_fname <- "edges.weights.f64.bin"
  sources_path <- file.path(out_dir, sources_fname)
  dests_path <- file.path(out_dir, dests_fname)
  weights_path <- file.path(out_dir, weights_fname)

  if (identical(connectivity_edges$index_dtype, "uint16")) {
    .write_uint16(
      sources_path,
      connectivity_edges$sources,
      compression = compression
    )
    .write_uint16(
      dests_path,
      connectivity_edges$destinations,
      compression = compression
    )
  } else {
    .write_uint32(
      sources_path,
      connectivity_edges$sources,
      compression = compression
    )
    .write_uint32(
      dests_path,
      connectivity_edges$destinations,
      compression = compression
    )
  }
  .write_float64_vector(
    weights_path,
    connectivity_edges$weights,
    compression = compression
  )
  invisible(NULL)
}

.connectivity_index_dtype <- function(n_cells) {
  if (
    !is.numeric(n_cells) ||
      length(n_cells) != 1L ||
      is.na(n_cells) ||
      !is.finite(n_cells) ||
      n_cells <= 0 ||
      n_cells != floor(n_cells)
  ) {
    stop("n_cells must be one positive integer.", call. = FALSE)
  }
  if (n_cells <= 65536) {
    list(index_dtype = "uint16", index_bytes = 2L)
  } else if (n_cells <= 4294967296) {
    list(index_dtype = "uint32", index_bytes = 4L)
  } else {
    stop(
      "Connectivity cannot exceed 4,294,967,296 cells because the current ",
      "browser contract supports at most the complete uint32 index domain.",
      call. = FALSE
    )
  }
}
