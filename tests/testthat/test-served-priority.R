# The store cellucid_view() and cellucid_serve() serve for an R object holds
# its cells in priority order (R/priority.R) under the served seed, as
# cellucid_prepare() writes them: the viewer's default subsample of a served
# object is the first rows of every array, so the bytes it reads follow the
# cells it opens. A cell's original position (what pull, push, highlights and
# sessions name) is its own: the one an object read from a file in priority
# order records, as cellucid_prepare() takes it, else the object's row.
# cellucid-python serves an AnnData the same way (_served_store.py), its
# arrays those Python's prepare() writes, which the cross-writer test holds
# to these.

served_store_of <- function(given) {
  viewer <- new.env()
  suppressMessages(cellucid:::.serve_object(viewer, NULL, given, "fixture", quiet = TRUE, cache_bytes = 268435456))
  viewer$store
}

test_that("an object is served in priority order, the arrays cellucid_prepare() writes", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(file.path(tempfile("cellucid_served_"), "data.zarr"))
  given <- args[setdiff(names(args), c("out_path", "format", "priority_seed"))]
  served <- dump_store(served_store_of(given))
  args$priority_seed <- cellucid:::.served_priority_seed
  folder <- suppressMessages(do.call(cellucid_prepare, args))
  expect_same_arrays(folder, served)

  n <- 6L
  order <- order(cellucid:::.subsample_rank(0, seq_len(n) - 1))
  expect_identical(store_array(served, "uns/cellucid/priority/seed"), 0)
  expect_equal(as.numeric(store_array(served, "obs/cellucid_original_position")), order - 1)
  expect_identical(store_array(served, "obs/_index"), paste0("cell_", order))
  # Each gene's rows ascend, and its prefix counts hold its values below each level.
  levels <- store_array(served, "uns/cellucid/priority/levels")
  indptr <- store_array(served, "X/indptr")
  rows <- store_array(served, "X/indices")
  counts <- store_array(served, "varm/cellucid_prefix_counts")
  for (gene in seq_len(length(indptr) - 1L)) {
    gene_rows <- rows[seq_len(indptr[[gene + 1L]] - indptr[[gene]]) + indptr[[gene]]]
    expect_false(is.unsorted(gene_rows, strictly = TRUE))
    for (j in seq_along(levels)) expect_equal(counts[gene, j], sum(gene_rows < levels[[j]]), info = paste(gene, j))
  }
})

test_that("a cell keeps the original position cellucid_prepare() gives it, by the one rule for a stored order", {
  # Served and prepared, an object takes its original positions by one rule
  # (.stored_order()): in a verified stored order (read from a file in
  # priority order, its record kept) it is put back first, and every cell
  # keeps its position; a position column with no record, a record that no
  # longer describes the cells (a subset of such an object, a column not the
  # seed's order) name nothing, the cells taken as they are, with the same
  # warnings. The viewer names cells by original position
  # (viewer$cell_names).
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  first <- suppressMessages(cellucid_prepare(file.path(tempfile(), "first.zarr"), "zarr", data = object,
                                             dataset_id = "s", dataset_name = "S", priority_seed = 3))
  original <- store_array(first, "obs/cellucid_original_position")
  record <- list(
    seed = store_array(first, "uns/cellucid/priority/seed"),
    specVersion = store_array(first, "uns/cellucid/priority/specVersion"),
    nObs = store_array(first, "uns/cellucid/priority/nObs"),
    levels = store_array(first, "uns/cellucid/priority/levels")
  )
  # The stored cells in their stored order, as an importer that keeps uns reads them.
  stored <- object[, original + 1]
  SummarizedExperiment::colData(stored)$cellucid_original_position <- original
  S4Vectors::metadata(stored)$cellucid <- list(priority = record)
  column_only <- stored
  S4Vectors::metadata(column_only)$cellucid <- NULL
  subset <- stored[, -1]
  swapped <- stored
  SummarizedExperiment::colData(swapped)$cellucid_original_position <- original[c(2, 1, seq.int(3, length(original)))]
  cases <- list(
    verified = list(stored, colnames(object), character()),
    `column only` = list(column_only, colnames(stored), c("has no uns", "is not written")),
    subset = list(subset, colnames(subset), c("written for 12 cells, and there are 11", "is not written")),
    `not the seed's order` = list(swapped, colnames(swapped), c("is not its seed's order", "is not written"))
  )
  for (name in names(cases)) {
    data <- cases[[name]][[1L]]
    caught <- function(expr) {
      said <- character()
      value <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
        said <<- c(said, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      list(value = value, said = said)
    }
    viewer <- new.env()
    serving <- caught(cellucid:::.serve_object(viewer, data, list(dataset_id = "s", dataset_name = "S"), "s",
                                               quiet = TRUE, cache_bytes = 268435456))
    preparing <- caught(cellucid_prepare(file.path(tempfile(), "again.zarr"), "zarr", data = data,
                                         dataset_id = "s", dataset_name = "S", priority_seed = cellucid:::.served_priority_seed))
    expect_identical(serving$said, preparing$said, info = name)
    for (part in cases[[name]][[3L]]) expect_true(any(grepl(part, serving$said, fixed = TRUE)), info = paste(name, part))
    if (length(cases[[name]][[3L]]) == 0L) expect_length(serving$said, 0L)
    expect_same_arrays(preparing$value, viewer$store)
    # Cells by original position: the object's original order when verified, else its own.
    expect_identical(viewer$cell_names, cases[[name]][[2L]], info = name)
    served <- dump_store(viewer$store)
    positions <- as.integer(store_array(served, "obs/cellucid_original_position"))
    expect_identical(viewer$cell_names[positions + 1L], store_array(served, "obs/_index"), info = name)
  }
})

# The expression matrix is never copied in priority order: its gene-major
# arrays are made from the object as it is, the cells visited (a genes x
# cells matrix) or renumbered (a cells x genes one) in that order, a dense
# one read through the order. They hold exactly the matrix with its rows
# reordered, and the prefix counts count them.
served_x <- function(folder) {
  shape <- store_json(folder, "X/.zattrs")$shape
  pointers <- store_array(folder, "X/indptr")
  rows <- store_array(folder, "X/indices")
  values <- store_array(folder, "X/data")
  dense <- matrix(0, shape[[1L]], shape[[2L]])
  for (gene in seq_len(shape[[2L]])) {
    span <- seq_len(pointers[[gene + 1L]] - pointers[[gene]]) + pointers[[gene]]
    dense[rows[span] + 1, gene] <- values[span]
  }
  dense
}

test_that("a served expression matrix is read in priority order, never copied in it", {
  skip_if_not_installed("Matrix")
  set.seed(3)
  n <- 300L
  cells_by_genes <- Matrix::rsparsematrix(n, 40L, density = 0.2, rand.x = function(k) round(runif(k, 1, 9), 1))
  dimnames(cells_by_genes) <- list(paste0("cell_", seq_len(n)), paste0("g", seq_len(40L)))
  order <- order(cellucid:::.subsample_rank(0, seq_len(n) - 1))
  base <- list(embeddings = list(X_umap = matrix(runif(n * 2), ncol = 2)), dataset_id = "x", dataset_name = "X")
  for (case in c("cells x genes", "genes x cells", "dense")) {
    given <- base
    given$gene_expression <- switch(case,
      "cells x genes" = cells_by_genes,
      "genes x cells" = cellucid:::.feature_major(Matrix::t(cells_by_genes)),
      "dense" = as.matrix(cells_by_genes)
    )
    served <- dump_store(served_store_of(given))
    if (case == "dense") {
      got <- store_array(served, "X")
    } else {
      got <- served_x(served)
      levels <- store_array(served, "uns/cellucid/priority/levels")
      counts <- store_array(served, "varm/cellucid_prefix_counts")
      want_counts <- vapply(levels, function(level) Matrix::colSums(cells_by_genes[order[seq_len(level)], , drop = FALSE] != 0), numeric(40L))
      expect_equal(unname(counts), unname(want_counts), info = case)
    }
    expect_equal(got, unname(as.matrix(cells_by_genes))[order, ], tolerance = 1e-6, info = case)
  }
})

test_that("a neighbour graph is served and written with both axes in priority order, never copied in it", {
  # Its row-compressed arrays are made a chunk at a time from the graph as
  # it is (src/chunks.c): row r is the graph's row order[r], its columns
  # renumbered and sorted. Rows straddle the chunks (2^18 entries) here.
  skip_if_not_installed("Matrix")
  set.seed(4)
  n <- 40000L
  from <- rep(seq_len(n), times = sample(0:14, n, replace = TRUE))
  to <- ((from - 1L + sample.int(n - 1L, length(from), replace = TRUE)) %% n) + 1L
  weight <- sample(1:255, length(from), replace = TRUE) / 256
  graph <- Matrix::sparseMatrix(i = c(from, to), j = c(to, from), x = c(weight, weight), dims = c(n, n))
  expect_gt(length(graph@x), 2 * 2^18)
  order <- order(cellucid:::.subsample_rank(0, seq_len(n) - 1))
  # The permuted graph is symmetric: its column-compressed arrays are the
  # row-compressed ones a store holds.
  want <- graph[order, order]
  given <- list(embeddings = list(X_umap = matrix(runif(n * 2), ncol = 2)), connectivities = graph,
                dataset_id = "g", dataset_name = "G")
  served <- dump_store(served_store_of(given))
  written <- file.path(tempfile(), "g.zarr")
  suppressMessages(do.call(cellucid_prepare, c(list(written, format = "zarr", priority_seed = 0), given)))
  for (store in list(served, written)) {
    expect_identical(as.integer(store_array(store, "obsp/connectivities/indptr")), want@p)
    expect_identical(as.integer(store_array(store, "obsp/connectivities/indices")), want@i)
    expect_identical(as.numeric(store_array(store, "obsp/connectivities/data")), want@x)
  }
})

test_that("a served object names its cells by its own order at the boundary", {
  # Pushes, selections and highlights name cells by position or name in the
  # object; the store holds them in priority order, each row's position in
  # obs['cellucid_original_position'].
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(file.path(tempfile("cellucid_served_"), "data.zarr"))
  given <- args[setdiff(names(args), c("out_path", "format", "priority_seed"))]
  viewer <- new.env()
  suppressMessages(cellucid:::.serve_object(viewer, NULL, given, "fixture", quiet = TRUE, cache_bytes = 268435456))
  expect_identical(viewer$cell_names, paste0("cell_", 1:6))
  served <- dump_store(viewer$store)
  positions <- as.integer(store_array(served, "obs/cellucid_original_position"))
  expect_false(identical(positions, 0:5))
  expect_identical(store_array(served, "obs/_index"), viewer$cell_names[positions + 1L])
})
