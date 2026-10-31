# A served expression matrix is laid out gene-major in the background, on
# every core (src/chunks.c, cellucid_layouts_start()). A chunk asked for before
# it ends is made from the layout's progress (cellucid_layout_block()): the
# values already laid out are copied and only the cells (or genes) not laid out
# yet are passed over, the bytes the layout will hold. These tests catch a
# layout halfway (`park`: every thread stands that many units into its share)
# and count what a request passed over and copied, never the time it took.

# A genes x cells matrix of `n_cells` cells, each holding `per_cell` values
# at genes spread over `n_genes`, every value distinct.
progress_matrix <- function(n_genes, n_cells, per_cell) {
  rows <- unlist(lapply(seq_len(n_cells), function(cell) sort((cell * 7L + seq_len(per_cell) * 13L) %% n_genes)))
  Matrix::sparseMatrix(
    i = rows, p = seq(0, by = per_cell, length.out = n_cells + 1L), x = seq_along(rows) / 7,
    dims = c(n_genes, n_cells), index1 = FALSE
  )
}

# Wait (at most 30 s) until layout `index` has `units` units laid out.
wait_units <- function(layouts, index, units) {
  deadline <- proc.time()[["elapsed"]] + 30
  repeat {
    state <- cellucid:::.native_layout_state(layouts, index)
    if (state[[2L]] >= units || state[[1L]] >= 2) return(state)
    if (proc.time()[["elapsed"]] > deadline) stop("The layout did not reach ", units, " units.")
    Sys.sleep(0.005)
  }
}

block_matches <- function(block, cells, values, range, info) {
  positions <- seq_len(range[[2]] - range[[1]]) + range[[1]]
  expect_identical(readBin(block[[1]], "integer", n = length(positions), size = 4L, endian = "little"), cells[positions], info = info)
  expect_identical(readBin(block[[2]], "double", n = length(positions), size = 8L, endian = "little"), values[positions], info = info)
}

test_that("a genes x cells matrix laid out halfway answers from its progress, passing over only the cells left", {
  skip_if_not_installed("Matrix")
  matrix <- progress_matrix(n_genes = 37L, n_cells = 3000L, per_cell = 9L)
  set.seed(3)
  visit <- sample.int(ncol(matrix)) - 1L
  map <- seq_len(nrow(matrix)) - 1L
  pointers <- cellucid:::.native_row_pointers(matrix@i, map, nrow(matrix))
  # The gene-major arrays, cells in priority order (output row c = column visit[c]).
  transposed <- Matrix::t(matrix[, visit + 1L])
  job <- list(0L, matrix@i, matrix@p, matrix@x, map, map, pointers, visit, 4L, "<f8")
  # Two threads, each standing 1,024 cells (four steps) into its 1,500:
  # cells 0..1023 and 1500..2523 laid out, 952 left.
  layouts <- cellucid:::.native_layouts_start(list(job), tempfile(c("i", "v")), threads = 2L, park = 1024)
  on.exit(cellucid:::.native_layouts_stop(layouts))
  expect_identical(wait_units(layouts, 0L, 2048)[[2L]], 2048)
  laid <- function(range) {
    cells <- transposed@i[seq_len(range[[2]] - range[[1]]) + range[[1]]]
    sum(cells < 1024 | (cells >= 1500 & cells < 2524))
  }
  n <- length(transposed@x)
  passed <- 0
  copied <- 0
  # Boundaries inside genes, on gene starts, the first and last values.
  for (range in list(c(0, n), c(0, 1), c(n - 1, n), c(17, 950), c(pointers[[4]], pointers[[5]]))) {
    block <- cellucid:::.native_layout_block(layouts, 0L, range[[1]], range[[2]])
    block_matches(block, transposed@i, transposed@x, range, paste(range, collapse = "-"))
    passed <- passed + 952
    copied <- copied + laid(range)
    state <- cellucid:::.native_layout_state(layouts, 0L)
    expect_identical(state[[1L]], 1)
    expect_identical(state[[5L]], passed, info = paste(range, collapse = "-"))
    expect_identical(state[[6L]], copied, info = paste(range, collapse = "-"))
  }
  # Let go, laid out: every block a copy, nothing passed over.
  cellucid:::.native_layouts_unpark(layouts)
  expect_identical(wait_units(layouts, 0L, 3000)[[1L]], 2)
  block <- cellucid:::.native_layout_block(layouts, 0L, 5, 800)
  block_matches(block, transposed@i, transposed@x, c(5, 800), "laid out")
  state <- cellucid:::.native_layout_state(layouts, 0L)
  expect_identical(state[[5L]], passed)
  expect_identical(state[[6L]], copied + 795)
  expect_identical(cellucid:::.native_layout_slice(layouts, 0L, 1L, 0, 4 * n, 4 * n),
                   writeBin(transposed@i, raw(), size = 4L, endian = "little"))
})

test_that("a matrix waiting behind another is made directly on every thread, the one laid out held", {
  skip_if_not_installed("Matrix")
  matrix <- progress_matrix(n_genes = 23L, n_cells = 1200L, per_cell = 5L)
  map <- seq_len(nrow(matrix)) - 1L
  pointers <- cellucid:::.native_row_pointers(matrix@i, map, nrow(matrix))
  transposed <- Matrix::t(matrix)
  job <- list(0L, matrix@i, matrix@p, matrix@x, map, map, pointers, integer(0), 4L, "<f8")
  # The first job parked: the second waits behind it.
  layouts <- cellucid:::.native_layouts_start(list(job, job), tempfile(c("i", "v", "i", "v")), threads = 2L, park = 100)
  on.exit(cellucid:::.native_layouts_stop(layouts))
  wait_units(layouts, 0L, 200)
  expect_identical(cellucid:::.native_layout_state(layouts, 1L)[[1L]], 0)
  block <- cellucid:::.native_layout_block(layouts, 1L, 3, 4000)
  block_matches(block, transposed@i, transposed@x, c(3, 4000), "waiting")
  # Two threads: every cell counted, then every cell placed.
  expect_identical(cellucid:::.native_layout_state(layouts, 1L)[5:6], c(2400, 0))
  cellucid:::.native_layouts_unpark(layouts)
  cellucid:::.store_laid_out(list(layouts = layouts, layout_labels = c("a", "b")))
  expect_identical(cellucid:::.native_layout_slice(layouts, 1L, 2L, 0, 8 * length(transposed@x), 8 * length(transposed@x)),
                   writeBin(transposed@x, raw(), size = 8L, endian = "little"))
})

test_that("a cells x genes matrix in priority order laid out halfway copies its columns laid out and sorts only the rest", {
  skip_if_not_installed("Matrix")
  # Cells x genes: a column a gene, its rows renumbered by priority.
  matrix <- Matrix::t(progress_matrix(n_genes = 40L, n_cells = 900L, per_cell = 6L))
  set.seed(5)
  new_row <- sample.int(nrow(matrix)) - 1L
  expected <- matrix[order(new_row), ]
  job <- list(1L, matrix@i, matrix@p, matrix@x, new_row, integer(0), matrix@p, 4L, "<f8")
  # Two threads, each standing 5 columns into its share (whole columns, about
  # half the values each).
  layouts <- cellucid:::.native_layouts_start(list(job), tempfile(c("i", "v")), threads = 2L, park = 5)
  on.exit(cellucid:::.native_layouts_stop(layouts))
  wait_units(layouts, 0L, 10)
  segment_start <- sum(matrix@p[-1] <= length(matrix@x) %/% 2)
  laid_columns <- c(0:4, segment_start + 0:4)
  # Every column's values: those of laid columns copied, the rest sorted here.
  for (range in list(c(0, length(matrix@x)), c(7, 1000), c(matrix@p[[3]], matrix@p[[23]]))) {
    before <- cellucid:::.native_layout_state(layouts, 0L)
    block <- cellucid:::.native_layout_block(layouts, 0L, range[[1]], range[[2]])
    block_matches(block, expected@i, expected@x, range, paste(range, collapse = "-"))
    after <- cellucid:::.native_layout_state(layouts, 0L)
    columns <- which(matrix@p[-1] > range[[1]] & matrix@p[-length(matrix@p)] < range[[2]]) - 1L
    sizes <- diff(matrix@p)[columns + 1L]
    inside <- pmin(matrix@p[columns + 2L], range[[2]]) - pmax(matrix@p[columns + 1L], range[[1]])
    laid <- columns %in% laid_columns
    expect_identical(after[[5L]] - before[[5L]], as.double(sum(sizes[!laid])), info = paste(range, collapse = "-"))
    expect_identical(after[[6L]] - before[[6L]], as.double(sum(inside[laid])), info = paste(range, collapse = "-"))
  }
})

test_that("a written store's gene-major block is the same bytes on any number of threads", {
  skip_if_not_installed("Matrix")
  matrix <- progress_matrix(n_genes = 31L, n_cells = 700L, per_cell = 7L)
  set.seed(2)
  visit <- sample.int(ncol(matrix)) - 1L
  for (kept in list(seq_len(nrow(matrix)), c(30L, 2L, 17L, 3L, 9L))) {
    map <- rep.int(-1L, nrow(matrix))
    map[kept] <- seq_along(kept) - 1L
    pointers <- cellucid:::.native_row_pointers(matrix@i, map, length(kept))
    n <- pointers[[length(pointers)]]
    one <- cellucid:::.native_gene_major_block(matrix, map, kept - 1L, pointers, 0, n, 4L, "<f4", visit, threads = 1L)
    # Two threads at most under R CMD check (CRAN's limit), more elsewhere.
    for (threads in if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2L else c(2L, 3L, 8L)) {
      expect_identical(cellucid:::.native_gene_major_block(matrix, map, kept - 1L, pointers, 0, n, 4L, "<f4", visit, threads = threads),
                       one, info = paste(length(kept), threads))
      expect_identical(cellucid:::.native_gene_major_block(matrix, map, kept - 1L, pointers, 11, n - 3, 4L, "<f4", visit, threads = threads),
                       cellucid:::.native_gene_major_block(matrix, map, kept - 1L, pointers, 11, n - 3, 4L, "<f4", visit, threads = 1L),
                       info = paste(length(kept), threads, "inside"))
    }
  }
})

test_that("a chunk's indices and values asked for during a layout are made by one block, its other ranges by none", {
  store <- new.env(parent = emptyenv())
  store$layout_jobs <- list()
  store$layout_labels <- character()
  blocks <- 0L
  made <- function(first, last) {
    n <- last - first
    list(as.raw(seq_len(n * 4) %% 256), as.raw(rev(seq_len(n * 4)) %% 256))
  }
  arrays <- cellucid:::.laid_out_arrays(store, "X", 100, c(4, 4), function() list())
  local_mocked_bindings(
    .layout_done = function(store, index) FALSE,
    .native_layout_block = function(layouts, index, first, last) {
      blocks <<- blocks + 1L
      made(first, last)
    }
  )
  # Chunks of 10 values: positions 11..15 lie in the second, 11..20.
  indices <- arrays$indices$chunk(11, 5, 5, chunk_first = 11, chunk_last = 20)
  values <- arrays$values$chunk(11, 5, 5, chunk_first = 11, chunk_last = 20)
  expect_identical(blocks, 1L)
  whole <- made(10, 20)
  expect_identical(indices, whole[[1L]][1:20])
  expect_identical(values, whole[[2L]][1:20])
  # Another range of the same chunk: no block; of the next chunk: one.
  expect_identical(arrays$values$chunk(16, 5, 5, chunk_first = 11, chunk_last = 20), whole[[2L]][21:40])
  expect_identical(blocks, 1L)
  arrays$values$chunk(21, 5, 5, chunk_first = 21, chunk_last = 30)
  expect_identical(blocks, 2L)
})

test_that("every size taken from the machine's memory starts from one rule, 8 GiB where the system does not say", {
  for (machine in list(NA_real_, 0, 16 * 2^30)) {
    local_mocked_bindings(.native_physical_memory = function() machine, .package = "cellucid")
    expected <- if (is.na(machine) || machine <= 0) 8 * 2^30 else machine
    expect_identical(cellucid:::.machine_memory_bytes(), expected)
    expect_identical(cellucid:::.default_cache_bytes(), floor(expected / 64))
    expect_identical(cellucid:::.backed_group_values(8), floor(expected / 4 / 8))
  }
})
