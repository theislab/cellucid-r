# The same object served by R (cellucid_view(), cellucid_serve()) and by
# Python (show_anndata(), serve_anndata(): cellucid-python's
# _served_store.py) holds the same bytes in every array both serve: the
# expression matrix in priority order (pointers, cells, values), the
# neighbour graph with both axes in priority order (obsp; each made chunk by
# chunk from the object as it is, on both sides), the original positions,
# the UMAP, the prefix counts and the priority record --
# each array's metadata and every chunk, byte for byte (both are
# uncompressed) -- on enough cells that every per-cell array and the
# matrix's entries span several chunks, and with the chunk sizes both
# languages' writers take from one place each (.chunk_elements_* in
# zarr-store.R; CELL/SLICED/WHOLE_CHUNK_ELEMENTS in
# cellucid.prepare_data._store) compared, so a drift of either fails here.
# Both store a matrix's values the one way each language's prepare does:
# float64 values as float32, stored zeros dropped, duplicate coordinates
# summed, a value that is not finite refused.
#
# The Python side runs in the interpreter named by CELLUCID_PYTHON (else
# python3) with cellucid and anndata importable (require_python(),
# helper-objects.R): skipped on CRAN and where there is none, an error where
# CELLUCID_PYTHON names one without them.

require_served_python <- function() {
  require_python("import cellucid._served_store, anndata", "CELLUCID_PYTHON", "cellucid and anndata")
}

# What Python serves for the object `spec` describes: each key's metadata or
# chunk digest (with the shared chunk sizes under __chunk_elements__).
python_served <- function(python, spec, arrays = served_arrays) {
  input <- tempfile(fileext = ".json")
  jsonlite::write_json(spec, input, digits = NA, auto_unbox = TRUE, na = "null")
  output <- system2(python, c(shQuote(testthat::test_path("python", "served_store.py")), shQuote(input), arrays),
                    stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  jsonlite::fromJSON(output[[length(output)]], simplifyVector = FALSE)
}

# Every key both serve holds the same metadata and bytes.
expect_same_served <- function(store, python_keys, arrays = served_arrays) {
  for (array in arrays) {
    r_keys <- grep(paste0("^", array, "/([0-9.]+|\\.zarray)$"), store$keys, value = TRUE)
    expect_gt(length(r_keys), 1L)
    py <- names(python_keys)[startsWith(names(python_keys), paste0(array, "/"))]
    expect_setequal(py, r_keys)
  }
  for (key in names(python_keys)) {
    if (endsWith(key, ".zarray")) {
      ours <- jsonlite::fromJSON(rawToChar(cellucid:::.store_get(store, key)), simplifyVector = FALSE)
      theirs <- python_keys[[key]]
      for (field in c("shape", "chunks", "dtype", "compressor", "filters", "order")) {
        expect_identical(ours[[field]], theirs[[field]], info = paste(key, field))
      }
      expect_equal(ours$fill_value, theirs$fill_value, info = key)
      next
    }
    expect_identical(cellucid:::.native_sha256(cellucid:::.store_get(store, key)), python_keys[[key]], info = key)
  }
}

r_served <- function(expression, umap, graph = NULL) {
  given <- list(gene_expression = expression, embeddings = list(X_umap = umap), connectivities = graph,
                dataset_id = "cross", dataset_name = "Cross")
  viewer <- new.env()
  suppressMessages(cellucid:::.serve_object(viewer, NULL, given, "cross", quiet = TRUE, cache_bytes = 268435456))
  viewer$store
}

served_arrays <- c(
  "X/indptr", "X/indices", "X/data", "obs/cellucid_original_position", "obsm/X_umap",
  "varm/cellucid_prefix_counts", "uns/cellucid/priority/seed", "uns/cellucid/priority/specVersion",
  "uns/cellucid/priority/nObs", "uns/cellucid/priority/levels"
)
graph_arrays <- c("obsp/connectivities/indptr", "obsp/connectivities/indices", "obsp/connectivities/data")

test_that("an object served by R and by Python holds the same bytes in every array both serve", {
  python <- require_served_python()
  skip_if_not_installed("Matrix")
  set.seed(11)
  # More cells than one per-cell chunk holds (2^16), more entries than one sliced chunk (2^14).
  n <- 70000L
  m <- 30L
  # Values float32 holds exactly, so both sides store the same numbers.
  x <- Matrix::rsparsematrix(n, m, density = 0.05, rand.x = function(k) sample(1:400, k, replace = TRUE) / 8)
  x <- methods::as(x, "CsparseMatrix")
  dimnames(x) <- list(sprintf("cell_%04d", seq_len(n)), sprintf("gene_%02d", seq_len(m)))
  umap <- matrix(sample(-800:800, 2 * n, replace = TRUE) / 16, ncol = 2)
  # A symmetric graph of more entries than one whole-read chunk (2^18) holds,
  # its weights exact in float32, its rows of different lengths.
  from <- rep(seq_len(n), times = sample(0:10, n, replace = TRUE))
  to <- ((from - 1L + sample.int(n - 1L, length(from), replace = TRUE)) %% n) + 1L
  weight <- sample(1:255, length(from), replace = TRUE) / 256
  graph <- Matrix::sparseMatrix(i = c(from, to), j = c(to, from), x = c(weight, weight), dims = c(n, n))
  arrays <- c(served_arrays, graph_arrays)
  python_keys <- python_served(python, list(
    cells = rownames(x), genes = colnames(x), p = x@p, i = x@i, x = x@x, umap = as.vector(t(umap)),
    graph = list(p = graph@p, i = graph@i, x = graph@x)
  ), arrays)
  sizes <- python_keys[["__chunk_elements__"]]
  python_keys[["__chunk_elements__"]] <- NULL
  expect_identical(
    lapply(sizes, as.numeric),
    list(cell = cellucid:::.chunk_elements_cell, sliced = cellucid:::.chunk_elements_sliced, whole = cellucid:::.chunk_elements_whole)
  )
  for (array in c("obs/cellucid_original_position", "obsm/X_umap", "X/indices", "X/data", graph_arrays[2:3])) {
    # Several chunks: their size is what is compared, not one chunk the whole length.
    expect_gt(sum(startsWith(names(python_keys), paste0(array, "/"))), 2L, label = paste(array, "chunks"))
  }
  # As a caller gives it (cells x genes, renumbered) and as an R object holds
  # it (genes x cells, its cells visited in priority order).
  for (expression in list(x, cellucid:::.feature_major(Matrix::t(x)))) {
    expect_same_served(r_served(expression, umap, graph), python_keys, arrays)
  }
  # cellucid_prepare() writes the same graph, compressed: its values are the
  # served ones.
  out <- file.path(tempfile(), "cross.zarr")
  suppressMessages(cellucid_prepare(out, gene_expression = x, embeddings = list(X_umap = umap), connectivities = graph,
                                    dataset_id = "cross", dataset_name = "Cross", format = "zarr", priority_seed = 0))
  served <- dump_store(r_served(x, umap, graph))
  for (array in graph_arrays) expect_identical(store_array(out, array), store_array(served, array), info = array)
})

test_that("float64 values, stored zeros and duplicate coordinates are stored alike by R and Python", {
  python <- require_served_python()
  skip_if_not_installed("Matrix")
  set.seed(12)
  n <- 3000L
  m <- 12L
  k <- 4000L
  rows <- sample.int(n, k, replace = TRUE)
  cols <- sample.int(m, k, replace = TRUE)
  # float64 values float32 rounds (six decimals: the same doubles once read
  # from JSON); a tenth stored zeros; some coordinates repeated.
  values <- round(stats::runif(k) * 10, 6)
  values[sample.int(k, k %/% 10)] <- 0
  rows <- c(rows, rows[1:200])
  cols <- c(cols, cols[1:200])
  values <- c(values, round(stats::runif(200), 6))
  umap <- matrix(sample(-800:800, 2 * n, replace = TRUE) / 16, ncol = 2)
  cells <- sprintf("cell_%04d", seq_len(n))
  genes <- sprintf("gene_%02d", seq_len(m))
  python_keys <- python_served(python, list(
    cells = cells, genes = genes, rows = rows - 1L, cols = cols - 1L, values = values, umap = as.vector(t(umap))
  ))
  python_keys[["__chunk_elements__"]] <- NULL
  # A dgTMatrix keeps duplicates and zeros as given; R's prepare sums and drops them.
  x <- methods::new("dgTMatrix", i = rows - 1L, j = cols - 1L, x = values, Dim = c(n, m), Dimnames = list(cells, genes))
  expect_identical(jsonlite::fromJSON(rawToChar(cellucid:::.store_get(r_served(x, umap), "X/data/.zarray")))$dtype, "<f4")
  expect_same_served(r_served(x, umap), python_keys)

  # A value that is not finite is refused by both.
  values[[5]] <- NaN
  refused <- python_served(python, list(
    cells = cells, genes = genes, rows = rows - 1L, cols = cols - 1L, values = values, umap = as.vector(t(umap))
  ))
  expect_match(refused[["__refused__"]], "finite")
  x@x[[5]] <- NaN
  expect_error(r_served(x, umap), "finite|NaN|missing")
})
