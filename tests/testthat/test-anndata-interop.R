# The store is read back by Python's anndata, the reference implementation of
# the format, and compared with what R put in, value for value. The test runs
# wherever a Python with anndata is found: the interpreter named by the
# CELLUCID_ANNDATA_PYTHON environment variable, else python3 on the PATH. It
# is skipped otherwise, and on CRAN.

anndata_python <- function() {
  candidates <- c(Sys.getenv("CELLUCID_ANNDATA_PYTHON"), Sys.which("python3"))
  for (candidate in candidates[nzchar(candidates)]) {
    status <- suppressWarnings(system2(
      candidate, c("-c", shQuote("import anndata, zarr")),
      stdout = FALSE, stderr = FALSE
    ))
    if (identical(status, 0L)) return(candidate)
  }
  NULL
}

# The outlier quantiles of the categorical obs columns with a category of 3
# cells (the rule itself is pinned by test-latent-outliers-conformance.R);
# NaN is null, as read_with_anndata.py reads it.
interop_outlier_fields <- function(args) {
  stored <- cellucid:::.stored_outlier_quantiles(args$obs, "X_pca", args$embeddings$X_pca, 3L)
  expect_identical(names(stored$fields), c("cluster", "batch", "flag"))
  lapply(stored$fields, function(field) {
    lapply(field$value, function(value) if (is.nan(value)) NULL else value)
  })
}

interop_expected <- function(args, x_dtype = "float32") {
  obs <- lapply(args$obs, function(column) {
    if (is.factor(column)) {
      codes <- as.integer(column) - 1L
      codes[is.na(codes)] <- -1L
      list(kind = "categorical", categories = levels(column), codes = codes, ordered = is.ordered(column))
    } else if (is.character(column) && anyDuplicated(column)) {
      categories <- sort(unique(column), method = "radix")
      list(kind = "categorical", categories = categories, codes = match(column, categories) - 1L, ordered = FALSE)
    } else {
      list(kind = "values", values = lapply(column, function(value) if (is.na(value)) NULL else value))
    }
  })
  rows <- function(matrix) lapply(seq_len(nrow(matrix)), function(row) as.numeric(matrix[row, ]))
  list(
    obs_names = rownames(args$gene_expression),
    var_names = rownames(args$var),
    obs = obs,
    X = rows(as.matrix(args$gene_expression)),
    X_format = "csc",
    X_dtype = x_dtype,
    X_nnz = sum(args$gene_expression != 0),
    obsm = lapply(args$embeddings, rows),
    connectivities = rows(as.matrix(args$connectivities)),
    uns = list(
      title = args$dataset_name,
      cellucid = list(
        schema = 1L,
        id = args$dataset_id,
        name = args$dataset_name,
        citation = list(text = args$citation_text),
        defaults = list(
          embedding2d = "X_umap",
          latentKey = "X_pca",
          obsKeys = list("cluster", "score"),
          centroidMinPoints = 3L
        ),
        outlierQuantiles = list(
          latentKey = "X_pca",
          minPoints = 3L,
          fields = interop_outlier_fields(args)
        )
      )
    )
  )
}

test_that("Python's anndata reads the folder and the .zarr.zip, Zarr format 2 and 3, without a warning", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  python <- anndata_python()
  skip_if(is.null(python), "no Python with anndata and zarr")
  script <- testthat::test_path("python", "read_with_anndata.py")

  args <- prepare_fixture(file.path(tempfile("cellucid_interop_"), "fixture.zarr"))
  # A stored zero is no value: anndata reads one fewer stored value.
  args$gene_expression@x[[3L]] <- 0
  args$citation_text <- "Doe et al."
  args$embedding_defaults <- c("2d" = "X_umap")
  args$latent_key <- "X_pca"
  args$obs_keys <- c("cluster", "score")
  args$centroid_min_points <- 3
  expected_path <- tempfile(fileext = ".json")
  jsonlite::write_json(
    interop_expected(args), expected_path,
    auto_unbox = TRUE, null = "null", digits = I(17)
  )

  folder <- dirname(args$out_path)
  for (zarr_format in c(2L, 3L)) {
    for (format in c("zarr", "zarr.zip")) {
      args$out_path <- file.path(folder, paste0("fixture-v", zarr_format, ".", format))
      args$format <- format
      args$zarr_format <- zarr_format
      store <- suppressMessages(do.call(cellucid_prepare, args))
      output <- suppressWarnings(system2(
        python, c(shQuote(script), shQuote(store), shQuote(expected_path)),
        stdout = TRUE, stderr = TRUE
      ))
      expect_identical(output[length(output)], "ok", info = paste(format, zarr_format, paste(output, collapse = "\n")))
    }
  }
})

test_that("Python's anndata reads integer counts above 2^24 back exactly, as float64", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  python <- anndata_python()
  skip_if(is.null(python), "no Python with anndata and zarr")
  script <- testthat::test_path("python", "read_with_anndata.py")

  args <- prepare_fixture(file.path(tempfile("cellucid_interop_"), "counts.zarr"))
  args$gene_expression@x <- c(1, 2, 3, 2^24 + 1, 5, 6, 7, 2^30, 9)
  args$citation_text <- "Doe et al."
  args$embedding_defaults <- c("2d" = "X_umap")
  args$latent_key <- "X_pca"
  args$obs_keys <- c("cluster", "score")
  args$centroid_min_points <- 3
  expected_path <- tempfile(fileext = ".json")
  jsonlite::write_json(
    interop_expected(args, x_dtype = "float64"), expected_path,
    auto_unbox = TRUE, null = "null", digits = I(17)
  )
  store <- suppressMessages(do.call(cellucid_prepare, args))
  output <- suppressWarnings(system2(
    python, c(shQuote(script), shQuote(store), shQuote(expected_path)),
    stdout = TRUE, stderr = TRUE
  ))
  expect_identical(output[length(output)], "ok", info = paste(output, collapse = "\n"))
})
