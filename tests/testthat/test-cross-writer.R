# The R and the Python writer produce the same AnnData from the same input.
# Each case is written here with cellucid_prepare(), in Zarr format 2 and in
# format 3, and, from a JSON description of the very same R objects, with
# Python's cellucid.prepare() in the same format;
# python/cross_writer.py then compares the two stores element by element:
# encodings, dtypes, shapes, chunks, codecs, fill values and every value.
# The test runs wherever a Python with cellucid, anndata and zarr is found:
# the interpreter named by the CELLUCID_PYTHON environment variable, else
# python3 on the PATH (require_python(), helper-objects.R). It is skipped
# otherwise, and on CRAN; where CELLUCID_PYTHON names an interpreter without
# them it fails.

cellucid_python <- function() {
  require_python("import cellucid, anndata, zarr, scipy", "CELLUCID_PYTHON", "cellucid, anndata and zarr")
}

# --- The R input, described in JSON -----------------------------------------

cross_writer_matrix <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    x <- methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix")
    return(list(kind = "sparse", shape = I(dim(x)), p = I(x@p), i = I(x@i), x = I(x@x)))
  }
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  list(kind = "dense", type = typeof(x), shape = I(dim(x)), values = I(as.vector(t(x))))
}

cross_writer_column <- function(name, values) {
  if (is.factor(values)) {
    return(list(
      name = name, type = "factor", values = I(as.character(values)),
      levels = I(levels(values)), ordered = is.ordered(values)
    ))
  }
  list(name = name, type = typeof(values), values = I(values))
}

cross_writer_frame <- function(frame, index) {
  list(
    index = if (is.null(index)) NULL else I(index),
    columns = I(unname(Map(cross_writer_column, names(frame), frame)))
  )
}

cross_writer_input <- function(args) {
  embeddings <- lapply(args$embeddings, cross_writer_matrix)
  supplied_names <- function(x) {
    if (is.null(x)) NULL else if (is.data.frame(x) && .row_names_info(x) < 0L) NULL else rownames(x)
  }
  cell_names <- Find(Negate(is.null), c(
    list(supplied_names(args$obs), supplied_names(args$gene_expression)),
    lapply(args$embeddings, rownames)
  ))
  obs <- if (is.null(args$obs)) data.frame(row = seq_len(nrow(args$embeddings[[1L]])))[, 0L, drop = FALSE] else args$obs
  spec <- list(
    arguments = args[intersect(names(args), c(
      "dataset_id", "dataset_name", "dataset_description", "citation_text",
      "citation_url", "source_name", "source_url", "latent_key",
      "var_gene_id_column", "vector_field_default", "normalize_embeddings",
      "centroid_min_points"
    ))],
    embeddings = embeddings,
    obs = cross_writer_frame(obs, cell_names)
  )
  if (!is.null(args$obs_keys)) spec$arguments$obs_keys <- I(args$obs_keys)
  if (!is.null(args$gene_identifiers)) spec$arguments$gene_identifiers <- I(args$gene_identifiers)
  if (!is.null(args$centroid_min_points)) spec$arguments$centroid_min_points <- as.integer(args$centroid_min_points)
  if (!is.null(args$gene_expression)) {
    spec$gene_expression <- cross_writer_matrix(args$gene_expression)
    var <- if (is.null(args$var)) data.frame(row = seq_len(ncol(args$gene_expression)))[, 0L, drop = FALSE] else args$var
    gene_names <- supplied_names(args$var)
    if (is.null(gene_names)) gene_names <- colnames(args$gene_expression)
    rownames(var) <- NULL
    spec$var <- cross_writer_frame(var, gene_names)
  }
  if (!is.null(args$connectivities)) spec$connectivities <- cross_writer_matrix(args$connectivities)
  if (!is.null(args$vector_fields)) {
    spec$vector_fields <- lapply(args$vector_fields, function(bindings) lapply(bindings, cross_writer_matrix))
  }
  if (!is.null(args$embedding_defaults)) {
    spec$embedding_defaults <- lapply(as.list(args$embedding_defaults), function(value) {
      if (is.list(value)) list(key = value$key, columns = I(as.integer(value$columns) - 1L)) else value
    })
  }
  spec
}

# --- The cases --------------------------------------------------------------

# Counts with stored zeros: the zeros are dropped by both writers.
with_stored_zeros <- function(x) {
  x@x[c(1L, 4L)] <- 0
  x
}

cross_writer_cases <- function() {
  base <- prepare_fixture("unused")
  base$out_path <- NULL
  counts <- with_stored_zeros(base$gene_expression)
  counts@x <- round(counts@x * 2)
  big_counts <- counts
  big_counts@x[[2L]] <- 2^24 + 1

  many_levels <- factor(rep(sprintf("type_%03d", 1:127), length.out = 6L), levels = sprintf("type_%03d", 1:127))
  full <- base
  full$gene_expression <- with_stored_zeros(base$gene_expression)
  full$obs$many <- many_levels
  full$obs$stage <- factor(c("early", "late", "mid", "early", NA, "late"), levels = c("early", "mid", "late"), ordered = TRUE)
  full$obs$whole <- c(1, 2, 3, 4, 5, 6)
  full$obs$note <- c("a", NA, "b", "c", "d", "e")
  # Outlier quantiles with centroid_min_points = 3 for a factor, a logical
  # and a character column, each with a category of 3 cells and with NA,
  # which belongs to no category.
  full$obs$grade <- factor(c("lo", "hi", "lo", NA, "lo", "hi"), levels = c("lo", "hi", "top"))
  full$obs$seen <- c(TRUE, TRUE, NA, TRUE, FALSE, NA)
  full$obs$tissue <- c("gut", NA, "gut", "gut", "lung", NA)
  full$vector_fields <- list(velocity = list(X_umap = base$embeddings$X_umap / 2))
  full$dataset_description <- "Six cells."
  full$citation_text <- "Doe et al."
  full$citation_url <- "https://example.org/doe"
  full$source_name <- "A lab"
  full$embedding_defaults <- list("2d" = "X_umap", "3d" = list(key = "X_pca", columns = c(3, 1, 2)))
  full$latent_key <- "X_pca"
  full$obs_keys <- c("cluster", "score")
  full$var_gene_id_column <- "symbol"
  full$vector_field_default <- "velocity"
  full$normalize_embeddings <- FALSE
  full$centroid_min_points <- 3

  dense_names <- list(NULL, c("g1", "g2", "g3"))
  embedding <- list(X_umap = matrix(seq_len(8) / 4, ncol = 2))
  latent_embeddings <- c(embedding, list(X_pca = matrix(c(0.5, 1, 2, 4, 3, 1, 0, 2, 7, 1, 1, 9), ncol = 3)))
  list(
    "every kind of input" = full,
    "integer counts float32 holds exactly" = modifyList(base, list(gene_expression = counts)),
    "integer counts above 2^24" = modifyList(base, list(gene_expression = big_counts)),
    "a count above 2^24 in a gene gene_identifiers leaves out" = modifyList(
      base, list(gene_expression = big_counts, gene_identifiers = c("ENSG2", "ENSG3"))
    ),
    "dense fractions" = list(
      embeddings = embedding, dataset_id = "dense", dataset_name = "Dense",
      gene_expression = matrix(c(0, 0.5, 1, 1.5, 0, 2.5, 3, 0, 0.25, 0, 0, 1) / 3, nrow = 4, dimnames = dense_names)
    ),
    "dense integer counts above 2^24" = list(
      embeddings = embedding, dataset_id = "dense", dataset_name = "Dense",
      gene_expression = matrix(c(0L, 1L, 2L, 16777217L, 0L, 3L, 4L, 0L, 5L, 0L, 0L, 6L), nrow = 4, dimnames = dense_names)
    ),
    "embeddings alone" = list(embeddings = embedding, dataset_id = "bare", dataset_name = "Bare"),
    # A latent key with no obs column, and with no category of 10 cells:
    # outlierQuantiles with an empty fields group.
    "a latent key and no obs" = list(
      embeddings = latent_embeddings, dataset_id = "latent", dataset_name = "Latent", latent_key = "X_pca"
    ),
    "a latent key and no category of 10 cells" = list(
      embeddings = latent_embeddings, dataset_id = "latent", dataset_name = "Latent", latent_key = "X_pca",
      obs = data.frame(group = factor(c("a", "b", "a", "b")), flag = c(TRUE, FALSE, NA, TRUE))
    )
  )
}

test_that("R and Python write the same AnnData from the same input", {
  skip_if_not_installed("Matrix")
  python <- cellucid_python()
  script <- testthat::test_path("python", "cross_writer.py")
  root <- tempfile("cellucid_cross_")
  dir.create(root)

  cases <- cross_writer_cases()
  pairs <- character()
  for (zarr_format in c(2L, 3L)) {
    for (index in seq_along(cases)) {
      args <- cases[[index]]
      input <- file.path(root, paste0("case", index, ".json"))
      jsonlite::write_json(
        cross_writer_input(args), input,
        auto_unbox = TRUE, null = "null", na = "null", digits = I(17)
      )
      args$out_path <- file.path(root, paste0("case", index, "-v", zarr_format, ".zarr"))
      args$format <- "zarr"
      args$zarr_format <- zarr_format
      # Both writers at their default: the cells in priority order.
      args$priority_seed <- NULL
      store <- suppressMessages(do.call(cellucid_prepare, args))
      pairs <- c(pairs, store, input)
    }
  }
  output <- suppressWarnings(system2(python, shQuote(c(script, pairs)), stdout = TRUE, stderr = TRUE))
  expect_identical(output[length(output)], "ok", info = paste(output, collapse = "\n"))

  # The cases exercise what they are named for.
  stores <- pairs[c(TRUE, FALSE)][seq_along(cases)]
  dtype <- function(case, key) array_meta(stores[[case]], key)$dtype
  expect_identical(dtype(1L, "obs/many/codes"), "<i2")
  expect_identical(dtype(2L, "X/data"), "<f4")
  expect_identical(dtype(3L, "X/data"), "<f8")
  expect_identical(dtype(4L, "X/data"), "<f4")
  expect_identical(dtype(5L, "X"), "<f4")
  expect_identical(dtype(6L, "X"), "<f8")
  expect_length(store_array(stores[[1L]], "X/data"), length(cases[[1L]]$gene_expression@x) - 2L)
  # Outlier quantiles for the columns with a category of 3 cells, and only
  # for them: not for the numbers, and not for maybe (2, 2, NA), stage (2, 1,
  # 2, NA), note and barcode (every value once) or many (127 levels, 6 cells).
  fields <- file.path(stores[[1L]], "uns/cellucid/outlierQuantiles/fields")
  expect_identical(
    list.dirs(fields, full.names = FALSE, recursive = FALSE),
    sort(c("cluster", "batch", "flag", "grade", "seen", "tissue"))
  )
  expect_identical(sum(is.nan(store_array(fields, "grade"))), 3L)
  expect_false(dir.exists(file.path(stores[[2L]], "uns/cellucid/outlierQuantiles")))
  for (case in c("a latent key and no obs", "a latent key and no category of 10 cells")) {
    outliers <- file.path(stores[[match(case, names(cases))]], "uns/cellucid/outlierQuantiles")
    expect_identical(store_array(outliers, "latentKey"), "X_pca", info = case)
    expect_identical(store_array(outliers, "minPoints"), 10, info = case)
    expect_identical(list.dirs(file.path(outliers, "fields"), full.names = FALSE, recursive = FALSE), character(0), info = case)
  }
})

# A multimodal object: R's Cellucid-ready MuData against Python's
# cellucid.prepare(data = mdata) of the same MuData. It runs wherever a Python
# with cellucid (with MuData support) and mudata is found:
# CELLUCID_MUDATA_PYTHON, else CELLUCID_PYTHON, else python3; it fails where
# CELLUCID_MUDATA_PYTHON names an interpreter without them.
cross_writer_mudata_python <- function() {
  require_python(
    "import cellucid, mudata, inspect, sys; sys.exit(0 if 'data' in inspect.signature(cellucid.prepare).parameters else 1)",
    c("CELLUCID_MUDATA_PYTHON", "CELLUCID_PYTHON"), "cellucid (MuData support) and mudata"
  )
}

test_that("R and Python write the same MuData from the same multimodal input", {
  skip_without_sce()
  python <- cross_writer_mudata_python()
  root <- tempfile("cellucid_cross_mudata_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  arguments <- character()
  for (name in names(mudata_cases())) {
    case <- mudata_cases()[[name]]
    write_mudata_case(case$object(), name, root, case$zarr_format)
    arguments <- c(arguments, file.path(root, paste0(name, c(".zarr", ".json"))))
  }
  output <- suppressWarnings(system2(
    python, shQuote(c(testthat::test_path("python", "cross_writer_mudata.py"), arguments)),
    stdout = TRUE, stderr = TRUE
  ))
  expect_identical(output[length(output)], "ok", info = paste(output, collapse = "\n"))
})

# A store's file as the writer means it: a metadata document byte for byte, a
# chunk by its bytes once inflated (a Blosc frame decoded). The deflate data
# is the platform's zlib's, which differs between zlib releases for some
# inputs (1.3.2 against 1.2.12): the bytes a writer stores are compared, not
# the zlib that compressed them (R/zarr-store.R).
stored_content <- function(path) {
  bytes <- readBin(path, "raw", file.size(path))
  if (basename(path) %in% c(".zarray", ".zattrs", ".zgroup", ".zmetadata", "zarr.json") || grepl("\\.json$", path)) {
    return(bytes)
  }
  if (length(bytes) >= 16L && bytes[[1L]] == as.raw(2L) && bytes[[2L]] == as.raw(1L)) return(blosc_decode(bytes))
  inflated <- tryCatch(memDecompress(bytes, type = "gzip"), error = function(error) NULL)
  if (is.null(inflated)) bytes else inflated
}

test_that("the MuData cases cellucid-python checks itself against are what this package writes", {
  skip_on_cran()
  skip_without_sce()
  fixtures <- testthat::test_path("fixtures", "cross_writer_mudata")
  # Kept in the repository for cellucid-python, not in the built package
  # (.Rbuildignore: Zarr's dot-files and long paths are not package files).
  skip_if_not(dir.exists(fixtures), "the MuData cross-writer cases are in the source repository only")
  root <- tempfile("cellucid_mudata_cases_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  for (name in names(mudata_cases())) {
    case <- mudata_cases()[[name]]
    write_mudata_case(case$object(), name, root, case$zarr_format)
    for (ending in c(".zarr", ".json")) {
      written <- list.files(file.path(root, paste0(name, ending)), recursive = TRUE, all.files = TRUE)
      kept <- list.files(file.path(fixtures, paste0(name, ending)), recursive = TRUE, all.files = TRUE)
      if (identical(ending, ".json")) {
        written <- kept <- ""
      }
      expect_identical(sort(written), sort(kept), info = name)
      paths <- if (identical(ending, ".json")) paste0(name, ending) else file.path(paste0(name, ending), written)
      for (path in paths) {
        expect_identical(
          stored_content(file.path(root, path)), stored_content(file.path(fixtures, path)),
          info = paste(name, path, "(regenerate with fixtures/cross_writer_mudata/generate.R)")
        )
      }
    }
  }
})

test_that("both repositories hold the same MuData cross-writer harness", {
  python <- file.path(cellucid_repository_root(), "..", "cellucid-python", "tests", "cross_writer_mudata.py")
  skip_if_not(file.exists(python), "cellucid-python is not checked out beside this repository")
  expect_identical(
    unname(tools::md5sum(python)),
    unname(tools::md5sum(testthat::test_path("python", "cross_writer_mudata.py")))
  )
})
