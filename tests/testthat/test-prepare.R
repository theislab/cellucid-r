store_path <- function(name = "store.zarr") {
  file.path(tempfile("cellucid_"), name)
}

prepare_fixture_store <- function(..., name = "store.zarr") {
  args <- prepare_fixture(store_path(name))
  extra <- list(...)
  args[names(extra)] <- extra
  suppressMessages(do.call(cellucid_prepare, args))
}

test_that("the root, obs and var are AnnData elements with their encodings", {
  skip_if_not_installed("Matrix")
  store <- prepare_fixture_store()
  expect_identical(store_json(store, ".zgroup"), list(zarr_format = 2L))
  expect_identical(
    store_json(store, ".zattrs"),
    list(`encoding-type` = "anndata", `encoding-version` = "0.1.0")
  )
  for (group in c("obsm", "obsp", "varm", "varp", "layers", "uns")) {
    expect_identical(
      store_json(store, file.path(group, ".zattrs")),
      list(`encoding-type` = "dict", `encoding-version` = "0.1.0"),
      info = group
    )
  }

  obs <- store_json(store, "obs/.zattrs")
  expect_identical(obs$`encoding-type`, "dataframe")
  expect_identical(obs$`encoding-version`, "0.2.0")
  expect_identical(obs$`_index`, "_index")
  expect_identical(
    unlist(obs$`column-order`),
    c("cluster", "batch", "barcode", "score", "flag", "maybe", "n_genes")
  )
  expect_identical(store_array(store, "obs/_index"), paste0("cell_", 1:6))
  expect_identical(
    store_json(store, "obs/_index/.zattrs")$`encoding-type`,
    "string-array"
  )

  # A factor keeps its levels in their order; a repeated string column is a
  # categorical in code-point order; a unique one stays a string array.
  expect_identical(
    store_json(store, "obs/cluster/.zattrs"),
    list(ordered = FALSE, `encoding-type` = "categorical", `encoding-version` = "0.2.0")
  )
  expect_identical(store_array(store, "obs/cluster/categories"), c("c", "b", "a"))
  expect_identical(store_array(store, "obs/cluster/codes"), c(1L, 2L, 1L, 0L, 2L, 1L))
  expect_identical(store_json(store, "obs/cluster/codes/.zarray")$dtype, "|i1")
  expect_identical(store_array(store, "obs/batch/categories"), c("x", "y"))
  expect_identical(store_array(store, "obs/barcode"), c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC"))

  score <- store_array(store, "obs/score")
  expect_identical(store_json(store, "obs/score/.zarray")$dtype, "<f8")
  expect_identical(score[-2], c(0.5, 1.25, -2, 3, 0))
  expect_true(is.nan(score[[2]]))
  expect_identical(store_array(store, "obs/flag"), c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_identical(store_json(store, "obs/maybe/.zattrs")$`encoding-type`, "nullable-boolean")
  expect_identical(store_array(store, "obs/maybe/mask"), c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_identical(store_array(store, "obs/maybe/values"), c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(store_json(store, "obs/n_genes/.zattrs")$`encoding-type`, "nullable-integer")
  expect_identical(store_array(store, "obs/n_genes/values"), c(10L, 20L, 0L, 40L, 50L, 60L))

  expect_identical(store_array(store, "var/_index"), c("ENSG1", "ENSG2", "ENSG3", "ENSG4"))
  expect_identical(unlist(store_json(store, "var/.zattrs")$`column-order`), "symbol")
  expect_identical(store_array(store, "var/symbol/categories"), c("CD3E", "CD8A", "MS4A1"))
})

test_that("X is gene-major float32 CSC with int32 indices", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  store <- suppressMessages(do.call(cellucid_prepare, args))
  expression <- args$gene_expression
  expect_identical(
    store_json(store, "X/.zattrs"),
    list(shape = list(6L, 4L), `encoding-type` = "csc_matrix", `encoding-version` = "0.1.0")
  )
  data <- store_json(store, "X/data/.zarray")
  expect_identical(data$dtype, "<f4")
  expect_identical(data$compressor, list(id = "zlib", level = 6L))
  expect_identical(data$zarr_format, 2L)
  expect_identical(data$order, "C")
  expect_identical(store_json(store, "X/indices/.zarray")$dtype, "<i4")
  expect_identical(store_json(store, "X/indptr/.zarray")$dtype, "<i4")
  expect_identical(store_array(store, "X/data"), as_float32(expression@x))
  expect_identical(store_array(store, "X/indices"), expression@i)
  expect_identical(store_array(store, "X/indptr"), expression@p)
})

test_that("embeddings, vector fields and the graph are stored as obsm and obsp", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  args$vector_fields <- list(velocity = list(
    X_umap = args$embeddings$X_umap / 2,
    X_pca = args$embeddings$X_pca / 4
  ))
  store <- suppressMessages(do.call(cellucid_prepare, args))
  for (key in c("X_umap", "X_pca")) {
    meta <- store_json(store, file.path("obsm", key, ".zarray"))
    expect_identical(meta$dtype, "<f4", info = key)
    expect_identical(unlist(meta$chunks), unlist(meta$shape), info = key)
    expect_identical(store_array(store, file.path("obsm", key)), args$embeddings[[key]], info = key)
  }
  # Each field is stored under scVelo's name for it, '<field>_<basis>' for
  # the embedding 'X_<basis>', UMAP or not, as wide as that embedding.
  expect_identical(store_array(store, "obsm/velocity_umap"), args$embeddings$X_umap / 2)
  expect_identical(store_array(store, "obsm/velocity_pca"), args$embeddings$X_pca / 4)
  expect_setequal(
    list.dirs(file.path(store, "obsm"), full.names = FALSE, recursive = FALSE),
    c("X_umap", "X_pca", "velocity_umap", "velocity_pca")
  )

  expect_identical(
    store_json(store, "obsp/connectivities/.zattrs"),
    list(shape = list(6L, 6L), `encoding-type` = "csr_matrix", `encoding-version` = "0.1.0")
  )
  graph <- Matrix::sparseMatrix(
    i = store_array(store, "obsp/connectivities/indices") + 1L,
    p = store_array(store, "obsp/connectivities/indptr"),
    x = store_array(store, "obsp/connectivities/data"),
    dims = c(6L, 6L)
  )
  expect_identical(as.matrix(Matrix::t(graph)), as.matrix(args$connectivities))
})

test_that("uns['cellucid'] carries the metadata and only the defaults that were set", {
  skip_if_not_installed("Matrix")
  minimal <- prepare_fixture_store()
  expect_identical(store_array(minimal, "uns/title"), "Fixture")
  expect_identical(
    sort(list.dirs(file.path(minimal, "uns/cellucid"), full.names = FALSE, recursive = FALSE)),
    c("id", "name", "schema")
  )
  expect_identical(store_json(minimal, "uns/cellucid/schema/.zarray")$dtype, "<i8")
  expect_identical(store_array(minimal, "uns/cellucid/schema"), 1)
  expect_identical(
    store_json(minimal, "uns/cellucid/name/.zattrs"),
    list(`encoding-type` = "string", `encoding-version` = "0.2.0")
  )

  full <- prepare_fixture_store(
    dataset_description = "Six cells.",
    citation_text = "Doe et al.",
    citation_url = "https://example.org/doe",
    source_name = "A lab",
    embedding_defaults = list("2d" = "X_umap", "3d" = list(key = "X_pca", columns = c(3, 1, 2))),
    latent_key = "X_pca",
    obs_keys = c("cluster", "score"),
    var_gene_id_column = "symbol",
    normalize_embeddings = FALSE,
    centroid_min_points = 3
  )
  cellucid <- file.path(full, "uns/cellucid")
  expect_identical(store_array(cellucid, "description"), "Six cells.")
  expect_identical(store_array(cellucid, "citation/text"), "Doe et al.")
  expect_identical(store_array(cellucid, "citation/url"), "https://example.org/doe")
  expect_identical(store_array(cellucid, "source/name"), "A lab")
  expect_false(dir.exists(file.path(cellucid, "source/url")))
  defaults <- file.path(cellucid, "defaults")
  expect_identical(
    sort(list.dirs(defaults, full.names = FALSE, recursive = FALSE)),
    sort(c(
      "embedding2d", "embedding3d", "geneIdColumn", "latentKey", "obsKeys",
      "normalizeEmbeddings", "centroidMinPoints"
    ))
  )
  expect_identical(store_array(defaults, "embedding2d"), "X_umap")
  expect_identical(store_array(defaults, "embedding3d/key"), "X_pca")
  expect_identical(store_array(defaults, "embedding3d/columns"), c(2, 0, 1))
  expect_identical(store_array(defaults, "obsKeys"), c("cluster", "score"))
  expect_identical(store_array(defaults, "normalizeEmbeddings"), FALSE)
  expect_identical(store_array(defaults, "centroidMinPoints"), 3)
  expect_identical(store_array(defaults, "geneIdColumn"), "symbol")
})

test_that("a latent key stores the outlier quantiles of every categorical obs column", {
  skip_if_not_installed("Matrix")
  store <- prepare_fixture_store(latent_key = "X_pca", centroid_min_points = 3)
  outliers <- file.path(store, "uns/cellucid/outlierQuantiles")
  expect_identical(
    store_json(store, "uns/cellucid/outlierQuantiles/.zattrs"),
    list(`encoding-type` = "dict", `encoding-version` = "0.1.0")
  )
  expect_identical(store_array(outliers, "latentKey"), "X_pca")
  expect_identical(store_json(outliers, "minPoints/.zarray")$dtype, "<i8")
  expect_identical(store_array(outliers, "minPoints"), 3)
  # cluster (b: 3 cells), batch (text with repeats: 3 and 3) and flag (3 and
  # 3) have a category of 3 cells; maybe (2, 2 and two NA) and barcode (all
  # different) have none; score and n_genes are numbers.
  expect_identical(
    sort(list.dirs(file.path(outliers, "fields"), full.names = FALSE, recursive = FALSE)),
    c("batch", "cluster", "flag")
  )
  meta <- array_meta(outliers, "fields/cluster")
  expect_identical(meta$dtype, "<f4")
  expect_identical(meta$shape, 6L)
  expect_identical(meta$shuffle, 4L)
  # The b cells (1, 3, 6) lie on a line through X_pca, at 7/3, 1/3 and 8/3
  # steps from their centre; a and c have fewer than 3 cells.
  expect_identical(
    store_array(outliers, "fields/cluster"),
    as_float32(c(2 / 3, NaN, 1 / 3, NaN, NaN, 1))
  )

  # Without centroid_min_points a category needs 10 cells: none has.
  plain <- prepare_fixture_store(latent_key = "X_pca", name = "plain.zarr")
  outliers <- file.path(plain, "uns/cellucid/outlierQuantiles")
  expect_identical(store_array(outliers, "minPoints"), 10)
  expect_identical(list.dirs(file.path(outliers, "fields"), full.names = FALSE, recursive = FALSE), character(0))
  expect_false(dir.exists(file.path(prepare_fixture_store(name = "none.zarr"), "uns/cellucid/outlierQuantiles")))
})

test_that("an obs column may have any storable name in the outlier fields", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  args$obs <- data.frame(.uns = factor(rep("a", 6L)), .unsorted = rep("b", 6L), check.names = FALSE)
  args$latent_key <- "X_pca"
  args$centroid_min_points <- 6
  store <- suppressMessages(do.call(cellucid_prepare, args))
  fields <- file.path(store, "uns/cellucid/outlierQuantiles/fields")
  expect_identical(store_array(fields, ".uns"), as_float32(c(6, 4, 2, 2, 4, 6) / 6))
  expect_identical(store_array(fields, ".unsorted"), store_array(fields, ".uns"))
})

test_that("consolidated metadata holds every metadata document of the store", {
  skip_if_not_installed("Matrix")
  store <- prepare_fixture_store()
  consolidated <- store_json(store, ".zmetadata")
  expect_identical(consolidated$zarr_consolidated_format, 1L)
  documents <- list.files(store, all.files = TRUE, recursive = TRUE, pattern = "^\\.z(group|attrs|array)$")
  expect_setequal(names(consolidated$metadata), documents)
  for (document in documents) {
    expect_identical(consolidated$metadata[[document]], store_json(store, document), info = document)
  }
})

test_that("a .zarr.zip is the same store, every entry stored", {
  skip_if_not_installed("Matrix")
  folder <- prepare_fixture_store()
  archive <- prepare_fixture_store(format = "zarr.zip", name = "store.zarr.zip")
  expect_true(file.exists(archive) && !dir.exists(archive))

  listing <- utils::unzip(archive, list = TRUE)
  files <- list.files(folder, all.files = TRUE, recursive = TRUE)
  expect_setequal(listing$Name, files)
  # Every chunk first, then the metadata documents, the consolidated metadata
  # last: the end of the archive, read when it opens, describes the store.
  documents <- basename(listing$Name) %in% c(".zarray", ".zattrs", ".zgroup", ".zmetadata")
  expect_false(any(documents[seq_len(sum(!documents))]))
  expect_identical(listing$Name[[nrow(listing)]], ".zmetadata")

  extracted <- tempfile("cellucid_unzip_")
  utils::unzip(archive, exdir = extracted)
  for (file in files) {
    expect_identical(
      readBin(file.path(extracted, file), "raw", file.size(file.path(extracted, file))),
      readBin(file.path(folder, file), "raw", file.size(file.path(folder, file))),
      info = file
    )
  }

  # Every central-directory record declares method 0 (stored), so a reader
  # can fetch an entry with one range request.
  bytes <- readBin(archive, "raw", file.size(archive))
  signature <- as.raw(c(0x50, 0x4b, 0x01, 0x02))
  starts <- which(bytes == signature[[1L]])
  starts <- starts[vapply(starts, function(at) identical(bytes[at + 0:3], signature), logical(1))]
  expect_length(starts, length(files))
  methods <- vapply(starts, function(at) readBin(bytes[at + 10:11], "integer", size = 2L, endian = "little"), integer(1))
  expect_true(all(methods == 0L))
})

test_that("chunks are padded to the full chunk shape at the edges, in both Zarr formats", {
  for (zarr_format in c(2L, 3L)) {
    root <- tempfile("cellucid_chunks_")
    store <- cellucid:::.store_open(root, zip = FALSE, zarr_format = zarr_format)
    cellucid:::.zarr_vector(store, "v", c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5), "<f4", target = 3)
    cellucid:::.zarr_vector(store, "s", c("été", "", "中文", "z"), "|O", target = 3)
    matrix <- matrix(seq_len(15) / 2, nrow = 5)
    cellucid:::.zarr_matrix(store, "m", matrix, "<f4", 2, 2)
    cellucid:::.store_close(store)

    chunk_files <- function(key) {
      folder <- if (zarr_format == 2L) file.path(root, key) else file.path(root, key, "c")
      files <- list.files(folder, recursive = TRUE)
      sort(files[!startsWith(basename(files), ".") & basename(files) != "zarr.json"])
    }
    v <- array_meta(root, "v")
    s <- array_meta(root, "s")
    m <- array_meta(root, "m")
    expect_identical(v$chunks, 3L)
    expect_identical(v$shuffle, 4L)
    expect_identical(chunk_files("v"), c("0", "1", "2"))
    expect_identical(read_chunk(root, "v", v, 2L, 3L), c(7.5, 0, 0))
    expect_identical(store_array(root, "v"), c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5))
    expect_identical(read_chunk(root, "s", s, 1L, 3L), c("z", "", ""))
    expect_identical(store_array(root, "s"), c("été", "", "中文", "z"))
    expect_identical(s$shuffle, 0L)
    expect_identical(
      chunk_files("m"),
      if (zarr_format == 2L) c("0.0", "0.1", "1.0", "1.1", "2.0", "2.1") else c("0/0", "0/1", "1/0", "1/1", "2/0", "2/1")
    )
    expect_identical(read_chunk(root, "m", m, c(2L, 1L), 4L), c(7.5, 0, 0, 0))
    expect_identical(store_array(root, "m"), matrix)
  }
  expect_identical(store_json(root, "s/zarr.json")$codecs, list(
    list(name = "vlen-utf8", configuration = structure(list(), names = character())),
    list(name = "gzip", configuration = list(level = 6L))
  ))
})

test_that("dense expression is stored gene-major in tall strips", {
  args <- list(
    out_path = store_path(),
    format = "zarr",
    embeddings = list(X_umap = matrix(seq_len(6), ncol = 2)),
    dataset_id = "dense",
    dataset_name = "Dense",
    gene_expression = matrix(c(0, 1, 2, 3, 4, 5), nrow = 3, dimnames = list(NULL, c("g1", "g2")))
  )
  store <- do.call(cellucid_prepare, args)
  meta <- store_json(store, "X/.zarray")
  expect_identical(meta$dtype, "<f4")
  expect_identical(unlist(meta$shape), c(3L, 2L))
  expect_identical(store_array(store, "X"), unname(args$gene_expression))
  expect_identical(store_array(store, "obs/_index"), c("0", "1", "2"))
  expect_identical(store_array(store, "var/_index"), c("g1", "g2"))
})

test_that("a stored zero of sparse expression is dropped, and the matrix stays sparse", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  expression <- args$gene_expression
  expression@x[c(2L, 5L)] <- 0
  args$gene_expression <- expression
  store <- suppressMessages(do.call(cellucid_prepare, args))
  kept <- Matrix::drop0(expression)
  expect_identical(store_json(store, "X/.zattrs")$`encoding-type`, "csc_matrix")
  expect_length(store_array(store, "X/data"), length(expression@x) - 2L)
  expect_identical(store_array(store, "X/data"), as_float32(kept@x))
  expect_identical(store_array(store, "X/indices"), kept@i)
  expect_identical(store_array(store, "X/indptr"), kept@p)
})

test_that("integer counts above 2^24 keep the whole X in float64, and nothing else does", {
  skip_if_not_installed("Matrix")
  sparse_dtype <- function(values, ...) {
    args <- prepare_fixture(store_path())
    args$gene_expression@x <- values
    extra <- list(...)
    args[names(extra)] <- extra
    store <- suppressMessages(do.call(cellucid_prepare, args))
    list(dtype = store_json(store, "X/data/.zarray")$dtype, data = store_array(store, "X/data"))
  }
  counts <- c(1, 2, 3, 2^24 + 1, 5, 6, 7, 8, -2^30)
  above <- sparse_dtype(counts)
  expect_identical(above$dtype, "<f8")
  expect_identical(above$data, counts)
  # 2^24 itself is a float32 integer; a fraction anywhere means the matrix
  # is not counts, and float32 is its type.
  expect_identical(sparse_dtype(c(1, 2, 3, 2^24, 5, 6, 7, 8, -2^24))$dtype, "<f4")
  expect_identical(sparse_dtype(c(1, 2, 3, 2^24 + 1, 5, 6, 7, 8, 0.5))$dtype, "<f4")
  # The type is decided on the genes stored, after gene_identifiers chose
  # them, as the Python writer decides it: ENSG1 holds the large count.
  expect_identical(sparse_dtype(counts, gene_identifiers = c("ENSG2", "ENSG3"))$dtype, "<f4")
  expect_identical(sparse_dtype(counts, gene_identifiers = c("ENSG1", "ENSG3"))$dtype, "<f8")

  dense <- function(values) {
    store <- cellucid_prepare(
      store_path(),
      format = "zarr",
      embeddings = list(X_umap = matrix(seq_len(4), ncol = 2)),
      dataset_id = "dense",
      dataset_name = "Dense",
      gene_expression = matrix(values, nrow = 2, dimnames = list(NULL, c("g1", "g2")))
    )
    list(dtype = store_json(store, "X/.zarray")$dtype, values = store_array(store, "X"))
  }
  integers <- dense(c(0L, 16777217L, 3L, 0L))
  expect_identical(integers$dtype, "<f8")
  expect_identical(integers$values, matrix(c(0, 16777217, 3, 0), nrow = 2))
  expect_identical(dense(c(0, 2^25, 3, 0.25))$dtype, "<f4")
  expect_identical(dense(c(TRUE, FALSE, TRUE, TRUE))$dtype, "<f4")
})

test_that("category codes take the type pandas gives them", {
  codes_dtype <- function(n_levels) {
    levels <- sprintf("t%05d", seq_len(n_levels))
    store <- cellucid_prepare(
      store_path(),
      format = "zarr",
      embeddings = list(X_umap = matrix(seq_len(4), ncol = 2)),
      dataset_id = "codes",
      dataset_name = "Codes",
      obs = data.frame(type = factor(levels[c(1L, n_levels)], levels = levels))
    )
    store_json(store, "obs/type/codes/.zarray")$dtype
  }
  # pandas compares the count with each type's maximum strictly.
  expect_identical(codes_dtype(126L), "|i1")
  expect_identical(codes_dtype(127L), "<i2")
  expect_identical(codes_dtype(32766L), "<i2")
  expect_identical(codes_dtype(32767L), "<i4")
  expect_identical(cellucid:::.ad_codes_dtype(2147483646), "<i4")
  expect_identical(cellucid:::.ad_codes_dtype(2147483647), "<i8")
})

test_that("expression is chunked for gene reads and the graph for whole reads", {
  skip_if_not_installed("Matrix")
  n <- 70000L
  ring <- c(seq(2L, n), 1L)
  graph <- Matrix::sparseMatrix(
    i = c(seq_len(n), ring), j = c(ring, seq_len(n)), x = 1, dims = c(n, n)
  )
  expression <- Matrix::sparseMatrix(i = seq_len(n), j = rep(1:2, length.out = n), x = 1, dims = c(n, 2L))
  store <- cellucid_prepare(
    store_path(),
    format = "zarr",
    embeddings = list(X_umap = matrix(as.double(seq_len(2L * n)), ncol = 2)),
    dataset_id = "chunks",
    dataset_name = "Chunks",
    gene_expression = expression,
    var = data.frame(row.names = c("g1", "g2")),
    connectivities = graph
  )
  chunks <- function(key) unlist(store_json(store, file.path(key, ".zarray"))$chunks)
  expect_identical(chunks("X/data"), 65536L)
  expect_identical(chunks("X/indices"), 65536L)
  expect_identical(chunks("X/indptr"), 3L)
  expect_identical(chunks("obsp/connectivities/data"), 2L * n)
  expect_identical(chunks("obsp/connectivities/indices"), 2L * n)
  expect_identical(chunks("obsp/connectivities/indptr"), n + 1L)
  expect_identical(chunks("obsm/X_umap"), c(n, 2L))
  expect_identical(chunks("obs/_index"), n)
})

test_that("repeated cell and gene names are made unique the anndata way", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  colnames(args$gene_expression) <- c("A", "B", "A", "C")
  rownames(args$var) <- NULL
  rownames(args$gene_expression)[[3L]] <- "cell_1"
  args$obs <- NULL
  expect_message(
    expect_message(store <- do.call(cellucid_prepare, args), "Cell names: 1 duplicate"),
    "Gene names: 1 duplicate name was made unique, e.g. 'A' → 'A-1'."
  )
  expect_identical(store_array(store, "var/_index"), c("A", "B", "A-1", "C"))
  expect_identical(store_array(store, "obs/_index")[[3L]], "cell_1-1")
})

test_that("gene_identifiers narrows X and var to the named genes", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  args$gene_identifiers <- c("ENSG4", "ENSG2")
  store <- suppressMessages(do.call(cellucid_prepare, args))
  expect_identical(store_array(store, "var/_index"), c("ENSG4", "ENSG2"))
  expect_identical(
    store_array(store, "X/indptr"),
    args$gene_expression[, c(4, 2)]@p
  )

  by_symbol <- prepare_fixture(store_path())
  by_symbol$var_gene_id_column <- "symbol"
  by_symbol$gene_identifiers <- c("CD3E-1", "MS4A1")
  store <- suppressMessages(do.call(cellucid_prepare, by_symbol))
  expect_identical(store_array(store, "var/_index"), c("ENSG4", "ENSG3"))

  by_symbol$out_path <- store_path()
  by_symbol$gene_identifiers <- c("CD3E", "NOPE")
  expect_error(
    do.call(cellucid_prepare, by_symbol),
    "gene_identifiers contains identifiers not found in var: c\\(\"NOPE\"\\)"
  )
})

test_that("a required argument that was never supplied is reported as missing", {
  for (argument in c("out_path", "format", "embeddings", "dataset_id", "dataset_name")) {
    expect_true(identical(formals(cellucid_prepare)[[argument]], quote(expr = )), info = argument)
  }
  expect_error(
    cellucid_prepare(),
    paste0(
      "^cellucid_prepare\\(\\) is missing 5 required arguments:\\n",
      "  - out_path: .*\\n  - format: .*\\n  - embeddings: .*\\n",
      "  - dataset_id: .*\\n  - dataset_name: .*$"
    )
  )
  expect_error(
    cellucid_prepare(store_path(), "zarr", list(X_umap = diag(2)), dataset_id = "a"),
    "^cellucid_prepare\\(\\) is missing 1 required argument:\\n  - dataset_name: one non-empty name"
  )
  # format has no default: the message names both choices and what each is for.
  expect_error(
    cellucid_prepare(store_path(), embeddings = list(X_umap = diag(2)), dataset_id = "a", dataset_name = "A"),
    paste0(
      "^cellucid_prepare\\(\\) is missing 1 required argument:\\n",
      "  - format: the kind of output to write, \"zarr\" \\(a \\.zarr folder of small ",
      "files: for GitHub.*\"zarr\\.zip\" \\(the same store in one ",
      "\\.zarr\\.zip file: for single-file hosts"
    )
  )
})

test_that("format chooses the output, and out_path must agree with it", {
  args <- list(
    embeddings = list(X_umap = diag(2)), dataset_id = "fmt", dataset_name = "Format"
  )
  write <- function(out_path, format) {
    do.call(cellucid_prepare, c(list(out_path = out_path, format = format), args))
  }
  folder <- write(store_path("a.zarr"), "zarr")
  expect_true(dir.exists(folder))
  archive <- write(store_path("a.zarr.zip"), "zarr.zip")
  expect_true(file.exists(archive) && !dir.exists(archive))
  # Matched as match.arg() matches: an unambiguous prefix names its choice.
  expect_true(file.exists(write(store_path("b.zarr.zip"), "zarr.")))

  root <- tempfile("cellucid_format_")
  choices <- paste0(
    "^format must be \"zarr\" or \"zarr\\.zip\":\\n",
    "  \"zarr\"      a \\.zarr folder of small files: for GitHub \\(no file over ",
    "100 MB, no LFS\\)\\.\\n",
    "  \"zarr\\.zip\"  the same store in one \\.zarr\\.zip file: for single-file .*\\n",
    "\\(\\.h5ad is written by Python's cellucid\\.prepare\\(format = \"h5ad\"\\); ",
    "this package writes no HDF5\\.\\)\\n"
  )
  expect_error(write(file.path(root, "a.zarr"), "h5ad"), paste0(choices, "Got \"h5ad\"\\.$"))
  expect_error(write(file.path(root, "a.zarr"), "z"), paste0(choices, "Got \"z\"\\.$"))
  expect_error(write(file.path(root, "a.zarr"), NA), paste0(choices, "Got NA\\.$"))
  expect_error(write(file.path(root, "a.zarr"), c("zarr", "zarr.zip")), paste0(choices, "Got c\\("))

  expect_error(
    write(file.path(root, "a.zarr.zip"), "zarr"),
    paste0(
      "format = \"zarr\" writes a .zarr folder, and out_path '.*a.zarr.zip' ends in ",
      "\".zarr.zip\", which format = \"zarr.zip\" writes: write '.*a.zarr' instead."
    ),
    fixed = FALSE
  )
  expect_error(
    write(file.path(root, "a.zarr"), "zarr.zip"),
    "format = \"zarr.zip\" writes one .zarr.zip file, and out_path '.*a.zarr' ends in \".zarr\", which format = \"zarr\" writes: write '.*a.zarr.zip' instead."
  )
  expect_error(
    write(file.path(root, "a.h5ad"), "zarr"),
    "ends in \".h5ad\", which this package does not write \\(Python's cellucid.prepare\\(format = \"h5ad\"\\) does\\): write '.*a.zarr' instead."
  )
  expect_error(
    write("dataset", "zarr.zip"),
    "out_path 'dataset' does not end in \".zarr.zip\": write 'dataset.zarr.zip' instead.",
    fixed = TRUE
  )
  expect_false(dir.exists(root))

  # A folder where one file is to be written is never replaced.
  taken <- store_path("taken.zarr.zip")
  dir.create(taken, recursive = TRUE)
  expect_error(
    do.call(cellucid_prepare, c(list(out_path = taken, format = "zarr.zip", force = TRUE), args)),
    "is a folder, and format = \"zarr.zip\" writes a file"
  )
  expect_true(dir.exists(taken))
})

test_that("each mistake is refused with its reason before anything is written", {
  skip_if_not_installed("Matrix")
  refused <- function(pattern, ..., name = "store.zarr") {
    args <- prepare_fixture(store_path(name))
    changes <- list(...)
    for (key in names(changes)) args[key] <- list(changes[[key]])
    expect_error(suppressMessages(do.call(cellucid_prepare, args)), pattern)
    expect_false(file.exists(args$out_path))
    expect_identical(list.files(dirname(args$out_path), all.files = TRUE, no.. = TRUE), character())
  }
  refused("out_path '.*store' does not end in \".zarr\": write '.*store.zarr' instead", name = "store")
  refused("ends in \".zarr\", which format = \"zarr\" writes", format = "zarr.zip")
  refused("^format must be \"zarr\" or \"zarr\\.zip\"", format = "folder")
  refused("dataset_id 'a b' is not a portable identifier", dataset_id = "a b")
  refused("^dataset_name is displayed verbatim", dataset_name = "Fixture ")
  refused("^citation_text is required whenever citation_url is supplied\\.$", citation_url = "https://x.org")
  refused("source_url must be an absolute http:// or https:// address", source_name = "Lab", source_url = "example.org")
  refused("embeddings must contain at least one embedding", embeddings = list())
  refused("Embedding 'X_umap' must contain only finite values", embeddings = list(X_umap = matrix(c(1, NA), ncol = 2)))
  refused("All embeddings must have one row per cell", embeddings = list(X_umap = diag(6), X_b = diag(5)))
  refused("Embedding 'X/umap' contains a character that cannot appear in a file name", embeddings = list(`X/umap` = diag(6)))
  refused("Embedding names 'X_umap' and 'x_UMAP' differ only in letter case", embeddings = list(X_umap = diag(6), x_UMAP = diag(6)))
  refused("obs column '_index' is reserved", obs = data.frame(`_index` = 1:6, check.names = FALSE))
  refused("obs has 5 rows, but the embeddings have 6 cells", obs = data.frame(a = 1:5))
  refused("Every input must list the same cells in the same order", obs = data.frame(a = 1:6, row.names = paste0("c", 1:6)))
  refused("gene_expression has 5 rows", gene_expression = Matrix::Matrix(0, 5, 4, sparse = TRUE))
  refused("rownames\\(var\\) and colnames\\(gene_expression\\) disagree at gene 1", var = data.frame(a = 1:4, row.names = paste0("g", 1:4)))
  refused("var describes the genes of gene_expression", gene_expression = NULL)
  refused("Connectivity matrix must be exactly symmetric", connectivities = Matrix::sparseMatrix(i = 1, j = 2, x = 1, dims = c(6L, 6L)))
  refused("Connectivity diagonal values must all be exactly zero", connectivities = Matrix::Diagonal(6))
  refused("connectivities must be 6 x 6", connectivities = Matrix::Diagonal(5))
  refused("has 2 columns, and 'X_pca' has 3: a vector field has one component per embedding column", vector_fields = list(velocity = list(X_pca = matrix(0, 6, 2))))
  refused(
    "the embedding needs an 'X_<basis>' key",
    embeddings = list(umap = diag(6)[, 1:2]),
    vector_fields = list(velocity = list(umap = matrix(0, 6, 2)))
  )
  refused(
    "obsm\\[\\['X_velocity_umap'\\]\\], which the viewer reads as an embedding",
    vector_fields = list(X_velocity = list(X_umap = matrix(0, 6, 2)))
  )
  refused(
    "reads as the field 'T' on 'X_fwd_umap'",
    embeddings = list(X_umap = diag(6)[, 1:2], X_fwd_umap = diag(6)[, 1:2]),
    vector_fields = list(T_fwd = list(X_umap = matrix(0, 6, 2)))
  )
  refused(
    "would both be stored as obsm\\[\\['a_b_umap'\\]\\]",
    embeddings = list(X_umap = diag(6)[, 1:2], X_b_umap = diag(6)[, 1:2]),
    vector_fields = list(a_b = list(X_umap = matrix(0, 6, 2)), a = list(X_b_umap = matrix(0, 6, 2)))
  )
  refused("is not one of the embeddings", vector_fields = list(velocity = list(X_tsne = matrix(0, 6, 2))))
  refused("embedding_defaults\\[\\[\"3d\"\\]\\] names 'X_umap', which has 2 columns", embedding_defaults = c("3d" = "X_umap"))
  refused("embedding_defaults\\[\\[\"2d\"\\]\\]\\$columns must name 2 different", embedding_defaults = list("2d" = list(key = "X_pca", columns = c(1, 1))))
  refused("latent_key must name one of the embeddings", latent_key = "X_scvi")
  refused("obs_keys contain columns not in obs: c\\(\"nope\"\\)", obs_keys = c("cluster", "nope"))
  refused("vector_field_default names a field, and no vector_fields were supplied", vector_field_default = "velocity")
  refused(
    "vector_field_default must name one of vector_fields: c\\(\"velocity\"\\)",
    vector_fields = list(velocity = list(X_umap = matrix(0, 6, 2))),
    vector_field_default = "drift"
  )
  # The latent space is not shown when the file opens, so a field drawn only
  # on it cannot be the one the overlay opens on.
  refused(
    "vector_field_default names 'drift', and the viewer shows none of its embeddings",
    vector_fields = list(drift = list(X_pca = matrix(0, 6, 3))),
    latent_key = "X_pca",
    vector_field_default = "drift"
  )
  refused("var_gene_id_column 'nope' not found in var", var_gene_id_column = "nope")
})

test_that("a column with no AnnData encoding fails the write and leaves nothing", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  args$obs$day <- as.Date("2026-01-01") + 0:5
  expect_error(
    suppressMessages(do.call(cellucid_prepare, args)),
    "obs column 'day' is a Date column, which has no AnnData encoding"
  )
  expect_identical(list.files(dirname(args$out_path), all.files = TRUE, no.. = TRUE), character())

  args$obs$day <- c(1, Inf, 2, 3, 4, 5)
  expect_error(
    suppressMessages(do.call(cellucid_prepare, args)),
    "obs column 'day' holds 1 infinite value \\(first at row 2\\)"
  )
})

test_that("force replaces a store only once the new one is complete", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  store <- suppressMessages(do.call(cellucid_prepare, args))
  expect_error(
    do.call(cellucid_prepare, args),
    "already exists; pass force = TRUE to replace it"
  )

  # A failed replacement leaves the existing store as it was.
  before <- tools::md5sum(list.files(store, recursive = TRUE, all.files = TRUE, full.names = TRUE))
  failing <- args
  failing$force <- TRUE
  failing$obs$day <- as.Date("2026-01-01") + 0:5
  expect_error(suppressMessages(do.call(cellucid_prepare, failing)), "Date column")
  expect_identical(
    tools::md5sum(list.files(store, recursive = TRUE, all.files = TRUE, full.names = TRUE)),
    before
  )
  expect_identical(list.files(dirname(store), all.files = TRUE, no.. = TRUE), "store.zarr")

  args$force <- TRUE
  args$dataset_name <- "Replaced"
  suppressMessages(do.call(cellucid_prepare, args))
  expect_identical(store_array(store, "uns/title"), "Replaced")
  expect_identical(list.files(dirname(store), all.files = TRUE, no.. = TRUE), "store.zarr")

  # force never deletes a folder that is not a Zarr store.
  other <- store_path("notes.zarr")
  dir.create(other, recursive = TRUE)
  writeLines("keep me", file.path(other, "notes.txt"))
  args$out_path <- other
  expect_error(do.call(cellucid_prepare, args), "exists and is not a Zarr store")
  expect_true(file.exists(file.path(other, "notes.txt")))
})

test_that("the gene-name column is stored as it is, with a note on what the viewer changes", {
  skip_if_not_installed("Matrix")
  args <- prepare_fixture(store_path())
  args$var_gene_id_column <- "symbol"
  expect_message(
    store <- do.call(cellucid_prepare, args),
    "Gene names from var column 'symbol': 1 duplicate name was made unique, e.g. 'CD3E' → 'CD3E-1'.",
    fixed = TRUE
  )
  expect_identical(store_array(store, "var/symbol/categories"), c("CD3E", "CD8A", "MS4A1"))
  expect_identical(store_array(store, "var/symbol/codes"), c(0L, 1L, 2L, 0L))
})

test_that("a default vector field is one drawn on an embedding shown when the file opens", {
  skip_if_not_installed("Matrix")
  # X_pca is the latent space, so it is shown only because a slot names it.
  store <- prepare_fixture_store(
    embedding_defaults = c("3d" = "X_pca"),
    latent_key = "X_pca",
    vector_fields = list(drift = list(X_pca = matrix(0, 6, 3))),
    vector_field_default = "drift"
  )
  expect_identical(store_array(file.path(store, "uns/cellucid/defaults"), "vectorFieldDefault"), "drift")
  expect_identical(
    cellucid:::.best_guess_order(c("X_pca", "X_spatial", "X_other", "X_tsne")),
    c("X_tsne", "X_spatial", "X_other", "X_pca")
  )
  expect_identical(
    cellucid:::.shown_embedding_keys(
      list(X_pca = matrix(0, 2, 50), X_draw_graph_fa = matrix(0, 2, 8)),
      list(),
      NULL
    ),
    "X_draw_graph_fa"
  )
})
