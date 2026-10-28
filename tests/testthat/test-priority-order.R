# The priority order (R/priority.R; docs: "Cellucid-ready AnnData
# specification"): cellucid_prepare() stores the cells by ascending uniform
# rank under the file's seed, records each row's original position and, beside
# a sparse X, each gene's values below each prefix level; and
# cellucid_restore_original_order() puts the cells back.

store_path <- function() file.path(tempfile("cellucid_priority_"), "store.zarr")

priority_case <- function(out_path, ...) {
  args <- prepare_fixture(out_path)
  args$priority_seed <- NULL
  utils::modifyList(args, list(...))
}

test_that("the cells are stored by ascending rank, every per-cell element with them", {
  for (seed in c(0, 7, 4294967295)) {
    args <- priority_case(store_path(), priority_seed = seed)
    store <- suppressMessages(do.call(cellucid_prepare, args))
    original <- store_array(store, "obs/cellucid_original_position")
    expect_identical(array_meta(store, "obs/cellucid_original_position")$dtype, "<u4")
    expect_setequal(original, 0:5)
    ranks <- cellucid:::.subsample_rank(seed, original)
    expect_true(all(diff(ranks) > 0))
    rows <- original + 1
    expect_identical(store_array(store, "obs/_index"), rownames(args$gene_expression)[rows])
    score <- args$obs$score[rows]
    score[is.na(score)] <- NaN
    expect_identical(store_array(store, "obs/score"), score)
    expect_identical(store_array(store, "obsm/X_pca"), as_float32_matrix(args$embeddings$X_pca[rows, ]))
    # obsp: both axes permuted.
    graph <- as.matrix(args$connectivities)[rows, rows]
    pointers <- store_array(store, "obsp/connectivities/indptr")
    indices <- store_array(store, "obsp/connectivities/indices")
    rebuilt <- matrix(0, 6, 6)
    for (r in seq_len(6)) {
      for (e in seq_len(pointers[[r + 1L]] - pointers[[r]])) {
        rebuilt[r, indices[[pointers[[r]] + e]] + 1] <- store_array(store, "obsp/connectivities/data")[[pointers[[r]] + e]]
      }
    }
    expect_identical(rebuilt, unname(graph))
    # X: each gene's rows ascending, the values of the cells now at them.
    x_pointers <- store_array(store, "X/indptr")
    x_rows <- store_array(store, "X/indices")
    x_values <- store_array(store, "X/data")
    dense <- as.matrix(args$gene_expression)[rows, ]
    for (gene in seq_len(4)) {
      range <- seq_len(x_pointers[[gene + 1L]] - x_pointers[[gene]]) + x_pointers[[gene]]
      expect_true(all(diff(x_rows[range]) > 0))
      expect_identical(x_values[range], as_float32(dense[x_rows[range] + 1, gene]))
    }
    # uns['cellucid']['priority'] and the prefix counts.
    expect_identical(store_array(store, "uns/cellucid/priority/seed"), seed)
    expect_identical(store_array(store, "uns/cellucid/priority/specVersion"), 1)
    levels <- store_array(store, "uns/cellucid/priority/levels")
    expect_identical(levels, c(1, 2, 3, 4))
    counts <- store_array(store, "varm/cellucid_prefix_counts")
    expect_identical(dim(counts), c(4L, 4L))
    for (gene in seq_len(4)) {
      range <- seq_len(x_pointers[[gene + 1L]] - x_pointers[[gene]]) + x_pointers[[gene]]
      expect_identical(counts[gene, ], vapply(levels, function(level) sum(x_rows[range] < level), numeric(1)))
    }
  }
})

test_that("the prefix levels step down from the cell count by 3/4, as cellucid-python's do", {
  expect_identical(cellucid:::.prefix_levels(1), numeric(0))
  expect_identical(cellucid:::.prefix_levels(4), c(1, 2, 3))
  expect_identical(cellucid:::.prefix_levels(10000)[c(1, 2, 29)], c(1, 2, 7500))
  big <- cellucid:::.prefix_levels(4294967295)
  expect_identical(big[[length(big)]], floor(4294967295 * 3 / 4))
  expect_true(all(diff(big) > 0))
})

test_that("a genes x cells assay counts its prefixes by cell, the genes kept in their order", {
  matrix <- Matrix::rsparsematrix(30, 50, 0.2)
  matrix@x <- abs(matrix@x) + 1
  dimnames(matrix) <- list(paste0("g", 1:30), paste0("c", 1:50))
  expression <- cellucid:::.expression_input(structure(list(matrix = matrix), class = "cellucid_feature_major"), "x")
  expression$genes <- c(5L, 2L, 30L)
  levels <- cellucid:::.prefix_levels(50)
  counts <- cellucid:::.expression_prefix_counts(expression, levels)
  for (g in seq_along(expression$genes)) {
    cells <- which(matrix[expression$genes[[g]], ] != 0) - 1
    expect_identical(counts[g, ], vapply(levels, function(level) as.double(sum(cells < level)), numeric(1)))
  }
  by_cell <- cellucid:::.expression_input(Matrix::t(matrix), "x")
  counts <- cellucid:::.expression_prefix_counts(by_cell, levels)
  expect_identical(dim(counts), c(30L, length(levels)))
})

test_that("priority_seed = NULL stores the cells as given; a seed outside uint32 is refused", {
  args <- priority_case(store_path())
  args["priority_seed"] <- list(NULL)
  store <- suppressMessages(do.call(cellucid_prepare, args))
  expect_identical(store_array(store, "obs/_index"), rownames(args$gene_expression))
  expect_false(file.exists(file.path(store, "varm", "cellucid_prefix_counts")))
  expect_false("priority" %in% list.files(file.path(store, "uns", "cellucid")))
  for (seed in list(-1, 2^32, 1.5, "0")) {
    expect_error(
      suppressMessages(do.call(cellucid_prepare, priority_case(store_path(), priority_seed = seed))),
      "priority_seed must be one whole number"
    )
  }
})

# The record cellucid_prepare() wrote to `store`, as an importer that keeps
# uns gives it (anndataR, zellkonverter: metadata(x)$cellucid$priority).
stored_record <- function(store) {
  list(
    seed = store_array(store, "uns/cellucid/priority/seed"),
    specVersion = store_array(store, "uns/cellucid/priority/specVersion"),
    nObs = store_array(store, "uns/cellucid/priority/nObs"),
    levels = store_array(store, "uns/cellucid/priority/levels")
  )
}

test_that("the record names its cell count and its levels", {
  store <- suppressMessages(do.call(cellucid_prepare, priority_case(store_path(), priority_seed = 7)))
  record <- stored_record(store)
  expect_equal(as.double(record$nObs), 6)
  expect_equal(as.double(record$levels), cellucid:::.prefix_levels(6))
  verified <- cellucid:::.stored_order(record, 6, store_array(store, "obs/cellucid_original_position"))
  expect_null(verified$notice)
  expect_equal(verified$positions, as.double(store_array(store, "obs/cellucid_original_position")))
  # The column is a check, not the source: a record without it verifies.
  expect_equal(cellucid:::.stored_order(record, 6, NULL)$positions, verified$positions)
})

test_that("a record that does not describe the cells is taken as cells Cellucid did not write", {
  store <- suppressMessages(do.call(cellucid_prepare, priority_case(store_path(), priority_seed = 7)))
  record <- stored_record(store)
  column <- store_array(store, "obs/cellucid_original_position")
  swapped <- column[c(1, 2, 4, 3, 5, 6)]
  cases <- list(
    list(record, 5, column[-6], "written for 6 cells, and there are 5"),
    list(record, 6, swapped, "is not its seed's order"),
    list(utils::modifyList(record, list(specVersion = 2)), 6, column, "version 2 of the subsampling specification"),
    list(record[setdiff(names(record), "nObs")], 6, column, "records no cell count"),
    list(utils::modifyList(record, list(levels = cellucid:::.prefix_levels(12))), 6, column, "levels are not those of 6 cells"),
    list(utils::modifyList(record, list(seed = -1)), 6, column, "seed is not a whole number"),
    list(NULL, 6, column, "has no uns\\['cellucid'\\]\\['priority'\\] beside it")
  )
  for (case in cases) {
    order <- cellucid:::.stored_order(case[[1L]], case[[2L]], case[[3L]])
    expect_null(order$positions)
    expect_match(order$notice, case[[4L]])
  }
  expect_identical(cellucid:::.stored_order(NULL, 6, NULL), list(positions = NULL, notice = NULL))
})

test_that("a file prepared again from an object in its stored order keeps its cells' original positions", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "s", dataset_name = "S", priority_seed = 3))
  original <- store_array(store, "obs/cellucid_original_position")
  # The stored cells, in their stored order, as an importer that keeps uns reads them.
  stored <- object[, original + 1]
  SummarizedExperiment::colData(stored)$cellucid_original_position <- original
  S4Vectors::metadata(stored)$cellucid <- list(priority = stored_record(store))
  second <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = stored, dataset_id = "s", dataset_name = "S", priority_seed = 3))
  expect_identical(store_array(second, "obs/cellucid_original_position"), original)
  expect_identical(store_array(second, "obs/_index"), store_array(store, "obs/_index"))
  # Without a seed: written in the original order, no stale column or record.
  given <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = stored, dataset_id = "s", dataset_name = "S", priority_seed = NULL))
  expect_identical(store_array(given, "obs/_index"), colnames(object))
  expect_false(file.exists(file.path(given, "obs", "cellucid_original_position")))
  expect_false("priority" %in% list.files(file.path(given, "uns", "cellucid")))
  # A column that is not the record's order: the cells are taken as given,
  # the column not written, and both warnings say why.
  broken <- stored
  SummarizedExperiment::colData(broken)$cellucid_original_position <- c(0, 0, seq_len(ncol(broken) - 2L))
  warnings <- character()
  rewritten <- withCallingHandlers(
    suppressMessages(cellucid_prepare(store_path(), "zarr", data = broken, dataset_id = "s", dataset_name = "S", priority_seed = NULL)),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  expect_true(any(grepl("is not its seed's order", warnings)))
  expect_true(any(grepl("is not written", warnings)))
  expect_identical(store_array(rewritten, "obs/_index"), colnames(stored))
  expect_false(file.exists(file.path(rewritten, "obs", "cellucid_original_position")))
})

# A priority record of `n` cells under `seed`, and the order its rows hold
# (row r: original position positions[r]).
priority_of <- function(n, seed = 7) {
  list(
    record = list(seed = seed, specVersion = 1L, nObs = n, levels = cellucid:::.prefix_levels(n)),
    positions = order(cellucid:::.subsample_rank(seed, seq_len(n) - 1), method = "radix") - 1
  )
}

test_that("cellucid_restore_original_order() puts a file's cells back in the order given", {
  priority <- priority_of(3)
  obs <- data.frame(cell_type = c("B", "T", "NK"), row.names = c("c1", "c2", "c3"))[priority$positions + 1, , drop = FALSE]
  obs$cellucid_original_position <- priority$positions
  attr(obs, "cellucid") <- list(priority = priority$record)
  restored <- cellucid_restore_original_order(obs)
  expect_identical(rownames(restored), c("c1", "c2", "c3"))
  expect_identical(restored$cell_type, c("B", "T", "NK"))
  expect_false("cellucid_original_position" %in% names(restored))
  expect_null(attr(restored, "cellucid"))
  plain <- data.frame(a = 1:3)
  expect_identical(cellucid_restore_original_order(plain), plain)
  # A column without a record names nothing; a record that no longer
  # describes the cells (one removed) leaves them as they are. Both warn.
  column_only <- obs
  attr(column_only, "cellucid") <- NULL
  expect_warning(expect_identical(cellucid_restore_original_order(column_only), column_only), "names no order")
  dropped <- obs[-2, , drop = FALSE]
  attr(dropped, "cellucid") <- list(priority = priority$record)
  expect_warning(expect_identical(cellucid_restore_original_order(dropped), dropped), "written for 3 cells, and there are 2")
  expect_error(cellucid_restore_original_order(1:3), "x must be a data.frame")
  skip_without_sce()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1:6, nrow = 2, dimnames = list(c("g1", "g2"), c("c1", "c2", "c3"))[c(1, 2)])),
    reducedDims = list(UMAP = matrix(as.double(1:6), ncol = 2))
  )
  stored <- sce[, priority$positions + 1]
  SummarizedExperiment::colData(stored)$cellucid_original_position <- priority$positions
  S4Vectors::metadata(stored)$cellucid <- list(priority = priority$record, id = "kept")
  back <- cellucid_restore_original_order(stored)
  expect_identical(colnames(back), c("c1", "c2", "c3"))
  expect_identical(SingleCellExperiment::reducedDim(back, "UMAP"), SingleCellExperiment::reducedDim(sce, "UMAP"))
  expect_false("cellucid_original_position" %in% names(SummarizedExperiment::colData(back)))
  expect_identical(S4Vectors::metadata(back)$cellucid, list(id = "kept"))
})

# Every matrix of a Seurat assay, by name: an Assay5's layers (SeuratObject
# 5), a version 3 Assay's counts, data and scale.data (SeuratObject 4 and 5;
# SeuratObject 4 has no layers).
assay_matrices <- function(assay) {
  if (methods::is(assay, "Assay5")) {
    names <- SeuratObject::Layers(assay)
    return(stats::setNames(lapply(names, function(layer) as.matrix(SeuratObject::LayerData(assay, layer))), names))
  }
  names <- c("counts", "data", "scale.data")
  stats::setNames(lapply(names, function(slot) as.matrix(methods::slot(assay, slot))), names)
}

test_that("cellucid_restore_original_order() reorders every part of a Seurat object", {
  skip_without_seurat()
  object <- pbmc_multimodal()
  # A version 3 assay beside the others: with SeuratObject 5 they are
  # version 5 assays, so both kinds are reordered.
  adt <- assay_matrices(object[["ADT"]])$counts
  object[["V3"]] <- SeuratObject::CreateAssayObject(counts = Matrix::Matrix(adt, sparse = TRUE))
  n <- ncol(object)
  priority <- priority_of(n, seed = 11)
  stored <- cellucid:::.seurat_in_cell_order(object, priority$positions + 1)
  expect_identical(colnames(stored), colnames(object)[priority$positions + 1])
  stored$cellucid_original_position <- priority$positions
  methods::slot(stored, "misc")$cellucid <- list(priority = priority$record)
  back <- cellucid_restore_original_order(stored)
  expect_identical(colnames(back), colnames(object))
  for (assay in SeuratObject::Assays(object)) {
    expect_identical(assay_matrices(back[[assay]]), assay_matrices(object[[assay]]), info = assay)
  }
  for (reduction in SeuratObject::Reductions(object)) {
    expect_identical(SeuratObject::Embeddings(back, reduction), SeuratObject::Embeddings(object, reduction), info = reduction)
  }
  for (graph in SeuratObject::Graphs(object)) {
    expect_identical(as.matrix(back[[graph]]), as.matrix(object[[graph]]), info = graph)
  }
  expect_identical(SeuratObject::Idents(back), SeuratObject::Idents(object))
  expect_identical(back[[]], object[[]])
  expect_null(methods::slot(back, "misc")$cellucid)
  expect_true(methods::validObject(back))
  # A neighbour index points at its cells by row: the rows move, and so do the indices.
  neighbors <- methods::new("Neighbor", nn.idx = cbind(seq_len(n), c(2:n, 1L)), nn.dist = cbind(rep(0, n), as.double(seq_len(n))),
                            alg.idx = NULL, alg.info = list(), cell.names = colnames(object))
  object[["RNA.nn"]] <- neighbors
  stored <- cellucid:::.seurat_in_cell_order(object, priority$positions + 1)
  stored$cellucid_original_position <- priority$positions
  methods::slot(stored, "misc")$cellucid <- list(priority = priority$record)
  back <- cellucid_restore_original_order(stored)
  expect_identical(SeuratObject::Indices(back[["RNA.nn"]]), SeuratObject::Indices(neighbors))
  expect_identical(SeuratObject::Distances(back[["RNA.nn"]]), SeuratObject::Distances(neighbors))
  # A session is refused on the stored order, applied after the restore.
  session <- testthat::test_path("fixtures", "session_format2", "session.cellucid-session")
  expect_error(cellucid_apply_session(stored, session, expected_dataset_id = "format2"), "cellucid_restore_original_order")
})

test_that("a MuData's modalities are ordered by their cells' priority, prefix counts in global levels", {
  skip_without_sce()
  object <- small_sce(alt = TRUE)
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "m", dataset_name = "M"))
  levels <- store_array(store, "uns/cellucid/priority/levels")
  for (name in store_json(store, "mod/.zattrs")$`mod-order`) {
    map <- store_array(store, file.path("obsmap", name))
    rows <- map[map > 0]
    expect_true(all(diff(rows) > 0))
    counts_path <- file.path("mod", name, "varm", "cellucid_prefix_counts")
    if (file.exists(file.path(store, counts_path, ".zarray"))) {
      counts <- store_array(store, counts_path)
      x_pointers <- store_array(store, file.path("mod", name, "X/indptr"))
      x_rows <- store_array(store, file.path("mod", name, "X/indices"))
      below <- vapply(levels, function(level) as.double(sum(which(map > 0) - 1 < level)), numeric(1))
      for (gene in seq_len(nrow(counts))) {
        range <- seq_len(x_pointers[[gene + 1L]] - x_pointers[[gene]]) + x_pointers[[gene]]
        expect_identical(counts[gene, ], vapply(below, function(rows_below) sum(x_rows[range] < rows_below), numeric(1)))
      }
    }
    own <- store_array(store, file.path("mod", name, "obs/cellucid_original_position"))
    expect_setequal(own, seq_along(own) - 1)
  }
})

test_that("a session is applied in original order, and an object in a file's stored order is refused", {
  priority <- priority_of(2)
  obs <- data.frame(row = c(1, 2), row.names = c("a", "b"))[priority$positions + 1, , drop = FALSE]
  obs$cellucid_original_position <- priority$positions
  attr(obs, "cellucid") <- list(priority = priority$record)
  session <- testthat::test_path("fixtures", "session_format2", "session.cellucid-session")
  expect_error(cellucid_apply_session(obs, session, expected_dataset_id = "format2"), "cellucid_restore_original_order")
  # A record that does not describe the cells: each row its own original
  # position, as the web app read the file; the warning says why.
  attr(obs, "cellucid") <- list(priority = utils::modifyList(priority$record, list(nObs = 3)))
  warned <- character()
  refused <- tryCatch(
    withCallingHandlers(cellucid_apply_session(obs, session, expected_dataset_id = "format2"),
                        warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) conditionMessage(e)
  )
  expect_true(any(grepl("written for 3 cells, and there are 2", warned)))
  expect_false(grepl("cellucid_restore_original_order", paste(refused, collapse = " ")))
})
