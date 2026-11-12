# The round trip: what a person made in the viewer, pulled onto the object in
# one call, and columns pushed into the viewer. The viewer is played by a
# WebSocket client in this process that sends what app/viewer-host-results.js
# sends: one result per result, then the results inventory.

# One result frame as the viewer sends it: header, then cells (0-based) and values.
result_frame <- function(kind, id, meta, cells, values = NULL, type = NULL,
                         revision = 1, n_obs, epoch = "0123456789abcdef") {
  header <- list(
    type = "result", epoch = epoch, revision = revision, kind = kind, id = id, nObs = n_obs,
    meta = meta, cells = if (is.null(cells)) NULL else length(cells), values = type
  )
  text <- as.character(jsonlite::toJSON(header, auto_unbox = TRUE, null = "null", digits = NA))
  if (is.list(meta) && length(meta) == 0L) text <- sub("\"meta\":\\[\\]", "\"meta\":{}", text)
  bytes <- charToRaw(text)
  c(cellucid:::.encode_uint32(length(bytes)), bytes,
    if (length(cells) > 0L) cellucid:::.encode_uint32(cells),
    if (!is.null(values)) cellucid:::.encode_wire_values(values, type))
}

viewer_results <- function(n) {
  list(
    list(kind = "highlight", id = "highlight_1", cells = c(0, 2),
         meta = list(label = "Lasso selection", tool = "lasso", pageId = "page_1", pageName = "Page 1",
                     color = "#e41a1c", enabled = TRUE)),
    list(kind = "categorical", id = "user_cat_1", cells = NULL, type = "uint8",
         values = c(0, 1, 255, rep(0, n - 3)),
         meta = list(key = "annotation", categories = I(c("B", "T")), colors = I(c("#e41a1c", "#377eb8")), origin = "user")),
    list(kind = "continuous", id = "signature:T score", cells = c(1, 3), type = "float32", values = c(0.5, -1),
         meta = list(key = "T score", origin = "analysis")),
    list(kind = "continuous", id = "user_obs_value_1", cells = NULL, type = "float32", values = rep(0, n),
         meta = list(key = "pushed", origin = "notebook")),
    list(kind = "geneSet", id = "gene_set_1", cells = numeric(0),
         meta = list(name = "T", description = "", genes = I("CD3E"), symbols = NULL, origin = "notebook")),
    list(kind = "visibility", id = "visible", cells = c(0, 1, 2), meta = structure(list(), names = character()))
  )
}

# What the viewer sends once its dataset opened: every result, then the
# inventory (app/viewer-host-results.js), over the client's socket.
send_snapshot <- function(viewer, client, held, n, epoch = "0123456789abcdef", dataset = viewer$dataset_id) {
  for (index in seq_along(held)) {
    result <- held[[index]]
    ws_send_binary(client, result_frame(result$kind, result$id, result$meta, result$cells, result$values,
                                        result$type, revision = index, n_obs = n, epoch = epoch))
  }
  items <- lapply(seq_along(held), function(index) list(kind = held[[index]]$kind, id = held[[index]]$id, revision = index))
  ws_send_json(client, list(type = "results", epoch = epoch, dataset = dataset, revision = length(held), nObs = n,
                            items = items))
  run_loop_until(function() identical(viewer$results_epoch, epoch) && !is.null(viewer$results_inventory))
}

test_that("one pull writes every highlight, label and score onto a Seurat object, from what the viewer sent", {
  skip_without_seurat()
  object <- pbmc_seurat()
  n <- ncol(object)
  viewer <- test_viewer(object, dataset_id = "p", dataset_name = "P")
  expect_error(cellucid_pull(object, viewer), "has not sent what it holds yet")
  client <- ws_welcomed(viewer)
  send_snapshot(viewer, client, viewer_results(n), n)
  pulled <- cellucid_pull(object, viewer, return_summary = TRUE)
  meta <- pulled$data@meta.data
  expect_identical(meta$cellucid_highlight__highlight_1[1:4], c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(levels(meta$annotation), c("B", "T"))
  expect_identical(as.character(meta$annotation[1:3]), c("B", "T", NA))
  # The colors the viewer shows it in go with it (RCC-1), where the object keeps uns.
  expect_identical(pulled$data@misc$annotation_colors, c("#e41a1c", "#377eb8"))
  expect_equal(meta[["T score"]][1:4], c(NA, 0.5, NA, -1))
  expect_false("pushed" %in% names(meta))
  expect_false("cellucid_visible" %in% names(meta))
  expect_identical(pulled$summary$skipped[["continuous:user_obs_value_1"]],
                   "pushed from a notebook (include_pushed = TRUE brings it back)")
  expect_identical(pulled$data@misc$cellucid$gene_sets$T$genes, "CD3E")
  expect_setequal(names(pulled$data@misc$cellucid$pulled), c("cellucid_highlight__highlight_1", "annotation", "T score"))

  # The person deletes the highlight: its column goes, the rest is replaced.
  ws_send_json(client, list(type = "result_removed", epoch = "0123456789abcdef", revision = 20, kind = "highlight",
                            id = "highlight_1"))
  run_loop_until(function() is.null(viewer$results[["highlight:highlight_1"]]))
  again <- cellucid_pull(pulled$data, viewer, return_summary = TRUE)
  expect_false("cellucid_highlight__highlight_1" %in% names(again$data@meta.data))
  expect_identical(again$summary$removed, "cellucid_highlight__highlight_1")

  # A column Cellucid did not write is a conflict.
  object$annotation <- "mine"
  expect_error(cellucid_pull(object, viewer), "no pull from the same source wrote them")
  suffixed <- cellucid_pull(object, viewer, on_conflict = "suffix")
  expect_true(all(suffixed$annotation == "mine"))
  expect_true("annotation_2" %in% names(suffixed@meta.data))
  # Pulled again, the same result is written where it was: no annotation_3.
  for (round in 1:3) suffixed <- cellucid_pull(suffixed, viewer, on_conflict = "suffix")
  expect_identical(grep("^annotation", names(suffixed@meta.data), value = TRUE), c("annotation", "annotation_2"))
  replaced <- cellucid_pull(suffixed, viewer, on_conflict = "replace", return_summary = TRUE)
  expect_identical(levels(replaced$data$annotation), c("B", "T"))
  expect_identical(replaced$summary$removed, "annotation_2")
  ws_close(client)
})

test_that("prune removes only what the viewer no longer holds", {
  skip_without_seurat()
  object <- pbmc_seurat()
  n <- ncol(object)
  viewer <- test_viewer(object, dataset_id = "p", dataset_name = "P")
  client <- ws_welcomed(viewer)
  send_snapshot(viewer, client, viewer_results(n), n)
  pulled <- cellucid_pull(object, viewer, include_pushed = TRUE, visibility = TRUE)
  expect_identical(pulled@misc$cellucid$pulled$annotation,
                   list(kind = "categorical", id = "user_cat_1", pulled_from = "viewer", key = "annotation", origin = "user"))
  expect_identical(pulled@misc$cellucid$gene_sets$T,
                   list(genes = "CD3E", description = "", origin = "notebook", pulled_from = "viewer",
                        id = "gene_set_1"))
  # Results left out of a pull keep their columns.
  again <- cellucid_pull(pulled, viewer, highlights = FALSE, return_summary = TRUE)
  expect_identical(again$summary$removed, character())
  expect_true(all(c("cellucid_highlight__highlight_1", "pushed", "cellucid_visible") %in% names(again$data@meta.data)))
  # A file's columns are not the viewer's to remove; a gene set gone from the
  # viewer goes from the record too.
  folder <- tempfile("exports_")
  dir.create(folder)
  writeLines(c("obs_name,my_groups", "kind,categorical", 'categories,"[""a""]"', paste0(colnames(object)[1], ",a")), file.path(folder, "my_groups.csv"))
  withfile <- cellucid_pull(again$data, file.path(folder, "my_groups.csv"))
  expect_identical(withfile@misc$cellucid$pulled$my_groups$pulled_from, "my_groups.csv")
  ws_send_json(client, list(type = "result_removed", epoch = "0123456789abcdef", revision = 30, kind = "geneSet",
                            id = "gene_set_1"))
  run_loop_until(function() is.null(viewer$results[["geneSet:gene_set_1"]]))
  last <- cellucid_pull(withfile, viewer, return_summary = TRUE)
  expect_true("my_groups" %in% names(last$data@meta.data))
  expect_identical(last$summary$removed, character())
  expect_identical(last$summary$removed_gene_sets, "T")
  expect_length(last$data@misc$cellucid$gene_sets, 0L)
  ws_close(client)
})

test_that("a new page replaces the old one's results; another dataset's are refused", {
  skip_without_seurat()
  object <- pbmc_seurat()
  n <- ncol(object)
  viewer <- test_viewer(object, dataset_id = "p", dataset_name = "P")
  client <- ws_welcomed(viewer)
  lasso <- function(cells) {
    list(kind = "highlight", id = "highlight_1", cells = cells,
         meta = list(label = "L", tool = "lasso", pageId = "page_1", pageName = "Page 1", color = "#e41a1c",
                     enabled = TRUE))
  }
  visible <- function(cells) list(kind = "visibility", id = "visible", cells = cells, meta = structure(list(), names = character()))
  # Page 1 at revisions 15 and 16; the page reloads and counts from 1 again.
  send_snapshot(viewer, client, list(lasso(c(1, 2, 3)), visible(c(0, 1))), n)
  ws_send_binary(client, result_frame("highlight", "highlight_1", lasso(c(1, 2, 3))$meta, c(1, 2, 3), revision = 15, n_obs = n))
  ws_send_binary(client, result_frame("visibility", "visible", visible(c(0, 1))$meta, c(0, 1), revision = 16, n_obs = n))
  send_snapshot(viewer, client, list(lasso(c(7, 8)), visible(0:8)), n, epoch = "fedcba9876543210")
  held <- cellucid_results(viewer)
  cells <- stats::setNames(lapply(held, `[[`, "cells"), vapply(held, `[[`, "", "kind"))
  expect_identical(cells$highlight, c(8, 9))
  expect_identical(cells$visibility, as.numeric(1:9))

  # The person opens another dataset in the viewer: nothing of it is pulled.
  ws_send_json(client, list(type = "results", epoch = "1111111111111111", dataset = NULL, revision = 0, nObs = 2700,
                            items = list()))
  run_loop_until(function() identical(viewer$results_epoch, "1111111111111111"))
  expect_error(cellucid_pull(object, viewer), "did not open: another one was opened in it. Open 'p'")
  expect_error(cellucid_gene_sets(viewer), "did not open")
  send_snapshot(viewer, client, list(lasso(c(0))), n, epoch = "2222222222222222", dataset = "other")
  expect_error(cellucid_subset(object, viewer, "highlight_1"), "describe the dataset 'other', and this viewer serves 'p'")
  send_snapshot(viewer, client, list(lasso(c(4))), n, epoch = "3333333333333333")
  expect_identical(cellucid_results(viewer)[[1L]]$cells, 5)
  ws_close(client)
})

test_that("push sends the field.push action at once; gene sets, subsets and selections read what the viewer sent", {
  skip_without_seurat()
  object <- pbmc_seurat()
  n <- ncol(object)
  viewer <- test_viewer(object, dataset_id = "p", dataset_name = "P")
  client <- ws_welcomed(viewer)
  send_snapshot(viewer, client, viewer_results(n), n)
  some <- c(2.5, -1)
  names(some) <- colnames(object)[c(5, 2)]
  pushed <- cellucid_push(viewer, some, key = "score")
  frame <- ws_receive(client)
  expect_identical(frame$opcode, 2L)
  expect_identical(frame$json$action, "field.push")
  expect_identical(frame$json$requestId, pushed$id)
  # Cells ascending by original position, each value with its cell, float32.
  expect_identical(frame$json$params, list(key = "score", kind = "continuous", categories = NULL, n_obs = as.integer(n),
                                           cells = c(1, 4), values = c(-1, 2.5)))
  expect_identical(vapply(frame$json$arrays, function(array) array$type, ""), c("uint32", "float32"))
  expect_error(cellucid_push(viewer, c(nope = 1), key = "x"), "not cells of the viewer's dataset")
  cellucid_gene_set_save(viewer, "T cells", c("CD3E", "CD3D"))
  expect_identical(ws_receive(client)$json$params, list(name = "T cells", genes = list("CD3E", "CD3D")))

  expect_identical(cellucid_gene_sets(viewer), list(T = "CD3E"))
  # Named by the var index, the viewer sends no symbols: the names shown are
  # the identifiers.
  expect_identical(cellucid_gene_sets(viewer, symbols = TRUE), list(T = "CD3E"))
  subset <- cellucid_subset(object, viewer, "Lasso selection")
  expect_identical(colnames(subset), colnames(object)[c(1, 3)])
  expect_error(cellucid_subset(object, viewer, "none"), "names 0 highlight groups")
  ws_send_binary(client, c(
    cellucid:::.encode_uint32(nchar("{\"type\":\"selection\",\"source\":\"lasso\"}")),
    charToRaw("{\"type\":\"selection\",\"source\":\"lasso\"}"),
    cellucid:::.encode_uint32(c(1, 4))
  ))
  run_loop_until(function() !is.null(cellucid_selection(viewer)))
  expect_identical(cellucid_selection(viewer)$cell_names, colnames(object)[c(2, 5)])
  expect_identical(colnames(cellucid_subset(object, viewer)), colnames(object)[c(2, 5)])
  # The same cells in another order are refused, never given the wrong cells
  # (a Seurat object keeps its cells' order when subset: its meta.data reordered).
  reordered <- object@meta.data[rev(seq_len(n)), , drop = FALSE]
  expect_error(cellucid_subset(reordered, viewer), "another order")
  expect_error(cellucid_pull(reordered, viewer), "another order")
  ws_close(client)
})

test_that("an object whose cells are not the viewer's in order is refused by pull and subset", {
  names <- paste0("cell-", 1:6)
  by_names <- new.env()
  by_names$cell_names <- names
  in_order <- data.frame(score = 1:6, row.names = names)
  expect_silent(cellucid:::.require_viewer_order(in_order, by_names))
  expect_error(cellucid:::.require_viewer_order(in_order[6:1, , drop = FALSE], by_names),
               "another order.*Reorder it first")
  other <- data.frame(score = 1:6, row.names = paste0("other-", 1:6))
  expect_error(cellucid:::.require_viewer_order(other, by_names), "6 of this object's cells are not cells of the viewer's dataset")
  # A file served from disk: its names are not read here, so an object in a
  # Cellucid-ready file's stored order is told by its priority record.
  from_disk <- new.env()
  from_disk$cell_names <- NULL
  positions <- order(cellucid:::.subsample_rank(0, 0:5), method = "radix") - 1
  stored <- data.frame(cellucid_original_position = positions, row.names = names)
  attr(stored, "cellucid") <- list(priority = list(seed = 0, specVersion = 1, nObs = 6, levels = cellucid:::.prefix_levels(6)))
  expect_error(cellucid:::.require_viewer_order(stored, from_disk), "cellucid_restore_original_order")
  # A position column with no record names no order: the rows are taken as they are.
  attr(stored, "cellucid") <- NULL
  expect_silent(cellucid:::.require_viewer_order(stored, from_disk))
  expect_silent(cellucid:::.require_viewer_order(in_order, from_disk))
  # The object the viewer serves keeps a stale position column out of the
  # store: in its own order, it is taken.
  expect_silent(cellucid:::.require_viewer_order(cbind(in_order, cellucid_original_position = 5:0), by_names))
})

test_that("a subsampled viewer leaves its unopened cells NA", {
  n <- 6
  open <- list(kind = "subset", id = "open", cells = c(0, 2, 3, 4), meta = list(digest = "0123456789abcdef"))
  results <- lapply(c(list(open), viewer_results(n)[1:2]), function(result) {
    frame <- result_frame(result$kind, result$id, result$meta, result$cells,
                          if (identical(result$kind, "categorical")) c(0, 1, 255, 1) else NULL,
                          result$type, n_obs = n,
                          revision = 1)
    if (identical(result$kind, "categorical")) {
      frame <- result_frame(result$kind, result$id, result$meta, open$cells, c(0, 1, 255, 1), "uint8", n_obs = n)
    }
    cellucid:::.result_record(cellucid:::.parse_viewer_binary_event(frame)$payload, NULL)
  })
  frame <- data.frame(row.names = paste0("cell-", 1:6))
  options <- cellucid:::.pull_options(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, "cellucid_highlight__", "", "error", TRUE)
  planned <- cellucid:::.result_candidates(frame, results, options)
  columns <- stats::setNames(lapply(planned$candidates, `[[`, "values"), vapply(planned$candidates, `[[`, "", "wanted"))
  expect_identical(columns$cellucid_highlight__highlight_1, c(TRUE, NA, TRUE, FALSE, FALSE, NA))
  expect_identical(as.character(columns$annotation), c("B", NA, "T", NA, "T", NA))
})

test_that("a session file pulls the state it shows, its columns named with every other one of the pull", {
  skip_without_sce()
  expected <- jsonlite::read_json(testthat::test_path("fixtures", "history_golden", "history_golden.expected.json"), simplifyVector = FALSE)
  dataset <- expected$dataset
  cells <- unlist(dataset$obs_names)
  genes <- unlist(dataset$var_names)
  obs <- S4Vectors::DataFrame(row.names = cells)
  for (name in names(dataset$obs)) {
    categories <- unlist(dataset$obs[[name]]$categories)
    obs[[name]] <- factor(categories[unlist(dataset$obs[[name]]$codes) + 1L], levels = categories)
  }
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = Matrix::Matrix(0, length(genes), length(cells), sparse = TRUE, dimnames = list(genes, cells))),
    colData = obs,
    reducedDims = list(umap = do.call(rbind, lapply(dataset$obsm$X_umap, function(row) as.double(unlist(row)))))
  )
  path <- testthat::test_path("fixtures", "history_golden", "history_golden.cellucid-session")
  out <- cellucid_pull(object, path, highlight_prefix = "hl_", return_summary = TRUE)
  added <- SummarizedExperiment::colData(out$data)
  expect_true(all(c("hl_highlight_1", "cell_type (merged)", "score") %in% names(added)))
  # A pull writes columns; the labels a session renames are cellucid_apply_session()'s.
  expect_identical(levels(added$cell_type), c("B cell", "T cell", "NK cell"))
  record <- S4Vectors::metadata(out$data)$cellucid
  expect_identical(record$pulled$hl_highlight_1$pulled_from, "history_golden.cellucid-session")
  expect_identical(record$gene_sets[["T cell markers"]]$genes, c("CD3E", "NKG7"))
  expect_identical(record$gene_sets[["T cell markers"]]$origin, "notebook")
  # Pulled again: the same columns, replaced; only the fields asked for.
  again <- cellucid_pull(out$data, path, highlight_prefix = "hl_")
  expect_identical(names(SummarizedExperiment::colData(again)), names(added))
  fields_only <- cellucid_pull(object, path, highlights = FALSE, return_summary = TRUE)
  expect_setequal(names(fields_only$summary$columns), c("score", "cell_type (merged)"))
})

test_that("the web app's exports pull onto a data frame in one call", {
  folder <- tempfile("exports_")
  dir.create(folder)
  writeLines(c("obs_name,cellucid_highlight__highlight_1", "kind,highlight", "categories,",
               paste0("cell-", 1:6, ",", c("true", "false", "true", "false", "false", "false"))),
             file.path(folder, "highlight.csv"))
  writeLines(c("obs_name,annotation", "kind,categorical", 'categories,"[""B"",""T""]"', "cell-1,T", "cell-2,B", "cell-6,T"),
             file.path(folder, "annotation.csv"))
  writeLines(c("obs_name,T score", "kind,continuous", "categories,", "cell-4,1.5", "cell-2,"), file.path(folder, "score.csv"))
  cellucid_write_gmt(list(T = "CD3E"), file.path(folder, "markers.gmt"))
  cells <- data.frame(n = 1:6, row.names = paste0("cell-", 1:6))
  pulled <- cellucid_pull(cells, folder, return_summary = TRUE)
  out <- pulled$data
  expect_identical(out$cellucid_highlight__highlight_1, c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(as.character(out$annotation), c("T", "B", NA, NA, NA, "T"))
  expect_identical(levels(out$annotation), c("B", "T"))
  expect_equal(out[["T score"]], c(NA, NA, NA, 1.5, NA, NA))
  expect_identical(pulled$summary$gene_sets, list(T = "CD3E"))
  expect_identical(attr(out, "cellucid")$gene_sets$T$genes, "CD3E")
  # Again: the columns and gene sets it wrote are replaced, not conflicts.
  again <- cellucid_pull(out, folder, on_conflict = "suffix", return_summary = TRUE)
  expect_setequal(names(again$summary$columns), names(pulled$summary$columns))
  expect_identical(names(attr(again$data, "cellucid")$gene_sets), "T")
  writeLines(c("obs_name,x", "kind,continuous", "categories,", "cell-99,1"), file.path(folder, "bad.csv"))
  expect_error(cellucid_pull(cells, file.path(folder, "bad.csv")), "not in this object")
  writeLines(c("obs_name,x", "cell-1,1"), file.path(folder, "bad.csv"))
  expect_error(cellucid_pull(cells, file.path(folder, "bad.csv")), "kind row")
  writeLines(c("obs_name,x", "kind,continuous", "categories,", "cell-1,one"), file.path(folder, "bad.csv"))
  expect_error(cellucid_pull(cells, file.path(folder, "bad.csv")), "continuous column 'x' holds 'one'")
  writeLines(c("obs_name,x", "kind,label", "categories,", "cell-1,one"), file.path(folder, "bad.csv"))
  expect_error(cellucid_pull(cells, file.path(folder, "bad.csv")), "of kind 'label'")
})

test_that("a field's CSV keeps the field's category order and its empty categories (F125)", {
  # Categories came back in the rows' order of first appearance, and a
  # category no cell held was lost.
  path <- tempfile(fileext = ".csv")
  writeLines(c("obs_name,stage", "kind,categorical",
               'categories,"[""Pre-endocrine"",""Progenitor"",""Endocrine"",""Unassigned""]"',
               "cell-1,Endocrine", "cell-2,Progenitor", "cell-3,Pre-endocrine", "cell-4,"), path)
  cells <- data.frame(n = 1:6, row.names = paste0("cell-", 1:6))
  out <- cellucid_pull(cells, path)
  expect_identical(levels(out$stage), c("Pre-endocrine", "Progenitor", "Endocrine", "Unassigned"))
  expect_identical(as.character(out$stage), c("Endocrine", "Progenitor", "Pre-endocrine", NA, NA, NA))
  # The categories row is the file's word on each column; anything else is refused.
  refusals <- list(
    list('categories,"[""a""]"', "cell-1,b", "categorical column 'x' holds 'b'"),
    list("categories,", "cell-1,a", "categorical column 'x' lists no categories"),
    list('categories,"[""a"",1]"', "cell-1,a", "not a JSON array of strings"),
    list('categories,"[""a"""', "cell-1,a", "not a JSON array of strings"),
    list('categories,"{""a"":1}"', "cell-1,a", "not a JSON array of strings"),
    list('categories,"[""a"",""a""]"', "cell-1,a", "lists the category 'a' twice"),
    # A label is the categories row's exactly: " a" is not "a".
    list('categories,"[""a""]"', "cell-1, a", "holds ' a'")
  )
  for (refusal in refusals) {
    writeLines(c("obs_name,x", "kind,categorical", refusal[[1]], refusal[[2]]), path)
    expect_error(cellucid_pull(cells, path), refusal[[3]], fixed = TRUE)
  }
  writeLines(c("obs_name,x", "kind,continuous", 'categories,"[""a""]"', "cell-1,1"), path)
  expect_error(cellucid_pull(cells, path), "continuous column 'x' lists categories", fixed = TRUE)
  # The two-row layout (no categories row) is not an obs-column CSV.
  writeLines(c("obs_name,x", "kind,categorical", "cell-1,a"), path)
  expect_error(cellucid_pull(cells, path), "categories row", fixed = TRUE)
})

shared_fixture <- function(name, sha256) {
  expect_pinned_fixture(name, sha256)
  cellucid_fixture_path(name)
}

test_that("the shared obs-column CSV fixture keeps every declared kind", {
  path <- shared_fixture("obs_columns.csv", "0353146eabba6a9fa76c565a5710dcce275471fcbfd7674f3a86d8c6e2799db2")
  cells <- data.frame(n = 1:6, row.names = c("c0", "c,1", "c2", "c3", "c4", "c5"))
  pulled <- cellucid_pull(cells, path, return_summary = TRUE)
  out <- pulled$data
  expect_identical(out$cellucid_highlight__highlight_1, c(TRUE, FALSE, TRUE, NA, FALSE, NA))
  # Labels that read as numbers, or as True/False, stay labels; the levels
  # are the categories row's, in its order (not the rows'), one no cell
  # holds among them.
  expect_identical(levels(out$leiden), c("10", "1", "0", 'unused, "quoted"'))
  expect_identical(as.character(out$leiden), c("0", "1", "10", "1", NA, NA))
  expect_identical(levels(out$flag), c("False", "True"))
  expect_identical(as.character(out$flag), c("True", "False", "True", NA, "False", NA))
  expect_equal(out$score, c(0.5, NA, -1.25, 3, 1000, NA))
})

test_that("GMT files round trip, by the rules the web app and Python share", {
  path <- tempfile(fileext = ".gmt")
  cellucid_write_gmt(list(`T cells` = c("CD3E", "CD3D", "CD3E"), `B\tcells` = "MS4A1", `B\r\ncells` = "CD79A", `a  b` = "X"),
                     path, descriptions = list(`T cells` = "markers"))
  expect_identical(readLines(path), c("T cells\tmarkers\tCD3E\tCD3D", "B cells\t\tMS4A1", "B cells (2)\t\tCD79A", "a  b\t\tX"))
  expect_identical(cellucid_read_gmt(path),
                   list(`T cells` = c("CD3E", "CD3D"), `B cells` = "MS4A1", `B cells (2)` = "CD79A", `a  b` = "X"))
  fixture <- shared_fixture("gene_sets.gmt", "071d3e007164ff4d9c311820ba3a1be3061b7fd4b80e59142df17690c37b8bfd")
  expect_identical(cellucid_read_gmt(fixture),
                   list(`T cells` = c("CD3E", "CD3D"), `Set 2` = "X", `T cells (2)` = "CD2", B = "MS4A1"))
})

test_that("gene sets of one name are named by one rule; a logical column is pushed as false/true", {
  sets <- cellucid:::.gene_sets_of(lapply(c("a", "b", "c"), function(id) {
    list(kind = "geneSet", id = id, meta = list(name = "T", description = "", genes = list(id), symbols = NULL, origin = "markers"))
  }))
  expect_identical(vapply(sets, `[[`, "", "name"), c("T", "T (2)", "T (3)"))
  # The names shown go beside the identifiers only when the viewer sent them
  # (a var column names the genes in its data mapping).
  named <- cellucid:::.gene_sets_of(list(list(kind = "geneSet", id = "g", meta = list(
    name = "S", description = "", genes = list("ENSG00000198851"), symbols = list("CD3E"), origin = "markers"
  ))))
  expect_identical(named[[1L]]$entry$genes, "ENSG00000198851")
  expect_identical(named[[1L]]$entry$symbols, "CD3E")
  expect_null(sets[[1L]]$entry$symbols)
  params <- cellucid:::.field_push_params(c(TRUE, FALSE, NA), "flag", 3, NULL)
  expect_identical(unlist(params$categories), c("false", "true"))
  expect_identical(params$values$values, c(1, 0, 255))
  expect_identical(params$values$type, "uint8")
})

test_that("a pushed column is the field.push action, encoded as the protocol's binary action", {
  params <- cellucid:::.field_push_params(factor(c("a", NA, "b")), "k", 3, NULL)
  expect_identical(params$values$values, c(0, 255, 1))
  frame <- cellucid:::.encode_viewer_action("r1", "field.push", params)
  command <- binary_command(frame)
  expect_identical(command$params, list(key = "k", kind = "categorical", categories = list("a", "b"), n_obs = 3L,
                                        cells = NULL, values = c(0, 255, 1)))
  expect_identical(vapply(command$arrays, function(array) array$type, ""), "uint8")
  # Past 255 categories the codes are 16-bit, their largest value no category.
  wide <- cellucid:::.field_push_params(factor(c(paste0("c", 1:300), NA)), "k", 301, NULL)
  expect_identical(wide$values$type, "uint16")
  expect_identical(wide$values$values[[301L]], 65535)
  expect_error(cellucid:::.field_push_params(c(1, Inf, 0), "k", 3, NULL), "infinity")
  expect_error(cellucid:::.field_push_params(1:2, "k", 3, NULL), "2 values for the viewer's 3 cells")
  # Before the dataset's cell count is known (a file only the browser reads,
  # before the viewer said), the column names every cell in order and the
  # viewer refuses one of another length.
  unknown <- cellucid:::.field_push_params(c(0.5, 1), "k", NULL, NULL)
  expect_identical(unknown$n_obs, 2L)
  expect_null(unknown$cells)
})

test_that("Run All: a push and a data mapping sent before the viewer opened go once it connects, answered by the viewer", {
  skip_without_seurat()
  object <- pbmc_seurat()
  n <- ncol(object)
  viewer <- test_viewer(object, dataset_id = "p", dataset_name = "P")
  # As for a format only the browser reads: the cell count is not known yet.
  viewer$n_cells <- NULL
  score <- cellucid_push(viewer, seq_len(n) / n, key = "score")
  short <- cellucid_push(viewer, 1:2, key = "short")
  mapping <- cellucid_mapping_set(viewer, expression = "X")
  expect_identical(score$state, "waiting for the viewer")
  client <- ws_welcomed(viewer)
  sent <- lapply(1:3, function(i) ws_receive(client)$json)
  expect_identical(vapply(sent, function(command) command$action, ""), c("field.push", "field.push", "mapping.set"))
  expect_identical(sent[[1L]]$params$n_obs, as.integer(n))
  expect_identical(sent[[2L]]$params$n_obs, 2L)
  expect_identical(sent[[3L]]$params, list(expression = "X"))
  # The viewer refuses the column of another dataset, by its request.
  ws_send_json(client, list(type = "command_error", command = "action", requestId = short$id,
                            reason = "The column \"short\" has 2 cells and this dataset 80 cells: it belongs to another dataset."))
  expect_warning(run_loop_until(function() short$done), "The viewer did not take cellucid_push\\(\\): The column \"short\" has 2 cells")
  expect_error(cellucid_action_result(short), "belongs to another dataset")
  ws_close(client)
})
