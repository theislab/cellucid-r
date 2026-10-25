# cellucid_prepare(data = ) reads Seurat (SeuratObject 4 and 5),
# SingleCellExperiment and SummarizedExperiment objects natively, with the
# mapping the web app applies to the same object saved with saveRDS().

store_path <- function(name = "data.zarr") file.path(tempfile("cellucid_object_"), name)

attrs_of <- function(store, key) store_json(store, file.path(key, ".zattrs"))

test_that("a Seurat object is written with its cells, embeddings, assay and graph", {
  skip_without_seurat()
  object <- pbmc_seurat()
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "pbmc", dataset_name = "PBMC"))
  expect_identical(store_array(store, "obs/_index"), colnames(object))
  expect_identical(store_array(store, "var/_index"), rownames(object))
  expect_identical(attrs_of(store, "")$`encoding-type`, "anndata")
  obsm <- list.dirs(file.path(store, "obsm"), full.names = FALSE, recursive = FALSE)
  expect_setequal(obsm, c("X_pca", "X_tsne"))
  tsne <- unname(SeuratObject::Embeddings(object, "tsne"))
  expect_equal(store_array(store, "obsm/X_tsne"), matrix(as_float32(tsne), nrow(tsne)))
  # meta.data columns, and the active identity.
  columns <- unlist(attrs_of(store, "obs")$`column-order`)
  expect_true(all(c("nCount_RNA", "groups", "RNA_snn_res.1", "active_ident") %in% columns))
  # X is the log-normalised data, gene-major CSC; one matrix per modality.
  data <- if (seurat_v5()) SeuratObject::LayerData(object, "data") else SeuratObject::GetAssayData(object, slot = "data")
  expect_equal(store_array(store, "X/indptr"), Matrix::t(data)@p)
  expect_equal(store_array(store, "X/data"), as_float32(Matrix::t(data)@x))
  expect_identical(list.files(file.path(store, "layers")), character(0))
  # var.features become highly_variable.
  expect_identical(
    sum(store_array(store, "var/highly_variable")),
    length(SeuratObject::VariableFeatures(object))
  )
  # The SNN graph, its diagonal removed, as the neighbour graph.
  expect_true(dir.exists(file.path(store, "obsp", "connectivities")))
})

test_that("a Seurat object's own record names it, and arguments override the object", {
  skip_without_seurat()
  object <- pbmc_seurat()
  object@misc$cellucid <- list(
    id = "pbmc-small", name = "PBMC small", citation = list(text = "Satija lab"),
    defaults = list(embedding2d = "X_tsne", latentKey = "X_pca")
  )
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object))
  expect_identical(store_array(store, "uns/cellucid/id"), "pbmc-small")
  expect_identical(store_array(store, "uns/cellucid/name"), "PBMC small")
  expect_identical(store_array(store, "uns/cellucid/citation/text"), "Satija lab")
  expect_identical(store_array(store, "uns/cellucid/defaults/embedding2d"), "X_tsne")
  expect_true(dir.exists(file.path(store, "uns/cellucid/outlierQuantiles")))

  # An argument replaces what the object holds; the launch defaults that
  # name its embeddings go with them.
  object@misc$cellucid$defaults <- NULL
  obs <- data.frame(kind = rep(c("x", "y"), length.out = ncol(object)), row.names = colnames(object))
  store <- suppressMessages(cellucid_prepare(
    store_path(), "zarr", data = object, dataset_name = "Renamed", obs = obs,
    embeddings = list(X_custom = matrix(1, ncol(object), 2, dimnames = list(colnames(object), NULL)))
  ))
  expect_identical(store_array(store, "uns/cellucid/name"), "Renamed")
  expect_identical(unlist(attrs_of(store, "obs")$`column-order`), "kind")
  expect_identical(list.dirs(file.path(store, "obsm"), full.names = FALSE, recursive = FALSE), "X_custom")
})

test_that("an object without an identity or an embedding says what to add", {
  skip_without_seurat()
  object <- pbmc_seurat()
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object),
    "dataset_id and dataset_name are required.*misc\\$cellucid"
  )
  object@reductions <- list()
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A"),
    "holds no embedding.*RunUMAP"
  )
})

test_that("an old or foreign object is refused with the line that fixes it", {
  skip_without_seurat()
  object <- pbmc_seurat()
  object@version <- package_version("2.3.4")
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A"),
    "UpdateSeuratObject"
  )
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = list(a = 1), dataset_id = "a", dataset_name = "A"),
    "must be a Seurat, SingleCellExperiment or SummarizedExperiment object; got a list object"
  )
  expect_error(
    cellucid_prepare(store_path(), "zarr", modalities = "RNA", embeddings = list(X_umap = diag(2)), dataset_id = "a", dataset_name = "A"),
    "modalities chooses among the modalities of data"
  )
})

test_that("Seurat v5 split layers are joined, on a copy", {
  skip_without_seurat()
  skip_if(!seurat_v5(), "SeuratObject 5 is not installed (split layers are a v5 feature)")
  object <- pbmc_seurat()
  object[["RNA"]] <- split(object[["RNA"]], f = object$groups)
  expect_true(any(grepl("^counts\\.", SeuratObject::Layers(object[["RNA"]]))))
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A"))
  expect_identical(store_array(store, "obs/_index"), colnames(object))
  joined <- SeuratObject::JoinLayers(object[["RNA"]])
  data <- SeuratObject::LayerData(joined, "data")[, colnames(object)]
  expect_equal(store_array(store, "X/data"), as_float32(Matrix::t(data)@x))
  # The caller's object still has its split layers.
  expect_true(any(grepl("^counts\\.", SeuratObject::Layers(object[["RNA"]]))))
})

test_that("several assays are written as a Cellucid-ready MuData, the active one first", {
  skip_without_seurat()
  object <- pbmc_multimodal()
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A"))
  root <- attrs_of(store, "")
  expect_identical(root$`encoding-type`, "MuData")
  expect_identical(root$axis, 0L)
  expect_identical(unlist(attrs_of(store, "mod")$`mod-order`), c("RNA", "ADT"))
  expect_identical(store_array(store, "mod/ADT/var/_index"), paste0("prot", 1:5))
  expect_identical(store_array(store, "obsmap/ADT"), as.numeric(seq_len(ncol(object))))
  expect_false(dir.exists(file.path(store, "uns/cellucid/defaults")))

  # modality names the modality the viewer opens with.
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", modality = "ADT"))
  expect_identical(store_array(store, "uns/cellucid/defaults/modality"), "ADT")
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", modality = "ATAC"),
    "modality must name one of the modalities: c\\(\"RNA\", \"ADT\"\\)"
  )

  # The genes of one feature space cannot be given for several.
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", gene_identifiers = "CD3E"),
    "c\\(\"gene_identifiers\"\\) describe the genes of one feature space.*modalities = \"RNA\""
  )

  # modalities chooses and orders; one modality is a plain AnnData.
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", modalities = "ADT"))
  expect_identical(attrs_of(store, "")$`encoding-type`, "anndata")
  expect_identical(store_array(store, "var/_index"), paste0("prot", 1:5))
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", modalities = "ATAC"),
    "modalities names c\\(\"ATAC\"\\), which data does not hold"
  )
})

test_that("a modality measured on some cells is placed by obsmap", {
  skip_if_not_installed("Matrix")
  cells <- paste0("c", 1:6)
  rna <- Matrix::Matrix(matrix(1:18, 3, 6, dimnames = list(paste0("g", 1:3), cells)), sparse = TRUE)
  adt <- Matrix::Matrix(matrix(1:4, 2, 2, dimnames = list(c("p1", "p2"), cells[c(2, 5)])), sparse = TRUE)
  components <- list(
    dataset_id = "part", dataset_name = "Partial",
    embeddings = list(X_umap = matrix(seq_len(12), 6, 2, dimnames = list(cells, NULL))),
    modalities = list(
      list(name = "RNA", gene_expression = cellucid:::.feature_major(rna), var = NULL, cells = NULL),
      list(name = "ADT", gene_expression = cellucid:::.feature_major(adt), var = NULL, cells = cells[c(2, 5)])
    ),
    notes = character()
  )
  folder <- store_path()
  plan <- cellucid:::.plan_dataset(components)
  target <- cellucid:::.resolve_output_path(folder, "zarr", FALSE)
  target$zarr_format <- 2L
  cellucid:::.write_output(target, plan$write)
  expect_identical(store_array(folder, "obsmap/ADT"), c(0, 1, 0, 0, 2, 0))
  expect_identical(store_array(folder, "obsm/ADT"), matrix(c(0, 1, 0, 0, 1, 0), ncol = 1L))
  expect_identical(store_array(folder, "obsmap/RNA"), as.numeric(1:6))
  expect_identical(store_array(folder, "mod/ADT/obs/_index"), cells[c(2, 5)])
  expect_identical(store_array(folder, "varmap/ADT"), c(0, 0, 0, 1, 2))
  expect_identical(store_array(folder, "var/_index"), c("g1", "g2", "g3", "p1", "p2"))
})

test_that("a SingleCellExperiment is written with its assays, reducedDims, altExps and colPairs", {
  skip_without_sce()
  object <- small_sce()
  object@metadata$cellucid <- list(id = "sce", name = "SCE")
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object))
  expect_identical(attrs_of(store, "")$`encoding-type`, "MuData")
  expect_identical(unlist(attrs_of(store, "mod")$`mod-order`), c("main", "ADT"))
  # reducedDims keep their case: X_UMAP, X_PCA; beside them mudata's masks of
  # measured cells, one per modality.
  expect_setequal(list.dirs(file.path(store, "obsm"), full.names = FALSE, recursive = FALSE), c("X_UMAP", "X_PCA", "main", "ADT"))
  expect_identical(store_array(store, "obsm/ADT"), matrix(1, ncol(object), 1L))
  columns <- unlist(attrs_of(store, "obs")$`column-order`)
  # Size factors: colData's sizeFactor where SingleCellExperiment keeps them
  # there, else obs['size_factor'] from int_colData.
  expect_identical(columns[1:3], c("cluster", "depth", "label"))
  factor_column <- intersect(columns, c("sizeFactor", "size_factor"))
  expect_length(factor_column, 1L)
  expect_equal(store_array(store, file.path("obs", factor_column)), SingleCellExperiment::sizeFactors(object))
  # logcounts is X (the other assays are not written); rowData is var.
  logcounts <- SummarizedExperiment::assay(object, "logcounts")
  expect_equal(store_array(store, "mod/main/X/data"), as_float32(Matrix::t(logcounts)@x))
  expect_false(dir.exists(file.path(store, "mod/main/layers/counts")))
  expect_identical(store_array(store, "mod/main/var/symbol"), paste0("S", 1:30))
  expect_true(dir.exists(file.path(store, "obsp/connectivities")))
})

test_that("a DelayedArray assay is read block by block", {
  skip_without_sce()
  skip_if_not_installed("DelayedArray")
  object <- small_sce(alt = FALSE)
  expected <- SummarizedExperiment::assay(object, "logcounts")
  SummarizedExperiment::assay(object, "logcounts") <- DelayedArray::DelayedArray(expected)
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "d", dataset_name = "D"))
  expect_equal(store_array(store, "X/data"), as_float32(Matrix::t(methods::as(expected, "CsparseMatrix"))@x))
  expect_message(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "d", dataset_name = "D"),
    "block by block"
  )
})

test_that("a SummarizedExperiment needs an embedding, and an old SCE is refused", {
  skip_without_sce()
  object <- small_sce(alt = FALSE)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = SummarizedExperiment::assay(object, "counts")),
    colData = SummarizedExperiment::colData(object)
  )
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = se, dataset_id = "a", dataset_name = "A"),
    "SummarizedExperiment holds no embedding.*scater::runUMAP"
  )
  umap <- matrix(seq_len(2L * ncol(se)), ncol = 2L, dimnames = list(colnames(se), NULL))
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = se, dataset_id = "a", dataset_name = "A", embeddings = list(X_umap = umap)))
  expect_identical(store_array(store, "obs/_index"), colnames(se))
  object@int_metadata$version <- package_version("1.6.0")
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A"),
    "updateObject"
  )
})
