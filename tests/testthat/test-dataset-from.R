# cellucid_prepare(data = ) reads Seurat (SeuratObject 4 and 5),
# SingleCellExperiment and SummarizedExperiment objects natively, with the
# mapping the web app applies to the same object saved with saveRDS().

store_path <- function(name = "data.zarr") file.path(tempfile("cellucid_object_"), name)

attrs_of <- function(store, key) store_json(store, file.path(key, ".zattrs"))

test_that("a Seurat object is written with its cells, embeddings, assay and graph", {
  skip_without_seurat()
  object <- pbmc_seurat()
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "pbmc", dataset_name = "PBMC", priority_seed = NULL))
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
  # X is the log-normalized data, gene-major CSC; one matrix per modality.
  data <- if (seurat_v5()) SeuratObject::LayerData(object, "data") else SeuratObject::GetAssayData(object, slot = "data")
  expect_equal(store_array(store, "X/indptr"), Matrix::t(data)@p)
  expect_equal(store_array(store, "X/data"), as_float32(Matrix::t(data)@x))
  expect_identical(list.files(file.path(store, "layers")), character(0))
  # var.features become highly_variable.
  expect_identical(
    sum(store_array(store, "var/highly_variable")),
    length(SeuratObject::VariableFeatures(object))
  )
  # The SNN graph, its diagonal removed, as the neighbor graph.
  expect_true(dir.exists(file.path(store, "obsp", "connectivities")))
})

test_that("a Seurat object's own record names it, and arguments override the object", {
  skip_without_seurat()
  object <- pbmc_seurat()
  object@misc$cellucid <- list(
    id = "pbmc-small", name = "PBMC small", citation = list(text = "Satija lab"),
    defaults = list(embedding2d = "X_tsne", latentKey = "X_pca")
  )
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, priority_seed = NULL))
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
    priority_seed = NULL,
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
    cellucid_prepare(store_path(), "zarr", data = object, priority_seed = NULL),
    "dataset_id and dataset_name are required.*misc\\$cellucid"
  )
  object@reductions <- list()
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", priority_seed = NULL),
    "holds no embedding.*RunUMAP"
  )
})

test_that("an old or foreign object is refused with the line that fixes it", {
  skip_without_seurat()
  object <- pbmc_seurat()
  object@version <- package_version("2.3.4")
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", priority_seed = NULL),
    "UpdateSeuratObject"
  )
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = list(a = 1), dataset_id = "a", dataset_name = "A", priority_seed = NULL),
    "must be a Seurat, SingleCellExperiment or SummarizedExperiment object; got a list object"
  )
  expect_error(
    cellucid_prepare(store_path(), "zarr", modalities = "RNA", embeddings = list(X_umap = diag(2)), dataset_id = "a", dataset_name = "A", priority_seed = NULL),
    "modalities chooses among the modalities of data"
  )
})

test_that("Seurat v5 split layers are joined, on a copy", {
  skip_without_seurat()
  skip_if(!seurat_v5(), "SeuratObject 5 is not installed (split layers are a v5 feature)")
  object <- pbmc_seurat()
  object[["RNA"]] <- split(object[["RNA"]], f = object$groups)
  expect_true(any(grepl("^counts\\.", SeuratObject::Layers(object[["RNA"]]))))
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", priority_seed = NULL))
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
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", priority_seed = NULL))
  root <- attrs_of(store, "")
  expect_identical(root$`encoding-type`, "MuData")
  expect_identical(root$axis, 0L)
  expect_identical(unlist(attrs_of(store, "mod")$`mod-order`), c("RNA", "ADT"))
  expect_identical(store_array(store, "mod/ADT/var/_index"), paste0("prot", 1:5))
  expect_identical(store_array(store, "obsmap/ADT"), as.numeric(seq_len(ncol(object))))
  expect_false(dir.exists(file.path(store, "uns/cellucid/defaults")))

  # modality names the modality the viewer opens with.
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", modality = "ADT", priority_seed = NULL))
  expect_identical(store_array(store, "uns/cellucid/defaults/modality"), "ADT")
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", modality = "ATAC", priority_seed = NULL),
    "modality must name one of the modalities: c\\(\"RNA\", \"ADT\"\\)"
  )

  # The genes of one feature space cannot be given for several.
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", gene_identifiers = "CD3E", priority_seed = NULL),
    "c\\(\"gene_identifiers\"\\) describe the genes of one feature space.*modalities = \"RNA\""
  )

  # modalities chooses and orders; one modality is a plain AnnData.
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", modalities = "ADT", priority_seed = NULL))
  expect_identical(attrs_of(store, "")$`encoding-type`, "anndata")
  expect_identical(store_array(store, "var/_index"), paste0("prot", 1:5))
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", modalities = "ATAC", priority_seed = NULL),
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
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, priority_seed = NULL))
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
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "d", dataset_name = "D", priority_seed = NULL))
  expect_equal(store_array(store, "X/data"), as_float32(Matrix::t(methods::as(expected, "CsparseMatrix"))@x))
  expect_message(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "d", dataset_name = "D", priority_seed = NULL),
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
    cellucid_prepare(store_path(), "zarr", data = se, dataset_id = "a", dataset_name = "A", priority_seed = NULL),
    "SummarizedExperiment holds no embedding.*scater::runUMAP"
  )
  umap <- matrix(seq_len(2L * ncol(se)), ncol = 2L, dimnames = list(colnames(se), NULL))
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = se, dataset_id = "a", dataset_name = "A", embeddings = list(X_umap = umap), priority_seed = NULL))
  expect_identical(store_array(store, "obs/_index"), colnames(se))
  object@int_metadata$version <- package_version("1.6.0")
  expect_error(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "a", dataset_name = "A", priority_seed = NULL),
    "updateObject"
  )
})

test_that("an object's record names knn and expression as launch defaults, left out as Python leaves them", {
  declared <- cellucid:::.declared_components(
    list(id = "d", name = "D", defaults = list(knn = "connectivities", expression = "X", latentKey = "X_pca", bogus = 1)),
    "misc$cellucid"
  )
  expect_identical(declared$components$latent_key, "X_pca")
  expect_null(declared$components$knn)
  # Only a key that is no launch default is named.
  expect_identical(declared$notes, "misc$cellucid$defaults: not written, no such launch default: c(\"bogus\").")
})

test_that("an object's record of the cells' totals over all genes keeps its column named", {
  anndata <- cellucid:::.declared_components(
    list(id = "d", name = "D", totalCounts = list(field = "nCount_RNA", over = "all genes")), "misc$cellucid"
  )
  expect_identical(anndata$components$total_counts, "nCount_RNA")
  stored <- cellucid:::.declared_components(
    list(id = "d", name = "D", totalCounts = list(field = "total_counts", over = "stored genes")), "misc$cellucid"
  )
  expect_null(stored$components$total_counts)
  modalities <- cellucid:::.declared_components(
    list(id = "d", name = "D", totalCounts = list(RNA = list(field = "nCount_RNA", over = "all genes"),
                                                  ADT = list(field = "ADT:total_counts", over = "stored genes"))),
    "misc$cellucid"
  )
  expect_identical(modalities$components$total_counts, c(RNA = "nCount_RNA"))
})

test_that("an object's Scanpy colors go with the columns written", {
  skip_without_seurat()
  object <- pbmc_seurat()
  groups <- sort(unique(as.character(object$groups)))
  object@misc$groups_colors <- c("#aa0000", "#00aa00")[seq_along(groups)]
  object@misc$nCount_RNA_colors <- "#000000"
  object@misc$absent_colors <- "#000000"
  said <- character()
  store <- withCallingHandlers(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "pbmc", dataset_name = "PBMC", priority_seed = NULL),
    message = function(message) {
      said <<- c(said, conditionMessage(message))
      invokeRestart("muffleMessage")
    }
  )
  expect_identical(store_array(store, "uns/groups_colors"), object@misc$groups_colors)
  expect_false(dir.exists(file.path(store, "uns/nCount_RNA_colors")))
  expect_false(dir.exists(file.path(store, "uns/absent_colors")))
  expect_true(any(grepl("uns['nCount_RNA_colors'] not written: 'nCount_RNA' is not written as a categorical", said, fixed = TRUE)))
  # Too few for its categories: left out, saying why.
  object@misc$groups_colors <- "#aa0000"
  said <- character()
  store <- withCallingHandlers(
    cellucid_prepare(store_path(), "zarr", data = object, dataset_id = "pbmc", dataset_name = "PBMC", priority_seed = NULL),
    message = function(message) {
      said <<- c(said, conditionMessage(message))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(dir.exists(file.path(store, "uns/groups_colors")))
  expect_true(any(grepl(paste0("uns['groups_colors'] not written: 1 color for the ", length(groups), " categories of 'groups'"), said, fixed = TRUE)))
  # obs_colors replaces them.
  store <- suppressMessages(cellucid_prepare(
    store_path(), "zarr", data = object, dataset_id = "pbmc", dataset_name = "PBMC", priority_seed = NULL,
    obs_colors = list(groups = rep("#123456", length(groups)))
  ))
  expect_identical(store_array(store, "uns/groups_colors"), rep("#123456", length(groups)))
})

test_that("a SingleCellExperiment's colors in metadata(x) are written at the MuData root", {
  skip_without_sce()
  object <- small_sce()
  object@metadata$cellucid <- list(id = "sce", name = "SCE")
  object@metadata$cluster_colors <- c("#010101", "#020202", "#030303")
  store <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = object, priority_seed = NULL))
  expect_identical(store_array(store, "uns/cluster_colors"), c("#010101", "#020202", "#030303"))
  expect_false(dir.exists(file.path(store, "mod/main/uns/cluster_colors")))
  # One feature space in priority order: the rows move, the categories and their colors do not.
  single <- small_sce(alt = FALSE)
  single@metadata$cluster_colors <- c("#010101", "#020202", "#030303")
  one <- suppressMessages(cellucid_prepare(store_path(), "zarr", data = single, dataset_id = "sce", dataset_name = "SCE", priority_seed = 7))
  expect_identical(store_array(one, "obs/cluster/categories"), c("a", "b", "c"))
  expect_identical(store_array(one, "uns/cluster_colors"), c("#010101", "#020202", "#030303"))
  # cellucid_view() and cellucid_serve() serve the same plan, and take obs_colors as prepare does.
  expect_true("obs_colors" %in% cellucid:::.view_plan_arguments)
  served <- suppressMessages(memory_store_of(cellucid:::.dataset_components(
    single, list(dataset_id = "sce", dataset_name = "SCE", obs_colors = list(cluster = c("#0a0a0a", "#0b0b0b", "#0c0c0c")))
  )))
  expect_false(is.null(cellucid:::.store_get(served, "uns/cluster_colors/.zattrs")))
  expect_error(
    cellucid:::.dataset_components(single, list(obs_colors = list(cluster = 1:3))),
    "obs_colors\\$cluster must be a character vector of colors"
  )
})

test_that("the layers of a Seurat assay not written are named, counts beside data included, in v3 and v5 assays", {
  skip_without_seurat()
  set.seed(9)
  counts <- Matrix::rsparsematrix(10L, 12L, 0.4, rand.x = function(n) as.double(rpois(n, 3) + 1))
  dimnames(counts) <- list(paste0("g", 1:10), paste0("c", 1:12))
  logged <- counts
  logged@x <- log1p(logged@x)
  scaled <- matrix(rnorm(36L), 3L, 12L, dimnames = list(paste0("g", 1:3), colnames(counts)))
  creators <- list(Assay = SeuratObject::CreateAssayObject)
  if (seurat_v5()) creators$Assay5 <- SeuratObject::CreateAssay5Object
  for (class in names(creators)) {
    object <- suppressWarnings(SeuratObject::CreateSeuratObject(creators[[class]](counts = counts)))
    object <- set_assay_layer(object, "data", logged)
    object <- set_assay_layer(object, "scale.data", scaled)
    object[["umap"]] <- SeuratObject::CreateDimReducObject(
      embeddings = matrix(rnorm(24L), 12L, 2L, dimnames = list(colnames(counts), c("umap_1", "umap_2"))),
      key = "umap_", assay = "RNA"
    )
    components <- cellucid:::.dataset_from(object)
    expect_match(components$notes, 'assays\\$RNA layers not read: c\\("counts", "scale.data"\\)\\.', all = FALSE, info = class)
  }
})
