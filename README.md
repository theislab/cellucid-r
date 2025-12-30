# cellucid

Minimal R package for exporting single-cell data to the [Cellucid](https://github.com/theislab/cellucid) visualization format.

## Aims

- **Minimal dependencies**: Only `jsonlite` as a hard dependency
- **Simple API**: `cellucid_prepare()` (also exported as `prepare()`) mirroring `cellucid-python`
- **Bioconductor submission**: Complete package skeleton (vignette, tests, biocViews, NEWS)
- **Raw R data structures**: No Seurat/SingleCellExperiment dependencies - accepts standard matrices and data.frames

## Installation

This package is designed for Bioconductor, but can be installed from GitHub during development:

```r
remotes::install_github("theislab/cellucid-r")
```

## Usage

```r
library(cellucid)

latent <- matrix(c(0, 0,
                   1, 1,
                   2, 2),
                 ncol = 2, byrow = TRUE)

obs <- data.frame(
  cluster = factor(c("A", "A", "B")),
  score = c(0.1, 0.2, 0.3)
)

umap2 <- matrix(c(0, 0,
                  1, 0,
                  0, 1),
                ncol = 2, byrow = TRUE)

expr <- matrix(c(0, 1,
                 2, 3,
                 4, 5),
               nrow = 3, ncol = 2, byrow = TRUE)

var <- data.frame(symbol = c("G1", "G2"))
rownames(var) <- var$symbol

out_dir <- "exports/my_dataset"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cellucid_prepare(
  latent_space = latent,
  obs = obs,
  var = var,
  gene_expression = expr,
  X_umap_2d = umap2,
  out_dir = out_dir,
  centroid_min_points = 1,
  force = TRUE
)
```

### Using with Seurat

`cellucid` intentionally does **not** depend on Seurat. This example shows how to
extract the required inputs from an existing Seurat object.

Important: Seurat expression matrices are typically **genes × cells**, but
`cellucid_prepare()` expects **cells × genes**. Transpose before passing.

```r
if (requireNamespace("Seurat", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE)) {

  # Replace with your object
  seu <- your_seurat_object

  # Embedding (cells × 2)
  umap2 <- Seurat::Embeddings(seu, "umap")[, 1:2, drop = FALSE]

  # Latent space (cells × dims) used for categorical outlier quantiles.
  # PCA is a reasonable default; fall back to UMAP if needed.
  latent <- tryCatch(
    Seurat::Embeddings(seu, "pca"),
    error = function(e) umap2
  )

  # Cell metadata (cells × fields)
  obs <- seu@meta.data

  # Gene expression: Seurat stores features × cells -> transpose to cells × features
  expr_genes_by_cells <- Seurat::GetAssayData(seu, slot = "data")
  expr_cells_by_genes <- Matrix::t(expr_genes_by_cells)

  # Feature metadata must have one row per gene
  var <- data.frame(symbol = rownames(expr_genes_by_cells))
  rownames(var) <- var$symbol

  cellucid_prepare(
    latent_space = latent,
    obs = obs,
    var = var,
    gene_expression = expr_cells_by_genes,
    X_umap_2d = umap2,
    out_dir = "exports/seurat_export",
    force = TRUE
  )
}
```

### Using with SingleCellExperiment

SingleCellExperiment expression assays are also typically **genes × cells**, so
transpose before passing to `cellucid_prepare()`.

```r
if (requireNamespace("SingleCellExperiment", quietly = TRUE) &&
    requireNamespace("SummarizedExperiment", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE)) {

  # Replace with your object
  sce <- your_sce_object

  # Embedding (cells × 2)
  umap2 <- SingleCellExperiment::reducedDim(sce, "UMAP")

  # Latent space (cells × dims)
  latent <- tryCatch(
    SingleCellExperiment::reducedDim(sce, "PCA"),
    error = function(e) umap2
  )

  obs <- as.data.frame(SummarizedExperiment::colData(sce))

  expr_genes_by_cells <- SummarizedExperiment::assay(sce, "logcounts")
  expr_cells_by_genes <- Matrix::t(expr_genes_by_cells)

  var <- as.data.frame(SummarizedExperiment::rowData(sce))
  rownames(var) <- rownames(sce)

  cellucid_prepare(
    latent_space = latent,
    obs = obs,
    var = var,
    gene_expression = expr_cells_by_genes,
    X_umap_2d = umap2,
    out_dir = "exports/sce_export",
    force = TRUE
  )
}
```

## Output

The export directory contains:

- Embeddings: `points_1d.bin`, `points_2d.bin`, `points_3d.bin` (float32; optionally `.gz`)
- Cell metadata: `obs_manifest.json` and `obs/` binaries
- Gene expression: `var_manifest.json` and `var/` binaries (one file per gene)
- Connectivity (optional): `connectivity_manifest.json` and `connectivity/` binaries
- Dataset metadata: `dataset_identity.json`

## Related

- [cellucid](https://github.com/theislab/cellucid) - WebGL viewer
- [cellucid-python](https://github.com/theislab/cellucid-python) - Python package

## Publishing

See `docs/publishing.md`.

## Contributing

See `CONTRIBUTING.md`.
