<p>
  <img src="https://raw.githubusercontent.com/theislab/cellucid-python/main/cellucid-logo.svg" alt="Cellucid logo" width="360">
</p>

[![CRAN status](https://www.r-pkg.org/badges/version/cellucid)](https://CRAN.R-project.org/package=cellucid)
[![R-CMD-check](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml)
[![BiocCheck](https://github.com/theislab/cellucid-r/actions/workflows/bioccheck.yaml/badge.svg)](https://github.com/theislab/cellucid-r/actions/workflows/bioccheck.yaml)
[![pkgdown](https://github.com/theislab/cellucid-r/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/theislab/cellucid-r/actions/workflows/pkgdown.yaml)
[![test-coverage](https://github.com/theislab/cellucid-r/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/theislab/cellucid-r/actions/workflows/test-coverage.yaml)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# Cellucid

R package for exporting single-cell data to the [Cellucid](https://github.com/theislab/cellucid) visualization format.

## Features

- **Minimal dependencies**: only `jsonlite` is required
- **Simple API**: `cellucid_prepare()` (also exported as `prepare()`) mirroring `cellucid-python`
- **Raw R data structures**: accepts standard matrices and `data.frame`s (no Seurat/SingleCellExperiment dependency)
- **Optional exports**: gene expression (`var/`), connectivity graphs (`connectivity/`), vector fields (`vectors/`)
- **Bioconductor-ready skeleton**: vignette, tests, `biocViews`, `NEWS.md`, CI workflows

## Installation

Install from GitHub (recommended for now):

```r
install.packages("remotes")
remotes::install_github("theislab/cellucid-r")
```

Optional but recommended (sparse matrices + connectivity export):

```r
install.packages("Matrix")
```

## Quickstart

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
               nrow = 3, ncol = 2, byrow = TRUE) # cells × genes

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

## Seurat / SingleCellExperiment

`cellucid` intentionally does **not** depend on Seurat or SingleCellExperiment. These examples show how to extract the required inputs.

Important: expression matrices are typically **genes × cells**, but `cellucid_prepare()` expects **cells × genes**. Transpose before passing.

### Seurat

```r
if (requireNamespace("Seurat", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE)) {

  seu <- your_seurat_object

  umap2 <- Seurat::Embeddings(seu, "umap")[, 1:2, drop = FALSE]
  latent <- tryCatch(Seurat::Embeddings(seu, "pca"), error = function(e) umap2)
  obs <- seu@meta.data

  expr_genes_by_cells <- Seurat::GetAssayData(seu, slot = "data")
  expr_cells_by_genes <- Matrix::t(expr_genes_by_cells)

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

### SingleCellExperiment

```r
if (requireNamespace("SingleCellExperiment", quietly = TRUE) &&
    requireNamespace("SummarizedExperiment", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE)) {

  sce <- your_sce_object

  umap2 <- SingleCellExperiment::reducedDim(sce, "UMAP")
  latent <- tryCatch(SingleCellExperiment::reducedDim(sce, "PCA"), error = function(e) umap2)
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
- Gene expression (optional): `var_manifest.json` and `var/` binaries (one file per gene)
- Connectivity (optional): `connectivity_manifest.json` and `connectivity/` binaries
- Vector fields (optional): `vectors/` binaries + metadata in `dataset_identity.json`
- Dataset metadata: `dataset_identity.json`

## Documentation

- Function docs: `?cellucid_prepare`
- Website: https://theislab.github.io/cellucid-r/
- Vignette: `browseVignettes("cellucid")` (source: `vignettes/cellucid.Rmd`)
- Publishing checklist: `docs/publishing.md`

## Citation

- In R: `citation("cellucid")`
- GitHub: `CITATION.cff`

## Related

- [cellucid](https://github.com/theislab/cellucid) - WebGL viewer
- [cellucid-python](https://github.com/theislab/cellucid-python) - Python package + docs

## Contributing

See `CONTRIBUTING.md`.
