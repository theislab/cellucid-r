<p>
  <img src="https://raw.githubusercontent.com/theislab/cellucid-python/main/cellucid-logo.svg" alt="Cellucid logo" width="360">
</p>

[![CRAN status](https://www.r-pkg.org/badges/version/cellucid)](https://CRAN.R-project.org/package=cellucid)
[![R-CMD-check](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml)
[![Documentation Status](https://readthedocs.org/projects/cellucid/badge/?version=latest)](https://cellucid.readthedocs.io/en/latest/user_guide/r_package/index.html)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# Cellucid (R)

**Export single-cell data from R to Cellucid — interactive, GPU-accelerated visualization in the browser.**

`cellucid` (repo: `cellucid-r`) writes your embeddings, metadata, gene expression, connectivities, and vector fields
to the on-disk format consumed by the Cellucid web app.

## Highlights

- **Exporter-only**: generate a shareable “export folder” and open it in the web app
- **Minimal dependencies**: hard dependency `jsonlite` (optional/recommended: `Matrix`)
- **Flexible inputs**: works with raw matrices/data.frames, with docs recipes for Seurat and SingleCellExperiment
- **Optional extras**: gene expression, connectivity graphs, and vector fields for overlays like velocity/drift

## Installation

Install from GitHub:

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

cellucid_prepare(
  latent_space = latent,   # cells × dims
  obs = obs,               # data.frame (cells × fields)
  var = var,               # data.frame (genes × fields)
  gene_expression = expr,  # optional: cells × genes
  X_umap_2d = umap2,       # optional
  out_dir = "exports/my_dataset",
  force = TRUE
)
```

## Links

- Web app: https://cellucid.com
- R package docs (installation + recipes): https://cellucid.readthedocs.io/en/latest/user_guide/r_package/index.html
- Seurat recipe: https://cellucid.readthedocs.io/en/latest/user_guide/r_package/e_integrations_recipes/01_seurat_recipe.html
- SingleCellExperiment recipe: https://cellucid.readthedocs.io/en/latest/user_guide/r_package/e_integrations_recipes/02_singlecellexperiment_recipe.html
- Source: https://github.com/theislab/cellucid-r
- Viewer: https://github.com/theislab/cellucid
- Citation: `citation("cellucid")` (or `CITATION.cff`)

## License

BSD-3-Clause
