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

Connectivity matrices may carry positive edge weights. They must be exactly
symmetric in topology and weight, have a zero diagonal, and contain no negative,
missing, or infinite values. Sparse inputs must omit stored zero entries.

Continuous observation fields, gene-expression matrices, and vector fields
must contain only finite values. For unquantized float32 output, every nonzero
value after required vector scaling must have magnitude from `2^-149` through
`(2 - 2^-23) * 2^127`; values that would encode as zero or infinity reject the
complete candidate. The reserved quantized missing marker is used only for
categorical outlier quantiles that Cellucid generates as `NaN` when a category
has fewer than `centroid_min_points` cells.

`gene_expression` is always interpreted as cells × genes. Cellucid validates
its row and column counts but never guesses or transposes orientation; for a
square matrix, shape alone cannot reveal a genes × cells input.

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
  dataset_id = "my-dataset",
  dataset_name = "My dataset",
  latent_space = latent,   # cells × dims
  obs = obs,               # data.frame (cells × fields)
  var = var,               # data.frame (genes × fields)
  gene_expression = expr,  # optional: cells × genes
  X_umap_2d = umap2,       # optional
  out_dir = "exports/my_dataset",
  force = TRUE,
  obs_categorical_dtype = "uint16"
)
```

## Links

- [Web app](https://cellucid.com)
- [Documentation](https://cellucid.readthedocs.io/en/latest/user_guide/r_package/index.html)
- Recipes: [Seurat](https://cellucid.readthedocs.io/en/latest/user_guide/r_package/e_integrations_recipes/01_seurat_recipe.html) · [SingleCellExperiment](https://cellucid.readthedocs.io/en/latest/user_guide/r_package/e_integrations_recipes/02_singlecellexperiment_recipe.html)
- Source: [cellucid-r](https://github.com/theislab/cellucid-r) · [cellucid](https://github.com/theislab/cellucid)
- Citation: `citation("cellucid")`

## Community

- Contributing: [CONTRIBUTING.md](CONTRIBUTING.md)
- Code of Conduct: [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)
- Security: [SECURITY.md](SECURITY.md)
- Support: [SUPPORT.md](SUPPORT.md)
- GitHub citation metadata: [CITATION.cff](CITATION.cff)

## License

BSD-3-Clause
