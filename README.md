<p>
  <img src="https://raw.githubusercontent.com/theislab/cellucid-python/main/cellucid-logo.svg" alt="Cellucid logo" width="360">
</p>

[![R-CMD-check](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml)
[![Documentation Status](https://readthedocs.org/projects/cellucid/badge/?version=latest)](https://cellucid.readthedocs.io/en/latest/user_guide/r_package/index.html)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# Cellucid (R)

**Export single-cell data from R to Cellucid — interactive, GPU-accelerated visualization in the browser.**

`cellucid` (repo: `cellucid-r`) writes your embeddings, metadata, gene expression, connectivities, and vector fields
to the on-disk format consumed by the Cellucid web app.

> **Active package version — 0.9.1**
>
> Version 0.9.1 is the active Cellucid for R source and documentation version,
> and it is the CRAN submission release. The CRAN package index is authoritative
> for registry availability; the Installation guide below gives both exact
> installation paths.

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
complete candidate. Quantized continuous fields and their manifest bounds use
the viewer's exact float32 values; a native-double range that collapses to one
float32 value rejects the complete candidate. Individual nonzero source values
may round to zero when the resulting float32 range remains non-collapsed. The
reserved quantized missing
marker is used only for
categorical outlier quantiles that Cellucid generates as `NaN` when a category
has fewer than `centroid_min_points` cells.

`gene_expression` is always interpreted as cells × genes. Cellucid validates
its row and column counts but never guesses or transposes orientation; for a
square matrix, shape alone cannot reveal a genes × cells input.

## Installation

Use the dedicated Installation guide to choose CRAN or GitHub based on current
registry availability. To install the active source directly from the official
GitHub repository:

```r
install.packages("remotes")
remotes::install_github("theislab/cellucid-r")
```

Source installation compiles a small native lock primitive and therefore needs
Rtools on Windows, the Xcode Command Line Tools on macOS, or a C toolchain plus
R development headers on Linux. A platform CRAN binary, when available, does
not require a local compiler.

Optional but recommended (sparse matrices + connectivity export):

```r
install.packages("Matrix")
```

See the dedicated [Installation guide](articles/installation.html)
for requirements, verification, upgrades, and the CRAN availability status.

Each export target has a persistent hidden sibling lock file. Cellucid never
places it inside the published dataset or removes it after success: keeping the
same filesystem identity is what lets independent Python and R exporters
reject concurrent writers safely.

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

For reproducible identity metadata, pass an exact UTC-seconds value such as
`created_at = "2026-07-29T12:34:56Z"`; otherwise the exporter records the
current UTC time. Optional provenance starts with `source_name`.
`source_url` and `source_citation` may then be supplied independently, but
neither is accepted without `source_name`.

## Documentation and ecosystem

- [R package guide](https://cellucid.readthedocs.io/en/latest/user_guide/r_package/index.html)
  with recipes for
  [Seurat](https://cellucid.readthedocs.io/en/latest/user_guide/r_package/e_integrations_recipes/01_seurat_recipe.html)
  and
  [SingleCellExperiment](https://cellucid.readthedocs.io/en/latest/user_guide/r_package/e_integrations_recipes/02_singlecellexperiment_recipe.html)
- [Complete Cellucid documentation](https://cellucid.readthedocs.io/en/latest/)
- [Live web application](https://www.cellucid.com) and
  [web viewer source](https://github.com/theislab/cellucid)
- [Python package](https://github.com/theislab/cellucid-python)
- [Official public demo datasets](https://github.com/theislab/cellucid-datasets)
- [Three custom dataset repository examples](https://github.com/theislab/cellucid-demo-custom-datasets)
- [Community annotation repository](https://github.com/theislab/cellucid-annotation)
- Citation: `citation("cellucid")`

## Community

- Contributing: [CONTRIBUTING.md](CONTRIBUTING.md)
- Code of Conduct: [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)
- Security: [SECURITY.md](SECURITY.md)
- Support: [SUPPORT.md](SUPPORT.md)
- GitHub citation metadata: [CITATION.cff](CITATION.cff)

## License

BSD-3-Clause
