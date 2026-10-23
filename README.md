<p>
  <img src="https://raw.githubusercontent.com/theislab/cellucid-python/main/cellucid-logo.svg" alt="Cellucid logo" width="360">
</p>

[![R-CMD-check](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml)
[![Documentation Status](https://readthedocs.org/projects/cellucid/badge/?version=latest)](https://docs.cellucid.com/en/latest/user_guide/r_package/index.html)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# Cellucid (R)

**Prepare single-cell data in R for Cellucid — interactive, GPU-accelerated visualization in the browser.**

`cellucid` (repo: `cellucid-r`) writes your embeddings, metadata, gene expression, neighbour graph, and vector fields
as a **Cellucid-ready AnnData**: a standard AnnData Zarr store that the Cellucid web app opens directly and Python's
`anndata` reads like any other.

> **Active package version — 0.9.5** <!-- CELLUCID_VERSION -->
>
> Version 0.9.5 is the active Cellucid for R source and documentation version,
> and it is the CRAN submission release. The CRAN package index is authoritative
> for registry availability; the Installation guide below gives both exact
> installation paths.

## Highlights

- **One standard format**: an AnnData Zarr store (format 2, or 3 with `zarr_format = 3`) with consolidated
  metadata — a `.zarr` folder for GitHub (`format = "zarr"`), or one `.zarr.zip` file (`format = "zarr.zip"`) for
  everywhere else: single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files,
  `cellucid serve` and Jupyter
- **Fast in the browser**: gene-major expression (CSC), chunked for range reads, byte-shuffled where that helps
  and compressed with deflate, which browsers decompress natively; a `.zarr.zip` keeps its metadata at the end
  of the archive, so it opens with one read
- **Self-describing**: the dataset's id, name, citation, source, and launch defaults live in `uns['cellucid']`,
  with each categorical `obs` column's outlier quantiles in the latent space named by `latent_key`, so the
  viewer's outlier filter works as soon as the file opens
- **Minimal dependencies**: hard dependency `jsonlite` (optional/recommended: `Matrix`)
- **Flexible inputs**: raw matrices and data frames, with recipes for Seurat and SingleCellExperiment

Every value is stored as float32 (embeddings, expression, vector fields, the
graph) and must be finite. A neighbour graph must be exactly symmetric, with a
zero diagonal and non-negative weights; a stored zero in a sparse input is not
an edge and is dropped.

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

Source installation compiles one small C routine (the CRC-32 checksum of the
`.zarr.zip` writer) and therefore needs Rtools on Windows, the Xcode Command
Line Tools on macOS, or a C toolchain plus R development headers on Linux. A
platform CRAN binary, when available, does not require a local compiler.

Optional but recommended (sparse expression matrices and neighbour graphs):

```r
install.packages("Matrix")
```

See the dedicated [Installation guide](https://theislab.github.io/cellucid-r/articles/installation.html)
for requirements, verification, upgrades, and the CRAN availability status.

## Quickstart

```r
library(cellucid)

cellucid_prepare(
  out_path = "my_dataset.zarr",       # or "my_dataset.zarr.zip" with format = "zarr.zip"
  format = "zarr",
  embeddings = list(                  # obsm: any number, every column kept
    X_umap = umap,
    X_pca = pca
  ),
  dataset_id = "my-dataset",
  dataset_name = "My dataset",
  obs = obs,                          # data.frame, one row per cell
  var = var,                          # data.frame, one row per gene
  gene_expression = expr,             # cells × genes, sparse or dense
  connectivities = knn,               # optional cells × cells graph
  embedding_defaults = c("2d" = "X_umap"),
  latent_key = "X_pca"
)
```

Open a `.zarr.zip` file from your computer in the web app at
<https://www.cellucid.com>, or a `.zarr` folder or `.zarr.zip` file published at
a web address with **Open from link**; read either in Python with
`anndata.read_zarr()`.

The store is written beside `out_path` under a temporary name and moved into
place only when it is complete; `force = TRUE` replaces an existing store.
Vector fields are stored under scVelo's names and drawn over the embedding they
are named for: `vector_fields = list(velocity = list(X_umap = v))` stores
`obsm['velocity_umap']`, and `X_branches = vb` stores `obsm['velocity_branches']`,
both listed by the viewer as the field `velocity`. `cellucid_transition_drift()` computes such a field
from a cell-cell transition matrix (CellRank's, for example), the computation of Python's
`compute_transition_drift()`. Repeated cell or gene names are made unique the way
anndata does it (`MATR3`, `MATR3-1`), and `var_gene_id_column` names the `var`
column, such as gene symbols beside Ensembl row names, that the viewer shows.

## Documentation and ecosystem

- [R package guide](https://docs.cellucid.com/en/latest/user_guide/r_package/index.html)
  with recipes for
  [Seurat](https://docs.cellucid.com/en/latest/user_guide/r_package/e_integrations_recipes/01_seurat_recipe.html)
  and
  [SingleCellExperiment](https://docs.cellucid.com/en/latest/user_guide/r_package/e_integrations_recipes/02_singlecellexperiment_recipe.html)
- [Complete Cellucid documentation](https://docs.cellucid.com/en/latest/)
- [Live web application](https://www.cellucid.com) and
  [web viewer source](https://github.com/theislab/cellucid)
- [Python package](https://github.com/theislab/cellucid-python)
- [Official public demo datasets](https://github.com/theislab/cellucid-datasets)
- [Three custom dataset repository examples](https://github.com/theislab/cellucid-demo-custom-datasets)
- [Community annotation repository](https://github.com/theislab/cellucid-annotation)
- Citation: `citation("cellucid")`

## Citation

<!-- BEGIN cellucid-citation (generated by cellucid-private/scripts/sync_citation.py from CITATION.cff; do not edit) source-sha256=2486463f2b7c8cc656ae41e2df0261c1eb79780c4d635be147215f409b89aa2a body-sha256=5186fc40b9eddfb904fab52bee5362bf27a8d6c3713f9b67b1defa1c68c54d3f -->
> **Placeholder citation.** The Cellucid preprint is in preparation for bioRxiv and has no DOI yet; the reference below is a placeholder and will be replaced when the preprint is published.

If you use Cellucid, please cite:

Inecik, K. PLACEHOLDER: Cellucid preprint, title to be announced (forthcoming on bioRxiv). *bioRxiv* (in preparation).

BibTeX and other formats: [Citing Cellucid](https://docs.cellucid.com/en/latest/citing.html).
<!-- END cellucid-citation -->

## Community

- Contributing: [CONTRIBUTING.md](https://github.com/theislab/cellucid-r/blob/main/CONTRIBUTING.md)
- Code of Conduct: [CODE_OF_CONDUCT.md](https://github.com/theislab/cellucid-r/blob/main/CODE_OF_CONDUCT.md)
- Security: [SECURITY.md](https://github.com/theislab/cellucid-r/blob/main/SECURITY.md)
- Support: [SUPPORT.md](https://github.com/theislab/cellucid-r/blob/main/SUPPORT.md)
- GitHub citation metadata: [CITATION.cff](https://github.com/theislab/cellucid-r/blob/main/CITATION.cff)

## License

BSD-3-Clause
