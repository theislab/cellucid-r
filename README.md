<p>
  <img src="https://raw.githubusercontent.com/theislab/cellucid-python/main/cellucid-logo.svg" alt="Cellucid logo" width="360">
</p>

[![R-CMD-check](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/theislab/cellucid-r/actions/workflows/R-CMD-check.yaml)
[![Documentation Status](https://readthedocs.org/projects/cellucid/badge/?version=latest)](https://docs.cellucid.com/en/latest/user_guide/r_package/index.html)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# Cellucid (R)

**View and prepare single-cell data in R for Cellucid — interactive, GPU-accelerated visualization in the browser.**

`cellucid` (repo: `cellucid-r`) opens a Seurat, SingleCellExperiment or SummarizedExperiment object in the Cellucid
viewer from R — in the RStudio or Positron Viewer pane or a browser, with no Python — and writes it (or your
embeddings, metadata, gene expression, neighbour graph and vector fields) as a **Cellucid-ready AnnData**: a standard
AnnData Zarr store that the Cellucid web app opens directly and Python's `anndata` reads like any other. Several assays
become a **Cellucid-ready MuData**.

> **Active package version — 0.9.5** <!-- CELLUCID_VERSION -->
>
> Version 0.9.5 is the active Cellucid for R source and documentation version,
> and it is the CRAN submission release. The CRAN package index is authoritative
> for registry availability; the Installation guide below gives both exact
> installation paths.

## Highlights

- **View an object from R**: `cellucid_view(seu)` serves the object from memory as the store `cellucid_prepare()`
  would write, uncompressed over loopback and making only what the viewer reads (a million-cell object's first
  gene in about a third of a second, never transposing its matrix), so a 100,000-cell object opens in under a
  second; talk to the
  viewer with `cellucid_highlight()`, `cellucid_color_by()`, `cellucid_on("selection", ...)` and
  `cellucid_session()`, and write what was made in it back onto the object with
  `cellucid_apply_session(seu, viewer)`; `cellucid_serve()` for a cluster behind an SSH tunnel
- **Objects in, directly**: `cellucid_prepare(data = seu)` reads Seurat (v3–v5), SingleCellExperiment and
  SummarizedExperiment objects as the web app reads them; any argument given replaces that component
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
- **Light dependencies**: `jsonlite`, and `httpuv`, `later` and `curl` for the viewer (optional: `Matrix`,
  `SeuratObject`, `SingleCellExperiment` / `SummarizedExperiment`, each needed only for its own objects)
- **Flexible inputs**: R objects, or raw matrices and data frames

Every value is stored as float32 (embeddings, expression, vector fields, the
graph; expression as float64 only for integer counts above 2^24, which float32
cannot hold exactly) and must be finite. A neighbour graph must be exactly symmetric, with a
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

Source installation compiles a few small C routines (checksums, the outlier
quantiles, one-pass matrix checks, the decoders of a saved session, linked
against zlib) and therefore needs Rtools on Windows, the Xcode Command Line
Tools on macOS, or a C toolchain plus the R and zlib development headers on
Linux. A
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

# View an object
viewer <- cellucid_view(seu)                 # a Seurat, SingleCellExperiment or SummarizedExperiment
cellucid_color_by(viewer, "seurat_clusters")
cellucid_selection(viewer)                   # the cells you selected in the viewer
cellucid_stop(viewer)

# Write it for sharing
cellucid_prepare("pbmc.zarr", "zarr", data = seu, dataset_id = "pbmc", dataset_name = "PBMC")

# Or write components you hold as matrices

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
a web address with the **Link or server** address bar; read either in Python with
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
  with pages for
  [Seurat](https://docs.cellucid.com/en/latest/user_guide/r_package/e_integrations_recipes/01_seurat_recipe.html)
  and
  [SingleCellExperiment](https://docs.cellucid.com/en/latest/user_guide/r_package/e_integrations_recipes/02_singlecellexperiment_recipe.html)
  objects, and for viewing them from R (Viewing and loading)
- [Complete Cellucid documentation](https://docs.cellucid.com/en/latest/)
- [Live web application](https://www.cellucid.com) and
  [web viewer source](https://github.com/theislab/cellucid)
- [Python package](https://github.com/theislab/cellucid-python)
- [Official public demo datasets](https://github.com/theislab/cellucid-datasets)
- [Custom dataset repository examples](https://github.com/theislab/cellucid-demo-custom-datasets)
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
