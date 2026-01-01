# Contributing to `cellucid` (R)

Contributions are welcome — code, docs, bug reports, and reproducible examples all help.

This file focuses on `cellucid-r` (the R package that exports data to the Cellucid viewer format), and also helps route issues/PRs to the correct repo in the Cellucid ecosystem.

By participating, you agree to follow the project’s Code of Conduct:
- `CODE_OF_CONDUCT.md`

---

## Which repo should I contribute to?

Cellucid is split by responsibility:

| Repo | What it is | Contribute here when you… |
|---|---|---|
| `cellucid` | Web app (UI + state + WebGL rendering) | are fixing UI bugs, rendering/performance, figure export, sessions, or community annotation frontend |
| `cellucid-python` | Python package + CLI (`prepare`, `serve`, `show_anndata`, hooks) + Sphinx docs | are fixing Python/CLI bugs, data prep/export, server endpoints, Jupyter hooks, or docs on ReadTheDocs |
| `cellucid-r` (this repo) | R package for exporting data to the Cellucid viewer format | are changing the R exporter (`cellucid_prepare()`), adding R-side tests/docs, or preparing for Bioconductor |
| `cellucid-annotation` | GitHub repo template for community annotation | are changing the repo schema/validation/workflows |

If you’re not sure where a bug belongs, open an issue in the repo you’re currently using and include:
- how you loaded data (exports vs h5ad/zarr vs remote server vs Jupyter),
- the UI environment (hosted app vs local app vs Jupyter iframe),
- and the smallest reproduction you can share.

---

## Fast paths (pick your contribution type)

### I want to report a bug

Please include:
- R version (`R.version.string`)
- OS (macOS/Windows/Linux)
- how you generated exports (which function call, which inputs)
- expected vs actual behavior
- the smallest dataset you can share (prefer synthetic/reproducible)
- any `R CMD check` output if relevant

If the bug is “the viewer looks wrong”, also include:
- which viewer environment you used (hosted app vs local app vs Jupyter iframe)
- browser + version
- whether you loaded the export via browser file picker vs remote server vs GitHub exports

### I want to contribute docs only

Docs live in:
- `cellucid-r/man/` (generated `.Rd` files)
- `cellucid-r/vignettes/` (BiocStyle vignette)
- `cellucid-r/README.md`

If you edit `.Rd` files directly, be aware they are usually generated from roxygen comments in `cellucid-r/R/`.
Prefer editing roxygen in `R/` and re-running `devtools::document()`.

### I want to add/modify R code

Fast workflow:
1) Set up your R dev environment (see “Development setup”)
2) Make a small, focused change
3) Add/adjust tests (`testthat`)
4) Run `devtools::check()` (and `BiocCheck::BiocCheck()` if relevant)
5) Submit a PR with a clear “what/why/how to verify”

---

## Development setup (cellucid-r)

### Prerequisites

- R `>= 4.3.0` (matches `DESCRIPTION`)
- Git

Recommended:
- RStudio (optional, but convenient for vignettes)
- Pandoc (needed for some vignette builds; RStudio bundles it)

### Clone

```bash
git clone https://github.com/theislab/cellucid-r.git
cd cellucid-r
```

### Install dev tools

In R:

```r
install.packages(c("devtools", "roxygen2", "testthat"))
```

For vignette builds and Bioconductor-style docs:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("BiocStyle", "knitr", "rmarkdown"))
```

Optional (useful during release prep):

```r
BiocManager::install("BiocCheck")
install.packages(c("pkgdown", "covr"))
```

Install package dependencies:

```r
devtools::install_deps(dependencies = TRUE)
```

Load the package from source (fast inner loop):

```r
devtools::load_all()
```

---

## Tests

Run tests:

```r
devtools::test()
```

Run a full check (recommended before PRs):

```r
devtools::check()
```

For Bioconductor submission readiness:

```r
BiocCheck::BiocCheck(".")
```

Guidelines:
- Add tests when behavior changes (especially edge cases like missing values, mismatched dimensions, sparse matrices).
- Prefer small synthetic inputs; avoid committing real datasets.

---

## Documentation (roxygen + vignettes)

### Regenerate `.Rd` files

The package ships its help page(s) under `cellucid-r/man/`.

If you change the public API (especially `cellucid_prepare()`), keep the help page in sync:
- update `cellucid-r/man/cellucid_prepare.Rd`, or
- regenerate it with `roxygen2` and commit the result.

```r
devtools::document()
```

### Build vignettes

```r
devtools::build_vignettes()
```

### pkgdown site (optional)

If you’re working on the website output (GitHub Pages):

```r
pkgdown::build_site(dest_dir = "docs/site")
```

---

## Design constraints (important for contributors)

This package is intentionally:

- **minimal-dependency** (only `jsonlite` is a hard dependency)
- **format-first** (exports must match what the web app expects)
- **Bioconductor-friendly** (checks, vignette style, and package structure matter)

If you propose adding a new dependency:
- prefer `Suggests` over `Imports` unless strictly required
- keep the core exporter usable on minimal installations

If you change the export format or schema:
- coordinate with `cellucid-python` and the web app (`cellucid`) so the ecosystem stays compatible
- add tests that validate the new behavior

---

## How to validate exports end-to-end

The most reliable validation is to export a tiny dataset and load it in the viewer.

Recommended workflow:
1) Run `cellucid_prepare()` on a tiny synthetic dataset.
2) Open the Cellucid web app (hosted or local).
3) Use “Browse local data…” and select the export folder.
4) Confirm:
   - the correct number of cells render
   - categorical and continuous fields appear
   - gene expression search works (if exported)

This catches many “format is technically written but semantically wrong” issues.

---

## PR guidelines

- Keep PRs small and focused (one feature/bugfix at a time).
- Include:
  - what changed
  - why it changed
  - how to verify (commands + expected outcome)
- If user-facing behavior changes:
  - update docs (README/vignette/man page as appropriate)
  - add/adjust tests
  - consider updating `NEWS.md` if it’s release-notable

---

## Troubleshooting (common contributor problems)

### `devtools::document()` changes lots of files unexpectedly

Cause:
- roxygen output depends on roxygen version and formatting.

Fix:
- ensure you’re using a recent `roxygen2`
- keep roxygen comments stable and avoid rewrapping large blocks unnecessarily

### `R CMD check` fails on vignettes

Common causes:
- missing suggested packages (`BiocStyle`, `knitr`, `rmarkdown`)
- missing Pandoc

Fix:
- install suggested packages (see setup above)
- use RStudio’s bundled Pandoc or install Pandoc system-wide

### Windows/macOS differences

Common causes:
- line ending differences
- path handling (`\` vs `/`)

Fix:
- use `file.path()` in R code
- avoid assuming writable directories; use `tempdir()` in tests/vignettes
