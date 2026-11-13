# Contributing to `cellucid` (R)

Contributions are welcome — code, docs, bug reports, and reproducible examples all help.

This file focuses on `cellucid-r` (the R package that writes Cellucid-ready AnnData stores for the Cellucid viewer), and also helps route issues/PRs to the correct repo in the Cellucid ecosystem.

By participating, you agree to follow the project’s Code of Conduct:
- `CODE_OF_CONDUCT.md`

If you’re reporting a security issue, please follow:
- `SECURITY.md`

---

## Which repo should I contribute to?

Cellucid is split by responsibility:

| Repo | What it is | Contribute here when you… |
|---|---|---|
| `cellucid` | Web app (UI + state + WebGL rendering) | are fixing UI bugs, rendering/performance, figure export, sessions, or community annotation frontend |
| `cellucid-python` | Python package + CLI (`prepare`, `serve`, `show_anndata`, hooks) + Sphinx docs | are fixing Python/CLI bugs, data preparation, server endpoints, Jupyter hooks, or docs on ReadTheDocs |
| `cellucid-r` (this repo) | R package that writes Cellucid-ready AnnData stores | are changing the R writer (`cellucid_prepare()`), adding R-side tests/docs, or preparing for CRAN |
| `cellucid-annotation` | GitHub repo template for community annotation | are changing the repo schema/validation/workflows |
| `cellucid-datasets` | The published demo catalog | are correcting a published generation or its catalog entry |
| `cellucid-demo-custom-datasets` | Worked example of publishing your own datasets | are changing that guide or its synthetic examples |

If you’re not sure where a bug belongs, open an issue in the repo you’re currently using and include:
- how you loaded data (a local `.zarr`/`.zarr.zip`/`.h5ad` file, a link, a server, or Jupyter),
- the UI environment (hosted app vs local app vs Jupyter iframe),
- and the smallest reproduction you can share.

---

## Fast paths (pick your contribution type)

### I want to report a bug

Please include:
- R version (`R.version.string`)
- OS (macOS/Windows/Linux)
- how you wrote the store (which function call, which inputs)
- expected vs actual behavior
- the smallest dataset you can share (prefer synthetic/reproducible)
- any `R CMD check` output if relevant

If the bug is “the viewer looks wrong”, also include:
- which viewer environment you used (hosted app vs local app vs Jupyter iframe)
- browser + version
- whether you opened the store with the browser file picker, from a link, or from a server

### I want to contribute docs only

Docs live in:
- `cellucid-r/man/` (hand-written `.Rd` files)
- `cellucid-r/vignettes/` (package vignette)
- `cellucid-r/README.md`

Edit the `.Rd` files directly. This package runs no documentation generator:
`man/` and `NAMESPACE` are hand-written and authoritative, and `R/` holds no
roxygen comments. See “Documentation (help page + vignettes)” below.

### I want to add/modify R code

Fast workflow:
1) Set up your R dev environment (see “Development setup”)
2) Make a small, focused change
3) Add/adjust tests (`testthat`)
4) Run `devtools::check()`
5) Submit a PR with a clear “what/why/how to verify”

---

## Development setup (cellucid-r)

### Prerequisites

- R `>= 4.3.0` (matches `DESCRIPTION`)
- Git
- A C toolchain for package builds: Rtools on Windows, the Xcode Command Line
  Tools on macOS, or a compiler plus R development headers on Linux

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
install.packages(c("devtools", "testthat"))
```

For vignette builds:

```r
install.packages(c("knitr", "rmarkdown"))
```

Optional (useful during release prep):

```r
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

Run the core suite after every change (eight files that cover each part of
the package once, about a minute on an idle machine; the list is in `tools/test-core.R`, and
`devtools::test()` stays the full run):

```bash
Rscript tools/test-core.R
```

Run a full check (recommended before PRs):

```r
devtools::check()
```

Guidelines:
- Add tests when behavior changes (especially edge cases like missing values, mismatched dimensions, sparse matrices).
- Prefer small synthetic inputs; avoid committing real datasets.

---

## Documentation (help page + vignettes)

### Edit `.Rd` files by hand

The package ships its help page(s) under `cellucid-r/man/`, and they are
hand-written. Do not run `devtools::document()` or `roxygen2::roxygenise()`
here — this package has no roxygen source. Two things depend on that:

- `NAMESPACE` carries `useDynLib(cellucid, .registration = TRUE, .fixes = "C_")`,
  which no generator would emit, so a regenerated `NAMESPACE` would break the
  `.Call()` into the CRC-32 routine of the `.zarr.zip` writer.
- `man/cellucid_prepare.Rd` carries the runnable `\examples{}` block and the
  store-layout and orientation notes that the help page needs.

If you change the public API (especially `cellucid_prepare()`), update
`cellucid-r/man/cellucid_prepare.Rd` in the same commit.
`tests/testthat/test-current-contract.R` reads the `\usage` block back out of
the `.Rd` file and compares it to the real formals, so the help page cannot
silently fall out of step with the code if an argument name, its position, or
its default value changes. `R CMD check` runs the same comparison through
`tools::codoc()`.

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
- **format-first** (every store is a standard AnnData store that Python's
  `anndata` reads without a warning and the web app opens directly)
- **CRAN-ready** (checks, documentation, and package structure matter)

If you propose adding a new dependency:
- prefer `Suggests` over `Imports` unless strictly required
- keep the core writer usable on minimal installations

If you change what a store holds (its layout or `uns['cellucid']`):
- coordinate with `cellucid-python` and the web app (`cellucid`) so every repository implements the same exact current contract
- add tests that validate the new behavior, and keep
  `tests/testthat/test-anndata-interop.R` passing: it reads the store back with
  Python's `anndata` (set `CELLUCID_ANNDATA_PYTHON` to a Python with `anndata`
  and `zarr`, or have one as `python3`)

---

## Comments say what the code does; history goes to the decision records

A comment says what the code does and why it must be so: the invariant it
keeps, the failure it prevents, the unit of a number, the function that owns
the other half. It is written in the present tense and stays true for as long
as the code below it is unchanged. This holds for `R/`, `src/`, `tools/` and
`tests/`. A help page under `man/` holds no history either: it is what a user
reads.

A comment never holds:

1. a CI run number or address;
2. a finding, audit or plan id (`F123`, `L25-3`, `R9W-5`, "round 9", "bug
   hunt", a section of a plan), or a path into a repository the reader cannot
   open;
3. a date of our own work (a date of an outside version the code depends on,
   such as a protocol revision, an R or a Bioconductor release, stays);
4. who decided, or their words;
5. what the code used to do, what was tried, what was found or fixed;
6. numbers measured on a machine, with the machine, engine or dataset.

Each of these goes, in the same change, to the page that records decisions
for that part of the package. Those pages are on the documentation site,
which is built from `cellucid-python/docs/`:

| Code | Where its history goes |
|---|---|
| The served store, the viewer channel and the notebook widget (`viewer-server.R`, `viewer-widget.R`, `zarr-store.R`) | `cellucid-python/docs/user_guide/web_app/p_developer_docs/26_decisions_served_store_and_widgets.md` |
| A measurement of the R package alone (matrices held outside R's memory, the writer) | The R architecture page, `cellucid-python/docs/user_guide/r_package/h_developer_docs/01_codebase_architecture.md`, under "Measurements" (the first one starts that section) |
| What CI showed that a local run does not | The R testing page, `cellucid-python/docs/user_guide/r_package/h_developer_docs/02_testing_and_ci.md`, "What CI sees that a local run does not" |
| A test that pins a behavior of the web app | The decision record of that subsystem; `cellucid-python/docs/user_guide/web_app/p_developer_docs/index.md` lists them |

In a decision record, a decision goes to "Decisions and why", numbers to
"Measurements that justify it", an attempt to "Tried and reverted". The
record tells the story in plain words, with its date and its numbers; it
holds no finding, audit or plan id either. Where a constant or a design rests
on one, the comment keeps the reason in one clause and names the page:
`(measured: 26_decisions_served_store_and_widgets.md)`.

What stays in a comment: the failure the code prevents, told as a condition;
pointers to tests, functions and pages that exist; facts about outside
versions.

Tests follow the same rule. A test's name says what it pins and holds no id.
A test file's header says what the file pins and what breaks without it; the
incident that led to it goes to the record.

`tests/testthat/test-comment-history.R` checks the comments of `R`, `tests`, `tools` and `src`, the comments of `.github/workflows` (in a source checkout) and the `test_that` descriptions for finding ids, rounds, CI run addresses and plan pointers.

**Before and after: a reason told through an incident**
(`R/backed-matrix.R`, `.delayed_untransposed()`).

```r
# Before
# ... up to Bioconductor 3.16 (DelayedArray 0.20 to 0.24, R 4.1
# and 4.2) a sparse block read through it subsets a SparseArraySeed with
# MatrixGenerics::rowAnyNAs(useNames = NA), which matrixStats 1.x, the one
# CRAN serves those R versions, refuses ("Argument 'useNames' must be either
# TRUE or FALSE"), so cellucid_prepare() of a cells x genes DelayedArray
# failed there (CI run 43106603588).

# After
# ... up to Bioconductor 3.16 (DelayedArray 0.20 to 0.24, R 4.1
# and 4.2) a sparse block read through it subsets a SparseArraySeed with
# MatrixGenerics::rowAnyNAs(useNames = NA), which matrixStats 1.x, the one
# CRAN serves those R versions, refuses ("Argument 'useNames' must be either
# TRUE or FALSE"), so cellucid_prepare() of a cells x genes DelayedArray
# would fail there.
```

The R testing page, "What CI sees that a local run does not", gains the run
and what failed in it. The facts about Bioconductor and matrixStats stay in
the comment: they are outside versions the code depends on.

**Before and after: a measurement that justifies a design**
(`R/viewer-widget.R`).

```r
# Before
# Bodies and the channel's binary frames go to the page as the binary buffers
# of a Jupyter message ... Measured on one 1 MiB body,
# base64 and its JSON took 17.6 ms of R's one thread (60 MB/s); its buffer,
# 0.3 ms.

# After
# Bodies and the channel's binary frames go to the page as the binary buffers
# of a Jupyter message, never as base64 inside its JSON, which costs R's one
# thread far more (measured: 26_decisions_served_store_and_widgets.md).
```

Record 26 already holds those numbers; a number it lacks is added there in
the same change.

---

## How to validate a store end-to-end

The most reliable validation is to write a tiny dataset and open it in the viewer.

Recommended workflow:
1) Run `cellucid_prepare()` on a tiny synthetic dataset with `format = "zarr.zip"`.
2) Open the Cellucid web app (hosted or local).
3) Open the `.zarr.zip` file from the data panel.
4) Confirm:
   - the correct number of cells render
   - categorical and continuous fields appear
   - gene expression search works (if written)

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

### `devtools::document()` rewrote `NAMESPACE` or `DESCRIPTION`

Cause:
- this package has no roxygen source, so running a documentation generator can
  only remove content. It drops `useDynLib(...)` from `NAMESPACE` and adds a
  `Config/roxygen2/version` field to `DESCRIPTION`.

Fix:
- `git checkout NAMESPACE DESCRIPTION man/` and edit `man/*.Rd` by hand instead.
- `tests/testthat/test-current-contract.R` fails while any of those generated
  changes are still present.

### `R CMD check` fails on vignettes

Common causes:
- missing suggested packages (`knitr`, `rmarkdown`)
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
