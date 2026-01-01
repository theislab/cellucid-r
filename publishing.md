# Publishing `cellucid` (R package) — very beginner-friendly

This guide explains how to publish the **R package** `cellucid` from the **GitHub repo**
`cellucid-r`, even if you’ve never published an R package before.

It focuses on a “CI-first” approach:
- GitHub Actions builds and checks the package
- you download the resulting `cellucid_<version>.tar.gz`
- you upload it to CRAN / reference it for Bioconductor / use it for conda packaging

## Naming (read this first)

You asked for:
- **R package name**: `cellucid`
- **conda package name only**: `cellucid-r`
- everywhere else: `cellucid`

That mapping is intentional and valid.

| Thing | Name | Where it appears |
|---|---|---|
| R package name | `cellucid` | `cellucid-r/DESCRIPTION` → `Package:` |
| What users type in R | `cellucid` | `library(cellucid)` |
| GitHub repo name | `cellucid-r` | `https://github.com/theislab/cellucid-r` |
| Source tarball name | `cellucid_<version>.tar.gz` | built by CI (`R CMD build`) |
| CRAN / Bioconductor | `cellucid` | indexed by R package name |
| conda | `cellucid-r` | `conda install ... cellucid-r` |

Notes:
- R package names cannot contain `-`, so the package itself must stay `cellucid`.
- conda package names can contain `-`, which is why the conda name can be `cellucid-r`.

## What you need (accounts + prerequisites)

### Accounts (what requires what)

- GitHub: required (to push changes + create releases)
- CRAN: no account required up front, but you must provide a working email and reply to CRAN emails
- Bioconductor: GitHub account required (submission happens via a GitHub issue)
- r-universe: GitHub account required (one-time registration)
- conda-forge: GitHub account required (submission happens via a PR)

### Local tooling (optional, but helpful)

You can ship releases using only GitHub Actions (no local R install), but for day-to-day development
you’ll eventually want:
- R (this package depends on `R (>= 4.3.0)`; see `cellucid-r/DESCRIPTION`)
- (Optional) RStudio
- Build tools (Windows: Rtools; macOS: Xcode Command Line Tools; Linux: build-essential)

If you don’t have these yet, you can still do the “CI-first” release flow below.

## Glossary (2-minute)

- **Tag**: a git reference like `v0.99.0` pointing at a commit <!-- CELLUCID_VERSION -->
- **GitHub Release**: a GitHub UI entry for a tag (with release notes + attached files)
- **Workflow**: a GitHub Actions automation (CI)
- **Artifact**: a file produced by a workflow that you can download (here: the tarball)
- **Source tarball**: `cellucid_<version>.tar.gz` — this is what CRAN expects you to upload

## How publishing works in this repo

GitHub Actions can reliably automate:
- running `R CMD check` on multiple OSes / R versions
- running `BiocCheck::BiocCheck()` for Bioconductor readiness
- building a reproducible source tarball (`cellucid_<version>.tar.gz`)
- building + publishing the docs website (pkgdown) to GitHub Pages

GitHub Actions cannot (by design) fully automate:
- CRAN submission (you must upload via a web form + respond to email)
- Bioconductor submission (you must open an issue in the Bioconductor Contributions tracker)
- conda-forge onboarding (you must open PR(s) to `conda-forge/staged-recipes`)

## Quick start: “I only want a correct tarball I can submit”

1. Merge a PR to `main` that bumps the version and updates `NEWS.md`.
2. Create a GitHub Release with tag `vX.Y.Z`.
3. Download `cellucid_<version>.tar.gz` from:
   - GitHub → Actions → `Build Release Tarball` → workflow run → Artifacts → `r-source-tarball`
4. Upload that tarball to CRAN, or reference it for Bioconductor/conda packaging.

If you want the full, beginner-friendly detail, keep reading.

---

## One-time setup (do once per GitHub repo)

### 1) Confirm workflows exist

In `cellucid-r/.github/workflows/` you should have:
- `R-CMD-check.yaml` (main CI gate)
- `bioccheck.yaml` (Bioconductor policy gate)
- `release.yaml` (workflow name: “Build Release Tarball”)
- `pkgdown.yaml` (website)

### 2) Enable GitHub Actions

In GitHub:
- Settings → Actions → General → ensure Actions are enabled

### 3) Enable GitHub Pages (for pkgdown)

This repo deploys a generated website to the `gh-pages` branch (root).

In GitHub:
- Settings → Pages
  - Source: “Deploy from a branch”
  - Branch: `gh-pages` / `(root)`

After the next successful `pkgdown` workflow run, your site should appear at:
- `https://theislab.github.io/cellucid-r/` (adjust org/repo if you fork)

---

## Every release (detailed checklist)

### Step 0 — Decide where you’re publishing

Pick your targets up front:
- **GitHub only**: easiest; good for early alpha
- **r-universe**: easiest “binary install” path for users
- **CRAN**: broad R audience; stricter checks and policy expectations
- **Bioconductor**: best fit for bioinformatics tooling; additional policies
- **conda**: for conda env installs; we want the conda name to be `cellucid-r`

A common, low-friction order is:
GitHub Release → r-universe (optional) → CRAN or Bioconductor → conda packaging

### Step 1 — Update the version number (and keep it consistent)

Update these files:
- `cellucid-r/DESCRIPTION` → `Version:`
- `cellucid-r/CITATION.cff` → `version:` (keep in sync)
- `cellucid-r/NEWS.md` → add a new section for the version

Versioning tips:
- Versions must be monotonically increasing. Never re-use a version number.
- Git tags/releases use a `v` prefix: tag `v0.99.0` corresponds to `Version: 0.99.0`. <!-- CELLUCID_VERSION -->

Bioconductor note (important for first-time submitters):
- New Bioconductor submissions commonly start at `0.99.0` (not `0.0.x`).
- If you plan to submit to Bioconductor soon, consider making a “Bioconductor submission release”
  where you bump `Version:` to `0.99.0` and tag `v0.99.0`.

### Step 2 — Update docs + generated files (optional, requires local R)

If you have R installed locally, it’s worth doing a local check before releasing.

From within `cellucid-r/`:

1) Run unit tests:

```r
install.packages("devtools")
devtools::test()
```

2) Run a full check:

```r
devtools::check()
```

If you don’t have R locally, skip this and rely on CI (GitHub Actions).

### Step 3 — Open a release PR and merge it

Make a PR containing:
- version bump(s)
- `NEWS.md` updates
- any code changes

Wait for CI to be green:
- `R-CMD-check` (required)
- `BiocCheck (Bioconductor release)` (required if you care about Bioconductor)
- (Recommended before Bioconductor submission) run `BiocCheck` manually with `run_devel=true`

Then merge to `main`.

### Step 4 — Create a GitHub Release (this triggers the tarball build)

In GitHub:
1. Go to Releases → “Draft a new release”
2. Tag: create `vX.Y.Z` (must match `cellucid-r/DESCRIPTION` version `X.Y.Z`)
3. Target: `main`
4. Title: `cellucid X.Y.Z`
5. Release notes: paste the matching section from `cellucid-r/NEWS.md`
6. Click “Publish release”

This triggers the `Build Release Tarball` workflow (defined in `cellucid-r/.github/workflows/release.yaml`).

Alternative (no Release yet): you can build a tarball via the Actions UI:
- Actions → `Build Release Tarball` → “Run workflow”
  - keep `run_as_cran=true` (recommended if CRAN is a target)

### Step 5 — Download the tarball artifact

From the successful `Build Release Tarball` workflow run:
- Download the artifact named `r-source-tarball`
  - it contains `cellucid_<version>.tar.gz`

If you created a GitHub Release, the same tarball is also attached to that Release.

### Step 6 — (Optional) sanity-check the tarball

If you have R installed locally, you can verify the exact tarball:

```sh
R CMD check --as-cran cellucid_<version>.tar.gz
```

If you don’t have R locally, you can still at least confirm the tarball unzips cleanly.

---

## Publishing targets (step-by-step)

### 1) GitHub (source of truth)

You’re effectively “published” as soon as:
- the repo is public, and
- you create a GitHub Release.

Install from GitHub (users still load as `cellucid`):

```r
remotes::install_github("theislab/cellucid-r")
library(cellucid)
```

### 2) Website (pkgdown → GitHub Pages)

The `pkgdown` workflow builds the site into `docs/site/` and deploys that folder to the `gh-pages` branch.

Once GitHub Pages is enabled (see one-time setup), the site will update on:
- pushes to `main`, and
- GitHub Releases.

### 3) r-universe (easy binaries)

`r-universe` builds your package automatically from GitHub and hosts binaries for Windows/macOS/Linux.

One-time setup:
1. Create a GitHub repo named `universe` under your org/user, e.g. `theislab/universe`.
2. Add a `packages.json` file to that repo:
   ```json
   [
     {
       "package": "cellucid",
       "url": "https://github.com/theislab/cellucid-r"
     }
   ]
   ```
3. Register at: https://r-universe.dev/add

Users can then install via:
```r
install.packages("cellucid", repos = "https://theislab.r-universe.dev")
```

### 4) CRAN (manual submission)

CRAN submission is manual, but CI can do the hard parts:
- build the `.tar.gz`
- run `R CMD check --as-cran`

What to expect (first time):
- you upload a `.tar.gz` via a web form
- you reply to a confirmation email
- CRAN may email you review requests
- if you need changes, you must submit a new version (never re-submit the same version)

CRAN checklist:
1. Create a GitHub Release (or run the `Build Release Tarball` workflow manually).
2. Confirm the workflow log shows `R CMD check --as-cran` passed with:
   - 0 ERROR
   - 0 WARNING
   - ideally 0 NOTE
3. Download `cellucid_<version>.tar.gz` from the workflow artifact.
4. Submit at: https://cran.r-project.org/submit.html
   - upload the `.tar.gz`
   - reply to the confirmation email

CRAN note specific to this repo:
- CRAN does not recognize Bioconductor-only `DESCRIPTION` fields like `biocViews`, which can produce a NOTE.
- If you want a 0-NOTE CRAN submission, consider making a CRAN-specific branch/release that drops Bioconductor-only fields.

### 5) Bioconductor (manual submission; primary target for bio tooling)

Bioconductor submission is also manual, but CI makes it much easier.

Bioconductor checklist:
1. Ensure `R-CMD-check` is green.
2. Run `BiocCheck`:
   - Actions → `BiocCheck` → “Run workflow”
   - set `run_devel=true` (recommended before submitting)
3. Read the contributor guide:
   - https://contributions.bioconductor.org/
4. Open a new issue at:
   - https://github.com/Bioconductor/Contributions/issues/new
5. In your issue, provide:
   - package name: `cellucid`
   - repo URL: `https://github.com/theislab/cellucid-r`
   - brief package summary (2–5 sentences)
   - confirmation that `BiocCheck` passes (link the Actions run)
   - any reviewer notes you already addressed

Versioning note (important):
- New Bioconductor submissions commonly start at `0.99.0`. If you’re currently at `0.0.x`,
  plan a version bump for submission.

### 6) conda (publish as `cellucid-r`)

Goal:
- conda users install `cellucid-r`
- R users load `cellucid`

```sh
# If you publish `cellucid-r` on conda-forge:
mamba install -c conda-forge cellucid-r

# If you publish `cellucid-r` on your own channel:
mamba install -c <your-channel> cellucid-r

R -q -e 'library(cellucid); packageVersion("cellucid")'
```

#### Option A: conda-forge (widest reach) + alias package

conda-forge conventions for CRAN R packages use `r-<pkgname>`, so the “canonical” build would be:
- `r-cellucid` (contains the actual R package)

To keep your desired conda install name (`cellucid-r`), publish an additional tiny “metapackage”:
- `cellucid-r` (depends on `r-cellucid`)

Important:
- This option still publishes `r-cellucid` to conda-forge (that’s where the actual R code lives).
  `cellucid-r` is just a beginner-friendly alias you document.
- If you truly need **only** `cellucid-r` to exist in conda (no `r-cellucid` at all), use Option B.

Step-by-step:
1. Get `cellucid` onto CRAN first (recommended), so conda-forge can build from the CRAN tarball.
2. Fork https://github.com/conda-forge/staged-recipes
3. Add the main recipe under `recipes/r-cellucid/` (builds the real R package).
4. Add an alias recipe under `recipes/cellucid-r/` that depends on `r-cellucid`.
5. Submit a PR and follow conda-forge reviewer feedback.

Minimal `meta.yaml` idea for the alias recipe:
```yaml
package:
  name: cellucid-r
  version: 0.99.0  # CELLUCID_VERSION

build:
  number: 0
  noarch: generic

requirements:
  run:
    - r-cellucid ==0.99.0  # CELLUCID_VERSION
```

After merge, conda-forge creates feedstock repo(s) for updates.

#### Option B (strict naming): publish `cellucid-r` on your own conda channel

If you want **only** `cellucid-r` to exist in conda (no `r-cellucid` at all), publish to your own
channel (e.g., Anaconda.org) where you control naming.

High-level steps:
1. Create an account on Anaconda.org and an upload token.
2. Create a conda recipe that installs the R package `cellucid`, but sets:
   - `package.name: cellucid-r`
3. Build and upload the package.

Trade-offs:
- Pros: exactly one conda package name (`cellucid-r`)
- Cons: users must add your channel; less discoverable than conda-forge

---

## Troubleshooting (common)

### Where to click in GitHub to find logs + artifacts

- Actions tab: list of workflow runs
- Click a workflow run: shows jobs + logs
- Artifacts are listed at the bottom of the workflow run page

### “CRAN/Bioc asked me to change something”

- Make the requested fix on a PR.
- Merge to `main`.
- Bump the version.
- Cut a new GitHub Release.
- Re-submit with the new tarball / link updated CI.

---

## Quick reference: workflows

| Workflow | File | Purpose |
|---|---|---|
| `R-CMD-check` | `cellucid-r/.github/workflows/R-CMD-check.yaml` | cross-platform correctness gate |
| `BiocCheck` | `cellucid-r/.github/workflows/bioccheck.yaml` | Bioconductor policy gate |
| `Build Release Tarball` | `cellucid-r/.github/workflows/release.yaml` | builds `cellucid_<version>.tar.gz` and runs `--as-cran` |
| `pkgdown` | `cellucid-r/.github/workflows/pkgdown.yaml` | builds + deploys documentation website |
