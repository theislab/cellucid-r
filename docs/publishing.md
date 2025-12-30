# Publishing `cellucid` (beginner-friendly)

This guide is a step-by-step checklist for publishing the `cellucid` R package
via GitHub Actions, with minimal local setup.

It covers:
- GitHub releases
- r-universe
- CRAN
- Bioconductor
- conda-forge (after CRAN)

## What GitHub Actions can and cannot do

GitHub Actions can reliably automate:
- Running `R CMD check` on multiple OSes / R versions
- Running `BiocCheck::BiocCheck()` for Bioconductor readiness
- Building a source tarball (`.tar.gz`) you can submit to CRAN/Bioconductor
- Building and publishing the pkgdown documentation website

GitHub Actions cannot (by design) fully automate:
- **CRAN submission** (you must upload via the CRAN web form and respond to email)
- **Bioconductor submission** (you must open a GitHub issue in the Bioconductor Contributions tracker)
- **conda-forge onboarding** (you must open a PR to `staged-recipes`)

This guide shows how to make those manual steps “download artifact → upload”.

## One-time setup (do this once per repo)

### 1) Confirm workflows exist

In `.github/workflows/` you should have:
- `R-CMD-check.yaml` (main CI gate)
- `BiocCheck` workflow (`bioccheck.yaml`)
- `Build Release Tarball` workflow (`release.yaml`)
- `pkgdown.yaml` (website)
- `test-coverage.yaml` (optional Codecov)

What each workflow does:

| Workflow | When it runs | What you use it for |
|---|---|---|
| `R-CMD-check` | push/PR/manual | Main correctness gate across OS + R versions |
| `BiocCheck` | push/PR/manual | Bioconductor readiness gate (policy + structure) |
| `Build Release Tarball` | release/manual | Produces `cellucid_<version>.tar.gz` for CRAN/Bioc |
| `pkgdown` | push/release/manual | Builds and publishes your docs website |
| `test-coverage` | push/PR/manual | Coverage reporting (optional) |

### 2) Enable GitHub Actions

In your GitHub repo:
- **Settings → Actions → General**
  - Allow Actions to run (default is fine for most repos)

### 2b) Learn where to click (Actions + Artifacts)

You will use these two pages constantly:
- **Actions tab**: shows workflow runs and logs
- **A workflow run page**: contains the **Artifacts** section (downloadable `.zip` files)

For CRAN/Bioconductor submissions, you typically download the artifact from:
- Actions → `Build Release Tarball` → latest successful run → Artifacts → `r-source-tarball`

### 3) Enable GitHub Pages (for pkgdown)

This repo uses `pkgdown` + the `gh-pages` branch.

In your GitHub repo:
- **Settings → Pages**
  - Source: “Deploy from a branch”
  - Branch: `gh-pages` / `(root)`

After the next successful `pkgdown` workflow run, your site should appear at:
`https://<github-username>.github.io/<repo>/`

### 4) (Optional) Configure Codecov

If you want coverage reporting:
- Create a Codecov project for your repo.
- If the repo is private, add a `CODECOV_TOKEN` secret in:
  - **Settings → Secrets and variables → Actions**

If you don’t want Codecov, you can delete `test-coverage.yaml`.

## The release flow (recommended)

Use this flow for every release, even “small” ones.

### Step A — Prepare the release PR

1. Update `DESCRIPTION`:
   - Bump `Version:`
2. Update `NEWS.md`:
   - Add a bullet list of user-visible changes
3. Ensure documentation is up to date:
   - `README.md`
   - vignette(s) in `vignettes/`

### Step B — Merge, then create a GitHub Release

1. Merge your PR to `main`.
2. Create a GitHub Release:
   - **Releases → Draft a new release**
   - Tag: `vX.Y.Z` (example: `v0.99.1`)
   - Target: `main`
   - Title: `cellucid X.Y.Z`
   - Release notes: paste the matching section from `NEWS.md`

This triggers the `Build Release Tarball` workflow.

### Step C — Wait for CI to go green

Before publishing anywhere:
- Ensure `R-CMD-check` is green
- Ensure `BiocCheck (Bioconductor release)` is green
- (Optional) run `BiocCheck (devel)` from the Actions tab (recommended before Bioconductor submission)

### Step D — Download the tarball artifact

From the successful `Build Release Tarball` workflow run:
- Download the artifact named `r-source-tarball`
  - It contains `cellucid_<version>.tar.gz`

If you created a GitHub Release, the same tarball is also attached to the release.

That `.tar.gz` is what you submit to CRAN and reference for Bioconductor.

## Publishing targets

### 1) GitHub (source of truth)

You’re effectively “published” as soon as:
- The repo is public, and
- You create a GitHub Release

Recommended:
- Keep release notes in GitHub Releases and in `NEWS.md`.

### 2) r-universe (easy binary builds)

`r-universe` builds your package automatically from GitHub and hosts binaries for
Windows/macOS/Linux.

One-time setup:
1. Create a GitHub repo named `universe` under your org/user, e.g.:
   - `theislab/universe`
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

Result:
- Users can install via:
  ```r
  install.packages("cellucid", repos = "https://theislab.r-universe.dev")
  ```

### 3) CRAN

CRAN submission is manual, but GitHub Actions can do the hard parts:
- Build the `.tar.gz`
- Run `R CMD check --as-cran`

Checklist:
1. Create a GitHub Release (or run the workflow manually):
   - Actions → `Build Release Tarball`
2. Confirm the workflow log shows `R CMD check --as-cran` passed with:
   - 0 ERROR
   - 0 WARNING
   - ideally 0 NOTE (some NOTEs can be acceptable, but CRAN will ask)
3. Download `cellucid_<version>.tar.gz` from the workflow artifact.
4. Submit at: https://cran.r-project.org/submit.html
   - Upload the `.tar.gz`
   - Reply to the confirmation email

If CRAN emails you about check failures:
- Fix the issue on `main`
- Cut a new patch release
- Re-submit the new tarball

### 4) Bioconductor (primary target)

Bioconductor submission is also manual, but CI makes it much easier.

Checklist:
1. Ensure `R-CMD-check` is green.
2. Run `BiocCheck`:
   - Actions → `BiocCheck`
   - Run workflow, set `run_devel=true` (recommended before submission)
3. Read the contributor guide:
   - https://contributions.bioconductor.org/
4. Open a new issue at:
   - https://github.com/Bioconductor/Contributions/issues/new
5. Provide:
   - Repo URL
   - Brief package summary
   - Confirmation that `BiocCheck` passes (link the Actions run)
    - Any reviewer notes you already addressed (helps a lot)

### 5) conda-forge (after CRAN)

The easiest conda-forge path is after you’re on CRAN:
1. Fork: https://github.com/conda-forge/staged-recipes
2. Add a recipe under `recipes/r-cellucid/` that builds from the CRAN tarball.
3. Submit a PR.

After merge, conda-forge creates a feedstock repo for updates.

## Quick reference: what each workflow is for

- `R-CMD-check`: cross-platform correctness gate; keep it green always.
- `BiocCheck`: Bioconductor policy gate; run devel before submission.
- `Build Release Tarball`: produces the `.tar.gz` you submit to CRAN/Bioconductor.
- `pkgdown`: publishes the documentation website.
- `test-coverage`: optional coverage reporting.
