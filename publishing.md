# Publishing `cellucid` (R package)

A beginner-friendly guide to publishing the R package `cellucid` to CRAN and beyond.

**Everything is automated via GitHub Actions** — you don't need R installed locally.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [One-Time Setup](#one-time-setup)
3. [Releasing a New Version](#releasing-a-new-version)
4. [Publishing to CRAN](#publishing-to-cran)
5. [Publishing to r-universe](#publishing-to-r-universe)
6. [Publishing to conda](#publishing-to-conda)
7. [Troubleshooting](#troubleshooting)
8. [Reference](#reference)

---

## Quick Start

**To release a new version:**

1. **Bump version**: Update `DESCRIPTION`, `CITATION.cff`, and `NEWS.md`
2. **Merge to main**: Create a PR, wait for CI, merge
3. **Create release**: Releases → "Draft new release" → tag `vX.Y.Z` → publish
4. **Download tarball**: Actions → "Build Release Tarball" → download `r-source-tarball`
5. **Submit to CRAN**: https://cran.r-project.org/submit.html → upload tarball → confirm email

That's it! Keep reading for detailed explanations.

---

## One-Time Setup

Do these steps once when setting up the repo.

### 1. Enable GitHub Actions

Go to: **Settings → Actions → General**

- Ensure "Allow all actions" is selected
- Under "Workflow permissions", select "Read and write permissions"

### 2. Enable GitHub Pages (for documentation website)

Go to: **Settings → Pages**

- Source: **GitHub Actions**

Your documentation will appear at: `https://theislab.github.io/cellucid-r/`

### 3. Verify workflows exist

Check that `.github/workflows/` contains:

| File | Purpose |
|------|---------|
| `R-CMD-check.yaml` | Runs tests on every push/PR (Windows, macOS, Linux) |
| `release.yaml` | Builds the CRAN-ready tarball when you create a release |
| `pkgdown.yaml` | Builds and deploys the documentation website |

---

## Releasing a New Version

### Step 1: Bump the version

Update these files with the new version number:

| File | What to change |
|------|----------------|
| `DESCRIPTION` | `Version:` field (e.g., `Version: 0.10.0`) |
| `CITATION.cff` | `version:` field (e.g., `version: 0.10.0`) |
| `NEWS.md` | Add a new section header (e.g., `# cellucid 0.10.0`) |

### Step 2: Create a PR and merge

1. Create a branch with your version bump
2. Push and create a PR
3. Wait for CI checks to pass (green checkmark)
4. Merge to `main`

### Step 3: Create a GitHub Release

1. Go to **Releases** → **"Draft a new release"**
2. Click **"Choose a tag"** → type `vX.Y.Z` (e.g., `v0.10.0`) → **"Create new tag"**
3. Target: `main`
4. Title: `cellucid X.Y.Z`
5. Description: copy from `NEWS.md` or write release notes
6. Click **"Publish release"**

This automatically triggers the "Build Release Tarball" workflow.
The workflow requires the tag to equal `v` plus the `DESCRIPTION` version and
requires the tagged commit to belong to `main`.

### Step 4: Download the tarball

1. Go to **Actions** → find the **"Build Release Tarball"** run
2. Scroll down to **"Artifacts"**
3. Download **`r-source-tarball`**
4. Unzip it — you'll get `cellucid_X.Y.Z.tar.gz` and its `.sha256` checksum

The tarball is also attached directly to your GitHub Release.

### Step 5: Next steps

- **For CRAN submission**: See [Publishing to CRAN](#publishing-to-cran) below
- **For r-universe**: It auto-updates from GitHub (after one-time setup)
- **For users**: They can install immediately:
  ```r
  remotes::install_github("theislab/cellucid-r")
  library(cellucid)
  ```

---

## Publishing to CRAN

CRAN is the main R package repository. Once accepted, users can install with just `install.packages("cellucid")`.

### Before Your First Submission

**Read these (they're short but important):**
- [CRAN policies](https://cran.r-project.org/web/packages/policies.html)
- [Submission checklist](https://cran.r-project.org/web/packages/submission_checklist.html)

**Your package must have:**
- 0 errors, 0 warnings on `R CMD check --as-cran`
- Working examples for all exported functions
- All URLs must work (no 404s, no redirects)

**The CI already checks most of this** — if GitHub Actions is green, you're likely ready.

### How to Submit

#### 1. Go to the submission form

Open: **https://cran.r-project.org/submit.html**

#### 2. Fill out the form

| Field | Value |
|-------|-------|
| Package tar.gz file | Upload `cellucid_X.Y.Z.tar.gz` |
| Package name | `cellucid` |
| Package version | `X.Y.Z` |
| Maintainer name | Your name |
| Maintainer email | Must match email in DESCRIPTION |
| Optional comment | Paste contents of `cran-comments.md` |

#### 3. Submit and confirm via email

1. Click **"Upload package"**
2. **Check your email immediately** (within minutes)
3. **Click the confirmation link** — this is required!

#### 4. Wait for results

- **Automatic checks**: 30 minutes to a few hours
- **Manual review** (first submission): A few days to 2 weeks
- **Track your submission**: [CRAN incoming dashboard](https://r-hub.github.io/cransays/articles/dashboard.html)

You'll receive an email when accepted or if changes are needed.

### If CRAN Asks for Changes

Don't panic — most first submissions get feedback.

1. **Fix the issues** in your code
2. **Bump the version** (e.g., `0.9.0` → `0.9.1`)
3. **Update `cran-comments.md`** with a "Resubmission" section:
   ```markdown
   ## Resubmission

   This is a resubmission. Changes made:
   - Fixed the URL that was returning 404
   - Reduced example runtime to under 5 seconds

   ## R CMD check results
   ...
   ```
4. **Create a new release** and **resubmit**

### CRAN Tips

| Rule | Details |
|------|---------|
| Submission frequency | Max once per month (unless responding to feedback) |
| Example runtime | Total < 5 seconds |
| Documentation | All exported functions need `@returns` tag |
| Title | < 65 characters, no period at end |
| Description | Must end with a period |

---

## Publishing to r-universe

[r-universe](https://r-universe.dev) builds your package automatically from GitHub and provides pre-built binaries (no compilation needed for users).

### Option A: Manual Setup (recommended — works before CRAN)

1. **Create a new GitHub repo** named exactly: `theislab.r-universe.dev`
   - Go to: https://github.com/new
   - Make it **public**

2. **Add `packages.json`** to the repo root:
   ```json
   [
     {
       "package": "cellucid",
       "url": "https://github.com/theislab/cellucid-r"
     }
   ]
   ```

3. **Install the r-universe GitHub App**:
   - Go to: https://github.com/apps/r-universe/installations/new
   - Select `theislab`
   - Grant access to all repositories

4. **Wait ~1 hour** for the first build

Your packages will appear at: https://theislab.r-universe.dev

### Option B: Auto-detection (after CRAN only)

Once your package is on CRAN with a GitHub URL in DESCRIPTION, r-universe can auto-detect it:

1. Just install the r-universe GitHub App (same link as above)
2. Wait up to 24 hours for auto-detection

### How Users Install from r-universe

```r
install.packages("cellucid", repos = "https://theislab.r-universe.dev")
```

---

## Publishing to conda

Goal: `mamba install -c conda-forge r-cellucid`

### Recommended: Get on CRAN first

conda-forge prefers to build R packages from CRAN tarballs. Once on CRAN:

1. Fork https://github.com/conda-forge/staged-recipes
2. Add a recipe under `recipes/r-cellucid/`
3. Submit a PR

conda-forge reviewers will help with the recipe format.

The conda-forge package name is `r-cellucid`, following conda-forge's
required `r-` prefix for R packages.

---

## Troubleshooting

### Where are the workflow artifacts?

1. Go to **Actions**
2. Click on the workflow run
3. Scroll to the bottom — **Artifacts** section
4. Download `r-source-tarball`

### CI fails with "pdflatex not available"

The release workflow installs TinyTeX and the declared LaTeX packages before building the vignette.

### CRAN rejected my package

Common issues and fixes:

| Issue | Fix |
|-------|-----|
| Examples take too long | Wrap slow examples in `\donttest{}` |
| Missing `@returns` | Add return value documentation to all exported functions |
| URLs don't work | Run `urlchecker::url_check()` locally and fix broken URLs |
| Title too long | Keep under 65 characters |
| Description doesn't end with period | Add a period |

### Need to update just the citation?

Edit `CITATION.cff` directly and push. For changes that need a new version, follow the full release process.

---

## Reference

### Package Naming

| Context | Name |
|---------|------|
| R package name | `cellucid` |
| GitHub repo | `cellucid-r` |
| What users type in R | `library(cellucid)` |
| CRAN | `cellucid` |
| conda-forge | `r-cellucid` |

R package names cannot contain `-`, so the package itself must stay `cellucid`.

### Files to Update for Each Release

| File | What changes |
|------|--------------|
| `DESCRIPTION` | `Version:` field |
| `CITATION.cff` | `version:` field |
| `NEWS.md` | New section header added |

### GitHub Actions Workflows

| Workflow | Trigger | Purpose |
|----------|---------|---------|
| R-CMD-check | Push, PR | Tests on Windows/macOS/Linux |
| Build Release Tarball | Release, manual | Creates CRAN-ready `.tar.gz` |
| pkgdown | Push to main, release | Deploys documentation website |

### Useful Links

| Resource | URL |
|----------|-----|
| Submit to CRAN | https://cran.r-project.org/submit.html |
| CRAN policies | https://cran.r-project.org/web/packages/policies.html |
| CRAN submission checklist | https://cran.r-project.org/web/packages/submission_checklist.html |
| Track CRAN submission | https://r-hub.github.io/cransays/articles/dashboard.html |
| r-universe setup | https://docs.r-universe.dev/publish/set-up.html |
| R Packages book | https://r-pkgs.org/ |
