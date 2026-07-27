current_contract_fixture <- function(out_dir) {
  list(
    latent_space = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE),
    obs = data.frame(group = factor(c("A", "B"))),
    var = data.frame(
      index = c("column_gene_1", "column_gene_2"),
      row.names = c("row_gene_1", "row_gene_2")
    ),
    gene_expression = matrix(c(1, 2, 3, 4), nrow = 2),
    X_umap_2d = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE),
    out_dir = out_dir,
    centroid_min_points = 1,
    dataset_id = "current-contract",
    dataset_name = "Current contract",
    force = TRUE
  )
}

test_that("shipped R repository surfaces contain only current contract language", {
  repository_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = TRUE
  )
  roots <- file.path(
    repository_root,
    c("R", "man", "vignettes", ".github")
  )
  files <- unlist(
    lapply(
      roots,
      list.files,
      recursive = TRUE,
      full.names = TRUE,
      all.files = TRUE
    ),
    use.names = FALSE
  )
  files <- c(
    files,
    file.path(
      repository_root,
      c(
        "README.md",
        "DESCRIPTION",
        "NAMESPACE",
        "_pkgdown.yml",
        "publishing.md"
      )
    )
  )
  text_suffix <- "\\.(R|Rd|Rmd|md|ya?ml)$"
  files <- files[
    file.exists(files) &
      (!grepl("\\.", basename(files), fixed = TRUE) |
        grepl(text_suffix, files, ignore.case = TRUE))
  ]
  retired_language <- paste0(
    "\\bfallback\\b|fall(?:ing)?[[:space:]-]+back|",
    "best[[:space:]_-]*effort|\\blegacy\\b|\\bdeprecated\\b|",
    "deprecation|backwards?[[:space:]_-]*compat|\\bshim\\b"
  )
  violations <- character()
  for (file in files) {
    lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
    matches <- grep(retired_language, lines, ignore.case = TRUE, perl = TRUE)
    if (length(matches) > 0L) {
      relative <- substring(
        normalizePath(file, mustWork = TRUE),
        nchar(repository_root) + 2L
      )
      violations <- c(
        violations,
        sprintf("%s:%d: %s", relative, matches, trimws(lines[matches]))
      )
    }
  }
  expect_identical(violations, character())

  publishing_path <- file.path(repository_root, "publishing.md")
  if (file.exists(publishing_path)) {
    publishing <- paste(
      readLines(
        publishing_path,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      collapse = "\n"
    )
    expect_false(
      grepl(
        "alias[[:space:]-]+metapackage|cellucid-r[[:space:]]*\\(alias\\)",
        publishing,
        ignore.case = TRUE,
        perl = TRUE
      )
    )
  }
})

test_that("README cross-references the current Cellucid ecosystem", {
  repository_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = TRUE
  )
  readme_candidates <- c(
    file.path(repository_root, "README.md"),
    file.path(repository_root, "00_pkg_src", "cellucid", "README.md")
  )
  readme_path <- readme_candidates[file.exists(readme_candidates)]
  expect_length(readme_path, 1L)
  readme <- paste(
    readLines(
      readme_path,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    collapse = "\n"
  )
  required_urls <- c(
    "https://cellucid.readthedocs.io/en/latest/user_guide/r_package/index.html",
    "https://github.com/theislab/cellucid",
    "https://github.com/theislab/cellucid-python",
    "https://github.com/theislab/cellucid-datasets",
    "https://github.com/theislab/cellucid-demo-custom-datasets",
    "https://github.com/theislab/cellucid-annotation"
  )
  expect_true(all(vapply(required_urls, grepl, logical(1), x = readme, fixed = TRUE)))
  expect_match(readme, "currently distributed from its GitHub repository", fixed = TRUE)
  expect_false(grepl("CRAN.R-project.org/package=cellucid", readme, fixed = TRUE))

  build_ignore_path <- file.path(repository_root, ".Rbuildignore")
  if (file.exists(build_ignore_path)) {
    build_ignore <- readLines(
      build_ignore_path,
      warn = FALSE,
      encoding = "UTF-8"
    )
    expect_false(any(build_ignore == "^README\\.md$"))
  }
})

test_that("main is the only workflow branch", {
  repository_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = TRUE
  )
  workflow_directory <- file.path(repository_root, ".github", "workflows")
  workflow_paths <- file.path(
    workflow_directory,
    c("R-CMD-check.yaml", "pkgdown.yaml")
  )
  workflow_presence <- file.exists(workflow_paths)
  expect_true(
    all(workflow_presence) ||
      (!any(workflow_presence) && !dir.exists(workflow_directory))
  )

  workflow_branches <- function(path) {
    lines <- readLines(
      path,
      warn = FALSE,
      encoding = "UTF-8"
    )
    trimws(lines[grepl("^[[:space:]]*branches:", lines)])
  }

  if (all(workflow_presence)) {
    expect_identical(
      workflow_branches(workflow_paths[[1L]]),
      rep("branches: [main]", 2L)
    )
    expect_identical(
      workflow_branches(workflow_paths[[2L]]),
      "branches: [main]"
    )
  }
})

test_that("release metadata and artifacts have one exact 0.9.1 identity", {
  repository_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = FALSE
  )
  description_path <- file.path(repository_root, "DESCRIPTION")
  version <- if (file.exists(description_path)) {
    unname(read.dcf(description_path)[1L, "Version"])
  } else {
    as.character(utils::packageVersion("cellucid"))
  }
  expect_identical(version, "0.9.1")

  news_path <- file.path(repository_root, "NEWS.md")
  if (file.exists(news_path)) {
    news <- readLines(news_path, warn = FALSE, encoding = "UTF-8")
    expect_identical(news[[1L]], "# cellucid 0.9.1 <!-- CELLUCID_VERSION -->")
  }

  citation_path <- file.path(repository_root, "CITATION.cff")
  if (file.exists(citation_path)) {
    citation <- readLines(citation_path, warn = FALSE, encoding = "UTF-8")
    expect_true(any(citation == "version: 0.9.1  # CELLUCID_VERSION"))
  }

  release_path <- file.path(
    repository_root,
    ".github",
    "workflows",
    "release.yaml"
  )
  if (file.exists(release_path)) {
    release <- paste(
      readLines(release_path, warn = FALSE, encoding = "UTF-8"),
      collapse = "\n"
    )
    expect_match(release, 'EXPECTED_TAG="v\\$\\{VERSION\\}"')
    expect_match(release, 'test "\\$\\{GITHUB_REF_TYPE\\}" = "tag"')
    expect_match(release, 'test "\\$\\{GITHUB_REF_NAME\\}" = "\\$\\{EXPECTED_TAG\\}"')
    expect_match(
      release,
      'git merge-base --is-ancestor "\\$\\{GITHUB_SHA\\}" origin/main'
    )
    expect_match(release, 'sha256sum "\\$\\{TARBALL\\}"')
    expect_lt(
      regexpr("R CMD check --as-cran", release, fixed = TRUE)[[1L]],
      regexpr("Attach tarball to GitHub Release", release, fixed = TRUE)[[1L]]
    )
  }
})

test_that("categorical code width is one required exact caller choice", {
  expect_true(
    identical(
      formals(cellucid_prepare)[["obs_categorical_dtype"]],
      quote(expr = )
    )
  )

  missing_out <- tempfile("cellucid_r_missing_categorical_dtype_")
  expect_error(
    do.call(cellucid_prepare, current_contract_fixture(missing_out)),
    "obs_categorical_dtype"
  )
  expect_false(dir.exists(missing_out))

  automatic_out <- tempfile("cellucid_r_automatic_categorical_dtype_")
  automatic_args <- current_contract_fixture(automatic_out)
  automatic_args$obs_categorical_dtype <- "auto"
  expect_error(
    do.call(cellucid_prepare, automatic_args),
    "exactly one of \"uint8\" or \"uint16\""
  )
  expect_false(dir.exists(automatic_out))
})

test_that("NULL alone selects var row names and every string names one exact column", {
  gene_id_default <- formals(cellucid_prepare)$var_gene_id_column
  expect_null(gene_id_default)

  rowname_out <- tempfile("cellucid_r_rowname_gene_ids_")
  rowname_args <- current_contract_fixture(rowname_out)
  rowname_args$obs_categorical_dtype <- "uint16"
  rowname_args$var_gene_id_column <- NULL
  do.call(cellucid_prepare, rowname_args)
  expect_true(file.exists(file.path(rowname_out, "var", "row_gene_1.values.f32")))
  expect_false(file.exists(file.path(rowname_out, "var", "column_gene_1.values.f32")))

  column_out <- tempfile("cellucid_r_literal_index_column_")
  column_args <- current_contract_fixture(column_out)
  column_args$obs_categorical_dtype <- "uint16"
  column_args$var_gene_id_column <- "index"
  do.call(cellucid_prepare, column_args)
  expect_true(file.exists(file.path(column_out, "var", "column_gene_1.values.f32")))
  expect_false(file.exists(file.path(column_out, "var", "row_gene_1.values.f32")))

  missing_out <- tempfile("cellucid_r_missing_gene_id_column_")
  missing_args <- current_contract_fixture(missing_out)
  missing_args$obs_categorical_dtype <- "uint16"
  missing_args$var_gene_id_column <- "missing"
  expect_error(
    do.call(cellucid_prepare, missing_args),
    "var_gene_id_column 'missing' not found in var"
  )
  expect_false(dir.exists(missing_out))
})
