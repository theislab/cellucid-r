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
  expect_match(
    readme,
    "choose CRAN or GitHub based on current",
    fixed = TRUE
  )
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

test_that("DESCRIPTION is the one in-repo source of the 0.9.1 release identity", {
  repository_root <- cellucid_repository_root()
  description_path <- file.path(repository_root, "DESCRIPTION")
  version <- cellucid_description_version()
  expect_identical(version, "0.9.1")

  r_dependency <- NULL
  if (file.exists(description_path)) {
    description <- readLines(description_path, warn = FALSE, encoding = "UTF-8")
    version_line <- which(startsWith(description, "Version:"))
    expect_length(version_line, 1L)
    expect_identical(description[[version_line]], paste0("Version: ", version))
    # DESCRIPTION is DCF and cannot carry a `#` comment, so the CELLUCID_VERSION
    # sweep marker that every other version site holds inline is a Config field
    # pinned to the line immediately after Version:.
    expect_identical(
      description[[version_line + 1L]],
      "Config/cellucid/version-marker: CELLUCID_VERSION"
    )

    depends <- unname(read.dcf(description_path)[1L, "Depends"])
    r_dependency <- regmatches(
      depends,
      regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", depends)
    )
    expect_length(r_dependency, 1L)
  }

  news_path <- file.path(repository_root, "NEWS.md")
  if (file.exists(news_path)) {
    news <- readLines(news_path, warn = FALSE, encoding = "UTF-8")
    expect_identical(
      news[[1L]],
      sprintf("# cellucid %s <!-- CELLUCID_VERSION -->", version)
    )
    expect_false(any(grepl("unreleased", tolower(news), fixed = TRUE)))
    expect_true(
      any(news == sprintf("Version %s is the CRAN submission release.", version))
    )
  }

  citation_path <- file.path(repository_root, "CITATION.cff")
  if (file.exists(citation_path)) {
    citation <- readLines(citation_path, warn = FALSE, encoding = "UTF-8")
    expect_true(
      any(citation == sprintf("version: %s  # CELLUCID_VERSION", version))
    )
  }

  # README.md and the installation vignette are the two in-repo pages that
  # quote the version to users. Nothing outside this repository compares them
  # to DESCRIPTION, so every version string they contain is checked here.
  readme_path <- file.path(repository_root, "README.md")
  if (file.exists(readme_path)) {
    expect_identical(cellucid_semantic_versions(readme_path), version)
  }

  installation_candidates <- c(
    file.path(repository_root, "vignettes", "installation.Rmd"),
    system.file("doc", "installation.Rmd", package = "cellucid")
  )
  installation_path <- installation_candidates[
    nzchar(installation_candidates) & file.exists(installation_candidates)
  ][1L]
  if (!is.na(installation_path) && length(r_dependency) == 1L) {
    expect_setequal(
      cellucid_semantic_versions(installation_path),
      c(version, r_dependency)
    )
  }

  package_citation_path <- file.path(
    repository_root,
    "inst",
    "CITATION"
  )
  if (file.exists(package_citation_path)) {
    package_citation <- paste(
      readLines(
        package_citation_path,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      collapse = "\n"
    )
    expect_match(package_citation, "meta$Version", fixed = TRUE)
    expect_false(
      grepl(
        'packageVersion("cellucid")',
        package_citation,
        fixed = TRUE
      )
    )
    citation_environment <- new.env(parent = globalenv())
    citation_environment$meta <- list(Version = version)
    sys.source(package_citation_path, envir = citation_environment)
    expect_identical(citation_environment$package_version, version)
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

test_that("pkgdown has an exact current installation and navigation contract", {
  repository_root <- cellucid_repository_root()
  version <- cellucid_description_version()
  read_text <- function(path) {
    paste(
      readLines(path, warn = FALSE, encoding = "UTF-8"),
      collapse = "\n"
    )
  }

  installation_candidates <- c(
    file.path(repository_root, "vignettes", "installation.Rmd"),
    system.file("doc", "installation.Rmd", package = "cellucid")
  )
  installation_path <- installation_candidates[
    nzchar(installation_candidates) & file.exists(installation_candidates)
  ][1L]
  expect_length(installation_path, 1L)
  installation <- read_text(installation_path)
  expect_match(
    installation,
    sprintf(
      "**The active Cellucid for R source and documentation version is %s.**",
      version
    ),
    fixed = TRUE
  )
  expect_match(
    installation,
    'install.packages("cellucid")',
    fixed = TRUE
  )
  expect_match(
    installation,
    'remotes::install_github("theislab/cellucid-r")',
    fixed = TRUE
  )
  expect_match(
    installation,
    'packageVersion("cellucid")',
    fixed = TRUE
  )
  expect_match(
    installation,
    "is authoritative for registry availability",
    fixed = TRUE
  )

  readme_path <- file.path(repository_root, "README.md")
  if (file.exists(readme_path)) {
    readme <- read_text(readme_path)
    expect_match(
      readme,
      sprintf("**Active package version — %s**", version),
      fixed = TRUE
    )
    expect_match(
      readme,
      "articles/installation.html",
      fixed = TRUE
    )
  }

  pkgdown_path <- file.path(repository_root, "_pkgdown.yml")
  if (file.exists(pkgdown_path)) {
    pkgdown <- read_text(pkgdown_path)
    expect_match(
      pkgdown,
      "left: [home, installation, reference, articles, news]",
      fixed = TRUE
    )
    expect_match(
      pkgdown,
      'in_header: >-\n      <link rel="icon" type="image/svg+xml"',
      fixed = TRUE
    )
    expect_match(pkgdown, "right: [search, webapp]", fixed = TRUE)
    expect_match(pkgdown, "href: https://www.cellucid.com", fixed = TRUE)
    expect_false(grepl("development", pkgdown, ignore.case = TRUE))
  }
})

test_that("man/ and NAMESPACE are hand-written and match the code exactly", {
  # This is the codoc check R CMD check performs, run against the installed Rd
  # database: the \usage block is turned back into a function and compared to
  # the real formals, argument names, order, and default expressions included.
  # It is what keeps a hand-written help page from drifting away from the
  # signature, so no generator is needed to hold the two together.
  #
  # tools::codoc() itself is not usable here: it loads and unloads the package
  # it inspects, and cellucid's .onUnload() calls library.dynam.unload(), which
  # would pull the native export lock out from under the running suite.
  repository_root <- cellucid_repository_root()
  man_directory <- file.path(repository_root, "man")
  rd_database <- if (dir.exists(man_directory)) {
    # test_local() loads the package from source, where no installed help
    # database exists yet, so the .Rd files themselves are the ones to read.
    tools::Rd_db(dir = repository_root)
  } else {
    tools::Rd_db("cellucid")
  }
  rd <- rd_database[["cellucid_prepare.Rd"]]
  expect_false(is.null(rd))
  section_tags <- vapply(rd, function(node) attr(node, "Rd_tag"), character(1))

  usage <- trimws(paste(unlist(rd[section_tags == "\\usage"]), collapse = ""))
  expect_true(startsWith(usage, "cellucid_prepare("))
  documented_signature <- eval(parse(
    text = paste0(sub("^cellucid_prepare", "function", usage), "\nNULL")
  ))
  expect_identical(formals(documented_signature), formals(cellucid_prepare))

  arguments <- rd[section_tags == "\\arguments"][[1L]]
  argument_tags <- vapply(
    arguments,
    function(node) attr(node, "Rd_tag"),
    character(1)
  )
  documented_arguments <- vapply(
    arguments[argument_tags == "\\item"],
    function(item) trimws(paste(unlist(item[[1L]]), collapse = "")),
    character(1)
  )
  expect_identical(
    unname(documented_arguments),
    names(formals(cellucid_prepare))
  )

  namespace_path <- file.path(repository_root, "NAMESPACE")
  if (file.exists(namespace_path)) {
    namespace <- readLines(namespace_path, warn = FALSE, encoding = "UTF-8")
    # No generator input for this line exists under R/, so a regenerated
    # NAMESPACE would drop it and every .Call() would stop resolving.
    expect_true(
      any(namespace == 'useDynLib(cellucid, .registration = TRUE, .fixes = "C_")')
    )
    expect_false(any(grepl("Generated by roxygen2", namespace, fixed = TRUE)))
  }

  man_directory <- file.path(repository_root, "man")
  if (dir.exists(man_directory)) {
    man_paths <- list.files(man_directory, pattern = "\\.Rd$", full.names = TRUE)
    expect_gt(length(man_paths), 0L)
    for (man_path in man_paths) {
      man <- readLines(man_path, warn = FALSE, encoding = "UTF-8")
      expect_false(any(grepl("Generated by roxygen2", man, fixed = TRUE)))
      expect_true(any(grepl("Hand-written and authoritative", man, fixed = TRUE)))
    }
    prepare_path <- file.path(man_directory, "cellucid_prepare.Rd")
    if (file.exists(prepare_path)) {
      prepare <- readLines(prepare_path, warn = FALSE, encoding = "UTF-8")
      expect_true(any(grepl("\\examples{", prepare, fixed = TRUE)))
    }
  }

  # A second copy of the same help text under R/ is what diverged before. The
  # help page is the only copy, so no roxygen comment may reappear.
  r_directory <- file.path(repository_root, "R")
  if (dir.exists(r_directory)) {
    roxygen_comments <- character()
    for (source_path in list.files(r_directory, pattern = "\\.R$", full.names = TRUE)) {
      source_lines <- readLines(source_path, warn = FALSE, encoding = "UTF-8")
      matches <- grep("^[[:space:]]*#'", source_lines)
      if (length(matches) > 0L) {
        roxygen_comments <- c(
          roxygen_comments,
          sprintf("R/%s:%d", basename(source_path), matches)
        )
      }
    }
    expect_identical(roxygen_comments, character())
  }

  description_path <- file.path(repository_root, "DESCRIPTION")
  if (file.exists(description_path)) {
    fields <- colnames(read.dcf(description_path))
    expect_false("RoxygenNote" %in% fields)
    expect_false("Roxygen" %in% fields)
    expect_false(any(startsWith(fields, "Config/roxygen2/")))
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
