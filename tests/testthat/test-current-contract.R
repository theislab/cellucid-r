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

test_that("every path the current contract reads is present in this layout", {
  # The assertions in this file are worth exactly as much as the files they
  # read, and every one of them used to sit inside an `if (file.exists(...))`
  # or `if (dir.exists(...))`, so deleting a page, a workflow, or a config
  # turned its checks into no-ops rather than into a failure. The declared
  # list is checked here instead, in both directions: everything exists in the
  # checkout, and the entries .Rbuildignore strips really are absent from the
  # built package.
  layout <- cellucid_test_layout()
  required <- cellucid_required_paths()
  # Pinned so that deleting a row is a failure and not a quieter suite.
  expect_identical(nrow(required), 18L)
  expect_setequal(unique(required$kind), c("file", "directory"))
  expect_setequal(unique(required$layout), c("both", "checkout"))

  for (index in seq_len(nrow(required))) {
    relative <- required$path[[index]]
    absolute <- file.path(layout$root, relative)
    present <- if (required$kind[[index]] == "directory") {
      dir.exists(absolute)
    } else {
      file.exists(absolute) && !dir.exists(absolute)
    }
    expect_identical(
      present,
      layout$checkout || required$layout[[index]] == "both",
      info = relative
    )
  }
})

test_that("shipped R repository surfaces contain only current contract language", {
  repository_root <- cellucid_repository_root()
  roots <- c("R", "man", "vignettes")
  scalars <- c("README.md", "DESCRIPTION", "NAMESPACE")
  if (cellucid_is_source_checkout()) {
    roots <- c(roots, ".github")
    scalars <- c(scalars, "_pkgdown.yml", "publishing.md")
  }
  files <- unlist(
    lapply(
      file.path(repository_root, roots),
      list.files,
      recursive = TRUE,
      full.names = TRUE,
      all.files = TRUE
    ),
    use.names = FALSE
  )
  files <- c(files, file.path(repository_root, scalars))
  text_suffix <- "\\.(R|Rd|Rmd|md|ya?ml)$"
  files <- files[
    !grepl("\\.", basename(files), fixed = TRUE) |
      grepl(text_suffix, files, ignore.case = TRUE)
  ]
  # A scan that scanned nothing reports no violations, which is what this test
  # did under R CMD check while the root resolved to the check directory. Name
  # the sources it has to have reached.
  expect_true(all(
    file.path(repository_root, scalars) %in% files
  ))
  expect_true(all(
    file.path(
      repository_root,
      c("R/cellucid_prepare.R", "man/cellucid_prepare.Rd")
    ) %in% files
  ))
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

  if (cellucid_is_source_checkout()) {
    publishing <- paste(
      readLines(
        file.path(repository_root, "publishing.md"),
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
  repository_root <- cellucid_repository_root()
  readme_path <- file.path(repository_root, "README.md")
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

  if (cellucid_is_source_checkout()) {
    build_ignore <- readLines(
      file.path(repository_root, ".Rbuildignore"),
      warn = FALSE,
      encoding = "UTF-8"
    )
    expect_false(any(build_ignore == "^README\\.md$"))
  }
})

test_that("main is the only workflow branch", {
  repository_root <- cellucid_repository_root()
  workflow_directory <- file.path(repository_root, ".github", "workflows")

  if (!cellucid_is_source_checkout()) {
    # .Rbuildignore strips .github, so its absence from the built package is
    # the contract here rather than a check that did not run.
    expect_false(dir.exists(file.path(repository_root, ".github")))
    return(invisible(NULL))
  }

  workflow_branches <- function(name) {
    lines <- readLines(
      file.path(workflow_directory, name),
      warn = FALSE,
      encoding = "UTF-8"
    )
    trimws(lines[grepl("^[[:space:]]*branches:", lines)])
  }

  expect_identical(
    workflow_branches("R-CMD-check.yaml"),
    rep("branches: [main]", 2L)
  )
  expect_identical(
    workflow_branches("pkgdown.yaml"),
    "branches: [main]"
  )
})

test_that("DESCRIPTION is the one in-repo source of the 0.9.1 release identity", {
  repository_root <- cellucid_repository_root()
  description_path <- file.path(repository_root, "DESCRIPTION")
  version <- cellucid_description_version()
  expect_identical(version, "0.9.1")

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

  news <- readLines(
    file.path(repository_root, "NEWS.md"),
    warn = FALSE,
    encoding = "UTF-8"
  )
  expect_identical(
    news[[1L]],
    sprintf("# cellucid %s <!-- CELLUCID_VERSION -->", version)
  )
  expect_false(any(grepl("unreleased", tolower(news), fixed = TRUE)))
  expect_true(
    any(news == sprintf("Version %s is the CRAN submission release.", version))
  )

  if (cellucid_is_source_checkout()) {
    citation <- readLines(
      file.path(repository_root, "CITATION.cff"),
      warn = FALSE,
      encoding = "UTF-8"
    )
    expect_true(
      any(citation == sprintf("version: %s  # CELLUCID_VERSION", version))
    )
  }

  # README.md and the installation vignette are the two in-repo pages that
  # quote the version to users. Nothing outside this repository compares them
  # to DESCRIPTION, so every version string they contain is checked here.
  expect_identical(
    cellucid_semantic_versions(file.path(repository_root, "README.md")),
    version
  )

  expect_setequal(
    cellucid_semantic_versions(
      file.path(repository_root, "vignettes", "installation.Rmd")
    ),
    c(version, r_dependency)
  )

  package_citation_path <- file.path(repository_root, "inst", "CITATION")
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

  if (cellucid_is_source_checkout()) {
    release <- paste(
      readLines(
        file.path(repository_root, ".github", "workflows", "release.yaml"),
        warn = FALSE,
        encoding = "UTF-8"
      ),
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
    # Anchored on the step headers, not on the bare phrases. The file opens
    # with a comment naming `R CMD check --as-cran`, and `regexpr()` reports
    # the first match, so a search for the phrase alone lands in that comment
    # whatever the workflow body does: the ordering was satisfied even with the
    # check step renamed out of existence. `expect_gt(., 0L)` closes the other
    # end of the same hole -- an absent step reports position -1, which is less
    # than any real position and so satisfies an ordering all by itself.
    check_at <- regexpr(
      "- name: R CMD check --as-cran",
      release,
      fixed = TRUE
    )[[1L]]
    attach_at <- regexpr(
      "- name: Attach tarball to GitHub Release",
      release,
      fixed = TRUE
    )[[1L]]
    expect_gt(check_at, 0L)
    expect_gt(attach_at, 0L)
    expect_lt(check_at, attach_at)
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

  installation <- read_text(
    file.path(repository_root, "vignettes", "installation.Rmd")
  )
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

  readme <- read_text(file.path(repository_root, "README.md"))
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

  if (cellucid_is_source_checkout()) {
    pkgdown <- read_text(file.path(repository_root, "_pkgdown.yml"))
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
  # The .Rd sources themselves, in both layouts: test_local() loads the package
  # from source, where no installed help database exists yet, and under
  # R CMD check the unpacked sources carry man/ as well.
  rd_database <- tools::Rd_db(dir = repository_root)
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

  namespace <- readLines(
    file.path(repository_root, "NAMESPACE"),
    warn = FALSE,
    encoding = "UTF-8"
  )
  # No generator input for this line exists under R/, so a regenerated
  # NAMESPACE would drop it and every .Call() would stop resolving.
  expect_true(
    any(namespace == 'useDynLib(cellucid, .registration = TRUE, .fixes = "C_")')
  )
  expect_false(any(grepl("Generated by roxygen2", namespace, fixed = TRUE)))

  man_paths <- list.files(man_directory, pattern = "\\.Rd$", full.names = TRUE)
  prepare_path <- file.path(man_directory, "cellucid_prepare.Rd")
  # An empty listing would satisfy the loop below without reading anything, so
  # name the page every other assertion in this test already depends on.
  expect_true(prepare_path %in% man_paths)
  for (man_path in man_paths) {
    man <- readLines(man_path, warn = FALSE, encoding = "UTF-8")
    expect_false(
      any(grepl("Generated by roxygen2", man, fixed = TRUE)),
      info = basename(man_path)
    )
    expect_true(
      any(grepl("Hand-written and authoritative", man, fixed = TRUE)),
      info = basename(man_path)
    )
  }
  prepare <- readLines(prepare_path, warn = FALSE, encoding = "UTF-8")
  expect_true(any(grepl("\\examples{", prepare, fixed = TRUE)))

  # A second copy of the same help text under R/ is what diverged before. The
  # help page is the only copy, so no roxygen comment may reappear.
  r_directory <- file.path(repository_root, "R")
  r_sources <- list.files(r_directory, pattern = "\\.R$", full.names = TRUE)
  expect_true(file.path(r_directory, "cellucid_prepare.R") %in% r_sources)
  roxygen_comments <- character()
  for (source_path in r_sources) {
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

  fields <- colnames(read.dcf(file.path(repository_root, "DESCRIPTION")))
  expect_false("RoxygenNote" %in% fields)
  expect_false("Roxygen" %in% fields)
  expect_false(any(startsWith(fields, "Config/roxygen2/")))
})

test_that("categorical code width is one required exact caller choice", {
  expect_true(
    identical(
      formals(cellucid_prepare)[["obs_categorical_dtype"]],
      quote(expr = )
    )
  )

  # An argument that was never supplied is reported as missing, not as a value
  # of the wrong type, and the message says what a valid one is.
  missing_out <- tempfile("cellucid_r_missing_categorical_dtype_")
  expect_error(
    do.call(cellucid_prepare, current_contract_fixture(missing_out)),
    paste0(
      "^cellucid_prepare\\(\\) is missing 1 required argument:\\n",
      "  - obs_categorical_dtype: exactly one of \"uint8\" or \"uint16\", ",
      "the storage width of every categorical code\\.$"
    )
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

test_that("a required argument that was never supplied is reported as missing", {
  # No value of any of these three is a reasonable guess about someone else's
  # dataset, so none of them carries a default. Leaving one out is a different
  # mistake from supplying a bad one, and the two must not read alike.
  for (argument in c("obs_categorical_dtype", "dataset_name", "dataset_id")) {
    expect_true(
      identical(formals(cellucid_prepare)[[argument]], quote(expr = )),
      info = argument
    )
  }

  guidance <- c(
    obs_categorical_dtype = paste0(
      "exactly one of \"uint8\" or \"uint16\", the storage width of every ",
      "categorical code"
    ),
    dataset_name = "one non-empty name for the dataset, shown to the reader verbatim",
    dataset_id = paste0(
      "one portable identifier of 1-180 ASCII letters, numbers, '.', '_', or ",
      "'-', beginning with a letter or number and not ending with '.'"
    )
  )

  for (argument in names(guidance)) {
    out <- tempfile(paste0("cellucid_r_absent_", argument, "_"))
    args <- current_contract_fixture(out)
    args$obs_categorical_dtype <- "uint16"
    args[[argument]] <- NULL
    failure <- tryCatch(
      do.call(cellucid_prepare, args),
      error = conditionMessage
    )
    expect_identical(
      failure,
      paste0(
        "cellucid_prepare() is missing 1 required argument:\n",
        "  - ", argument, ": ", guidance[[argument]], "."
      )
    )
    expect_false(dir.exists(out))
  }

  # Every argument left out is named at once, in signature order, exactly as
  # the cellucid Python package names them in one TypeError.
  all_absent <- tempfile("cellucid_r_absent_all_")
  args <- current_contract_fixture(all_absent)
  args$dataset_name <- NULL
  args$dataset_id <- NULL
  failure <- tryCatch(
    do.call(cellucid_prepare, args),
    error = conditionMessage
  )
  expect_identical(
    failure,
    paste0(
      "cellucid_prepare() is missing 3 required arguments:\n",
      "  - obs_categorical_dtype: ", guidance[["obs_categorical_dtype"]], ".\n",
      "  - dataset_name: ", guidance[["dataset_name"]], ".\n",
      "  - dataset_id: ", guidance[["dataset_id"]], "."
    )
  )
  expect_false(dir.exists(all_absent))

  # A value that was supplied and is wrong stays a message about that value.
  supplied_out <- tempfile("cellucid_r_supplied_null_id_")
  supplied <- current_contract_fixture(supplied_out)
  supplied$obs_categorical_dtype <- "uint16"
  expect_error(
    do.call(
      cellucid_prepare,
      c(supplied[names(supplied) != "dataset_id"], list(dataset_id = NULL))
    ),
    "^dataset_id must be exactly one string\\.$"
  )
  expect_false(dir.exists(supplied_out))
})

test_that("dataset_description is optional and NULL and \"\" publish the same", {
  expect_null(formals(cellucid_prepare)[["dataset_description"]])

  described <- function(value, suffix) {
    out <- tempfile(paste0("cellucid_r_description_", suffix, "_"))
    args <- current_contract_fixture(out)
    args$obs_categorical_dtype <- "uint16"
    args$dataset_description <- value
    do.call(cellucid_prepare, args)
    identity <- jsonlite::read_json(
      file.path(out, "dataset_identity.json"),
      simplifyVector = FALSE
    )
    identity$description
  }

  expect_identical(described(NULL, "null"), "")
  expect_identical(described("", "empty"), "")
  expect_identical(described("An exact description", "text"), "An exact description")

  # Nothing else about the argument is relaxed: text the viewer cannot draw as
  # written is still rejected.
  padded_out <- tempfile("cellucid_r_description_padded_")
  padded <- current_contract_fixture(padded_out)
  padded$obs_categorical_dtype <- "uint16"
  padded$dataset_description <- " padded "
  expect_error(
    do.call(cellucid_prepare, padded),
    "^dataset_description is displayed verbatim"
  )
  expect_false(dir.exists(padded_out))

  wrong_type_out <- tempfile("cellucid_r_description_wrong_type_")
  wrong_type <- current_contract_fixture(wrong_type_out)
  wrong_type$obs_categorical_dtype <- "uint16"
  wrong_type$dataset_description <- 5
  expect_error(
    do.call(cellucid_prepare, wrong_type),
    "^dataset_description must be exactly one string\\.$"
  )
  expect_false(dir.exists(wrong_type_out))
})

current_contract_gene_names <- function(out_dir) {
  manifest <- jsonlite::read_json(
    file.path(out_dir, "var_manifest.json"),
    simplifyVector = FALSE
  )
  vapply(manifest$fields, function(field) field[[2L]], character(1))
}

test_that("NULL alone selects var row names and every string names one exact column", {
  gene_id_default <- formals(cellucid_prepare)$var_gene_id_column
  expect_null(gene_id_default)

  # The selected axis is now visible only in the manifest: the payload files
  # are named by index on every dataset, so no filename can report which
  # column was chosen.
  payload_names <- c("0.values.f32", "1.values.f32")

  rowname_out <- tempfile("cellucid_r_rowname_gene_ids_")
  rowname_args <- current_contract_fixture(rowname_out)
  rowname_args$obs_categorical_dtype <- "uint16"
  rowname_args$var_gene_id_column <- NULL
  do.call(cellucid_prepare, rowname_args)
  expect_identical(
    current_contract_gene_names(rowname_out),
    c("row_gene_1", "row_gene_2")
  )
  expect_identical(
    sort(list.files(file.path(rowname_out, "var"))),
    payload_names
  )

  column_out <- tempfile("cellucid_r_literal_index_column_")
  column_args <- current_contract_fixture(column_out)
  column_args$obs_categorical_dtype <- "uint16"
  column_args$var_gene_id_column <- "index"
  do.call(cellucid_prepare, column_args)
  expect_identical(
    current_contract_gene_names(column_out),
    c("column_gene_1", "column_gene_2")
  )
  expect_identical(
    sort(list.files(file.path(column_out, "var"))),
    payload_names
  )

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
