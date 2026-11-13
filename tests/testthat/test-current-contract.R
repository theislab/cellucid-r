test_that("every path the current contract reads is present in this layout", {
  # The assertions in this file are worth exactly as much as the files they
  # read: inside an `if (file.exists(...))` or `if (dir.exists(...))`, a
  # deleted page, workflow, or config would turn its checks into no-ops
  # rather than into a failure. The declared list is checked here instead, in
  # both directions: everything exists in the checkout, and the entries
  # .Rbuildignore strips really are absent from the built package.
  layout <- cellucid_test_layout()
  required <- cellucid_required_paths()
  # Pinned so that deleting a row is a failure and not a quieter suite.
  expect_identical(nrow(required), 19L)
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
  # The github-pages environment accepts only main: a release-triggered run
  # deploys from the tag ref and is refused.
  pkgdown_lines <- readLines(
    file.path(workflow_directory, "pkgdown.yaml"),
    warn = FALSE,
    encoding = "UTF-8"
  )
  expect_false(any(grepl("^[[:space:]]*release:", pkgdown_lines)))
})

test_that("DESCRIPTION is the one in-repo source of the release identity", {
  repository_root <- cellucid_repository_root()
  description_path <- file.path(repository_root, "DESCRIPTION")
  version <- cellucid_description_version()

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
  citation_environment <- new.env(parent = globalenv())
  citation_environment$meta <- list(Version = version)
  sys.source(package_citation_path, envir = citation_environment)
  expect_identical(citation_environment$package_version, version)
})

# The release workflow, parsed, and its steps by what they do.
release_workflow <- function() {
  workflow <- yaml::read_yaml(file.path(cellucid_repository_root(), ".github", "workflows", "release.yaml"))
  steps <- workflow$jobs$build$steps
  runs <- vapply(steps, function(step) if (is.character(step$run)) step$run else "", character(1L))
  uses <- vapply(steps, function(step) if (is.character(step$uses)) step$uses else "", character(1L))
  # YAML 1.1 reads the key `on` as TRUE.
  list(workflow = workflow, steps = steps, runs = runs, uses = uses,
       on = if (is.null(workflow[["on"]])) workflow[["TRUE"]] else workflow[["on"]])
}

test_that("the release workflow runs on a published release only, and checks before it attaches and submits", {
  skip_if_not(cellucid_is_source_checkout(), "the workflows are in the checkout only")
  skip_if_not_installed("yaml")
  release <- release_workflow()
  # One trigger: a published GitHub Release (no tag push, no dispatch).
  expect_identical(names(release$on), "release")
  expect_identical(unlist(release$on$release$types), "published")

  # Each step is found by what it does; the order is the guarantee.
  position <- function(hits) {
    expect_length(which(hits), 1L)
    which(hits)
  }
  verify <- position(grepl("GITHUB_REF_TYPE", release$runs, fixed = TRUE))
  build <- position(grepl("R CMD build", release$runs, fixed = TRUE))
  check <- position(grepl("(^|\n)\\s*R CMD check --as-cran", release$runs))
  attach <- position(startsWith(release$uses, "softprops/action-gh-release@"))
  refuse <- position(grepl("00check.log", release$runs, fixed = TRUE))
  submit <- position(grepl("upload_cran", release$runs, fixed = TRUE))
  expect_true(all(diff(c(verify, build, check, attach, refuse, submit)) > 0L))

  # Every step output a step reads is written by the step it names.
  text <- paste(unlist(release$steps), collapse = "\n")
  referenced <- regmatches(text, gregexpr("steps\\.[a-z]+\\.outputs\\.[a-z]+", text))[[1L]]
  expect_gt(length(referenced), 0L)
  for (reference in unique(referenced)) {
    parts <- strsplit(reference, ".", fixed = TRUE)[[1L]]
    writer <- Filter(function(step) identical(step$id, parts[[2L]]), release$steps)
    expect_length(writer, 1L)
    expect_true(grepl(paste0(parts[[4L]], "="), writer[[1L]]$run, fixed = TRUE), info = reference)
  }
  # The release assets are the tarball and its checksum, which the build step makes.
  expect_true(grepl("steps.build.outputs.checksum", release$steps[[attach]]$with$files, fixed = TRUE))
  expect_true(grepl("sha256sum", release$runs[[build]], fixed = TRUE))
})

test_that("the release verification passes a final version on its tag, from main, and nothing else", {
  skip_if_not(cellucid_is_source_checkout(), "the workflows are in the checkout only")
  skip_if_not_installed("yaml")
  skip_if(!nzchar(Sys.which("git")) || !nzchar(Sys.which("bash")), "git and bash are needed to run the step")
  release <- release_workflow()
  step <- release$steps[[which(grepl("GITHUB_REF_TYPE", release$runs, fixed = TRUE))]]

  git <- function(directory, ...) {
    output <- system2("git", c("-C", shQuote(directory), "-c", "user.name=test", "-c", "user.email=test@example.invalid", ...),
                      stdout = TRUE, stderr = TRUE)
    expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
    output
  }
  work <- tempfile("release-work-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  repository <- file.path(work, "package")
  dir.create(repository)
  git(repository, "init", "-q", "-b", "main")
  writeLines("Version: 1.2.3", file.path(repository, "DESCRIPTION"))
  git(repository, "add", "DESCRIPTION")
  git(repository, "commit", "-q", "-m", "release")
  # The remote is a bare copy of this scratch repository, never a real one.
  git(work, "clone", "-q", "--bare", shQuote(repository), shQuote(file.path(work, "origin.git")))
  git(repository, "remote", "add", "origin", shQuote(file.path(work, "origin.git")))
  head <- git(repository, "rev-parse", "HEAD")
  git(repository, "checkout", "-q", "-b", "side")
  writeLines("x", file.path(repository, "side.txt"))
  git(repository, "add", "side.txt")
  git(repository, "commit", "-q", "-m", shQuote("not on main"))
  side <- git(repository, "rev-parse", "HEAD")
  git(repository, "checkout", "-q", "main")

  passes <- function(..., commit = head, description = NULL) {
    environment <- c(
      RELEASE_TAG = "v1.2.3", RELEASE_DRAFT = "false", RELEASE_PRERELEASE = "false",
      GITHUB_REF_TYPE = "tag", GITHUB_REF_NAME = "v1.2.3", GITHUB_SHA = commit
    )
    changes <- list(...)
    environment[names(changes)] <- unlist(changes)
    # What `shell: bash` runs: bash -e -o pipefail on the step's script.
    script <- tempfile("release-step-", work)
    writeLines(step$run, script)
    if (!is.null(description)) writeLines(sub("^Version: .*", description, readLines(file.path(repository, "DESCRIPTION"))), file.path(repository, "DESCRIPTION"))
    on.exit(writeLines("Version: 1.2.3", file.path(repository, "DESCRIPTION")), add = TRUE)
    output <- suppressWarnings(system2("env", c(paste0(names(environment), "=", environment), "bash", "-e", "-o", "pipefail", shQuote(script)),
                                       stdout = TRUE, stderr = TRUE))
    is.null(attr(output, "status"))
  }
  old <- setwd(repository)
  on.exit(setwd(old), add = TRUE, after = FALSE)

  expect_true(passes())
  expect_false(passes(RELEASE_DRAFT = "true"))
  expect_false(passes(RELEASE_PRERELEASE = "true"))
  expect_false(passes(GITHUB_REF_TYPE = "branch"))
  expect_false(passes(RELEASE_TAG = "v1.2.4"))
  expect_false(passes(GITHUB_REF_NAME = "v1.2.4", RELEASE_TAG = "v1.2.4"))
  expect_false(passes(description = "Version: 1.2.3.9000"))
  expect_false(passes(commit = side))
})

# Every path of a sibling repository the tests name, as list(repository, path,
# whole): the calls of expect_shared_file() (its `theirs`),
# cellucid_sibling_file() and cellucid_rendered_by_web(), read from the parsed
# test files. `whole` is FALSE where only the leading folders are literal
# (file.path("tests", "fixtures", name)), and `path` is then those folders.
sibling_paths_the_tests_name <- function(test_dir) {
  found <- list()
  add <- function(repository, parts) {
    if (!is.character(repository) || length(repository) != 1L) return(invisible(NULL))
    literal <- vapply(parts, function(part) is.character(part) && length(part) == 1L, logical(1L))
    whole <- all(literal)
    text <- paste(unlist(parts[cumprod(literal) == 1L]), collapse = "/")
    if (!nzchar(text)) return(invisible(NULL))
    found[[length(found) + 1L]] <<- list(repository = repository, path = if (whole) text else paste0(text, "/"), whole = whole)
  }
  parts_of <- function(expression) {
    if (is.call(expression) && identical(expression[[1L]], quote(file.path))) as.list(expression)[-1L] else list(expression)
  }
  visit <- function(x) {
    if (is.call(x)) {
      name <- if (is.name(x[[1L]])) as.character(x[[1L]]) else ""
      if (identical(name, "expect_shared_file")) {
        matched <- match.call(function(ours, theirs, repository = "cellucid", then = "") NULL, x)
        add(if (is.null(matched$repository)) "cellucid" else matched$repository, parts_of(matched$theirs))
      } else if (identical(name, "cellucid_sibling_file") && length(x) > 2L) {
        add(x[[2L]], as.list(x)[-(1:2)])
      } else if (identical(name, "cellucid_rendered_by_web") && length(x) > 1L) {
        add("cellucid", list(x[[2L]]))
      }
      for (part in as.list(x)) if (!missing(part)) visit(part)
    } else if (is.pairlist(x) || is.expression(x)) {
      for (part in as.list(x)) if (!missing(part)) visit(part)
    }
  }
  for (file in list.files(test_dir, "^test-.*\\.R$", full.names = TRUE)) {
    for (expression in parse(file, keep.source = FALSE)) visit(expression)
  }
  found
}

# What the sparse checkout holds: its patterns hold `path` (or, where only the
# leading folders `path` are known, something under them).
sparse_pattern_holds <- function(pattern, path, whole) {
  pattern <- sub("^/", "", pattern)
  if (endsWith(pattern, "/")) return(startsWith(path, pattern) || (!whole && startsWith(pattern, path)))
  identical(pattern, path) || (!whole && startsWith(pattern, path))
}

test_that("CI fetches every sibling the cross-repository tests compare with, at main, and requires them", {
  # Without them those tests skip (helper-repository.R); CI fetches both at
  # main and sets CELLUCID_REQUIRE_SIBLINGS=1, so there they always run.
  skip_if_not(cellucid_is_source_checkout(), "the workflows are in the checkout only")
  skip_if_not_installed("yaml")
  root <- cellucid_repository_root()
  workflow <- yaml::read_yaml(file.path(root, ".github", "workflows", "R-CMD-check.yaml"))
  job <- workflow$jobs[["R-CMD-check"]]
  scripts <- Filter(is.character, lapply(job$steps, function(step) step$run))
  script <- paste(unlist(scripts), collapse = "\n")
  uses <- unlist(lapply(job$steps, function(step) step$uses))
  expect_true(any(startsWith(uses, "actions/setup-node@")))
  expect_match(script, "CELLUCID_REQUIRE_SIBLINGS=1", fixed = TRUE)

  # The `fetch <repository> <ref> <directory> <patterns...>` calls of the step
  # that checks the siblings out, each with its continuation lines joined.
  calls <- grep("^\\s*fetch (cellucid|cellucid-python) ", strsplit(gsub("\\\\\n", " ", script), "\n")[[1L]], value = TRUE)
  fetched <- lapply(strsplit(trimws(calls), "[[:space:]]+"), function(tokens) {
    list(repository = tokens[[2L]], ref = tokens[[3L]], directory = gsub('"', "", tokens[[4L]]), patterns = tokens[-(1:4)])
  })
  names(fetched) <- vapply(fetched, function(call) call$repository, character(1L))
  expect_setequal(names(fetched), c("cellucid", "cellucid-python"))

  named <- sibling_paths_the_tests_name(test_path())
  for (repository in names(fetched)) {
    checkout <- fetched[[repository]]
    # At main (the Python one through the job's variable), named for the tests by its root variable.
    ref <- if (identical(checkout$ref, "\"$CELLUCID_PYTHON_BRANCH\"")) job$env$CELLUCID_PYTHON_BRANCH else checkout$ref
    expect_identical(ref, "main", info = repository)
    variable <- cellucid_sibling_variable(repository)
    expect_true(grepl(paste0(variable, "=", checkout$directory), script, fixed = TRUE), info = variable)
    # Its sparse checkout holds every path the tests name of it.
    wanted <- Filter(function(path) identical(path$repository, repository), named)
    expect_gt(length(wanted), 0L)
    left_out <- vapply(wanted, function(path) {
      !any(vapply(checkout$patterns, sparse_pattern_holds, logical(1L), path = path$path, whole = path$whole))
    }, logical(1L))
    expect_identical(vapply(wanted[left_out], function(path) path$path, character(1L)), character(0L), info = repository)
  }
})

test_that("the pkgdown site contract holds: every vignette and help page indexed once", {
  # .github/scripts/validate-pkgdown.R gates the site's build (pkgdown.yaml);
  # it runs here too, so a help page added or removed fails the suite and
  # R-CMD-check, not only the site's deploy.
  skip_if_not(cellucid_is_source_checkout(), "the site contract is in the checkout only")
  root <- cellucid_repository_root()
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", ".github/scripts/validate-pkgdown.R"),
    stdout = TRUE, stderr = TRUE
  ))
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
})

test_that("man/ and NAMESPACE are hand-written and match the code exactly", {
  # This is the codoc check R CMD check performs, run against the installed Rd
  # database: the \usage block is turned back into a function and compared to
  # the real formals, argument names, order, and default expressions included.
  # It is what keeps a hand-written help page from drifting away from the
  # signature, so no generator is needed to hold the two together.
  #
  # tools::codoc() itself is not usable here: it loads and unloads the package
  # it inspects, which would unload the namespace the running suite uses.
  repository_root <- cellucid_repository_root()
  man_directory <- file.path(repository_root, "man")
  # The .Rd sources themselves, in both layouts: test_local() loads the package
  # from source, where no installed help database exists yet, and under
  # R CMD check the unpacked sources carry man/ as well.
  rd_database <- tools::Rd_db(dir = repository_root)
  # Every exported function is documented by exactly one page's \usage, and
  # each call there matches the function's formals; a page's \arguments list
  # every argument of its functions once, in order. The pages of the actions'
  # functions are generated from the action catalog (test-actions.R holds
  # them to it); every other page is hand-written.
  hand_written <- c(
    "cellucid_prepare", "cellucid_transition_drift", "cellucid_view", "cellucid_serve",
    "cellucid_stop", "cellucid_stop_all", "cellucid_viewers", "cellucid_act", "cellucid_action_result",
    "cellucid_on", "cellucid_off", "cellucid_selection", "cellucid_wait_event", "cellucid_wait_ready",
    "cellucid_clear_web_cache", "cellucid_debug_connection", "cellucid_apply_session", "cellucid_read_session",
    "cellucid_apply_plan", "cellucid_history_applied", "cellucid_trackers", "cellucid_tracker", "cellucid_branches",
    "cellucid_branch", "cellucid_groups", "cellucid_group", "cellucid_steps", "cellucid_step", "cellucid_step_state", "cellucid_step_diff",
    "cellucid_stale", "cellucid_notes", "cellucid_turns", "cellucid_journal", "cellucid_methods", "cellucid_pseudocode",
    "cellucid_conversations", "cellucid_turn", "cellucid_load_session", "cellucid_load_session_result",
    "cellucid_last_history_change", "cellucid_mcp",
    "cellucid_subsample_cells", "cellucid_restore_original_order", "cellucid_pull", "cellucid_results", "cellucid_subset",
    "cellucid_gene_sets", "cellucid_push", "cellucid_write_gmt",
    "cellucid_read_gmt", "cellucid_request_session", "cellucid_session_result"
  )
  catalog <- jsonlite::read_json(file.path(repository_root, "inst", "action_catalog.json"), simplifyVector = FALSE)
  generated <- vapply(catalog$actions, function(action) paste0("cellucid_", action$tool), "")
  expect_length(intersect(generated, hand_written), 0L)
  exported <- c(hand_written, generated)
  documented <- character()
  for (page in names(rd_database)) {
    rd <- rd_database[[page]]
    section_tags <- vapply(rd, function(node) attr(node, "Rd_tag"), character(1))
    usage <- paste(unlist(rd[section_tags == "\\usage"]), collapse = "")
    calls <- as.list(parse(text = usage))
    page_arguments <- character()
    for (call in calls) {
      function_name <- as.character(call[[1L]])
      documented <- c(documented, function_name)
      implementation <- get(function_name, envir = asNamespace("cellucid"))
      given <- as.list(call)[-1L]
      labels <- names(given)
      if (is.null(labels)) labels <- rep("", length(given))
      formal_names <- ifelse(nzchar(labels), labels, vapply(given, deparse, character(1)))
      documented_formals <- stats::setNames(
        lapply(seq_along(given), function(index) if (nzchar(labels[[index]])) given[[index]] else quote(expr = )),
        formal_names
      )
      expect_identical(as.pairlist(documented_formals), formals(implementation), info = function_name)
      page_arguments <- c(page_arguments, setdiff(names(formals(implementation)), page_arguments))
    }
    arguments <- rd[section_tags == "\\arguments"]
    if (length(page_arguments) == 0L) {
      expect_length(arguments, 0L)
      next
    }
    arguments <- arguments[[1L]]
    argument_tags <- vapply(arguments, function(node) attr(node, "Rd_tag"), character(1))
    documented_arguments <- vapply(
      arguments[argument_tags == "\\item"],
      function(item) trimws(paste(unlist(item[[1L]]), collapse = "")),
      character(1)
    )
    expect_identical(unname(documented_arguments), page_arguments, info = page)
  }
  expect_setequal(documented, exported)
  expect_false(anyDuplicated(documented) > 0L)

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
  expect_identical(
    sort(grep("^export\\(", namespace, value = TRUE)),
    paste0("export(", sort(exported), ")")
  )

  man_paths <- list.files(man_directory, pattern = "\\.Rd$", full.names = TRUE)
  prepare_path <- file.path(man_directory, "cellucid_prepare.Rd")
  # An empty listing would satisfy the loop below without reading anything, so
  # name the page every other assertion in this test already depends on.
  expect_true(prepare_path %in% man_paths)
  generated_pages <- file.path(man_directory, paste0(generated, ".Rd"))
  expect_true(all(generated_pages %in% man_paths))
  for (man_path in man_paths) {
    man <- readLines(man_path, warn = FALSE, encoding = "UTF-8")
    expect_false(
      any(grepl("Generated by roxygen2", man, fixed = TRUE)),
      info = basename(man_path)
    )
    if (man_path %in% generated_pages) {
      expect_true(startsWith(man[[1L]], "% Generated by tools/generate-actions.R from inst/action_catalog.json"),
                  info = basename(man_path))
    } else {
      expect_true(
        any(grepl("Hand-written and authoritative", man, fixed = TRUE)),
        info = basename(man_path)
      )
    }
  }
  prepare <- readLines(prepare_path, warn = FALSE, encoding = "UTF-8")
  expect_true(any(grepl("\\examples{", prepare, fixed = TRUE)))

  # A second copy of the same help text under R/ would diverge from the page.
  # The help page is the only copy, so no roxygen comment may appear.
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

test_that("R code is ASCII outside comments, as R CMD check --as-cran requires", {
  # CRAN's portability check warns on any non-ASCII character in code (a
  # string included), and a WARNING stops the release workflow. Characters
  # such as the middle dot and the ellipsis are written as \u escapes.
  # Tokens are read from the parse, so a `#` inside a string is not a comment.
  files <- list.files(file.path(cellucid_repository_root(), "R"), pattern = "[.]R$", full.names = TRUE)
  expect_gt(length(files), 0L)
  offenders <- character()
  for (file in files) {
    data <- utils::getParseData(parse(file, keep.source = TRUE, encoding = "UTF-8"), includeText = TRUE)
    code <- data[data$terminal & data$token != "COMMENT", , drop = FALSE]
    bad <- code[grepl("[^\x01-\x7f]", code$text, useBytes = TRUE), , drop = FALSE]
    if (nrow(bad) > 0L) offenders <- c(offenders, sprintf("%s:%d", basename(file), bad$line1))
  }
  expect_identical(offenders, character())
})
