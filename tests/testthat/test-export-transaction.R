.transaction_test_parent <- function(label) {
  parent <- tempfile(paste0("cellucid_r_transaction_", label, "_"))
  dir.create(parent)
  parent
}

.transaction_write_journal <- function(
    target,
    transaction_id,
    had_target
) {
  paths <- cellucid:::.export_transaction_paths(
    target,
    transaction_id
  )
  writeBin(
    charToRaw(cellucid:::.serialize_export_transaction(
      transaction_id,
      had_target
    )),
    paths$journal
  )
  paths
}

.transaction_create_directory <- function(path, owner) {
  expect_true(dir.create(path))
  writeLines(owner, file.path(path, "owner"), useBytes = TRUE)
}

.transaction_snapshot <- function(path) {
  info <- cellucid:::.export_path_info(path)
  if (info$kind == 0L) {
    return(NULL)
  }
  if (info$kind == 2L) {
    files <- list.files(
      path,
      recursive = TRUE,
      all.files = TRUE,
      no.. = TRUE,
      full.names = TRUE
    )
    files <- files[!dir.exists(files)]
    values <- lapply(files, function(file) {
      readBin(file, what = "raw", n = file.info(file)$size)
    })
    names(values) <- substring(files, nchar(path) + 2L)
    return(values)
  }
  readBin(path, what = "raw", n = file.info(path)$size)
}

.transaction_prepare <- function(
    out,
    dataset_id,
    force,
    created_at = NULL
) {
  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  obs <- data.frame(group = factor(c("A", "B")))
  umap1 <- matrix(c(0, 1), ncol = 1)
  cellucid_prepare(
    dataset_id = dataset_id,
    dataset_name = dataset_id,
    latent_space = latent,
    obs = obs,
    X_umap_1d = umap1,
    out_dir = out,
    centroid_min_points = 1,
    force = force,
    obs_categorical_dtype = "uint16",
    created_at = created_at
  )
}

.transaction_dataset_id <- function(target) {
  jsonlite::read_json(
    file.path(target, "dataset_identity.json"),
    simplifyVector = FALSE
  )$id
}

.expect_transaction_controls_cleared <- function(target) {
  parent <- normalizePath(
    dirname(target),
    winslash = "/",
    mustWork = TRUE
  )
  siblings <- list.files(parent, all.files = TRUE)
  expect_false(any(startsWith(
    siblings,
    paste0(".", basename(target), ".cellucid-")
  )))
}

.expect_transaction_reparse_node <- function(path) {
  expect_identical(cellucid:::.export_path_info(path)$kind, 3L)
  if (.Platform$OS.type != "windows") {
    link_target <- Sys.readlink(path)
    expect_identical(length(link_target), 1L)
    expect_false(is.na(link_target))
    expect_true(nzchar(link_target))
  }
  invisible(path)
}

.expect_transaction_cleanup <- function(
    parent,
    file_reparse_nodes = character(),
    directory_reparse_nodes = character(),
    after_reparse_cleanup = function() invisible(NULL)
) {
  stopifnot(is.function(after_reparse_cleanup))
  reparse_nodes <- unique(c(
    file_reparse_nodes,
    directory_reparse_nodes
  ))
  for (path in reparse_nodes) {
    if (cellucid:::.export_path_info(path)$kind == 0L) {
      next
    }
    .expect_transaction_reparse_node(path)
    status <- NULL
    expect_warning(
      status <- .remove_test_reparse(
        path,
        directory = path %in% directory_reparse_nodes
      ),
      regexp = NA
    )
    expect_identical(status, 0L)
    expect_identical(cellucid:::.export_path_info(path)$kind, 0L)
  }
  after_reparse_cleanup()

  if (cellucid:::.export_path_info(parent)$kind != 0L) {
    status <- NULL
    expect_warning(
      status <- unlink(parent, recursive = TRUE, force = TRUE),
      regexp = NA
    )
    expect_identical(status, 0L)
    expect_identical(cellucid:::.export_path_info(parent)$kind, 0L)
  }
  invisible(parent)
}

test_that("transaction descriptor bytes and reserved paths match Python exactly", {
  parent <- .transaction_test_parent("canonical_bytes")
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  target <- file.path(parent, "generation")
  transaction_id <- "0123456789abcdef0123456789abcdef"
  paths <- cellucid:::.export_transaction_paths(
    target,
    transaction_id
  )

  expect_identical(
    basename(paths$journal),
    ".generation.cellucid-transaction.json"
  )
  expect_identical(
    basename(paths$journal_temp),
    ".generation.cellucid-transaction.json.tmp"
  )
  expect_identical(
    basename(paths$stage),
    paste0(".generation.cellucid-stage-", transaction_id)
  )
  expect_identical(
    basename(paths$backup),
    paste0(".generation.cellucid-backup-", transaction_id)
  )
  expect_identical(
    charToRaw(cellucid:::.serialize_export_transaction(
      transaction_id,
      TRUE
    )),
    charToRaw(paste0(
      "{\"format\":\"cellucid-export-transaction\",\"version\":1,",
      "\"transaction_id\":\"",
      transaction_id,
      "\",\"had_target\":true}\n"
    ))
  )
  expect_identical(
    charToRaw(cellucid:::.serialize_export_transaction(
      transaction_id,
      FALSE
    )),
    charToRaw(paste0(
      "{\"format\":\"cellucid-export-transaction\",\"version\":1,",
      "\"transaction_id\":\"",
      transaction_id,
      "\",\"had_target\":false}\n"
    ))
  )

  cellucid:::.write_export_transaction(
    target,
    transaction_id,
    FALSE
  )
  expect_identical(
    readBin(paths$journal, what = "raw", n = 512L),
    charToRaw(cellucid:::.serialize_export_transaction(
      transaction_id,
      FALSE
    ))
  )
  expect_false(file.exists(paths$journal_temp))
  expect_identical(
    cellucid:::.read_export_transaction(paths$journal),
    list(transaction_id = transaction_id, had_target = FALSE)
  )
})

test_that("retired random stage and backup helpers have no compatibility aliases", {
  namespace <- asNamespace("cellucid")
  for (name in c(
    ".discard_staged_output",
    ".begin_staged_output",
    ".publish_staged_output"
  )) {
    expect_false(
      exists(name, envir = namespace, inherits = FALSE),
      info = name
    )
  }
})

test_that("native transaction identities are exact and collision-free", {
  identities <- replicate(
    256L,
    cellucid:::.native_export_transaction_id()
  )
  expect_true(all(grepl(
    "^[0-9a-f]{32}$",
    identities,
    perl = TRUE,
    useBytes = TRUE
  )))
  expect_identical(length(unique(identities)), length(identities))
})

test_that("every recoverable Python and R transaction state converges", {
  cases <- list(
    list(TRUE, c(TRUE, FALSE, FALSE), "target"),
    list(TRUE, c(TRUE, TRUE, FALSE), "target"),
    list(TRUE, c(FALSE, TRUE, TRUE), "backup"),
    list(TRUE, c(TRUE, FALSE, TRUE), "target"),
    list(FALSE, c(FALSE, FALSE, FALSE), NULL),
    list(FALSE, c(FALSE, TRUE, FALSE), NULL),
    list(FALSE, c(TRUE, FALSE, FALSE), "target")
  )
  transaction_id <- "0123456789abcdef0123456789abcdef"

  for (index in seq_along(cases)) {
    case <- cases[[index]]
    parent <- .transaction_test_parent(paste0("valid_", index))
    target <- file.path(parent, "generation")
    paths <- .transaction_write_journal(
      target,
      transaction_id,
      case[[1L]]
    )
    state <- case[[2L]]
    for (node_index in seq_along(state)) {
      if (state[[node_index]]) {
        path <- c(target, paths$stage, paths$backup)[[node_index]]
        owner <- c("target", "stage", "backup")[[node_index]]
        .transaction_create_directory(path, owner)
      }
    }

    expect_no_error(
      cellucid:::.recover_export_transaction(target)
    )
    expect_false(file.exists(paths$journal))
    expect_false(file.exists(paths$journal_temp))
    expect_false(file.exists(paths$stage))
    expect_false(file.exists(paths$backup))
    if (is.null(case[[3L]])) {
      expect_false(dir.exists(target))
    } else {
      expect_identical(
        readLines(file.path(target, "owner"), warn = FALSE),
        case[[3L]]
      )
    }
    unlink(parent, recursive = TRUE, force = TRUE)
  }
})

test_that("all nine unspecified transaction states fail closed", {
  valid <- list(
    `TRUE` = c("100", "110", "011", "101"),
    `FALSE` = c("000", "010", "100")
  )
  transaction_id <- "abcdef0123456789abcdef0123456789"

  for (had_target in c(TRUE, FALSE)) {
    for (target_exists in c(FALSE, TRUE)) {
      for (stage_exists in c(FALSE, TRUE)) {
        for (backup_exists in c(FALSE, TRUE)) {
          state <- c(target_exists, stage_exists, backup_exists)
          state_key <- paste0(as.integer(state), collapse = "")
          if (state_key %in% valid[[as.character(had_target)]]) {
            next
          }
          parent <- .transaction_test_parent(
            paste0(
              "invalid_",
              tolower(as.character(had_target)),
              "_",
              state_key
            )
          )
          target <- file.path(parent, "generation")
          paths <- .transaction_write_journal(
            target,
            transaction_id,
            had_target
          )
          nodes <- c(target, paths$stage, paths$backup)
          for (node_index in seq_along(state)) {
            if (state[[node_index]]) {
              .transaction_create_directory(
                nodes[[node_index]],
                c("target", "stage", "backup")[[node_index]]
              )
            }
          }
          before <- lapply(nodes, .transaction_snapshot)
          journal_before <- readBin(
            paths$journal,
            what = "raw",
            n = 512L
          )

          expect_error(
            cellucid:::.recover_export_transaction(target),
            "cannot be recovered without guessing",
            info = paste(had_target, state_key)
          )
          expect_identical(
            readBin(paths$journal, what = "raw", n = 512L),
            journal_before
          )
          expect_identical(
            lapply(nodes, .transaction_snapshot),
            before
          )
          unlink(parent, recursive = TRUE, force = TRUE)
        }
      }
    }
  }
})

test_that("malformed and noncanonical journals never mutate a generation", {
  transaction_id <- "0123456789abcdef0123456789abcdef"
  canonical <- cellucid:::.serialize_export_transaction(
    transaction_id,
    TRUE
  )
  malformed <- list(
    trailing_space = charToRaw(sub("\n$", " \n", canonical)),
    missing_newline = charToRaw(sub("\n$", "", canonical)),
    leading_newline = charToRaw(paste0("\n", canonical)),
    decimal_version = charToRaw(sub(
      "\"version\":1",
      "\"version\":1.0",
      canonical,
      fixed = TRUE
    )),
    uppercase_identity = charToRaw(sub(
      transaction_id,
      toupper(transaction_id),
      canonical,
      fixed = TRUE
    )),
    extra_member = charToRaw(sub(
      "}\n$",
      ",\"extra\":null}\n",
      canonical
    )),
    duplicate_member = charToRaw(sub(
      "\"version\":1,",
      "\"version\":1,\"version\":1,",
      canonical,
      fixed = TRUE
    )),
    wrong_order = charToRaw(paste0(
      "{\"version\":1,\"format\":\"cellucid-export-transaction\",",
      "\"transaction_id\":\"",
      transaction_id,
      "\",\"had_target\":true}\n"
    )),
    string_boolean = charToRaw(sub(
      "\"had_target\":true",
      "\"had_target\":\"true\"",
      canonical,
      fixed = TRUE
    )),
    non_ascii = c(charToRaw(canonical), as.raw(255L)),
    embedded_nul = c(charToRaw(canonical), as.raw(0L)),
    oversized = as.raw(rep.int(65L, 513L))
  )

  for (name in names(malformed)) {
    parent <- .transaction_test_parent(paste0("malformed_", name))
    target <- file.path(parent, "generation")
    .transaction_create_directory(target, "prior")
    paths <- cellucid:::.export_transaction_paths(
      target,
      transaction_id
    )
    writeBin(malformed[[name]], paths$journal)
    before <- .transaction_snapshot(target)
    journal_before <- readBin(
      paths$journal,
      what = "raw",
      n = 1024L
    )

    expect_error(
      cellucid:::.recover_export_transaction(target),
      "malformed|canonical|large",
      info = name
    )
    expect_identical(.transaction_snapshot(target), before)
    expect_identical(
      readBin(paths$journal, what = "raw", n = 1024L),
      journal_before
    )
    unlink(parent, recursive = TRUE, force = TRUE)
  }
})

test_that("unsafe journal and temporary control nodes fail without mutation", {
  transaction_id <- "0123456789abcdef0123456789abcdef"

  parent <- .transaction_test_parent("journal_hard_link")
  target <- file.path(parent, "generation")
  .transaction_create_directory(target, "prior")
  paths <- cellucid:::.export_transaction_paths(target, transaction_id)
  source <- file.path(parent, "journal-source")
  writeBin(
    charToRaw(cellucid:::.serialize_export_transaction(
      transaction_id,
      TRUE
    )),
    source
  )
  skip_if_not(file.link(source, paths$journal), "hard links unavailable")
  before <- .transaction_snapshot(target)
  expect_error(
    cellucid:::.recover_export_transaction(target),
    "non-linked"
  )
  expect_identical(.transaction_snapshot(target), before)
  expect_true(file.exists(source))
  expect_true(file.exists(paths$journal))
  unlink(parent, recursive = TRUE, force = TRUE)

  parent <- .transaction_test_parent("journal_symlink")
  target <- file.path(parent, "generation")
  .transaction_create_directory(target, "prior")
  paths <- cellucid:::.export_transaction_paths(target, transaction_id)
  source <- file.path(parent, "journal-source")
  writeBin(
    charToRaw(cellucid:::.serialize_export_transaction(
      transaction_id,
      TRUE
    )),
    source
  )
  link_created <- isTRUE(file.symlink(source, paths$journal))
  if (!link_created) {
    .expect_transaction_cleanup(parent)
  }
  skip_if_not(link_created, "symbolic links unavailable")
  expected_journal_source <- charToRaw(
    cellucid:::.serialize_export_transaction(
      transaction_id,
      TRUE
    )
  )
  expect_error(
    cellucid:::.recover_export_transaction(target),
    "non-linked|non-symbolic"
  )
  expect_identical(
    readBin(source, what = "raw", n = 512L),
    expected_journal_source
  )
  .expect_transaction_reparse_node(paths$journal)
  .expect_transaction_cleanup(
    parent,
    file_reparse_nodes = paths$journal,
    after_reparse_cleanup = function() {
      expect_identical(
        readBin(source, what = "raw", n = 512L),
        expected_journal_source
      )
    }
  )

  parent <- .transaction_test_parent("temporary_directory")
  target <- file.path(parent, "generation")
  .transaction_create_directory(target, "prior")
  controls <- cellucid:::.export_transaction_control_paths(target)
  expect_true(dir.create(controls$journal_temp))
  expect_error(
    cellucid:::.recover_export_transaction(target),
    "non-linked|regular file"
  )
  expect_true(dir.exists(controls$journal_temp))
  expect_identical(
    readLines(file.path(target, "owner"), warn = FALSE),
    "prior"
  )
  unlink(parent, recursive = TRUE, force = TRUE)
})

test_that("unsafe target, stage, and backup nodes fail closed", {
  transaction_id <- "0123456789abcdef0123456789abcdef"
  cases <- c("target-symlink", "stage-file", "backup-symlink")

  for (case in cases) {
    parent <- .transaction_test_parent(case)
    target <- file.path(parent, "generation")
    reparse_nodes <- character()
    paths <- .transaction_write_journal(
      target,
      transaction_id,
      TRUE
    )
    outside <- file.path(parent, "outside")
    .transaction_create_directory(outside, "outside")

    if (case == "target-symlink") {
      link_created <- .create_test_directory_reparse(outside, target)
      if (!link_created) {
        .expect_transaction_cleanup(parent)
      }
      skip_if_not(link_created, "symbolic links unavailable")
      reparse_nodes <- c(reparse_nodes, target)
    } else {
      .transaction_create_directory(target, "prior")
    }
    if (case == "stage-file") {
      writeLines("unsafe stage", paths$stage)
    }
    if (case == "backup-symlink") {
      link_created <- .create_test_directory_reparse(
        outside,
        paths$backup
      )
      if (!link_created) {
        .expect_transaction_cleanup(parent)
      }
      skip_if_not(link_created, "symbolic links unavailable")
      reparse_nodes <- c(reparse_nodes, paths$backup)
    }
    outside_before <- .transaction_snapshot(outside)
    journal_before <- readBin(
      paths$journal,
      what = "raw",
      n = 512L
    )

    expect_error(
      cellucid:::.recover_export_transaction(target),
      "ordinary non-symbolic directory",
      info = case
    )
    expect_identical(.transaction_snapshot(outside), outside_before)
    expect_identical(
      readBin(paths$journal, what = "raw", n = 512L),
      journal_before
    )
    for (path in reparse_nodes) {
      .expect_transaction_reparse_node(path)
    }
    .expect_transaction_cleanup(
      parent,
      directory_reparse_nodes = reparse_nodes,
      after_reparse_cleanup = function() {
        expect_identical(
          .transaction_snapshot(outside),
          outside_before
        )
      }
    )
  }
})

test_that("the public writer never follows the final output symlink", {
  parent <- .transaction_test_parent("public_target_symlink")
  outside <- file.path(parent, "outside")
  .transaction_create_directory(outside, "outside")
  target <- file.path(parent, "generation")
  link_created <- .create_test_directory_reparse(outside, target)
  if (!link_created) {
    .expect_transaction_cleanup(parent)
  }
  skip_if_not(link_created, "directory reparse points unavailable")
  before <- .transaction_snapshot(outside)
  on.exit(
    .expect_transaction_cleanup(
      parent,
      directory_reparse_nodes = target,
      after_reparse_cleanup = function() {
        expect_identical(.transaction_snapshot(outside), before)
      }
    ),
    add = TRUE
  )

  expect_error(
    .transaction_prepare(target, "replacement", TRUE),
    "ordinary non-symbolic directory"
  )
  .expect_transaction_reparse_node(target)
  expect_identical(.transaction_snapshot(outside), before)
  .expect_transaction_controls_cleared(target)
})

test_that("publication validates journal and temporary-node ownership before rename", {
  parent <- .transaction_test_parent("pre_publish_validation")
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  target <- file.path(parent, "generation")
  .transaction_create_directory(target, "prior")
  transaction <- cellucid:::.begin_export_transaction(target, TRUE)
  writeLines("replacement", file.path(transaction$stage, "owner"))
  paths <- cellucid:::.export_transaction_paths(
    target,
    transaction$transaction_id
  )

  writeBin(
    charToRaw('{"format":"replaced"}\n'),
    paths$journal
  )
  expect_error(
    cellucid:::.publish_export_generation(
      transaction$stage,
      target,
      transaction$transaction_id,
      transaction$had_target
    ),
    "canonical|active transaction"
  )
  expect_identical(
    readLines(file.path(target, "owner"), warn = FALSE),
    "prior"
  )
  expect_true(dir.exists(transaction$stage))
  expect_false(file.exists(transaction$backup))

  writeBin(
    charToRaw(cellucid:::.serialize_export_transaction(
      transaction$transaction_id,
      TRUE
    )),
    paths$journal
  )
  writeLines("reappeared", paths$journal_temp)
  expect_error(
    cellucid:::.publish_export_generation(
      transaction$stage,
      target,
      transaction$transaction_id,
      transaction$had_target
    ),
    "temporary path reappeared"
  )
  expect_identical(
    readLines(file.path(target, "owner"), warn = FALSE),
    "prior"
  )
  expect_true(dir.exists(transaction$stage))
  expect_false(file.exists(transaction$backup))
})

test_that("recovery precedes the public force policy check", {
  parent <- .transaction_test_parent("recover_before_force")
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  target <- file.path(parent, "generation")
  transaction_id <- "0123456789abcdef0123456789abcdef"
  paths <- .transaction_write_journal(
    target,
    transaction_id,
    TRUE
  )
  .transaction_create_directory(paths$stage, "rejected-stage")
  .transaction_create_directory(paths$backup, "prior")

  expect_error(
    .transaction_prepare(target, "replacement", FALSE),
    "out_dir already exists.*force = TRUE"
  )
  expect_identical(
    readLines(file.path(target, "owner"), warn = FALSE),
    "prior"
  )
  expect_false(file.exists(paths$journal))
  expect_false(file.exists(paths$journal_temp))
  expect_false(file.exists(paths$stage))
  expect_false(file.exists(paths$backup))
})

test_that("ordinary failures recover every publication rename boundary", {
  phases <- c(
    "journal-published",
    "before-stage-publish",
    "prior-moved",
    "stage-published"
  )

  for (phase in phases) {
    parent <- .transaction_test_parent(paste0("ordinary_", phase))
    target <- file.path(parent, "generation")
    .transaction_prepare(target, "prior-generation", TRUE)
    canonical_target <- normalizePath(
      target,
      winslash = "/",
      mustWork = TRUE
    )
    real_rename <- cellucid:::.rename_export_path
    crashing_rename <- function(source, destination) {
      is_journal <- endsWith(
        destination,
        ".cellucid-transaction.json"
      )
      is_stage_publish <- startsWith(
        basename(source),
        paste0(
          ".",
          basename(canonical_target),
          ".cellucid-stage-"
        )
      ) && identical(destination, canonical_target)
      if (phase == "before-stage-publish" && is_stage_publish) {
        stop("synthetic before-stage-publish failure")
      }
      real_rename(source, destination)
      is_prior_move <- identical(source, canonical_target) &&
        startsWith(
          basename(destination),
          paste0(
            ".",
            basename(canonical_target),
            ".cellucid-backup-"
          )
        )
      if (phase == "journal-published" && is_journal) {
        stop("synthetic journal-published failure")
      }
      if (phase == "prior-moved" && is_prior_move) {
        stop("synthetic prior-moved failure")
      }
      if (phase == "stage-published" && is_stage_publish) {
        stop("synthetic stage-published failure")
      }
      invisible(destination)
    }

    testthat::with_mocked_bindings(
      expect_error(
        .transaction_prepare(
          target,
          "replacement-generation",
          TRUE
        ),
        paste0("synthetic ", phase, " failure")
      ),
      .rename_export_path = crashing_rename,
      .package = "cellucid"
    )
    expected <- if (phase == "stage-published") {
      "replacement-generation"
    } else {
      "prior-generation"
    }
    expect_identical(
      .transaction_dataset_id(target),
      expected,
      info = phase
    )
    .expect_transaction_controls_cleared(target)
    unlink(parent, recursive = TRUE, force = TRUE)
  }
})

test_that("abrupt process death recovers every publication rename boundary", {
  phases <- c(
    "journal-published",
    "before-stage-publish",
    "prior-moved",
    "stage-published"
  )
  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  library_expression <- paste(deparse(.libPaths()), collapse = " ")

  for (phase in phases) {
    parent <- .transaction_test_parent(paste0("abrupt_", phase))
    target <- file.path(parent, "generation")
    .transaction_prepare(target, "prior-generation", TRUE)
    canonical_target <- normalizePath(
      target,
      winslash = "/",
      mustWork = TRUE
    )
    script <- file.path(parent, "crash-export.R")
    writeLines(
      c(
        paste0(".libPaths(", library_expression, ")"),
        "library(cellucid)",
        "arguments <- commandArgs(trailingOnly = TRUE)",
        "target <- arguments[[1L]]",
        "phase <- arguments[[2L]]",
        "namespace <- asNamespace('cellucid')",
        "real_rename <- get('.rename_export_path', envir = namespace)",
        "crashing_rename <- function(source, destination) {",
        "  is_journal <- endsWith(destination, '.cellucid-transaction.json')",
        paste0(
          "  is_stage_publish <- startsWith(basename(source), ",
          "paste0('.', basename(target), '.cellucid-stage-')) && ",
          "identical(destination, target)"
        ),
        "  if (phase == 'before-stage-publish' && is_stage_publish) {",
        "    quit(save = 'no', status = 73L, runLast = FALSE)",
        "  }",
        "  real_rename(source, destination)",
        paste0(
          "  is_prior_move <- identical(source, target) && ",
          "startsWith(basename(destination), ",
          "paste0('.', basename(target), '.cellucid-backup-'))"
        ),
        "  if (phase == 'journal-published' && is_journal) {",
        "    quit(save = 'no', status = 73L, runLast = FALSE)",
        "  }",
        "  if (phase == 'prior-moved' && is_prior_move) {",
        "    quit(save = 'no', status = 73L, runLast = FALSE)",
        "  }",
        "  if (phase == 'stage-published' && is_stage_publish) {",
        "    quit(save = 'no', status = 73L, runLast = FALSE)",
        "  }",
        "  invisible(destination)",
        "}",
        "unlockBinding('.rename_export_path', namespace)",
        paste0(
          "assign('.rename_export_path', crashing_rename, ",
          "envir = namespace)"
        ),
        "lockBinding('.rename_export_path', namespace)",
        ".transaction_prepare <- function() {",
        "  latent <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)",
        "  obs <- data.frame(group = factor(c('A', 'B')))",
        "  umap1 <- matrix(c(0, 1), ncol = 1)",
        "  cellucid_prepare(",
        "    dataset_id = 'replacement-generation',",
        "    dataset_name = 'replacement-generation',",
        "    latent_space = latent,",
        "    obs = obs,",
        "    X_umap_1d = umap1,",
        "    out_dir = target,",
        "    centroid_min_points = 1,",
        "    force = TRUE,",
        "    obs_categorical_dtype = 'uint16'",
        "  )",
        "}",
        ".transaction_prepare()",
        "quit(save = 'no', status = 4L, runLast = FALSE)"
      ),
      script,
      useBytes = TRUE
    )
    output <- suppressWarnings(system2(
      rscript,
      c(
        "--vanilla",
        shQuote(script),
        shQuote(canonical_target),
        shQuote(phase)
      ),
      stdout = TRUE,
      stderr = TRUE,
      timeout = 30
    ))
    status <- attr(output, "status", exact = TRUE)
    if (is.null(status)) {
      status <- 0L
    }
    expect_identical(
      as.integer(status),
      73L,
      info = paste(phase, paste(output, collapse = "\n"))
    )

    expect_error(
      .transaction_prepare(target, "retry-generation", FALSE),
      "out_dir already exists.*force = TRUE",
      info = phase
    )
    expected <- if (phase == "stage-published") {
      "replacement-generation"
    } else {
      "prior-generation"
    }
    expect_identical(
      .transaction_dataset_id(target),
      expected,
      info = phase
    )
    .expect_transaction_controls_cleared(target)
    unlink(parent, recursive = TRUE, force = TRUE)
  }
})
