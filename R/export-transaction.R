# Publishing a generation is one step or none. A journal beside the output
# directory names the transaction and records whether a prior generation
# existed, so an interrupted export is resolved on the next call from what
# is on disk rather than from a guess: stage, backup, and target together
# decide whether the generation is committed or rolled back.

.export_transaction_format <- "cellucid-export-transaction"
.export_transaction_version <- 1L
.export_transaction_id_pattern <- "^[0-9a-f]{32}$"

.validate_export_transaction_id <- function(transaction_id) {
  if (
    !is.character(transaction_id) ||
      length(transaction_id) != 1L ||
      is.na(transaction_id) ||
      is.object(transaction_id) ||
      !grepl(
        .export_transaction_id_pattern,
        transaction_id,
        perl = TRUE,
        useBytes = TRUE
      )
  ) {
    stop("Export transaction identity is not canonical.", call. = FALSE)
  }
  transaction_id
}

.export_transaction_paths <- function(final_path, transaction_id) {
  transaction_id <- .validate_export_transaction_id(transaction_id)
  parent <- dirname(final_path)
  stem <- paste0(".", basename(final_path))
  list(
    journal = file.path(
      parent,
      paste0(stem, ".cellucid-transaction.json")
    ),
    journal_temp = file.path(
      parent,
      paste0(stem, ".cellucid-transaction.json.tmp")
    ),
    stage = file.path(
      parent,
      paste0(stem, ".cellucid-stage-", transaction_id)
    ),
    backup = file.path(
      parent,
      paste0(stem, ".cellucid-backup-", transaction_id)
    )
  )
}

.export_transaction_control_paths <- function(final_path) {
  paths <- .export_transaction_paths(final_path, strrep("0", 32L))
  paths[c("journal", "journal_temp")]
}

.serialize_export_transaction <- function(transaction_id, had_target) {
  transaction_id <- .validate_export_transaction_id(transaction_id)
  if (
    !is.logical(had_target) ||
      length(had_target) != 1L ||
      is.na(had_target) ||
      is.object(had_target)
  ) {
    stop(
      "Export transaction had_target must be exactly boolean.",
      call. = FALSE
    )
  }
  paste0(
    "{\"format\":\"",
    .export_transaction_format,
    "\",\"version\":",
    .export_transaction_version,
    ",\"transaction_id\":\"",
    transaction_id,
    "\",\"had_target\":",
    if (had_target) "true" else "false",
    "}\n"
  )
}

.read_export_transaction <- function(journal_path) {
  .require_export_regular_file(
    journal_path,
    "Export transaction journal"
  )
  journal_bytes <- readBin(
    journal_path,
    what = "raw",
    n = 513L
  )
  if (length(journal_bytes) > 512L) {
    stop(
      "Export transaction journal is unexpectedly large: ",
      journal_path,
      call. = FALSE
    )
  }
  if (
    length(journal_bytes) == 0L ||
      any(as.integer(journal_bytes) > 127L)
  ) {
    stop(
      "Export transaction journal is malformed: ",
      journal_path,
      call. = FALSE
    )
  }
  journal <- tryCatch(
    rawToChar(journal_bytes),
    error = function(error) NULL
  )
  if (is.null(journal)) {
    stop(
      "Export transaction journal is malformed: ",
      journal_path,
      call. = FALSE
    )
  }

  prefix <- paste0(
    "{\"format\":\"",
    .export_transaction_format,
    "\",\"version\":",
    .export_transaction_version,
    ",\"transaction_id\":\""
  )
  for (had_target in c(TRUE, FALSE)) {
    suffix <- paste0(
      "\",\"had_target\":",
      if (had_target) "true" else "false",
      "}\n"
    )
    expected_length <- nchar(prefix, type = "bytes") +
      32L +
      nchar(suffix, type = "bytes")
    if (
      nchar(journal, type = "bytes") == expected_length &&
        startsWith(journal, prefix) &&
        endsWith(journal, suffix)
    ) {
      transaction_id <- substr(
        journal,
        nchar(prefix, type = "chars") + 1L,
        nchar(prefix, type = "chars") + 32L
      )
      if (
        grepl(
          .export_transaction_id_pattern,
          transaction_id,
          perl = TRUE,
          useBytes = TRUE
        ) &&
          identical(
            journal_bytes,
            charToRaw(.serialize_export_transaction(
              transaction_id,
              had_target
            ))
          )
      ) {
        return(list(
          transaction_id = transaction_id,
          had_target = had_target
        ))
      }
    }
  }
  stop(
    "Export transaction journal is not canonical: ",
    journal_path,
    call. = FALSE
  )
}

.require_active_export_transaction <- function(
    journal_path,
    transaction_id,
    had_target
) {
  transaction <- .read_export_transaction(journal_path)
  if (
    !identical(transaction$transaction_id, transaction_id) ||
      !identical(transaction$had_target, had_target)
  ) {
    stop(
      "Export transaction journal does not describe the active transaction.",
      call. = FALSE
    )
  }
  invisible(transaction)
}

.discard_export_transaction_temp <- function(journal_temp) {
  if (.export_path_info(journal_temp)$kind == 0L) {
    return(invisible(journal_temp))
  }
  .remove_export_control_file(journal_temp)
}

.recover_export_transaction <- function(final_path) {
  controls <- .export_transaction_control_paths(final_path)
  .discard_export_transaction_temp(controls$journal_temp)
  if (.export_path_info(controls$journal)$kind == 0L) {
    return(invisible(final_path))
  }

  transaction <- .read_export_transaction(controls$journal)
  paths <- .export_transaction_paths(
    final_path,
    transaction$transaction_id
  )
  target_exists <- .require_export_directory_or_absent(
    final_path,
    "Export target"
  )
  stage_exists <- .require_export_directory_or_absent(
    paths$stage,
    "Staged export generation"
  )
  backup_exists <- .require_export_directory_or_absent(
    paths$backup,
    "Prior export generation"
  )
  state <- c(target_exists, stage_exists, backup_exists)

  if (transaction$had_target) {
    if (identical(state, c(TRUE, FALSE, FALSE))) {
      NULL
    } else if (identical(state, c(TRUE, TRUE, FALSE))) {
      .remove_export_tree(paths$stage)
    } else if (identical(state, c(FALSE, TRUE, TRUE))) {
      .rename_export_path(paths$backup, final_path)
      .remove_export_tree(paths$stage)
    } else if (identical(state, c(TRUE, FALSE, TRUE))) {
      .remove_export_tree(paths$backup)
    } else {
      stop(
        "Export transaction cannot be recovered without guessing whether ",
        "to commit or roll back: target/stage/backup state is ",
        paste(state, collapse = "/"),
        ".",
        call. = FALSE
      )
    }
  } else {
    if (identical(state, c(FALSE, FALSE, FALSE))) {
      NULL
    } else if (identical(state, c(FALSE, TRUE, FALSE))) {
      .remove_export_tree(paths$stage)
    } else if (identical(state, c(TRUE, FALSE, FALSE))) {
      NULL
    } else {
      stop(
        "Initial export transaction cannot be recovered without guessing ",
        "whether to commit or roll back: target/stage/backup state is ",
        paste(state, collapse = "/"),
        ".",
        call. = FALSE
      )
    }
  }

  .require_active_export_transaction(
    paths$journal,
    transaction$transaction_id,
    transaction$had_target
  )
  .remove_export_control_file(paths$journal)
  invisible(final_path)
}

.new_export_transaction_id <- function(final_path) {
  for (attempt in seq_len(128L)) {
    transaction_id <- .native_export_transaction_id()
    transaction_id <- .validate_export_transaction_id(transaction_id)
    paths <- .export_transaction_paths(final_path, transaction_id)
    if (
      .export_path_info(paths$stage)$kind == 0L &&
        .export_path_info(paths$backup)$kind == 0L
    ) {
      return(transaction_id)
    }
  }
  stop(
    "Could not allocate a unique export transaction identity.",
    call. = FALSE
  )
}

.write_export_transaction <- function(
    final_path,
    transaction_id,
    had_target
) {
  paths <- .export_transaction_paths(final_path, transaction_id)
  if (
    .export_path_info(paths$journal)$kind != 0L ||
      .export_path_info(paths$journal_temp)$kind != 0L
  ) {
    stop(
      "Export transaction control path is already occupied for ",
      final_path,
      ".",
      call. = FALSE
    )
  }
  journal_bytes <- charToRaw(.serialize_export_transaction(
    transaction_id,
    had_target
  ))
  status <- .native_export_write_journal(
    paths$journal_temp,
    journal_bytes
  )
  if (!identical(status, TRUE)) {
    stop(
      "Native export transaction journal write returned invalid state.",
      call. = FALSE
    )
  }
  .rename_export_path(paths$journal_temp, paths$journal)
  invisible(paths$journal)
}

.begin_export_transaction <- function(final_path, force) {
  had_target <- .require_export_directory_or_absent(
    final_path,
    "Export target"
  )
  if (had_target && !force) {
    stop(
      "out_dir already exists; set force = TRUE to replace the complete ",
      "generation.",
      call. = FALSE
    )
  }
  transaction_id <- .new_export_transaction_id(final_path)
  .write_export_transaction(
    final_path,
    transaction_id,
    had_target
  )
  paths <- .export_transaction_paths(final_path, transaction_id)
  if (
    !dir.create(
      paths$stage,
      recursive = FALSE,
      showWarnings = FALSE
    )
  ) {
    stop("Could not create the staged output directory.", call. = FALSE)
  }
  .fsync_export_directory(dirname(final_path))
  list(
    transaction_id = transaction_id,
    had_target = had_target,
    stage = paths$stage,
    backup = paths$backup
  )
}

.publish_export_generation <- function(
    stage,
    final_path,
    transaction_id,
    had_target
) {
  paths <- .export_transaction_paths(final_path, transaction_id)
  if (!identical(stage, paths$stage)) {
    stop(
      "Staged export path does not belong to the active transaction.",
      call. = FALSE
    )
  }
  .require_active_export_transaction(
    paths$journal,
    transaction_id,
    had_target
  )
  if (.export_path_info(paths$journal_temp)$kind != 0L) {
    stop(
      "Export transaction journal temporary path reappeared during publication.",
      call. = FALSE
    )
  }
  if (!.require_export_directory_or_absent(
    stage,
    "Staged export generation"
  )) {
    stop("Staged export directory is missing: ", stage, call. = FALSE)
  }

  if (had_target) {
    if (!.require_export_directory_or_absent(
      final_path,
      "Export target"
    )) {
      stop(
        "Prior export generation is missing: ",
        final_path,
        call. = FALSE
      )
    }
    if (.export_path_info(paths$backup)$kind != 0L) {
      stop(
        "Export backup path is already occupied: ",
        paths$backup,
        call. = FALSE
      )
    }
    .rename_export_path(final_path, paths$backup)
  } else if (.export_path_info(final_path)$kind != 0L) {
    stop(
      "Initial export target appeared during publication: ",
      final_path,
      call. = FALSE
    )
  }

  .rename_export_path(stage, final_path)
  if (had_target) {
    .remove_export_tree(paths$backup)
  }
  .require_active_export_transaction(
    paths$journal,
    transaction_id,
    had_target
  )
  .remove_export_control_file(paths$journal)
  invisible(final_path)
}
