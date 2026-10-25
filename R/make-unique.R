# Making a list of names unique: the one rule every Cellucid reader and writer
# applies to gene names. It is anndata's make_index_unique(), the rule behind
# var_names_make_unique(): an empty name first takes the index name at the
# same position, the first occurrence of a name keeps it, and each later
# occurrence becomes "<name>-<k>", where k is a per-name counter that only
# counts up and skips any candidate already present or already produced.
# cellucid-python's make_unique_names() and the web app's makeUniqueNames()
# implement the same rule, and the shared conformance fixture
# tests/testthat/fixtures/data_mapping_conformance.json pins all three.

# A failure of the rule carries a stable reason code -- missing_name,
# missing_index_name, or name_not_string -- as the condition's `code` field,
# the same codes the conformance fixture records.
.data_mapping_error <- function(code, message) {
  structure(
    class = c("cellucid_data_mapping_error", "error", "condition"),
    list(message = message, call = NULL, code = code)
  )
}

.is_missing_name <- function(value) {
  is.null(value) ||
    (length(value) == 1L && is.atomic(value) && is.na(value)) ||
    (is.character(value) && length(value) == 1L && !nzchar(value))
}

# Returns list(names, renamed, filled_from_index, filled, from_index):
# `renamed` holds c(from, to) for every renamed position in input order,
# `filled` each position's name after empty names took their index name and before
# uniqueness, and `from_index` which positions took their index name, so a caller can
# report on a subset of the positions.
.make_unique_names <- function(names, index = NULL, join = "-") {
  if (!is.character(join) || length(join) != 1L || is.na(join) || !nzchar(join)) {
    .abort("join must be one non-empty string.")
  }
  if (is.factor(names)) {
    names <- as.character(names)
  }
  # Distinct, non-empty strings are already unique: the rule changes nothing,
  # and a million cell names are checked in one hashed pass.
  if (is.character(names) && !anyNA(names) && all(nzchar(names)) && !anyDuplicated(names) &&
      (is.null(index) || length(index) == length(names))) {
    return(list(
      names = names, renamed = list(), filled_from_index = 0L,
      filled = names, from_index = logical(length(names))
    ))
  }
  values <- if (is.list(names)) names else as.list(names)
  n <- length(values)
  if (!is.null(index) && length(index) != n) {
    .abort(
      "index holds ", length(index), " names but the name list holds ", n, "."
    )
  }
  filled <- character(n)
  from_index <- logical(n)
  for (position in seq_len(n)) {
    value <- values[[position]]
    if (.is_missing_name(value)) {
      from_index[[position]] <- TRUE
      if (is.null(index)) {
        stop(.data_mapping_error(
          "missing_name",
          paste0(
            "Name at position ", position - 1L,
            " is empty and no index was supplied to take its name from."
          )
        ))
      }
      index_name <- index[[position]]
      if (
        !is.character(index_name) ||
          length(index_name) != 1L ||
          is.na(index_name) ||
          !nzchar(index_name)
      ) {
        stop(.data_mapping_error(
          "missing_index_name",
          paste0(
            "Name at position ", position - 1L,
            " is empty, and its index name is not a ",
            "non-empty string: ", .describe_index_name(index_name), "."
          )
        ))
      }
      filled[[position]] <- index_name
      next
    }
    if (!is.character(value) || length(value) != 1L) {
      stop(.data_mapping_error(
        "name_not_string",
        paste0(
          "Name at position ", position - 1L, " must be a string, got ",
          class(value)[[1L]], "."
        )
      ))
    }
    filled[[position]] <- value
  }

  # An environment is a hash set keyed by exact string; names are compared
  # code point for code point, exactly as Python compares str. Its keys are
  # the names' UTF-8 bytes (.symbol_keys()): in an R session that is not
  # UTF-8, a name holding a Greek alpha would otherwise become the symbol
  # "<U+03B1>", with a warning each time, and equal the name "<U+03B1>".
  keys <- .symbol_keys(filled)
  join_key <- .symbol_keys(join)
  taken <- new.env(hash = TRUE, parent = emptyenv(), size = max(1L, 2L * n))
  for (key in keys) {
    assign(key, TRUE, envir = taken)
  }
  seen <- new.env(hash = TRUE, parent = emptyenv(), size = max(1L, n))
  counters <- new.env(hash = TRUE, parent = emptyenv())
  result <- character(n)
  renamed <- list()
  for (position in seq_len(n)) {
    value <- filled[[position]]
    key <- keys[[position]]
    if (!exists(key, envir = seen, inherits = FALSE)) {
      assign(key, TRUE, envir = seen)
      result[[position]] <- value
      next
    }
    counter <- if (exists(key, envir = counters, inherits = FALSE)) {
      get(key, envir = counters, inherits = FALSE)
    } else {
      0L
    }
    repeat {
      counter <- counter + 1L
      candidate_key <- paste0(key, join_key, counter)
      if (!exists(candidate_key, envir = taken, inherits = FALSE)) {
        break
      }
    }
    assign(key, counter, envir = counters)
    assign(candidate_key, TRUE, envir = taken)
    candidate <- paste0(value, join, counter)
    result[[position]] <- candidate
    renamed[[length(renamed) + 1L]] <- c(value, candidate)
  }
  list(
    names = result,
    renamed = renamed,
    filled_from_index = sum(from_index),
    filled = filled,
    from_index = from_index
  )
}

# How an index name that cannot be taken is quoted back: a string in single
# quotes, anything else as R prints it.
.describe_index_name <- function(value) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    return(paste0("'", value, "'"))
  }
  paste(format(value), collapse = ", ")
}

# The one-line note that reports what making names unique changed, shared
# word for word with cellucid-python and the web app. Names are quoted with
# plain single quotes and never escaped, and counts carry thousands
# separators. NULL when nothing changed.
.unique_names_note <- function(filled_from_index, renamed, source) {
  parts <- character(0)
  if (filled_from_index > 0) {
    if (is.null(source)) {
      .abort("The var index cannot fill names from itself.")
    }
    parts <- c(parts, paste0(
      .format_count(filled_from_index),
      " ",
      .plural(
        filled_from_index,
        paste0("gene has no '", source, "' value and uses its var index name."),
        paste0("genes have no '", source, "' value and use their var index name.")
      )
    ))
  }
  if (length(renamed) > 0L) {
    first <- renamed[[1L]]
    parts <- c(parts, paste0(
      .format_count(length(renamed)),
      " ",
      .plural(
        length(renamed),
        "duplicate name was made unique",
        "duplicate names were made unique"
      ),
      ", e.g. '", first[[1L]], "' \u2192 '", first[[2L]], "'."
    ))
  }
  if (length(parts) == 0L) {
    return(NULL)
  }
  enc2utf8(paste(parts, collapse = " "))
}
