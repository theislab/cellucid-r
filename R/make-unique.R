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
    stop("join must be one non-empty string.", call. = FALSE)
  }
  if (is.factor(names)) {
    names <- as.character(names)
  }
  values <- if (is.list(names)) names else as.list(names)
  n <- length(values)
  if (!is.null(index) && length(index) != n) {
    stop(
      "index holds ", length(index), " names but the name list holds ", n, ".",
      call. = FALSE
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
  # byte for byte, exactly as Python compares str.
  taken <- new.env(hash = TRUE, parent = emptyenv(), size = max(1L, 2L * n))
  for (value in filled) {
    assign(value, TRUE, envir = taken)
  }
  seen <- new.env(hash = TRUE, parent = emptyenv(), size = max(1L, n))
  counters <- new.env(hash = TRUE, parent = emptyenv())
  result <- character(n)
  renamed <- list()
  for (position in seq_len(n)) {
    value <- filled[[position]]
    if (!exists(value, envir = seen, inherits = FALSE)) {
      assign(value, TRUE, envir = seen)
      result[[position]] <- value
      next
    }
    counter <- if (exists(value, envir = counters, inherits = FALSE)) {
      get(value, envir = counters, inherits = FALSE)
    } else {
      0L
    }
    repeat {
      counter <- counter + 1L
      candidate <- paste0(value, join, counter)
      if (!exists(candidate, envir = taken, inherits = FALSE)) {
        break
      }
    }
    assign(value, counter, envir = counters)
    assign(candidate, TRUE, envir = taken)
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

.plural <- function(count, singular, plural) {
  if (count == 1) singular else plural
}

# The one-line note that reports what making names unique changed, shared
# word for word with cellucid-python and the web app. Names are quoted with
# plain single quotes and never escaped, and counts carry thousands
# separators. NULL when nothing changed.
.unique_names_note <- function(filled_from_index, renamed, source) {
  parts <- character(0)
  if (filled_from_index > 0) {
    if (is.null(source)) {
      stop("The var index cannot fill names from itself.", call. = FALSE)
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
