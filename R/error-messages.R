# How this package shows a set of values back to the caller. Every message
# that prints one goes through here, on every axis and in every subsystem,
# so one reader learns one list syntax.

# A message that shows the reader a set of values shows it as a list, not as a
# sentence that happens to continue after a comma: without a boundary,
# "columns not in obs: a, b. Available columns: x" cannot be read back as two
# lists. The list is written the way the caller writes one in R -- c("a", "b"),
# c(0, 1) -- so the message quotes back a value that can be pasted straight
# into the failing call. cellucid-python has the same reason to print a set as
# a list and prints its own language's, ['a', 'b'] and [0, 1]; the boundary is
# shared, the syntax inside it belongs to the reader's language. Every message
# in this package that shows a set of values uses this one function.
.format_value_list <- function(values) {
  if (is.character(values)) {
    items <- sprintf('"%s"', values)
  } else {
    items <- format(values, trim = TRUE)
  }
  paste0("c(", paste(items, collapse = ", "), ")")
}

# Why one continuous payload cannot be published, counted rather than named.
#
# `x must contain only finite values.` is true and stops there. One NaN in
# eighteen million cells and every value NaN produce that same sentence, and
# they want opposite responses: the first is a cell to look at, the second is
# the wrong matrix. cellucid-python's `continuous_payload_diagnosis` module
# counts the offenders and names the first few positions; this is the same
# report in R, so a dataset refused by either writer is refused with the same
# information.
#
# Counting is only reached once the cheap `anyNA(x) || any(!is.finite(x))` scan
# has already decided the payload is unusable, so a healthy column pays nothing
# and the working set is one the caller already had in memory.
.non_finite_example_limit <- 5L

.format_count <- function(count) {
  format(count, big.mark = ",", trim = TRUE, scientific = FALSE)
}

# Counting and locating in one bounded pass.
#
# The obvious form -- `is.nan(flat)`, `is.na(flat)`, `is.infinite(flat)` and
# `which(offending)` over the whole input -- allocates six logical vectors as
# long as the input plus one holding every offending position. The input that
# most often reaches here is a dense matrix, so explaining a failure would cost
# several times the memory of the thing that failed, and the diagnostic would
# run out of memory instead of printing the diagnosis.
#
# `sum()` is the same problem at the other end. On a logical vector it returns
# an integer, so past 2^31 - 1 offenders it returns NA with a warning -- and
# `if (sum(offending) == 1L)` then aborts on the missing value, replacing a
# useful message with `argument is not interpretable as logical`. A dense
# matrix with that many entries is ordinary: three million cells by a thousand
# genes is already past it.
#
# One pass over bounded chunks answers both. Counts accumulate in double, exact
# to 2^53, and only the first few positions are ever kept.
.non_finite_scan_chunk_size <- 1048576L

.scan_non_finite <- function(
    values,
    limit,
    .chunk_size = .non_finite_scan_chunk_size
) {
  n <- length(values)
  counts <- c(nan = 0, missing = 0, positive = 0, negative = 0)
  positions <- numeric(0)
  start <- 1
  while (start <= n) {
    end <- min(n, start + .chunk_size - 1)
    # Single-index selection on a matrix is linear indexing and yields a plain
    # vector, so no caller has to flatten its input to be diagnosed.
    chunk <- values[start:end]
    chunk_nan <- is.nan(chunk)
    chunk_missing <- is.na(chunk) & !chunk_nan
    chunk_infinite <- is.infinite(chunk)
    chunk_positive <- chunk_infinite & chunk > 0
    chunk_negative <- chunk_infinite & !chunk_positive
    counts[["nan"]] <- counts[["nan"]] + sum(chunk_nan)
    counts[["missing"]] <- counts[["missing"]] + sum(chunk_missing)
    counts[["positive"]] <- counts[["positive"]] + sum(chunk_positive)
    counts[["negative"]] <- counts[["negative"]] + sum(chunk_negative)
    if (length(positions) < limit) {
      found <- which(chunk_nan | chunk_missing | chunk_infinite)
      if (length(found) > 0L) {
        keep <- seq_len(min(length(found), limit - length(positions)))
        positions <- c(positions, found[keep] + (start - 1))
      }
    }
    start <- end + 1
  }
  list(counts = counts, positions = positions, length = n)
}

# The positions of the first offending values, as the caller indexes them: plain
# indices for a vector, `[row, column]` for a matrix. `shown` is already limited
# by the scan; `total` is how many there were, so the list can say it is one.
.format_non_finite_positions <- function(shown, total, dims) {
  if (is.null(dims)) {
    labels <- .format_count(shown)
  } else {
    rows <- ((shown - 1) %% dims[1L]) + 1
    columns <- ((shown - 1) %/% dims[1L]) + 1
    labels <- sprintf("[%s, %s]", .format_count(rows), .format_count(columns))
  }
  paste0(
    paste(labels, collapse = ", "),
    if (total > length(shown)) ", ..." else ""
  )
}

# `subject` is the name the caller used. `positions` is FALSE for a sparse
# matrix, whose stored non-zero entries have no index the caller would
# recognise; the counts are still exactly what was found.
.describe_non_finite <- function(values, subject, positions = TRUE) {
  dims <- if (is.matrix(values)) dim(values) else NULL
  scan <- .scan_non_finite(values, .non_finite_example_limit)
  counts <- scan$counts
  offending <- sum(counts)
  counted <- c(
    if (counts[["nan"]] > 0) paste(.format_count(counts[["nan"]]), "NaN"),
    if (counts[["missing"]] > 0) paste(.format_count(counts[["missing"]]), "NA"),
    if (counts[["positive"]] > 0) paste(.format_count(counts[["positive"]]), "+Inf"),
    if (counts[["negative"]] > 0) paste(.format_count(counts[["negative"]]), "-Inf")
  )
  scope <- if (positions) "values" else "stored non-zero values"
  where <- if (positions) {
    paste0(
      " First affected position",
      if (offending == 1) ": " else "s: ",
      .format_non_finite_positions(scan$positions, offending, dims),
      "."
    )
  } else {
    ""
  }
  paste0(
    subject, " cannot be published: of ", .format_count(scan$length), " ",
    scope, ", ", paste(counted, collapse = ", "), ".", where,
    " Cellucid publishes finite values only, because a colour scale has no",
    " position for an infinity and no value for a missing one."
  )
}
