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

# The positions of the first offending values, as the caller indexes them: plain
# indices for a vector, `[row, column]` for a matrix.
.format_non_finite_positions <- function(offending, dims) {
  positions <- which(offending)
  keep <- seq_len(min(length(positions), .non_finite_example_limit))
  shown <- positions[keep]
  if (is.null(dims)) {
    labels <- .format_count(shown)
  } else {
    rows <- ((shown - 1L) %% dims[1L]) + 1L
    columns <- ((shown - 1L) %/% dims[1L]) + 1L
    labels <- sprintf("[%s, %s]", .format_count(rows), .format_count(columns))
  }
  paste0(
    paste(labels, collapse = ", "),
    if (length(positions) > length(shown)) ", ..." else ""
  )
}

# `subject` is the name the caller used. `positions` is FALSE for a sparse
# matrix, whose stored non-zero entries have no index the caller would
# recognise; the counts are still exactly what was found.
.describe_non_finite <- function(values, subject, positions = TRUE) {
  flat <- as.vector(values)
  dims <- if (is.matrix(values)) dim(values) else NULL
  nan <- is.nan(flat)
  missing <- is.na(flat) & !nan
  infinite <- is.infinite(flat)
  positive <- infinite & !is.na(flat) & flat > 0
  negative <- infinite & !positive
  offending <- nan | missing | infinite
  counted <- c(
    if (any(nan)) paste(.format_count(sum(nan)), "NaN"),
    if (any(missing)) paste(.format_count(sum(missing)), "NA"),
    if (any(positive)) paste(.format_count(sum(positive)), "+Inf"),
    if (any(negative)) paste(.format_count(sum(negative)), "-Inf")
  )
  scope <- if (positions) "values" else "stored non-zero values"
  where <- if (positions) {
    paste0(
      " First affected position",
      if (sum(offending) == 1L) ": " else "s: ",
      .format_non_finite_positions(offending, dims),
      "."
    )
  } else {
    ""
  }
  paste0(
    subject, " cannot be published: of ", .format_count(length(flat)), " ",
    scope, ", ", paste(counted, collapse = ", "), ".", where,
    " Cellucid publishes finite values only, because a colour scale has no",
    " position for an infinity and no value for a missing one."
  )
}
