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
