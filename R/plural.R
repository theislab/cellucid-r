# The noun a count takes: .plural(1, "cell") is "cell", .plural(2, "cell") is
# "cells". Every message that writes a count before a noun writes the noun with
# this, never as paste0(n, " cells"), which reads "1 cells" the day n is one
# (tests/testthat/test-count-noun-agreement.R). The count is written by the
# caller, so its own formatting stays as it was. A count that is not known
# (NA) takes the plural.
.plural <- function(count, singular, plural = paste0(singular, "s")) {
  if (isTRUE(count == 1)) singular else plural
}
