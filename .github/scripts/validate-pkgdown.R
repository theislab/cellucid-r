expected <- c(
  "url: https://theislab.github.io/cellucid-r/",
  "template:",
  "  bootstrap: 5"
)
actual <- readLines("_pkgdown.yml", warn = FALSE, encoding = "UTF-8")
if (!identical(actual, expected)) {
  stop(
    "_pkgdown.yml must declare the canonical Cellucid R site URL and Bootstrap 5.",
    call. = FALSE
  )
}
