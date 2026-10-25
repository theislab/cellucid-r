# The pkgdown site contract, checked before the site is built.
#
# The site header is fixed line for line: the canonical URL, Bootstrap 5, the
# inline favicon, and a navigation bar with the Installation page and, on the
# right, only search and the Web app. The article and reference indexes after
# it follow the repository: every vignette is an article (Installation, then
# the package overview, first) and every help page is in the reference, each
# exactly once, so adding a vignette or a help page means listing it here and
# nothing else.

contract <- c(
  "url: https://theislab.github.io/cellucid-r/",
  "template:",
  "  bootstrap: 5",
  "  includes:",
  "    in_header: >-",
  paste0(
    "      <link rel=\"icon\" type=\"image/svg+xml\" ",
    "href=\"data:image/svg+xml;base64,",
    "PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9",
    "IjAgMCA2NCA2NCI+PGNpcmNsZSBjeD0iMzIiIGN5PSIzMiIgcj0iMjciIGZpbGw9",
    "IiNmOGY5ZmEiIHN0cm9rZT0iIzM3NDE1MSIgc3Ryb2tlLXdpZHRoPSI0Ii8+PGNp",
    "cmNsZSBjeD0iMzIiIGN5PSIzMiIgcj0iMTAiIGZpbGw9IiMzNzQxNTEiLz48cGF0",
    "aCBkPSJNMTYgMzVjOCAxMyAyNyAxMyAzNCAwIiBmaWxsPSJub25lIiBzdHJva2U9",
    "IiM2YjcyODAiIHN0cm9rZS13aWR0aD0iNCIgc3Ryb2tlLWxpbmVjYXA9InJvdW5k",
    "Ii8+PC9zdmc+\">"
  ),
  "navbar:",
  "  structure:",
  "    left: [home, installation, reference, articles, news]",
  "    right: [search, webapp]",
  "  components:",
  "    home:",
  "      text: Home",
  "      href: index.html",
  "    installation:",
  "      text: Installation",
  "      href: articles/installation.html",
  "    webapp:",
  "      text: Web app",
  "      href: https://www.cellucid.com"
)

problems <- character()
actual <- readLines("_pkgdown.yml", warn = FALSE, encoding = "UTF-8")
if (!identical(utils::head(actual, length(contract)), contract)) {
  problems <- c(
    problems,
    paste(
      "_pkgdown.yml must declare the canonical Bootstrap 5 site,",
      "its inline favicon, Installation page, and navigation with only search",
      "and the Web app."
    )
  )
}
indexes <- actual[-seq_along(contract)]

# One index: its key, then groups of "  - title: <title>", "    contents:"
# and "      - <entry>" lines. Returns the entries in order, or NULL when a
# line is not of that shape.
index_entries <- function(lines, key) {
  if (length(lines) < 4L || !identical(lines[[1L]], paste0(key, ":"))) return(NULL)
  shape <- c("^  - title: \\S.*$", "^    contents:$", "^      - \\S+$")
  body <- lines[-1L]
  kind <- vapply(body, function(line) {
    match <- which(vapply(shape, grepl, logical(1L), x = line))
    if (length(match) == 1L) match else NA_integer_
  }, integer(1L), USE.NAMES = FALSE)
  previous <- c(3L, utils::head(kind, -1L))
  follows <- (kind == 1L & previous == 3L) | (kind == 2L & previous == 1L) |
    (kind == 3L & previous %in% c(2L, 3L))
  if (anyNA(kind) || !all(follows) || kind[[length(kind)]] != 3L) return(NULL)
  sub("^      - ", "", body[kind == 3L])
}

split_at <- match("reference:", indexes)
articles <- if (!is.na(split_at)) index_entries(indexes[seq_len(split_at - 1L)], "articles")
reference <- if (!is.na(split_at)) index_entries(indexes[split_at:length(indexes)], "reference")
if (is.null(articles) || is.null(reference)) {
  problems <- c(
    problems,
    paste(
      "_pkgdown.yml must end with an articles index and then a reference index,",
      "each a list of '- title:' groups with their 'contents:'."
    )
  )
} else {
  vignettes <- sub("\\.Rmd$", "", list.files("vignettes", pattern = "\\.Rmd$"))
  first <- c("installation", "cellucid")
  if (!identical(utils::head(articles, 2L), first) || anyDuplicated(articles) > 0L ||
      !setequal(articles, vignettes)) {
    problems <- c(
      problems,
      paste0(
        "The articles index must list installation, then cellucid, then every other ",
        "vignette once (it lists ", paste(articles, collapse = ", "),
        "; vignettes/ holds ", paste(vignettes, collapse = ", "), ")."
      )
    )
  }
  topics <- vapply(list.files("man", pattern = "\\.Rd$", full.names = TRUE), function(path) {
    lines <- readLines(path, warn = FALSE)
    name <- regmatches(lines, regexpr("(?<=^\\\\name\\{)[^}]+", lines, perl = TRUE))
    if (length(name) != 1L) stop(path, " has no single \\name{}.", call. = FALSE)
    name
  }, character(1L), USE.NAMES = FALSE)
  wrong <- list(
    "not listed" = setdiff(topics, reference),
    "not a help page" = setdiff(reference, topics),
    "listed twice" = unique(reference[duplicated(reference)])
  )
  wrong <- wrong[lengths(wrong) > 0L]
  if (length(wrong) > 0L) {
    problems <- c(
      problems,
      paste0(
        "The reference index must list every help page in man/ once (",
        paste(names(wrong), vapply(wrong, paste, character(1L), collapse = ", "), sep = ": ", collapse = "; "),
        ")."
      )
    )
  }
}

if (length(problems) > 0L) stop(paste(problems, collapse = "\n"), call. = FALSE)
