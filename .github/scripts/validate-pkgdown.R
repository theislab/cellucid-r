expected <- c(
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
  "      href: https://www.cellucid.com",
  "articles:",
  "  - title: Start here",
  "    contents:",
  "      - installation",
  "      - cellucid"
)
actual <- readLines("_pkgdown.yml", warn = FALSE, encoding = "UTF-8")
if (!identical(actual, expected)) {
  stop(
    paste(
      "_pkgdown.yml must declare the canonical Bootstrap 5 site,",
      "its inline favicon, Installation page, and navigation with only search",
      "and the Web app."
    ),
    call. = FALSE
  )
}
