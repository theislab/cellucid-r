# Every line of R this package ships or runs under R CMD check works on the
# oldest R that DESCRIPTION's Depends admits.
#
# The suite usually runs on a newer R, where a base function added since (such
# as tools::sha256sum() or %||%) or newer syntax (the pipe placeholder _) works
# and passes here, and then fails on the oldest supported R, in CI or for a
# user. The code is read token by token, and every base-package function it
# calls must be one that R exported: fixtures/r-<major>.<minor>-api.txt lists
# them, written by that R. A function the package or its tests define is
# theirs, whatever its name.

cellucid_minimum_r <- function() {
  depends <- read.dcf(file.path(cellucid_repository_root(), "DESCRIPTION"))[1L, "Depends"]
  version <- regmatches(depends, regexec("\\bR \\(>= ([0-9]+\\.[0-9]+)(\\.[0-9]+)?\\)", depends))[[1L]]
  if (length(version) < 2L) stop("DESCRIPTION's Depends names no oldest R.", call. = FALSE)
  version[[2L]]
}

# The R code under check: R/, the tests, inst/CITATION, and the vignettes'
# code chunks. Each element is the code with its name.
cellucid_r_sources <- function() {
  root <- cellucid_repository_root()
  files <- c(
    list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE),
    file.path(root, "tests", "testthat.R"),
    list.files(testthat::test_path(), pattern = "\\.R$", full.names = TRUE),
    file.path(root, "inst", "CITATION")
  )
  sources <- lapply(files, function(path) readLines(path, warn = FALSE, encoding = "UTF-8"))
  names(sources) <- substring(normalizePath(files, winslash = "/"), nchar(root) + 2L)
  for (path in list.files(file.path(root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE)) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    # Lines outside R chunks are blanked, so line numbers stay the file's.
    code <- character(length(lines))
    block <- FALSE
    r_chunk <- FALSE
    for (i in seq_along(lines)) {
      if (grepl("^\\s*```", lines[[i]])) {
        r_chunk <- !block && grepl("^\\s*```+\\s*\\{r[ ,}]", lines[[i]])
        block <- !block
      } else if (r_chunk) {
        code[[i]] <- lines[[i]]
      }
    }
    sources[[paste0("vignettes/", basename(path))]] <- code
  }
  sources
}

cellucid_parse_tokens <- function(lines, name) {
  data <- utils::getParseData(parse(text = lines, keep.source = TRUE, srcfile = srcfilecopy(name, lines)))
  # A file with no code at all (a vignette with no R chunk) has no parse data.
  if (is.null(data)) return(data.frame(line1 = integer(), col1 = integer(), token = character(), text = character()))
  data <- data[data$terminal, c("line1", "col1", "token", "text")]
  data[order(data$line1, data$col1), ]
}

test_that("the R API this package uses existed in the oldest R it supports", {
  minimum <- cellucid_minimum_r()
  fixture <- testthat::test_path("fixtures", paste0("r-", minimum, "-api.txt"))
  if (!file.exists(fixture)) {
    stop(
      "DESCRIPTION's Depends names R ", minimum, " but fixtures/r-", minimum,
      "-api.txt is missing; write it in R ", minimum, " as fixtures/r-4.1-api.txt says.",
      call. = FALSE
    )
  }
  api <- grep("^#", readLines(fixture, warn = FALSE, encoding = "UTF-8"), value = TRUE, invert = TRUE)
  base_packages <- unique(sub("::.*", "", api))
  expect_true(all(c("base", "stats", "tools", "utils") %in% base_packages))
  known <- sub("^[^:]+::", "", api)

  sources <- cellucid_r_sources()
  expect_true(all(c("R/cellucid_prepare.R", "inst/CITATION", "vignettes/cellucid.Rmd") %in% names(sources)))
  tokens <- Map(cellucid_parse_tokens, sources, names(sources))

  # Names the package and its tests define: assigned, function arguments, or
  # in the package namespace.
  own <- ls(asNamespace("cellucid"), all.names = TRUE)
  for (data in tokens) {
    target <- data$token %in% c("SYMBOL", "SPECIAL", "STR_CONST") &
      c(data$token[-1L], "") %in% c("LEFT_ASSIGN", "EQ_ASSIGN")
    own <- c(own, data$text[target], data$text[data$token == "SYMBOL_FORMALS"])
  }
  own <- unique(gsub("^[`'\"]|[`'\"]$", "", own))

  # What the running R offers from its base packages, by name.
  offered <- unique(unlist(lapply(base_packages, function(package) {
    if (package == "base") ls(baseenv(), all.names = TRUE) else getNamespaceExports(package)
  })))

  newer <- character()
  for (name in names(tokens)) {
    data <- tokens[[name]]
    found <- function(rows, what) sprintf("%s:%d %s", name, data$line1[rows], what)
    # The pipe placeholder _ is R 4.2 syntax.
    if (package_version(minimum) < "4.2") {
      placeholder <- which(data$token == "PLACEHOLDER")
      newer <- c(newer, found(placeholder, rep("the pipe placeholder _", length(placeholder))))
    }
    # package::name, for the base packages.
    qualified <- which(data$token == "NS_GET")
    qualified <- qualified[data$text[qualified - 1L] %in% base_packages]
    called <- sprintf("%s::%s", data$text[qualified - 1L], data$text[qualified + 1L])
    missing <- !called %in% api
    newer <- c(newer, found(qualified[missing] + 1L, called[missing]))
    # A bare call of a base-package function the package does not define (not
    # a namespaced one, and not an element called through $ or @).
    bare <- which(data$token %in% c("SYMBOL_FUNCTION_CALL", "SPECIAL") &
      !c("", data$token[-nrow(data)]) %in% c("NS_GET", "NS_GET_INT", "'$'", "'@'"))
    text <- data$text[bare]
    missing <- text %in% offered & !text %in% own & !text %in% known
    newer <- c(newer, found(bare[missing], text[missing]))
  }
  expect_identical(newer, character(), info = paste0("needs a newer R than ", minimum))
})
