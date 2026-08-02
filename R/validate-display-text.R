# The rule for text the viewer draws verbatim: dataset names, descriptions,
# source lines, category labels, and every field identifier. It is a
# cross-writer contract -- the cellucid Python package enforces the same
# classes in _display_text_defect() -- so the character classes, the
# character names, and the escaping all live in one place.

# Text the viewer prints verbatim -- a category label, a dataset name, a
# description, a source line -- is drawn with the browser's default white-space
# handling, in the legend, the field selector, and every exported figure.
# Characters that carry no glyph therefore change the value without changing
# the picture: "Liver " and "Liver" are two distinct categories that look
# identical on screen and in an exported SVG.
#
# All three classes below are rejected rather than repaired. Trimming would
# rewrite a scientific label the caller never asked to change, and would merge
# two distinct categories into one, silently moving cells between them. The
# Python exporter enforces the identical rule in _display_text_defect().
#
# R character vectors cannot hold U+0000, so the control class starts at U+0001.
.control_character_pattern <- "[\u0001-\u001f\u007f-\u009f]"
# Zero-width characters with no meaning of their own: ZERO WIDTH SPACE, WORD
# JOINER, and the byte-order mark a spreadsheet leaves at the front of a UTF-8
# CSV. U+200C and U+200D are deliberately absent: they join Indic, Persian, and
# emoji sequences, so banning them would reject real text.
.invisible_character_pattern <- "[\u200b\u2060\ufeff]"
# Every Unicode whitespace character that is not already a control character.
.edge_whitespace_class <- "[\u0020\u00a0\u1680\u2000-\u200a\u2028\u2029\u202f\u205f\u3000]"
.edge_whitespace_pattern <- paste0(
  "^",
  .edge_whitespace_class,
  "|",
  .edge_whitespace_class,
  "$"
)
.whitespace_run_pattern <- paste0(.edge_whitespace_class, "+")
# Every character the two non-control rejected classes can hold, so an
# invisible defect is named on screen exactly as the Python exporter names it.
.display_character_names <- c(
  "0020" = "U+0020 SPACE",
  "00A0" = "U+00A0 NO-BREAK SPACE",
  "1680" = "U+1680 OGHAM SPACE MARK",
  "2000" = "U+2000 EN QUAD",
  "2001" = "U+2001 EM QUAD",
  "2002" = "U+2002 EN SPACE",
  "2003" = "U+2003 EM SPACE",
  "2004" = "U+2004 THREE-PER-EM SPACE",
  "2005" = "U+2005 FOUR-PER-EM SPACE",
  "2006" = "U+2006 SIX-PER-EM SPACE",
  "2007" = "U+2007 FIGURE SPACE",
  "2008" = "U+2008 PUNCTUATION SPACE",
  "2009" = "U+2009 THIN SPACE",
  "200A" = "U+200A HAIR SPACE",
  "2028" = "U+2028 LINE SEPARATOR",
  "2029" = "U+2029 PARAGRAPH SEPARATOR",
  "202F" = "U+202F NARROW NO-BREAK SPACE",
  "205F" = "U+205F MEDIUM MATHEMATICAL SPACE",
  "3000" = "U+3000 IDEOGRAPHIC SPACE",
  "200B" = "U+200B ZERO WIDTH SPACE",
  "2060" = "U+2060 WORD JOINER",
  "FEFF" = "U+FEFF ZERO WIDTH NO-BREAK SPACE"
)

.describe_character <- function(character) {
  key <- sprintf("%04X", utf8ToInt(character))
  named <- unname(.display_character_names[key])
  if (!is.na(named)) {
    return(named)
  }
  # The C0 and C1 ranges carry no Unicode name at all, so say what they are.
  paste0("U+", key, " (control character)")
}

.escape_display_codepoint <- function(codepoint) {
  if (codepoint == 92L) {
    return("\\\\")
  }
  if (codepoint == 39L) {
    return("\\'")
  }
  if (codepoint == 9L) {
    return("\\t")
  }
  if (codepoint == 10L) {
    return("\\n")
  }
  if (codepoint == 13L) {
    return("\\r")
  }
  character <- intToUtf8(codepoint)
  if (
    codepoint == 32L ||
      !(
        grepl(.control_character_pattern, character, perl = TRUE) ||
          grepl(.invisible_character_pattern, character, perl = TRUE) ||
          grepl(.edge_whitespace_class, character, perl = TRUE)
      )
  ) {
    return(character)
  }
  if (codepoint <= 255L) {
    return(sprintf("\\x%02x", codepoint))
  }
  if (codepoint <= 65535L) {
    return(sprintf("\\u%04x", codepoint))
  }
  sprintf("\\U%08x", codepoint)
}

.escape_display_text <- function(text) {
  codepoints <- utf8ToInt(text)
  if (length(codepoints) == 0L) {
    return("''")
  }
  paste0(
    "'",
    paste(
      vapply(codepoints, .escape_display_codepoint, character(1)),
      collapse = ""
    ),
    "'"
  )
}

.display_text_defect <- function(text) {
  control <- regexpr(.control_character_pattern, text, perl = TRUE)
  if (control[[1L]] > 0L) {
    position <- control[[1L]]
    return(paste0(
      "contains ",
      .describe_character(substring(text, position, position))
    ))
  }
  hidden <- regexpr(.invisible_character_pattern, text, perl = TRUE)
  if (hidden[[1L]] > 0L) {
    position <- hidden[[1L]]
    return(paste0(
      "contains ",
      .describe_character(substring(text, position, position))
    ))
  }
  edge <- regexpr(.edge_whitespace_pattern, text, perl = TRUE)
  if (edge[[1L]] > 0L) {
    position <- edge[[1L]]
    return(paste0(
      if (position == 1L) "starts with" else "ends with",
      " ",
      .describe_character(substring(text, position, position))
    ))
  }
  NULL
}

.validate_display_text <- function(value, arg, allow_empty = FALSE) {
  value <- .validate_string(value, arg)
  if (!nzchar(value)) {
    if (allow_empty) {
      return(value)
    }
    stop(arg, " must be a non-empty string.", call. = FALSE)
  }
  defect <- .display_text_defect(value)
  if (!is.null(defect)) {
    stop(
      arg,
      " is displayed verbatim, so it must not carry characters that have no ",
      "glyph: ",
      .escape_display_text(value),
      " ",
      defect,
      ". Cellucid does not remove them for you, because that would change ",
      "text you did not ask to change. Pass the exact text you want shown.",
      call. = FALSE
    )
  }
  value
}

.validate_category_labels <- function(categories, key) {
  # Logical and numeric category labels carry no display text to inspect.
  if (!is.character(categories) || length(categories) == 0L) {
    return(invisible(categories))
  }
  suspect <- !nzchar(categories) |
    grepl(.control_character_pattern, categories, perl = TRUE) |
    grepl(.invisible_character_pattern, categories, perl = TRUE) |
    grepl(.edge_whitespace_pattern, categories, perl = TRUE)
  if (any(suspect)) {
    defects <- vapply(
      categories[suspect],
      function(label) {
        if (!nzchar(label)) {
          return("'' is empty")
        }
        paste0(.escape_display_text(label), " ", .display_text_defect(label))
      },
      character(1),
      USE.NAMES = FALSE
    )
    shown <- utils::head(defects, 10L)
    listed <- paste0("  - ", shown, collapse = "\n")
    remainder <- length(defects) - length(shown)
    if (remainder > 0L) {
      listed <- paste0(listed, "\n  - ... and ", remainder, " more")
    }
    counted <- if (length(defects) == 1L) {
      "1 label"
    } else {
      paste0(length(defects), " labels")
    }
    stop(
      "Categorical field '",
      key,
      "' has ",
      counted,
      " the viewer cannot show as written:\n",
      listed,
      "\nLabels are drawn verbatim in the legend, the field selector, and ",
      "exported figures, so these characters are invisible on screen while ",
      "still making the label a different value from the one it looks like. ",
      "Cellucid does not clean them for you: trimming a label rewrites your ",
      "annotation and can merge two categories into one, moving cells between ",
      "them. Clean the column and export again, for example:\n",
      "    obs[['",
      key,
      "']] <- factor(trimws(as.character(obs[['",
      key,
      "']]), whitespace = ",
      '"[\\\\h\\\\v]"',
      "))",
      call. = FALSE
    )
  }

  collapsed <- gsub(.whitespace_run_pattern, " ", categories, perl = TRUE)
  collision <- which(duplicated(collapsed))
  if (length(collision) > 0L) {
    later <- collision[[1L]]
    first <- match(collapsed[[later]], collapsed)
    stop(
      "Categorical field '",
      key,
      "' labels ",
      .escape_display_text(categories[[first]]),
      " and ",
      .escape_display_text(categories[[later]]),
      " are stored as different categories but are drawn identically: a run ",
      "of whitespace collapses to a single space in the legend, the field ",
      "selector, and exported figures. Rename one of them so the two ",
      "categories can be told apart, then export again.",
      call. = FALSE
    )
  }
  invisible(categories)
}
