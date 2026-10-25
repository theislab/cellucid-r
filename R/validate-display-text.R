# The rule for text the viewer draws verbatim: the dataset's name, its
# description, and its source and citation lines. The character classes, the
# character names, and the escaping all live in one place.

# The dataset's name, description, source and citation are drawn with the
# browser's default white-space handling, in the dataset picker and the
# dataset information panel. Characters that carry no glyph therefore change
# the value without changing the picture: "Atlas " and "Atlas" look identical
# on screen.
#
# All three classes below are rejected rather than repaired: trimming would
# rewrite text the caller never asked to change.
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
# Every character the two non-control rejected classes can hold, so an
# invisible defect is named on screen by its Unicode name.
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
  key <- sprintf("%04X", .code_points(character))
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
  codepoints <- .code_points(text)
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
    .abort(arg, " must be a non-empty string.")
  }
  defect <- .display_text_defect(value)
  if (!is.null(defect)) {
    .abort(
      arg,
      " is displayed verbatim, so it must not carry characters that have no ",
      "glyph: ",
      .escape_display_text(value),
      " ",
      defect,
      ". Cellucid does not remove them for you, because that would change ",
      "text you did not ask to change. Pass the exact text you want shown."
    )
  }
  value
}
