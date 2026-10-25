# Text that is the same in every R session, whatever its native encoding.
#
# R keeps each string with its encoding (UTF-8, latin1, or the session's
# native one) and makes it native wherever the C level takes a string as
# plain bytes: a condition's text, a symbol, a string read from bytes. A
# session whose native encoding is UTF-8 -- Linux and macOS, and Windows from
# R 4.2 -- loses nothing there. A session whose native encoding is a code page
# -- Windows before R 4.2 (CP1252 in Western Europe and the Americas), or a
# Unix locale such as C or ISO-8859-1 -- writes every character the code page
# lacks as an escape such as "<U+2192>", and takes UTF-8 bytes held as native
# text for characters of the code page (an e-acute, two bytes in UTF-8, is
# read as two characters). The package supports R 4.1 on Windows, so text
# crosses each of those borders here.

# --- Conditions ---------------------------------------------------------------

# message(), warning() and stop() given strings pass them through gettext(),
# and a warning's or an error's text through C as well, and each makes it
# native: the note "1 duplicate name was made unique, e.g. 'A' -> 'A-1'." (with
# U+2192 for the arrow) reached every handler, and a test, as
# "'A' <U+2192> 'A-1'". A condition object is signalled as it is, so its text
# is exactly the package's, in UTF-8; only a console that cannot show a
# character shows its escape. The arguments are joined as message() joins
# them: each made character, then all pasted with nothing between.
.condition_text <- function(...) {
  enc2utf8(paste(unlist(lapply(list(...), as.character)), collapse = ""))
}

# message(...): a message the caller can muffle with suppressMessages().
.inform <- function(...) {
  message(simpleMessage(paste0(.condition_text(...), "\n")))
}

# warning(..., call. = FALSE).
.warn <- function(...) {
  warning(simpleWarning(.condition_text(...)))
}

# stop(..., call. = FALSE).
.abort <- function(...) {
  stop(simpleError(.condition_text(...)))
}

# --- Bytes and text -----------------------------------------------------------

# The UTF-8 bytes of a string, whatever its encoding.
.utf8_bytes <- function(text) {
  charToRaw(enc2utf8(text))
}

# The string UTF-8 `bytes` hold, marked UTF-8, or NA when they are not UTF-8
# or hold a nul (which no R string holds; rawToChar() would drop a trailing
# one). rawToChar() marks its result native, so a JSON document or a frame
# from the browser would otherwise be read in the session's code page.
.utf8_text <- function(bytes) {
  if (any(bytes == as.raw(0L))) return(NA_character_)
  .as_utf8(rawToChar(bytes))
}

# A string whose bytes are UTF-8 but which R holds as native text (what
# rawToChar() and URLdecode() return, httpuv's text frames), marked UTF-8; NA
# when the bytes are not UTF-8.
.as_utf8 <- function(text) {
  if (is.na(text) || !validUTF8(text)) return(NA_character_)
  Encoding(text) <- "UTF-8"
  text
}

# The Unicode code points of one string, whatever its encoding: utf8ToInt()
# reads the bytes as UTF-8, and a latin1 or native e-acute is one byte that
# is not UTF-8 (NA).
.code_points <- function(text) {
  utf8ToInt(enc2utf8(text))
}

# Strings in code-point order, as Python's sorted() and the web app order
# them, never by the locale's collation. A radix sort compares the bytes as
# R holds them, which is code-point order only for UTF-8: a latin1 e-acute
# (the byte 0xE9) sorted after a UTF-8 Greek alpha (0xCE 0xB1).
.sort_code_points <- function(text) {
  sort(enc2utf8(text), method = "radix")
}

# --- Environment keys ---------------------------------------------------------

# An environment is R's hash table, and its names are symbols, which R makes
# in the native encoding: in a session that is not UTF-8 a name holding a
# character the code page lacks becomes an escape such as "<U+03B1>", with a
# warning for every lookup, and then equals the name that is that escape. A
# key is the string's UTF-8 bytes taken as they are (marked native, so that
# nothing translates them), so two keys are equal exactly when the strings
# are, in every session.
.symbol_keys <- function(text) {
  keys <- enc2utf8(as.character(text))
  Encoding(keys) <- "unknown"
  keys
}
