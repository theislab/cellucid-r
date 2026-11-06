# JSON as the web app writes it. Two things are compared with the viewer's own
# bytes, so their text must be exactly JavaScript's:
#
# - an action command (R/viewer-protocol.R): JSON.stringify, compact, keys in
#   the order given -- the protocol fixture's binary actions are matched byte
#   for byte;
# - a session value's content hash (R/session-file.R): the SHA-256 of its
#   canonical text (app/history/value-pool.js canonical()), keys sorted by
#   their UTF-16 code units, a typed array as {"$array":"<digest>"}, a number
#   JSON cannot hold as {"$number":"NaN"|"Infinity"|"-Infinity"|"-0"}.
#
# Numbers are written as JavaScript's Number::toString writes them (shortest
# round-trip digits, "1e+21", "1.5e-7", "5"; src/number_text.c), and strings
# as JSON.stringify escapes them: a quote, a backslash and the control
# characters, nothing else (no "/", no non-ASCII).
#
# A JSON value in R: NULL is null, a named list an object (its names the keys,
# in order; an empty one is structure(list(), names = character())), an
# unnamed list an array, a length-one logical, number or string a scalar, and
# a vector of any other length (or one wrapped in I()) an array of scalars.

# Strings as JSON string literals, as JSON.stringify writes them.
.js_strings <- function(text) {
  if (length(text) == 0L) return(character())
  text <- enc2utf8(as.character(text))
  if (anyNA(text)) .abort("JSON holds no NA text.")
  plain <- !grepl("[\"\\\\\001-\037]", text, useBytes = TRUE)
  out <- paste0("\"", text, "\"")
  if (!all(plain)) {
    out[!plain] <- vapply(text[!plain], .js_escaped_string, "", USE.NAMES = FALSE)
  }
  Encoding(out) <- "UTF-8"
  out
}

.js_escapes <- c(`8` = "\\b", `9` = "\\t", `10` = "\\n", `12` = "\\f", `13` = "\\r", `34` = "\\\"", `92` = "\\\\")

.js_escaped_string <- function(text) {
  points <- utf8ToInt(text)
  pieces <- intToUtf8(points, multiple = TRUE)
  special <- points < 32L | points == 34L | points == 92L
  for (at in which(special)) {
    escape <- .js_escapes[as.character(points[[at]])]
    pieces[[at]] <- if (is.na(escape)) sprintf("\\u%04x", points[[at]]) else escape
  }
  paste0("\"", paste(pieces, collapse = ""), "\"")
}

# Numbers as JSON.stringify writes them; JSON holds finite ones only.
.js_numbers <- function(numbers) {
  if (length(numbers) == 0L) return(character())
  if (is.integer(numbers)) {
    if (anyNA(numbers)) .abort("JSON holds no NA.")
    return(as.character(numbers))
  }
  numbers <- as.double(numbers)
  if (any(!is.finite(numbers))) .abort("JSON holds finite numbers only, not NA, NaN or an infinity.")
  .native_number_text(numbers)
}

# The scalars of an atomic vector as JSON texts.
.js_scalars <- function(value) {
  if (is.factor(value)) value <- as.character(value)
  if (is.character(value)) return(.js_strings(value))
  if (is.logical(value)) {
    if (anyNA(value)) .abort("JSON holds no NA.")
    return(ifelse(value, "true", "false"))
  }
  if (is.numeric(value)) return(.js_numbers(value))
  .abort("JSON holds text, numbers, TRUE/FALSE, lists and NULL; not ", class(value)[[1L]], ".")
}

# The texts of a list's items: its scalars written together, one call a
# type (a list of thousands of numbers is one call to src/number_text.c),
# every other item by `write`.
.js_texts <- function(items, write) {
  if (length(items) == 0L) return(character())
  scalar <- vapply(items, function(item) is.atomic(item) && length(item) == 1L && !inherits(item, "AsIs") &&
                     !inherits(item, "cellucid_typed_array"), TRUE, USE.NAMES = FALSE)
  out <- character(length(items))
  if (any(!scalar)) out[!scalar] <- vapply(items[!scalar], write, "", USE.NAMES = FALSE)
  if (any(scalar)) {
    at <- which(scalar)
    kinds <- vapply(items[at], function(item) if (is.factor(item)) "factor" else typeof(item), "", USE.NAMES = FALSE)
    for (kind in unique(kinds)) {
      same <- at[kinds == kind]
      values <- if (identical(kind, "factor")) {
        vapply(items[same], as.character, "", USE.NAMES = FALSE)
      } else {
        unlist(items[same], use.names = FALSE)
      }
      out[same] <- .js_scalars(values)
    }
  }
  out
}

# One JSON value as compact text, keys in the order given (JSON.stringify).
.js_json <- function(value) {
  if (is.null(value)) return("null")
  if (inherits(value, "cellucid_typed_array")) {
    .abort("A typed array is carried in a binary frame's body, not in JSON.")
  }
  if (is.list(value) && !is.data.frame(value)) {
    keys <- names(value)
    if (is.null(keys)) {
      return(paste0("[", paste(.js_texts(value, .js_json), collapse = ","), "]"))
    }
    if (anyDuplicated(keys)) .abort("A JSON object names a key twice.")
    if (length(value) == 0L) return("{}")
    parts <- paste0(.js_strings(keys), ":", .js_texts(value, .js_json))
    return(paste0("{", paste(parts, collapse = ","), "}"))
  }
  if (is.atomic(value) && !is.null(value)) {
    scalars <- .js_scalars(value)
    if (length(value) == 1L && !inherits(value, "AsIs")) return(scalars)
    return(paste0("[", paste(scalars, collapse = ","), "]"))
  }
  .abort("JSON holds text, numbers, TRUE/FALSE, lists and NULL; not ", class(value)[[1L]], ".")
}

# --- Canonical text of a stored value --------------------------------------------

# Strings in JavaScript's order (Array.prototype.sort: UTF-16 code units).
# UTF-8 bytes order code points, which is UTF-16's order except that a
# character past U+FFFF (a surrogate pair, D800-DFFF) sorts before U+E000 to
# U+FFFF: those are ordered by their units written as hexadecimal.
.utf16_order <- function(keys) {
  keys <- enc2utf8(keys)
  # A character from U+E000 on is the only one whose UTF-8 lead byte is 0xEE
  # or more.
  wide <- vapply(keys, function(key) any(charToRaw(key) >= as.raw(0xEE)), TRUE, USE.NAMES = FALSE)
  if (!any(wide)) return(order(keys, method = "radix"))
  units <- vapply(keys, function(key) {
    points <- utf8ToInt(key)
    astral <- points > 0xFFFF
    if (any(astral)) {
      offset <- points - 0x10000
      high <- 0xD800 + offset %/% 1024
      low <- 0xDC00 + offset %% 1024
      points <- unlist(lapply(seq_along(points), function(at) {
        if (astral[[at]]) c(high[[at]], low[[at]]) else points[[at]]
      }))
    }
    paste(sprintf("%04x", points), collapse = "")
  }, "", USE.NAMES = FALSE)
  order(units, method = "radix")
}

# The canonical text of a value's skeleton as a session file stores it (its
# placeholders already in their canonical form): keys sorted, compact. A part
# a skeleton names ({"$part": <name>}) is the canonical text of the value it
# stands for, `part_text(name)` (the reader's, memoized: one part's text is
# the same wherever it is named).
.canonical_text <- function(node, part_text = NULL) {
  if (is.null(node)) return("null")
  if (is.list(node)) {
    keys <- names(node)
    write <- function(child) .canonical_text(child, part_text)
    if (is.null(keys)) {
      return(paste0("[", paste(.js_texts(node, write), collapse = ","), "]"))
    }
    if (length(node) == 0L) return("{}")
    if (!is.null(part_text) && identical(keys, "$part")) return(part_text(node[[1L]]))
    ordered <- .utf16_order(keys)
    parts <- paste0(.js_strings(keys[ordered]), ":", .js_texts(node[ordered], write))
    return(paste0("{", paste(parts, collapse = ","), "}"))
  }
  .js_scalars(node)
}

# The content hash of a stored value: the SHA-256 of its canonical text.
.value_hash <- function(skeleton, part_text = NULL) {
  .native_sha256(charToRaw(enc2utf8(.canonical_text(skeleton, part_text))))
}
