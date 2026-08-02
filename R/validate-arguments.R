# Value rules for the arguments of cellucid_prepare(). Each function answers
# one question about one argument and returns the value it accepted, so a
# caller's mistake is named before any file is opened. The reporter for
# arguments that were never supplied at all lives here too, because it
# quotes the same rules these functions enforce.

.validate_force <- function(force) {
  if (
    !is.logical(force) ||
      length(force) != 1L ||
      is.na(force) ||
      is.object(force)
  ) {
    stop("force must be exactly TRUE or FALSE.", call. = FALSE)
  }
  force
}

.validate_string <- function(value, arg) {
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value)
  ) {
    stop(arg, " must be exactly one string.", call. = FALSE)
  }
  value
}

.validate_required_string <- function(value, arg) {
  value <- .validate_string(value, arg)
  if (!nzchar(value) || !identical(value, trimws(value))) {
    stop(
      arg,
      " must be one non-empty string without leading or trailing whitespace.",
      call. = FALSE
    )
  }
  value
}

# What a valid value of each required argument is, phrased so it can be read as
# the completion of "<argument>: ...". Each line repeats the rule its own
# validator enforces, because a caller who never supplied the argument has no
# other way to learn it.
.required_argument_guidance <- c(
  obs_categorical_dtype = paste0(
    "exactly one of \"uint8\" or \"uint16\", the storage width of every ",
    "categorical code"
  ),
  dataset_name = "one non-empty name for the dataset, shown to the reader verbatim",
  dataset_id = paste0(
    "one portable identifier of 1-180 ASCII letters, numbers, '.', '_', or ",
    "'-', beginning with a letter or number and not ending with '.'"
  )
)

# R's own message for an argument with no default -- 'argument "dataset_id" is
# missing, with no default' -- is translated into the session language, so it
# does not read the same twice, and it names no valid value. This one is fixed
# text, names every argument that was left out at once, and says what each one
# accepts. The cellucid Python package reports the same set in one TypeError.
.require_supplied_arguments <- function(absent) {
  if (length(absent) == 0L) {
    return(invisible(NULL))
  }
  counted <- if (length(absent) == 1L) {
    "1 required argument"
  } else {
    paste0(length(absent), " required arguments")
  }
  listed <- paste0(
    "  - ",
    absent,
    ": ",
    unname(.required_argument_guidance[absent]),
    ".",
    collapse = "\n"
  )
  stop(
    "cellucid_prepare() is missing ",
    counted,
    ":\n",
    listed,
    call. = FALSE
  )
}

.validate_obs_categorical_dtype <- function(value) {
  supported <- c("uint8", "uint16")
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value) ||
      !value %in% supported
  ) {
    stop(
      "obs_categorical_dtype must be exactly one of \"uint8\" or \"uint16\".",
      call. = FALSE
    )
  }
  value
}

# Every numeric argument below carries the same class term the string and
# logical ones do, and it is load-bearing rather than symmetric for its own
# sake. These four validators return the caller's own value, and each of those
# values is published in a manifest: compression and quantization in the obs
# and var schemas, centroid_outlier_quantile in the obs manifest. jsonlite has
# no asJSON method for an arbitrary S3 class, so a classed number that gets
# past here does not fail as a bad argument -- it fails much later, inside the
# manifest writer, with `No method asJSON S3 class: <class>`, which names no
# argument and tells the caller nothing they can act on.
#
# .validate_positive_integer() is the exception that shows the rule: it
# returns as.integer(value), so the class never survives it. The term is there
# too, because a caller who passed a classed count should be told so rather
# than have the attribute silently dropped -- the same reason a classed column
# is refused instead of unclassed in .is_numeric_data_frame().
.normalize_compression <- function(compression) {
  if (is.null(compression)) {
    return(NULL)
  }
  if (
    !is.numeric(compression) ||
      length(compression) != 1L ||
      is.na(compression) ||
      is.object(compression) ||
      !is.finite(compression) ||
      compression != floor(compression) ||
      compression < 1 ||
      compression > 9
  ) {
    stop(
      "compression must be NULL or one integer from 1 to 9.",
      call. = FALSE
    )
  }
  compression
}

.validate_quantization_bits <- function(bits, arg) {
  if (
    !is.numeric(bits) ||
      length(bits) != 1L ||
      is.na(bits) ||
      is.object(bits) ||
      !is.finite(bits) ||
      !(bits %in% c(8, 16))
  ) {
    stop(arg, " must be exactly 8 or 16.", call. = FALSE)
  }
  bits
}

.validate_centroid_outlier_quantile <- function(value) {
  if (is.null(value)) {
    return(NULL)
  }
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value) ||
      !is.finite(value) ||
      value <= 0.5 ||
      value >= 1
  ) {
    stop(
      "centroid_outlier_quantile must be NULL or one finite number ",
      "strictly between 0.5 and 1.",
      call. = FALSE
    )
  }
  value
}

.validate_positive_integer <- function(value, arg) {
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value) ||
      !is.finite(value) ||
      value != floor(value) ||
      value < 1 ||
      value > .Machine$integer.max
  ) {
    stop(arg, " must be one positive integer.", call. = FALSE)
  }
  as.integer(value)
}

# The plural of .validate_string(): a set of identifiers rather than one
# string, and held to the same three things -- character storage, no class,
# no dimensions. The class term is not decoration. Identity here is decided
# with base semantics throughout -- duplicated(), setdiff(), %in%, setNames(),
# and the identical() that checks the written var manifest against the staged
# gene names -- and identical() compares attributes, so a classed character
# vector that survives this check reaches that assertion and fails it. The
# caller then gets an internal manifest message naming a staging path instead
# of being told which argument was wrong. Every data axis in this package
# already refuses a classed value; the identifier axes refuse one here.
.validate_character_vector <- function(values, what) {
  if (
    !is.character(values) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    stop(what, " must be a native character vector.", call. = FALSE)
  }
  if (anyNA(values) || any(!nzchar(values))) {
    stop(
      what,
      " must contain only non-missing, non-empty strings.",
      call. = FALSE
    )
  }
  invisible(values)
}
