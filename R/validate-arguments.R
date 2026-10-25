# Value rules for the arguments of cellucid_prepare(). Each function answers
# one question about one argument and returns the value it accepted, so a
# caller's mistake is named before anything is written. The reporter for
# arguments that were never supplied at all lives here too, because it
# quotes the same rules these functions enforce.

.validate_flag <- function(value, arg) {
  if (
    !is.logical(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value)
  ) {
    .abort(arg, " must be exactly TRUE or FALSE.")
  }
  value
}

.validate_string <- function(value, arg) {
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      is.object(value)
  ) {
    .abort(arg, " must be exactly one string.")
  }
  value
}

.validate_required_string <- function(value, arg) {
  value <- .validate_string(value, arg)
  if (!nzchar(value) || !identical(value, trimws(value))) {
    .abort(
      arg,
      " must be one non-empty string without leading or trailing whitespace."
    )
  }
  value
}

# What a valid value of each required argument is, phrased so it can be read as
# the completion of "<argument>: ...". Each line repeats the rule its own
# validator enforces, because a caller who never supplied the argument has no
# other way to learn it.
.required_argument_guidance <- c(
  out_path = paste0(
    "the path to write: a folder ending in \".zarr\" for format = \"zarr\", ",
    "or a file ending in \".zarr.zip\" for format = \"zarr.zip\""
  ),
  format = paste0(
    "the kind of output to write, \"zarr\" (a .zarr folder of small files: ",
    "for GitHub, no file over 100 MB, no LFS) or \"zarr.zip\" (the same store ",
    "in one .zarr.zip file: for single-file hosts such as Zenodo, Figshare, ",
    "Hugging Face, S3 or a lab server, local files, cellucid serve and Jupyter)"
  ),
  embeddings = paste0(
    "a named list mapping each obsm key to a matrix with one row per cell, ",
    "for example list(X_umap = umap)"
  ),
  dataset_id = paste0(
    "one portable identifier of 1-180 ASCII letters, numbers, '.', '_', or ",
    "'-', beginning with a letter or number and not ending with '.'"
  ),
  dataset_name = "one non-empty name for the dataset, shown to the reader verbatim"
)

# R's own message for an argument with no default -- 'argument "dataset_id" is
# missing, with no default' -- is translated into the session language, so it
# does not read the same twice, and it names no valid value. This one is fixed
# text, names every argument that was left out at once, and says what each one
# accepts.
# The output formats, the suffix out_path ends in for each, and what each one
# is for. There is no default: a folder of many files and one archive suit
# different places, and only the caller knows where the store goes. (An .h5ad
# file is written by Python's cellucid.prepare(format = "h5ad"); this package
# writes no HDF5, so it depends on no HDF5 library.)
.output_formats <- c(zarr = ".zarr", zarr.zip = ".zarr.zip")

.format_choices <- paste0(
  "format must be \"zarr\" or \"zarr.zip\":\n",
  "  \"zarr\"      a .zarr folder of small files: for GitHub (no file over ",
  "100 MB, no LFS).\n",
  "  \"zarr.zip\"  the same store in one .zarr.zip file: for single-file ",
  "hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files, ",
  "cellucid serve and Jupyter.\n",
  "(.h5ad is written by Python's cellucid.prepare(format = \"h5ad\"); ",
  "this package writes no HDF5.)"
)

# format is matched as match.arg() matches: exactly, or by an unambiguous
# prefix ("zarr." is "zarr.zip"). Anything else is refused with the choices.
.validate_format <- function(format) {
  if (!is.character(format) || length(format) != 1L || is.na(format) || is.object(format)) {
    .abort(.format_choices, "\nGot ", .format_value(format), ".")
  }
  position <- pmatch(format, names(.output_formats), nomatch = 0L)
  if (position == 0L) {
    .abort(.format_choices, "\nGot \"", format, "\".")
  }
  names(.output_formats)[[position]]
}

# The Zarr format written: 2 when NULL (the default, read by every anndata
# release since 0.8 and by zarr-python 2 and 3), or 3 (anndata 0.11 and later,
# zarr-python 3). Python's cellucid.prepare(zarr_format=) takes the same.
.validate_zarr_format <- function(zarr_format) {
  if (is.null(zarr_format)) {
    return(2L)
  }
  if (is.numeric(zarr_format) && length(zarr_format) == 1L && !is.object(zarr_format) &&
    !is.na(zarr_format) && zarr_format %in% c(2, 3)) {
    return(as.integer(zarr_format))
  }
  .abort("zarr_format must be 2 or 3, got ", .format_value(zarr_format), ".")
}

.format_value <- function(value) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    return(paste0("\"", value, "\""))
  }
  paste(utils::capture.output(dput(value)), collapse = " ")
}

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
  .abort(
    "cellucid_prepare() is missing ",
    counted,
    ":\n",
    listed
  )
}

# One plain positive whole number; a classed number is refused rather than
# unclassed, so the caller learns that its attribute would have been dropped.
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
    .abort(arg, " must be one positive integer.")
  }
  as.integer(value)
}

# A set of names: character storage, no class, no dimensions, nothing
# missing or empty.
.validate_character_vector <- function(values, what) {
  if (
    !is.character(values) ||
      is.object(values) ||
      !is.null(dim(values))
  ) {
    .abort(what, " must be a native character vector.")
  }
  if (anyNA(values) || any(!nzchar(values))) {
    .abort(
      what,
      " must contain only non-missing, non-empty strings."
    )
  }
  invisible(values)
}
