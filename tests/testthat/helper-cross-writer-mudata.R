# The MuData cross-writer cases: a multimodal object written by
# cellucid_prepare() as a MuData Zarr store, beside a JSON description of the
# very same input in the form cellucid-python's tests/cross_writer_mudata.py
# reads (python/cross_writer_mudata.py here, the same file byte for byte).
# fixtures/cross_writer_mudata/ holds the cases cellucid-python checks its own
# writer against (generate.R writes them; test-cross-writer.R checks they are
# what this package writes now). The folder is in the source repository only
# (.Rbuildignore): Zarr's dot-files and long paths are not package files.

# The JSON form of one matrix (cells x features for X) or data frame.
mudata_case_matrix <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    x <- methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix")
    return(list(kind = "sparse", shape = I(dim(x)), p = I(x@p), i = I(x@i), x = I(x@x)))
  }
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  list(kind = "dense", type = typeof(x), shape = I(dim(x)), values = I(as.vector(t(x))))
}

mudata_case_frame <- function(frame, index) {
  column <- function(name, values) {
    if (is.factor(values)) {
      return(list(name = name, type = "factor", values = I(as.character(values)),
                  levels = I(levels(values)), ordered = is.ordered(values)))
    }
    list(name = name, type = typeof(values), values = I(values))
  }
  list(index = I(index), columns = I(unname(Map(column, names(frame), frame))))
}

# The input cellucid_prepare(data = object) reads, as the harness describes
# it: the object's own components, through the one extraction the writer uses.
mudata_case_spec <- function(object, dataset_id, dataset_name) {
  components <- cellucid:::.dataset_components(object, list())
  cells <- rownames(components$obs)
  no_columns <- function(n) data.frame(row = seq_len(n))[, 0L, drop = FALSE]
  modalities <- lapply(components$modalities, function(item) {
    features <- item$gene_expression
    if (inherits(features, "cellucid_feature_major")) features <- features$matrix
    own <- if (is.null(item$cells)) cells else item$cells
    var <- if (is.null(item$var)) no_columns(nrow(features)) else item$var
    list(
      name = item$name,
      X = mudata_case_matrix(Matrix::t(features)),
      obs = mudata_case_frame(no_columns(length(own)), own),
      var = mudata_case_frame(var, rownames(features))
    )
  })
  spec <- list(
    arguments = list(dataset_id = dataset_id, dataset_name = dataset_name),
    obs = mudata_case_frame(components$obs, cells),
    obsm = lapply(components$embeddings, mudata_case_matrix),
    modalities = modalities
  )
  if (!is.null(components$connectivities)) spec$obsp <- list(connectivities = mudata_case_matrix(components$connectivities))
  spec
}

# Write one case: <folder>/<name>.zarr (cellucid_prepare) and <name>.json.
write_mudata_case <- function(object, name, folder, zarr_format, dataset_id = name, dataset_name = name) {
  store <- file.path(folder, paste0(name, ".zarr"))
  suppressMessages(cellucid_prepare(
    store, "zarr", data = object, dataset_id = dataset_id, dataset_name = dataset_name,
    zarr_format = zarr_format, force = TRUE
  ))
  jsonlite::write_json(
    mudata_case_spec(object, dataset_id, dataset_name), file.path(folder, paste0(name, ".json")),
    auto_unbox = TRUE, null = "null", na = "null", digits = I(17), pretty = FALSE
  )
  invisible(store)
}

# The cases of fixtures/cross_writer_mudata/: an SCE with an alternative
# experiment, in Zarr format 2 and 3.
mudata_cases <- function() {
  list(
    `sce-altexp-v2` = list(object = function() small_sce(), zarr_format = 2L),
    `sce-altexp-v3` = list(object = function() small_sce(), zarr_format = 3L)
  )
}
