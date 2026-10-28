# A Cellucid-ready MuData written by R is a standard MuData: Python's mudata
# reads it, modality by modality, with the values R put in. The test runs
# wherever a Python with mudata is found: the interpreter named by
# CELLUCID_MUDATA_PYTHON, else CELLUCID_PYTHON, else python3 on the PATH
# (require_python(), helper-objects.R). It is skipped otherwise, and on CRAN;
# where CELLUCID_MUDATA_PYTHON names an interpreter without mudata it fails.

mudata_python <- function() {
  require_python("import mudata", c("CELLUCID_MUDATA_PYTHON", "CELLUCID_PYTHON"), "mudata")
}

test_that("mudata reads the MuData store of a multimodal object", {
  skip_without_sce()
  python <- mudata_python()
  object <- small_sce()
  store <- suppressMessages(cellucid_prepare(file.path(tempfile(), "m.zarr"), "zarr", data = object, dataset_id = "m", dataset_name = "M", priority_seed = NULL))
  output <- suppressWarnings(system2(python, shQuote(c(testthat::test_path("python", "read_with_mudata.py"), store)), stdout = TRUE, stderr = FALSE))
  read <- jsonlite::fromJSON(output[length(output)], simplifyVector = FALSE)
  expect_identical(read$n_obs, ncol(object))
  expect_identical(unlist(read$obs_names), colnames(object))
  expect_identical(unlist(read$mod), c("main", "ADT"))
  # mudata adds one boolean mask per modality to obsm when it reads.
  expect_true(all(c("X_PCA", "X_UMAP") %in% unlist(read$obsm)))
  expect_identical(read$cellucid_id, "m")
  main <- read$modalities$main
  expect_identical(unlist(main$var_names), rownames(object))
  expect_length(main$layers, 0L)
  logcounts <- as.matrix(SummarizedExperiment::assay(object, "logcounts"))
  x <- do.call(rbind, lapply(main$X, unlist))
  expect_equal(x, as_float32_matrix(t(logcounts)), ignore_attr = TRUE)
  adt <- read$modalities$ADT
  expect_identical(unlist(adt$var_names), paste0("ab", 1:4))
  expect_identical(unlist(read$obsmap$ADT), seq_len(ncol(object)))
})
