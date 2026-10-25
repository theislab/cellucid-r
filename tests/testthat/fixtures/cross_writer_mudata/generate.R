# Write the MuData cross-writer cases of this folder: for each case of
# mudata_cases() (tests/testthat/helper-cross-writer-mudata.R), <case>.zarr,
# the store cellucid_prepare() writes, and <case>.json, the same input in the
# form cellucid-python's tests/cross_writer_mudata.py reads. cellucid-python
# checks its own writer against them; test-cross-writer.R checks they are
# still what this package writes. Run from the repository root, with
# cellucid, SingleCellExperiment and jsonlite installed:
#
#   Rscript tests/testthat/fixtures/cross_writer_mudata/generate.R
suppressMessages(library(cellucid))
for (helper in c("helper-objects.R", "helper-cross-writer-mudata.R")) {
  source(file.path("tests", "testthat", helper))
}
folder <- file.path("tests", "testthat", "fixtures", "cross_writer_mudata")
for (name in names(mudata_cases())) {
  case <- mudata_cases()[[name]]
  write_mudata_case(case$object(), name, folder, case$zarr_format)
  cat("wrote", name, "\n")
}
