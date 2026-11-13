# The fast run: each part of the package once, in about a minute on an idle
# machine. Run it from the package's root after every change (Rscript tools/test-core.R);
# testthat::test_local() is the full run.
core <- c(
  "dataset-from",              # prepare and the writer, from the objects people hold
  "viewer-protocol",           # the viewer channel's messages
  "viewer-widget",             # the served store and the channel, end to end without a port
  "sessions",                  # sessions and the apply
  "data-mapping-conformance",  # the data-mapping rules against the shared fixture
  "subsample-conformance",     # the subsampling rules against the shared fixture
  "actions",                   # the generated actions
  "comment-history"            # comments and test names hold no history
)
files <- sub("^test-(.*)\\.R$", "\\1", list.files("tests/testthat", "^test-.*\\.R$"))
stopifnot(all(core %in% files))  # a renamed file fails here; it does not drop out
testthat::test_local(filter = paste0("^(", paste(core, collapse = "|"), ")$"))
