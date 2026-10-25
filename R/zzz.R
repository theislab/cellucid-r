.onUnload <- function(libpath) {
  .stop_viewers_on_unload()
  library.dynam.unload("cellucid", libpath)
}

# Every viewer's server is stopped before the package unloads: its handlers
# call this package's native routines, and a server left bound would keep its
# port (8765, the one Google sign-in accepts) from the next viewer.
.stop_viewers_on_unload <- function() {
  tryCatch(cellucid_stop_all(), error = function(error) {
    .warn("A Cellucid viewer did not stop as the package unloaded: ", conditionMessage(error))
  })
}
