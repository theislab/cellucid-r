# The complete surface of src/export_lock.c. Every .Call() in this package
# is one of these, so the R side of the native contract is readable in one
# screen and the callers never spell a symbol name themselves.

.native_export_lock_acquire <- function(lock_path) {
  .Call(C_cellucid_export_lock_acquire, lock_path)
}

.native_export_lock_release <- function(lock_handle) {
  .Call(C_cellucid_export_lock_release, lock_handle)
}

.native_process_handle_count <- function() {
  .Call(C_cellucid_process_handle_count)
}

.native_export_path_info <- function(path) {
  .Call(C_cellucid_export_path_info, path)
}

.native_export_transaction_id <- function() {
  .Call(C_cellucid_export_transaction_id)
}

.native_export_write_journal <- function(path, contents) {
  .Call(C_cellucid_export_write_journal, path, contents)
}

.native_export_sync_directory <- function(path) {
  .Call(C_cellucid_export_sync_directory, path)
}
