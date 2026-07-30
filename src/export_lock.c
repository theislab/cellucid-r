#include <R.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>

#include <limits.h>
#include <string.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

typedef struct cellucid_export_lock {
#ifdef _WIN32
  HANDLE handle;
#else
  int descriptor;
  pid_t owner_pid;
#endif
  int locked;
  int is_registered;
  SEXP finalizer_weakref;
  struct cellucid_export_lock *previous;
  struct cellucid_export_lock *next;
} cellucid_export_lock;

static SEXP cellucid_export_lock_tag = NULL;
static cellucid_export_lock *cellucid_export_lock_head = NULL;
#ifndef _WIN32
static pid_t cellucid_export_lock_registry_pid = (pid_t)0;
#endif

static void cellucid_register_export_lock(cellucid_export_lock *lock) {
  lock->previous = NULL;
  lock->next = cellucid_export_lock_head;
  if (cellucid_export_lock_head != NULL) {
    cellucid_export_lock_head->previous = lock;
  }
  cellucid_export_lock_head = lock;
  lock->is_registered = 1;
}

static void cellucid_unregister_export_lock(cellucid_export_lock *lock) {
  if (!lock->is_registered) {
    return;
  }
  if (lock->previous == NULL) {
    cellucid_export_lock_head = lock->next;
  } else {
    lock->previous->next = lock->next;
  }
  if (lock->next != NULL) {
    lock->next->previous = lock->previous;
  }
  lock->previous = NULL;
  lock->next = NULL;
  lock->is_registered = 0;
}

#ifdef _WIN32

static WCHAR *cellucid_utf8_to_extended_path(const char *path) {
  int converted_length;
  int prefix_length;
  int source_offset;
  int index;
  WCHAR *converted;
  WCHAR *extended;
  const WCHAR drive_prefix[] = L"\\\\?\\";
  const WCHAR unc_prefix[] = L"\\\\?\\UNC\\";

  converted_length =
      MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, NULL, 0);
  if (converted_length <= 0) {
    Rf_error("Could not convert the export lock path to UTF-16 (Windows error %lu).",
             (unsigned long)GetLastError());
  }
  converted = (WCHAR *)R_alloc((size_t)converted_length, sizeof(WCHAR));
  if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, converted,
                          converted_length) != converted_length) {
    Rf_error("Could not materialize the export lock path as UTF-16.");
  }
  for (index = 0; index < converted_length - 1; ++index) {
    if (converted[index] == L'/') {
      converted[index] = L'\\';
    }
  }

  if (converted_length >= 5 && converted[0] == L'\\' &&
      converted[1] == L'\\' && converted[2] == L'.' &&
      converted[3] == L'\\') {
    Rf_error("Windows device paths cannot own an export lock.");
  }
  if (converted_length >= 5 && converted[0] == L'\\' &&
      converted[1] == L'\\' && converted[2] == L'?' &&
      converted[3] == L'\\') {
    int is_prefixed_drive =
        converted_length >= 8 &&
        ((converted[4] >= L'A' && converted[4] <= L'Z') ||
         (converted[4] >= L'a' && converted[4] <= L'z')) &&
        converted[5] == L':' && converted[6] == L'\\';
    int is_prefixed_unc =
        converted_length >= 9 &&
        (converted[4] == L'U' || converted[4] == L'u') &&
        (converted[5] == L'N' || converted[5] == L'n') &&
        (converted[6] == L'C' || converted[6] == L'c') &&
        converted[7] == L'\\';
    if (!is_prefixed_drive && !is_prefixed_unc) {
      Rf_error("Unsupported extended Windows export lock namespace.");
    }
    if (converted_length > 32767) {
      Rf_error("Export lock path exceeds the Windows UTF-16 path limit.");
    }
    return converted;
  }

  if (converted_length >= 3 && converted[0] == L'\\' &&
      converted[1] == L'\\') {
    prefix_length = 8;
    source_offset = 2;
    if (converted_length > 32767 - prefix_length + source_offset) {
      Rf_error("Export lock path exceeds the Windows UTF-16 path limit.");
    }
    extended = (WCHAR *)R_alloc(
        (size_t)(prefix_length + converted_length - source_offset),
        sizeof(WCHAR));
    memcpy(extended, unc_prefix, (size_t)prefix_length * sizeof(WCHAR));
  } else if (
      converted_length >= 4 &&
      ((converted[0] >= L'A' && converted[0] <= L'Z') ||
       (converted[0] >= L'a' && converted[0] <= L'z')) &&
      converted[1] == L':' && converted[2] == L'\\') {
    prefix_length = 4;
    source_offset = 0;
    if (converted_length > 32767 - prefix_length) {
      Rf_error("Export lock path exceeds the Windows UTF-16 path limit.");
    }
    extended = (WCHAR *)R_alloc(
        (size_t)(prefix_length + converted_length), sizeof(WCHAR));
    memcpy(extended, drive_prefix, (size_t)prefix_length * sizeof(WCHAR));
  } else {
    Rf_error("Export lock path must be absolute on Windows.");
  }
  memcpy(
      extended + prefix_length, converted + source_offset,
      (size_t)(converted_length - source_offset) * sizeof(WCHAR));
  return extended;
}

static int cellucid_unlock_native(cellucid_export_lock *lock) {
  OVERLAPPED overlap;
  if (!lock->locked) {
    return 0;
  }
  memset(&overlap, 0, sizeof(overlap));
  if (!UnlockFileEx(lock->handle, 0, 1, 0, &overlap)) {
    return (int)GetLastError();
  }
  lock->locked = 0;
  return 0;
}

static int cellucid_close_native(cellucid_export_lock *lock) {
  if (lock->handle == INVALID_HANDLE_VALUE) {
    return 0;
  }
  if (!CloseHandle(lock->handle)) {
    return (int)GetLastError();
  }
  lock->handle = INVALID_HANDLE_VALUE;
  return 0;
}

#else

static int cellucid_lstat_native(const char *path, struct stat *information) {
  int result;
  do {
    result = lstat(path, information);
  } while (result == -1 && errno == EINTR);
  return result;
}

static int cellucid_fstat_native(int descriptor, struct stat *information) {
  int result;
  do {
    result = fstat(descriptor, information);
  } while (result == -1 && errno == EINTR);
  return result;
}

static int cellucid_open_native(const char *path, int flags) {
  int descriptor;
  do {
    descriptor = open(path, flags, S_IRUSR | S_IWUSR);
  } while (descriptor == -1 && errno == EINTR);
  return descriptor;
}

static int cellucid_record_lock_native(
    int descriptor, struct flock *operation) {
  int result;
  do {
    result = fcntl(descriptor, F_SETLK, operation);
  } while (result == -1 && errno == EINTR);
  return result;
}

static int cellucid_unlock_native(cellucid_export_lock *lock) {
  struct flock operation;
  if (!lock->locked || lock->owner_pid != getpid()) {
    return 0;
  }
  memset(&operation, 0, sizeof(operation));
  operation.l_type = F_UNLCK;
  operation.l_whence = SEEK_SET;
  operation.l_start = 0;
  operation.l_len = 0;
  if (cellucid_record_lock_native(lock->descriptor, &operation) == -1) {
    return errno;
  }
  lock->locked = 0;
  return 0;
}

static int cellucid_close_native(cellucid_export_lock *lock) {
  int result;
  if (lock->descriptor < 0) {
    return 0;
  }
  result = close(lock->descriptor);
  if (result == -1) {
    int close_error = errno;
    /*
     * POSIX leaves the descriptor state unspecified after EINTR. Retrying can
     * close an unrelated reused descriptor, so this handle is retired.
     */
    lock->descriptor = -1;
    return close_error;
  }
  lock->descriptor = -1;
  return 0;
}

#endif

static void cellucid_export_lock_finalizer(SEXP external_pointer) {
  cellucid_export_lock *lock =
      (cellucid_export_lock *)R_ExternalPtrAddr(external_pointer);
  if (lock == NULL) {
    return;
  }
  (void)cellucid_unlock_native(lock);
  (void)cellucid_close_native(lock);
  cellucid_unregister_export_lock(lock);
  R_ClearExternalPtr(external_pointer);
  R_Free(lock);
}

static void cellucid_finalize_export_lock(cellucid_export_lock *lock) {
  SEXP finalizer_weakref = lock->finalizer_weakref;
  R_RunWeakRefFinalizer(finalizer_weakref);
}

#ifndef _WIN32
static void cellucid_refresh_export_lock_process(void) {
  pid_t process_id = getpid();
  cellucid_export_lock *lock;
  int cleanup_failure_count = 0;
  int first_unlock_error = 0;
  int first_close_error = 0;

  if (cellucid_export_lock_registry_pid == process_id) {
    return;
  }

  /*
   * POSIX record locks are not inherited by a forked child, but the open
   * descriptors and R external pointers are. Retire every copied handle
   * before the child performs any native lock operation. Set the new PID
   * first so finalization cannot recurse into this refresh.
   */
  cellucid_export_lock_registry_pid = process_id;
  lock = cellucid_export_lock_head;
  while (lock != NULL) {
    cellucid_export_lock *next = lock->next;
    int unlock_error = cellucid_unlock_native(lock);
    int close_error = cellucid_close_native(lock);
    if (unlock_error != 0 || close_error != 0) {
      if (cleanup_failure_count == 0) {
        first_unlock_error = unlock_error;
        first_close_error = close_error;
      }
      cleanup_failure_count += 1;
    }
    /*
     * close() retires the descriptor even when it reports an error because
     * POSIX leaves descriptor state unspecified after EINTR.
     */
    cellucid_finalize_export_lock(lock);
    lock = next;
  }
  if (cleanup_failure_count != 0) {
    Rf_error(
        "Could not retire %d fork-inherited export lock handle(s) "
        "(first unlock/close errors %d/%d).",
        cleanup_failure_count, first_unlock_error, first_close_error);
  }
}
#else
static void cellucid_refresh_export_lock_process(void) {}
#endif

static SEXP cellucid_export_lock_acquire(SEXP path) {
  SEXP external_pointer;
  SEXP finalizer_weakref;
  cellucid_export_lock *lock;
  const char *native_path;
  SEXP path_string;

  cellucid_refresh_export_lock_process();
  if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 ||
      STRING_ELT(path, 0) == NA_STRING) {
    Rf_error("Export lock path must be one non-missing string.");
  }
  path_string = STRING_ELT(path, 0);
  if (LENGTH(path_string) == 0 ||
      memchr(CHAR(path_string), '\0', (size_t)LENGTH(path_string)) != NULL) {
    Rf_error("Export lock path must be non-empty and contain no NUL byte.");
  }
#ifdef _WIN32
  if (Rf_getCharCE(path_string) == CE_BYTES) {
    Rf_error("A byte-encoded export lock path is unsupported on Windows.");
  }
  native_path = Rf_translateCharUTF8(path_string);
#else
  native_path = Rf_getCharCE(path_string) == CE_BYTES
                    ? CHAR(path_string)
                    : Rf_translateChar(path_string);
#endif
  external_pointer =
      PROTECT(R_MakeExternalPtr(NULL, cellucid_export_lock_tag, R_NilValue));
  finalizer_weakref = PROTECT(R_MakeWeakRefC(
      external_pointer, R_NilValue, cellucid_export_lock_finalizer, TRUE));
  lock = R_Calloc(1, cellucid_export_lock);
#ifdef _WIN32
  lock->handle = INVALID_HANDLE_VALUE;
#else
  lock->descriptor = -1;
#endif
  lock->finalizer_weakref = finalizer_weakref;
  R_SetExternalPtrAddr(external_pointer, lock);
  cellucid_register_export_lock(lock);
  UNPROTECT(1);

#ifdef _WIN32
  {
    WCHAR *wide_path = cellucid_utf8_to_extended_path(native_path);
    BY_HANDLE_FILE_INFORMATION information;
    BY_HANDLE_FILE_INFORMATION locked_information;
    OVERLAPPED overlap;
    DWORD lock_error;
    DWORD validation_error;
    int unlock_error;
    int close_error;

    lock->handle = CreateFileW(
        wide_path, GENERIC_READ | GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_ALWAYS,
        FILE_ATTRIBUTE_NORMAL | FILE_FLAG_OPEN_REPARSE_POINT, NULL);
    if (lock->handle == INVALID_HANDLE_VALUE) {
      lock_error = GetLastError();
      cellucid_finalize_export_lock(lock);
      UNPROTECT(1);
      Rf_error("Could not open the export lock file (Windows error %lu).",
               (unsigned long)lock_error);
    }
    if (!GetFileInformationByHandle(lock->handle, &information)) {
      lock_error = GetLastError();
      close_error = cellucid_close_native(lock);
      if (close_error == 0) {
        cellucid_finalize_export_lock(lock);
      }
      UNPROTECT(1);
      if (close_error != 0) {
        Rf_error(
            "Could not inspect or close the export lock file "
            "(Windows errors %lu/%d).",
            (unsigned long)lock_error, close_error);
      }
      Rf_error("Could not inspect the export lock file (Windows error %lu).",
               (unsigned long)lock_error);
    }
    if ((information.dwFileAttributes &
         (FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_REPARSE_POINT)) != 0 ||
        information.nNumberOfLinks != 1) {
      close_error = cellucid_close_native(lock);
      if (close_error == 0) {
        cellucid_finalize_export_lock(lock);
      }
      UNPROTECT(1);
      if (close_error != 0) {
        Rf_error(
            "Export lock path is not a regular non-reparse file and its "
            "handle could not close (Windows error %d).",
            close_error);
      }
      Rf_error("Export lock path must identify a regular non-reparse file.");
    }

    memset(&overlap, 0, sizeof(overlap));
    if (!LockFileEx(
            lock->handle,
            LOCKFILE_EXCLUSIVE_LOCK | LOCKFILE_FAIL_IMMEDIATELY,
            0, 1, 0, &overlap)) {
      lock_error = GetLastError();
      close_error = cellucid_close_native(lock);
      if (close_error == 0) {
        cellucid_finalize_export_lock(lock);
      }
      UNPROTECT(1);
      if (close_error != 0) {
        Rf_error(
            "A failed export lock attempt also failed to close its handle "
            "(Windows errors %lu/%d).",
            (unsigned long)lock_error, close_error);
      }
      if (lock_error == ERROR_LOCK_VIOLATION) {
        return R_NilValue;
      }
      Rf_error("Could not acquire the export lock (Windows error %lu).",
               (unsigned long)lock_error);
    }
    lock->locked = 1;
    if (!GetFileInformationByHandle(lock->handle, &locked_information)) {
      validation_error = GetLastError();
      unlock_error = cellucid_unlock_native(lock);
      close_error = cellucid_close_native(lock);
      if (close_error == 0) {
        cellucid_finalize_export_lock(lock);
      }
      UNPROTECT(1);
      Rf_error(
          "Could not inspect the acquired export lock "
          "(validation/unlock/close Windows errors %lu/%d/%d).",
          (unsigned long)validation_error, unlock_error, close_error);
    }
    if ((locked_information.dwFileAttributes &
         (FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_REPARSE_POINT)) != 0 ||
        locked_information.nNumberOfLinks != 1 ||
        locked_information.dwVolumeSerialNumber !=
            information.dwVolumeSerialNumber ||
        locked_information.nFileIndexHigh != information.nFileIndexHigh ||
        locked_information.nFileIndexLow != information.nFileIndexLow) {
      unlock_error = cellucid_unlock_native(lock);
      close_error = cellucid_close_native(lock);
      if (close_error == 0) {
        cellucid_finalize_export_lock(lock);
      }
      UNPROTECT(1);
      Rf_error(
          "The acquired export lock changed identity or link safety "
          "(unlock/close Windows errors %d/%d).",
          unlock_error, close_error);
    }
  }
#else
  {
    struct stat before_open;
    struct stat descriptor_info;
    struct stat locked_descriptor_info;
    struct stat locked_path_info;
    struct stat path_info;
    struct flock operation;
    int flags = O_RDWR;
    int lock_error;
    int validation_error;
    int unlock_error;
    int close_error;
    int existed_before_open = 1;
    int descriptor_must_match_before_open = 0;

    lock->owner_pid = getpid();
#ifdef O_CLOEXEC
    flags |= O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
    flags |= O_NOFOLLOW;
#endif

    if (cellucid_lstat_native(native_path, &before_open) == -1) {
      if (errno == ENOENT) {
        existed_before_open = 0;
      } else {
        lock_error = errno;
        cellucid_finalize_export_lock(lock);
        UNPROTECT(1);
        Rf_error("Could not inspect the export lock path: %s.",
                 strerror(lock_error));
      }
    } else {
      if (!S_ISREG(before_open.st_mode) || before_open.st_nlink != 1) {
        cellucid_finalize_export_lock(lock);
        UNPROTECT(1);
        Rf_error(
            "Export lock path must identify one non-linked regular file.");
      }
      descriptor_must_match_before_open = 1;
    }

    if (!existed_before_open) {
      lock->descriptor = cellucid_open_native(
          native_path, flags | O_CREAT | O_EXCL);
    } else {
      lock->descriptor = cellucid_open_native(native_path, flags);
    }
    if (lock->descriptor == -1 && !existed_before_open && errno == EEXIST) {
      if (cellucid_lstat_native(native_path, &before_open) == -1) {
        lock_error = errno;
        cellucid_finalize_export_lock(lock);
        UNPROTECT(1);
        Rf_error("Could not validate the concurrently created export lock: %s.",
                 strerror(lock_error));
      }
      if (!S_ISREG(before_open.st_mode) || before_open.st_nlink != 1) {
        cellucid_finalize_export_lock(lock);
        UNPROTECT(1);
        Rf_error(
            "Export lock path must identify one non-linked regular file.");
      }
      descriptor_must_match_before_open = 1;
      lock->descriptor = cellucid_open_native(native_path, flags);
    }
    if (lock->descriptor == -1) {
      lock_error = errno;
      cellucid_finalize_export_lock(lock);
      UNPROTECT(1);
      Rf_error("Could not open the export lock file: %s.",
               strerror(lock_error));
    }
#ifndef O_CLOEXEC
    {
      int descriptor_flags = fcntl(lock->descriptor, F_GETFD);
      if (descriptor_flags == -1 ||
          fcntl(lock->descriptor, F_SETFD,
                descriptor_flags | FD_CLOEXEC) == -1) {
        lock_error = errno;
        close_error = cellucid_close_native(lock);
        cellucid_finalize_export_lock(lock);
        UNPROTECT(1);
        if (close_error != 0) {
          Rf_error(
              "Could not protect or close the export lock descriptor: %s.",
              strerror(close_error));
        }
        Rf_error("Could not protect the export lock descriptor: %s.",
                 strerror(lock_error));
      }
    }
#endif

    if (cellucid_fstat_native(lock->descriptor, &descriptor_info) == -1 ||
        cellucid_lstat_native(native_path, &path_info) == -1) {
      lock_error = errno;
      close_error = cellucid_close_native(lock);
      cellucid_finalize_export_lock(lock);
      UNPROTECT(1);
      if (close_error != 0) {
        Rf_error("Could not inspect or close the export lock descriptor: %s.",
                 strerror(close_error));
      }
      Rf_error("Could not inspect the export lock descriptor: %s.",
               strerror(lock_error));
    }
    if (!S_ISREG(descriptor_info.st_mode) ||
        !S_ISREG(path_info.st_mode) || descriptor_info.st_nlink != 1 ||
        path_info.st_nlink != 1 ||
        descriptor_info.st_dev != path_info.st_dev ||
        descriptor_info.st_ino != path_info.st_ino ||
        (descriptor_must_match_before_open &&
         (descriptor_info.st_dev != before_open.st_dev ||
          descriptor_info.st_ino != before_open.st_ino))) {
      close_error = cellucid_close_native(lock);
      cellucid_finalize_export_lock(lock);
      UNPROTECT(1);
      if (close_error != 0) {
        Rf_error(
            "The unsafe export lock descriptor also failed to close: %s.",
            strerror(close_error));
      }
      Rf_error(
          "Export lock path must identify one non-linked regular file.");
    }

    memset(&operation, 0, sizeof(operation));
    operation.l_type = F_WRLCK;
    operation.l_whence = SEEK_SET;
    operation.l_start = 0;
    operation.l_len = 0;
    if (cellucid_record_lock_native(lock->descriptor, &operation) == -1) {
      lock_error = errno;
      close_error = cellucid_close_native(lock);
      cellucid_finalize_export_lock(lock);
      UNPROTECT(1);
      if (close_error != 0) {
        Rf_error(
            "A failed export lock attempt also failed to close its "
            "descriptor: %s.",
            strerror(close_error));
      }
      if (lock_error == EACCES || lock_error == EAGAIN) {
        return R_NilValue;
      }
      Rf_error("Could not acquire the export lock: %s.",
               strerror(lock_error));
    }
    lock->locked = 1;
    if (cellucid_fstat_native(
            lock->descriptor, &locked_descriptor_info) == -1 ||
        cellucid_lstat_native(native_path, &locked_path_info) == -1) {
      validation_error = errno;
      unlock_error = cellucid_unlock_native(lock);
      close_error = cellucid_close_native(lock);
      cellucid_finalize_export_lock(lock);
      UNPROTECT(1);
      Rf_error(
          "Could not inspect the acquired export lock "
          "(validation error %d: %s; unlock/close errors %d/%d).",
          validation_error, strerror(validation_error),
          unlock_error, close_error);
    }
    if (!S_ISREG(locked_descriptor_info.st_mode) ||
        !S_ISREG(locked_path_info.st_mode) ||
        locked_descriptor_info.st_nlink != 1 ||
        locked_path_info.st_nlink != 1 ||
        locked_descriptor_info.st_dev != locked_path_info.st_dev ||
        locked_descriptor_info.st_ino != locked_path_info.st_ino ||
        locked_descriptor_info.st_dev != descriptor_info.st_dev ||
        locked_descriptor_info.st_ino != descriptor_info.st_ino ||
        locked_path_info.st_dev != path_info.st_dev ||
        locked_path_info.st_ino != path_info.st_ino ||
        (descriptor_must_match_before_open &&
         (locked_descriptor_info.st_dev != before_open.st_dev ||
          locked_descriptor_info.st_ino != before_open.st_ino))) {
      unlock_error = cellucid_unlock_native(lock);
      close_error = cellucid_close_native(lock);
      cellucid_finalize_export_lock(lock);
      UNPROTECT(1);
      Rf_error(
          "The acquired export lock changed identity or link safety "
          "(unlock/close errors %d/%d).",
          unlock_error, close_error);
    }
  }
#endif

  UNPROTECT(1);
  return external_pointer;
}

static SEXP cellucid_export_lock_release(SEXP external_pointer) {
  SEXP status;
  cellucid_export_lock *lock;
  int unlock_error;
  int close_error;

  if (TYPEOF(external_pointer) != EXTPTRSXP ||
      R_ExternalPtrTag(external_pointer) != cellucid_export_lock_tag) {
    Rf_error("Export lock handle is not owned by Cellucid.");
  }
  cellucid_refresh_export_lock_process();
  status = PROTECT(Rf_allocVector(INTSXP, 3));
  lock = (cellucid_export_lock *)R_ExternalPtrAddr(external_pointer);
  if (lock == NULL) {
    INTEGER(status)[0] = 0;
    INTEGER(status)[1] = 0;
    INTEGER(status)[2] = 1;
    UNPROTECT(1);
    return status;
  }

  unlock_error = cellucid_unlock_native(lock);
  close_error = cellucid_close_native(lock);
  INTEGER(status)[0] = unlock_error;
  INTEGER(status)[1] = close_error;
#ifdef _WIN32
  INTEGER(status)[2] = close_error == 0;
#else
  INTEGER(status)[2] = 1;
#endif
  if (INTEGER(status)[2]) {
    cellucid_finalize_export_lock(lock);
  }
  UNPROTECT(1);
  return status;
}

static SEXP cellucid_export_lock_drain(void) {
  SEXP status;
  cellucid_export_lock *lock;
  int lock_count = 0;
  int row = 0;

  cellucid_refresh_export_lock_process();
  for (lock = cellucid_export_lock_head; lock != NULL; lock = lock->next) {
    if (lock_count >= INT_MAX / 3) {
      Rf_error("Too many live export lock handles to drain safely.");
    }
    lock_count += 1;
  }

  status = PROTECT(Rf_allocMatrix(INTSXP, lock_count, 3));
  for (row = 0; row < lock_count; ++row) {
    INTEGER(status)[row] = 0;
    INTEGER(status)[row + lock_count] = 0;
    INTEGER(status)[row + (2 * lock_count)] = 1;
  }
  row = 0;
  lock = cellucid_export_lock_head;
  while (lock != NULL) {
    cellucid_export_lock *next = lock->next;
    int unlock_error = cellucid_unlock_native(lock);
    int close_error = cellucid_close_native(lock);
    int resource_closed;
#ifdef _WIN32
    resource_closed = close_error == 0;
#else
    resource_closed = 1;
#endif
    INTEGER(status)[row] = unlock_error;
    INTEGER(status)[row + lock_count] = close_error;
    INTEGER(status)[row + (2 * lock_count)] = resource_closed;
    if (resource_closed) {
      /*
       * Running the registered weak-reference finalizer now both frees the
       * native state and clears the finalizer function pointer before the
       * package DLL is unloaded. A later gc() of an escaped external pointer
       * is therefore inert instead of jumping into unloaded code.
       */
      cellucid_finalize_export_lock(lock);
    }
    row += 1;
    lock = next;
  }

  UNPROTECT(1);
  return status;
}

static SEXP cellucid_process_handle_count(void) {
#ifdef _WIN32
  DWORD count;
  if (!GetProcessHandleCount(GetCurrentProcess(), &count)) {
    return ScalarInteger(NA_INTEGER);
  }
  return ScalarReal((double)count);
#else
  return ScalarInteger(NA_INTEGER);
#endif
}

static SEXP cellucid_export_path_info(SEXP path) {
  SEXP result;
  SEXP path_string;
  const char *native_path;
  int kind = 0;
  double link_count = 0;

  if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 ||
      STRING_ELT(path, 0) == NA_STRING) {
    Rf_error("Export transaction path must be one non-missing string.");
  }
  path_string = STRING_ELT(path, 0);
  if (LENGTH(path_string) == 0 ||
      memchr(CHAR(path_string), '\0', (size_t)LENGTH(path_string)) != NULL) {
    Rf_error("Export transaction path must be non-empty and contain no NUL byte.");
  }
#ifdef _WIN32
  if (Rf_getCharCE(path_string) == CE_BYTES) {
    Rf_error("A byte-encoded export transaction path is unsupported on Windows.");
  }
  native_path = Rf_translateCharUTF8(path_string);
  {
    WCHAR *wide_path = cellucid_utf8_to_extended_path(native_path);
    WIN32_FILE_ATTRIBUTE_DATA attributes;
    BY_HANDLE_FILE_INFORMATION information;
    HANDLE handle;
    DWORD path_error;

    if (!GetFileAttributesExW(
            wide_path, GetFileExInfoStandard, &attributes)) {
      path_error = GetLastError();
      if (path_error != ERROR_FILE_NOT_FOUND &&
          path_error != ERROR_PATH_NOT_FOUND) {
        Rf_error(
            "Could not inspect the export transaction path "
            "(Windows error %lu).",
            (unsigned long)path_error);
      }
    } else if (
        (attributes.dwFileAttributes & FILE_ATTRIBUTE_REPARSE_POINT) != 0) {
      kind = 3;
    } else {
      handle = CreateFileW(
          wide_path, FILE_READ_ATTRIBUTES,
          FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
          NULL, OPEN_EXISTING,
          FILE_FLAG_BACKUP_SEMANTICS | FILE_FLAG_OPEN_REPARSE_POINT,
          NULL);
      if (handle == INVALID_HANDLE_VALUE) {
        Rf_error(
            "Could not open the export transaction path for inspection "
            "(Windows error %lu).",
            (unsigned long)GetLastError());
      }
      if (!GetFileInformationByHandle(handle, &information)) {
        path_error = GetLastError();
        (void)CloseHandle(handle);
        Rf_error(
            "Could not inspect the export transaction path handle "
            "(Windows error %lu).",
            (unsigned long)path_error);
      }
      if (!CloseHandle(handle)) {
        Rf_error(
            "Could not close the export transaction path handle "
            "(Windows error %lu).",
            (unsigned long)GetLastError());
      }
      link_count = (double)information.nNumberOfLinks;
      if ((information.dwFileAttributes & FILE_ATTRIBUTE_REPARSE_POINT) != 0) {
        kind = 3;
      } else if (
          (information.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0) {
        kind = 2;
      } else if (
          (information.dwFileAttributes & FILE_ATTRIBUTE_DEVICE) == 0) {
        kind = 1;
      } else {
        kind = 4;
      }
    }
  }
#else
  native_path = Rf_getCharCE(path_string) == CE_BYTES
                    ? CHAR(path_string)
                    : Rf_translateChar(path_string);
  {
    struct stat information;
    if (cellucid_lstat_native(native_path, &information) == -1) {
      if (errno != ENOENT) {
        Rf_error("Could not inspect the export transaction path: %s.",
                 strerror(errno));
      }
    } else {
      link_count = (double)information.st_nlink;
      if (S_ISREG(information.st_mode)) {
        kind = 1;
      } else if (S_ISDIR(information.st_mode)) {
        kind = 2;
      } else if (S_ISLNK(information.st_mode)) {
        kind = 3;
      } else {
        kind = 4;
      }
    }
  }
#endif

  result = PROTECT(Rf_allocVector(REALSXP, 2));
  REAL(result)[0] = (double)kind;
  REAL(result)[1] = link_count;
  UNPROTECT(1);
  return result;
}

static SEXP cellucid_export_transaction_id(void) {
  unsigned char bytes[16];
  char encoded[33];
  static const char hexadecimal[] = "0123456789abcdef";
  int index;

#ifdef _WIN32
  {
    typedef BOOLEAN(APIENTRY * cellucid_rtl_gen_random)(PVOID, ULONG);
    HMODULE provider = LoadLibraryW(L"advapi32.dll");
    cellucid_rtl_gen_random generate;
    if (provider == NULL) {
      Rf_error(
          "Could not load the Windows system random provider "
          "(Windows error %lu).",
          (unsigned long)GetLastError());
    }
    generate = (cellucid_rtl_gen_random)GetProcAddress(
        provider, "SystemFunction036");
    if (generate == NULL || !generate(bytes, (ULONG)sizeof(bytes))) {
      DWORD random_error = GetLastError();
      FreeLibrary(provider);
      Rf_error(
          "Could not generate an export transaction identity "
          "(Windows error %lu).",
          (unsigned long)random_error);
    }
    if (!FreeLibrary(provider)) {
      Rf_error(
          "Could not release the Windows system random provider "
          "(Windows error %lu).",
          (unsigned long)GetLastError());
    }
  }
#else
  {
    int descriptor;
    size_t offset = 0;
    int flags = O_RDONLY;
#ifdef O_CLOEXEC
    flags |= O_CLOEXEC;
#endif
    descriptor = cellucid_open_native("/dev/urandom", flags);
    if (descriptor == -1) {
      Rf_error("Could not open the system random provider: %s.",
               strerror(errno));
    }
    while (offset < sizeof(bytes)) {
      ssize_t count = read(
          descriptor, bytes + offset, sizeof(bytes) - offset);
      if (count == -1 && errno == EINTR) {
        continue;
      }
      if (count <= 0) {
        int random_error = count == -1 ? errno : EIO;
        (void)close(descriptor);
        Rf_error("Could not read the system random provider: %s.",
                 strerror(random_error));
      }
      offset += (size_t)count;
    }
    if (close(descriptor) == -1) {
      Rf_error("Could not close the system random provider: %s.",
               strerror(errno));
    }
  }
#endif

  for (index = 0; index < 16; ++index) {
    encoded[index * 2] = hexadecimal[bytes[index] >> 4];
    encoded[index * 2 + 1] = hexadecimal[bytes[index] & 0x0f];
  }
  encoded[32] = '\0';
  return Rf_mkString(encoded);
}

static SEXP cellucid_export_write_journal(SEXP path, SEXP contents) {
  SEXP path_string;
  const char *native_path;
  R_xlen_t content_length;

  if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 ||
      STRING_ELT(path, 0) == NA_STRING) {
    Rf_error("Export transaction journal path must be one non-missing string.");
  }
  if (TYPEOF(contents) != RAWSXP) {
    Rf_error("Export transaction journal contents must be raw bytes.");
  }
  content_length = XLENGTH(contents);
  if (content_length <= 0 || content_length > 512) {
    Rf_error("Export transaction journal contents have an invalid size.");
  }
  path_string = STRING_ELT(path, 0);
  if (LENGTH(path_string) == 0 ||
      memchr(CHAR(path_string), '\0', (size_t)LENGTH(path_string)) != NULL) {
    Rf_error(
        "Export transaction journal path must be non-empty and contain no NUL byte.");
  }

#ifdef _WIN32
  if (Rf_getCharCE(path_string) == CE_BYTES) {
    Rf_error("A byte-encoded export transaction path is unsupported on Windows.");
  }
  native_path = Rf_translateCharUTF8(path_string);
  {
    WCHAR *wide_path = cellucid_utf8_to_extended_path(native_path);
    HANDLE handle = CreateFileW(
        wide_path, GENERIC_WRITE, FILE_SHARE_READ,
        NULL, CREATE_NEW, FILE_ATTRIBUTE_NORMAL, NULL);
    DWORD write_error = ERROR_SUCCESS;
    DWORD close_error = ERROR_SUCCESS;
    DWORD written = 0;
    if (handle == INVALID_HANDLE_VALUE) {
      Rf_error(
          "Could not exclusively create the export transaction journal "
          "(Windows error %lu).",
          (unsigned long)GetLastError());
    }
    if (!WriteFile(
            handle, RAW(contents), (DWORD)content_length, &written, NULL) ||
        written != (DWORD)content_length) {
      write_error = GetLastError();
      if (write_error == ERROR_SUCCESS) {
        write_error = ERROR_WRITE_FAULT;
      }
    } else if (!FlushFileBuffers(handle)) {
      write_error = GetLastError();
    }
    if (!CloseHandle(handle)) {
      close_error = GetLastError();
    }
    if (write_error != ERROR_SUCCESS || close_error != ERROR_SUCCESS) {
      Rf_error(
          "Could not durably write the export transaction journal "
          "(Windows write/close errors %lu/%lu).",
          (unsigned long)write_error, (unsigned long)close_error);
    }
  }
#else
  native_path = Rf_getCharCE(path_string) == CE_BYTES
                    ? CHAR(path_string)
                    : Rf_translateChar(path_string);
  {
    int descriptor;
    int flags = O_WRONLY | O_CREAT | O_EXCL;
    int write_error = 0;
    int close_error = 0;
    size_t offset = 0;
#ifdef O_CLOEXEC
    flags |= O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
    flags |= O_NOFOLLOW;
#endif
    descriptor = cellucid_open_native(native_path, flags);
    if (descriptor == -1) {
      Rf_error(
          "Could not exclusively create the export transaction journal: %s.",
          strerror(errno));
    }
    while (offset < (size_t)content_length) {
      ssize_t count = write(
          descriptor, RAW(contents) + offset,
          (size_t)content_length - offset);
      if (count == -1 && errno == EINTR) {
        continue;
      }
      if (count <= 0) {
        write_error = count == -1 ? errno : EIO;
        break;
      }
      offset += (size_t)count;
    }
    if (write_error == 0) {
      int sync_result;
      do {
        sync_result = fsync(descriptor);
      } while (sync_result == -1 && errno == EINTR);
      if (sync_result == -1) {
        write_error = errno;
      }
    }
    if (close(descriptor) == -1) {
      close_error = errno;
    }
    if (write_error != 0 || close_error != 0) {
      Rf_error(
          "Could not durably write the export transaction journal "
          "(write error %d: %s; close error %d: %s).",
          write_error,
          write_error == 0 ? "none" : strerror(write_error),
          close_error,
          close_error == 0 ? "none" : strerror(close_error));
    }
  }
#endif

  return ScalarLogical(1);
}

static SEXP cellucid_export_sync_directory(SEXP path) {
#ifdef _WIN32
  (void)path;
  return ScalarLogical(0);
#else
  SEXP path_string;
  const char *native_path;
  int descriptor;
  int flags = O_RDONLY;
  int sync_result;
  int sync_error;
  int close_error;
  struct stat information;

  if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 ||
      STRING_ELT(path, 0) == NA_STRING) {
    Rf_error("Export transaction directory path must be one non-missing string.");
  }
  path_string = STRING_ELT(path, 0);
  if (LENGTH(path_string) == 0 ||
      memchr(CHAR(path_string), '\0', (size_t)LENGTH(path_string)) != NULL) {
    Rf_error(
        "Export transaction directory path must be non-empty and contain no NUL byte.");
  }
  native_path = Rf_getCharCE(path_string) == CE_BYTES
                    ? CHAR(path_string)
                    : Rf_translateChar(path_string);
#ifdef O_CLOEXEC
  flags |= O_CLOEXEC;
#endif
#ifdef O_DIRECTORY
  flags |= O_DIRECTORY;
#endif
#ifdef O_NOFOLLOW
  flags |= O_NOFOLLOW;
#endif
  descriptor = cellucid_open_native(native_path, flags);
  if (descriptor == -1) {
    Rf_error("Could not open the export transaction directory: %s.",
             strerror(errno));
  }
  if (cellucid_fstat_native(descriptor, &information) == -1 ||
      !S_ISDIR(information.st_mode)) {
    sync_error = errno == 0 ? ENOTDIR : errno;
    (void)close(descriptor);
    Rf_error("Export transaction sync path is not an ordinary directory: %s.",
             strerror(sync_error));
  }
  do {
    sync_result = fsync(descriptor);
  } while (sync_result == -1 && errno == EINTR);
  sync_error = sync_result == -1 ? errno : 0;
  close_error = close(descriptor) == -1 ? errno : 0;
  if (close_error != 0) {
    Rf_error("Could not close the export transaction directory: %s.",
             strerror(close_error));
  }
  if (sync_error == EINVAL || sync_error == ENOTSUP ||
      sync_error == EOPNOTSUPP) {
    return ScalarLogical(0);
  }
  if (sync_error != 0) {
    Rf_error("Could not synchronize the export transaction directory: %s.",
             strerror(sync_error));
  }
  return ScalarLogical(1);
#endif
}

static const R_CallMethodDef cellucid_call_methods[] = {
    {"cellucid_export_lock_acquire",
     (DL_FUNC)&cellucid_export_lock_acquire, 1},
    {"cellucid_export_lock_release",
     (DL_FUNC)&cellucid_export_lock_release, 1},
    {"cellucid_export_lock_drain",
     (DL_FUNC)&cellucid_export_lock_drain, 0},
    {"cellucid_process_handle_count",
     (DL_FUNC)&cellucid_process_handle_count, 0},
    {"cellucid_export_path_info",
     (DL_FUNC)&cellucid_export_path_info, 1},
    {"cellucid_export_transaction_id",
     (DL_FUNC)&cellucid_export_transaction_id, 0},
    {"cellucid_export_write_journal",
     (DL_FUNC)&cellucid_export_write_journal, 2},
    {"cellucid_export_sync_directory",
     (DL_FUNC)&cellucid_export_sync_directory, 1},
    {NULL, NULL, 0}};

void attribute_visible R_init_cellucid(DllInfo *dll) {
  cellucid_export_lock_tag = Rf_install("cellucid_export_lock_handle");
#ifndef _WIN32
  cellucid_export_lock_registry_pid = getpid();
#endif
  R_registerRoutines(dll, NULL, cellucid_call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
