/*
 * The version of a file served as it is (R/viewer-server.R, .file_version()):
 * its size and a strong validator, read from the file opened, never from its
 * path, as cellucid-python's _file_version() reads them.
 *
 * The validator is the size, the modification time, the status-change time
 * and the file's identity (inode; on Windows the volume's file index), each
 * as hex: on Linux and macOS the very tag cellucid-python states. A file
 * rewritten in place gets a new modification and status-change time, even
 * when its modification time is put back after; one replaced by another
 * moved over it (even one whose modification time was copied with it: cp -p,
 * rsync -t) a new identity and status-change time. Base R cannot read these:
 * file.info() has no inode, and on Windows its ctime is the creation time,
 * which a rewrite in place keeps. Windows' status-change time is NTFS's
 * ChangeTime, read from the handle (GetFileInformationByHandleEx(
 * FileBasicInfo)), the call cellucid-python makes.
 */

/* GetFileInformationByHandleEx: Windows Vista on. Set before any header,
 * as the C library's own headers set a default. */
#if defined(_WIN32) && (!defined(_WIN32_WINNT) || _WIN32_WINNT < 0x0600)
#undef _WIN32_WINNT
#define _WIN32_WINNT 0x0600
#endif

#include <stdint.h>

#include <R.h>
#include <Rinternals.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

/* `value` as lower-case hex (with a leading '-' when negative), as Python
 * formats an integer with :x, appended at `out`; returns the end. */
static char *append_hex(char *out, int64_t value) {
  uint64_t magnitude = value < 0 ? (uint64_t)0 - (uint64_t)value : (uint64_t)value;
  char digits[16];
  int count = 0;
  if (value < 0) *out++ = '-';
  do {
    digits[count++] = "0123456789abcdef"[magnitude & 0xf];
    magnitude >>= 4;
  } while (magnitude != 0);
  while (count > 0) *out++ = digits[--count];
  return out;
}

/* list(size, etag): the size, and the validator of the four parts. */
static SEXP version_of(int64_t size, int64_t modified, int64_t changed, int64_t identity) {
  char tag[4 * 17 + 4 + 1];
  char *at = tag;
  *at++ = '"';
  at = append_hex(at, size);
  *at++ = '-';
  at = append_hex(at, modified);
  *at++ = '-';
  at = append_hex(at, changed);
  *at++ = '-';
  at = append_hex(at, identity);
  *at++ = '"';
  *at = '\0';
  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("size"));
  SET_STRING_ELT(names, 1, mkChar("etag"));
  SET_VECTOR_ELT(out, 0, ScalarReal((double)size));
  SET_VECTOR_ELT(out, 1, mkString(tag));
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

/* The version of the regular file at `path` (one string), or NULL when it
 * cannot be opened as one (gone, unreadable, a folder). */
SEXP cellucid_file_version(SEXP path) {
  if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 || STRING_ELT(path, 0) == NA_STRING) {
    error("file version: path must be one string");
  }
#if defined(_WIN32)
  /* The name in UTF-16, from UTF-8 whatever the session's code page (a path
   * with characters outside it would otherwise name no file). */
  const char *name = translateCharUTF8(STRING_ELT(path, 0));
  int wide_length = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, name, -1, NULL, 0);
  if (wide_length <= 0) return R_NilValue;
  wchar_t *wide = (wchar_t *)R_alloc((size_t)wide_length, sizeof(wchar_t));
  if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, name, -1, wide, wide_length) != wide_length) {
    return R_NilValue;
  }
  /* Opened as Python's open() opens it: shared with writers, renames and
   * deletes, so serving never locks the file. A folder cannot be opened
   * without FILE_FLAG_BACKUP_SEMANTICS, so it is not a file here either. */
  HANDLE handle = CreateFileW(wide, GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL,
                              OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
  if (handle == INVALID_HANDLE_VALUE) return R_NilValue;
  BY_HANDLE_FILE_INFORMATION file;
  FILE_BASIC_INFO basic;
  BOOL answered = GetFileInformationByHandle(handle, &file) &&
                  GetFileInformationByHandleEx(handle, FileBasicInfo, &basic, (DWORD)sizeof(basic));
  CloseHandle(handle);
  if (!answered || (file.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) return R_NilValue;
  int64_t size = (int64_t)(((uint64_t)file.nFileSizeHigh << 32) | file.nFileSizeLow);
  /* The modification time in nanoseconds since 1970, as Python's st_mtime_ns
   * states it; the change time in NTFS's 100 ns ticks since 1601, as
   * cellucid-python states it. */
  int64_t modified = ((int64_t)basic.LastWriteTime.QuadPart - INT64_C(116444736000000000)) * 100;
  int64_t changed = (int64_t)basic.ChangeTime.QuadPart;
  int64_t identity = (int64_t)(((uint64_t)file.nFileIndexHigh << 32) | file.nFileIndexLow);
  return version_of(size, modified, changed, identity);
#else
  /* The path in the session's native encoding, as the system takes it. */
  const char *native = translateChar(STRING_ELT(path, 0));
  int descriptor = open(native, O_RDONLY);
  if (descriptor < 0) return R_NilValue;
  struct stat status;
  int failed = fstat(descriptor, &status);
  close(descriptor);
  if (failed != 0 || !S_ISREG(status.st_mode)) return R_NilValue;
#if defined(__APPLE__)
  int64_t modified = (int64_t)status.st_mtimespec.tv_sec * 1000000000 + status.st_mtimespec.tv_nsec;
  int64_t changed = (int64_t)status.st_ctimespec.tv_sec * 1000000000 + status.st_ctimespec.tv_nsec;
#else
  int64_t modified = (int64_t)status.st_mtim.tv_sec * 1000000000 + status.st_mtim.tv_nsec;
  int64_t changed = (int64_t)status.st_ctim.tv_sec * 1000000000 + status.st_ctim.tv_nsec;
#endif
  return version_of((int64_t)status.st_size, modified, changed, (int64_t)status.st_ino);
#endif
}
