/*
 * The version of a file served as it is (R/viewer-server.R, .file_version()):
 * its size and a strong validator, read from the file opened, never from its
 * path, as cellucid-python's _file_version() reads them; and a range of the
 * file read with its version from the one file opened (.file_response()), as
 * cellucid-python answers a range from the descriptor it reads the tag from.
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
/* One file opened for reading: its descriptor (handle on Windows), and its
 * version read from it. */
typedef struct {
#if defined(_WIN32)
  HANDLE handle;
#else
  int descriptor;
#endif
  int64_t size, modified, changed, identity;
} opened_file;

/* Open the regular file at `path` and read its version from what was opened;
 * 0 when it is not a readable regular file (nothing left open). */
static int open_file(SEXP path, opened_file *file) {
  if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 || STRING_ELT(path, 0) == NA_STRING) {
    error("file version: path must be one string");
  }
#if defined(_WIN32)
  /* The name in UTF-16, from UTF-8 whatever the session's code page (a path
   * with characters outside it would otherwise name no file). */
  const char *name = translateCharUTF8(STRING_ELT(path, 0));
  int wide_length = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, name, -1, NULL, 0);
  if (wide_length <= 0) return 0;
  wchar_t *wide = (wchar_t *)R_alloc((size_t)wide_length, sizeof(wchar_t));
  if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, name, -1, wide, wide_length) != wide_length) {
    return 0;
  }
  /* Opened as Python's open() opens it: shared with writers, renames and
   * deletes, so serving never locks the file. A folder cannot be opened
   * without FILE_FLAG_BACKUP_SEMANTICS, so it is not a file here either. */
  HANDLE handle = CreateFileW(wide, GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL,
                              OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
  if (handle == INVALID_HANDLE_VALUE) return 0;
  BY_HANDLE_FILE_INFORMATION information;
  FILE_BASIC_INFO basic;
  BOOL answered = GetFileInformationByHandle(handle, &information) &&
                  GetFileInformationByHandleEx(handle, FileBasicInfo, &basic, (DWORD)sizeof(basic));
  if (!answered || (information.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
    CloseHandle(handle);
    return 0;
  }
  file->handle = handle;
  file->size = (int64_t)(((uint64_t)information.nFileSizeHigh << 32) | information.nFileSizeLow);
  /* The modification time in nanoseconds since 1970, as Python's st_mtime_ns
   * states it; the change time in NTFS's 100 ns ticks since 1601, as
   * cellucid-python states it. */
  file->modified = ((int64_t)basic.LastWriteTime.QuadPart - INT64_C(116444736000000000)) * 100;
  file->changed = (int64_t)basic.ChangeTime.QuadPart;
  file->identity = (int64_t)(((uint64_t)information.nFileIndexHigh << 32) | information.nFileIndexLow);
  return 1;
#else
  /* The path in the session's native encoding, as the system takes it. */
  const char *native = translateChar(STRING_ELT(path, 0));
  int descriptor = open(native, O_RDONLY);
  if (descriptor < 0) return 0;
  struct stat status;
  if (fstat(descriptor, &status) != 0 || !S_ISREG(status.st_mode)) {
    close(descriptor);
    return 0;
  }
  file->descriptor = descriptor;
  file->size = (int64_t)status.st_size;
#if defined(__APPLE__)
  file->modified = (int64_t)status.st_mtimespec.tv_sec * 1000000000 + status.st_mtimespec.tv_nsec;
  file->changed = (int64_t)status.st_ctimespec.tv_sec * 1000000000 + status.st_ctimespec.tv_nsec;
#else
  file->modified = (int64_t)status.st_mtim.tv_sec * 1000000000 + status.st_mtim.tv_nsec;
  file->changed = (int64_t)status.st_ctim.tv_sec * 1000000000 + status.st_ctim.tv_nsec;
#endif
  file->identity = (int64_t)status.st_ino;
  return 1;
#endif
}

static void close_file(opened_file *file) {
#if defined(_WIN32)
  CloseHandle(file->handle);
#else
  close(file->descriptor);
#endif
}

/* `length` bytes from `offset` of an opened file into `out`; the bytes read
 * (fewer when the file ends first: it shrank since its size was read). */
static int64_t read_file(opened_file *file, int64_t offset, int64_t length, unsigned char *out) {
  int64_t done = 0;
  while (done < length) {
#if defined(_WIN32)
    LARGE_INTEGER at;
    at.QuadPart = offset + done;
    DWORD want = (DWORD)((length - done) > 0x40000000 ? 0x40000000 : (length - done));
    DWORD got = 0;
    if (!SetFilePointerEx(file->handle, at, NULL, FILE_BEGIN) || !ReadFile(file->handle, out + done, want, &got, NULL)) break;
    if (got == 0) break;
    done += (int64_t)got;
#else
    size_t want = (size_t)((length - done) > 0x40000000 ? 0x40000000 : (length - done));
    ssize_t got = pread(file->descriptor, out + done, want, (off_t)(offset + done));
    if (got <= 0) break;
    done += (int64_t)got;
#endif
  }
  return done;
}

/* list(size, etag) of the regular file at `path`, or NULL. */
SEXP cellucid_file_version(SEXP path) {
  opened_file file;
  if (!open_file(path, &file)) return R_NilValue;
  close_file(&file);
  return version_of(file.size, file.modified, file.changed, file.identity);
}

/*
 * list(size, etag, bytes) of the regular file at `path`, or NULL: its version
 * and `length` bytes from `offset` (0-based), read from the one file opened,
 * so a file replaced or rewritten meanwhile never sends its new bytes under
 * the old version. `bytes` is shorter than asked when the file ended first.
 */
SEXP cellucid_file_range(SEXP path, SEXP offset, SEXP length) {
  double first = asReal(offset);
  double count = asReal(length);
  if (!R_FINITE(first) || !R_FINITE(count) || first < 0 || count < 0 || first + count > 9007199254740992.0) {
    error("file range: offset and length must be whole numbers from 0");
  }
  /* Made before the file is opened: an allocation that fails leaves nothing open. */
  PROTECT_INDEX at;
  SEXP bytes;
  PROTECT_WITH_INDEX(bytes = allocVector(RAWSXP, (R_xlen_t)count), &at);
  opened_file file;
  if (!open_file(path, &file)) {
    UNPROTECT(1);
    return R_NilValue;
  }
  int64_t got = read_file(&file, (int64_t)first, (int64_t)count, RAW(bytes));
  close_file(&file);
  if (got < (int64_t)count) REPROTECT(bytes = xlengthgets(bytes, (R_xlen_t)got), at);
  SEXP version = PROTECT(version_of(file.size, file.modified, file.changed, file.identity));
  SEXP out = PROTECT(allocVector(VECSXP, 3));
  SEXP names = PROTECT(allocVector(STRSXP, 3));
  SET_VECTOR_ELT(out, 0, VECTOR_ELT(version, 0));
  SET_VECTOR_ELT(out, 1, VECTOR_ELT(version, 1));
  SET_VECTOR_ELT(out, 2, bytes);
  SET_STRING_ELT(names, 0, mkChar("size"));
  SET_STRING_ELT(names, 1, mkChar("etag"));
  SET_STRING_ELT(names, 2, mkChar("bytes"));
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(4);
  return out;
}
