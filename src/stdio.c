/*
 * The process's own standard input and output, for cellucid_mcp() (R/mcp.R):
 * a local MCP server speaks JSON-RPC on them while R's event loop keeps
 * running a viewer's server and the agents' requests (later::run_now()), so
 * reading may never block for longer than it is asked to wait, and nothing R
 * prints may reach the protocol's stream.
 *
 * - cellucid_stdin_read(timeout): the bytes standard input holds now, waiting
 *   at most `timeout` seconds for some; NULL when none came in that time,
 *   raw(0) at the end of the input (the app closed it). Linux and macOS ask
 *   poll() and read() what is there; on Windows an app's pipe is asked with
 *   PeekNamedPipe() (a pipe has no wait of its own: it is asked every 10 ms
 *   until the time is up) and read with ReadFile() for exactly what it holds,
 *   and a file read at once. A Windows console is refused before reading
 *   (cellucid_stdin_console()): its input is key events, which no app sends.
 * - cellucid_stdout_write(bytes): writes all the bytes to standard output, by
 *   its descriptor (Windows: its handle), past R's console and any sink();
 *   FALSE when it is closed.
 */

#include <R.h>
#include <Rinternals.h>
#include <string.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <errno.h>
#include <poll.h>
#include <unistd.h>
#endif

/* What one read takes at most; the rest is read by the next call. */
#define STDIN_CHUNK_BYTES 65536

static unsigned char stdin_chunk[STDIN_CHUNK_BYTES];

static SEXP bytes_of(const unsigned char *bytes, size_t count) {
  SEXP out = PROTECT(allocVector(RAWSXP, (R_xlen_t)count));
  if (count > 0) memcpy(RAW(out), bytes, count);
  UNPROTECT(1);
  return out;
}

SEXP cellucid_stdin_read(SEXP timeout) {
  double seconds = asReal(timeout);
  if (!R_FINITE(seconds) || seconds < 0) seconds = 0;
  int milliseconds = (int)(seconds * 1000.0);
#if defined(_WIN32)
  HANDLE handle = GetStdHandle(STD_INPUT_HANDLE);
  if (handle == NULL || handle == INVALID_HANDLE_VALUE) return bytes_of(NULL, 0);
  DWORD type = GetFileType(handle);
  DWORD wanted = STDIN_CHUNK_BYTES;
  if (type == FILE_TYPE_PIPE) {
    int waited = 0;
    for (;;) {
      DWORD held = 0;
      /* A pipe whose writer closed it fails here: the end of the input. */
      if (!PeekNamedPipe(handle, NULL, 0, NULL, &held, NULL)) return bytes_of(NULL, 0);
      if (held > 0) {
        wanted = held < STDIN_CHUNK_BYTES ? held : STDIN_CHUNK_BYTES;
        break;
      }
      if (waited >= milliseconds) return R_NilValue;
      Sleep(10);
      waited += 10;
    }
  } else if (type == FILE_TYPE_CHAR) {
    /* A console: cellucid_mcp() refuses it before reading. */
    return bytes_of(NULL, 0);
  }
  DWORD count = 0;
  if (!ReadFile(handle, stdin_chunk, wanted, &count, NULL)) return bytes_of(NULL, 0);
  return bytes_of(stdin_chunk, (size_t)count);
#else
  struct pollfd input;
  input.fd = 0;
  input.events = POLLIN;
  input.revents = 0;
  int ready = poll(&input, 1, milliseconds);
  if (ready < 0) return errno == EINTR ? R_NilValue : bytes_of(NULL, 0);
  if (ready == 0) return R_NilValue;
  ssize_t count = read(0, stdin_chunk, STDIN_CHUNK_BYTES);
  if (count < 0) return (errno == EINTR || errno == EAGAIN) ? R_NilValue : bytes_of(NULL, 0);
  return bytes_of(stdin_chunk, (size_t)count);
#endif
}

/* Whether standard input is a Windows console (FALSE elsewhere: a terminal
 * on Linux and macOS is polled like a pipe). */
SEXP cellucid_stdin_console(void) {
#if defined(_WIN32)
  HANDLE handle = GetStdHandle(STD_INPUT_HANDLE);
  if (handle == NULL || handle == INVALID_HANDLE_VALUE) return ScalarLogical(0);
  return ScalarLogical(GetFileType(handle) == FILE_TYPE_CHAR);
#else
  return ScalarLogical(0);
#endif
}

SEXP cellucid_stdout_write(SEXP bytes) {
  if (TYPEOF(bytes) != RAWSXP) error("bytes must be a raw vector");
  const unsigned char *data = RAW(bytes);
  R_xlen_t left = XLENGTH(bytes);
#if defined(_WIN32)
  HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
  if (handle == NULL || handle == INVALID_HANDLE_VALUE) return ScalarLogical(0);
  while (left > 0) {
    DWORD piece = left > 1073741824 ? 1073741824 : (DWORD)left;
    DWORD written = 0;
    if (!WriteFile(handle, data, piece, &written, NULL)) return ScalarLogical(0);
    data += written;
    left -= (R_xlen_t)written;
  }
#else
  while (left > 0) {
    ssize_t written = write(1, data, (size_t)left);
    if (written < 0) {
      if (errno == EINTR) continue;
      if (errno == EAGAIN) {
        struct pollfd output;
        output.fd = 1;
        output.events = POLLOUT;
        output.revents = 0;
        poll(&output, 1, 100);
        continue;
      }
      return ScalarLogical(0);
    }
    data += written;
    left -= (R_xlen_t)written;
  }
#endif
  return ScalarLogical(1);
}
