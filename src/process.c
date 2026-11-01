/*
 * Whether a process runs on this machine, never signalling it (R/output.R:
 * the staging folders a killed cellucid_prepare() left are removed by the
 * next one), as cellucid-python's _scratch._process_alive() asks: kill(pid,
 * 0) on Linux and macOS, OpenProcess() and its exit code on Windows (where a
 * signal ends a process). A process that runs but may not be asked (another
 * user's) counts as running.
 */

#include <R.h>
#include <Rinternals.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <errno.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#endif

SEXP cellucid_process_alive(SEXP pid) {
  double id = asReal(pid);
  if (!R_FINITE(id) || id <= 0 || id != (double)(long)id) return ScalarLogical(0);
#if defined(_WIN32)
  /* 0x1000: the access that asks a process's exit code and nothing else. */
  HANDLE handle = OpenProcess(0x1000, FALSE, (DWORD)id);
  if (handle == NULL) return ScalarLogical(GetLastError() == ERROR_ACCESS_DENIED);
  DWORD code = 0;
  BOOL answered = GetExitCodeProcess(handle, &code);
  CloseHandle(handle);
  return ScalarLogical(answered && code == 259 /* STILL_ACTIVE */);
#else
  if (kill((pid_t)id, 0) == 0) return ScalarLogical(1);
  return ScalarLogical(errno == EPERM);
#endif
}

/*
 * The physical memory of this machine in bytes (R/viewer-api.R: what a
 * served store keeps follows it, as cellucid-python's
 * _served_store.machine_memory_bytes() does): sysconf() on Linux and macOS,
 * GlobalMemoryStatusEx() on Windows; NA where neither answers.
 */
SEXP cellucid_physical_memory(void) {
#if defined(_WIN32)
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if (!GlobalMemoryStatusEx(&status)) return ScalarReal(NA_REAL);
  return ScalarReal((double)status.ullTotalPhys);
#else
  long pages = sysconf(_SC_PHYS_PAGES);
  long size = sysconf(_SC_PAGE_SIZE);
  if (pages <= 0 || size <= 0) return ScalarReal(NA_REAL);
  return ScalarReal((double)pages * (double)size);
#endif
}
