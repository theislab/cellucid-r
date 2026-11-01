/*
 * DEFLATE of many chunks at once, on every core (R/zarr-store.R): the chunks
 * of a written store, compressed a batch at a time. Each chunk is one zlib
 * stream (RFC 1950) at the store's level, made by zlib's compress2() exactly
 * as base R's memCompress(type = "gzip") makes one; the bytes are each
 * chunk's own whatever the number of threads, since every chunk is
 * compressed alone. The worker threads touch no R object: every input's
 * bytes and length are read, and every output's buffer is made into a raw
 * vector, on R's thread, before and after them (R's API is not thread-safe).
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

typedef struct {
  const Bytef **in;
  uLong *in_length;
  Bytef **out;
  uLongf *out_length;
  int *status;
  R_xlen_t n;
  int level;
  R_xlen_t next;
#ifdef _WIN32
  CRITICAL_SECTION lock;
#else
  pthread_mutex_t lock;
#endif
} deflate_job;

static R_xlen_t take_next(deflate_job *job) {
  R_xlen_t index;
#ifdef _WIN32
  EnterCriticalSection(&job->lock);
  index = job->next++;
  LeaveCriticalSection(&job->lock);
#else
  pthread_mutex_lock(&job->lock);
  index = job->next++;
  pthread_mutex_unlock(&job->lock);
#endif
  return index;
}

static void deflate_some(deflate_job *job) {
  for (;;) {
    R_xlen_t index = take_next(job);
    if (index >= job->n) return;
    uLongf room = compressBound(job->in_length[index]);
    Bytef *buffer = (Bytef *)malloc(room > 0 ? room : 1);
    if (buffer == NULL) {
      job->status[index] = Z_MEM_ERROR;
      continue;
    }
    job->status[index] = compress2(buffer, &room, job->in[index], job->in_length[index], job->level);
    job->out[index] = buffer;
    job->out_length[index] = room;
  }
}

#ifdef _WIN32
static DWORD WINAPI deflate_worker(LPVOID argument) {
  deflate_some((deflate_job *)argument);
  return 0;
}
#else
static void *deflate_worker(void *argument) {
  deflate_some((deflate_job *)argument);
  return NULL;
}
#endif

/*
 * `chunks`: a list of raw vectors; `level`: the deflate level; `threads`: how
 * many threads compress (the caller's count of cores; one compresses on R's
 * thread alone). Returns a list of raw vectors, each chunk's zlib stream.
 */
SEXP cellucid_deflate_chunks(SEXP chunks, SEXP level, SEXP threads) {
  if (TYPEOF(chunks) != VECSXP) error("chunks must be a list of raw vectors");
  R_xlen_t n = XLENGTH(chunks);
  deflate_job job;
  job.n = n;
  job.level = asInteger(level);
  job.next = 0;
  job.in = (const Bytef **)R_alloc(n > 0 ? n : 1, sizeof(Bytef *));
  job.in_length = (uLong *)R_alloc(n > 0 ? n : 1, sizeof(uLong));
  job.out = (Bytef **)R_alloc(n > 0 ? n : 1, sizeof(Bytef *));
  job.out_length = (uLongf *)R_alloc(n > 0 ? n : 1, sizeof(uLongf));
  job.status = (int *)R_alloc(n > 0 ? n : 1, sizeof(int));
  for (R_xlen_t index = 0; index < n; index++) {
    SEXP chunk = VECTOR_ELT(chunks, index);
    if (TYPEOF(chunk) != RAWSXP) error("chunk %lld is not a raw vector", (long long)index + 1);
    /* zlib counts a buffer in uLong: 32 bits on Windows. A chunk is far smaller. */
    if ((double)XLENGTH(chunk) > (double)((uLong)-1) / 2) error("chunk %lld is too large for zlib", (long long)index + 1);
    job.in[index] = (const Bytef *)RAW(chunk);
    job.in_length[index] = (uLong)XLENGTH(chunk);
    job.out[index] = NULL;
    job.out_length[index] = 0;
    job.status[index] = Z_OK;
  }
  int count = asInteger(threads);
  if (count == NA_INTEGER || count < 1) count = 1;
  if ((R_xlen_t)count > n) count = n > 0 ? (int)n : 1;
#ifdef _WIN32
  InitializeCriticalSection(&job.lock);
#else
  pthread_mutex_init(&job.lock, NULL);
#endif
  int started = 0;
  if (count > 1) {
#ifdef _WIN32
    HANDLE *handles = (HANDLE *)R_alloc(count - 1, sizeof(HANDLE));
    for (int t = 0; t < count - 1; t++) {
      handles[t] = CreateThread(NULL, 0, deflate_worker, &job, 0, NULL);
      if (handles[t] == NULL) break;
      started++;
    }
    deflate_some(&job);
    for (int t = 0; t < started; t++) {
      WaitForSingleObject(handles[t], INFINITE);
      CloseHandle(handles[t]);
    }
#else
    pthread_t *handles = (pthread_t *)R_alloc(count - 1, sizeof(pthread_t));
    for (int t = 0; t < count - 1; t++) {
      if (pthread_create(&handles[t], NULL, deflate_worker, &job) != 0) break;
      started++;
    }
    deflate_some(&job);
    for (int t = 0; t < started; t++) pthread_join(handles[t], NULL);
#endif
  } else {
    deflate_some(&job);
  }
#ifdef _WIN32
  DeleteCriticalSection(&job.lock);
#else
  pthread_mutex_destroy(&job.lock);
#endif
  int failed = 0;
  for (R_xlen_t index = 0; index < n; index++) {
    if (job.status[index] != Z_OK) failed = 1;
  }
  if (failed) {
    for (R_xlen_t index = 0; index < n; index++) free(job.out[index]);
    error("deflate: zlib could not compress a chunk");
  }
  SEXP out = PROTECT(allocVector(VECSXP, n));
  for (R_xlen_t index = 0; index < n; index++) {
    SEXP packed = allocVector(RAWSXP, (R_xlen_t)job.out_length[index]);
    SET_VECTOR_ELT(out, index, packed);
    if (job.out_length[index] > 0) memcpy(RAW(packed), job.out[index], job.out_length[index]);
    free(job.out[index]);
    job.out[index] = NULL;
  }
  UNPROTECT(1);
  return out;
}
