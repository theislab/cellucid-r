/*
 * Random bytes from the operating system's cryptographic generator: the
 * viewer's id and token (R/viewer-api.R). R's own generator is seedable
 * (set.seed() makes a token predictable) and drawing from it moves the
 * caller's .Random.seed, so an analysis would give other results depending
 * on whether a viewer was opened; this never touches either.
 */

#include <R.h>
#include <Rinternals.h>

#if defined(_WIN32)
#include <windows.h>
#include <bcrypt.h>
#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
#include <stdlib.h>
#else
#include <stdio.h>
#endif

SEXP cellucid_random_bytes(SEXP n) {
  double wanted = asReal(n);
  if (!R_FINITE(wanted) || wanted < 0 || wanted > (double)R_XLEN_T_MAX ||
      wanted != (double)(R_xlen_t)wanted) {
    error("random bytes: n must be a whole number of bytes");
  }
  R_xlen_t count = (R_xlen_t)wanted;
  SEXP out = PROTECT(allocVector(RAWSXP, count));
  if (count > 0) {
    unsigned char *bytes = RAW(out);
#if defined(_WIN32)
    R_xlen_t done = 0;
    while (done < count) {
      ULONG part = (ULONG)((count - done) > 0x40000000 ? 0x40000000 : (count - done));
      if (BCryptGenRandom(NULL, bytes + done, part, BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0) {
        UNPROTECT(1);
        error("random bytes: the system random generator (BCryptGenRandom) failed");
      }
      done += part;
    }
#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    arc4random_buf(bytes, (size_t)count);
#else
    FILE *source = fopen("/dev/urandom", "rb");
    size_t got = 0;
    if (source != NULL) {
      got = fread(bytes, 1, (size_t)count, source);
      fclose(source);
    }
    if (got != (size_t)count) {
      UNPROTECT(1);
      error("random bytes: /dev/urandom could not be read");
    }
#endif
  }
  UNPROTECT(1);
  return out;
}
