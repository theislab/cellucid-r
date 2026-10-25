/*
 * CRC-32 (ISO 3309, the ZIP checksum) of one raw vector, for the stored
 * entries of a .zarr.zip written by cellucid_prepare(format = "zarr.zip").
 * Base R has no CRC-32, and a per-byte loop in R is far too slow for
 * expression chunks.
 *
 * `crc` is the running value of the bytes before `bytes` (0 for none), so a
 * long entry can be checksummed in pieces.
 */

#include <R.h>
#include <Rinternals.h>

#include <stdint.h>

static uint32_t cellucid_crc32_table[256];
static int cellucid_crc32_ready = 0;

static void cellucid_crc32_init(void) {
  for (uint32_t n = 0; n < 256; n++) {
    uint32_t c = n;
    for (int k = 0; k < 8; k++) {
      c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
    }
    cellucid_crc32_table[n] = c;
  }
  cellucid_crc32_ready = 1;
}

SEXP cellucid_crc32(SEXP bytes, SEXP previous) {
  if (TYPEOF(bytes) != RAWSXP) {
    Rf_error("A CRC-32 is computed over a raw vector.");
  }
  if (!Rf_isReal(previous) || XLENGTH(previous) != 1) {
    Rf_error("The running CRC-32 must be one number.");
  }
  if (!cellucid_crc32_ready) cellucid_crc32_init();
  uint32_t crc = ((uint32_t)REAL(previous)[0]) ^ 0xFFFFFFFFu;
  const unsigned char *data = RAW(bytes);
  R_xlen_t count = XLENGTH(bytes);
  for (R_xlen_t i = 0; i < count; i++) {
    crc = cellucid_crc32_table[(crc ^ data[i]) & 0xFFu] ^ (crc >> 8);
  }
  return Rf_ScalarReal((double)(crc ^ 0xFFFFFFFFu));
}
