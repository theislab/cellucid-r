/*
 * What the gene-major passes share (src/chunks.c: a matrix in R's memory;
 * src/backed.c: one held outside it): the little-endian stores of a chunk's
 * fixed-width values, the sort of a gene's (row, value) pairs by row, and
 * memory-mapped scratch files.
 */

#ifndef CELLUCID_GENE_MAJOR_H
#define CELLUCID_GENE_MAJOR_H

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

/* Little-endian stores and loads that do not depend on the host's byte order. */
static inline void put32(unsigned char *at, uint32_t v) {
  at[0] = (unsigned char)v;
  at[1] = (unsigned char)(v >> 8);
  at[2] = (unsigned char)(v >> 16);
  at[3] = (unsigned char)(v >> 24);
}

static inline void put64(unsigned char *at, uint64_t v) {
  put32(at, (uint32_t)v);
  put32(at + 4, (uint32_t)(v >> 32));
}

static inline uint32_t get32(const unsigned char *at) {
  return (uint32_t)at[0] | ((uint32_t)at[1] << 8) | ((uint32_t)at[2] << 16) | ((uint32_t)at[3] << 24);
}

static inline uint64_t get64(const unsigned char *at) {
  return (uint64_t)get32(at) | ((uint64_t)get32(at + 4) << 32);
}

/*
 * R's NA and NaN are NaNs with platform-dependent payloads and sign: every one
 * is written as numpy's canonical quiet NaN (0x7FC00000, 0x7FF8000000000000),
 * so the bytes do not depend on the platform.
 */
static inline void put_f4(unsigned char *at, double v) {
  if (isnan(v)) {
    put32(at, 0x7FC00000u);
    return;
  }
  float f = (float)v;
  uint32_t bits;
  memcpy(&bits, &f, 4);
  put32(at, bits);
}

static inline void put_f8(unsigned char *at, double v) {
  if (isnan(v)) {
    put64(at, 0x7FF8000000000000ull);
    return;
  }
  uint64_t bits;
  memcpy(&bits, &v, 8);
  put64(at, bits);
}

static inline double get_f4(const unsigned char *at) {
  uint32_t bits = get32(at);
  float f;
  memcpy(&f, &bits, 4);
  return (double)f;
}

static inline double get_f8(const unsigned char *at) {
  uint64_t bits = get64(at);
  double v;
  memcpy(&v, &bits, 8);
  return v;
}

/* One (row, value) pair of a gene (or a graph's row) being put in row order. */
struct renumbered { int64_t row; double value; };

/*
 * The `n` pairs in ascending row (rows are distinct and below 2^`bits`), into
 * `pairs`; `spare` is as long (unused below 64 pairs).
 */
void cellucid_sort_by_row(struct renumbered *pairs, struct renumbered *spare, R_xlen_t n, int bits);

/* A scratch file mapped into memory, removed as soon as it is mapped (Windows: when closed). */
struct scratch_map {
  unsigned char *bytes;
  size_t length;
#ifdef _WIN32
  HANDLE file, map;
#else
  int fd;
#endif
};

/* An unmapped scratch_map (unmapping it is then a no-op). */
void cellucid_scratch_init(struct scratch_map *m);
/* Create, size and map the scratch file at `path`; NULL, or what failed. */
const char *cellucid_map_scratch(struct scratch_map *m, const char *path, size_t length);
void cellucid_unmap_scratch(struct scratch_map *m);

#endif
