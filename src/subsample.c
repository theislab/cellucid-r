/*
 * Subsampling: the cells a Cellucid subsample keeps, by the subsampling
 * specification (docs: "Subsampling specification"), exactly as the web app
 * (assets/js/data/subsample.js) and Python (cellucid.subsample_cells) keep
 * them; tests/testthat/test-subsample-conformance.R holds all three to one
 * fixture. Hashes are uint32 arithmetic, allocations exact 64-bit integer
 * arithmetic, and every float operation of the spatial grid one IEEE 754
 * binary64 operation in the specified order (never contracted into a fused
 * multiply-add). Positions cross to R 0-based, as doubles (a dataset holds up
 * to 2^32 - 1 cells, more than an R integer counts).
 */

#pragma STDC FP_CONTRACT OFF

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static uint32_t fmix32(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6bu;
  h ^= h >> 13;
  h *= 0xc2b2ae35u;
  h ^= h >> 16;
  return h;
}

static uint32_t seed_key(uint32_t seed) { return fmix32(seed ^ 0x9e3779b9u); }

static uint32_t as_uint32(SEXP value, const char *label) {
  double x = asReal(value);
  if (!(x >= 0 && x <= 4294967295.0) || x != floor(x)) {
    error("%s must be an integer in [0, 2^32 - 1]", label);
  }
  return (uint32_t)x;
}

static int compare_u64(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a;
  uint64_t y = *(const uint64_t *)b;
  return x < y ? -1 : x > y;
}

static int compare_u32(const void *a, const void *b) {
  uint32_t x = *(const uint32_t *)a;
  uint32_t y = *(const uint32_t *)b;
  return x < y ? -1 : x > y;
}

/* The k items of `items` (count m) of lowest rank, into `out` ascending. */
static void lowest(const uint32_t *items, size_t m, size_t k, uint32_t key, uint32_t *out) {
  if (k >= m) {
    memcpy(out, items, m * sizeof(uint32_t));
    qsort(out, m, sizeof(uint32_t), compare_u32);
    return;
  }
  uint64_t *keys = (uint64_t *)R_alloc(m, sizeof(uint64_t));
  for (size_t i = 0; i < m; i++) keys[i] = ((uint64_t)fmix32(items[i] ^ key) << 32) | items[i];
  qsort(keys, m, sizeof(uint64_t), compare_u64);
  for (size_t i = 0; i < k; i++) out[i] = (uint32_t)(keys[i] & 0xffffffffu);
  qsort(out, k, sizeof(uint32_t), compare_u32);
}

static SEXP positions_vector(const uint32_t *values, size_t n) {
  SEXP out = PROTECT(allocVector(REALSXP, (R_xlen_t)n));
  double *p = REAL(out);
  for (size_t i = 0; i < n; i++) p[i] = (double)values[i];
  UNPROTECT(1);
  return out;
}

SEXP cellucid_subsample_rank(SEXP seed, SEXP positions) {
  uint32_t key = seed_key(as_uint32(seed, "seed"));
  R_xlen_t n = XLENGTH(positions);
  SEXP out = PROTECT(allocVector(REALSXP, n));
  for (R_xlen_t i = 0; i < n; i++) {
    double x = REAL(positions)[i];
    if (!(x >= 0 && x <= 4294967295.0) || x != floor(x)) error("positions must be integers in [0, 2^32 - 1]");
    REAL(out)[i] = (double)fmix32((uint32_t)x ^ key);
  }
  UNPROTECT(1);
  return out;
}

SEXP cellucid_subsample_uniform(SEXP n_obs, SEXP k, SEXP seed) {
  uint32_t n = as_uint32(n_obs, "n_obs");
  uint32_t kk = as_uint32(k, "k");
  uint32_t *items = (uint32_t *)R_alloc(n, sizeof(uint32_t));
  for (uint32_t i = 0; i < n; i++) items[i] = i;
  uint32_t *out = (uint32_t *)R_alloc(kk, sizeof(uint32_t));
  lowest(items, n, kk, seed_key(as_uint32(seed, "seed")), out);
  return positions_vector(out, kk);
}

SEXP cellucid_subsample_blocks(SEXP n_obs, SEXP k, SEXP seed, SEXP block_size) {
  uint32_t n = as_uint32(n_obs, "n_obs");
  uint32_t kk = as_uint32(k, "k");
  uint64_t size = as_uint32(block_size, "block_size");
  if (size == 0) error("block_size must be at least 1");
  uint32_t key = seed_key(as_uint32(seed, "seed"));
  uint64_t count = ((uint64_t)n + size - 1) / size;
  uint64_t *keys = (uint64_t *)R_alloc(count, sizeof(uint64_t));
  for (uint64_t j = 0; j < count; j++) keys[j] = ((uint64_t)fmix32((uint32_t)j ^ key) << 32) | j;
  qsort(keys, count, sizeof(uint64_t), compare_u64);
  uint32_t *out = (uint32_t *)R_alloc(kk, sizeof(uint32_t));
  size_t at = 0;
  for (uint64_t r = 0; r < count && at < kk; r++) {
    uint64_t start = (keys[r] & 0xffffffffu) * size;
    uint64_t length = size;
    if (start + length > n) length = n - start;
    for (uint64_t i = 0; i < length && at < kk; i++) out[at++] = (uint32_t)(start + i);
  }
  qsort(out, kk, sizeof(uint32_t), compare_u32);
  return positions_vector(out, kk);
}

typedef struct {
  uint64_t remainder;
  size_t stratum;
} remainder_entry;

static int compare_remainder(const void *a, const void *b) {
  const remainder_entry *x = (const remainder_entry *)a;
  const remainder_entry *y = (const remainder_entry *)b;
  if (x->remainder != y->remainder) return x->remainder > y->remainder ? -1 : 1;
  return x->stratum < y->stratum ? -1 : x->stratum > y->stratum;
}

/* Split `total` over strata of `capacity` each (none zero), as specified. */
static void allocate(const uint64_t *capacity, size_t m, uint64_t total, int balanced, uint64_t *out) {
  uint64_t sum = 0;
  for (size_t s = 0; s < m; s++) sum += capacity[s];
  for (size_t s = 0; s < m; s++) out[s] = 0;
  if (total == 0 || m == 0) return;
  if (!balanced) {
    remainder_entry *entries = (remainder_entry *)R_alloc(m, sizeof(remainder_entry));
    uint64_t given = 0;
    for (size_t s = 0; s < m; s++) {
      /* total and capacity are below 2^32: the product fits in 64 bits. */
      uint64_t product = total * capacity[s];
      out[s] = product / sum;
      entries[s].remainder = product % sum;
      entries[s].stratum = s;
      given += out[s];
    }
    qsort(entries, m, sizeof(remainder_entry), compare_remainder);
    for (uint64_t i = 0; i < total - given; i++) out[entries[i].stratum]++;
    return;
  }
  uint64_t *sorted = (uint64_t *)R_alloc(m, sizeof(uint64_t));
  memcpy(sorted, capacity, m * sizeof(uint64_t));
  qsort(sorted, m, sizeof(uint64_t), compare_u64);
  uint64_t level = 0, used = 0, remaining = m;
  for (size_t i = 0; i < m; i++) {
    uint64_t step = sorted[i] - level;
    if (used + step * remaining > total) break;
    used += step * remaining;
    level = sorted[i];
    remaining--;
  }
  if (remaining > 0) level += (total - used) / remaining;
  uint64_t given = 0;
  for (size_t s = 0; s < m; s++) {
    out[s] = capacity[s] < level ? capacity[s] : level;
    given += out[s];
  }
  for (size_t s = 0; s < m && given < total; s++) {
    if (capacity[s] > level) {
      out[s]++;
      given++;
    }
  }
}

/*
 * The kept cells of a split over strata: `strata` (0-based, one per cell,
 * each below n_strata), k cells, `minimum` first to each occupied stratum,
 * the rest by `allocation` ("proportional" or "balanced"), the lowest-ranked
 * cells of each stratum.
 */
SEXP cellucid_subsample_split(SEXP strata, SEXP n_strata, SEXP k, SEXP seed, SEXP allocation, SEXP minimum) {
  R_xlen_t n = XLENGTH(strata);
  uint32_t m = as_uint32(n_strata, "n_strata");
  uint64_t kk = as_uint32(k, "k");
  uint64_t floor_count = as_uint32(minimum, "minimum");
  int balanced = strcmp(CHAR(STRING_ELT(allocation, 0)), "balanced") == 0;
  uint32_t key = seed_key(as_uint32(seed, "seed"));
  uint64_t *sizes = (uint64_t *)R_alloc(m, sizeof(uint64_t));
  memset(sizes, 0, m * sizeof(uint64_t));
  const double *values = REAL(strata);
  for (R_xlen_t i = 0; i < n; i++) {
    double s = values[i];
    if (!(s >= 0 && s < m) || s != floor(s)) error("stratum %g of cell %lld is outside [0, %u)", s, (long long)i, m);
    sizes[(uint32_t)s]++;
  }
  uint64_t *base = (uint64_t *)R_alloc(m, sizeof(uint64_t));
  uint64_t *capacity = (uint64_t *)R_alloc(m, sizeof(uint64_t));
  size_t *index = (size_t *)R_alloc(m, sizeof(size_t));
  size_t positive = 0, occupied = 0;
  uint64_t floor_total = 0;
  for (uint32_t s = 0; s < m; s++) {
    if (sizes[s] == 0) continue;
    occupied++;
    base[s] = sizes[s] < floor_count ? sizes[s] : floor_count;
    floor_total += base[s];
    if (sizes[s] - base[s] > 0) {
      capacity[positive] = sizes[s] - base[s];
      index[positive++] = s;
    }
  }
  if (floor_total > kk) {
    error("a subsample of %llu cells cannot hold %llu cells of each of %llu categories",
          (unsigned long long)kk, (unsigned long long)floor_count, (unsigned long long)occupied);
  }
  uint64_t *extra = (uint64_t *)R_alloc(positive > 0 ? positive : 1, sizeof(uint64_t));
  allocate(capacity, positive, kk - floor_total, balanced, extra);
  uint64_t *quota = (uint64_t *)R_alloc(m, sizeof(uint64_t));
  for (uint32_t s = 0; s < m; s++) quota[s] = sizes[s] == 0 ? 0 : base[s];
  for (size_t i = 0; i < positive; i++) quota[index[i]] += extra[i];
  /* Cells grouped by stratum, ascending within each. */
  uint64_t *start = (uint64_t *)R_alloc((size_t)m + 1, sizeof(uint64_t));
  start[0] = 0;
  for (uint32_t s = 0; s < m; s++) start[s + 1] = start[s] + sizes[s];
  uint64_t *fill = (uint64_t *)R_alloc(m, sizeof(uint64_t));
  memcpy(fill, start, m * sizeof(uint64_t));
  uint32_t *grouped = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  for (R_xlen_t i = 0; i < n; i++) grouped[fill[(uint32_t)values[i]]++] = (uint32_t)i;
  uint32_t *out = (uint32_t *)R_alloc(kk, sizeof(uint32_t));
  size_t at = 0;
  for (uint32_t s = 0; s < m; s++) {
    if (quota[s] == 0) continue;
    lowest(grouped + start[s], sizes[s], quota[s], key, out + at);
    at += quota[s];
  }
  qsort(out, at, sizeof(uint32_t), compare_u32);
  return positions_vector(out, at);
}

/*
 * The grid box of each cell (`spatial`): coordinates (a row-major double
 * matrix of n rows and `dimensions` columns) taken as float32, each axis'
 * box floor((x - min) * grid / (max - min)), grid - 1 at the maximum, 0 on
 * an axis of one value; a cell with a non-finite coordinate in box
 * grid^dimensions.
 */
SEXP cellucid_subsample_boxes(SEXP coordinates, SEXP dimensions, SEXP grid_size) {
  int d = asInteger(dimensions);
  double grid = (double)as_uint32(grid_size, "grid_size");
  if (d < 1 || d > 3) error("a spatial subsample needs a 1-, 2- or 3-dimensional embedding");
  if (grid < 1) error("grid_size must be at least 1");
  R_xlen_t n = XLENGTH(coordinates) / d;
  double *values = (double *)R_alloc((size_t)XLENGTH(coordinates), sizeof(double));
  for (R_xlen_t i = 0; i < XLENGTH(coordinates); i++) {
    volatile float narrow = (float)REAL(coordinates)[i];
    values[i] = (double)narrow;
  }
  double low[3] = {INFINITY, INFINITY, INFINITY};
  double high[3] = {-INFINITY, -INFINITY, -INFINITY};
  for (R_xlen_t i = 0; i < n; i++) {
    int finite = 1;
    for (int a = 0; a < d; a++) if (!isfinite(values[i * d + a])) finite = 0;
    if (!finite) continue;
    for (int a = 0; a < d; a++) {
      double x = values[i * d + a];
      if (x < low[a]) low[a] = x;
      if (x > high[a]) high[a] = x;
    }
  }
  SEXP out = PROTECT(allocVector(REALSXP, n));
  double unplaced = pow(grid, d);
  for (R_xlen_t i = 0; i < n; i++) {
    int finite = 1;
    for (int a = 0; a < d; a++) if (!isfinite(values[i * d + a])) finite = 0;
    if (!finite) {
      REAL(out)[i] = unplaced;
      continue;
    }
    double box = 0, scale = 1;
    for (int a = 0; a < d; a++) {
      volatile double width = high[a] - low[a];
      double b = 0;
      if (width > 0) {
        volatile double offset = values[i * d + a] - low[a];
        volatile double scaled = offset * grid;
        volatile double ratio = scaled / width;
        b = floor(ratio);
        if (b > grid - 1) b = grid - 1;
      }
      box += b * scale;
      scale *= grid;
    }
    REAL(out)[i] = box;
  }
  UNPROTECT(1);
  return out;
}

/* The cell-set digest of the specification, as 16 lowercase hex digits. */
SEXP cellucid_subsample_digest(SEXP n_obs, SEXP kept) {
  uint32_t n = as_uint32(n_obs, "n_obs");
  int all = isNull(kept);
  uint64_t k = all ? n : (uint64_t)XLENGTH(kept);
  uint32_t a = 0x811c9dc5u ^ n;
  uint32_t b = 0x9e3779b9u ^ (uint32_t)k;
  for (uint64_t i = 0; i < k; i++) {
    uint32_t p = all ? (uint32_t)i : (uint32_t)REAL(kept)[i];
    a = (a ^ p) * 0x01000193u;
    uint32_t x = b ^ p;
    b = ((x << 13) | (x >> 19)) * 5u + 0xe6546b64u;
  }
  char text[17];
  snprintf(text, sizeof text, "%08x%08x", fmix32(a ^ (uint32_t)k), fmix32(b ^ n));
  return mkString(text);
}

/* A cell set as a session stores it: gaps minus one, LEB128. */
SEXP cellucid_subsample_encode(SEXP kept) {
  R_xlen_t k = XLENGTH(kept);
  unsigned char *bytes = (unsigned char *)R_alloc((size_t)k * 5 + 1, 1);
  size_t at = 0;
  double previous = -1;
  for (R_xlen_t i = 0; i < k; i++) {
    double p = REAL(kept)[i];
    if (!(p > previous) || p != floor(p) || p > 4294967295.0) error("a cell set lists positions strictly ascending");
    uint64_t gap = (uint64_t)(p - previous - 1);
    previous = p;
    while (gap >= 0x80) {
      bytes[at++] = (unsigned char)((gap & 0x7f) | 0x80);
      gap >>= 7;
    }
    bytes[at++] = (unsigned char)gap;
  }
  SEXP out = PROTECT(allocVector(RAWSXP, (R_xlen_t)at));
  if (at > 0) memcpy(RAW(out), bytes, at);
  UNPROTECT(1);
  return out;
}

SEXP cellucid_subsample_decode(SEXP payload, SEXP n_obs) {
  uint32_t n = as_uint32(n_obs, "n_obs");
  R_xlen_t length = XLENGTH(payload);
  const unsigned char *bytes = RAW(payload);
  double *positions = (double *)R_alloc(length > 0 ? (size_t)length : 1, sizeof(double));
  R_xlen_t count = 0;
  double previous = -1;
  uint64_t gap = 0;
  int shift = 0;
  for (R_xlen_t i = 0; i < length; i++) {
    gap |= (uint64_t)(bytes[i] & 0x7f) << shift;
    if (bytes[i] & 0x80) {
      shift += 7;
      if (shift > 35) error("a stored cell set is malformed");
      continue;
    }
    previous = previous + (double)gap + 1;
    if (previous >= n) error("a stored cell set names a cell beyond the dataset's %u cells", n);
    positions[count++] = previous;
    gap = 0;
    shift = 0;
  }
  if (shift != 0) error("a stored cell set is truncated");
  SEXP out = PROTECT(allocVector(REALSXP, count));
  if (count > 0) memcpy(REAL(out), positions, (size_t)count * sizeof(double));
  UNPROTECT(1);
  return out;
}

/*
 * The shortest decimal that names a share x in (0, 1]: the fewest
 * significant digits whose text is read back as x, and among those the
 * nearest to x -- the text Python's repr and the web app's String(x) give,
 * which the share's rational is taken from (R/subsample.R). R's own
 * conversion of text to a double is not correctly rounded on every build
 * (it reads "0.09946616017259657" one ulp low), so the round trip is tested
 * with the C library's strtod, which is. At each length the correctly
 * rounded text is tried first; when it does not read back, its neighbour on
 * x's side may (x a power of two: its rounding interval is narrower below
 * it, so "5.960464477539062e-08" is another double and
 * "5.960464477539063e-08" is 2^-24).
 */
SEXP cellucid_shortest_decimal(SEXP x) {
  if (TYPEOF(x) != REALSXP || XLENGTH(x) != 1) error("a share must be one number");
  double value = REAL(x)[0];
  if (!(value > 0 && value <= 1)) error("a share must be in (0, 1]");
  char text[64];
  for (int digits = 1; digits <= 17; digits++) {
    snprintf(text, sizeof text, "%.*e", digits - 1, value);
    if (strtod(text, NULL) == value) return mkString(text);
    /* The neighbour on x's side: one unit of the last digit. */
    char *e = strchr(text, 'e');
    int exponent = atoi(e + 1);
    long long mantissa = 0;
    for (const char *p = text; p < e; p++) {
      if (*p >= '0' && *p <= '9') mantissa = mantissa * 10 + (*p - '0');
    }
    long long neighbour = strtod(text, NULL) < value ? mantissa + 1 : mantissa - 1;
    long long lowest = 1;
    for (int i = 1; i < digits; i++) lowest *= 10;
    if (neighbour >= lowest * 10) {
      /* 9.99 up one unit is 1.00 of the next power of ten: fewer digits, tried already. */
      continue;
    }
    if (neighbour < lowest && digits > 1) continue; /* 1.00 down one unit: fewer digits, tried already */
    char shown[32];
    snprintf(shown, sizeof shown, "%lld", neighbour);
    char candidate[64];
    if (digits > 1) {
      snprintf(candidate, sizeof candidate, "%c.%se%+03d", shown[0], shown + 1, exponent);
    } else {
      snprintf(candidate, sizeof candidate, "%ce%+03d", shown[0], exponent);
    }
    if (strtod(candidate, NULL) == value) return mkString(candidate);
  }
  snprintf(text, sizeof text, "%.16e", value);
  return mkString(text);
}
