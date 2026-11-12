/*
 * The text JavaScript gives a number (ECMAScript Number::toString, which
 * JSON.stringify writes): the fewest significant digits whose text is read
 * back as the number, and among those the nearest to it, the even one when
 * two are as near; positional from 1e-6 up to 1e21, else with an exponent
 * ("1e+21", "1.5e-7"). Three places need these exact characters: a session
 * value's content hash (the SHA-256 of its canonical text, R/json.R), an
 * action command's JSON (R/viewer-protocol.R), and a subsample share's
 * decimal (R/subsample.R), where Python's repr gives the same digits.
 *
 * It is found in exact integer arithmetic and its digits are written by hand,
 * never through the C library: R's own conversion of text to a double is not
 * correctly rounded on every build (it reads "0.09946616017259657" one ulp
 * low), the C libraries' printf breaks a tie as each likes (Windows' writes
 * 2^-25 to 17 digits as "2.9802322387695313e-08", where the specification's
 * text, the even one, is "...312"), and the compiler of R 4.1's Rtools on
 * Windows checks printf formats against msvcrt's, warns on %llu, and R CMD
 * check counts that warning.
 *
 * x = m * 2^e exactly, and every text that reads back as x lies in its
 * rounding interval: from halfway down to the double below to halfway up to
 * the double above, both ends included when m is even (a text exactly
 * halfway is read as the even neighbor). In units of 2^(e - 2) x is 4m, the
 * interval [4m - 2, 4m + 2] (from 4m - 1 for a power of two, whose double
 * below is nearer). At each length of `digits` digits the decimals are
 * c * 10^p; with A / B = x / 10^p as two big integers, c = A / B rounded half
 * to even is the nearest, and it, or its neighbor on x's other side when it
 * is outside, is compared with the interval's ends exactly.
 */

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

/* An unsigned big integer, 32-bit limbs least significant first. The largest
 * number formed is below 2^1260 (x at least 2^-1074 scaled by 10^341, or up
 * to 2^1024 by 10^-p and a 64-bit digit count); 64 limbs hold 2^2048. */
#define BIG_LIMBS 64

typedef struct {
  uint32_t limb[BIG_LIMBS];
  int n; /* limbs in use; 0 is zero */
} big_number;

static void big_trim(big_number *a) {
  while (a->n > 0 && a->limb[a->n - 1] == 0) a->n--;
}

static void big_set(big_number *a, uint64_t value) {
  a->n = 0;
  while (value > 0) {
    a->limb[a->n++] = (uint32_t)(value & 0xffffffffu);
    value >>= 32;
  }
}

static void big_multiply(big_number *a, uint32_t factor) {
  uint64_t carry = 0;
  for (int i = 0; i < a->n; i++) {
    uint64_t product = (uint64_t)a->limb[i] * factor + carry;
    a->limb[i] = (uint32_t)(product & 0xffffffffu);
    carry = product >> 32;
  }
  if (carry > 0) {
    if (a->n == BIG_LIMBS) error("a number's decimal does not fit the exact arithmetic");
    a->limb[a->n++] = (uint32_t)carry;
  }
  big_trim(a);
}

static void big_times_ten_to(big_number *a, int power) {
  static const uint32_t small[9] = {1u, 10u, 100u, 1000u, 10000u, 100000u, 1000000u, 10000000u, 100000000u};
  for (; power >= 9; power -= 9) big_multiply(a, 1000000000u);
  if (power > 0) big_multiply(a, small[power]);
}

static void big_times_two_to(big_number *a, int power) {
  if (a->n == 0 || power == 0) return;
  int words = power / 32, bits = power % 32;
  if (a->n + words + 1 > BIG_LIMBS) error("a number's decimal does not fit the exact arithmetic");
  uint32_t spill = 0;
  if (bits > 0) {
    for (int i = 0; i < a->n; i++) {
      uint32_t limb = a->limb[i];
      a->limb[i] = (limb << bits) | spill;
      spill = limb >> (32 - bits);
    }
  }
  for (int i = a->n - 1; i >= 0; i--) a->limb[i + words] = a->limb[i];
  for (int i = 0; i < words; i++) a->limb[i] = 0;
  a->n += words;
  if (spill > 0) a->limb[a->n++] = spill;
}

/* out = a * factor (factor below 2^64). */
static void big_times_u64(const big_number *a, uint64_t factor, big_number *out) {
  uint32_t parts[2] = {(uint32_t)(factor & 0xffffffffu), (uint32_t)(factor >> 32)};
  int n_parts = parts[1] != 0 ? 2 : (parts[0] != 0 ? 1 : 0);
  if (n_parts == 0 || a->n == 0) {
    out->n = 0;
    return;
  }
  int n = a->n + n_parts;
  if (n > BIG_LIMBS) error("a number's decimal does not fit the exact arithmetic");
  for (int i = 0; i < n; i++) out->limb[i] = 0;
  for (int j = 0; j < n_parts; j++) {
    uint64_t carry = 0;
    for (int i = 0; i < a->n; i++) {
      uint64_t sum = (uint64_t)a->limb[i] * parts[j] + out->limb[i + j] + carry;
      out->limb[i + j] = (uint32_t)(sum & 0xffffffffu);
      carry = sum >> 32;
    }
    for (int k = a->n + j; carry > 0 && k < n; k++) {
      uint64_t sum = (uint64_t)out->limb[k] + carry;
      out->limb[k] = (uint32_t)(sum & 0xffffffffu);
      carry = sum >> 32;
    }
  }
  out->n = n;
  big_trim(out);
}

static int big_compare(const big_number *a, const big_number *b) {
  if (a->n != b->n) return a->n < b->n ? -1 : 1;
  for (int i = a->n - 1; i >= 0; i--) {
    if (a->limb[i] != b->limb[i]) return a->limb[i] < b->limb[i] ? -1 : 1;
  }
  return 0;
}

/* The sign of c * B - target. */
static int compare_scaled(uint64_t c, const big_number *b, const big_number *target) {
  big_number product;
  big_times_u64(b, c, &product);
  return big_compare(&product, target);
}

/* u * 2^twos * 10^tens, both powers at least 0. */
static void big_scaled(big_number *a, uint64_t u, int twos, int tens) {
  big_set(a, u);
  big_times_ten_to(a, tens);
  big_times_two_to(a, twos);
}

static uint64_t ten_to(int power) {
  uint64_t value = 1;
  while (power-- > 0) value *= 10;
  return value;
}

/* The decimal digits of `value` (a 64-bit number has at most 20) into `out`,
 * padded on the left with zeros to at least `places`, with its terminator;
 * returns how many digits it wrote. */
static int write_digits(uint64_t value, int places, char *out) {
  int length = 0;
  for (uint64_t left = value; left > 0; left /= 10) length++;
  if (length < places) length = places;
  if (length < 1) length = 1;
  out[length] = '\0';
  for (int i = length - 1; i >= 0; i--) {
    out[i] = (char)('0' + (int)(value % 10));
    value /= 10;
  }
  return length;
}

/*
 * The shortest decimal of a positive finite x: its digits (no trailing zero)
 * into `digits` and the position of its decimal point, `point`: x is
 * 0.<digits> * 10^point (ECMAScript's n).
 */
static void shortest_digits(double value, char *digits, int *point) {
  /* An integer below 2^53 is its own shortest decimal: the doubles near it
   * are at most 1 apart, so no other decimal (one ending in fewer digits is
   * another integer) reads back as it. */
  if (value < 9007199254740992.0 && value == floor(value)) {
    char text[24];
    int length = write_digits((uint64_t)value, 1, text);
    *point = length;
    while (length > 1 && text[length - 1] == '0') length--;
    memcpy(digits, text, (size_t)length);
    digits[length] = '\0';
    return;
  }
  /* value = m * 2^e: m below 2^53, e at least -1074 (a subnormal's m is
   * smaller, its e -1074). */
  int e;
  double fraction = frexp(value, &e);
  uint64_t m = (uint64_t)ldexp(fraction, 53);
  e -= 53;
  if (e < -1074) {
    m >>= -1074 - e;
    e = -1074;
  }
  int s = e - 2; /* value is 4m * 2^s */
  int even = (m % 2) == 0;
  uint64_t here = 4 * m, high = 4 * m + 2;
  uint64_t low = (m == (UINT64_C(1) << 52) && e > -1074) ? 4 * m - 1 : 4 * m - 2;
  /* k: 10^k <= value < 10^(k + 1), found exactly. */
  int k = (int)floor(log10(value));
  for (;;) {
    /* 10^k against value: 10^k * 2^max(-s,0) vs here * 2^max(s,0), with 10^-k on value's side for k < 0. */
    big_number left, right;
    big_scaled(&left, 1, s < 0 ? -s : 0, k > 0 ? k : 0);
    big_scaled(&right, here, s > 0 ? s : 0, k < 0 ? -k : 0);
    if (big_compare(&left, &right) > 0) {
      k--;
      continue;
    }
    big_scaled(&left, 1, s < 0 ? -s : 0, k + 1 > 0 ? k + 1 : 0);
    big_scaled(&right, here, s > 0 ? s : 0, k + 1 < 0 ? -(k + 1) : 0);
    if (big_compare(&left, &right) <= 0) {
      k++;
      continue;
    }
    break;
  }
  for (int count = 1; count <= 17; count++) {
    int p = k - count + 1;
    /* x / 10^p = A / B; the interval's ends scaled alike. */
    int twos_a = s > 0 ? s : 0, tens_a = p < 0 ? -p : 0;
    big_number b, a, a_twice, a_low, a_high;
    big_scaled(&b, 1, s < 0 ? -s : 0, p > 0 ? p : 0);
    big_scaled(&a, here, twos_a, tens_a);
    big_scaled(&a_twice, here, twos_a + 1, tens_a);
    big_scaled(&a_low, low, twos_a, tens_a);
    big_scaled(&a_high, high, twos_a, tens_a);
    /* c = floor(A / B), from an estimate within a few units, made exact. */
    double estimate = value;
    int scale = -p;
    while (scale > 200) {
      estimate *= 1e200;
      scale -= 200;
    }
    while (scale < -200) {
      estimate *= 1e-200;
      scale += 200;
    }
    estimate *= pow(10.0, (double)scale);
    uint64_t limit = ten_to(count);
    uint64_t c = estimate <= 0 ? 0 : (estimate >= (double)limit ? limit : (uint64_t)estimate);
    while (c > 0 && compare_scaled(c, &b, &a) > 0) c--;
    while (compare_scaled(c + 1, &b, &a) <= 0) c++;
    /* Rounded half to even: (2c + 1) * B against 2A. */
    int half = compare_scaled(2 * c + 1, &b, &a_twice);
    if (half < 0 || (half == 0 && (c % 2) == 1)) c++;
    uint64_t lowest = ten_to(count - 1);
    uint64_t candidates[2] = {c, 0};
    int n_candidates = 1;
    /* The neighbor on value's other side, when it has `count` digits too
     * (9.99 up one unit is 1.00 of the next power of ten, and 1.00 down one
     * 9.99 of the one below: fewer digits, tried already). */
    int side = compare_scaled(c, &b, &a);
    if (side < 0 && c + 1 < limit) candidates[n_candidates++] = c + 1;
    if (side > 0 && c > lowest) candidates[n_candidates++] = c - 1;
    for (int i = 0; i < n_candidates; i++) {
      int from_low = compare_scaled(candidates[i], &b, &a_low);
      int to_high = compare_scaled(candidates[i], &b, &a_high);
      if ((even ? from_low >= 0 : from_low > 0) && (even ? to_high <= 0 : to_high < 0)) {
        uint64_t q = candidates[i];
        int places = count;
        int top = k + 1;
        if (q == limit) {
          /* Rounded up to the next power of ten: one digit, a place higher. */
          q = 1;
          places = 1;
          top = k + 2;
        }
        char text[24];
        int length = write_digits(q, places, text);
        while (length > 1 && text[length - 1] == '0') length--;
        memcpy(digits, text, (size_t)length);
        digits[length] = '\0';
        *point = top;
        return;
      }
    }
  }
  /* 17 significant digits always name a double: never reached. */
  error("no decimal of 17 digits names the number");
}

/* ECMAScript Number::toString of a finite number, into `out` (32 bytes). */
static void number_text(double value, char *out) {
  if (value == 0) {
    strcpy(out, "0"); /* -0 too, as JavaScript writes it */
    return;
  }
  char *at = out;
  if (value < 0) {
    *at++ = '-';
    value = -value;
  }
  char digits[24];
  int n;
  shortest_digits(value, digits, &n);
  int k = (int)strlen(digits);
  if (k <= n && n <= 21) {
    memcpy(at, digits, (size_t)k);
    at += k;
    for (int i = 0; i < n - k; i++) *at++ = '0';
    *at = '\0';
  } else if (0 < n && n <= 21) {
    memcpy(at, digits, (size_t)n);
    at += n;
    *at++ = '.';
    memcpy(at, digits + n, (size_t)(k - n));
    at += k - n;
    *at = '\0';
  } else if (-6 < n && n <= 0) {
    *at++ = '0';
    *at++ = '.';
    for (int i = 0; i < -n; i++) *at++ = '0';
    memcpy(at, digits, (size_t)k);
    at += k;
    *at = '\0';
  } else {
    *at++ = digits[0];
    if (k > 1) {
      *at++ = '.';
      memcpy(at, digits + 1, (size_t)(k - 1));
      at += k - 1;
    }
    /* A double's decimal exponent has at most three digits. */
    int exponent = n - 1 >= 0 ? n - 1 : 1 - n;
    *at++ = 'e';
    *at++ = n - 1 >= 0 ? '+' : '-';
    if (exponent >= 100) *at++ = (char)('0' + exponent / 100);
    if (exponent >= 10) *at++ = (char)('0' + exponent / 10 % 10);
    *at++ = (char)('0' + exponent % 10);
    *at = '\0';
  }
}

/* Each number of a double or integer vector as JavaScript writes it; every
 * one must be finite (JSON holds no other). */
SEXP cellucid_number_text(SEXP x) {
  R_xlen_t n = XLENGTH(x);
  SEXP out = PROTECT(allocVector(STRSXP, n));
  char text[32];
  if (TYPEOF(x) == INTSXP) {
    const int *values = INTEGER(x);
    for (R_xlen_t i = 0; i < n; i++) {
      if (values[i] == NA_INTEGER) error("JSON holds no NA");
      int64_t value = values[i];
      char *at = text;
      if (value < 0) *at++ = '-';
      write_digits((uint64_t)(value < 0 ? -value : value), 1, at);
      SET_STRING_ELT(out, i, mkChar(text));
    }
  } else if (TYPEOF(x) == REALSXP) {
    const double *values = REAL(x);
    for (R_xlen_t i = 0; i < n; i++) {
      if (!R_FINITE(values[i])) error("JSON holds finite numbers only");
      number_text(values[i], text);
      SET_STRING_ELT(out, i, mkChar(text));
    }
  } else {
    error("number_text: numbers must be a double or integer vector");
  }
  UNPROTECT(1);
  return out;
}
