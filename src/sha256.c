/*
 * SHA-256 (FIPS 180-4) of one raw vector, as lower-case hex: the digest the
 * Cellucid web build's inventory (cellucid-web-assets.json) declares for every
 * file, which the viewer's web cache checks before it serves a byte
 * (R/web-cache.R). Base R computes SHA-256 only from R 4.5
 * (tools::sha256sum, and only of files), and this package supports R 4.1.
 */

#include <R.h>
#include <Rinternals.h>

#include <stdint.h>
#include <string.h>

static const uint32_t cellucid_sha256_k[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
    0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
    0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
    0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
    0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
    0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};

#define CELLUCID_ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))

static void cellucid_sha256_block(uint32_t state[8], const unsigned char block[64]) {
  uint32_t w[64];
  for (int i = 0; i < 16; i++) {
    w[i] = ((uint32_t)block[4 * i] << 24) | ((uint32_t)block[4 * i + 1] << 16) |
           ((uint32_t)block[4 * i + 2] << 8) | (uint32_t)block[4 * i + 3];
  }
  for (int i = 16; i < 64; i++) {
    uint32_t s0 = CELLUCID_ROTR(w[i - 15], 7) ^ CELLUCID_ROTR(w[i - 15], 18) ^ (w[i - 15] >> 3);
    uint32_t s1 = CELLUCID_ROTR(w[i - 2], 17) ^ CELLUCID_ROTR(w[i - 2], 19) ^ (w[i - 2] >> 10);
    w[i] = w[i - 16] + s0 + w[i - 7] + s1;
  }
  uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
  uint32_t e = state[4], f = state[5], g = state[6], h = state[7];
  for (int i = 0; i < 64; i++) {
    uint32_t s1 = CELLUCID_ROTR(e, 6) ^ CELLUCID_ROTR(e, 11) ^ CELLUCID_ROTR(e, 25);
    uint32_t ch = (e & f) ^ (~e & g);
    uint32_t t1 = h + s1 + ch + cellucid_sha256_k[i] + w[i];
    uint32_t s0 = CELLUCID_ROTR(a, 2) ^ CELLUCID_ROTR(a, 13) ^ CELLUCID_ROTR(a, 22);
    uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
    uint32_t t2 = s0 + maj;
    h = g;
    g = f;
    f = e;
    e = d + t1;
    d = c;
    c = b;
    b = a;
    a = t1 + t2;
  }
  state[0] += a;
  state[1] += b;
  state[2] += c;
  state[3] += d;
  state[4] += e;
  state[5] += f;
  state[6] += g;
  state[7] += h;
}

SEXP cellucid_sha256(SEXP bytes) {
  if (TYPEOF(bytes) != RAWSXP) {
    Rf_error("A SHA-256 is computed over a raw vector.");
  }
  uint32_t state[8] = {0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
                       0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
  const unsigned char *data = RAW(bytes);
  uint64_t count = (uint64_t)XLENGTH(bytes);
  uint64_t whole = count / 64;
  for (uint64_t i = 0; i < whole; i++) {
    cellucid_sha256_block(state, data + 64 * i);
  }
  unsigned char tail[128];
  memset(tail, 0, sizeof(tail));
  uint64_t rest = count - 64 * whole;
  if (rest > 0) memcpy(tail, data + 64 * whole, (size_t)rest);
  tail[rest] = 0x80;
  int blocks = rest < 56 ? 1 : 2;
  uint64_t bits = count * 8;
  for (int i = 0; i < 8; i++) {
    tail[64 * blocks - 1 - i] = (unsigned char)(bits >> (8 * i));
  }
  for (int i = 0; i < blocks; i++) {
    cellucid_sha256_block(state, tail + 64 * i);
  }
  static const char digits[] = "0123456789abcdef";
  char hex[65];
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 4; j++) {
      unsigned char byte = (unsigned char)(state[i] >> (24 - 8 * j));
      hex[8 * i + 2 * j] = digits[byte >> 4];
      hex[8 * i + 2 * j + 1] = digits[byte & 0x0f];
    }
  }
  hex[64] = '\0';
  return Rf_mkString(hex);
}
