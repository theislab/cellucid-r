# Every compressed payload carries the same ten header bytes on every machine.
#
# `gzfile()` leaves the gzip member header to the zlib R was built against,
# which stamps its own build platform into the `OS` field: `0x03` on Unix,
# `0x0b` on Windows. The same input would then produce different export bytes
# depending on the operating system that ran the export, and would never agree
# with `cellucid.prepare()`, which pins the header explicitly. These tests
# assert the ten bytes directly, and assert that the member still decompresses
# to exactly the payload that went in -- a canonical header on a corrupted
# stream would be worse than a platform-dependent one on a sound stream.

.gzip_header_bytes <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  readBin(connection, what = "raw", n = 10L)
}

# RFC 1952 section 2.3.1. `XFL` records which end of the speed/ratio range the
# caller asked for; `cellucid.prepare()` writes the same three values.
.gzip_expected_header <- function(compression) {
  extra_flags <- if (compression >= 9L) {
    0x02
  } else if (compression <= 1L) {
    0x04
  } else {
    0x00
  }
  as.raw(c(
    0x1f, 0x8b, # magic
    0x08,       # deflate
    0x00,       # no FTEXT, FHCRC, FEXTRA, FNAME or FCOMMENT
    0x00, 0x00, 0x00, 0x00, # MTIME, fixed at the Unix epoch
    extra_flags,
    0xff        # OS, the "unknown" code, identical on every platform
  ))
}

.gzip_read_raw <- function(path, n) {
  connection <- gzfile(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  readBin(connection, what = "raw", n = n)
}

# `expect_identical()` reports a difference by building an element-by-element
# diff, which does not finish on two raw vectors of the size these payloads
# reach when most of their positions differ -- exactly the shape a byte-order
# or framing regression produces. Comparing the length, then the first
# differing position, then a bounded window around it keeps a real failure
# readable and bounded instead of hanging the run that found it.
.gzip_expect_payload <- function(observed, expected) {
  expect_identical(length(observed), length(expected))
  first_difference <- which(observed[seq_along(expected)] != expected)[1L]
  expect_identical(first_difference, NA_integer_)
  if (!is.na(first_difference)) {
    window <- seq.int(
      first_difference,
      min(length(expected), first_difference + 7L)
    )
    expect_identical(observed[window], expected[window])
  }
  invisible(NULL)
}

.gzip_fixture <- function(n_cells = 96L) {
  set.seed(20260731L)
  list(
    n_cells = n_cells,
    latent = matrix(stats::rnorm(n_cells * 3L), ncol = 3L),
    umap2 = matrix(stats::rnorm(n_cells * 2L), ncol = 2L),
    obs = data.frame(
      cluster = factor(rep(c("A", "B", "C", "D"), length.out = n_cells)),
      score = stats::runif(n_cells)
    )
  )
}

test_that("every compressed payload of an export carries the canonical header", {
  fixture <- .gzip_fixture()
  out <- file.path(tempdir(), "cellucid_gzip_header_export")
  unlink(out, recursive = TRUE, force = TRUE)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)

  cellucid_prepare(
    dataset_id = "gzip-header",
    dataset_name = "Gzip header",
    latent_space = fixture$latent,
    obs = fixture$obs,
    X_umap_2d = fixture$umap2,
    out_dir = out,
    centroid_min_points = 4L,
    force = TRUE,
    obs_categorical_dtype = "uint16",
    compression = 6L
  )

  # The exact set, not its size. A payload family that stopped being
  # compressed drops out of a `.gz` listing entirely, so a count-based guard
  # leaves the loop below asserting nothing whatever about it. points_%dd.bin
  # is the case that matters: it sits at the export root, where the
  # manifest-versus-disk check does not reach either, so nothing else in the
  # suite would notice it losing its .gz while dataset_identity.json went on
  # naming one.
  written <- list.files(out, pattern = "\\.gz$", recursive = TRUE)
  expect_identical(
    written,
    c(
      "obs/0.codes.u16.gz",
      "obs/0.outliers.f32.gz",
      "obs/1.values.f32.gz",
      "points_2d.bin.gz"
    )
  )
  expected <- .gzip_expected_header(6L)
  for (relative in written) {
    expect_identical(
      .gzip_header_bytes(file.path(out, relative)),
      expected,
      info = relative
    )
  }
})

test_that("the canonical header is written at every compression level", {
  values <- as.numeric(seq_len(4096L)) / 4096
  for (compression in 1:9) {
    path <- tempfile("cellucid-gzip-header-level-")
    on.exit(unlink(c(path, paste0(path, ".gz")), force = TRUE), add = TRUE)
    written <- cellucid:::.write_float32_vector(
      path,
      values,
      compression = compression
    )
    expect_identical(
      .gzip_header_bytes(written),
      .gzip_expected_header(compression)
    )
  }
})

test_that("a canonical header still decompresses to exactly the payload", {
  # One case per binary writer, so no writer can acquire the header without
  # the payload surviving the round trip. Every expectation is written out
  # rather than produced by the writer's own packer: an expectation that calls
  # `.pack_uint16()` to describe what `.write_uint16()` should have produced
  # moves with the packer and agrees with a reversed one. The exact byte layout
  # of each width is pinned in test-byte-order.R; these literals repeat only as
  # much of it as this test needs to be evidence.
  cases <- list(
    list(
      write = function(path, compression) {
        cellucid:::.write_float32_vector(
          path,
          as.numeric(seq_len(2048L)) / 64,
          compression = compression
        )
      },
      expected = writeBin(
        as.numeric(seq_len(2048L)) / 64,
        raw(),
        size = 4L,
        endian = "little"
      )
    ),
    list(
      write = function(path, compression) {
        cellucid:::.write_float64_vector(
          path,
          c(0.5, 1 + (2^-30), 2),
          compression = compression
        )
      },
      expected = writeBin(
        c(0.5, 1 + (2^-30), 2),
        raw(),
        size = 8L,
        endian = "little"
      )
    ),
    list(
      write = function(path, compression) {
        cellucid:::.write_uint8(path, 0:255, compression = compression)
      },
      expected = as.raw(0:255)
    ),
    list(
      write = function(path, compression) {
        cellucid:::.write_uint16(
          path,
          c(0L, 258L, 65535L),
          compression = compression
        )
      },
      # 258 is 0x0102, so its two bytes are not a palindrome and a reversed
      # packer would decompress to 0x0201 here.
      expected = as.raw(c(0x00, 0x00, 0x02, 0x01, 0xff, 0xff))
    ),
    list(
      write = function(path, compression) {
        cellucid:::.write_uint32(
          path,
          c(0, 16909060, 4294967295),
          compression = compression
        )
      },
      # 16909060 is 0x01020304: four distinct bytes.
      expected = as.raw(c(
        0x00, 0x00, 0x00, 0x00,
        0x04, 0x03, 0x02, 0x01,
        0xff, 0xff, 0xff, 0xff
      ))
    )
  )

  for (case in cases) {
    path <- tempfile("cellucid-gzip-roundtrip-")
    on.exit(unlink(c(path, paste0(path, ".gz")), force = TRUE), add = TRUE)
    written <- case$write(path, 6L)
    expect_identical(.gzip_header_bytes(written), .gzip_expected_header(6L))
    .gzip_expect_payload(
      .gzip_read_raw(written, length(case$expected) + 1L),
      case$expected
    )
  }
})

test_that("a payload larger than one deflate buffer keeps header and payload", {
  # The header is written once at the head of the file while the deflate
  # stream is still being produced, so a payload that spans many internal
  # buffers is the case where a rewritten header could land on the wrong
  # bytes or be overwritten again.
  set.seed(11L)
  values <- stats::rnorm(1024L * 1024L)
  expected <- writeBin(values, raw(), size = 4L, endian = "little")
  path <- tempfile("cellucid-gzip-large-")
  on.exit(unlink(c(path, paste0(path, ".gz")), force = TRUE), add = TRUE)

  written <- cellucid:::.write_float32_vector(path, values, compression = 6L)
  expect_identical(.gzip_header_bytes(written), .gzip_expected_header(6L))
  .gzip_expect_payload(.gzip_read_raw(written, length(expected) + 1L), expected)
})

test_that("an empty payload is a complete canonical member", {
  # `connectivities` accepts an exact empty graph, so an empty compressed
  # payload is a case the format reaches.
  path <- tempfile("cellucid-gzip-empty-")
  on.exit(unlink(c(path, paste0(path, ".gz")), force = TRUE), add = TRUE)
  written <- cellucid:::.write_uint32(path, numeric(0), compression = 6L)

  expect_identical(.gzip_header_bytes(written), .gzip_expected_header(6L))
  expect_identical(.gzip_read_raw(written, 1L), raw(0))
})

test_that("a member that is not the fixed ten-byte form stops the export", {
  # The header is replaced in place, so the writer has to prove it is looking
  # at the header it was verified against before it overwrites anything. A
  # member carrying an FNAME field would put the deflate stream elsewhere.
  path <- tempfile("cellucid-gzip-foreign-", fileext = ".gz")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  connection <- file(path, open = "wb")
  writeBin(
    c(
      as.raw(c(0x1f, 0x8b, 0x08, 0x08)), # FNAME set
      as.raw(c(0x00, 0x00, 0x00, 0x00)),
      as.raw(c(0x00, 0xff)),
      charToRaw("x"),
      as.raw(0x00)
    ),
    connection
  )
  close(connection)

  expect_error(
    cellucid:::.write_canonical_gzip_header(path, 6L),
    "does not begin with the fixed ten-byte gzip header"
  )
  # The refused member is left exactly as it was found.
  expect_identical(.gzip_header_bytes(path)[[4L]], as.raw(0x08))
})

test_that("an uncompressed payload is written with no gzip framing at all", {
  path <- tempfile("cellucid-gzip-absent-")
  on.exit(unlink(c(path, paste0(path, ".gz")), force = TRUE), add = TRUE)
  written <- cellucid:::.write_uint8(path, 0:255, compression = NULL)
  expect_identical(written, path)
  expect_false(file.exists(paste0(path, ".gz")))
  expect_identical(.gzip_header_bytes(written), as.raw(0:9))
})
