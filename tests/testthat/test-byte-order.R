# One input produces one sequence of bytes, and that sequence is
# little-endian, on every machine that runs this writer.
#
# The viewer reads every payload as a little-endian typed array, so a payload
# written in the host's own byte order would be read backwards on a big-endian
# machine -- silently, as plausible numbers rather than as an error. Two
# separate mechanisms stop that, and this file pins both:
#
#   * the float widths pass `endian = "little"` to `writeBin()`;
#   * `.pack_uint16()` and `.pack_uint32()` lay their bytes out by arithmetic,
#     least significant first, and hand `writeBin()` a `raw` vector.
#
# Every expectation below is a `raw` literal derived by hand from IEEE 754 and
# from base-2 place value -- never produced by `writeBin()`, never by the
# packers -- and every payload is read back as bytes, with no width and no byte
# order named on the reading side. Both halves carry weight. An expectation the
# code under test computed for itself moves with that code and agrees with a
# broken writer; a readback that names a width and a byte order re-encodes what
# it read, so a writer and a reader that are wrong in the same direction agree
# with each other.
#
# The values are chosen so that no expected byte sequence is a palindrome:
# reversing any literal below produces a different literal. A packer with its
# byte lanes swapped is the defect these catch that a fixture of 0 and 65535
# cannot.
#
# One thing a byte comparison cannot do, here or anywhere: tell `endian =
# "little"` apart from the host's own order while the host is little-endian.
# The structural rules at the end of this file are what hold that.

# Read the payload back as bytes, with no width and no byte order anywhere on
# the reading side.
.byte_order_bytes <- function(path, compression, n) {
  connection <- if (is.null(compression)) {
    file(path, open = "rb")
  } else {
    gzfile(path, open = "rb")
  }
  on.exit(close(connection), add = TRUE)
  # One byte more than expected, so a payload longer than the literal fails
  # rather than matching its own prefix.
  readBin(connection, what = "raw", n = n + 1L)
}

test_that("a float32 payload is the IEEE 754 bytes, least significant first", {
  # value                        IEEE 754 binary32   little-endian bytes
  # 1                            0x3F800000          00 00 80 3F
  # -2                           0xC0000000          00 00 00 C0
  # 0.5                          0x3F000000          00 00 00 3F
  # 16 * (1 + 0x38C4D2 / 2^23)   0x41B8C4D2          D2 C4 B8 41
  #
  # The last value is written as its own fields -- sign 0, exponent 0x83 so
  # 2^4, mantissa 0x38C4D2 over 2^23 -- so the expectation is readable as
  # arithmetic rather than as a magic decimal. All four of its bytes differ,
  # which is what makes a reversed writer visible in every byte position and
  # not only in the two the smaller constants exercise.
  values <- c(1, -2, 0.5, 16 * (1 + 0x38C4D2 / 2^23))
  expected <- as.raw(c(
    0x00, 0x00, 0x80, 0x3f,
    0x00, 0x00, 0x00, 0xc0,
    0x00, 0x00, 0x00, 0x3f,
    0xd2, 0xc4, 0xb8, 0x41
  ))

  for (compression in list(NULL, 6L)) {
    path <- tempfile("cellucid-byte-order-f32-")
    on.exit(unlink(c(path, paste0(path, ".gz")), force = TRUE), add = TRUE)
    written <- cellucid:::.write_float32_vector(
      path,
      values,
      compression = compression
    )
    expect_identical(
      .byte_order_bytes(written, compression, length(expected)),
      expected
    )
  }
})

test_that("a float64 payload is the IEEE 754 bytes, least significant first", {
  # value                              IEEE 754 binary64      little-endian
  # 1                                  0x3FF0000000000000     00 .. 00 F0 3F
  # -2                                 0xC000000000000000     00 .. 00 00 C0
  # 1 + 2^-52                          0x3FF0000000000001     01 00 .. F0 3F
  # 2048 * (1 + 0x123456789ABCD/2^52)  0x40A123456789ABCD     CD AB .. A1 40
  #
  # `1 + 2^-52` is the smallest double above 1, so its one set mantissa bit is
  # the first byte written when the order is little-endian and the last when it
  # is not. The final value has eight distinct bytes.
  values <- c(1, -2, 1 + 2^-52, 2048 * (1 + 0x123456789ABCD / 2^52))
  expected <- as.raw(c(
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc0,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,
    0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0xa1, 0x40
  ))

  for (compression in list(NULL, 6L)) {
    path <- tempfile("cellucid-byte-order-f64-")
    on.exit(unlink(c(path, paste0(path, ".gz")), force = TRUE), add = TRUE)
    written <- cellucid:::.write_float64_vector(
      path,
      values,
      compression = compression
    )
    expect_identical(
      .byte_order_bytes(written, compression, length(expected)),
      expected
    )
  }
})

test_that("the unsigned packers put the least significant byte first", {
  # 258 is 0x0102, so its bytes are 0x02 then 0x01 and never the other way
  # round; 65280 is 0xFF00, so a reversed packer would publish 65535's
  # neighbour instead. 0 and 65535 are byte palindromes and cannot tell the two
  # apart on their own, so they are here for the range boundary only.
  expect_identical(
    cellucid:::.pack_uint16(c(0, 1, 258, 65280, 65535)),
    as.raw(c(
      0x00, 0x00,
      0x01, 0x00,
      0x02, 0x01,
      0x00, 0xff,
      0xff, 0xff
    ))
  )

  # 16909060 is 0x01020304: four distinct bytes, so every one of the four
  # positions is pinned.
  expect_identical(
    cellucid:::.pack_uint32(c(0, 1, 16909060, 4278190080, 4294967295)),
    as.raw(c(
      0x00, 0x00, 0x00, 0x00,
      0x01, 0x00, 0x00, 0x00,
      0x04, 0x03, 0x02, 0x01,
      0x00, 0x00, 0x00, 0xff,
      0xff, 0xff, 0xff, 0xff
    ))
  )

  # The packers take integer and double input by the same rule, and the two
  # spellings of one value must produce one payload.
  expect_identical(
    cellucid:::.pack_uint16(c(258L, 65280L)),
    cellucid:::.pack_uint16(c(258, 65280))
  )
  expect_identical(
    cellucid:::.pack_uint32(c(1L, 16909060L)),
    cellucid:::.pack_uint32(c(1, 16909060))
  )
})

test_that("an unsigned payload reaches the export as the packed bytes", {
  # The packers own the byte order; these writers own passing it through
  # untouched, plain and compressed.
  cases <- list(
    list(
      write = cellucid:::.write_uint8,
      values = c(0, 1, 127, 128, 254, 255),
      expected = as.raw(c(0x00, 0x01, 0x7f, 0x80, 0xfe, 0xff))
    ),
    list(
      write = cellucid:::.write_uint16,
      values = c(0, 1, 258, 65280, 65535),
      expected = as.raw(c(
        0x00, 0x00,
        0x01, 0x00,
        0x02, 0x01,
        0x00, 0xff,
        0xff, 0xff
      ))
    ),
    list(
      write = cellucid:::.write_uint32,
      values = c(0, 1, 16909060, 4278190080, 4294967295),
      expected = as.raw(c(
        0x00, 0x00, 0x00, 0x00,
        0x01, 0x00, 0x00, 0x00,
        0x04, 0x03, 0x02, 0x01,
        0x00, 0x00, 0x00, 0xff,
        0xff, 0xff, 0xff, 0xff
      ))
    )
  )

  for (case in cases) {
    for (compression in list(NULL, 6L)) {
      path <- tempfile("cellucid-byte-order-uint-")
      on.exit(unlink(c(path, paste0(path, ".gz")), force = TRUE), add = TRUE)
      written <- case$write(path, case$values, compression = compression)
      expect_identical(
        .byte_order_bytes(written, compression, length(case$expected)),
        case$expected
      )
    }
  }
})

test_that("an exported embedding carries the same bytes as the writer", {
  # The literals above pin the writers; this pins the path from a real
  # `cellucid_prepare()` call to those writers, so a payload that reached the
  # export some other way would be visible here.
  coords <- matrix(
    c(
      -1, -1,
      1, 1,
      0.5, -0.5,
      -0.5, 0.5
    ),
    ncol = 2L,
    byrow = TRUE
  )
  out <- file.path(tempdir(), "cellucid_byte_order_export")
  unlink(out, recursive = TRUE, force = TRUE)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)

  cellucid_prepare(
    dataset_id = "byte-order",
    dataset_name = "Byte order",
    latent_space = coords,
    obs = data.frame(cluster = factor(c("A", "B", "A", "B"))),
    X_umap_2d = coords,
    out_dir = out,
    centroid_min_points = 2L,
    force = TRUE,
    obs_categorical_dtype = "uint16"
  )

  # The coordinates are already centred on the origin and already span
  # [-1, 1] on the widest axis, so the published values are the input values:
  # centre 0, scale 2 / 2 = 1. Row major, x then y.
  expect_identical(
    .byte_order_bytes(file.path(out, "points_2d.bin"), NULL, 32L),
    as.raw(c(
      0x00, 0x00, 0x80, 0xbf, 0x00, 0x00, 0x80, 0xbf, # -1, -1
      0x00, 0x00, 0x80, 0x3f, 0x00, 0x00, 0x80, 0x3f, #  1,  1
      0x00, 0x00, 0x00, 0x3f, 0x00, 0x00, 0x00, 0xbf, #  0.5, -0.5
      0x00, 0x00, 0x00, 0xbf, 0x00, 0x00, 0x00, 0x3f  # -0.5,  0.5
    ))
  )

  # Two categories alternating, stored as uint16: codes 0, 1, 0, 1.
  expect_identical(
    .byte_order_bytes(file.path(out, "obs", "0.codes.u16"), NULL, 8L),
    as.raw(c(0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00))
  )
})

# One mutation is invisible to every assertion above, and to any assertion that
# reaches its verdict by looking at bytes: delete `endian = "little"` from a
# float writer and `writeBin()` falls back to `.Platform$endian`, which is
# "little" on this machine and on all three operating systems the workflow
# matrix runs. The published bytes do not move, so nothing observes the loss.
# It would surface only on a big-endian host -- the one host where the argument
# is load bearing, and the one host nobody here runs. Measured rather than
# argued: with `endian = "little"` removed from both float writers, every byte
# comparison in this suite still passed.
#
# So the argument is checked where it is written rather than by what it does.
# The rules below read the package's own code as code -- every call it makes,
# from the parsed bodies in its namespace -- which is what lets them fail on
# the mutation a byte comparison cannot see. They are the peer of the
# structural test in cellucid-python that refuses a `tobytes`/`tofile` outside
# its byte-order rule.

# Read from the installed namespace rather than from R/*.R. `body()` keeps the
# parsed call with its arguments whether or not source references were kept, so
# the rule is read off the code that will actually run, and it is read the same
# way here, under R CMD check, and against an installed binary package. Reading
# the source files instead would leave the rule unexamined wherever they are
# absent -- which is every R CMD check, and so every run this project's CI
# performs.
.byte_order_package_functions <- function() {
  namespace <- asNamespace("cellucid")
  objects <- mget(
    ls(namespace, all.names = TRUE),
    envir = namespace,
    ifnotfound = list(NULL)
  )
  objects[vapply(objects, is.function, logical(1))]
}

# Every call to `name` in one parse tree, found by walking the tree rather than
# by matching text, so a call split across lines, nested in another call, or
# written inside a callback is still seen.
.byte_order_find_calls <- function(node, name) {
  found <- list()
  walk <- function(current) {
    if (is.call(current)) {
      target <- current[[1L]]
      if (is.name(target) && identical(as.character(target), name)) {
        found[[length(found) + 1L]] <<- current
      }
    }
    if (is.recursive(current)) {
      for (index in seq_along(current)) {
        # A formal without a default, and a missing argument as in `x[, 1]`,
        # parse to the empty symbol. Binding it to a name and then reading that
        # name is what raises "argument is missing", so it is recognised in
        # place and stepped over instead.
        if (identical(current[[index]], quote(expr = ))) {
          next
        }
        walk(current[[index]])
      }
    }
  }
  walk(node)
  found
}

test_that("every sized writeBin and readBin in the package names little-endian", {
  functions <- .byte_order_package_functions()
  sized <- character()
  offenders <- character()
  for (name in names(functions)) {
    for (target in c("writeBin", "readBin")) {
      definition <- get(target, envir = baseenv())
      for (call_node in .byte_order_find_calls(body(functions[[name]]), target)) {
        arguments <- as.list(match.call(definition, call_node))[-1L]
        if (!("size" %in% names(arguments))) {
          # One byte per element: there is no order to name, and naming one
          # would claim a guarantee the call does not provide.
          next
        }
        sized <- c(sized, name)
        if (!identical(arguments$endian, "little")) {
          offenders <- c(
            offenders,
            sprintf("%s(): %s", name, deparse1(call_node))
          )
        }
      }
    }
  }
  expect_identical(offenders, character())

  # The verdict above is only evidence if the walk reached the calls it is
  # meant to judge, so the walk is measured before it is believed: a walker
  # that found nothing would report no offenders. Naming the functions rather
  # than counting them also makes this the register of every place in the
  # package that reads or writes a multi-byte number -- a new one has to be
  # added here, deliberately, with its byte order already correct.
  expect_setequal(
    sized,
    c(
      ".write_float32_vector",
      ".write_float64_vector",
      ".roundtrip_finite_float32_chunk"
    )
  )
  expect_length(
    .byte_order_find_calls(
      quote(g(function(x) writeBin(x, con, size = 4L))),
      "writeBin"
    ),
    1L
  )
})

test_that("binary payloads are written in one place, through one function", {
  functions <- .byte_order_package_functions()
  writing <- character()
  detached <- character()
  for (name in names(functions)) {
    definition <- body(functions[[name]])
    if (length(.byte_order_find_calls(definition, "writeBin")) == 0L) {
      next
    }
    writing <- c(writing, name)
    if (
      length(.byte_order_find_calls(definition, ".write_binary_payload")) == 0L
    ) {
      detached <- c(detached, name)
    }
  }

  # Five payload families, one write path. A sixth added anywhere in the
  # package fails here, which is the point: a second write path is how a
  # payload would acquire a second byte order, or reach the export without the
  # canonical gzip header.
  expect_setequal(
    setdiff(writing, detached),
    c(
      ".write_float32_vector",
      ".write_float64_vector",
      ".write_uint8",
      ".write_uint16",
      ".write_uint32"
    )
  )

  # The two functions that legitimately call writeBin() without going through
  # it, named so that a third cannot appear unnoticed:
  # .write_canonical_gzip_header() rewrites the ten header bytes of a member
  # .write_binary_payload() has already closed and verified, and
  # .roundtrip_finite_float32_chunk() writes into a raw() vector in memory to
  # round a value into the viewer's float32 domain, publishing nothing.
  expect_setequal(
    detached,
    c(".write_canonical_gzip_header", ".roundtrip_finite_float32_chunk")
  )
})
