# Cellucid's own category colors, as the web app gives them (R/category-colors.R).
# A user field a session applies without colors of its own is written with
# the colors the viewer showed it in: Cellucid's category colors by position
# (data/palettes.js getCategoryColor). The shared fixture, written by the web
# app's scripts/generate-category-default-colors.mjs, pins the first 1,024
# byte for byte in cellucid, cellucid-python and cellucid-r.

.category_default_colors_sha256 <- "66e652c560afa8aa5bbe149ed87379fd9b937e5a3ca0585d669fdc09092576e0"

test_that("the colors are the web app's", {
  expect_pinned_fixture("category_default_colors.json", .category_default_colors_sha256)
  expected <- cellucid_read_fixture("category_default_colors.json")$colors
  expect_identical(cellucid:::.default_category_colors(length(expected)), unlist(expected))
  expect_identical(cellucid:::.default_category_colors(0), character())
})

test_that("a count is a non-negative whole number", {
  for (count in list(-1, 1.5, NA, "2", c(1, 2))) {
    expect_error(cellucid:::.default_category_colors(count), "non-negative whole number")
  }
})
