test_that("polar and Cartesian conversions round trip", {
  d <- c(1, 5, 12.3); a <- c(0, 90, 237)
  xy <- polar_to_xy(d, a)
  back <- xy_to_polar(xy$x, xy$y)
  expect_equal(back$distance, d, tolerance = 1e-8)
  expect_equal(back$azimuth, a, tolerance = 1e-6)
})

test_that("north is +y and east is +x", {
  expect_equal(polar_to_xy(10, 0), data.frame(x = 0, y = 10), tolerance = 1e-9)
  expect_equal(polar_to_xy(10, 90), data.frame(x = 10, y = 0), tolerance = 1e-9)
})

test_that("polygon area is exact for a square and a regular polygon", {
  sq <- cbind(c(0, 2, 2, 0), c(0, 0, 2, 2))
  expect_equal(forestCI:::polygon_area(sq), 4)
  circ <- forestCI:::regular_polygon(0, 0, 1, 3600)
  expect_equal(forestCI:::polygon_area(circ), pi, tolerance = 1e-5)
})

test_that("half plane clipping halves a square", {
  sq <- cbind(c(-1, 1, 1, -1), c(-1, -1, 1, 1))
  half <- forestCI:::clip_halfplane(sq, 1, 0, 0)
  expect_equal(forestCI:::polygon_area(half), 2, tolerance = 1e-9)
})
