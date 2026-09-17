test_that("geometric crown dimensions match closed form solids", {
  # cone of base radius R and height L: lateral area = pi*R*sqrt(R^2+L^2),
  # volume = pi*R^2*L/3. Put the widest point at the crown base.
  p <- crown_profile("geometric", solid = "cone")
  R <- 2; L <- 8
  d <- crown_dimensions(p, lcw = 2 * R, height = 20, hcb = 20 - L,
                        widest = 1, n_sub = 4000)
  expect_equal(d$csa, pi * R * sqrt(R^2 + L^2), tolerance = 1e-3)
  expect_equal(d$cv, pi * R^2 * L / 3, tolerance = 1e-3)
  expect_equal(d$cpa, pi * R^2, tolerance = 1e-9)
})

test_that("paraboloid volume is half the enclosing cylinder", {
  p <- crown_profile("geometric", solid = "paraboloid")
  R <- 3; L <- 10
  d <- crown_dimensions(p, lcw = 2 * R, height = 20, hcb = 10,
                        widest = 1, n_sub = 4000)
  expect_equal(d$cv, pi * R^2 * L / 2, tolerance = 1e-3)
})

test_that("crown radius is zero at the apex and largest at the widest point", {
  p <- crown_profile("dual_exponent")
  h <- 20; hcb <- 10; lcw <- 4; w <- 0.6
  expect_equal(crown_radius_at(p, h, lcw, h, hcb, widest = w), 0,
               tolerance = 1e-9)
  zw <- h - w * (h - hcb)
  expect_equal(crown_radius_at(p, zw, lcw, h, hcb, widest = w), lcw / 2,
               tolerance = 1e-6)
  expect_equal(crown_radius_at(p, hcb - 1, lcw, h, hcb, widest = w), 0)
})

test_that("every profile family runs and gives positive finite dimensions", {
  for (fam in c("dual_exponent", "variable_exponent", "geometric")) {
    p <- crown_profile(fam)
    d <- crown_dimensions(p, lcw = 4, height = 20, hcb = 10, species = "RS")
    expect_true(all(is.finite(unlist(d))), info = fam)
    expect_true(all(unlist(d) > 0), info = fam)
    expect_lt(d$csa_upper + d$csa_lower - d$csa, 1e-9)
  }
})

test_that("largest crown width never exceeds maximum crown width", {
  dbh <- seq(5, 80, by = 5)
  for (sp in c("RS", "BF", "SM", "RO", "WP", "ZZ")) {
    m <- max_crown_width(dbh, sp)
    l <- largest_crown_width(dbh, sp, mcw = m)
    expect_true(all(l <= m + 1e-9), info = sp)
  }
})

test_that("adaptive integration beats Simpson on a sharply pointed crown", {
  # a gothic dome with a high exponent turns vertical at the apex, which is the
  # case a fixed rule under-integrates
  p <- crown_profile("dual_exponent")
  a <- crown_dimensions(p, lcw = 4, height = 20, hcb = 10, widest = 0.8,
                        up_exp = 5, lo_exp = 5, shape = "g",
                        method = "adaptive")
  s1 <- crown_dimensions(p, lcw = 4, height = 20, hcb = 10, widest = 0.8,
                         up_exp = 5, lo_exp = 5, shape = "g",
                         method = "simpson", n_sub = 200L)
  s2 <- crown_dimensions(p, lcw = 4, height = 20, hcb = 10, widest = 0.8,
                         up_exp = 5, lo_exp = 5, shape = "g",
                         method = "simpson", n_sub = 4000L)
  expect_true(all(is.finite(unlist(a))))
  # Simpson converges up toward the adaptive value as it is refined
  expect_lt(abs(s2$csa - a$csa), abs(s1$csa - a$csa))
  expect_equal(s2$cv, a$cv, tolerance = 1e-4)
})

test_that("adaptive integration still reproduces the closed form solids", {
  p <- crown_profile("geometric", solid = "cone")
  R <- 2; L <- 8
  d <- crown_dimensions(p, lcw = 2 * R, height = 20, hcb = 20 - L,
                        widest = 1, method = "adaptive")
  expect_equal(d$csa, pi * R * sqrt(R^2 + L^2), tolerance = 1e-5)
  expect_equal(d$cv, pi * R^2 * L / 3, tolerance = 1e-5)
})
