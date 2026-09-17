st <- pef_stand()

test_that("the example stand builds and is spatial with complete crowns", {
  expect_s3_class(st, "stand")
  expect_true(is_spatial(st))
  expect_true(isTRUE(attr(st, "has_crown")))
  expect_equal(nrow(st$trees), 74L)
  expect_equal(nrow(st$plots), 2L)
})

test_that("basal area of larger trees is monotone and bounded by plot BA", {
  di <- ci_distance_independent(st)
  d <- merge(di, st$trees[, c("plot_id", "tree_id", "dbh")],
             by = c("plot_id", "tree_id"))
  for (p in unique(d$plot_id)) {
    s <- d[d$plot_id == p, ]
    s <- s[order(-s$dbh), ]
    expect_true(all(diff(s$bal) >= -1e-8), info = p)
    expect_true(all(s$bal <= s$ba + 1e-8), info = p)
    expect_equal(s$bal[1], 0, tolerance = 1e-8)   # largest tree has no larger
  }
  expect_true(all(abs(d$bal - d$bal_sw - d$bal_hw) < 1e-8))
  expect_true(all(abs(d$ccfl - d$ccfl_sw - d$ccfl_hw) < 1e-8))
})

test_that("BAL is symmetric under tied diameters", {
  tst <- data.frame(plot = "a", tree = as.character(1:3), species = "RS",
                    dbh = c(20, 20, 30), expf = 10)
  s <- as_stand(tst, plot = "plot", tree = "tree", species = "species",
                dbh = "dbh", expf = "expf", quiet = TRUE)
  di <- ci_distance_independent(s, sdi_max = NULL)
  expect_equal(di$bal[di$tree_id == "1"], di$bal[di$tree_id == "2"])
})

test_that("distance-dependent indices behave as competition should", {
  nb <- neighbors(st, method = "radius", radius = 6)
  dd <- ci_distance_dependent(st, nb = nb)
  expect_true(all(dd$hegyi >= 0, na.rm = TRUE))
  expect_true(all(dd$n_comp >= 0, na.rm = TRUE))
  # a one-sided filter can only remove competitors
  nb1 <- neighbors(st, method = "radius", radius = 6, filter = "larger_dbh")
  expect_lte(nrow(nb1), nrow(nb))
  dd1 <- ci_distance_dependent(st, nb = nb1, indices = "hegyi")
  expect_true(all(dd1$hegyi <= dd$hegyi + 1e-8, na.rm = TRUE))
  # a wider radius can only add competitors
  dd2 <- ci_distance_dependent(st, radius = 10, indices = "hegyi")
  expect_true(all(dd2$hegyi >= dd$hegyi - 1e-8, na.rm = TRUE))
})

test_that("Hegyi for a known two tree configuration is exact", {
  tst <- data.frame(plot = "a", tree = c("1", "2"), species = "RS",
                    dbh = c(20, 40), x = c(0, 3), y = c(0, 4), expf = 10)
  s <- as_stand(tst, plot = "plot", tree = "tree", species = "species",
                dbh = "dbh", x = "x", y = "y", expf = "expf", quiet = TRUE)
  dd <- ci_distance_dependent(s, radius = 10, indices = "hegyi")
  expect_equal(dd$hegyi[dd$tree_id == "1"], (40 / 20) / 5, tolerance = 1e-10)
  expect_equal(dd$hegyi[dd$tree_id == "2"], (20 / 40) / 5, tolerance = 1e-10)
})

test_that("APA sums to plot area and rAPA agrees with APA", {
  a <- apa(st, weight = "none")
  tot <- tapply(a$apa, a$plot_id, sum)
  area <- st$plots$plot_area * 10000
  expect_equal(as.numeric(tot), area, tolerance = 0.01)
  r <- rapa(st, weight = "none", resolution = 0.25)
  m <- merge(a, r, by = c("plot_id", "tree_id"))
  expect_gt(stats::cor(m$apa, m$rapa), 0.95)
})

test_that("weighting APA moves growing space from small trees to large ones", {
  a0 <- apa(st, weight = "none")
  a1 <- apa(st, weight = "dbh", exponent = 2)
  m <- merge(merge(a0, a1, by = c("plot_id", "tree_id"),
                   suffixes = c("_none", "_dbh")),
             st$trees[, c("plot_id", "tree_id", "dbh")],
             by = c("plot_id", "tree_id"))
  delta <- m$apa_dbh - m$apa_none
  expect_gt(stats::cor(delta, m$dbh), 0.3)
})

test_that("crown exposure is bounded by crown surface area", {
  ce <- crown_exposure(st, n_depth = 8L, n_azimuth = 24L, osv = TRUE,
                       ray_steps = 10L)
  expect_true(all(ce$csax <= ce$csa + 1e-6))
  expect_true(all(ce$osv <= ce$csax + 1e-6))
  expect_true(all(ce$csax_rel >= 0 & ce$csax_rel <= 1 + 1e-9))
  expect_true(all(ce$cpax <= ce$cpa + 1e-6))
})

test_that("an isolated tree is fully exposed", {
  tst <- data.frame(plot = "a", tree = c("1", "2"), species = "RS",
                    dbh = c(20, 20), x = c(0, 60), y = c(0, 0),
                    height = 15, hcb = 8, expf = 10)
  s <- as_stand(tst, plot = "plot", tree = "tree", species = "species",
                dbh = "dbh", height = "height", hcb = "hcb",
                x = "x", y = "y", expf = "expf", quiet = TRUE)
  ce <- crown_exposure(s, n_depth = 8L, n_azimuth = 24L, osv = FALSE)
  expect_equal(ce$csax_rel, c(1, 1), tolerance = 0.02)
})

test_that("Clark-Evans is near 1 for a random pattern", {
  set.seed(42)
  n <- 400; r <- 30
  th <- stats::runif(n, 0, 2 * pi); rad <- r * sqrt(stats::runif(n))
  tst <- data.frame(plot = "a", tree = as.character(seq_len(n)),
                    species = "RS", dbh = 20,
                    x = rad * cos(th), y = rad * sin(th))
  s <- as_stand(tst, plot = "plot", tree = "tree", species = "species",
                dbh = "dbh", x = "x", y = "y", plot_radius = r, quiet = TRUE)
  sp <- spatial_pattern(s)
  expect_gt(sp$clark_evans, 0.9)
  expect_lt(sp$clark_evans, 1.15)
})

test_that("the wrapper runs, and degrades gracefully without coordinates", {
  ci <- competition_indices(st, indices = c("independent", "dependent", "apa"),
                            edge = "flag", quiet = TRUE)
  expect_equal(nrow(ci), 74L)
  expect_true(all(c("bal", "ccfl", "hegyi", "apa", "is_edge") %in% names(ci)))

  flat <- as.data.frame(st$trees)
  s2 <- as_stand(flat, plot = "plot_id", tree = "tree_id", species = "species",
                 dbh = "dbh", expf = "expf", quiet = TRUE)
  expect_false(is_spatial(s2))
  ci2 <- competition_indices(s2, quiet = TRUE)
  expect_true("bal" %in% names(ci2))
  expect_false("hegyi" %in% names(ci2))
})

test_that("edge exclusion blanks distance-dependent values only", {
  ci <- competition_indices(st, indices = c("independent", "dependent"),
                            edge = "exclude", quiet = TRUE)
  expect_true(all(!is.na(ci$bal)))
  if (any(ci$is_edge)) expect_true(all(is.na(ci$hegyi[ci$is_edge])))
})

test_that("unknown species fall back rather than fail", {
  tst <- data.frame(plot = "a", tree = c("1", "2"), species = c("ZZZ", "RS"),
                    dbh = c(20, 30), expf = 10)
  expect_message(
    s <- as_stand(tst, plot = "plot", tree = "tree", species = "species",
                  dbh = "dbh", expf = "expf"),
    "not in the trait table")
  expect_true(all(is.finite(s$trees$mcw)))
})
