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

test_that("node aligned rAPA covers the inclusive grid", {
  # node alignment samples -extent to +extent inclusive, so it has one more
  # point per axis than centre alignment and counts a larger nominal area
  rn <- rapa(st, weight = "none", resolution = 0.5, boundary = "square",
             extent = 25, align = "node")
  rc <- rapa(st, weight = "none", resolution = 0.5, boundary = "square",
             extent = 25, align = "center")
  tot_n <- tapply(rn$rapa, rn$plot_id, sum)
  tot_c <- tapply(rc$rapa, rc$plot_id, sum)
  expect_equal(as.numeric(tot_c), rep(50^2, length(tot_c)), tolerance = 1e-6)
  expect_equal(as.numeric(tot_n), rep(101^2 * 0.25, length(tot_n)),
               tolerance = 1e-6)
  m <- merge(rn, rc, by = c("plot_id", "tree_id"), suffixes = c("_n", "_c"))
  expect_gt(stats::cor(m$rapa_n, m$rapa_c), 0.999)
})

test_that("square and circular APA domains differ only by the clipped corners", {
  ac <- apa(st, weight = "none", boundary = "circle")
  aq <- apa(st, weight = "none", boundary = "square", extent = 25)
  expect_true(all(aq$apa >= ac$apa - 1e-6))
})

test_that("apa_bounded marks the trees whose polygon no neighbour closes", {
  g <- expand.grid(x = c(-5, 0, 5), y = c(-5, 0, 5))
  d <- data.frame(plot = "a", tree = as.character(seq_len(9)), species = "RS",
                  dbh = 20, x = g$x, y = g$y)
  s <- as_stand(d, plot = "plot", tree = "tree", species = "species",
                dbh = "dbh", x = "x", y = "y", plot_radius = 12, quiet = TRUE)
  a <- apa(s, weight = "none")
  expect_type(a$apa_bounded, "logical")
  # on a 3 by 3 grid only the centre is closed by neighbours; edge midpoints
  # have neighbours spanning exactly a half turn, which leaves the cell open
  expect_identical(a$tree_id[a$apa_bounded], "5")
  expect_equal(a$apa[a$tree_id == "5"], 25, tolerance = 1e-10)

  a0 <- apa(st, weight = "none")
  a1 <- apa(st, weight = "dbh", exponent = 2, boundary = "square", extent = 25)
  m <- merge(a0, a1, by = c("plot_id", "tree_id"))
  expect_identical(m$apa_bounded.x, m$apa_bounded.y)
  tr <- merge(st$trees[, c("plot_id", "tree_id", "x", "y")], a0,
              by = c("plot_id", "tree_id"))
  for (p in unique(tr$plot_id)) {
    k <- tr[tr$plot_id == p, ]
    hull <- k$tree_id[grDevices::chull(k$x, k$y)]
    expect_false(any(k$apa_bounded[k$tree_id %in% hull]), info = p)
    expect_true(any(k$apa_bounded), info = p)
  }
})

test_that("the APA early exit does not change any answer", {
  set.seed(5)
  n <- 60; r <- 20
  th <- stats::runif(n, 0, 2 * pi); rad <- r * 0.9 * sqrt(stats::runif(n))
  d <- data.frame(plot = "a", tree = as.character(seq_len(n)), species = "RS",
                  dbh = stats::runif(n, 10, 50),
                  x = rad * cos(th), y = rad * sin(th))
  s <- as_stand(d, plot = "plot", tree = "tree", species = "species",
                dbh = "dbh", x = "x", y = "y", plot_radius = r, quiet = TRUE)
  a0 <- apa(s, weight = "none")
  a1 <- apa(s, weight = "dbh", exponent = 2)
  expect_equal(sum(a0$apa), pi * r^2, tolerance = 0.01)
  expect_true(all(a0$apa > 0))
  expect_true(all(a1$apa > 0))
  expect_true(all(a0$apa <= pi * r^2 + 1e-6))
})

test_that("crown exposure is unchanged by the neighbour restriction", {
  far <- data.frame(plot = "a", tree = c("1", "2"), species = "RS",
                    dbh = c(25, 25), x = c(0, 80), y = c(0, 0),
                    height = 18, hcb = 9)
  sf <- as_stand(far, plot = "plot", tree = "tree", species = "species",
                 dbh = "dbh", height = "height", hcb = "hcb", x = "x", y = "y",
                 plot_radius = 60, quiet = TRUE)
  cf <- crown_exposure(sf, n_depth = 8L, n_azimuth = 24L, osv = TRUE,
                       ray_steps = 10L)
  expect_equal(cf$csax_rel, c(1, 1), tolerance = 0.02)
  expect_true(all(cf$osv_rel > 0.4))

  tight <- data.frame(plot = "a", tree = as.character(1:4), species = "RS",
                      dbh = c(45, 15, 15, 15),
                      x = c(0, 1.5, -1.5, 0), y = c(0, 0, 0, 1.5),
                      height = c(24, 12, 12, 12), hcb = c(10, 6, 6, 6))
  stt <- as_stand(tight, plot = "plot", tree = "tree", species = "species",
                  dbh = "dbh", height = "height", hcb = "hcb", x = "x", y = "y",
                  plot_radius = 12, quiet = TRUE)
  ct <- crown_exposure(stt, n_depth = 8L, n_azimuth = 24L, osv = FALSE)
  expect_gt(ct$csax_rel[ct$tree_id == "1"], 0.9)
  expect_lt(min(ct$csax_rel[ct$tree_id != "1"]), 0.999)
})

test_that("a caller supplied sp_type steers the fallback and survives it", {
  d <- data.frame(plot = "a", tree = c("1", "2"), species = c("ZZZ", "ZZZ"),
                  typ = c("SW", "HW"), dbh = c(20, 20), expf = 10)
  s <- as_stand(d, plot = "plot", tree = "tree", species = "species",
                dbh = "dbh", sp_type = "typ", expf = "expf", quiet = TRUE)
  expect_identical(s$trees$sp_type, c("SW", "HW"))
  expect_equal(length(unique(s$trees$widest)), 2L)

  # a recognised code takes its type from the trait table, not from the caller
  d2 <- data.frame(plot = "a", tree = c("1", "2"), species = c("RS", "SM"),
                   typ = c("HW", "SW"), dbh = c(20, 20), expf = 10)
  s2 <- as_stand(d2, plot = "plot", tree = "tree", species = "species",
                 dbh = "dbh", sp_type = "typ", expf = "expf", quiet = TRUE)
  expect_identical(s2$trees$sp_type, c("SW", "HW"))

  expect_error(
    as_stand(data.frame(plot = "a", tree = "1", species = "RS",
                        typ = "conifer", dbh = 20, expf = 10),
             plot = "plot", tree = "tree", species = "species", dbh = "dbh",
             sp_type = "typ", expf = "expf", quiet = TRUE),
    "SW")
})

test_that("an empty neighbourhood means zero competition, not unknown", {
  # six trees far enough apart that none tallies another
  d <- data.frame(plot = "a", tree = as.character(1:6), species = "RS",
                  dbh = c(12, 18, 24, 30, 36, 42),
                  x = c(-40, -20, 0, 20, 40, 0),
                  y = c(-40, 20, 0, -20, 40, 45),
                  height = 18, hcb = 9)
  s <- as_stand(d, plot = "plot", tree = "tree", species = "species",
                dbh = "dbh", height = "height", hcb = "hcb", x = "x", y = "y",
                plot_radius = 60, quiet = TRUE)
  nb <- neighbors(s, method = "radius", radius = 0.5)
  expect_equal(nrow(nb), 0L)
  v <- ci_distance_dependent(s, nb = nb,
                             indices = c("hegyi", "martin_ek", "spurr",
                                         "local_ba", "local_bal", "n_comp",
                                         "mean_dist"))
  expect_equal(v$hegyi, rep(0, 6))
  expect_equal(v$martin_ek, rep(0, 6))
  expect_equal(v$n_comp, rep(0, 6))
  expect_true(all(is.na(v$mean_dist)))
  expect_true(all(v$local_ba > 0))          # a tree still occupies its own space
  expect_equal(v$local_bal, rep(0, 6))

  # and it must agree with the per-tree empty case in a stand that has both
  d2 <- rbind(d, data.frame(plot = "a", tree = "7", species = "RS", dbh = 20,
                            x = 0.4, y = 0, height = 18, hcb = 9))
  s2 <- as_stand(d2, plot = "plot", tree = "tree", species = "species",
                 dbh = "dbh", height = "height", hcb = "hcb", x = "x", y = "y",
                 plot_radius = 60, quiet = TRUE)
  v2 <- ci_distance_dependent(s2, radius = 0.5, indices = c("hegyi", "n_comp"))
  lone <- v2[v2$tree_id %in% c("1", "2", "4", "5", "6"), ]
  expect_equal(lone$hegyi, rep(0, 5))
  expect_equal(lone$n_comp, rep(0, 5))
})
