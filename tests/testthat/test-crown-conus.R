# Tests for the CONUS crown width source. The Acadian path is covered by
# test-crown.R and by inst/scripts/regression_vs_original.R; these probes are
# about the second source and about the fallback defect the 0.1.2 sp_type fix
# did not reach.

test_that("the acadian source is unchanged by the arrival of a second source", {
  tr <- species_traits()
  expect_true(all(tr$cw_form == "acadian"))
  d <- c(5, 12.7, 25, 40, 80)
  for (sp in c("RS", "BF", "SM", "RO", "WP", "YB")) {
    p <- tr[tr$code == sp, ]
    expect_equal(max_crown_width(d, sp), p$mcw_a1 * d^p$mcw_a2)
    expect_equal(largest_crown_width(d, sp),
                 pmin((p$mcw_a1 * d^p$mcw_a2) / (p$lcw_b1 * d^p$lcw_b2),
                      p$mcw_a1 * d^p$mcw_a2))
  }
  # cr is accepted and ignored under the acadian form
  expect_equal(largest_crown_width(d, "RS"),
               largest_crown_width(d, "RS", cr = 0.3))
})

test_that("conus traits carry the whole crosswalk and nothing is silently lost", {
  tr <- species_traits(source = "conus")
  ac <- species_traits()
  expect_identical(sort(tr$code), sort(ac$code))
  expect_identical(nrow(tr), nrow(ac))
  expect_equal(sum(tr$cw_form == "conus"), 38L)
  expect_setequal(tr$code[tr$cw_form == "acadian"],
                  c("99", "AS", "HI", "OH", "OS"))
  cf <- tr[tr$cw_form == "conus", ]
  expect_true(all(!is.na(cf$spcd)))
  expect_true(all(cf$cw_provenance == "fitted"))
  expect_true(all(is.finite(cf$mcw_dbh_max) & cf$mcw_dbh_max > 0))
  expect_true(all(cf$mcw_ceiling %in% round(c(65.2, 75) * 0.3048, 10)))
  # every mapped species agrees with the shipped library on its coefficients
  L <- conus_crown_width
  i <- match(cf$spcd, L$spcd)
  expect_equal(cf$mcw_a1, L$mcw_a1[i])
  expect_equal(cf$mcw_a2, L$mcw_a2[i])
  expect_equal(cf$lcw_r0, L$lcw_r0[i])
})

test_that("the conus recipe holds the domain and then applies the ceiling", {
  tr <- species_traits(source = "conus")
  p <- tr[tr$code == "RS", ]
  raw <- function(d) p$mcw_a1 * d^p$mcw_a2
  dm <- p$mcw_dbh_max
  # below the hold the curve is the raw curve
  expect_equal(max_crown_width(c(10, 30), "RS", tr), raw(c(10, 30)))
  # at and above the hold it is flat at the held value
  expect_equal(max_crown_width(c(dm, dm + 1, 200), "RS", tr),
               rep(min(raw(dm), p$mcw_ceiling), 3))
  # the ceiling is never breached for any species at any diameter
  d <- c(1, 5, 20, 50, 100, 300)
  for (sp in tr$code[tr$cw_form == "conus"]) {
    m <- max_crown_width(d, sp, tr)
    expect_true(all(m <= tr$mcw_ceiling[tr$code == sp] + 1e-9))
  }
})

test_that("conus LCW is below MCW everywhere and responds to crown ratio", {
  tr <- species_traits(source = "conus")
  d <- c(5, 10, 20, 40, 60, 100)
  for (sp in c("RS", "BF", "SM", "RO", "WP", "YB", "PB", "NC")) {
    m <- max_crown_width(d, sp, tr)
    for (cr in c(0.05, 0.3, 0.6, 1)) {
      l <- largest_crown_width(d, sp, tr, mcw = m, cr = cr)
      expect_true(all(l < m))
      expect_true(all(l > 0))
    }
    # monotone increasing in crown ratio
    l1 <- largest_crown_width(d, sp, tr, mcw = m, cr = 0.2)
    l2 <- largest_crown_width(d, sp, tr, mcw = m, cr = 0.9)
    expect_true(all(l2 > l1))
  }
  # flat above the hold, exactly where MCW is flat
  dm <- tr$mcw_dbh_max[tr$code == "RS"]
  l <- largest_crown_width(c(dm, dm + 5, dm + 50), "RS", tr, cr = 0.5)
  expect_equal(l, rep(l[1], 3))
  # NA crown ratio is treated as the open grown level
  expect_equal(largest_crown_width(20, "RS", tr, cr = NA_real_),
               largest_crown_width(20, "RS", tr, cr = 1))
  # a percentage is rejected rather than silently used
  expect_error(largest_crown_width(20, "RS", tr, cr = 50), "proportion")
})

test_that("an unknown species falls back on its own type, not always hardwood", {
  tr <- species_traits()
  dflt <- forestCI:::type_defaults(tr)
  sw <- dflt[dflt$sp_type == "SW", ]
  hw <- dflt[dflt$sp_type == "HW", ]
  d <- 20

  # this is the 0.1.2 defect: without sp_type an unknown code took hardwood
  expect_equal(max_crown_width(d, "ZZZ", tr), hw$mcw_a1 * d^hw$mcw_a2)
  # and with sp_type it now takes the right one
  expect_equal(max_crown_width(d, "ZZZ", tr, sp_type = "SW"),
               sw$mcw_a1 * d^sw$mcw_a2)
  expect_equal(max_crown_width(d, "ZZZ", tr, sp_type = "HW"),
               hw$mcw_a1 * d^hw$mcw_a2)
  expect_false(isTRUE(all.equal(sw$mcw_a1 * d^sw$mcw_a2,
                                hw$mcw_a1 * d^hw$mcw_a2)))
})

test_that("as_stand routes sp_type into the crown width coefficients", {
  dat <- data.frame(plot = 1, tree = 1:2, species = c("ZZZ", "ZZZ"),
                    dbh = c(20, 20), type = c("SW", "HW"))
  s <- suppressMessages(as_stand(dat, sp_type = "type", plot_radius = 10))
  # the two trees must now differ, which is the whole point of the fix
  expect_false(isTRUE(all.equal(s$trees$mcw[1], s$trees$mcw[2])))
  # and each must equal the coefficients its own row carries
  expect_equal(s$trees$mcw, s$trees$mcw_a1 * s$trees$dbh^s$trees$mcw_a2)
})

test_that("as_stand carries the conus source through to the indices", {
  dat <- data.frame(plot = 1, tree = 1:6,
                    species = c("RS", "BF", "SM", "RO", "WP", "YB"),
                    dbh = c(10, 18, 26, 34, 42, 50),
                    height = c(9, 13, 17, 21, 24, 26),
                    hcb = c(4, 6, 8, 9, 10, 11))
  a <- as_stand(dat, height = "height", hcb = "hcb", plot_radius = 11.28)
  k <- as_stand(dat, traits = species_traits(source = "conus"),
                height = "height", hcb = "hcb", plot_radius = 11.28)
  expect_true(all(k$trees$cw_form == "conus"))
  expect_false(isTRUE(all.equal(a$trees$mcw, k$trees$mcw)))
  # LCW must respond to the stand's own crown ratios under the conus source
  expect_true(all(k$trees$lcw < k$trees$mcw))
  expect_equal(k$trees$cr, (dat$height - dat$hcb) / dat$height)
  open <- largest_crown_width(k$trees$dbh, k$trees$species, k$traits,
                              mcw = k$trees$mcw, cr = 1)
  expect_true(all(k$trees$lcw < open))
  # crown competition factor still computes and is finite on both sources
  for (st in list(a, k)) {
    ci <- competition_indices(st, rule = "fixed_radius", radius = 6)
    expect_true(all(is.finite(ci$trees$mca)))
  }
})

test_that("max_crown_area and as_stand use one definition of crown area", {
  dat <- data.frame(plot = 1, tree = 1:3, species = c("RS", "SM", "WP"),
                    dbh = c(12, 28, 44))
  s <- as_stand(dat, plot_radius = 10)
  expect_equal(s$trees$mca,
               max_crown_area(s$trees$dbh, s$trees$species, s$trees$expf))
})

test_that("the shipped conus library is internally consistent", {
  L <- conus_crown_width
  expect_equal(nrow(L), 466L)
  expect_true(all(!is.na(L$mcw_a1) & L$mcw_a1 > 0))
  expect_true(all(L$mcw_a2 > 0 & L$mcw_a2 < 1.5))
  expect_true(all(L$dbh_min_fit > 0 & L$dbh_max_fit >= L$dbh_min_fit))
  # 17 of 466 rows carry a single observation, so the fitted range collapses
  # to a point. All are donor carried or borrowed and none is an FIA code
  # that forestCI maps to, so no shipped prediction is pinned by them.
  degen <- L$dbh_max_fit <= L$dbh_min_fit
  expect_equal(sum(degen), 17L)
  expect_true(all(L$n_obs[degen] == 1))
  expect_equal(L$cw_shipped_max,
               pmin(L$mcw_a1 * L$dbh_max_fit^L$mcw_a2, L$cw_ceiling))
  expect_setequal(unique(L$clade), c("S", "H"))
  expect_setequal(unique(L$provenance),
                  c("fitted", "genus_donor", "spgrp_donor", "clade_population"))
  # the ratio slopes are global, which is a documented limitation, so assert it
  expect_equal(length(unique(L$lcw_rcr)), 1L)
  expect_equal(length(unique(L$lcw_rdbh)), 1L)
  expect_gt(length(unique(L$lcw_r0)), 200L)
  # no plot level or coordinate field may ever appear here
  expect_false(any(grepl("lat|lon|coord|plt|plot", names(L), ignore.case = TRUE)))
})
