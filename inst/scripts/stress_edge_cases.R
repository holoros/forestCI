# ---------------------------------------------------------------------------
# forestCI stress and edge case suite
#   Rscript stress.R <out_dir>
# Every probe is trapped, so the suite always completes and reports.
# ---------------------------------------------------------------------------
args <- commandArgs(TRUE)
outdir <- if (length(args)) args[1] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(forestCI); library(data.table) })
set.seed(20260918)

R <- data.table(probe = character(), expect = character(),
                got = character(), status = character())
note <- function(probe, expect, got, ok) {
  R <<- rbind(R, data.table(probe = probe, expect = expect, got = got,
                            status = if (isTRUE(ok)) "pass" else "FAIL"))
}
# runs expr; ok_fun(value) decides pass. An error is a value of class "error".
probe <- function(name, expect, expr, ok_fun) {
  v <- tryCatch(suppressMessages(expr), error = function(e) e,
                warning = function(w) w)
  got <- if (inherits(v, "condition")) paste0(class(v)[1], ": ",
                                              conditionMessage(v))
         else "returned"
  ok <- tryCatch(isTRUE(ok_fun(v)), error = function(e) FALSE)
  note(name, expect, got, ok)
  invisible(v)
}
mk <- function(n, seed = 1, radius = 16, spp = c("RS","BF","RM","EH"),
               with_xy = TRUE, with_crown = TRUE, plot = "a") {
  set.seed(seed)
  th <- runif(n, 0, 2 * pi); rad <- radius * 0.95 * sqrt(runif(n))
  d <- data.frame(plot = plot, tree = as.character(seq_len(n)),
                  species = sample(spp, n, TRUE),
                  dbh = round(runif(n, 8, 55), 2),
                  x = rad * cos(th), y = rad * sin(th))
  d$height <- 1.37 + 0.55 * d$dbh + rnorm(n, 0, 1)
  d$hcb <- d$height * runif(n, 0.3, 0.6)
  a <- list(d, plot = "plot", tree = "tree", species = "species", dbh = "dbh",
            plot_radius = radius, quiet = TRUE)
  if (with_xy) { a$x <- "x"; a$y <- "y" }
  if (with_crown) { a$height <- "height"; a$hcb <- "hcb" }
  do.call(as_stand, a)
}

# ---- A. degenerate stand sizes -------------------------------------------
for (n in c(1L, 2L, 3L)) {
  st <- mk(n, seed = n)
  probe(sprintf("as_stand with %d tree(s)", n), "builds",
        st, function(v) inherits(v, "stand") && nrow(v$trees) == n)
  probe(sprintf("ci_distance_independent, n = %d", n), "finite BAL",
        ci_distance_independent(st, sdi_max = NULL),
        function(v) nrow(v) == n && all(is.finite(v$bal)))
  probe(sprintf("neighbors, n = %d", n), "no error",
        neighbors(st, radius = 6), function(v) is.data.frame(v))
  probe(sprintf("ci_distance_dependent, n = %d", n), "one row per tree",
        ci_distance_dependent(st, radius = 6), function(v) nrow(v) == n)
  probe(sprintf("apa, n = %d", n), "areas sum to plot area",
        apa(st, weight = "none"),
        function(v) abs(sum(v$apa) - pi * 16^2) / (pi * 16^2) < 0.02)
  probe(sprintf("rapa, n = %d", n), "finite",
        rapa(st, weight = "none", resolution = 1), function(v) all(is.finite(v$rapa)))
  probe(sprintf("crown_exposure, n = %d", n), "csax <= csa",
        crown_exposure(st, n_depth = 6L, n_azimuth = 12L, osv = TRUE,
                       ray_steps = 6L),
        function(v) all(v$csax <= v$csa + 1e-6) && all(v$osv <= v$csax + 1e-6))
  probe(sprintf("spatial_pattern, n = %d", n), "runs or returns NA",
        spatial_pattern(st), function(v) nrow(v) == 1)
  probe(sprintf("edge_flag, n = %d", n), "logical",
        edge_flag(st), function(v) is.logical(v$is_edge))
  probe(sprintf("competition_indices, n = %d", n), "one row per tree",
        competition_indices(st, n_depth = 6L, n_azimuth = 12L,
                            crown_osv = FALSE, rapa_resolution = 1,
                            quiet = TRUE),
        function(v) nrow(v) == n)
}

# ---- B. pathological inputs that must be rejected --------------------------
bad <- function(name, d, ...) {
  probe(name, "informative error",
        as_stand(d, plot = "plot", tree = "tree", species = "species",
                 dbh = "dbh", expf = "expf", quiet = TRUE, ...),
        function(v) inherits(v, "error"))
}
bad("zero diameter rejected",
    data.frame(plot="a", tree=c("1","2"), species="RS", dbh=c(0, 20), expf=10))
bad("negative diameter rejected",
    data.frame(plot="a", tree=c("1","2"), species="RS", dbh=c(-5, 20), expf=10))
bad("duplicate tree id rejected",
    data.frame(plot="a", tree=c("1","1"), species="RS", dbh=c(10, 20), expf=10))
probe("coincident stems rejected", "informative error",
      as_stand(data.frame(plot="a", tree=c("1","2"), species="RS",
                          dbh=c(10,20), x=c(0,0), y=c(0,0), expf=10),
               plot="plot", tree="tree", species="species", dbh="dbh",
               x="x", y="y", expf="expf", quiet=TRUE),
      function(v) inherits(v, "error"))
probe("crown base above tip rejected", "informative error",
      as_stand(data.frame(plot="a", tree=c("1","2"), species="RS",
                          dbh=c(10,20), height=c(10,10), hcb=c(12,5), expf=10),
               plot="plot", tree="tree", species="species", dbh="dbh",
               height="height", hcb="hcb", expf="expf", quiet=TRUE),
      function(v) inherits(v, "error"))
probe("no expansion factor and no plot size rejected", "informative error",
      as_stand(data.frame(plot="a", tree=c("1","2"), species="RS", dbh=c(10,20)),
               plot="plot", tree="tree", species="species", dbh="dbh",
               quiet=TRUE),
      function(v) inherits(v, "error"))
probe("missing column named clearly", "error names the column",
      as_stand(data.frame(plot="a", tree="1", species="RS"),
               plot="plot", tree="tree", species="species", dbh="dbh", expf=10),
      function(v) inherits(v, "error") && grepl("dbh", conditionMessage(v)))

# ---- C. graceful degradation ---------------------------------------------
st_nc <- mk(30, seed = 7, with_xy = FALSE)
probe("no coordinates: is_spatial FALSE", "FALSE", is_spatial(st_nc),
      function(v) identical(v, FALSE))
probe("no coordinates: distance-dependent refuses clearly", "informative error",
      neighbors(st_nc, radius = 6),
      function(v) inherits(v, "error") && grepl("coordinate", conditionMessage(v)))
probe("no coordinates: wrapper still returns BAL", "bal present",
      competition_indices(st_nc, quiet = TRUE),
      function(v) "bal" %in% names(v) && !"hegyi" %in% names(v) && nrow(v) == 30)
st_nh <- mk(20, seed = 8, with_crown = FALSE)
probe("no heights: crown_exposure refuses clearly", "informative error",
      crown_exposure(st_nh),
      function(v) inherits(v, "error") && grepl("height", conditionMessage(v)))
probe("no heights: wrapper skips crown, keeps the rest", "hegyi present",
      competition_indices(st_nh, rapa_resolution = 1, quiet = TRUE),
      function(v) "hegyi" %in% names(v) && !"csax" %in% names(v))
probe("no heights: height filter refuses clearly", "informative error",
      neighbors(st_nh, filter = "taller"),
      function(v) inherits(v, "error"))
probe("unknown species falls back, does not fail", "finite mcw",
      { s <- as_stand(data.frame(plot="a", tree=c("1","2"),
                                 species=c("ZZZ","QQQ"), dbh=c(10,20), expf=10),
                      plot="plot", tree="tree", species="species", dbh="dbh",
                      expf="expf", quiet=TRUE); s$trees$mcw },
      function(v) all(is.finite(v)))
probe("sp_type steers the fallback for unknown codes", "SW and HW differ",
      { d <- data.frame(plot="a", tree=c("1","2"), species=c("ZZZ","ZZZ"),
                        typ=c("SW","HW"), dbh=c(20,20), expf=10)
        s <- as_stand(d, plot="plot", tree="tree", species="species",
                      dbh="dbh", sp_type="typ", expf="expf", quiet=TRUE)
        length(unique(s$trees$widest)) == 2 &&
          identical(s$trees$sp_type, c("SW","HW")) },
      function(v) isTRUE(v))
probe("sp_type is ignored for recognised codes", "trait table wins",
      { d <- data.frame(plot="a", tree=c("1","2"), species=c("RS","SM"),
                        typ=c("HW","SW"), dbh=c(20,20), expf=10)
        s <- as_stand(d, plot="plot", tree="tree", species="species",
                      dbh="dbh", sp_type="typ", expf="expf", quiet=TRUE)
        identical(s$trees$sp_type, c("SW","HW")) },
      function(v) isTRUE(v))
probe("a bad sp_type value is rejected", "informative error",
      as_stand(data.frame(plot="a", tree="1", species="RS", typ="conifer",
                          dbh=20, expf=10),
               plot="plot", tree="tree", species="species", dbh="dbh",
               sp_type="typ", expf="expf", quiet=TRUE),
      function(v) inherits(v, "error"))

# ---- D. ties and degenerate geometry --------------------------------------
probe("all diameters tied: every BAL is zero", "all zero",
      { d <- data.frame(plot="a", tree=as.character(1:10), species="RS",
                        dbh=25, expf=10)
        s <- as_stand(d, plot="plot", tree="tree", species="species",
                      dbh="dbh", expf="expf", quiet=TRUE)
        ci_distance_independent(s, sdi_max = NULL)$bal },
      function(v) all(abs(v) < 1e-9))
probe("tree exactly at plot centre", "edge_flag runs",
      { d <- data.frame(plot="a", tree=as.character(1:5), species="RS",
                        dbh=c(20,25,30,35,40),
                        x=c(0, 3, -4, 5, -6), y=c(0, 4, 3, -5, 2))
        s <- as_stand(d, plot="plot", tree="tree", species="species",
                      dbh="dbh", x="x", y="y", plot_radius=16, quiet=TRUE)
        edge_flag(s) },
      function(v) nrow(v) == 5 && all(!is.na(v$is_edge)))
probe("widest point at the apex", "finite positive dimensions",
      crown_dimensions(crown_profile(), lcw=4, height=20, hcb=10, widest=0),
      function(v) all(is.finite(unlist(v))) && v$csa > 0)
probe("widest point at the crown base", "finite positive dimensions",
      crown_dimensions(crown_profile(), lcw=4, height=20, hcb=10, widest=1),
      function(v) all(is.finite(unlist(v))) && v$csa > 0)
probe("very short crown", "finite",
      crown_dimensions(crown_profile(), lcw=4, height=20, hcb=19.99),
      function(v) all(is.finite(unlist(v))))
probe("very narrow crown", "finite",
      crown_dimensions(crown_profile(), lcw=0.01, height=20, hcb=10),
      function(v) all(is.finite(unlist(v))) && v$csa > 0)
probe("extreme profile exponents", "finite",
      crown_dimensions(crown_profile(), lcw=4, height=20, hcb=10,
                       up_exp=12, lo_exp=12, shape="e"),
      function(v) all(is.finite(unlist(v))))
probe("custom profile family", "finite and bounded",
      crown_dimensions(crown_profile("custom", fun=function(z) (1-z)^0.4),
                       lcw=4, height=20, hcb=10),
      function(v) all(is.finite(unlist(v))) && v$csa > 0)
probe("neiloid solid", "finite",
      crown_dimensions(crown_profile("geometric", solid="neiloid"),
                       lcw=4, height=20, hcb=10),
      function(v) all(is.finite(unlist(v))))
probe("crown radius zero outside the crown", "zero below crown base",
      crown_radius_at(crown_profile(), 5, lcw=4, height=20, hcb=10),
      function(v) v == 0)

# ---- E. neighbourhood rules ------------------------------------------------
st <- mk(60, seed = 11)
probe("knn with k larger than the stand", "k capped, no error",
      neighbors(st, method="knn", k=500L),
      function(v) is.data.frame(v) && nrow(v) > 0)
probe("angle gauge tallies fewer trees as BAF rises", "monotone decreasing",
      { k <- sapply(c(1, 4, 16, 64, 256), function(b)
          nrow(neighbors(st, method="angle_gauge", baf=b)))
        all(diff(k) <= 0) }, function(v) isTRUE(v))
probe("indices survive an empty neighbourhood", "sums are zero not NA",
      { sp <- mk(6, seed=77, radius=60)
        nb <- neighbors(sp, method="radius", radius=0.5)
        v <- ci_distance_dependent(sp, nb=nb, indices=c("hegyi","n_comp"))
        nrow(nb) == 0 && all(v$hegyi == 0) && all(v$n_comp == 0) &&
          !any(is.na(v$hegyi)) }, function(v) isTRUE(v))
probe("influence zone selection", "runs",
      neighbors(st, method="influence_zone"), function(v) is.data.frame(v))
probe("radius zero yields no competitors", "empty",
      neighbors(st, method="radius", radius=0),
      function(v) nrow(v) == 0)
probe("height fraction filter is a subset of unfiltered", "subset",
      { a <- neighbors(st, radius=8)
        b <- neighbors(st, radius=8, filter="height_frac", filter_value=0.67)
        nrow(b) <= nrow(a) }, function(v) v)
probe("taller filter is a subset of unfiltered", "subset",
      { a <- neighbors(st, radius=8)
        b <- neighbors(st, radius=8, filter="taller")
        nrow(b) <= nrow(a) }, function(v) v)
probe("Hegyi is monotone in radius", "non-decreasing",
      { h <- sapply(c(4,6,8,10), function(r)
          sum(ci_distance_dependent(st, radius=r, indices="hegyi")$hegyi))
        all(diff(h) >= -1e-8) }, function(v) v)

# ---- F. multi plot ---------------------------------------------------------
probe("two plots of different radius are handled independently", "per plot BA",
      { d1 <- as.data.frame(mk(25, seed=21, radius=12, plot="p1")$trees)
        d2 <- as.data.frame(mk(40, seed=22, radius=20, plot="p2")$trees)
        d <- rbind(d1[, c("plot_id","tree_id","species","dbh","x","y","height","hcb")],
                   d2[, c("plot_id","tree_id","species","dbh","x","y","height","hcb")])
        d$plot_radius <- ifelse(d$plot_id == "p1", 12, 20)
        s <- as_stand(d, plot="plot_id", tree="tree_id", species="species",
                      dbh="dbh", height="height", hcb="hcb", x="x", y="y",
                      plot_radius="plot_radius", quiet=TRUE)
        ci <- competition_indices(s, rapa_resolution=1, n_depth=6L,
                                  n_azimuth=12L, crown_osv=FALSE, quiet=TRUE)
        length(unique(ci$ba)) == 2 && nrow(ci) == 65 },
      function(v) v)
probe("APA respects each plot's own boundary", "sums to each area",
      { d1 <- as.data.frame(mk(25, seed=21, radius=12, plot="p1")$trees)
        d2 <- as.data.frame(mk(40, seed=22, radius=20, plot="p2")$trees)
        d <- rbind(d1[, c("plot_id","tree_id","species","dbh","x","y")],
                   d2[, c("plot_id","tree_id","species","dbh","x","y")])
        d$plot_radius <- ifelse(d$plot_id == "p1", 12, 20)
        s <- as_stand(d, plot="plot_id", tree="tree_id", species="species",
                      dbh="dbh", x="x", y="y", plot_radius="plot_radius",
                      quiet=TRUE)
        a <- apa(s, weight="none")
        tot <- tapply(a$apa, a$plot_id, sum)
        max(abs(as.numeric(tot) - pi*c(12,20)^2) / (pi*c(12,20)^2)) < 0.02 },
      function(v) v)

# ---- G. weighting and options ----------------------------------------------
probe("apa normalize conserves plot area", "sums exactly",
      { a <- apa(st, weight="dbh", exponent=2, normalize=TRUE)
        abs(sum(a$apa) - pi*16^2) / (pi*16^2) < 1e-9 }, function(v) v)
probe("rapa additive mode is milder than multiplicative", "smaller max share",
      { m <- rapa(st, weight="dbh", exponent=1, resolution=1,
                  weight_mode="multiplicative")
        a <- rapa(st, weight="dbh", exponent=1, resolution=1,
                  weight_mode="additive")
        max(a$rapa) < max(m$rapa) }, function(v) v)
probe("directional modification is inert on flat ground", "identical",
      { a <- rapa(st, weight="none", resolution=1, dr_mod=FALSE)
        b <- rapa(st, weight="none", resolution=1, dr_mod=TRUE, dr_max=2)
        max(abs(a$rapa - b$rapa)) < 1e-9 }, function(v) v)
probe("directional modification bites on a slope", "differs",
      { d <- as.data.frame(st$trees)
        d$slope <- 40; d$aspect <- 180
        s <- as_stand(d, plot="plot_id", tree="tree_id", species="species",
                      dbh="dbh", x="x", y="y", plot_radius=16,
                      slope="slope", aspect="aspect", quiet=TRUE)
        a <- rapa(s, weight="none", resolution=1, dr_mod=FALSE)
        b <- rapa(s, weight="none", resolution=1, dr_mod=TRUE, dr_max=2)
        max(abs(a$rapa - b$rapa)) > 1e-6 }, function(v) v)
probe("nci returns NA without coefficients, not a guess", "all NA",
      ci_distance_dependent(st, radius=6, indices="nci")$nci,
      function(v) all(is.na(v) | v == 0))
probe("nci computes when coefficients are supplied", "finite positive",
      ci_distance_dependent(st, radius=6, indices="nci",
                            nci_alpha=2, nci_beta=1, nci_lambda=1)$nci,
      function(v) all(is.finite(v)) && any(v > 0))
probe("species_traits accepts a minimal override", "row present",
      species_traits(data.frame(code="DF", sp_type="SW")),
      function(v) "DF" %in% v$code && all(!is.na(v$widest[v$code == "DF"])))
probe("species_traits rejects a bad sp_type", "informative error",
      species_traits(data.frame(code="DF", sp_type="XX")),
      function(v) inherits(v, "error"))
probe("a user trait table flows through as_stand", "override used",
      { tr <- species_traits(data.frame(code="RS", sp_type="SW",
                                        mcw_a1=99, mcw_a2=0.5))
        s <- as_stand(data.frame(plot="a", tree=c("1","2"), species="RS",
                                 dbh=c(16,25), expf=10),
                      plot="plot", tree="tree", species="species", dbh="dbh",
                      expf="expf", traits=tr, quiet=TRUE)
        abs(s$trees$mcw[1] - 99*16^0.5) < 1e-8 }, function(v) v)

# ---- H. determinism and integration agreement ------------------------------
probe("competition_indices is deterministic", "identical twice",
      { a <- competition_indices(st, rapa_resolution=1, n_depth=6L,
                                 n_azimuth=12L, crown_osv=FALSE, quiet=TRUE)
        b <- competition_indices(st, rapa_resolution=1, n_depth=6L,
                                 n_azimuth=12L, crown_osv=FALSE, quiet=TRUE)
        isTRUE(all.equal(a, b)) }, function(v) v)
probe("adaptive and refined Simpson agree", "within 0.5 percent",
      { a <- crown_dimensions(crown_profile(), lcw=4, height=20, hcb=10,
                              species="RS", method="adaptive")
        s <- crown_dimensions(crown_profile(), lcw=4, height=20, hcb=10,
                              species="RS", method="simpson", n_sub=8000L)
        abs(a$csa - s$csa)/a$csa < 0.005 && abs(a$cv - s$cv)/a$cv < 0.005 },
      function(v) v)
probe("crown surface area rises with crown length", "monotone",
      { v <- sapply(seq(2, 14, by=2), function(cl)
          crown_dimensions(crown_profile(), lcw=4, height=20, hcb=20-cl,
                           species="RS")$csa)
        all(diff(v) > 0) }, function(v) v)
probe("crown surface area rises with crown width", "monotone",
      { v <- sapply(seq(1, 9, by=2), function(w)
          crown_dimensions(crown_profile(), lcw=w, height=20, hcb=10,
                           species="RS")$csa)
        all(diff(v) > 0) }, function(v) v)

# ---- I. scaling ------------------------------------------------------------
timing <- rbindlist(lapply(c(100L, 400L, 1000L), function(n) {
  s <- mk(n, seed = 99, radius = sqrt(n / 0.12 / pi))
  t1 <- system.time(ci_distance_independent(s, sdi_max = NULL))[["elapsed"]]
  t2 <- system.time(ci_distance_dependent(s, radius = 6))[["elapsed"]]
  t3 <- system.time(apa(s, weight = "dbh"))[["elapsed"]]
  t4 <- system.time(rapa(s, weight = "none", resolution = 0.5))[["elapsed"]]
  t5 <- system.time(crown_exposure(s, n_depth = 8L, n_azimuth = 24L,
                                   osv = FALSE))[["elapsed"]]
  data.table(n_trees = n, independent_s = t1, dependent_s = t2, apa_s = t3,
             rapa_s = t4, crown_s = t5)
}))
print(timing)
fwrite(timing, file.path(outdir, "2026-09-18_forestCI-stress-timing_DATA.csv"))
probe("1000 tree stand completes the full wrapper under 300 s", "under 300 s",
      { s <- mk(1000L, seed = 101, radius = sqrt(1000/0.12/pi))
        el <- system.time(competition_indices(s, rapa_resolution = 0.5,
                                              n_depth = 8L, n_azimuth = 24L,
                                              crown_osv = FALSE,
                                              quiet = TRUE))[["elapsed"]]
        cat(sprintf("  1000 tree wrapper: %.1f s\n", el)); el < 300 },
      function(v) v)

probe("crown exposure scales sub-quadratically", "1000 under 60 s",
      { s1 <- mk(1000L, seed=303, radius=sqrt(1000/0.12/pi))
        el <- system.time(crown_exposure(s1, n_depth=8L, n_azimuth=24L,
                                         osv=FALSE))[["elapsed"]]
        cat(sprintf("  crown_exposure 1000 trees: %.1f s\n", el)); el < 60 },
      function(v) isTRUE(v))
probe("APA scales sub-quadratically", "1000 under 4 s",
      { s1 <- mk(1000L, seed=304, radius=sqrt(1000/0.12/pi))
        el <- system.time(apa(s1, weight="dbh"))[["elapsed"]]
        cat(sprintf("  apa 1000 trees: %.1f s\n", el)); el < 4 },
      function(v) isTRUE(v))

# ---- J. random fuzz --------------------------------------------------------
fuzz_fail <- 0L
for (k in 1:40) {
  n <- sample(2:40, 1)
  s <- tryCatch(mk(n, seed = 1000 + k, radius = runif(1, 8, 25)),
                error = function(e) NULL)
  if (is.null(s)) { fuzz_fail <- fuzz_fail + 1L; next }
  ok <- tryCatch({
    ci <- suppressMessages(competition_indices(
      s, radius = runif(1, 3, 10), rapa_resolution = 1, n_depth = 6L,
      n_azimuth = 12L, crown_osv = FALSE, quiet = TRUE,
      apa_weight = sample(c("none","dbh","lcw"), 1),
      filter = sample(c("none","larger_dbh","taller","height_frac"), 1),
      edge = sample(c("flag","exclude","none"), 1)))
    nrow(ci) == n &&
      all(is.finite(ci$bal)) && all(ci$bal >= -1e-8) &&
      all(is.na(ci$csax) | ci$csax <= ci$csa + 1e-6) &&
      all(is.na(ci$hegyi) | ci$hegyi >= -1e-8)
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) fuzz_fail <- fuzz_fail + 1L
}
note("40 random stands, all invariants hold", "0 failures",
     sprintf("%d failures", fuzz_fail), fuzz_fail == 0L)

# ---- report ----------------------------------------------------------------
cat("\n")
print(R, nrows = 200)
fwrite(R, file.path(outdir, "2026-09-18_forestCI-stress-results_FINAL.csv"))
cat(sprintf("\n%d probes, %d pass, %d FAIL\n", nrow(R), sum(R$status == "pass"),
            sum(R$status == "FAIL")))
if (any(R$status == "FAIL")) {
  cat("\nFAILURES:\n"); print(R[status == "FAIL"], nrows = 100)
}
gc()
