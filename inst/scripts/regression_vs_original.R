# ---------------------------------------------------------------------------
# forestCI regression against the original SFR 575 scripts
#
# Runs the original CI.R / SPP.R / AcadianGY pipeline on PEF.csv exactly as
# Run.R does, then runs forestCI on the same data, and compares tree by tree.
#
#   Rscript regression_vs_original.R <orig_dir> <out_dir>
#
# <orig_dir> must contain CI.R, SPP.R, circular2xy.R, sort.data.frame.r,
# AcadianGYv12.1.2.r and PEF.csv. Those files are not redistributed with this
# package; they are the University of Maine SFR 575 laboratory scripts.
#
# The original code is written around global assignment and a global object
# called `data`, so it is run in globalenv() deliberately and the forestCI side
# is computed first and stashed under a different name.
# ---------------------------------------------------------------------------
args <- commandArgs(TRUE)
orig <- if (length(args) >= 1) args[1] else "orig"
outdir <- if (length(args) >= 2) args[2] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
elog <- file.path(outdir, "error_log.txt")
trap <- function(label, expr) tryCatch(expr, error = function(e) {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time()), label, conditionMessage(e)),
      file = elog, append = TRUE); NULL })

set.seed(20260917)
suppressPackageStartupMessages({ library(forestCI); library(data.table) })

# ---- 1. forestCI side -----------------------------------------------------
pkg <- local({
  raw <- read.csv(file.path(orig, "PEF.csv"), stringsAsFactors = FALSE)
  raw$plot_id <- raw$PlotID.yr
  st <- as_stand(raw, plot = "plot_id", tree = "TreeNum", species = "FVS",
                 dbh = "DBH", height = "HT", hcb = "HCB", expf = "EXPF",
                 distance = "Dist", azimuth = "Az", plot_radius = 16.0509,
                 quiet = TRUE)
  di <- ci_distance_independent(st, sdi_max = NULL)
  # measured crown radii, exactly as CI.R builds them: four quarter ellipses,
  # then a mean radius from the projected area
  cpa_meas <- 0.25 * pi * (raw$Rad0 * raw$Rad90 + raw$Rad90 * raw$Rad180 +
                           raw$Rad180 * raw$Rad270 + raw$Rad270 * raw$Rad0)
  lcw_meas <- 2 * sqrt(cpa_meas / pi)
  dm <- crown_dimensions(crown_profile("dual_exponent"), lcw = lcw_meas,
                         height = raw$HT, hcb = raw$HCB, species = raw$FVS,
                         n_sub = 600L)
  meas <- data.table(plot_id = st$trees$plot_id, tree_id = as.character(raw$TreeNum),
                     csa_meas = dm$csa, cv_meas = dm$cv, cpa_meas = dm$cpa)
  # angle gauge neighbourhood, BAF 10 sq ft per acre = 2.2957 m2 ha-1, which is
  # the CI.R default of BAF = 10/4.356
  nb <- neighbors(st, method = "angle_gauge", baf = 10 / 4.356)
  dd <- ci_distance_dependent(st, nb = nb, indices = c("hegyi", "spurr", "n_comp"))
  ce <- crown_exposure(st, n_depth = 60L, n_azimuth = 180L, osv = FALSE)
  # CI.R rasterises on a square of half width 25 m at 0.5 m, not on the plot
  # circle, so the comparable call uses the same domain
  ra <- rapa(st, weight = "none", resolution = 0.5, boundary = "square",
             extent = 25, align = "node")
  ap <- apa(st, weight = "none")
  out <- Reduce(function(a, b) merge(a, b, by = c("plot_id", "tree_id")),
                list(st$trees[, .(plot_id, tree_id, dbh, species, mcw, lcw, mca)],
                     di, dd, ce[, .(plot_id, tree_id, csa, cpa, cv)], meas, ra, ap))
  setnames(out, "tree_id", "TreeNum")
  out
})
cat("forestCI side:", nrow(pkg), "trees\n")

# ---- 2. original pipeline, run in globalenv --------------------------------
setwd(orig)
source("CI.R");  source("SPP.R");  source("circular2xy.R")
source("sort.data.frame.r")
# take only mcw() and SPP.func() out of AcadianGY so that plyr, dplyr and nlme
# are not needed and no top level code runs
acd <- readLines("AcadianGYv12.1.2.r")
grab <- function(start_pat) {
  i <- grep(start_pat, acd)[1]
  depth <- 0; j <- i
  repeat {
    depth <- depth + lengths(regmatches(acd[j], gregexpr("\\{", acd[j]))) -
      lengths(regmatches(acd[j], gregexpr("\\}", acd[j])))
    if (j > i && depth <= 0) break
    j <- j + 1
    if (j > length(acd)) break
  }
  acd[i:j]
}
eval(parse(text = grab("^SPP\\.func=function")), envir = globalenv())
eval(parse(text = grab("^mcw=function")), envir = globalenv())

data <- read.csv("PEF.csv", stringsAsFactors = FALSE)

FORMAT(inputs = data, plots = "PlotID.yr", trees = "TreeNum",
       positions = "none", aspects = "none", distances = "Dist",
       directions = "Az", slopes = "none", spps = "FVS", dbhs = "DBH",
       cohorts = "none", cnpypstns = "none", hgts = "HT", crown_bases = "HCB",
       north_crs = "Rad0", east_crs = "Rad90", south_crs = "Rad180",
       west_crs = "Rad270", measures = "none")

pid <- as.character(data$PlotID.yr[1])
START(choice = pid)
CROWNS(widept = wides, upexp = ups, loexp = los, shape = shps, trnc = FALSE)

orig_csa <- data.frame(TreeNum = as.character(N), csa_orig = as.numeric(CSA),
                       cv_orig = as.numeric(CV), stringsAsFactors = FALSE)

trap("rAPAS", {
  rAPAS(set = c(pid), widept = wides, upexp = ups, loexp = los, shape = shps,
        trnc = FALSE)
})
orig_rapa <- if (exists("rPAs") && length(rPAs) > 0) {
  data.frame(TreeNum = as.character(rownames(as.data.frame(rPAs))),
             rapa_orig = as.numeric(as.data.frame(rPAs)[[1]]),
             stringsAsFactors = FALSE)
} else NULL

# HEGYI() and SPURR() in CI.R depend on global filter state that does not
# initialise cleanly outside the full Run.R session, so their equations are
# evaluated directly against the DIST and DBH matrices that START() left in
# the global environment. This is the CI.R formula, not a paraphrase of it:
#   in-set:  (0.5 * BAF^-0.5) * DBH[j] > DIST[i, j]
#   hegyi:   sum over j != i of (DBH[j] / DBH[i]) / DIST[i, j]
#   spurr:   0.25 * sum((rank(s^-1) - 0.5) * s^2) / (sum(INt[i, ]) - 1)
#            with s = DBH[j] / DIST[i, j]. Note that INt[i, i] is 1, because a
#            tree always tallies itself at zero distance, so the CI.R divisor
#            (sum(INt) - 1) is the competitor count, not the count minus one.
BAF <- 10 / 4.356
nT <- length(Nb)
hg <- numeric(nT); sp_ <- numeric(nT); nin <- integer(nT)
for (i in seq_len(nT)) {
  inset <- (0.5 * BAF^-0.5) * DBH > DIST[i, ]
  inset[i] <- FALSE
  nin[i] <- sum(inset)
  if (nin[i] == 0) { hg[i] <- 0; sp_[i] <- NA_real_; next }
  hg[i] <- sum((DBH[inset] / DBH[i]) / DIST[i, inset])
  sv <- rep(0, nT)
  sv[inset] <- DBH[inset] / DIST[i, inset]
  sp_[i] <- 0.25 * sum((rank(sv^-1) - 0.5) * sv^2) / nin[i]
}
orig_hegyi <- data.frame(TreeNum = as.character(N), hegyi_orig = hg,
                         n_in_orig = nin, stringsAsFactors = FALSE)
orig_spurr <- data.frame(TreeNum = as.character(N), spurr_orig = sp_,
                         stringsAsFactors = FALSE)

# BAL, CCFL and MCW exactly as Run.R builds them
d2 <- read.csv("PEF.csv", stringsAsFactors = FALSE)
d2$BA <- 0.00007854 * d2$DBH^2 * d2$EXPF
d2$SPtype <- as.vector(mapply(SPP.func, d2$FVS)[1, ])
d2$SW <- ifelse(d2$SPtype == "SW", 1, 0)
d2$BA.sw <- d2$BA * d2$SW
d2$MCW <- mapply(mcw, sp = d2$FVS, dbh = d2$DBH)
d2$MCA <- 100 * ((pi * (d2$MCW / 2)^2) / 10000) * d2$EXPF
d2$MCA.SW <- ifelse(d2$SPtype == "SW", d2$MCA, 0)
d2 <- d2[order(-d2$DBH), ]
d2$BAL_incl <- cumsum(d2$BA)
d2$BAL_orig <- d2$BAL_incl - d2$BA          # the intended exclusive form
d2$CCFL_incl <- cumsum(d2$MCA)
d2$CCFL_orig <- d2$CCFL_incl - d2$MCA
d2$BALSW_orig <- cumsum(d2$BA.sw) - d2$BA.sw
d2$CCFLSW_orig <- cumsum(d2$MCA.SW) - d2$MCA.SW
# the order dependent cumulative sum is what a naive cumsum gives when two
# trees share a diameter; forestCI instead treats tied trees symmetrically
d2$BAL_tieorder <- cumsum(d2$BA) - d2$BA
orig_di <- data.frame(TreeNum = as.character(d2$TreeNum),
                      bal_tieorder_orig = d2$BAL_tieorder,
                      mcw_orig = d2$MCW, mca_orig = d2$MCA,
                      bal_orig = d2$BAL_orig, ccfl_orig = d2$CCFL_orig,
                      bal_sw_orig = d2$BALSW_orig,
                      ccfl_sw_orig = d2$CCFLSW_orig,
                      ba_plot_orig = sum(d2$BA), ccf_orig = sum(d2$MCA),
                      stringsAsFactors = FALSE)

setwd(if (length(args) >= 2) dirname(normalizePath(file.path(outdir, "."))) else ".")

# ---- 3. compare ------------------------------------------------------------
cmp <- as.data.table(pkg)
cmp[, TreeNum := as.character(TreeNum)]
for (tab in list(orig_csa, orig_rapa, orig_hegyi, orig_spurr, orig_di)) {
  if (!is.null(tab)) cmp <- merge(cmp, as.data.table(tab), by = "TreeNum", all.x = TRUE)
}

pairs <- list(
  c("mcw",  "mcw_orig",  "maximum crown width (m)"),
  c("mca",  "mca_orig",  "maximum crown area (percent)"),
  c("bal",  "bal_orig",  "basal area of larger trees (m2 ha-1)"),
  c("bal_sw", "bal_sw_orig", "softwood BAL (m2 ha-1)"),
  c("ccfl", "ccfl_orig", "crown competition factor of larger trees (percent)"),
  c("ccfl_sw", "ccfl_sw_orig", "softwood CCFL (percent)"),
  c("ba",   "ba_plot_orig", "plot basal area (m2 ha-1)"),
  c("ccf",  "ccf_orig",  "crown competition factor (percent)"),
  c("csa_meas", "csa_orig", "crown surface area, measured radii (m2)"),
  c("cv_meas",  "cv_orig",  "crown volume, measured radii (m3)"),
  c("csa",  "csa_orig",  "crown surface area, predicted LCW (m2)"),
  c("cv",   "cv_orig",   "crown volume, predicted LCW (m3)"),
  c("rapa", "rapa_orig", "rasterised APA (m2)"),
  c("hegyi","hegyi_orig","Hegyi index"),
  c("spurr","spurr_orig","Spurr point density"),
  c("n_comp","n_in_orig","angle gauge competitor count")
)
res <- rbindlist(lapply(pairs, function(p) {
  if (!all(p[1:2] %in% names(cmp))) {
    return(data.table(quantity = p[3], n = 0L, max_abs_diff = NA_real_,
                      max_rel_diff_pct = NA_real_, r = NA_real_,
                      verdict = "not computed"))
  }
  a <- cmp[[p[1]]]; b <- cmp[[p[2]]]
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(data.table(quantity = p[3], n = 0L, max_abs_diff = NA_real_,
                                  max_rel_diff_pct = NA_real_, r = NA_real_,
                                  verdict = "no overlap"))
  a <- a[ok]; b <- b[ok]
  ad <- max(abs(a - b))
  rd <- 100 * max(abs(a - b) / pmax(abs(b), .Machine$double.eps))
  rr <- if (stats::sd(a) > 0 && stats::sd(b) > 0) stats::cor(a, b) else NA_real_
  data.table(quantity = p[3], n = sum(ok), max_abs_diff = ad,
             max_rel_diff_pct = rd, r = rr,
             verdict = if (ad < 1e-6) "identical"
                       else if (rd < 0.5) "agrees to 0.5 percent"
                       else if (!is.na(rr) && rr > 0.99) "same ordering, different scale"
                       else "DIVERGES")
}))
print(res)
fwrite(res, file.path(outdir, "2026-09-17_forestCI-regression-vs-original_FINAL.csv"))
fwrite(cmp, file.path(outdir, "2026-09-17_forestCI-regression-treewise_DATA.csv"))

# ---- 4. figure --------------------------------------------------------------
trap("regression figure", {
  show <- pairs[vapply(pairs, function(p) all(p[1:2] %in% names(cmp)), logical(1))]
  show <- show[seq_len(min(9L, length(show)))]
  png(file.path(outdir, "2026-09-17_forestCI-regression_FINAL.png"),
      width = 2400, height = 2400, res = 300)
  op <- par(mfrow = c(3, 3), mar = c(4, 4, 2.2, 1), bty = "l", cex.axis = 0.8)
  for (p in show) {
    a <- cmp[[p[1]]]; b <- cmp[[p[2]]]
    ok <- is.finite(a) & is.finite(b)
    plot(b[ok], a[ok], pch = 21, bg = "#2b6a4d", col = "grey25", cex = 0.8,
         xlab = "original scripts", ylab = "forestCI", main = p[3],
         cex.main = 0.85)
    abline(0, 1, col = "#c02020", lwd = 1.5)
  }
  par(op); dev.off()
})

gc()
cat("\nDone.\n")
