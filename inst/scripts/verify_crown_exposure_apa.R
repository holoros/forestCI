# Regression of forestCI crown exposure and polygon APA against CI.R.
# Writes CSVs to the directory given as arg 2. No graphics, no interactive device.
#
# Polygon APA is compared unclipped (boundary = "none") and classified with
# apa()'s apa_bounded flag. On the Penobscot example plot the 15 bounded trees
# agree to a relative difference below 1e-10. The 10 unbounded trees are the
# trees on the stem convex hull, and the original has no answer for any of
# them: ENGINE's elimination step assumes the neighbours surround the subject,
# so across an angular gap of more than a half turn it discards genuine
# bounding neighbours. With two neighbours left the two corners coincide and
# the area is exactly 0 (6 trees), with one it is NA (2 trees), and with three
# or more it closes a polygon that does not contain the subject stem (2 trees,
# TreeNum 15000 and 24600). The orig_poly_has_stem column tests the last case
# directly on ENGINE's own corner coordinates.
args <- commandArgs(trailingOnly = TRUE)
orig <- if (length(args) >= 1) args[1] else "orig"
outd <- if (length(args) >= 2) args[2] else "."
dir.create(outd, showWarnings = FALSE, recursive = TRUE)
outd <- normalizePath(outd); orig <- normalizePath(orig)
elog <- file.path(outd, "error_log.txt")
trap <- function(tag, expr) tryCatch(expr, error = function(e) {
  cat(sprintf("[%s] %s\n", tag, conditionMessage(e)), file = elog, append = TRUE)
  cat("FAILED:", tag, ":", conditionMessage(e), "\n"); NULL })
set.seed(20260923)
suppressPackageStartupMessages({ library(forestCI); library(data.table) })

## ---- 1. forestCI side ------------------------------------------------------
raw <- read.csv(file.path(orig, "PEF.csv"), stringsAsFactors = FALSE)
pid <- as.character(raw$PlotID.yr[1])
d1  <- raw[raw$PlotID.yr == pid, ]
st <- as_stand(d1, plot = "PlotID.yr", tree = "TreeNum", species = "FVS",
               dbh = "DBH", height = "HT", hcb = "HCB",
               distance = "Dist", azimuth = "Az", expf = "EXPF")
ce <- crown_exposure(st, osv = FALSE)
ap <- apa(st, weight = "none")
ap_open <- as.data.table(apa(st, weight = "none", boundary = "none"))[
  , .(tree_id, apa_unclipped = apa, apa_bounded)]
pkg <- merge(as.data.table(ce)[, .(tree_id, csa, csax, csax_rel, cpa, cpax, cpax_rel)],
             as.data.table(ap)[, .(tree_id, apa)], by = "tree_id")
pkg <- merge(pkg, ap_open, by = "tree_id")
setnames(pkg, "tree_id", "TreeNum")
cat("forestCI side:", nrow(pkg), "trees\n")

## ---- 2. original side ------------------------------------------------------
setwd(orig)
source("CI.R"); source("SPP.R"); source("circular2xy.R"); source("sort.data.frame.r")
data <- read.csv("PEF.csv", stringsAsFactors = TRUE)
FORMAT(inputs = data, plots = "PlotID.yr", trees = "TreeNum",
       positions = "none", aspects = "none", distances = "Dist",
       directions = "Az", slopes = "none", spps = "FVS", dbhs = "DBH",
       cohorts = "none", cnpypstns = "none", hgts = "HT", crown_bases = "HCB",
       north_crs = "Rad0", east_crs = "Rad90", south_crs = "Rad180",
       west_crs = "Rad270", measures = "none")
CPARS <- list(widept = wides, upexp = ups, loexp = los, shape = shps, trnc = FALSE)

# 2a. polygon APA. APAr wraps START, CROWNS, FLTnWGTdfs, ENGINE.
# even odd point in polygon test on ENGINE's absolute corner coordinates
stem_in_corners <- function(k) {
  if (k > length(CornersXk) || k > length(CornersYk)) return(NA)
  cx <- CornersXk[[k]]; cy <- CornersYk[[k]]
  if (is.null(cx) || length(cx) < 3L || any(!is.finite(c(cx, cy)))) return(NA)
  px <- X[k]; py <- Y[k]; L <- length(cx); inside <- FALSE
  for (q in seq_len(L)) {
    r <- if (q == L) 1L else q + 1L
    if ((cy[q] > py) != (cy[r] > py) &&
        px < (cx[r] - cx[q]) * (py - cy[q]) / (cy[r] - cy[q]) + cx[q])
      inside <- !inside
  }
  inside
}
orig_apa <- trap("APAr", {
  do.call(APAr, c(list(choice = pid, crowncenters = FALSE), CPARS))
  data.frame(TreeNum = as.character(names(APAk)), apa_orig = as.numeric(APAk),
             orig_poly_has_stem = vapply(seq_along(APAk), stem_in_corners, NA),
             stringsAsFactors = FALSE)
})

# 2b. crown exposure. OVERLAP hard codes crowncenters = TRUE at CI.R 3607, so the
# chain is called directly instead, once stem centred to match forestCI and once
# crown centred to measure what that choice is worth.
run_motor <- function(cc, res, rad) {
  do.call(START, list(choice = pid, crowncenters = cc, xy = FALSE))
  do.call(CROWNS, CPARS)
  BOUNDS(FALSE, 25, res)
  MOTOR("NONE", FALSE, FALSE, 24, 16, 8, rad)
  data.frame(TreeNum = as.character(N),
             csax_orig = as.numeric(CSAx), cpax_orig = as.numeric(CPAx),
             upcsa_orig = as.numeric(upCSA), locsa_orig = as.numeric(loCSA),
             csa_orig = as.numeric(CSA), cpa_orig = as.numeric(CPA),
             stringsAsFactors = FALSE)
}
t0 <- proc.time()[3]
orig_ce_stem <- trap("MOTOR stem centred", run_motor(FALSE, 0.25, 16.05))
cat("MOTOR stem centred:", round(proc.time()[3] - t0, 1), "s\n")
t0 <- proc.time()[3]
orig_ce_crwn <- trap("MOTOR crown centred", run_motor(TRUE, 0.25, 16.05))
cat("MOTOR crown centred:", round(proc.time()[3] - t0, 1), "s\n")
# grid sensitivity: the rasterised exposure ratio must be stable in the cell size
orig_ce_fine <- trap("MOTOR stem centred res 0.125", run_motor(FALSE, 0.125, 16.05))

## ---- 3. compare ------------------------------------------------------------
cmp <- as.data.table(pkg)
addc <- function(dt, x, sfx) {
  if (is.null(x)) return(dt)
  x <- as.data.table(x); setnames(x, setdiff(names(x), "TreeNum"),
                                  paste0(setdiff(names(x), "TreeNum"), sfx))
  merge(dt, x, by = "TreeNum", all.x = TRUE)
}
cmp <- addc(cmp, orig_apa, "")
cmp <- addc(cmp, orig_ce_stem, "_s")
cmp <- addc(cmp, orig_ce_crwn, "_c")
cmp <- addc(cmp, orig_ce_fine, "_f")
# the original normalises exposed surface by the UPPER crown only (CI.R 3342),
# so its exposure ratio has a different denominator than csax_rel
if ("csax_orig_s" %in% names(cmp))
  cmp[, `:=`(ratio_orig_s = csax_orig_s / upcsa_orig_s,
             ratio_orig_f = csax_orig_f / upcsa_orig_f,
             ratio_orig_c = csax_orig_c / upcsa_orig_c,
             cpa_ratio_orig_s = cpax_orig_s / cpa_orig_s)]
if ("apa_orig" %in% names(cmp)) {
  cmp[, apa_rel_diff := abs(apa_unclipped - apa_orig) / pmax(abs(apa_orig), 1e-12)]
  cmp[, apa_class := fifelse(apa_bounded,
    fifelse(!is.na(apa_rel_diff) & apa_rel_diff < 1e-10, "bounded, exact",
            "bounded, DIFFERS"),
    fifelse(is.na(apa_orig), "unbounded, original NA",
    fifelse(apa_orig == 0, "unbounded, original zero",
    fifelse(orig_poly_has_stem %in% FALSE,
            "unbounded, original polygon excludes the stem",
            "unbounded, other"))))]
  cls <- cmp[, .N, by = apa_class][order(apa_class)]
  fwrite(cls, file.path(outd, "2026-09-23_forestCI-apa-class_DATA.csv"))
  print(cls)
}
fwrite(cmp, file.path(outd, "2026-09-23_forestCI-ce-apa-treewise_DATA.csv"))

pd <- function(a, b) { k <- is.finite(a) & is.finite(b) & (abs(a) + abs(b)) > 0
  if (!any(k)) return(NA_real_); max(abs(a[k] - b[k]) / pmax(abs(b[k]), 1e-12)) * 100 }
bnd <- if ("apa_orig" %in% names(cmp)) cmp[apa_bounded == TRUE] else cmp[0]
res <- data.table(
  quantity = c("polygon APA, unweighted",
               "polygon APA, unweighted, unclipped, bounded cells",
               "analytic crown surface area",
               "analytic crown projection area",
               "exposed CSA ratio, pkg vs orig stem centred",
               "exposed CSA ratio, orig res 0.25 vs 0.125",
               "exposed CSA ratio, orig stem vs crown centred",
               "exposed CPA ratio, pkg vs orig stem centred"),
  max_pct_diff = c(
    if (!is.null(orig_apa)) pd(cmp$apa, cmp$apa_orig) else NA_real_,
    if (nrow(bnd)) pd(bnd$apa_unclipped, bnd$apa_orig) else NA_real_,
    pd(cmp$csa, cmp$csa_orig_s), pd(cmp$cpa, cmp$cpa_orig_s),
    pd(cmp$csax_rel, cmp$ratio_orig_s),
    pd(cmp$ratio_orig_s, cmp$ratio_orig_f),
    pd(cmp$ratio_orig_s, cmp$ratio_orig_c),
    pd(cmp$cpax_rel, cmp$cpa_ratio_orig_s)))
res[, correlation := c(
  if (!is.null(orig_apa)) cor(cmp$apa, cmp$apa_orig, use = "complete") else NA_real_,
  if (nrow(bnd) > 2L) cor(bnd$apa_unclipped, bnd$apa_orig, use = "complete") else NA_real_,
  cor(cmp$csa, cmp$csa_orig_s, use = "complete"),
  cor(cmp$cpa, cmp$cpa_orig_s, use = "complete"),
  cor(cmp$csax_rel, cmp$ratio_orig_s, use = "complete"),
  cor(cmp$ratio_orig_s, cmp$ratio_orig_f, use = "complete"),
  cor(cmp$ratio_orig_s, cmp$ratio_orig_c, use = "complete"),
  cor(cmp$cpax_rel, cmp$cpa_ratio_orig_s, use = "complete"))]
fwrite(res, file.path(outd, "2026-09-23_forestCI-ce-apa-summary_DATA.csv"))
print(res)
cat("\nmeans: pkg csax_rel", round(mean(cmp$csax_rel, na.rm = TRUE), 4),
    " orig ratio", round(mean(cmp$ratio_orig_s, na.rm = TRUE), 4), "\n")
cat("means: pkg cpax_rel", round(mean(cmp$cpax_rel, na.rm = TRUE), 4),
    " orig ratio", round(mean(cmp$cpa_ratio_orig_s, na.rm = TRUE), 4), "\n")
gc()
cat("=== DONE ===\n")
