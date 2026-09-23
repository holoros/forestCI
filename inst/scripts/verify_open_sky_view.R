# Open sky view against CI.R POVall(fly = TRUE)
#
# The comparator column is potencover, which counts only rays with potential to
# face sky, not covertypes, which includes downward facing rays scored as
# floor. Result on the Penobscot example plot: forestCI reads a mean sky
# fraction of 0.5465 against 0.4814, a 13.5 percent high bias in the expected
# direction, correlating at 0.834 on the ratio and 0.956 on the area.
# Convergence, same plot: res1 4, 8, 16, 32 at res2 16 give original means of
# 0.4509, 0.4814, 0.4727 and 0.4760, and res2 32 gives 0.4732 at res1 8 and
# 0.4705 at res1 16, so above the default the reference stays within 0.011
# and forestCI's bias is 13.5 to 16.1 percent. The convergence is not
# monotone, so no Richardson limit is reported. res1 = 32 takes 835 s.
#
# Usage: Rscript verify_open_sky_view.R [orig_dir] [out_dir] [res1] [res2]
#   orig_dir  directory holding CI.R, SPP.R, circular2xy.R,
#             sort.data.frame.r and PEF.csv
#   out_dir   where the CSVs are written, default the working directory
#   res1      POVall facet resolution on the subject crown, default 8, which
#             gives 32 facets per tree; a comma list runs a convergence study
#   res2      POVall ray resolution of the hemisphere cast from each facet,
#             default 16; a comma list of the same length as res1 is paired
#             with it element by element, any other length is crossed
# Each resolution pair writes its own treewise and summary CSV, tagged -resA-B
# unless it is the default 8, 16, and a pair whose treewise CSV already exists
# in out_dir is read back rather than recomputed. That lets the pairs run as
# separate SLURM jobs and a final call with the full lists only summarize.
# With more than one pair a convergence table is written as well.
# Headless. No graphics device is opened.

args <- commandArgs(trailingOnly = TRUE)
orig <- if (length(args) >= 1) args[1] else "orig"
outd <- if (length(args) >= 2) args[2] else "."
res1s <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else 8L
res2s <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else 16L
dir.create(outd, showWarnings = FALSE, recursive = TRUE)
outd <- normalizePath(outd); orig <- normalizePath(orig)
elog <- file.path(outd, "error_log.txt")
pairs <- if (length(res1s) == length(res2s)) data.frame(r1 = res1s, r2 = res2s) else
  expand.grid(r1 = res1s, r2 = res2s)
set.seed(20260923)
suppressPackageStartupMessages({ library(forestCI); library(data.table) })
raw <- read.csv(file.path(orig, "PEF.csv"), stringsAsFactors = FALSE)
pid <- as.character(raw$PlotID.yr[1]); d1 <- raw[raw$PlotID.yr == pid, ]
st <- as_stand(d1, plot="PlotID.yr", tree="TreeNum", species="FVS", dbh="DBH",
               height="HT", hcb="HCB", distance="Dist", azimuth="Az", expf="EXPF")
# match the crown radius input to CI.R, which builds it from the four measured
# radii (CI.R 194, 196) rather than from predicted LCW
cpa_meas <- (pi/4) * (d1$Rad0*d1$Rad90 + d1$Rad90*d1$Rad180 +
                      d1$Rad180*d1$Rad270 + d1$Rad270*d1$Rad0)
mm <- data.table(tree_id = as.character(d1$TreeNum),
                 lcw_meas = 2 * sqrt(cpa_meas / pi))
tr <- merge(as.data.table(st$trees), mm, by = "tree_id", sort = FALSE)
tr[, lcw := lcw_meas]; st$trees <- tr
t0 <- proc.time()[3]
ce <- as.data.table(crown_exposure(st, osv = TRUE))
cat("forestCI crown_exposure with osv:", round(proc.time()[3]-t0,1), "s\n")
pkg <- ce[, .(TreeNum = as.integer(tree_id), csa, csax, csax_rel, osv, osv_rel)]

pd <- function(a,b){k <- is.finite(a)&is.finite(b)&abs(b)>1e-9
  if(!any(k)) return(NA_real_); max(abs(a[k]-b[k])/abs(b[k]))*100}
tagof <- function(r1, r2) if (r1 == 8L && r2 == 16L) "" else sprintf("-res%d-%d", r1, r2)

setwd(orig)
source("CI.R"); source("SPP.R"); source("circular2xy.R"); source("sort.data.frame.r")
run_pair <- function(r1, r2) {
  tag <- tagof(r1, r2)
  f_tw <- file.path(outd, sprintf("2026-09-23_forestCI-osv-matched%s-treewise_DATA.csv", tag))
  if (file.exists(f_tw) && file.size(f_tw) > 0) {
    cat("res1", r1, "res2", r2, ": read back from", basename(f_tw), "\n")
    return(fread(f_tw))
  }
  data <- read.csv("PEF.csv", stringsAsFactors = TRUE)
  FORMAT(inputs = data, plots = "PlotID.yr", trees = "TreeNum",
         positions = "none", aspects = "none", distances = "Dist",
         directions = "Az", slopes = "none", spps = "FVS", dbhs = "DBH",
         cohorts = "none", cnpypstns = "none", hgts = "HT", crown_bases = "HCB",
         north_crs = "Rad0", east_crs = "Rad90", south_crs = "Rad180",
         west_crs = "Rad270", measures = "none")
  t0 <- proc.time()[3]
  ok <- tryCatch({
    POVall(set = pid, widept = wides, upexp = ups, loexp = los, shape = shps,
           trnc = FALSE, fly = TRUE, res1 = r1, res2 = r2,
           draw = FALSE, mercator = FALSE, spherical = FALSE, solid = TRUE)
    TRUE }, error = function(e) {
    cat("[POVall res", r1, r2, "]", conditionMessage(e), "\n", file = elog, append = TRUE)
    cat("POVall FAILED:", conditionMessage(e), "\n"); FALSE })
  el <- round(proc.time()[3]-t0, 1)
  cat("POVall fly res1 =", r1, "res2 =", r2, ", 25 trees:", el, "s\n")
  if (!ok) return(NULL)
  pc <- as.data.table(potencover); ct <- as.data.table(covertypes)
  o <- data.table(TreeNum = as.integer(pc$TREE), sky_poten = as.numeric(pc$SKY),
                  tot_poten = as.numeric(pc$TOTAL), sky_cover = as.numeric(ct$SKY),
                  tot_cover = as.numeric(ct$TOTAL))
  m <- merge(pkg, o, by = "TreeNum")
  m[, `:=`(osv_rel_orig_poten = sky_poten / tot_poten,
           osv_rel_orig_cover = sky_cover / tot_cover,
           res1 = r1, res2 = r2, elapsed_s = el)]
  fwrite(m, f_tw)
  s <- data.table(
    quantity = c("open sky view ratio, pkg vs orig potencover",
                 "open sky view ratio, pkg vs orig covertypes",
                 "open sky view area m2, pkg vs orig potencover",
                 "orig TOTAL vs analytic CSA"),
    max_pct_diff = c(pd(m$osv_rel, m$osv_rel_orig_poten),
                     pd(m$osv_rel, m$osv_rel_orig_cover),
                     pd(m$osv, m$sky_poten), pd(m$tot_poten, m$csa)),
    correlation = c(cor(m$osv_rel, m$osv_rel_orig_poten, use="complete"),
                    cor(m$osv_rel, m$osv_rel_orig_cover, use="complete"),
                    cor(m$osv, m$sky_poten, use="complete"),
                    cor(m$tot_poten, m$csa, use="complete")))
  fwrite(s, file.path(outd, sprintf("2026-09-23_forestCI-osv-matched%s-summary_DATA.csv", tag)))
  print(s)
  cat("\nmean osv_rel pkg", round(mean(m$osv_rel, na.rm=TRUE),4),
      " orig poten", round(mean(m$osv_rel_orig_poten, na.rm=TRUE),4),
      " orig cover", round(mean(m$osv_rel_orig_cover, na.rm=TRUE),4), "\n")
  m
}
runs <- Filter(Negate(is.null), Map(run_pair, pairs$r1, pairs$r2))

# convergence of the reference in its own resolution, against the finest pair
if (length(runs) > 1L) {
  allr <- rbindlist(runs, fill = TRUE)
  fin <- allr[res1 == max(res1)][res2 == max(res2)]
  ref <- fin[, .(TreeNum, ref = osv_rel_orig_poten)]
  cv <- merge(allr, ref, by = "TreeNum")[, .(
    mean_orig = mean(osv_rel_orig_poten, na.rm = TRUE),
    mean_pkg = mean(osv_rel, na.rm = TRUE),
    pkg_bias_pct = (mean(osv_rel, na.rm = TRUE) /
                    mean(osv_rel_orig_poten, na.rm = TRUE) - 1) * 100,
    cor_pkg = cor(osv_rel, osv_rel_orig_poten, use = "complete"),
    max_abs_vs_finest = max(abs(osv_rel_orig_poten - ref), na.rm = TRUE),
    max_pct_vs_finest = pd(osv_rel_orig_poten, ref),
    cor_vs_finest = cor(osv_rel_orig_poten, ref, use = "complete"),
    tot_vs_csa = mean(tot_poten / csa, na.rm = TRUE),
    elapsed_s = elapsed_s[1]), by = .(res1, res2)]
  setorder(cv, res2, res1)
  # Richardson estimate of the res1 limit of the stand mean, where three res1
  # values in ratio 2 share a res2 and the differences shrink monotonically
  cv[, `:=`(order_p = NA_real_, mean_orig_limit = NA_real_)]
  for (r2 in unique(cv$res2)) {
    k <- cv[res2 == r2]
    if (nrow(k) < 3L) next
    k <- tail(k, 3L); f <- k$mean_orig
    if (all(diff(k$res1 / head(k$res1, 1)) > 0) && all(k$res1[-1] / head(k$res1, -1) == 2)) {
      q <- (f[1] - f[2]) / (f[2] - f[3])
      if (is.finite(q) && q > 1) {
        p <- log2(q)
        cv[res2 == r2 & res1 == max(k$res1),
           `:=`(order_p = p, mean_orig_limit = f[3] + (f[3] - f[2]) / (2^p - 1))]
      }
    }
  }
  fwrite(cv, file.path(outd, "2026-09-23_forestCI-osv-convergence_DATA.csv"))
  print(cv)
}
gc(); cat("=== DONE ===\n")
