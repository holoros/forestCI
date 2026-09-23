# Open sky view against CI.R POVall(fly = TRUE)
#
# The comparator column is potencover, which counts only rays with potential to
# face sky, not covertypes, which includes downward facing rays scored as
# floor. Result on the Penobscot example plot: forestCI reads a mean sky
# fraction of 0.5465 against 0.4814, a 13.5 percent high bias in the expected
# direction, correlating at 0.834 on the ratio and 0.956 on the area.
#
# Usage: Rscript verify_open_sky_view.R [orig_dir] [out_dir]
#   orig_dir  directory holding CI.R, SPP.R, circular2xy.R,
#             sort.data.frame.r and PEF.csv
#   out_dir   where the CSVs are written, default the working directory
# Headless. No graphics device is opened.

# Open sky view: forestCI crown_exposure() against CI.R POVall(fly = TRUE).
args <- commandArgs(trailingOnly = TRUE)
orig <- if (length(args) >= 1) args[1] else "orig"
outd <- if (length(args) >= 2) args[2] else "."
dir.create(outd, showWarnings = FALSE, recursive = TRUE)
elog <- file.path(outd, "error_log.txt")
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

setwd(orig)
source("CI.R"); source("SPP.R"); source("circular2xy.R"); source("sort.data.frame.r")
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
         trnc = FALSE, fly = TRUE, res1 = 8L, res2 = 16L,
         draw = FALSE, mercator = FALSE, spherical = FALSE, solid = TRUE)
  TRUE }, error = function(e) {
  cat("[POVall]", conditionMessage(e), "\n", file = elog, append = TRUE)
  cat("POVall FAILED:", conditionMessage(e), "\n"); FALSE })
cat("POVall fly res1=8, 25 trees:", round(proc.time()[3]-t0,1), "s\n")
if (ok) {
  pc <- as.data.table(potencover); ct <- as.data.table(covertypes)
  o <- data.table(TreeNum = as.integer(pc$TREE), sky_poten = as.numeric(pc$SKY),
                  tot_poten = as.numeric(pc$TOTAL), sky_cover = as.numeric(ct$SKY),
                  tot_cover = as.numeric(ct$TOTAL))
  m <- merge(pkg, o, by = "TreeNum")
  m[, `:=`(osv_rel_orig_poten = sky_poten / tot_poten,
           osv_rel_orig_cover = sky_cover / tot_cover)]
  fwrite(m, file.path(outd, "2026-09-23_forestCI-osv-matched-treewise_DATA.csv"))
  pd <- function(a,b){k <- is.finite(a)&is.finite(b)&abs(b)>1e-9
    if(!any(k)) return(NA_real_); max(abs(a[k]-b[k])/abs(b[k]))*100}
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
  fwrite(s, file.path(outd, "2026-09-23_forestCI-osv-matched-summary_DATA.csv"))
  print(s)
  cat("\nmean osv_rel pkg", round(mean(m$osv_rel, na.rm=TRUE),4),
      " orig poten", round(mean(m$osv_rel_orig_poten, na.rm=TRUE),4),
      " orig cover", round(mean(m$osv_rel_orig_cover, na.rm=TRUE),4), "\n")
}
gc(); cat("=== DONE ===\n")
