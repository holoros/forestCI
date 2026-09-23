# Exposed crown surface and projection against CI.R, crown radius matched
#
# CI.R builds crown radius from the four measured radii while crown_exposure()
# reads predicted largest crown width. Substituting the measured equivalent
# isolates the exposure algorithms. Result on the Penobscot example plot:
# analytic crown projection area exact, crown surface area 0.264 percent,
# exposed projection ratio r = 0.996, exposed surface ratio r = 0.984.
# Run verify_crown_exposure_apa.R first; this script reads its treewise CSV.
#
# Usage: Rscript verify_crown_exposure_matched.R [orig_dir] [out_dir]
#   orig_dir  directory holding CI.R, SPP.R, circular2xy.R,
#             sort.data.frame.r and PEF.csv
#   out_dir   where the CSVs are written, default the working directory
# Headless. No graphics device is opened.

# Exposed crown surface and projection, with the crown RADIUS input matched.
# CI.R builds crown radius from the four measured radii (CPA at 194, CR at 196);
# forestCI's crown_exposure reads lcw, which is predicted. Overwriting lcw with
# the measured equivalent removes the one remaining input difference, so what is
# left is a comparison of the two exposure algorithms and nothing else.
args <- commandArgs(trailingOnly = TRUE)
orig <- if (length(args) >= 1) args[1] else "orig"
outd <- if (length(args) >= 2) args[2] else "."
dir.create(outd, showWarnings = FALSE, recursive = TRUE)
set.seed(20260923)
suppressPackageStartupMessages({ library(forestCI); library(data.table) })
raw <- read.csv(file.path(orig, "PEF.csv"), stringsAsFactors = FALSE)
pid <- as.character(raw$PlotID.yr[1]); d1 <- raw[raw$PlotID.yr == pid, ]
st <- as_stand(d1, plot="PlotID.yr", tree="TreeNum", species="FVS", dbh="DBH",
               height="HT", hcb="HCB", distance="Dist", azimuth="Az", expf="EXPF")
# CI.R 194: four quadrant ellipse sector area; 196: its mean circle radius
cpa_meas <- (pi/4) * (d1$Rad0*d1$Rad90 + d1$Rad90*d1$Rad180 +
                      d1$Rad180*d1$Rad270 + d1$Rad270*d1$Rad0)
cr_meas  <- sqrt(cpa_meas / pi)
mm <- data.table(tree_id = as.character(d1$TreeNum), lcw_meas = 2 * cr_meas)
tr <- as.data.table(st$trees)
tr <- merge(tr, mm, by = "tree_id", sort = FALSE)
tr[, lcw := lcw_meas]
st$trees <- tr
ce <- as.data.table(crown_exposure(st, osv = FALSE))
o  <- fread(file.path(outd, "2026-09-23_forestCI-ce-apa-treewise_DATA.csv"))
m  <- merge(ce[, .(TreeNum = as.integer(tree_id), csa, csax, csax_rel,
                   cpa, cpax, cpax_rel)], o[, .(TreeNum, csa_orig_s, cpa_orig_s,
                   csax_orig_s, upcsa_orig_s, ratio_orig_s, cpa_ratio_orig_s,
                   ratio_orig_f)], by = "TreeNum")
fwrite(m, file.path(outd, "2026-09-23_forestCI-ce-matched-treewise_DATA.csv"))
pd <- function(a,b){k <- is.finite(a)&is.finite(b)&abs(b)>1e-9
  if(!any(k)) return(NA_real_); max(abs(a[k]-b[k])/abs(b[k]))*100}
s <- data.table(
  quantity = c("analytic crown surface area, matched radius",
               "analytic crown projection area, matched radius",
               "exposed crown surface ratio",
               "exposed crown projection ratio"),
  max_pct_diff = c(pd(m$csa, m$csa_orig_s), pd(m$cpa, m$cpa_orig_s),
                   pd(m$csax_rel, m$ratio_orig_s), pd(m$cpax_rel, m$cpa_ratio_orig_s)),
  correlation = c(cor(m$csa, m$csa_orig_s), cor(m$cpa, m$cpa_orig_s),
                  cor(m$csax_rel, m$ratio_orig_s), cor(m$cpax_rel, m$cpa_ratio_orig_s)))
fwrite(s, file.path(outd, "2026-09-23_forestCI-ce-matched-summary_DATA.csv"))
print(s)
cat("\nmean exposed CSA ratio  pkg", round(mean(m$csax_rel),4),
    " orig", round(mean(m$ratio_orig_s),4), "\n")
cat("mean exposed CPA ratio  pkg", round(mean(m$cpax_rel),4),
    " orig", round(mean(m$cpa_ratio_orig_s),4), "\n")
cat("trees whose exposed CSA ratio differs by more than 0.05:",
    sum(abs(m$csax_rel - m$ratio_orig_s) > 0.05), "of", nrow(m), "\n")
cat("=== DONE ===\n")
