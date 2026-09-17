# ----------------------------------------------------------------------------
# forestCI demonstration and internal validation on the Penobscot example plots
#
# Runs the full index set on the two stem-mapped plots shipped with the
# package, checks the identities that must hold, writes a summary table and
# two headless figures.
#
#   Rscript inst/scripts/pef_competition_demo.R [outdir]
#
# Headless throughout: png() and dev.off(), never a screen device. Every read,
# fit and plot is wrapped so a failure appends to error_log.txt rather than
# killing a batch job.
# ----------------------------------------------------------------------------

set.seed(20260917)
outdir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
elog <- file.path(outdir, "error_log.txt")

trap <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    cat(sprintf("[%s] %s: %s\n", format(Sys.time()), label,
                conditionMessage(e)), file = elog, append = TRUE)
    NULL
  })
}

suppressPackageStartupMessages({
  library(forestCI)
  library(data.table)
})

st <- trap("build stand", pef_stand())
stopifnot(!is.null(st))

ci <- trap("competition indices", competition_indices(
  st, radius = 6, edge = "flag", apa_weight = "dbh", apa_exponent = 2,
  rapa_exponent = 1, rapa_resolution = 0.25, n_depth = 12L, n_azimuth = 36L,
  crown_osv = TRUE, quiet = TRUE))
stopifnot(!is.null(ci))

# --- identity checks that must hold regardless of the data ------------------
checks <- data.table(check = character(), result = character())
add <- function(name, ok, detail = "") {
  checks <<- rbind(checks, data.table(
    check = name, result = paste0(if (ok) "pass" else "FAIL",
                                  if (nzchar(detail)) paste0(" (", detail, ")") else "")))
}

add("BAL of the largest tree on each plot is zero",
    all(abs(ci[, .(m = min(bal)), by = plot_id]$m) < 1e-8))
add("BAL never exceeds plot basal area", all(ci$bal <= ci$ba + 1e-8))
add("softwood and hardwood BAL partitions sum to BAL",
    max(abs(ci$bal - ci$bal_sw - ci$bal_hw)) < 1e-8)
add("softwood and hardwood CCFL partitions sum to CCFL",
    max(abs(ci$ccfl - ci$ccfl_sw - ci$ccfl_hw)) < 1e-8)
# APA tiles the plot exactly only when unweighted: displaced bisectors are an
# approximation, not a partition, so the weighted version is checked for
# plausibility rather than for exact tiling.
apa0 <- apa(st, weight = "none")
apa_tot <- apa0[, .(s = sum(apa)), by = plot_id]
plot_m2 <- st$plots$plot_area * 10000
add("unweighted APA polygons tile the plot",
    max(abs(apa_tot$s - plot_m2) / plot_m2) < 0.01,
    sprintf("max error %.2f percent",
            100 * max(abs(apa_tot$s - plot_m2) / plot_m2)))
# Diagnostic, not a pass or fail: displaced bisectors are an approximation to a
# weighted tessellation, so the weighted areas are not expected to conserve.
apaw_tot <- ci[, .(s = sum(apa)), by = plot_id]
checks <- rbind(checks, data.table(
  check = "weighted APA total departure from plot area (diagnostic)",
  result = sprintf("%.1f percent",
                   100 * max(abs(apaw_tot$s - plot_m2) / plot_m2))))
add("exposed crown surface area never exceeds crown surface area",
    all(ci$csax <= ci$csa + 1e-6))
add("open sky view never exceeds exposed crown surface area",
    all(ci$osv <= ci$csax + 1e-6))
# apa() and rapa() weight by different rules (bounded bisector displacement
# versus multiplicative Voronoi), so they are comparable only unweighted.
rapa0 <- rapa(st, weight = "none", resolution = 0.25)
cmp <- merge(apa0, rapa0, by = c("plot_id", "tree_id"))
add("unweighted rasterised and polygon APA agree",
    stats::cor(cmp$apa, cmp$rapa) > 0.95,
    sprintf("r = %.3f", stats::cor(cmp$apa, cmp$rapa)))
add("plot basal area is on a per hectare scale",
    all(ci$ba > 5 & ci$ba < 100),
    sprintf("range %.1f to %.1f m2 ha-1", min(ci$ba), max(ci$ba)))
add("local basal area is on the same scale as stand basal area",
    abs(mean(ci$local_ba) - mean(ci$ba)) / mean(ci$ba) < 0.6,
    sprintf("local %.1f vs stand %.1f m2 ha-1",
            mean(ci$local_ba), mean(ci$ba)))

fwrite(checks, file.path(outdir, "2026-09-17_forestCI-identity-checks_FINAL.csv"))
print(checks)

# --- summary table ----------------------------------------------------------
idx <- c("ba", "bal", "ccf", "ccfl", "sdi", "rd", "hegyi", "martin_ek",
         "lorimer", "spurr", "local_ba", "local_bal", "mean_dist",
         "apa", "rapa", "csa", "csax_rel", "osv_rel")
idx <- intersect(idx, names(ci))
summ <- rbindlist(lapply(idx, function(v) {
  x <- ci[[v]]
  data.table(index = v, n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE),
             sd = stats::sd(x, na.rm = TRUE), min = min(x, na.rm = TRUE),
             max = max(x, na.rm = TRUE))
}))
fwrite(summ, file.path(outdir, "2026-09-17_forestCI-pef-index-summary_FINAL.csv"))
print(summ)

# --- Figure 1: stem maps ----------------------------------------------------
trap("figure 1", {
  png(file.path(outdir, "2026-09-17_pef-stem-maps_FINAL.png"),
      width = 2400, height = 1200, res = 300)
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 0.8), bty = "l",
            cex.axis = 0.8, cex.lab = 0.9)
  for (p in unique(ci$plot_id)) {
    d <- ci[plot_id == p]
    plot(d$x, d$y, cex = 0.35 + d$dbh / 15, pch = 21, lwd = 0.6,
         bg = ifelse(d$sp_type == "SW", "#2b6a4d", "#c08a2e"),
         col = "grey25", asp = 1, xlim = c(-17, 17), ylim = c(-17, 17),
         xlab = "Easting (m)", ylab = "Northing (m)", main = p)
    symbols(0, 0, circles = 16.05, inches = FALSE, add = TRUE, fg = "grey55")
  }
  par(op); dev.off()
})

# --- Figure 2: index correlation --------------------------------------------
trap("figure 2", {
  keep <- intersect(c("bal", "ccfl", "local_bal", "hegyi", "martin_ek",
                      "spurr", "apa", "rapa", "csax_rel", "osv_rel"), names(ci))
  m <- stats::cor(as.matrix(ci[, ..keep]), use = "pairwise.complete.obs")
  png(file.path(outdir, "2026-09-17_pef-index-correlation_FINAL.png"),
      width = 1800, height = 1700, res = 300)
  op <- par(mar = c(6.5, 6.5, 1.5, 1.5))
  image(seq_len(ncol(m)), seq_len(ncol(m)), t(m[nrow(m):1, ]), axes = FALSE,
        xlab = "", ylab = "", zlim = c(-1, 1),
        col = grDevices::hcl.colors(41, "Blue-Red 3", rev = TRUE))
  axis(1, seq_len(ncol(m)), colnames(m), las = 2, cex.axis = 0.75)
  axis(2, seq_len(ncol(m)), rev(rownames(m)), las = 2, cex.axis = 0.75)
  for (i in seq_len(ncol(m))) for (j in seq_len(ncol(m))) {
    text(i, ncol(m) - j + 1, sprintf("%.2f", m[j, i]), cex = 0.5,
         col = ifelse(abs(m[j, i]) > 0.6, "white", "grey15"))
  }
  box(); par(op); dev.off()
})

# --- Figure 3: crown profile families ---------------------------------------
trap("figure 3", {
  h <- 20; hcb <- 10; lcw <- 4
  zz <- seq(hcb, h, length.out = 400)
  fams <- c("dual_exponent", "variable_exponent", "geometric")
  cols <- c("#2b6a4d", "#c08a2e", "#4a4a4a")
  png(file.path(outdir, "2026-09-17_crown-profile-families_FINAL.png"),
      width = 1600, height = 1600, res = 300)
  op <- par(mar = c(4.2, 4.2, 1.5, 1.5), bty = "l")
  plot(NA, xlim = c(-2.4, 2.4), ylim = c(hcb - 0.4, h + 0.4),
       xlab = "Crown radius (m)", ylab = "Height above ground (m)")
  for (k in seq_along(fams)) {
    pr <- crown_profile(fams[k])
    r <- crown_radius_at(pr, zz, lcw, h, hcb, widest = 0.8, up_exp = 3,
                         lo_exp = 3, shape = "g", cr = 0.5, hw = 0)
    lines(r, zz, col = cols[k], lwd = 2)
    lines(-r, zz, col = cols[k], lwd = 2)
  }
  abline(h = h - 0.8 * (h - hcb), lty = 3, col = "grey60")
  legend("topright", fams, col = cols, lwd = 2, bty = "n", cex = 0.8)
  par(op); dev.off()
})

gc()
cat("\nDone. Outputs in ", normalizePath(outdir), "\n", sep = "")
if (file.exists(elog)) cat("Errors were logged to ", elog, "\n", sep = "")
