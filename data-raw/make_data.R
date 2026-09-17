# Builds every .rda in data/. Run from the package root with
#   Rscript data-raw/make_data.R
# Requires only base R plus data.table.
library(data.table)

# 1. Penobscot example plots -------------------------------------------
pef <- read.csv("inst/extdata/pef.csv", stringsAsFactors = FALSE)
pef$plot_id <- as.character(pef$plot_id)
pef$tree_id <- as.character(pef$tree_id)
save(pef, file = "data/pef.rda", compress = "xz", version = 2)

# 2. Species trait defaults --------------------------------------------
species_traits_default <- as.data.table(
  read.csv("inst/extdata/species_traits.csv", stringsAsFactors = FALSE,
           na.strings = c("NA", ""))
)
setkeyv(species_traits_default, "code")
save(species_traits_default, file = "data/species_traits_default.rda",
     compress = "xz", version = 2)

# 3. Variable exponent calibration -------------------------------------
# Fits r(z) = (1 - z)^(1/k), k = a0 + a1*z + a2*CR + a3*HW, to the
# dual exponent softwood and hardwood defaults over a grid of crown ratios.
# This is a calibration to the package's own default geometry, recorded as
# such in crown_profile_info(); it is not a fit to measured crown profiles.
dual_rel <- function(z, exp_, shape) {
  xx <- if (shape == "p") 1 else exp_
  yy <- if (shape == "g") 1 else exp_
  pmax(1 - z^xx, 0)^(1 / yy)
}

fit_section <- function(section) {
  grid <- expand.grid(z = seq(0.01, 0.99, by = 0.01),
                      cr = seq(0.2, 0.9, by = 0.1),
                      hw = c(0, 1))
  # softwood default: up_exp 3, lo_exp 3, shape "g"
  # hardwood default: up_exp 3, lo_exp 4, shape "e"
  grid$target <- mapply(function(z, hw) {
    if (hw == 0) dual_rel(z, 3, "g")
    else dual_rel(z, if (section == "upper") 3 else 4, "e")
  }, grid$z, grid$hw)
  obj <- function(p) {
    k <- pmax(p[1] + p[2] * grid$z + p[3] * grid$cr + p[4] * grid$hw, 0.05)
    pred <- pmax(1 - grid$z, 0)^(1 / k)
    sum((pred - grid$target)^2)
  }
  opt <- optim(c(2, 0, 0, 0), obj, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-12))
  stopifnot(opt$convergence == 0)
  list(par = opt$par, rmse = sqrt(opt$value / nrow(grid)))
}

up <- fit_section("upper")
lo <- fit_section("lower")
variable_exponent_defaults <- c(
  c0 = up$par[1], c1 = up$par[2], c2 = up$par[3], c3 = up$par[4],
  d0 = lo$par[1], d1 = lo$par[2], d2 = lo$par[3], d3 = lo$par[4]
)
cat(sprintf("variable exponent calibration RMSE: upper %.5f, lower %.5f\n",
            up$rmse, lo$rmse))
save(variable_exponent_defaults, file = "data/variable_exponent_defaults.rda",
     compress = "xz", version = 2)
