#' Distance-independent competition indices
#'
#' Computes the non-spatial family: basal area, basal area of larger trees,
#' stand density index, relative density, crown competition factor and crown
#' competition factor of larger trees. All are per plot, and all require only
#' species, diameter and an expansion factor, so they are available for any
#' stand whether or not stem coordinates were measured.
#'
#' Softwood and hardwood partitions are returned for the one-sided indices
#' (`bal_sw`, `bal_hw`, `ccfl_sw`, `ccfl_hw`), because the Acadian diameter and
#' height increment equations take them separately and because the partition is
#' the cheapest available proxy for the shade tolerance of the competitors.
#'
#' Stand density index is summed over trees, SDI = sum of (DBH_i / 25.4)^1.605
#' times the expansion factor, rather than computed from quadratic mean
#' diameter, since the summation form is unbiased under irregular diameter
#' distributions. The quadratic mean diameter form is returned alongside it as
#' `sdi_reineke` for comparison.
#'
#' @param stand A `stand` from [as_stand()].
#' @param reineke_slope Slope of the self-thinning line. Default 1.605.
#' @param sdi_max Maximum stand density index used for relative density. Either
#'   a single number, or `"acadian"` to use the mixed species prediction below,
#'   or `"woodall"` for SDImax = 2098.6 * SG + 2057.3 where SG is the basal area
#'   weighted specific gravity, or `NULL` to skip relative density.
#' @param stand_age Stand age in years, used only by the `"acadian"` maximum
#'   SDI prediction. Default 50.
#'
#' @section Provenance of the maximum SDI predictions:
#' The `"acadian"` option reproduces, coefficient for coefficient, the mixed
#' species maximum SDI equation as implemented in Weiskittel's Acadian analysis
#' scripts. It predicts maximum SDI from the hardwood proportion of basal area,
#' basal area weighted specific gravity, the diameter range, the species count
#' and stand age. The `"woodall"` option is the specific gravity relation of
#' Woodall et al. (2005). Neither is refitted here; both are ported.
#'
#' @return A `data.table` of tree level indices with columns `plot_id`,
#'   `tree_id`, `bal`, `bal_sw`, `bal_hw`, `ccfl`, `ccfl_sw`, `ccfl_hw`,
#'   `sdi_tree`, and the plot level `ba`, `ccf`, `tph`, `qmd`, `sdi`,
#'   `sdi_reineke`, `sdi_max` and `rd` repeated on every tree of the plot.
#' @references
#' Reineke, L.H. (1933) Perfecting a stand-density index for even-aged forests.
#' Journal of Agricultural Research 46: 627-638.
#'
#' Woodall, C.W., Miles, P.D., Vissage, J.S. (2005) Determining maximum stand
#' density index in mixed species stands for strategic-scale stocking
#' assessments. Forest Ecology and Management 216: 367-377.
#' @export
#' @examples
#' data(pef)
#' st <- pef_stand()
#' head(ci_distance_independent(st))
ci_distance_independent <- function(stand, reineke_slope = 1.605,
                                    sdi_max = "acadian", stand_age = 50) {
  stopifnot(inherits(stand, "stand"))
  tr <- data.table::copy(stand$trees)
  tr[, "is_sw" := as.integer(tr$sp_type == "SW")]
  tr[, "ba_sw" := tr$ba_ha * tr$is_sw]
  tr[, "mca_sw" := tr$mca * tr$is_sw]
  tr[, "sdi_tree" := (tr$dbh / 25.4)^reineke_slope * tr$expf]

  # One-sided sums: everything strictly larger than the subject. Ties in DBH
  # are excluded from each other's "larger" set, which keeps the index
  # symmetric under tied diameters rather than depending on row order.
  one_sided <- function(size, value) {
    o <- order(-size)
    s <- size[o]; v <- value[o]
    cs <- cumsum(v) - v                     # sum of strictly earlier rows
    # pull tied groups back to the sum above the tie
    grp <- match(s, unique(s))
    first_of_grp <- !duplicated(grp)
    base <- cs[first_of_grp][grp]
    out <- numeric(length(size))
    out[o] <- base
    out
  }

  tr[, "bal"     := one_sided(dbh, ba_ha),  by = "plot_id"]
  tr[, "bal_sw"  := one_sided(dbh, ba_sw),  by = "plot_id"]
  tr[, "bal_hw"  := bal - bal_sw]
  tr[, "ccfl"    := one_sided(dbh, mca),    by = "plot_id"]
  tr[, "ccfl_sw" := one_sided(dbh, mca_sw), by = "plot_id"]
  tr[, "ccfl_hw" := ccfl - ccfl_sw]

  ps <- tr[, list(
    ba_plot  = sum(ba_ha, na.rm = TRUE),
    ccf      = sum(mca, na.rm = TRUE),
    tph      = sum(expf, na.rm = TRUE),
    sdi      = sum(sdi_tree, na.rm = TRUE),
    qmd      = qmd_from(dbh, expf),
    sg_ba    = stats::weighted.mean(sg, ba_ha, na.rm = TRUE),
    p_hw_ba  = sum(ba_ha * (sp_type == "HW"), na.rm = TRUE) /
               sum(ba_ha, na.rm = TRUE),
    dbh_min  = min(dbh, na.rm = TRUE),
    dbh_max  = max(dbh, na.rm = TRUE),
    n_spp    = length(unique(species))
  ), by = "plot_id"]
  ps[, "sdi_reineke" := (qmd / 25.4)^reineke_slope * tph]

  if (!is.null(sdi_max)) {
    ps[, "sdi_max" := sdi_max_predict(sdi_max, ps, stand_age)]
    ps[, "rd" := sdi / sdi_max]
  }

  keep_plot <- intersect(c("plot_id", "ba_plot", "ccf", "tph", "sdi",
                           "sdi_reineke", "qmd", "sdi_max", "rd"), names(ps))
  out <- merge(tr[, c("plot_id", "tree_id", "bal", "bal_sw", "bal_hw",
                      "ccfl", "ccfl_sw", "ccfl_hw", "sdi_tree"), with = FALSE],
               ps[, keep_plot, with = FALSE], by = "plot_id", all.x = TRUE)
  data.table::setnames(out, "ba_plot", "ba", skip_absent = TRUE)
  out[]
}

sdi_max_predict <- function(spec, ps, stand_age) {
  if (is.numeric(spec)) return(rep_len(spec, nrow(ps)))
  spec <- match.arg(as.character(spec), c("acadian", "woodall"))
  if (spec == "woodall") {
    return(2098.6 * ps$sg_ba + 2057.3)
  }
  # Acadian mixed species maximum SDI, ported verbatim from the Acadian
  # analysis script. Coefficients are not refitted here.
  483.2448 -
    1.4563 * ps$p_hw_ba -
    212.705 * log(ps$sg_ba) +
    45.351 * sqrt(pmax(ps$dbh_max - ps$dbh_min, 0)) +
    14.811 * ps$n_spp -
    0.0848 * stand_age +
    0.0001 * stand_age^2 +
    331.3714 * (1 / 12)
}
