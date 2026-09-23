#' Penobscot Experimental Forest stem-mapped example plots
#'
#' Two fully stem-mapped permanent plots from the Penobscot Experimental
#' Forest, Bradley, Maine, in the Acadian Forest Ecoregion, used throughout
#' the examples and vignette. Each plot is 0.0809 ha (16.05 m radius) and every
#' tree carries a distance and bearing from plot centre, diameter, total
#' height, height to crown base, and crown radii measured in the four cardinal
#' directions.
#'
#' The species mix is the usual Acadian one: red spruce (*Picea rubens*
#' Sarg.), balsam fir (*Abies balsamea* (L.) Mill.), eastern hemlock (*Tsuga
#' canadensis* (L.) Carriere), eastern white pine (*Pinus strobus* L.), white
#' spruce (*P. glauca* (Moench) Voss), red maple (*Acer rubrum* L.) and paper
#' birch (*Betula papyrifera* Marshall).
#'
#' @format A `data.frame` with 74 rows and 16 columns:
#' \describe{
#'   \item{plot_id}{Plot and year identifier.}
#'   \item{tree_id}{Tree number, unique within a plot.}
#'   \item{mu, plot, year}{Management unit, plot number, measurement year.}
#'   \item{species}{Two letter Acadian variant species code.}
#'   \item{distance}{Horizontal distance from plot centre, m.}
#'   \item{azimuth}{Bearing from plot centre, degrees.}
#'   \item{expf}{Expansion factor, trees ha-1.}
#'   \item{dbh}{Diameter at breast height, cm.}
#'   \item{height}{Total height, m.}
#'   \item{hcb}{Height to crown base, m.}
#'   \item{cr_north, cr_east, cr_south, cr_west}{Crown radius in each cardinal
#'     direction, m.}
#' }
#' @source USDA Forest Service Penobscot Experimental Forest permanent plot
#'   network, as distributed with the University of Maine SFR 575 competition
#'   index laboratory.
#' @examples
#' data(pef)
#' str(pef)
"pef"

#' Built-in species trait table
#'
#' The trait defaults shipped with the package. Use [species_traits()] rather
#' than this object directly, since that function handles user overrides and
#' softwood or hardwood fallback.
#'
#' @format A `data.table` with 43 rows. See [species_traits()] for the columns.
"species_traits_default"

#' Shipped coefficients for the variable exponent crown profile
#'
#' A named numeric vector of eight coefficients, `c0` to `c3` for the upper
#' crown and `d0` to `d3` for the lower crown. These are a least squares
#' calibration of the variable exponent form to the dual exponent softwood and
#' hardwood defaults, produced by `data-raw/variable_exponent.R`. They are not
#' an independent fit to measured crown profiles and should not be cited as
#' one; [crown_profile_info()] reports this.
#'
#' @format Named numeric vector of length 8.
"variable_exponent_defaults"

#' The Penobscot example as a ready-made stand
#'
#' Convenience wrapper that loads [pef] and passes it through [as_stand()] with
#' the correct column mapping and plot radius. Used by the examples so that
#' each one does not have to repeat the mapping.
#'
#' @param quiet Suppress the species fallback message.
#' @return A `stand`.
#' @export
#' @examples
#' pef_stand()
pef_stand <- function(quiet = TRUE) {
  pef <- NULL
  utils::data("pef", package = "forestCI", envir = environment())
  as_stand(pef, plot = "plot_id", tree = "tree_id", species = "species",
           dbh = "dbh", height = "height", hcb = "hcb", expf = "expf",
           distance = "distance", azimuth = "azimuth",
           plot_radius = 16.0509, quiet = quiet)
}

#' CONUS maximum and largest crown width coefficients
#'
#' A metric conversion of the CONUS maximum and largest crown width library,
#' one row per FIA species code, used when [species_traits()] is called with
#' `source = "conus"`. Maximum crown width is the tau 0.95 quantile of crown
#' width on diameter at crown ratio 1, fit as a Bayesian hierarchical quantile
#' regression on the power form MCW = `mcw_a1` * DBH^`mcw_a2` with species
#' nested in family and crossed with EPA Level III ecoregion, the ecoregion
#' effect marginalized out. Largest crown width is a bounded fraction of that
#' envelope, logit(LCW/MCW) = `lcw_r0` + `lcw_rcr` * CR + `lcw_rdbh` * DBH.
#'
#' The English unit library was converted analytically:
#' `mcw_a1` = 0.3048 * a1_en * 2.54^(-a2), `mcw_a2` unchanged, diameters
#' multiplied by 2.54, crown widths by 0.3048, and `lcw_rdbh` divided by 2.54
#' because it multiplies diameter inside a logit. The conversion was checked
#' against the source project's own independent metric export and agrees on
#' all 221 fitted species to a maximum relative difference of 5.0e-15 in
#' `mcw_a1`, exactly in `mcw_a2`, and 1.4e-13 cm in `dbh_max_fit`.
#'
#' Read `provenance` and `confidence` before relying on a species. Only 221 of
#' the 466 codes are fitted; the rest are carried by a genus or species group
#' donor or by the clade population equation, and `source_spcd` names the
#' species whose coefficients were actually used. See [max_crown_width()] for
#' the three step prediction recipe and for the caveats that apply to this
#' library as a whole.
#'
#' @format A `data.table` with 466 rows and 22 columns. Units are metric: DBH
#'   cm, crown width m.
#' \describe{
#'   \item{spcd}{FIA species code.}
#'   \item{common_name, genus, family, clade}{Identity. `clade` is S or H.}
#'   \item{provenance}{fitted, genus_donor, spgrp_donor or clade_population.}
#'   \item{source_spcd}{The species whose coefficients were used.}
#'   \item{confidence, source_confidence}{Five level shrinkage based class.}
#'   \item{n_obs}{Crown measured trees behind the fit.}
#'   \item{mcw_a1, mcw_a2}{Power coefficients, trait corrected level.}
#'   \item{dbh_min_fit, dbh_max_fit}{Fitted diameter domain, cm. In 17 of
#'     the 466 rows the two are equal, because the row carries a single
#'     crown measured tree; all 17 are donor carried or borrowed and none
#'     of them is an FIA species code that the forestCI trait table maps
#'     to, so no shipped prediction is held at a point by them. Every
#'     species has `dbh_min_fit` of at least 12.7 cm and the recipe
#'     applies no lower guard, so a small tree inventory extrapolates
#'     below the fitted range without warning.}
#'   \item{cw_ceiling}{Widest crown observed in the clade, m.}
#'   \item{cw_shipped_max}{What the recipe returns at and above `dbh_max_fit`,
#'     the lesser of the curve there and the ceiling.}
#'   \item{ceiling_binds_in_fit}{The curve breaches the ceiling inside its own
#'     fitted range.}
#'   \item{lcw_r0, lcw_rcr, lcw_rdbh}{Logit ratio terms. `lcw_rcr` and
#'     `lcw_rdbh` are global constants, not species specific.}
#'   \item{library_version}{The library export these coefficients came from.}
#' }
#' @source CONUS crown width library, FIA Phase 3 and Forest Health Monitoring
#'   crown width data after Bechtold (2003). No plot coordinates or per plot
#'   records are used or carried.
"conus_crown_width"
