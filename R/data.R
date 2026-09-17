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
