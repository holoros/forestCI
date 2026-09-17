#' Maximum and largest crown width
#'
#' `max_crown_width()` returns maximum crown width, the crown width an
#' open-grown tree of that diameter would carry, and is the quantity crown
#' competition factor is built from. `largest_crown_width()` returns the
#' widest crown width of the tree as it actually stands, which is the radius
#' the crown profile is scaled to.
#'
#' Both use the power forms of the Acadian variant,
#' MCW = a1 * DBH^a2 and LCW = MCW / (b1 * DBH^b2), with MCW and LCW in m and
#' DBH in cm. These are Eq. 1 and Eq. 2 of Russell and Weiskittel (2011).
#' Coefficients come from the species trait table, so a species with no local
#' coefficients falls back to the softwood or hardwood default rather than
#' failing.
#'
#' Eq. 3 of the same paper adds crown ratio to the denominator,
#' LCW = MCW / (c1 * DBH^c2 * CR^c3), and fits better for 13 of the 15 species.
#' It is not shipped here. The Acadian implementation that this package ports
#' and the published coefficient table as machine-extracted disagree on the
#' sign of the diameter exponent for several species, and on the value for
#' yellow birch, and that discrepancy is unresolved against the rendered paper.
#' Supplying the published coefficients through the `extra` argument of
#' [species_traits()] is the intended route once a reader has checked them.
#'
#' `largest_crown_width()` caps its result at maximum crown width. With the
#' shipped coefficients the ratio form stays below one over any realistic
#' diameter, so the cap binds only below about 1 cm DBH; it is a guard against
#' a user supplied coefficient set, not a correction to these.
#'
#' @param dbh Numeric vector of diameter at breast height, cm.
#' @param species Character vector of species codes, recycled against `dbh`.
#' @param traits Trait table, normally from [species_traits()].
#' @param mcw Maximum crown width, m. Supplied to `largest_crown_width()` to
#'   avoid recomputing it.
#' @return Numeric vector of crown widths, m.
#' @references
#' Russell, M.B., Weiskittel, A.R. (2011) Maximum and largest crown width
#' equations for 15 tree species in Maine. Northern Journal of Applied
#' Forestry 28: 84-91.
#' @export
#' @examples
#' max_crown_width(c(10, 25, 40), c("RS", "RS", "SM"))
max_crown_width <- function(dbh, species, traits = species_traits()) {
  p <- lookup_traits(species, traits, c("mcw_a1", "mcw_a2"))
  out <- p$mcw_a1 * dbh^p$mcw_a2
  out[is.na(dbh) | dbh <= 0] <- NA_real_
  out
}

#' @rdname max_crown_width
#' @export
largest_crown_width <- function(dbh, species, traits = species_traits(),
                                mcw = NULL) {
  if (is.null(mcw)) mcw <- max_crown_width(dbh, species, traits)
  p <- lookup_traits(species, traits, c("lcw_b1", "lcw_b2"))
  denom <- p$lcw_b1 * dbh^p$lcw_b2
  out <- mcw / denom
  # LCW is bounded above by MCW by definition; the ratio form can exceed it
  # for small diameters in some species, so cap rather than return a crown
  # wider than the open-grown maximum.
  pmin(out, mcw, na.rm = FALSE)
}

#' Maximum crown area and crown competition factor contribution
#'
#' Maximum crown area is the ground area a single open-grown tree of that
#' diameter would occupy, expressed as a percentage of one hectare and scaled
#' by the tree expansion factor. Summed over a plot it gives crown competition
#' factor.
#'
#' @inheritParams max_crown_width
#' @param expf Expansion factor, trees ha-1 represented by each record.
#' @return Numeric vector, percent of a hectare.
#' @export
max_crown_area <- function(dbh, species, expf, traits = species_traits()) {
  mcw <- max_crown_width(dbh, species, traits)
  100 * ((pi * (mcw / 2)^2) / 10000) * expf
}

# Vectorised trait lookup with softwood or hardwood fallback
lookup_traits <- function(species, traits, cols) {
  traits <- data.table::as.data.table(traits)
  idx <- match(as.character(species), traits$code)
  defaults <- type_defaults(traits)
  hw <- which(defaults$sp_type == "HW")
  out <- lapply(cols, function(cl) {
    v <- traits[[cl]][idx]
    v[is.na(idx)] <- defaults[[cl]][hw]
    v
  })
  stats::setNames(out, cols)
}
