#' Maximum and largest crown width
#'
#' `max_crown_width()` returns maximum crown width, the crown width an
#' open-grown tree of that diameter would carry, and is the quantity crown
#' competition factor is built from. `largest_crown_width()` returns the
#' widest crown width of the tree as it actually stands, which is the radius
#' the crown profile is scaled to.
#'
#' Two coefficient sources are supported, selected through the trait table
#' rather than through these functions, so that [species_traits()] remains the
#' single place a coefficient set is chosen. Which one a row carries is
#' recorded in its `cw_form` column.
#'
#' @section The Acadian form (`cw_form == "acadian"`, the default):
#' The power forms of the Acadian variant, MCW = a1 * DBH^a2 and
#' LCW = MCW / (b1 * DBH^b2), with MCW and LCW in m and DBH in cm. These are
#' Eq. 1 and Eq. 2 of Russell and Weiskittel (2011), 15 species in Maine.
#' `cr` is not used. The curve has no asymptote and no domain limit is
#' applied, which is how the original Acadian implementation behaves and is
#' what the regression against it certifies.
#'
#' Eq. 3 of the same paper adds crown ratio to the denominator,
#' LCW = MCW / (c1 * DBH^c2 * CR^c3), and fits better for 13 of the 15
#' species. It is still not shipped. The Acadian implementation that this
#' package ports and the published coefficient table as machine-extracted
#' disagree on the sign of the diameter exponent for several species, and on
#' the value for yellow birch, and that discrepancy remains unresolved against
#' the rendered paper. The CONUS source below does not settle it: it is a
#' different model of a different quantity, not a reading of Table 3.
#'
#' @section The CONUS form (`cw_form == "conus"`):
#' A metric conversion of the CONUS maximum and largest crown width library,
#' fit on FIA Phase 3 and Forest Health Monitoring crown width data. MCW is the
#' tau 0.95 quantile of crown width on diameter at crown ratio 1, evaluated in
#' three steps whose order is load bearing. The curve MCW = a1 * DBH^a2 is
#' evaluated on the trait corrected level; any diameter above `mcw_dbh_max` is
#' held at `mcw_dbh_max`, so the curve is never evaluated outside the range it
#' was fitted over; and the result is capped at `mcw_ceiling`, the widest crown
#' observed in that clade, 22.86 m hardwood and 19.87 m softwood. The cap is
#' not implied by the hold: for some species the curve breaches the ceiling
#' inside its own fitted range.
#'
#' LCW is not a second allometry under this form. It is a bounded fraction of
#' the MCW envelope, LCW = MCW * ratio with
#' logit(ratio) = `lcw_r0` + `lcw_rcr` * CR + `lcw_rdbh` * DBH, so LCW lies
#' below MCW at every diameter by construction. The ratio is evaluated on the
#' held diameter, so LCW goes flat exactly where MCW does. `cr` defaults to 1,
#' the crown ratio the envelope is delivered at.
#'
#' @section Caveats on the CONUS source:
#' Three, all of which the source project records and none of which this
#' package can repair. The envelope is a 0.95 quantile of stand grown trees
#' rather than a conditional mean of open grown trees, so crown competition
#' factor built from it is not on the scale Krajicek et al. (1961) calibrated.
#' The domain hold goes flat inside a normal inventory for most Northeastern
#' species, which understates competition in larger stands. And the LCW ratio
#' carries a single global crown ratio slope and a single global diameter
#' slope, so the whole between-species difference in LCW/MCW is one intercept.
#' Below `dbh_min_fit` no guard is applied and the curve is extrapolated.
#' Use [conus_crown_width] to inspect the fitted domain, provenance and
#' confidence class of any species before relying on it.
#'
#' `largest_crown_width()` caps its result at maximum crown width under either
#' form. Under the CONUS form the bound already holds by construction, so the
#' cap is a guard against a user supplied coefficient set rather than a
#' correction.
#'
#' @param dbh Numeric vector of diameter at breast height, cm.
#' @param species Character vector of species codes, recycled against `dbh`.
#' @param traits Trait table, normally from [species_traits()].
#' @param sp_type Optional character vector, `"SW"` or `"HW"`, used to pick the
#'   softwood or hardwood default for a species code that is not in the trait
#'   table. Recycled against `species`. When it is not supplied the hardwood
#'   default is used, which is the historical behaviour.
#' @param mcw Maximum crown width, m. Supplied to `largest_crown_width()` to
#'   avoid recomputing it.
#' @param cr Crown ratio as a proportion in (0, 1], default 1. Used by the
#'   CONUS form and ignored by the Acadian form. A missing value is treated as
#'   1, the open-grown level the envelope is delivered at.
#' @return Numeric vector of crown widths, m.
#' @references
#' Russell, M.B., Weiskittel, A.R. (2011) Maximum and largest crown width
#' equations for 15 tree species in Maine. Northern Journal of Applied
#' Forestry 28: 84-91.
#'
#' Krajicek, J.E., Brinkman, K.A., Gingrich, S.F. (1961) Crown competition, a
#' measure of density. Forest Science 7: 35-42.
#'
#' Bechtold, W.A. (2003) Crown-diameter prediction models for 87 species of
#' stand-grown trees in the eastern United States. Southern Journal of Applied
#' Forestry 27: 269-278.
#' @export
#' @examples
#' max_crown_width(c(10, 25, 40), c("RS", "RS", "SM"))
#' # the same trees on the CONUS library
#' max_crown_width(c(10, 25, 40), c("RS", "RS", "SM"),
#'                 traits = species_traits(source = "conus"))
max_crown_width <- function(dbh, species, traits = species_traits(),
                            sp_type = NULL) {
  cols <- c("mcw_a1", "mcw_a2")
  extra <- c("cw_form", "mcw_dbh_max", "mcw_ceiling")
  has_conus <- all(extra %in% names(traits))
  p <- lookup_traits(species, traits, if (has_conus) c(cols, extra) else cols,
                     sp_type = sp_type)
  d <- as.numeric(dbh)
  # recycle the coefficient columns to the length of the result before any
  # logical subsetting: a single species code returns length one columns, and
  # indexing those with a longer logical silently yields NA
  n <- max(length(d), length(p[[1L]]))
  d <- rep_len(d, n)
  p <- lapply(p, rep_len, n)
  dd <- d
  if (has_conus) {
    conus <- !is.na(p$cw_form) & p$cw_form == "conus"
    hold <- conus & !is.na(p$mcw_dbh_max) & is.finite(d) & d > p$mcw_dbh_max
    dd[hold] <- p$mcw_dbh_max[hold]
  }
  out <- p$mcw_a1 * dd^p$mcw_a2
  if (has_conus) {
    cap <- conus & !is.na(p$mcw_ceiling) & is.finite(out) & out > p$mcw_ceiling
    out[cap] <- p$mcw_ceiling[cap]
  }
  out[is.na(d) | d <= 0] <- NA_real_
  out
}

#' @rdname max_crown_width
#' @export
largest_crown_width <- function(dbh, species, traits = species_traits(),
                                mcw = NULL, cr = 1, sp_type = NULL) {
  if (is.null(mcw)) mcw <- max_crown_width(dbh, species, traits, sp_type)
  d <- as.numeric(dbh)
  cols <- c("lcw_b1", "lcw_b2")
  extra <- c("cw_form", "mcw_dbh_max", "lcw_r0", "lcw_rcr", "lcw_rdbh")
  has_conus <- all(extra %in% names(traits))
  p <- lookup_traits(species, traits, if (has_conus) c(cols, extra) else cols,
                     sp_type = sp_type)
  # as in max_crown_width(), recycle before any logical subsetting
  n <- max(length(mcw), length(d), length(p[[1L]]))
  mcw <- rep_len(mcw, n)
  d <- rep_len(d, n)
  p <- lapply(p, rep_len, n)

  out <- mcw / (p$lcw_b1 * d^p$lcw_b2)

  if (has_conus) {
    conus <- !is.na(p$cw_form) & p$cw_form == "conus"
    if (any(conus)) {
      crv <- rep_len(as.numeric(cr), n)
      bad <- !is.na(crv) & (crv <= 0 | crv > 1)
      if (any(bad)) {
        stop("`cr` must lie in (0, 1]; it is a proportion, not a percentage. ",
             "Offending range ", signif(min(crv[bad]), 4), " to ",
             signif(max(crv[bad]), 4), ".", call. = FALSE)
      }
      crv[is.na(crv)] <- 1               # the level the envelope is delivered at
      dd <- d
      hold <- !is.na(p$mcw_dbh_max) & is.finite(d) & d > p$mcw_dbh_max
      dd[hold] <- p$mcw_dbh_max[hold]
      eta <- p$lcw_r0 + p$lcw_rcr * crv + p$lcw_rdbh * dd
      out[conus] <- (mcw * (1 / (1 + exp(-eta))))[conus]
    }
  }
  # LCW is bounded above by MCW by definition; the Acadian ratio form can
  # exceed it for small diameters in some species, so cap rather than return a
  # crown wider than the open-grown maximum.
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
max_crown_area <- function(dbh, species, expf, traits = species_traits(),
                           sp_type = NULL) {
  mcw <- max_crown_width(dbh, species, traits, sp_type)
  crown_area_pct(mcw, expf)
}

# Percent of a hectare occupied by an open-grown crown of width `mcw`, scaled
# by the expansion factor. Defined once so that as_stand() and max_crown_area()
# cannot drift apart.
crown_area_pct <- function(mcw, expf) {
  100 * ((pi * (mcw / 2)^2) / 10000) * expf
}

# Vectorised trait lookup with softwood or hardwood fallback. The fallback row
# is chosen per element from `sp_type` when the caller supplies one, mirroring
# attach_traits(); with no `sp_type` it is the hardwood default, which is the
# historical behaviour and what the direct-call tests assume.
lookup_traits <- function(species, traits, cols, sp_type = NULL) {
  traits <- data.table::as.data.table(traits)
  idx <- match(as.character(species), traits$code)
  defaults <- type_defaults(traits)
  ft <- if (is.null(sp_type)) {
    rep("HW", length(idx))
  } else {
    rep_len(as.character(sp_type), length(idx))
  }
  ft[is.na(ft) | !ft %in% c("SW", "HW")] <- "HW"
  fb <- match(ft, defaults$sp_type)
  miss <- is.na(idx)
  out <- lapply(cols, function(cl) {
    v <- traits[[cl]][idx]
    v[miss] <- defaults[[cl]][fb][miss]
    v
  })
  stats::setNames(out, cols)
}
