#' Crown profile families
#'
#' A crown profile describes crown radius as a function of depth below the
#' crown apex. Every crown based competition index in the package
#' (crown surface area, exposed crown surface area, open sky view, crown
#' volume, crown projection area) is computed by integrating or sampling a
#' crown profile, so swapping the profile swaps the geometry everywhere at
#' once. Three families ship with the package and a fourth slot takes an
#' arbitrary user function.
#'
#' @section Families:
#' \describe{
#'   \item{`"dual_exponent"`}{The default. Radius at depth d below the apex is
#'     r(d) = R * (1 - (d/L)^x)^(1/y), applied separately to the section above
#'     the widest point and the section below it, with L the length of the
#'     section concerned and R the largest crown width divided by two. The
#'     `shape` trait overrides one exponent: `"p"` (paraboloid) forces x = 1,
#'     `"g"` (gothic dome) forces y = 1, `"e"` (ellipsoid) leaves both free.
#'     This is the form of Waskiewicz (2011) and Kuehne et al. (2019, Eq. 3).
#'     It needs four species constants and is the most flexible of the three.}
#'   \item{`"variable_exponent"`}{A single continuous exponent that varies with
#'     relative depth rather than two fixed exponents,
#'     r(z) = R * (1 - z)^(c0 + c1 * z + c2 * CR + c3 * HW), with z relative
#'     depth in the section, CR crown ratio and HW a hardwood indicator. It
#'     needs no species table at all, only crown ratio and a softwood or
#'     hardwood flag, so it is the family to use where species specific crown
#'     constants do not exist. The shipped coefficients are a least squares
#'     calibration of this form to the `"dual_exponent"` softwood and hardwood
#'     defaults, not an independent field fit, and are reported as such by
#'     [crown_profile_info()].}
#'   \item{`"geometric"`}{Classical solids of revolution: cone, paraboloid,
#'     ellipsoid, neiloid. One constant, no species table, closed form surface
#'     area and volume. Use it as a transparent baseline or where crown data
#'     are too thin to justify anything else.}
#'   \item{`"custom"`}{Any function `f(z, ...)` returning relative radius
#'     (0 to 1) at relative depth `z` (0 at apex, 1 at the section end).}
#' }
#'
#' @param family One of `"dual_exponent"`, `"variable_exponent"`,
#'   `"geometric"`, `"custom"`.
#' @param solid For `"geometric"`, one of `"cone"`, `"paraboloid"`,
#'   `"ellipsoid"`, `"neiloid"`.
#' @param coef Named numeric vector of coefficients for
#'   `"variable_exponent"`: `c0`, `c1`, `c2`, `c3` for the upper crown and
#'   `d0`, `d1`, `d2`, `d3` for the lower crown. Defaults to the shipped
#'   calibration.
#' @param fun For `"custom"`, a vectorised function of relative depth `z`
#'   returning relative radius.
#' @param widest_default Position of the widest point as a proportion of crown
#'   length from the apex, used when the trait table supplies none. Given as a
#'   length two vector, softwood then hardwood.
#' @param ... Additional arguments stored on the object and passed to `fun`.
#'
#' @return An object of class `crown_profile`.
#' @references
#' Kuehne, C., Weiskittel, A.R., Waskiewicz, J. (2019) Comparing performance of
#' contrasting distance-independent and distance-dependent competition metrics
#' in predicting individual tree diameter increment and survival within
#' structurally-heterogeneous, mixed-species forests of Northeastern United
#' States. Forest Ecology and Management 433: 205-216.
#' @export
#' @examples
#' p <- crown_profile("dual_exponent")
#' crown_dimensions(p, lcw = 4, height = 20, hcb = 10, species = "RS")
crown_profile <- function(family = c("dual_exponent", "variable_exponent",
                                     "geometric", "custom"),
                          solid = c("paraboloid", "cone", "ellipsoid", "neiloid"),
                          coef = NULL, fun = NULL,
                          widest_default = c(SW = 0.80, HW = 0.60), ...) {
  family <- match.arg(family)
  solid <- match.arg(solid)
  if (family == "custom" && !is.function(fun)) {
    stop("family = \"custom\" requires `fun` to be a function of relative depth.",
         call. = FALSE)
  }
  if (family == "variable_exponent") {
    coef <- utils::modifyList(as.list(variable_exponent_coef()), as.list(coef %||% list()))
    coef <- unlist(coef)
  }
  structure(
    list(family = family, solid = solid, coef = coef, fun = fun,
         widest_default = widest_default, args = list(...)),
    class = "crown_profile"
  )
}

#' @export
print.crown_profile <- function(x, ...) {
  cat("<crown_profile>\n")
  cat("  family:", x$family, "\n")
  if (x$family == "geometric") cat("  solid: ", x$solid, "\n")
  if (x$family == "variable_exponent") {
    cat("  coef:  ", paste(sprintf("%s=%.4f", names(x$coef), x$coef),
                           collapse = ", "), "\n")
  }
  cat("  widest point default (SW, HW):",
      paste(x$widest_default, collapse = ", "), "\n")
  invisible(x)
}

#' Provenance of a crown profile
#'
#' Reports where each constant of a crown profile came from, so that a figure
#' or table built on it can state its own basis. Constants calibrated inside
#' this package are labelled as such and are not presented as field fits.
#'
#' @param profile A `crown_profile`.
#' @return A `data.frame` with columns `constant`, `value` and `source`.
#' @export
crown_profile_info <- function(profile) {
  stopifnot(inherits(profile, "crown_profile"))
  switch(profile$family,
    dual_exponent = data.frame(
      constant = c("widest", "up_exp", "lo_exp", "shape"),
      value = NA_character_,
      source = "species trait table; Waskiewicz (2011) parameterisation",
      stringsAsFactors = FALSE),
    variable_exponent = data.frame(
      constant = names(profile$coef),
      value = as.character(round(profile$coef, 6)),
      source = paste("least squares calibration of this form to the",
                     "dual_exponent softwood and hardwood defaults,",
                     "computed in data-raw/variable_exponent.R;",
                     "not an independent field fit"),
      stringsAsFactors = FALSE),
    geometric = data.frame(
      constant = "solid", value = profile$solid,
      source = "assumed solid of revolution, no fitted constants",
      stringsAsFactors = FALSE),
    custom = data.frame(
      constant = "fun", value = "user supplied",
      source = "user supplied", stringsAsFactors = FALSE)
  )
}

# Relative radius at relative position z within one crown section, where
# z = 0 is the widest point of the crown and z = 1 is the tip of that section
# (the apex for the upper section, the crown base for the lower section).
# Returns values in [0, 1], decreasing in z; multiply by the largest crown
# radius to get metres.
rel_radius <- function(profile, z, section = c("upper", "lower"),
                       up_exp = 3, lo_exp = 3, shape = "e",
                       cr = 0.5, hw = 0) {
  section <- match.arg(section)
  z <- pmin(pmax(z, 0), 1)
  switch(profile$family,
    dual_exponent = {
      e <- if (section == "upper") up_exp else lo_exp
      # shape overrides one exponent, following Waskiewicz (2011)
      xx <- ifelse(shape == "p", 1, e)
      yy <- ifelse(shape == "g", 1, e)
      pmax(1 - z^xx, 0)^(1 / yy)
    },
    variable_exponent = {
      k <- if (section == "upper") {
        profile$coef[["c0"]] + profile$coef[["c1"]] * z +
          profile$coef[["c2"]] * cr + profile$coef[["c3"]] * hw
      } else {
        profile$coef[["d0"]] + profile$coef[["d1"]] * z +
          profile$coef[["d2"]] * cr + profile$coef[["d3"]] * hw
      }
      k <- pmax(k, 0.05)
      pmax(1 - z, 0)^(1 / k)
    },
    geometric = {
      p <- switch(profile$solid,
                  cone = 1, paraboloid = 2, ellipsoid = 2, neiloid = 3)
      if (profile$solid == "ellipsoid") {
        sqrt(pmax(1 - z^2, 0))
      } else {
        pmax(1 - z, 0)^(1 / p)
      }
    },
    custom = do.call(profile$fun, c(list(z), profile$args))
  )
}

# First derivative of relative radius with respect to z, by central difference.
# Analytic where cheap, numeric otherwise; the numeric branch is stable because
# rel_radius is smooth on the open interval and the endpoints are handled by
# the integration limits.
rel_radius_deriv <- function(profile, z, ...) {
  h <- 1e-5
  zl <- pmax(z - h, 0)
  zu <- pmin(z + h, 1)
  (rel_radius(profile, zu, ...) - rel_radius(profile, zl, ...)) / (zu - zl)
}

#' Crown radius at a given height
#'
#' @param profile A `crown_profile`.
#' @param height_above_ground Numeric vector of heights, m, at which to
#'   evaluate crown radius.
#' @param lcw Largest crown width, m.
#' @param height Total tree height, m.
#' @param hcb Height to crown base, m.
#' @param widest Position of the widest point as a proportion of crown length
#'   from the apex.
#' @param up_exp,lo_exp,shape Crown shape constants, used by the
#'   `"dual_exponent"` family.
#' @param cr Crown ratio, used by the `"variable_exponent"` family.
#' @param hw Hardwood indicator, 0 or 1, used by the `"variable_exponent"`
#'   family.
#' @return Numeric vector of crown radii, m. Zero outside the crown.
#' @export
crown_radius_at <- function(profile, height_above_ground, lcw, height, hcb,
                            widest = 0.7, up_exp = 3, lo_exp = 3, shape = "e",
                            cr = NULL, hw = 0) {
  cl <- height - hcb
  if (any(cl <= 0, na.rm = TRUE)) {
    stop("crown length (height - hcb) must be positive for every tree.",
         call. = FALSE)
  }
  cr <- cr %||% (cl / height)
  rmax <- lcw / 2
  depth <- height - height_above_ground          # depth below apex, m
  ucl <- widest * cl                             # upper section length
  lcl <- cl - ucl                                # lower section length
  out <- numeric(length(depth))
  up <- depth >= 0 & depth <= ucl
  lo <- depth > ucl & depth <= cl
  # relative position within a section runs 0 at the widest point to 1 at the
  # tip of that section, matching rel_radius(), which decreases from 1 to 0.
  if (any(up)) {
    z <- ifelse(ucl > 0, (ucl - depth[up]) / ucl, 0)
    out[up] <- rmax * rel_radius(profile, z, "upper", up_exp, lo_exp, shape,
                                 cr, hw)
  }
  if (any(lo)) {
    z <- ifelse(lcl > 0, (depth[lo] - ucl) / lcl, 1)
    out[lo] <- rmax * rel_radius(profile, z, "lower", up_exp, lo_exp, shape,
                                 cr, hw)
  }
  out
}

#' Crown surface area, crown volume and crown projection area
#'
#' Integrates a crown profile to give the lateral surface area of the crown,
#' the crown volume, and the vertically projected crown area. The crown is
#' treated as a solid of revolution about the stem axis, made of an upper
#' section from the apex to the widest point and a lower section from the
#' widest point to the crown base.
#'
#' Crown projection area is pi times the square of the largest crown radius by
#' construction, since the widest point is the widest point.
#'
#' @inheritParams crown_radius_at
#' @param species Optional species codes. When supplied, `widest`, `up_exp`,
#'   `lo_exp` and `shape` are taken from `traits` and any values passed
#'   directly are ignored.
#' @param traits Trait table, normally from [species_traits()].
#' @param sp_type Optional "SW" or "HW", used to pick the softwood or hardwood
#'   default for a species code that is not in the trait table. With none the
#'   hardwood default is used, which is the historical behaviour.
#' @param n_sub Number of subdivisions used by `method = "simpson"`.
#' @param method Integration rule. `"adaptive"` calls [stats::integrate()] once
#'   per tree and per crown section, which handles the vertical tangent that the
#'   profile develops at the tip of the crown and is accurate to the default
#'   tolerance of that routine. `"simpson"` uses a fixed composite Simpson rule
#'   over `n_sub` subdivisions, which is several times faster and vectorised
#'   over trees, but converges slowly near that tangent and can be a few percent
#'   low on the most sharply pointed crowns. `"adaptive"` reproduces the
#'   integration used by the original Acadian competition index scripts.
#' @return A `data.frame` with columns `csa` (crown surface area, m2),
#'   `csa_upper`, `csa_lower`, `cv` (crown volume, m3) and `cpa` (crown
#'   projection area, m2).
#' @export
#' @examples
#' crown_dimensions(crown_profile(), lcw = 4, height = 20, hcb = 10,
#'                  species = "RS")
crown_dimensions <- function(profile, lcw, height, hcb,
                             species = NULL, traits = species_traits(),
                             widest = 0.7, up_exp = 3, lo_exp = 3, shape = "e",
                             cr = NULL, hw = 0, sp_type = NULL, n_sub = 200L,
                             method = c("adaptive", "simpson")) {
  method <- match.arg(method)
  stopifnot(inherits(profile, "crown_profile"))
  n <- max(length(lcw), length(height), length(hcb))
  lcw <- rep_len(lcw, n); height <- rep_len(height, n); hcb <- rep_len(hcb, n)

  if (!is.null(species)) {
    species <- rep_len(as.character(species), n)
    p <- lookup_traits(species, traits,
                       c("widest", "up_exp", "lo_exp", "shape", "sp_type"),
                       sp_type = sp_type)
    widest <- as.numeric(p$widest)
    up_exp <- as.numeric(p$up_exp)
    lo_exp <- as.numeric(p$lo_exp)
    shape  <- as.character(p$shape)
    hw     <- as.numeric(p$sp_type == "HW")
  } else {
    widest <- rep_len(widest, n); up_exp <- rep_len(up_exp, n)
    lo_exp <- rep_len(lo_exp, n); shape <- rep_len(shape, n)
    hw <- rep_len(hw, n)
  }
  cl <- height - hcb
  cr <- rep_len(cr %||% (cl / height), n)
  rmax <- lcw / 2
  ucl <- widest * cl
  lcl <- cl - ucl

  if (method == "adaptive") {
    sec_adaptive <- function(section, seclen, i) {
      if (!is.finite(seclen[i]) || seclen[i] <= 0) return(c(0, 0))
      rfun <- function(z) {
        rmax[i] * rel_radius(profile, z, section, up_exp[i], lo_exp[i],
                             shape[i], cr[i], hw[i])
      }
      dfun <- function(z) {
        rmax[i] * rel_radius_deriv(profile, z, section = section,
                                   up_exp = up_exp[i], lo_exp = lo_exp[i],
                                   shape = shape[i], cr = cr[i], hw = hw[i])
      }
      # integrate over length along the section axis, l = z * seclen
      area_f <- function(l) {
        z <- l / seclen[i]
        rfun(z) * sqrt(1 + (dfun(z) / seclen[i])^2)
      }
      vol_f <- function(l) rfun(l / seclen[i])^2
      a <- tryCatch(stats::integrate(area_f, 0, seclen[i],
                                     subdivisions = 1000L,
                                     rel.tol = 1e-8, stop.on.error = FALSE)$value,
                    error = function(e) NA_real_)
      v <- tryCatch(stats::integrate(vol_f, 0, seclen[i],
                                     subdivisions = 1000L,
                                     rel.tol = 1e-8, stop.on.error = FALSE)$value,
                    error = function(e) NA_real_)
      c(2 * pi * a, pi * v)
    }
    up_a <- up_v <- lo_a <- lo_v <- numeric(n)
    for (i in seq_len(n)) {
      u <- sec_adaptive("upper", ucl, i); l <- sec_adaptive("lower", lcl, i)
      up_a[i] <- u[1]; up_v[i] <- u[2]; lo_a[i] <- l[1]; lo_v[i] <- l[2]
    }
    return(data.frame(
      csa       = up_a + lo_a,
      csa_upper = up_a,
      csa_lower = lo_a,
      cv        = up_v + lo_v,
      cpa       = pi * rmax^2
    ))
  }

  # composite Simpson on z in (0, 1), vectorised over trees
  m <- if (n_sub %% 2L == 0L) n_sub else n_sub + 1L
  z <- seq(0, 1, length.out = m + 1L)
  w <- c(1, rep(c(4, 2), length.out = m - 1L), 1) / 3 * (1 / m)

  sec_dims <- function(section, seclen, exps_up, exps_lo) {
    # returns list(area, vol) per tree for one section
    area <- numeric(n); vol <- numeric(n)
    for (k in seq_along(z)) {
      zz <- rep_len(z[k], n)
      rr <- rel_radius(profile, zz, section, exps_up, exps_lo, shape, cr, hw)
      dr <- rel_radius_deriv(profile, zz, section = section,
                             up_exp = exps_up, lo_exp = exps_lo,
                             shape = shape, cr = cr, hw = hw)
      # r(t) = rmax * rr, dr/dl = rmax * dr / seclen
      slope <- ifelse(seclen > 0, rmax * dr / seclen, 0)
      area <- area + w[k] * (rmax * rr) * sqrt(1 + slope^2) * seclen
      vol  <- vol  + w[k] * (rmax * rr)^2 * seclen
    }
    list(area = 2 * pi * area, vol = pi * vol)
  }

  up <- sec_dims("upper", ucl, up_exp, lo_exp)
  lo <- sec_dims("lower", lcl, up_exp, lo_exp)

  data.frame(
    csa       = up$area + lo$area,
    csa_upper = up$area,
    csa_lower = lo$area,
    cv        = up$vol + lo$vol,
    cpa       = pi * rmax^2
  )
}

# Shipped calibration for the variable_exponent family. Values are produced by
# data-raw/variable_exponent.R, which fits this form to the dual_exponent
# softwood and hardwood defaults. They are a calibration, not a field fit.
variable_exponent_coef <- function() {
  forestCI::variable_exponent_defaults
}
