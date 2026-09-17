# Crown surface elevation of a set of trees above a set of horizontal points.
# Returns a matrix, points by trees, of the elevation (m) of the upper crown
# surface, or NA where the point lies outside that tree's crown.
crown_surface_z <- function(profile, px, py, trees, n_grid = 128L) {
  n <- nrow(trees)
  out <- matrix(NA_real_, nrow = length(px), ncol = n)
  z <- seq(0, 1, length.out = n_grid)
  for (i in seq_len(n)) {
    tt <- trees[i, ]
    cl <- tt$height - tt$hcb
    ucl <- tt$widest * cl
    if (ucl <= 0) next
    rmax <- tt$lcw / 2
    rr <- rmax * rel_radius(profile, z, "upper", tt$up_exp, tt$lo_exp,
                            tt$shape, tt$cr, as.numeric(tt$sp_type == "HW"))
    rho <- sqrt((px - tt$x)^2 + (py - tt$y)^2)
    inside <- rho <= max(rr)
    if (!any(inside)) next
    # rr decreases in z, and elevation rises from the widest point at z = 0 to
    # the apex at z = 1; invert on radius to get the surface elevation.
    zz <- stats::approx(x = rev(rr), y = rev(z), xout = rho[inside],
                        rule = 2, ties = "ordered")$y
    out[inside, i] <- (tt$height - ucl) + zz * ucl
  }
  out
}

#' Exposed crown surface area, exposed crown projection area and open sky view
#'
#' Three one-sided crown based indices, all obtained by sampling the crown
#' profile rather than by assuming a fixed solid, so they follow whatever
#' family [crown_profile()] was given.
#'
#' Crown surface area is the lateral surface of the crown as a solid of
#' revolution. Exposed crown surface area is the part of that surface not
#' overtopped by a neighbouring crown, computed by dividing the crown into
#' `n_depth` by `n_azimuth` facets, placing each facet centre in three
#' dimensions, and asking whether any other crown surface stands above it at
#' that horizontal position. Exposed crown projection area is the same test
#' applied to the vertical projection. Open sky view is stricter: a facet
#' counts only when a ray leaving it along the outward normal escapes without
#' entering any neighbouring crown, so a facet under an open gap counts and a
#' facet under a leaning neighbour does not.
#'
#' Kuehne et al. (2019) found exposed crown surface area the single most
#' useful distance-dependent metric for diameter increment, and Waskiewicz
#' (2011) found open sky view the best of everything he tested, so these are
#' the indices the crown profile machinery exists to serve.
#'
#' The relative forms `csax_rel` and `osv_rel`, each divided by crown surface
#' area, are returned alongside the absolute ones because the absolute forms
#' are strongly correlated with tree size and the relative forms are not.
#'
#' @param stand A `stand` with stem coordinates, height and height to crown
#'   base for every tree.
#' @param profile A `crown_profile`.
#' @param n_depth,n_azimuth Facet grid. Default 18 by 72, which is the 1296
#'   facet grid of Waskiewicz (2011). Reduce both for large plots.
#' @param osv Compute open sky view. It is the expensive part; set `FALSE` to
#'   skip it.
#' @param ray_steps Sample points along each open sky view ray.
#' @param ray_length Maximum ray length, m.
#' @return A `data.table` with `plot_id`, `tree_id`, `csa`, `csax`,
#'   `csax_rel`, `cpa`, `cpax`, `cpax_rel`, `cv`, and, when requested, `osv`
#'   and `osv_rel`.
#' @export
#' @examples
#' st <- pef_stand()
#' head(crown_exposure(st, n_depth = 8, n_azimuth = 24, osv = FALSE))
crown_exposure <- function(stand, profile = crown_profile(),
                           n_depth = 18L, n_azimuth = 72L, osv = TRUE,
                           ray_steps = 40L, ray_length = 30) {
  stopifnot(inherits(stand, "stand"), inherits(profile, "crown_profile"))
  if (!is_spatial(stand)) stop("this stand has no stem coordinates.", call. = FALSE)
  if (!isTRUE(attr(stand, "has_crown"))) {
    stop("crown exposure needs `height` and `hcb` for every tree.", call. = FALSE)
  }
  tr <- stand$trees

  res <- lapply(split(seq_len(nrow(tr)), tr$plot_id), function(rows) {
    p <- tr[rows, ]
    n <- nrow(p)
    dims <- crown_dimensions(profile, p$lcw, p$height, p$hcb,
                             species = p$species, traits = stand$traits)
    csa <- dims$csa; cpa <- dims$cpa; cv <- dims$cv
    csax <- numeric(n); cpax <- numeric(n); osvv <- rep(NA_real_, n)

    az <- seq(0, 2 * pi, length.out = n_azimuth + 1L)[-(n_azimuth + 1L)]
    daz <- 2 * pi / n_azimuth

    for (i in seq_len(n)) {
      tt <- p[i, ]
      cl <- tt$height - tt$hcb
      ucl <- tt$widest * cl
      lcl <- cl - ucl
      rmax <- tt$lcw / 2
      hw <- as.numeric(tt$sp_type == "HW")

      # facet centres in (section, relative depth)
      zc <- (seq_len(n_depth) - 0.5) / n_depth
      up_r <- rmax * rel_radius(profile, zc, "upper", tt$up_exp, tt$lo_exp,
                                tt$shape, tt$cr, hw)
      lo_r <- rmax * rel_radius(profile, zc, "lower", tt$up_exp, tt$lo_exp,
                                tt$shape, tt$cr, hw)
      up_dr <- rmax * rel_radius_deriv(profile, zc, section = "upper",
                                       up_exp = tt$up_exp, lo_exp = tt$lo_exp,
                                       shape = tt$shape, cr = tt$cr, hw = hw)
      lo_dr <- rmax * rel_radius_deriv(profile, zc, section = "lower",
                                       up_exp = tt$up_exp, lo_exp = tt$lo_exp,
                                       shape = tt$shape, cr = tt$cr, hw = hw)
      up_z <- (tt$height - ucl) + zc * ucl
      lo_z <- (tt$height - ucl) - zc * lcl
      up_a <- up_r * sqrt(1 + (up_dr / max(ucl, 1e-9))^2) * (ucl / n_depth) * daz
      lo_a <- lo_r * sqrt(1 + (lo_dr / max(lcl, 1e-9))^2) * (lcl / n_depth) * daz

      rvec <- c(up_r, lo_r); zvec <- c(up_z, lo_z); avec <- c(up_a, lo_a)
      # d(radius)/d(height): negative on the upper surface (radius shrinks as
      # height rises toward the apex), positive on the lower surface.
      slope <- c(up_dr / max(ucl, 1e-9), -lo_dr / max(lcl, 1e-9))
      # facet areas are normalised to the analytic crown surface area, so the
      # discretisation can never make exposed area exceed total area
      avec <- avec * (csa[i] / max(sum(avec) * n_azimuth, .Machine$double.eps))

      fx <- as.vector(outer(rvec, cos(az))) + tt$x
      fy <- as.vector(outer(rvec, sin(az))) + tt$y
      fz <- rep(zvec, times = n_azimuth)
      fa <- rep(avec, times = n_azimuth)
      fs <- rep(slope, times = n_azimuth)
      fth <- rep(az, each = length(rvec))

      others <- p[-i, ]
      if (nrow(others) == 0L) {
        csax[i] <- sum(fa); cpax[i] <- cpa[i]
        if (osv) osvv[i] <- sum(fa)
        next
      }
      zmat <- crown_surface_z(profile, fx, fy, others)
      top <- suppressWarnings(apply(zmat, 1L, max, na.rm = TRUE))
      top[!is.finite(top)] <- -Inf
      exposed <- fz >= top
      csax[i] <- sum(fa[exposed])

      # exposed projection: sample the horizontal disc of the crown
      cpax[i] <- cpa[i] * exposed_projection_fraction(profile, tt, others,
                                                      n_depth, n_azimuth)

      if (osv) {
        ok <- ray_escapes(profile, fx[exposed], fy[exposed], fz[exposed],
                          fth[exposed], fs[exposed], others,
                          ray_steps, ray_length)
        osvv[i] <- sum(fa[exposed][ok])
      }
    }
    data.table::data.table(
      plot_id = p$plot_id, tree_id = p$tree_id,
      csa = csa, csax = csax, csax_rel = sdiv(csax, csa),
      cpa = cpa, cpax = cpax, cpax_rel = sdiv(cpax, cpa),
      cv = cv,
      osv = osvv, osv_rel = sdiv(osvv, csa)
    )
  })
  out <- data.table::rbindlist(res)
  if (!osv) out[, c("osv", "osv_rel") := NULL]
  out[]
}

exposed_projection_fraction <- function(profile, tt, others, n_depth, n_azimuth) {
  rmax <- tt$lcw / 2
  rr <- rmax * sqrt((seq_len(n_depth) - 0.5) / n_depth)   # equal-area rings
  az <- seq(0, 2 * pi, length.out = n_azimuth + 1L)[-(n_azimuth + 1L)]
  px <- as.vector(outer(rr, cos(az))) + tt$x
  py <- as.vector(outer(rr, sin(az))) + tt$y
  hw <- as.numeric(tt$sp_type == "HW")
  cl <- tt$height - tt$hcb
  ucl <- tt$widest * cl
  self <- crown_surface_z(profile, px, py, tt)
  zmat <- crown_surface_z(profile, px, py, others)
  top <- suppressWarnings(apply(zmat, 1L, max, na.rm = TRUE))
  top[!is.finite(top)] <- -Inf
  mean(self[, 1L] >= top, na.rm = TRUE)
}

ray_escapes <- function(profile, fx, fy, fz, fth, fslope, others,
                        ray_steps, ray_length) {
  n <- length(fx)
  if (n == 0L) return(logical(0))
  # outward normal of a surface of revolution: radial component 1, vertical
  # component the negative reciprocal of the profile slope, normalised
  vert <- 1 / sqrt(1 + fslope^2)
  radial <- fslope / sqrt(1 + fslope^2)
  dx <- radial * cos(fth)
  dy <- radial * sin(fth)
  dz <- pmax(vert, 0.05)
  ok <- rep(TRUE, n)
  tt <- seq(ray_length / ray_steps, ray_length, length.out = ray_steps)
  for (s in tt) {
    still <- which(ok)
    if (!length(still)) break
    qx <- fx[still] + s * dx[still]
    qy <- fy[still] + s * dy[still]
    qz <- fz[still] + s * dz[still]
    zmat <- crown_surface_z(profile, qx, qy, others)
    hcb <- rep(others$hcb, each = length(still))
    blocked <- rowSums(!is.na(zmat) &
                         zmat >= matrix(qz, nrow = length(still),
                                        ncol = nrow(others)) &
                         matrix(rep(others$hcb, each = length(still)),
                                nrow = length(still)) <=
                           matrix(qz, nrow = length(still),
                                  ncol = nrow(others))) > 0
    ok[still][blocked] <- FALSE
  }
  ok
}
