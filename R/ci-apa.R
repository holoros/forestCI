# Convex polygon clipping (Sutherland-Hodgman) against a half plane
# a*x + b*y <= c. Polygon given as a two column matrix, vertices ordered.
clip_halfplane <- function(poly, a, b, c) {
  n <- nrow(poly)
  if (n == 0L) return(poly)
  inside <- a * poly[, 1L] + b * poly[, 2L] <= c + 1e-12
  if (all(inside)) return(poly)
  if (!any(inside)) return(poly[0L, , drop = FALSE])
  out <- matrix(NA_real_, nrow = 0L, ncol = 2L)
  for (i in seq_len(n)) {
    j <- if (i == n) 1L else i + 1L
    pi_ <- poly[i, ]; pj <- poly[j, ]
    if (inside[i]) out <- rbind(out, pi_)
    if (inside[i] != inside[j]) {
      di <- a * pi_[1L] + b * pi_[2L] - c
      dj <- a * pj[1L] + b * pj[2L] - c
      t <- di / (di - dj)
      out <- rbind(out, pi_ + t * (pj - pi_))
    }
  }
  out
}

polygon_area <- function(poly) {
  n <- nrow(poly)
  if (n < 3L) return(0)
  x <- poly[, 1L]; y <- poly[, 2L]
  abs(sum(x * c(y[-1L], y[1L]) - c(x[-1L], x[1L]) * y)) / 2
}

regular_polygon <- function(cx, cy, r, n = 72L) {
  th <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)]
  cbind(cx + r * cos(th), cy + r * sin(th))
}

#' Area potentially available
#'
#' Growing space allocated to each tree by a weighted Voronoi tessellation of
#' the plot. Each pair of trees is separated by a line perpendicular to the
#' segment joining them; with `weight = "none"` the line is the perpendicular
#' bisector and the result is the ordinary Voronoi polygon, and with a weight
#' it is displaced toward the smaller tree in proportion to relative size,
#' which is the construction of Moore, Budelsky and Schlesinger (1973).
#'
#' The displacement is at a distance from the subject of
#' `d_ij * w_i / (w_i + w_j)`, with `w` the weighting variable raised to
#' `exponent`. With `weight = "dbh"` and `exponent = 2` this is the original
#' Moore et al. form, which assumes partial size asymmetry with respect to
#' diameter and perfect size symmetry with respect to basal area.
#'
#' @section Weighted polygons do not tile exactly:
#' With `weight = "none"` the polygons partition the plot and their areas sum
#' to the plot area to within rounding. With a weight they do not. Displacing
#' each bisector independently is an approximation to a weighted tessellation,
#' not a tessellation, so neighbouring polygons can overlap slightly or leave
#' slivers, and on the example plots the total departs from the plot area by
#' about 24 percent at `exponent = 2`. This is a property of the construction,
#' not of the implementation. The returned `apa_share` makes the departure
#' visible, and `normalize = TRUE` rescales each plot's polygons to sum to the
#' plot area if a conserved quantity is wanted, at the cost of making a tree's
#' value depend on its neighbours' totals.
#'
#' Polygons are clipped to the plot boundary, so no tree is allocated growing
#' space outside the plot. Trees on the plot edge have no competitor on their
#' outer side and their polygon is therefore bounded by the plot boundary
#' rather than by a neighbour; `edge_flag()` marks them and
#' `competition_indices()` sets their APA to `NA` when `edge = "exclude"`.
#'
#' @param stand A `stand` carrying stem coordinates.
#' @param weight Weighting variable: `"none"`, `"dbh"`, `"lcw"`, `"csa"`, or
#'   the name of any numeric column of `stand$trees`.
#' @param exponent Exponent applied to the weighting variable. Waskiewicz
#'   (2011) found an optimum of 1 for diameter weighting and 2 for crown
#'   surface width weighting.
#' @param boundary Plot boundary. `"circle"` uses `plot_radius`, `"convex"`
#'   uses the convex hull of the stems, `"square"` uses a square of half width
#'   `extent` centred on plot centre, and `"none"` uses a square sized to
#'   contain every stem.
#' @param extent Half width, m, of the square domain when
#'   `boundary = "square"`.
#' @param n_vertices Vertices used to approximate a circular boundary.
#' @param normalize Rescale each plot's polygon areas so that they sum to the
#'   plot area. Default `FALSE`.
#'
#' @return A `data.table` with `plot_id`, `tree_id`, `apa` (m2) and
#'   `apa_share`, the polygon area as a proportion of plot area.
#' @references
#' Moore, J.A., Budelsky, C.A., Schlesinger, R.C. (1973) A new index
#' representing individual tree competitive status. Canadian Journal of Forest
#' Research 3: 495-500.
#' @export
#' @examples
#' st <- pef_stand()
#' head(apa(st, weight = "dbh", exponent = 2))
apa <- function(stand, weight = c("none", "dbh", "lcw", "csa"), exponent = 1,
                boundary = c("circle", "convex", "square", "none"),
                n_vertices = 180L, extent = 25, normalize = FALSE) {
  stopifnot(inherits(stand, "stand"))
  if (!is_spatial(stand)) stop("this stand has no stem coordinates.", call. = FALSE)
  boundary <- match.arg(boundary)
  wname <- if (length(weight) > 1L) "none" else as.character(weight)
  tr <- stand$trees

  tr[, {
    n <- .N
    w <- if (wname == "none") rep(1, n) else {
      if (!wname %in% names(.SD) && !wname %in% names(tr)) {
        stop("weight column `", wname, "` not found in the stand.", call. = FALSE)
      }
      v <- if (wname %in% names(.SD)) .SD[[wname]] else tr[[wname]][.I]
      as.numeric(v)^exponent
    }
    base_poly <- switch(boundary,
      circle = {
        r <- plot_radius[1L]
        if (is.na(r)) regular_polygon(mean(x), mean(y),
                                      max(sqrt(x^2 + y^2)) + 1, n_vertices)
        else regular_polygon(0, 0, r, n_vertices)
      },
      convex = {
        h <- grDevices::chull(x, y)
        cbind(x[h], y[h])
      },
      square = cbind(c(-extent, extent, extent, -extent),
                     c(-extent, -extent, extent, extent)),
      none = {
        r <- max(abs(c(x, y))) * 2
        cbind(c(-r, r, r, -r), c(-r, -r, r, r))
      }
    )
    area <- numeric(n)
    wmax <- max(w)
    for (i in seq_len(n)) {
      poly <- base_poly
      dall <- sqrt((x - x[i])^2 + (y - y[i])^2)
      ord <- order(dall); ord <- ord[ord != i]
      # The bisector against neighbour j sits at frac_j * d_ij from the subject,
      # and frac_j is at least w_i / (w_i + max(w)). Working outward from the
      # nearest neighbour, once that lower bound exceeds the polygon's own
      # circumradius no remaining neighbour can cut it, so the loop stops. This
      # is an exact early exit, not an approximation, and it is what keeps the
      # cost near linear instead of quadratic in stand size.
      fmin <- w[i] / (w[i] + wmax)
      for (j in ord) {
        dij <- dall[j]
        if (dij == 0) next
        rpoly <- if (nrow(poly) >= 1L) {
          max(sqrt((poly[, 1L] - x[i])^2 + (poly[, 2L] - y[i])^2))
        } else 0
        if (fmin * dij > rpoly) break
        dx <- x[j] - x[i]; dy <- y[j] - y[i]
        frac <- w[i] / (w[i] + w[j])
        px <- x[i] + frac * dx
        py <- y[i] + frac * dy
        # half plane keeping the side containing tree i
        a <- dx / dij; b <- dy / dij
        cc <- a * px + b * py
        poly <- clip_halfplane(poly, a, b, cc)
        if (nrow(poly) < 3L) break
      }
      area[i] <- polygon_area(poly)
    }
    plot_m2 <- plot_area[1L] * 10000
    if (normalize && sum(area) > 0) area <- area * plot_m2 / sum(area)
    list(tree_id = tree_id, apa = area, apa_share = area / plot_m2)
  }, by = "plot_id"]
}

#' Rasterised area potentially available
#'
#' The rasterised counterpart of [apa()]. The plot is discretised on a square
#' grid and every cell is assigned to the tree minimising `d / w`, the
#' multiplicatively weighted Voronoi assignment of Mu (2004). Unlike the
#' polygon construction this allows a polygon to bend around and to nest
#' completely inside another, which is what happens when a small tree stands
#' under a large one, and it is the reason Kuehne et al. (2019) preferred the
#' rasterised form.
#'
#' @section Weighting is not the same operation as in [apa()]:
#' [apa()] displaces a straight bisector to a point `w_i / (w_i + w_j)` of the
#' way along the segment joining two trees, so the weight enters as a bounded
#' ratio and the boundary can never move past the competitor. The rasterised
#' construction is multiplicatively weighted: the boundary is the locus where
#' `d_i / w_i` equals `d_j / w_j`, so a tree with twice the weight reaches twice
#' as far and a large enough tree can swallow a neighbour entirely. That
#' nesting is the point of the rasterised form, and it is why Kuehne et al.
#' (2019) used it, but it means the two functions agree only when unweighted
#' and that an exponent above 1 is far more aggressive here than in [apa()].
#' `weight_mode = "additive"` gives the milder alternative in which the
#' boundary is where `d_i - w_i` equals `d_j - w_j`.
#'
#' `dr_mod = TRUE` adds the directional modification: each tree's weight is
#' scaled between 1 and `dr_max` as a function of bearing, reaching `dr_max`
#' toward `dr_bearing`, so that a tree competes more strongly on one side. The
#' default bearing is downhill, that is plot aspect plus 180 degrees, scaled by
#' plot slope. This is the package's implementation of the idea described by
#' Kuehne et al. (2019), not a ported equation, and it reduces exactly to the
#' unmodified assignment when `dr_max` is 1 or slope is 0.
#'
#' @inheritParams apa
#' @param resolution Cell size, m. Default 0.25. Halving it quadruples the
#'   cost and changes APA in the third significant figure.
#' @param align Grid alignment. `"center"` samples cell centres, so the sampled
#'   points tile the domain exactly and the areas sum to it. `"node"` samples
#'   grid nodes from `-extent` to `extent` inclusive, which is the convention of
#'   the original Acadian rasterised APA code and counts the boundary ring of
#'   half cells at full weight. Use `"node"` only to reproduce that code.
#' @param dr_mod Apply the directional modification.
#' @param dr_max Maximum directional multiplier on the weight. Default 2.
#' @param dr_bearing Bearing, degrees, of maximum competitive reach. `NULL`
#'   uses downhill.
#' @param weight_mode `"multiplicative"` assigns each cell to the tree
#'   minimising `d / w`, which is the weighted Voronoi diagram of Mu (2004).
#'   `"additive"` minimises `d - w`, which is milder and never lets one tree
#'   nest inside another.
#' @return A `data.table` with `plot_id`, `tree_id`, `rapa` (m2),
#'   `rapa_share` and `n_cells`.
#' @references
#' Mu, L. (2004) Polygon characterization with the multiplicatively weighted
#' Voronoi diagram. The Professional Geographer 56: 223-239.
#' @export
rapa <- function(stand, weight = c("none", "dbh", "lcw", "csa"), exponent = 1,
                 resolution = 0.25,
                 boundary = c("circle", "convex", "square", "none"),
                 extent = 25, align = c("center", "node"),
                 dr_mod = FALSE, dr_max = 2, dr_bearing = NULL,
                 weight_mode = c("multiplicative", "additive")) {
  align <- match.arg(align)
  stopifnot(inherits(stand, "stand"))
  if (!is_spatial(stand)) stop("this stand has no stem coordinates.", call. = FALSE)
  boundary <- match.arg(boundary)
  weight_mode <- match.arg(weight_mode)
  wname <- if (length(weight) > 1L) "none" else as.character(weight)
  tr <- stand$trees

  tr[, {
    n <- .N
    w <- if (wname == "none") rep(1, n) else as.numeric(.SD[[wname]])^exponent
    # scale weights to a mean of one so that the additive mode is on the same
    # length scale as the grid, and so the multiplicative mode is invariant to
    # the units of the weighting variable
    w <- w / mean(w)
    r <- if (boundary == "square") extent else plot_radius[1L]
    if (is.na(r)) r <- max(sqrt(x^2 + y^2)) + 1
    gx <- if (align == "node") {
      seq(-r, r, by = resolution)
    } else {
      seq(-r + resolution / 2, r - resolution / 2, by = resolution)
    }
    gy <- gx
    g <- expand.grid(gx = gx, gy = gy)
    keep <- switch(boundary,
      circle = sqrt(g$gx^2 + g$gy^2) <= r,
      convex = {
        h <- grDevices::chull(x, y)
        point_in_poly(g$gx, g$gy, cbind(x[h], y[h]))
      },
      square = rep(TRUE, nrow(g)),
      none = rep(TRUE, nrow(g))
    )
    g <- g[keep, , drop = FALSE]
    if (nrow(g) == 0L) {
      return(list(tree_id = tree_id, rapa = rep(NA_real_, n),
                  rapa_share = rep(NA_real_, n), n_cells = rep(0L, n)))
    }
    bearing <- if (!is.null(dr_bearing)) {
      rep_len(dr_bearing, n)
    } else {
      rep_len((aspect[1L] + 180) %% 360, n)
    }
    slope_scale <- if (is.null(dr_bearing)) min(slope[1L] / 100, 1) else 1

    best <- rep(Inf, nrow(g)); owner <- rep(NA_integer_, nrow(g))
    for (i in seq_len(n)) {
      dx <- g$gx - x[i]; dy <- g$gy - y[i]
      d <- sqrt(dx^2 + dy^2)
      wi <- w[i]
      if (dr_mod && dr_max > 1 && slope_scale > 0) {
        th <- (atan2(dx, dy) * 180 / pi) %% 360
        m <- 1 + (dr_max - 1) * slope_scale *
          (1 + cos((th - bearing[i]) * pi / 180)) / 2
        wi <- wi * m
      }
      s <- if (weight_mode == "multiplicative") d / wi else d - wi
      upd <- s < best
      best[upd] <- s[upd]; owner[upd] <- i
    }
    cnt <- tabulate(owner, nbins = n)
    area <- cnt * resolution^2
    list(tree_id = tree_id, rapa = area,
         rapa_share = area / (plot_area[1L] * 10000),
         n_cells = cnt)
  }, by = "plot_id"]
}

# Ray casting point in polygon, vectorised over points
point_in_poly <- function(px, py, poly) {
  n <- nrow(poly)
  inside <- rep(FALSE, length(px))
  j <- n
  for (i in seq_len(n)) {
    xi <- poly[i, 1L]; yi <- poly[i, 2L]
    xj <- poly[j, 1L]; yj <- poly[j, 2L]
    cross <- ((yi > py) != (yj > py)) &
      (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
    inside <- xor(inside, cross)
    j <- i
  }
  inside
}
