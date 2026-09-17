# Internal helpers -------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Convert polar field measurements to Cartesian coordinates
#'
#' Converts a distance and a compass bearing measured from plot centre into
#' x (east) and y (north) coordinates. This is the standard field to stem map
#' conversion in forest inventory.
#'
#' @param distance Numeric vector of horizontal distances from plot centre, m.
#' @param azimuth Numeric vector of compass bearings in degrees, 0 to 360,
#'   measured clockwise from north.
#' @param center_x,center_y Numeric scalars or vectors giving the plot centre
#'   coordinates. Default 0.
#'
#' @return A `data.frame` with columns `x` and `y`.
#' @export
#' @examples
#' polar_to_xy(c(5, 10), c(0, 90))
polar_to_xy <- function(distance, azimuth, center_x = 0, center_y = 0) {
  if (!is.numeric(distance) || !is.numeric(azimuth)) {
    stop("`distance` and `azimuth` must be numeric.", call. = FALSE)
  }
  if (any(distance < 0, na.rm = TRUE)) {
    stop("`distance` must be non-negative.", call. = FALSE)
  }
  az <- azimuth %% 360
  rad <- az * pi / 180
  data.frame(
    x = center_x + distance * sin(rad),
    y = center_y + distance * cos(rad)
  )
}

#' Convert Cartesian coordinates to distance and azimuth
#'
#' @param x,y Numeric vectors of coordinates.
#' @param center_x,center_y Plot centre coordinates.
#' @return A `data.frame` with columns `distance` and `azimuth` (degrees).
#' @export
xy_to_polar <- function(x, y, center_x = 0, center_y = 0) {
  dx <- x - center_x
  dy <- y - center_y
  az <- (atan2(dx, dy) * 180 / pi) %% 360
  data.frame(distance = sqrt(dx^2 + dy^2), azimuth = az)
}

# Basal area of a single tree, m2, from DBH in cm
tree_ba <- function(dbh) (pi / 40000) * dbh^2

# Quadratic mean diameter, cm, from DBH (cm) and expansion factor (trees ha-1)
qmd_from <- function(dbh, expf) {
  sqrt(sum(expf * dbh^2, na.rm = TRUE) / sum(expf, na.rm = TRUE))
}

# Pairwise Euclidean distance between two coordinate sets, no matrix blowup
# for the common case of one plot at a time.
pairwise_dist <- function(x, y) {
  as.matrix(stats::dist(cbind(x, y)))
}

# Assert that required columns exist
need_cols <- function(dt, cols, what = "data") {
  missing <- setdiff(cols, names(dt))
  if (length(missing)) {
    stop(sprintf("%s is missing required column(s): %s",
                 what, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

# Safe division: 0 / 0 returns 0 rather than NaN
sdiv <- function(a, b) ifelse(b == 0 | is.na(b), 0, a / b)
