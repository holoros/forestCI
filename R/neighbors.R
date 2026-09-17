#' Select competitor neighbourhoods
#'
#' Builds the subject-competitor pair table every distance-dependent index is
#' computed from. Competitor selection and competitor filtering are kept
#' separate: selection decides which trees are near enough to matter, filtering
#' decides which of those are large enough to matter.
#'
#' Waskiewicz (2011) found that height based filters outperformed canopy
#' position filters, which in turn outperformed diameter based filters, with a
#' threshold near 67 percent of subject crown length. `filter = "height_frac"`
#' with `filter_value = 0.67` reproduces that rule.
#'
#' @param stand A `stand` from [as_stand()] carrying stem coordinates.
#' @param method Competitor selection rule.
#'   `"radius"` takes every tree within `radius` m.
#'   `"knn"` takes the `k` nearest.
#'   `"angle_gauge"` takes trees whose diameter subtends more than the critical
#'   angle of the basal area factor `baf`, that is DBH_j / dist_ij greater than
#'   the gauge constant, which is horizontal point sampling from the subject.
#'   `"influence_zone"` takes trees whose crown influence circle, of radius
#'   `iz_scale` times largest crown radius, overlaps the subject's.
#' @param radius Search radius, m, for `method = "radius"`. Default 6, the
#'   radius used by Kuehne et al. (2019).
#' @param k Number of neighbours for `method = "knn"`. Default 4, following
#'   Motz et al. (2010).
#' @param baf Basal area factor, m2 ha-1, for `method = "angle_gauge"`.
#'   Default 2.3, which is 10 sq ft per acre.
#' @param iz_scale Multiplier on largest crown radius for
#'   `method = "influence_zone"`. Default 1.
#' @param filter Competitor filter.
#'   `"none"` keeps every selected competitor (two-sided index).
#'   `"larger_dbh"` keeps competitors with greater diameter (one-sided).
#'   `"taller"` keeps competitors taller than the subject.
#'   `"height_frac"` keeps competitors whose total height reaches at least
#'   `filter_value` of the way up the subject's crown from its crown base.
#' @param filter_value Threshold for `"height_frac"`. Default 0.67.
#' @param include_self Keep the subject in its own neighbourhood. Needed by
#'   local basal area, not by the ratio indices. Default `FALSE`.
#'
#' @return A `data.table` with one row per subject-competitor pair and columns
#'   `plot_id`, `tree_id`, `comp_id`, `d_ij`, plus the subject and competitor
#'   attributes each index needs.
#' @references
#' Motz, K., Sterba, H., Pommerening, A. (2010) Sampling measures of tree
#' diversity. Forest Ecology and Management 260: 1985-1996.
#' @export
neighbors <- function(stand,
                      method = c("radius", "knn", "angle_gauge", "influence_zone"),
                      radius = 6, k = 4L, baf = 2.3, iz_scale = 1,
                      filter = c("none", "larger_dbh", "taller", "height_frac"),
                      filter_value = 0.67, include_self = FALSE) {
  stopifnot(inherits(stand, "stand"))
  method <- match.arg(method)
  filter <- match.arg(filter)
  if (!is_spatial(stand)) {
    stop("this stand has no stem coordinates, so distance-dependent indices ",
         "cannot be computed. Supply `x` and `y`, or `distance` and `azimuth`, ",
         "to as_stand().", call. = FALSE)
  }
  if (filter %in% c("taller", "height_frac") && !isTRUE(attr(stand, "has_crown"))) {
    stop("filter = \"", filter, "\" needs `height`", 
         if (filter == "height_frac") " and `hcb`" else "",
         " for every tree.", call. = FALSE)
  }
  tr <- stand$trees

  pairs <- tr[, {
    n <- .N
    if (n < 2L) {
      data.table::data.table(tree_id = character(0), comp_id = character(0),
                             d_ij = numeric(0))
    } else {
      dm <- as.matrix(stats::dist(cbind(x, y)))
      diag(dm) <- Inf
      keep <- switch(method,
        radius = dm <= radius,
        knn = {
          kk <- min(k, n - 1L)
          t(apply(dm, 1L, function(r) rank(r, ties.method = "first") <= kk))
        },
        angle_gauge = {
          # gauge constant: DBH (cm) / distance (m) limit for a given BAF
          # metric gauge constant: a tree is tallied when DBH (cm) divided by
          # horizontal distance (m) exceeds 2 * sqrt(BAF)
          gk <- 2 * sqrt(baf)
          outer(rep(1, n), dbh) / dm > gk
        },
        influence_zone = {
          iz <- iz_scale * lcw / 2
          dm < outer(iz, iz, "+")
        }
      )
      idx <- which(keep, arr.ind = TRUE)
      data.table::data.table(
        tree_id = tree_id[idx[, 1L]],
        comp_id = tree_id[idx[, 2L]],
        d_ij    = dm[idx]
      )
    }
  }, by = "plot_id"]

  if (nrow(pairs) == 0L) return(pairs)

  sub <- tr[, c("plot_id", "tree_id", "dbh", "height", "hcb", "ba_m2", "ba_ha",
                "lcw", "mca", "species", "sp_type", "x", "y", "expf",
                "shade"), with = FALSE]
  cmp <- data.table::copy(sub)
  data.table::setnames(cmp, setdiff(names(cmp), "plot_id"),
                       paste0(setdiff(names(cmp), "plot_id"), "_j"))
  data.table::setnames(cmp, "tree_id_j", "comp_id")

  pairs <- merge(pairs, sub, by = c("plot_id", "tree_id"), all.x = TRUE)
  pairs <- merge(pairs, cmp, by = c("plot_id", "comp_id"), all.x = TRUE)

  pairs <- switch(filter,
    none        = pairs,
    larger_dbh  = pairs[pairs$dbh_j > pairs$dbh, ],
    taller      = pairs[pairs$height_j > pairs$height, ],
    height_frac = pairs[pairs$height_j >=
                          pairs$hcb + filter_value * (pairs$height - pairs$hcb), ]
  )

  if (include_self) {
    self <- data.table::copy(sub)
    self[, "comp_id" := self$tree_id]
    self[, "d_ij" := 0]
    for (nm in setdiff(names(sub), c("plot_id", "tree_id"))) {
      self[[paste0(nm, "_j")]] <- self[[nm]]
    }
    pairs <- data.table::rbindlist(list(pairs, self), use.names = TRUE, fill = TRUE)
  }
  data.table::setorderv(pairs, c("plot_id", "tree_id", "d_ij"))
  # the neighbourhood area is what local basal area is expressed per hectare of,
  # and it is only defined for a fixed radius selection
  attr(pairs, "nb_area_ha") <- if (method == "radius") {
    pi * radius^2 / 10000
  } else {
    NA_real_
  }
  attr(pairs, "nb_method") <- method
  pairs[]
}
