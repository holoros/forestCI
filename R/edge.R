#' Flag plot edge trees
#'
#' Distance-dependent indices are biased low for trees near the plot boundary,
#' because part of their competitive neighbourhood was never measured. Two
#' rules are offered.
#'
#' `"buffer"` marks every tree within `buffer` m of the plot boundary. It is
#' simple, and with `buffer` equal to the search radius it guarantees that no
#' retained tree has an unmeasured competitor inside its neighbourhood.
#'
#' `"exterior_angle"` is the rule of Waskiewicz (2011), also used by Kuehne et
#' al. (2019): a sector of `angle` degrees is opened from the subject tree
#' facing directly away from plot centre, and the tree is an edge tree when no
#' measured neighbour falls inside it. It retains more trees than a buffer of
#' the same nominal width, which matters on the small plots typical of
#' stem-mapped inventories, at the cost of a weaker guarantee.
#'
#' @param stand A `stand` with stem coordinates.
#' @param method `"buffer"` or `"exterior_angle"`.
#' @param buffer Buffer width, m, for `method = "buffer"`.
#' @param angle Sector width, degrees, for `method = "exterior_angle"`.
#'   Default 120.
#' @return A `data.table` with `plot_id`, `tree_id`, `is_edge` and, for the
#'   angle rule, `n_exterior`, the count of neighbours in the sector.
#' @export
#' @examples
#' st <- pef_stand()
#' table(edge_flag(st)$is_edge)
edge_flag <- function(stand, method = c("exterior_angle", "buffer"),
                      buffer = 6, angle = 120) {
  stopifnot(inherits(stand, "stand"))
  if (!is_spatial(stand)) stop("this stand has no stem coordinates.", call. = FALSE)
  method <- match.arg(method)
  tr <- stand$trees

  tr[, {
    n <- .N
    if (method == "buffer") {
      r <- plot_radius[1L]
      if (is.na(r)) {
        stop("`plot_radius` is needed for the buffer edge rule.", call. = FALSE)
      }
      list(tree_id = tree_id,
           is_edge = sqrt(x^2 + y^2) > (r - buffer),
           n_exterior = NA_integer_)
    } else {
      half <- angle / 2
      cnt <- integer(n)
      for (i in seq_len(n)) {
        out_bear <- (atan2(x[i], y[i]) * 180 / pi) %% 360
        if (x[i] == 0 && y[i] == 0) out_bear <- 0
        th <- (atan2(x - x[i], y - y[i]) * 180 / pi) %% 360
        dev <- abs(((th - out_bear + 180) %% 360) - 180)
        cnt[i] <- sum(dev <= half & seq_len(n) != i)
      }
      list(tree_id = tree_id, is_edge = cnt == 0L, n_exterior = cnt)
    }
  }, by = "plot_id"]
}
