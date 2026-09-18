#' Distance-dependent competition indices
#'
#' Computes the spatially explicit pairwise family from a neighbourhood built
#' by [neighbors()]. Every index here is a sum over competitors of a function
#' of relative size and distance, so they differ only in how sharply they
#' discount distance and how they scale size.
#'
#' @param stand A `stand` carrying stem coordinates.
#' @param nb A pair table from [neighbors()]. Built with the defaults if
#'   absent.
#' @param indices Character vector naming the indices to compute. See Details.
#' @param nci_alpha,nci_beta,nci_lambda Coefficients of the Canham
#'   neighbourhood competition index. `nci_lambda` is either a single number
#'   applied to every competitor or a named numeric vector indexed by
#'   competitor species code. Supplying values is the caller's responsibility;
#'   the package ships no fitted NCI parameters and will not invent any.
#' @param ... Passed to [neighbors()] when `nb` is not supplied.
#'
#' @section Indices:
#' \describe{
#'   \item{`hegyi`}{sum over j of (DBH_j / DBH_i) / d_ij. Hegyi (1974).}
#'   \item{`hegyi_ba`}{The basal area ratio variant used by Kuehne et al.
#'     (2019), sum of (BA_j / BA_i) / d_ij.}
#'   \item{`martin_ek`}{sum of (DBH_j / DBH_i) * exp(-16 * d_ij /
#'     (DBH_i + DBH_j)). Martin and Ek (1984).}
#'   \item{`daniels`}{sum of DBH_j^2 / (DBH_i^2 * d_ij). Daniels (1976).}
#'   \item{`lorimer`}{sum of DBH_j / DBH_i, no distance weight, over a fixed
#'     radius. Lorimer (1983).}
#'   \item{`rouvinen`}{sum of arctan(DBH_j / (100 * d_ij)), the angular size of
#'     each competitor. Rouvinen and Kuuluvainen (1997).}
#'   \item{`spurr`}{Point density. Competitors are ranked on DBH_j / d_ij, and
#'     the index is the mean over ranks k of 0.25 * (k - 0.5) *
#'     (DBH_k / d_k)^2 expressed per hectare. Spurr (1962).}
#'   \item{`nci`}{sum of lambda * DBH_j^alpha / d_ij^beta. Canham et al.
#'     (2004). Returns `NA` unless coefficients are supplied.}
#'   \item{`local_ba`}{Basal area of the neighbourhood including the subject,
#'     expressed per hectare of the neighbourhood itself rather than per
#'     hectare of the plot, so that it is on the same scale as stand basal
#'     area and comparable across search radii. Kuehne et al. (2019). Defined
#'     for `method = "radius"`; with any other selection rule it falls back to
#'     the plot expansion factor and says so.}
#'   \item{`local_bal`}{The same restricted to competitors larger than the
#'     subject.}
#'   \item{`mean_dist`}{Mean distance to the selected competitors. With
#'     `method = "knn"` and `k = 4` this is the meanDIST of Motz et al. (2010).}
#'   \item{`n_comp`}{Competitor count, the simplest possible index and a useful
#'     diagnostic of whether the neighbourhood is empty.}
#' }
#'
#' @return A `data.table` with `plot_id`, `tree_id` and one column per index.
#' @references
#' Daniels, R.F. (1976) Simple competition indices and their correlation with
#' annual loblolly pine tree growth. Forest Science 22: 454-456.
#'
#' Hegyi, F. (1974) A simulation model for managing jack pine stands. In:
#' Fries, J. (ed.) Growth Models for Tree and Stand Simulation. Royal College
#' of Forestry, Stockholm, pp. 74-90.
#'
#' Lorimer, C.G. (1983) Tests of age-independent competition indices for
#' individual trees in natural hardwood stands. Forest Ecology and Management
#' 6: 343-360.
#'
#' Martin, G.L., Ek, A.R. (1984) A comparison of competition measures and
#' growth models for predicting plantation red pine diameter and height
#' growth. Forest Science 30: 731-743.
#'
#' Rouvinen, S., Kuuluvainen, T. (1997) Structure and asymmetry of tree crowns
#' in relation to local competition in a natural mature Scots pine forest.
#' Canadian Journal of Forest Research 27: 890-902.
#'
#' Spurr, S.H. (1962) A measure of point density. Forest Science 8: 85-96.
#' @export
#' @examples
#' st <- pef_stand()
#' head(ci_distance_dependent(st, radius = 6))
ci_distance_dependent <- function(stand, nb = NULL,
                                  indices = c("hegyi", "hegyi_ba", "martin_ek",
                                              "daniels", "lorimer", "rouvinen",
                                              "spurr", "local_ba", "local_bal",
                                              "mean_dist", "n_comp"),
                                  nci_alpha = NULL, nci_beta = NULL,
                                  nci_lambda = NULL, ...) {
  stopifnot(inherits(stand, "stand"))
  indices <- match.arg(indices, choices = c(
    "hegyi", "hegyi_ba", "martin_ek", "daniels", "lorimer", "rouvinen",
    "spurr", "nci", "local_ba", "local_bal", "mean_dist", "n_comp"),
    several.ok = TRUE)
  if (is.null(nb)) nb <- neighbors(stand, ...)
  nb_area <- attr(nb, "nb_area_ha")

  base <- stand$trees[, c("plot_id", "tree_id"), with = FALSE]
  if (nrow(nb) == 0L) {
    # No tree in the stand has a competitor. That is zero competition, not
    # unknown competition, and it must agree with what a single tree with an
    # empty neighbourhood gets when other trees do have one. Only mean distance
    # to a competitor is genuinely undefined.
    for (i in indices) {
      base[, (i) := if (identical(i, "mean_dist")) NA_real_ else 0]
    }
    if ("local_ba" %in% indices) {
      nb_area <- attr(nb, "nb_area_ha")
      own <- stand$trees$ba_m2
      base[, "local_ba" := if (is.null(nb_area) || is.na(nb_area)) {
        own * stand$trees$expf
      } else own / nb_area]
    }
    return(base[])
  }
  nb_area <- attr(nb, "nb_area_ha")
  if (is.null(nb_area)) nb_area <- NA_real_
  if (is.na(nb_area) && any(c("local_ba", "local_bal") %in% indices)) {
    message("forestCI: local basal area is expressed per hectare of the ",
            "neighbourhood, which is only defined for a fixed radius. With ",
            "method \"", attr(nb, "nb_method") %||% "unknown",
            "\" it falls back to the plot expansion factor.")
  }
  nb <- nb[nb$comp_id != nb$tree_id, ]
  nb[, "d_ij" := pmax(d_ij, .Machine$double.eps)]

  agg <- nb[, {
    dbh_i <- dbh[1L]; ba_i <- ba_m2[1L]
    res <- list()
    if ("hegyi"     %in% indices) res$hegyi     <- sum((dbh_j / dbh_i) / d_ij)
    if ("hegyi_ba"  %in% indices) res$hegyi_ba  <- sum((ba_m2_j / ba_i) / d_ij)
    if ("martin_ek" %in% indices) res$martin_ek <-
      sum((dbh_j / dbh_i) * exp(-16 * d_ij / (dbh_i + dbh_j)))
    if ("daniels"   %in% indices) res$daniels   <- sum(dbh_j^2 / (dbh_i^2 * d_ij))
    if ("lorimer"   %in% indices) res$lorimer   <- sum(dbh_j / dbh_i)
    if ("rouvinen"  %in% indices) res$rouvinen  <- sum(atan(dbh_j / (100 * d_ij)))
    if ("spurr"     %in% indices) {
      o <- order(-(dbh_j / d_ij))
      kk <- seq_along(o)
      res$spurr <- sum(0.25 * (kk - 0.5) * (dbh_j[o] / d_ij[o])^2) / length(o)
    }
    if ("nci" %in% indices) {
      res$nci <- if (is.null(nci_alpha) || is.null(nci_beta)) {
        NA_real_
      } else {
        lam <- if (is.null(nci_lambda)) {
          1
        } else if (length(nci_lambda) == 1L) {
          as.numeric(nci_lambda)
        } else {
          l <- nci_lambda[as.character(species_j)]
          l[is.na(l)] <- 0
          as.numeric(l)
        }
        sum(lam * dbh_j^nci_alpha / d_ij^nci_beta)
      }
    }
    if ("local_ba" %in% indices) {
      res$local_ba <- if (is.na(nb_area)) sum(ba_ha_j) + ba_ha[1L]
                      else (sum(ba_m2_j) + ba_m2[1L]) / nb_area
    }
    if ("local_bal" %in% indices) {
      res$local_bal <- if (is.na(nb_area)) sum(ba_ha_j[dbh_j > dbh_i])
                       else sum(ba_m2_j[dbh_j > dbh_i]) / nb_area
    }
    if ("mean_dist" %in% indices) res$mean_dist <- mean(d_ij)
    if ("n_comp"    %in% indices) res$n_comp    <- as.numeric(length(d_ij))
    res
  }, by = c("plot_id", "tree_id")]

  out <- merge(base, agg, by = c("plot_id", "tree_id"), all.x = TRUE)
  # a tree with an empty neighbourhood has zero competition, not missing
  # competition, for every index that is a sum. mean_dist stays NA.
  sums <- setdiff(intersect(indices, names(out)), "mean_dist")
  for (i in sums) out[is.na(get(i)), (i) := 0]
  if ("local_ba" %in% names(out)) {
    # a tree with no competitor still occupies its own neighbourhood
    empty <- out$local_ba == 0
    if (any(empty)) {
      own <- stand$trees$ba_m2[match(paste(out$plot_id[empty], out$tree_id[empty]),
                                     paste(stand$trees$plot_id, stand$trees$tree_id))]
      out[empty, "local_ba" := if (is.na(nb_area)) {
        own * stand$trees$expf[1L]
      } else own / nb_area]
    }
  }
  out[]
}

#' Plot level spatial pattern statistics
#'
#' Clark and Evans (1954) aggregation index R and the mean directional index of
#' Corral-Rivas et al. (2010), the two plot level descriptors Kuehne et al.
#' (2019) used to stratify their results.
#'
#' R is the mean nearest neighbour distance divided by the expectation under
#' complete spatial randomness. R of 1 is random, above 1 regular, below 1
#' clustered. The Donnelly edge correction is applied for circular plots when
#' `plot_radius` is available.
#'
#' The mean directional index summarises the angular evenness of the `k`
#' nearest neighbours around each subject: it is the length of the resultant of
#' the k unit vectors pointing at those neighbours, averaged over trees, and
#' runs from 0 for perfectly even spacing to `k` for all neighbours on one side.
#'
#' @param stand A `stand` carrying stem coordinates.
#' @param k Neighbours used by the mean directional index. Default 4.
#' @param edge_correction Apply the Donnelly correction to R. Default `TRUE`.
#' @return A `data.table` with `plot_id`, `clark_evans`, `mdi`,
#'   `mean_nn_dist` and `n_trees`.
#' @references
#' Clark, P.J., Evans, F.C. (1954) Distance to nearest neighbor as a measure of
#' spatial relationships in populations. Ecology 35: 445-453.
#'
#' Corral-Rivas, J.J., Wehenkel, C., Castellanos-Bocaz, H.A., Vargas-Larreta,
#' B., Dieguez-Aranda, U. (2010) A permutation test of spatial randomness.
#' Journal of Forest Research 15: 218-225.
#' @export
spatial_pattern <- function(stand, k = 4L, edge_correction = TRUE) {
  stopifnot(inherits(stand, "stand"))
  if (!is_spatial(stand)) {
    stop("this stand has no stem coordinates.", call. = FALSE)
  }
  tr <- stand$trees
  tr[, {
    n <- .N
    if (n < 3L) {
      list(clark_evans = NA_real_, mdi = NA_real_,
           mean_nn_dist = NA_real_, n_trees = n)
    } else {
      dm <- as.matrix(stats::dist(cbind(x, y)))
      diag(dm) <- Inf
      nn <- apply(dm, 1L, min)
      area_m2 <- plot_area[1L] * 10000
      lambda <- n / area_m2
      r_exp <- 1 / (2 * sqrt(lambda))
      if (edge_correction && !is.na(plot_radius[1L])) {
        # Donnelly (1978) correction for a circular plot
        per <- 2 * pi * plot_radius[1L]
        r_exp <- r_exp + (0.0514 + 0.041 / sqrt(n)) * per / n
      }
      kk <- min(k, n - 1L)
      res <- vapply(seq_len(n), function(i) {
        j <- order(dm[i, ])[seq_len(kk)]
        ang <- atan2(y[j] - y[i], x[j] - x[i])
        sqrt(sum(cos(ang))^2 + sum(sin(ang))^2)
      }, numeric(1))
      list(clark_evans = mean(nn) / r_exp,
           mdi = mean(res),
           mean_nn_dist = mean(nn),
           n_trees = n)
    }
  }, by = "plot_id"]
}
