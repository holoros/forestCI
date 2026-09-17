#' Compute competition indices
#'
#' The single entry point. Given a stand, it computes every index the data
#' will support and returns one row per tree. Distance-independent indices are
#' always computed. Distance-dependent indices are computed when stem
#' coordinates are present, crown based indices when height and height to
#' crown base are present, and the rest are returned as `NA` columns with one
#' message saying which were skipped and why, rather than silently omitted.
#'
#' @param stand A `stand` from [as_stand()].
#' @param indices Character vector of index families to compute, any of
#'   `"independent"`, `"dependent"`, `"apa"`, `"rapa"`, `"crown"`,
#'   `"pattern"`. Default all.
#' @param profile A `crown_profile` for the crown based indices.
#' @param radius Search radius, m, for the pairwise indices. Default 6.
#' @param filter Competitor filter passed to [neighbors()].
#' @param edge Edge handling. `"flag"` adds the `is_edge` column and leaves the
#'   values in place. `"exclude"` additionally sets every distance-dependent
#'   value to `NA` for edge trees, which is the conservative choice and the one
#'   Kuehne et al. (2019) made. `"none"` skips the edge computation.
#' @param edge_method,edge_buffer,edge_angle Passed to [edge_flag()].
#' @param apa_weight,apa_exponent Weighting variable and exponent for [apa()].
#' @param rapa_exponent Exponent for [rapa()]. Kept separate because the
#'   rasterised construction is multiplicatively weighted, so an exponent that
#'   is mild in [apa()] is aggressive here. Default 1.
#' @param rapa_resolution Cell size for [rapa()], m.
#' @param crown_osv Compute open sky view, the expensive crown index.
#' @param n_depth,n_azimuth Crown facet grid, passed to [crown_exposure()].
#' @param quiet Suppress the message listing skipped families.
#' @param ... Passed to [neighbors()].
#'
#' @return A `data.table` with one row per tree, carrying the stand columns and
#'   every computed index.
#' @export
#' @examples
#' st <- pef_stand()
#' ci <- competition_indices(st, indices = c("independent", "dependent"))
#' str(ci)
competition_indices <- function(stand,
                                indices = c("independent", "dependent", "apa",
                                            "rapa", "crown", "pattern"),
                                profile = crown_profile(),
                                radius = 6,
                                filter = "none",
                                edge = c("flag", "exclude", "none"),
                                edge_method = "exterior_angle",
                                edge_buffer = radius, edge_angle = 120,
                                apa_weight = "dbh", apa_exponent = 1,
                                rapa_exponent = 1, rapa_resolution = 0.25,
                                crown_osv = TRUE,
                                n_depth = 18L, n_azimuth = 72L,
                                quiet = FALSE, ...) {
  stopifnot(inherits(stand, "stand"))
  indices <- match.arg(indices, several.ok = TRUE)
  edge <- match.arg(edge)

  spatial <- is_spatial(stand)
  crowned <- isTRUE(attr(stand, "has_crown"))
  skipped <- character(0)

  out <- data.table::copy(stand$trees)
  dep_cols <- character(0)

  if ("independent" %in% indices) {
    di <- ci_distance_independent(stand)
    out <- merge(out, di, by = c("plot_id", "tree_id"), all.x = TRUE,
                 suffixes = c("", ".di"))
  }

  if ("dependent" %in% indices) {
    if (spatial) {
      nb <- neighbors(stand, method = "radius", radius = radius,
                      filter = filter, ...)
      dd <- ci_distance_dependent(stand, nb = nb)
      dep_cols <- c(dep_cols, setdiff(names(dd), c("plot_id", "tree_id")))
      out <- merge(out, dd, by = c("plot_id", "tree_id"), all.x = TRUE)
    } else {
      skipped <- c(skipped, "dependent (no stem coordinates)")
    }
  }

  if ("apa" %in% indices) {
    if (spatial) {
      a <- apa(stand, weight = apa_weight, exponent = apa_exponent)
      dep_cols <- c(dep_cols, "apa", "apa_share")
      out <- merge(out, a, by = c("plot_id", "tree_id"), all.x = TRUE)
    } else {
      skipped <- c(skipped, "apa (no stem coordinates)")
    }
  }

  if ("rapa" %in% indices) {
    if (spatial) {
      a <- rapa(stand, weight = apa_weight, exponent = rapa_exponent,
                resolution = rapa_resolution)
      dep_cols <- c(dep_cols, "rapa", "rapa_share", "n_cells")
      out <- merge(out, a, by = c("plot_id", "tree_id"), all.x = TRUE)
    } else {
      skipped <- c(skipped, "rapa (no stem coordinates)")
    }
  }

  if ("crown" %in% indices) {
    if (spatial && crowned) {
      ce <- crown_exposure(stand, profile = profile, n_depth = n_depth,
                           n_azimuth = n_azimuth, osv = crown_osv)
      dep_cols <- c(dep_cols, setdiff(names(ce), c("plot_id", "tree_id")))
      out <- merge(out, ce, by = c("plot_id", "tree_id"), all.x = TRUE)
    } else {
      skipped <- c(skipped, paste0("crown (",
        if (!spatial) "no stem coordinates" else "no height or hcb", ")"))
    }
  }

  if ("pattern" %in% indices) {
    if (spatial) {
      sp <- spatial_pattern(stand)
      out <- merge(out, sp, by = "plot_id", all.x = TRUE)
    } else {
      skipped <- c(skipped, "pattern (no stem coordinates)")
    }
  }

  if (edge != "none" && spatial) {
    ef <- edge_flag(stand, method = edge_method, buffer = edge_buffer,
                    angle = edge_angle)
    out <- merge(out, ef, by = c("plot_id", "tree_id"), all.x = TRUE)
    if (edge == "exclude" && length(dep_cols)) {
      dep_cols <- intersect(unique(dep_cols), names(out))
      for (cl in dep_cols) out[out$is_edge %in% TRUE, (cl) := NA]
    }
  }

  if (length(skipped) && !quiet) {
    message("forestCI: skipped ", paste(skipped, collapse = "; "), ".")
  }
  data.table::setorderv(out, c("plot_id", "tree_id"))
  out[]
}
