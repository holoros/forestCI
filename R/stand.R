#' Build a stand object
#'
#' Standardises an inventory table into the object every `forestCI` function
#' consumes. Column names are supplied once here and nowhere else. Stem
#' coordinates may be given directly as `x` and `y` or as `distance` and
#' `azimuth` from plot centre, and if neither is available the object is still
#' valid: the distance-independent indices work without coordinates, and the
#' distance-dependent ones report that they were not computable rather than
#' failing.
#'
#' @param data A `data.frame` of tree records, one row per tree per
#'   measurement.
#' @param plot,tree,species,dbh Column names for plot identifier, tree
#'   identifier, species code and diameter at breast height in cm. Required.
#' @param sp_type Optional column name holding `"SW"` or `"HW"` per tree. It is
#'   used only for species codes the trait table does not recognise, where it
#'   decides whether the softwood or the hardwood default applies. Without it
#'   an unrecognised code falls back to the hardwood default. Recognised codes
#'   take their type from the trait table and ignore this column.
#' @param height,hcb Column names for total height and height to crown base,
#'   both m. Optional. Crown based indices need them; supply them or supply
#'   `crown_ratio`.
#' @param crown_ratio Column name for crown ratio, m m-1. Used to derive `hcb`
#'   when `hcb` is absent.
#' @param expf Column name for the expansion factor, trees ha-1 represented by
#'   each record. If absent it is derived from `plot_area`.
#' @param x,y Column names for stem coordinates, m.
#' @param distance,azimuth Column names for horizontal distance, m, and
#'   compass bearing, degrees, from plot centre. Used when `x` and `y` are
#'   absent.
#' @param plot_area Column name for plot area in ha, or a single numeric value
#'   applied to every plot.
#' @param plot_radius Column name for plot radius in m, or a single numeric
#'   value. Used for edge correction and, for circular plots, to derive
#'   `plot_area`.
#' @param slope,aspect Column names for plot slope in percent and aspect in
#'   degrees, or single numeric values. Default 0.
#' @param year Optional column name for measurement year, so that a
#'   remeasured plot is treated as separate stands.
#' @param traits Species trait table, from [species_traits()].
#' @param quiet Suppress the message about unmatched species codes.
#'
#' @return An object of class `stand`: a list with element `trees`
#'   (a `data.table`) and element `plots` (a `data.table`), plus attributes
#'   recording whether coordinates and crown data are present.
#' @export
#' @examples
#' data(pef)
#' st <- as_stand(pef, plot = "plot_id", tree = "tree_id", species = "species",
#'                dbh = "dbh", height = "height", hcb = "hcb", expf = "expf",
#'                distance = "distance", azimuth = "azimuth",
#'                plot_radius = 16.1)
#' st
as_stand <- function(data,
                     plot = "plot", tree = "tree", species = "species",
                     dbh = "dbh", sp_type = NULL, height = NULL, hcb = NULL,
                     crown_ratio = NULL, expf = NULL,
                     x = NULL, y = NULL,
                     distance = NULL, azimuth = NULL,
                     plot_area = NULL, plot_radius = NULL,
                     slope = 0, aspect = 0, year = NULL,
                     traits = species_traits(), quiet = FALSE) {

  d <- data.table::as.data.table(data)
  need_cols(d, c(plot, tree, species, dbh), "`data`")

  tr <- data.table::data.table(
    plot_id = as.character(d[[plot]]),
    tree_id = as.character(d[[tree]]),
    species = as.character(d[[species]]),
    dbh     = as.numeric(d[[dbh]])
  )
  if (!is.null(sp_type)) {
    need_cols(d, sp_type, "`data`")
    st <- as.character(d[[sp_type]])
    if (!all(st %in% c("SW", "HW"))) {
      stop("`sp_type` must be \"SW\" or \"HW\" for every record.", call. = FALSE)
    }
    tr[, "sp_type" := st]
  }
  if (!is.null(year)) {
    tr[, "plot_id" := paste(tr$plot_id, as.character(d[[year]]), sep = ".")]
    tr[, "year" := d[[year]]]
  }
  tr[, "height" := if (!is.null(height)) as.numeric(d[[height]]) else NA_real_]
  tr[, "hcb"    := if (!is.null(hcb))    as.numeric(d[[hcb]])    else NA_real_]
  if (!is.null(crown_ratio) && is.null(hcb)) {
    tr[, "hcb" := tr$height * (1 - as.numeric(d[[crown_ratio]]))]
  }

  # Plot geometry -------------------------------------------------------
  scalar_or_col <- function(arg, default = NA_real_) {
    if (is.null(arg)) return(rep(default, nrow(d)))
    if (is.character(arg) && length(arg) == 1L && arg %in% names(d)) {
      return(as.numeric(d[[arg]]))
    }
    rep_len(as.numeric(arg), nrow(d))
  }
  tr[, "plot_radius" := scalar_or_col(plot_radius)]
  tr[, "plot_area"   := scalar_or_col(plot_area)]
  tr[, "slope"       := scalar_or_col(slope, 0)]
  tr[, "aspect"      := scalar_or_col(aspect, 0)]
  # circular plot area from radius where area was not given
  need_area <- is.na(tr$plot_area) & !is.na(tr$plot_radius)
  tr[need_area, "plot_area" := pi * tr$plot_radius[need_area]^2 / 10000]

  # Expansion factor ----------------------------------------------------
  if (!is.null(expf)) {
    tr[, "expf" := as.numeric(d[[expf]])]
  } else if (!all(is.na(tr$plot_area))) {
    tr[, "expf" := 1 / tr$plot_area]
  } else {
    stop("Supply `expf`, or `plot_area`, or `plot_radius`, so that per hectare ",
         "quantities can be scaled.", call. = FALSE)
  }

  # Coordinates ---------------------------------------------------------
  has_xy <- !is.null(x) && !is.null(y)
  has_polar <- !is.null(distance) && !is.null(azimuth)
  if (has_xy) {
    tr[, "x" := as.numeric(d[[x]])]
    tr[, "y" := as.numeric(d[[y]])]
    pol <- xy_to_polar(tr$x, tr$y)
    tr[, "distance" := pol$distance]
    tr[, "azimuth"  := pol$azimuth]
  } else if (has_polar) {
    tr[, "distance" := as.numeric(d[[distance]])]
    tr[, "azimuth"  := as.numeric(d[[azimuth]])]
    xy <- polar_to_xy(tr$distance, tr$azimuth)
    tr[, "x" := xy$x]
    tr[, "y" := xy$y]
  } else {
    tr[, "x" := NA_real_]; tr[, "y" := NA_real_]
    tr[, "distance" := NA_real_]; tr[, "azimuth" := NA_real_]
  }
  spatial <- has_xy || has_polar

  # Validity ------------------------------------------------------------
  if (any(tr$dbh <= 0, na.rm = TRUE)) {
    stop("`dbh` must be positive for every record.", call. = FALSE)
  }
  bad_cl <- !is.na(tr$height) & !is.na(tr$hcb) & (tr$height - tr$hcb) <= 0
  if (any(bad_cl)) {
    stop(sprintf("%d tree(s) have height <= hcb, so crown length is not positive.",
                 sum(bad_cl)), call. = FALSE)
  }
  dup <- duplicated(tr[, c("plot_id", "tree_id")])
  if (any(dup)) {
    stop("tree identifiers must be unique within a plot; ",
         sum(dup), " duplicate(s) found.", call. = FALSE)
  }
  if (spatial) {
    dupxy <- duplicated(tr[, c("plot_id", "x", "y")])
    if (any(dupxy)) {
      stop("trees within a plot must have distinct coordinates; ",
           sum(dupxy), " coincident stem location(s) found.", call. = FALSE)
    }
  }

  tr <- attach_traits(tr, traits, quiet = quiet)
  tr[, "ba_m2" := tree_ba(tr$dbh)]
  tr[, "ba_ha" := tr$ba_m2 * tr$expf]
  # crown ratio is computed first because the CONUS largest crown width form
  # reads it; the Acadian form ignores it, so this reorder changes no number
  # under the default source.
  tr[, "cr"  := (tr$height - tr$hcb) / tr$height]
  st <- if ("sp_type" %in% names(tr)) tr$sp_type else NULL
  tr[, "mcw" := max_crown_width(tr$dbh, tr$species, traits, sp_type = st)]
  tr[, "lcw" := largest_crown_width(tr$dbh, tr$species, traits, mcw = tr$mcw,
                                    cr = tr$cr, sp_type = st)]
  tr[, "mca" := crown_area_pct(tr$mcw, tr$expf)]

  plots <- tr[, list(
    n_trees     = .N,
    plot_area   = plot_area[1],
    plot_radius = plot_radius[1],
    slope       = slope[1],
    aspect      = aspect[1]
  ), by = "plot_id"]

  structure(
    list(trees = tr[], plots = plots[], traits = data.table::as.data.table(traits)),
    class = "stand",
    spatial = spatial,
    has_crown = all(!is.na(tr$height)) && all(!is.na(tr$hcb))
  )
}

#' @export
print.stand <- function(x, ...) {
  cat("<stand>\n")
  cat(sprintf("  %d plot(s), %d tree(s), %d species\n",
              nrow(x$plots), nrow(x$trees), length(unique(x$trees$species))))
  cat("  stem coordinates:", if (attr(x, "spatial")) "yes" else "no", "\n")
  cat("  crown data (height and hcb):",
      if (attr(x, "has_crown")) "complete" else "incomplete", "\n")
  cat(sprintf("  DBH range: %.1f to %.1f cm\n",
              min(x$trees$dbh, na.rm = TRUE), max(x$trees$dbh, na.rm = TRUE)))
  invisible(x)
}

#' Does this stand carry stem coordinates
#' @param x A `stand`.
#' @return Logical scalar.
#' @export
is_spatial <- function(x) isTRUE(attr(x, "spatial"))
