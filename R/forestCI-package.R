#' forestCI: competition indices for forest stands
#'
#' `forestCI` computes individual tree competition indices from stand
#' inventory data. Distance-independent indices require only species,
#' diameter and an expansion factor. Distance-dependent indices additionally
#' require stem coordinates, supplied either directly as x and y or as
#' distance and azimuth from plot centre. Crown based indices require a
#' crown profile, which the package supplies through a pluggable interface
#' with softwood and hardwood defaults so that species with no local
#' parameterisation are still handled.
#'
#' All internal computation is metric. Diameter at breast height is in cm,
#' heights and distances in m, basal area in m2 ha-1, crown competition
#' factor in percent, and stand density index in trees ha-1.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom data.table := .N .SD .I
#' @importFrom data.table setDT setorder setorderv setnames setkeyv
#' @importFrom data.table as.data.table data.table rbindlist copy
#' @importFrom grDevices chull
#' @importFrom stats integrate setNames weighted.mean approx optim
#' @importFrom utils head modifyList data
NULL

# data.table non standard evaluation guards
utils::globalVariables(c(
  ".", ".I", ".N", ".SD", "..keep",
  "apa", "apa_share", "aspect",
  "ba", "ba_ha", "ba_ha_j", "ba_m2", "ba_m2_j", "ba_sw", "bal", "bal_sw",
  "ccfl", "ccfl_sw", "cr", "csa", "csax",
  "d_ij", "dbh", "dbh_j", "expf", "hcb", "height", "is_edge",
  "lcw", "lo_exp", "mca", "mca_sw", "n_cells", "n_comp",
  "plot_area", "plot_id", "plot_radius", "qmd",
  "rapa", "rapa_share", "sdi", "sdi_tree", "sg", "shade", "shape",
  "slope", "sp_type", "species", "species_j", "tph", "tree_id",
  "up_exp", "widest", "x", "y"
))

# data.table uses this to confirm the package intends its `:=` semantics
.datatable.aware <- TRUE
