#' Species trait defaults used by forestCI
#'
#' Returns the built-in species trait table, optionally overridden or extended
#' by a user supplied table. Every crown based index in the package reads its
#' species constants from here, so supplying a table for a new region is the
#' single step needed to move the package outside the Acadian species pool.
#'
#' The built-in table covers 43 species codes of the northeastern United States
#' and eastern Canada, using the two letter Forest Vegetation Simulator Acadian
#' variant codes. Three entries were checked against the Forest Vegetation
#' Simulator species crosswalk and Northeast variant overview rather than
#' against the source code: `BK` is black locust (*Robinia pseudoacacia* L.),
#' `SC` is Scots pine (*Pinus sylvestris* L.), and `PB` is paper birch
#' (*Betula papyrifera* Marshall), which is a hardwood. The Acadian source
#' classifies `PB` as a softwood, which is an error, and this table corrects
#' it. Specific gravity and shade tolerance follow the Acadian
#' variant parameter set, shade tolerance itself on the scale of Niinemets and
#' Valladares (2006). Crown shape constants (`widest`, `up_exp`, `lo_exp`,
#' `shape`) follow the parameterisation of Waskiewicz (2011). Maximum and
#' largest crown width coefficients follow Russell and Weiskittel (2011) as
#' implemented in the Acadian variant.
#'
#' @section Columns:
#' \describe{
#'   \item{code}{Species code, matched against the `species` column of a stand.}
#'   \item{scientific_name, authority, common_name}{Nomenclature. Authority is
#'     held separately so that the binomial can be italicised without it.}
#'   \item{sp_type}{`"SW"` softwood or `"HW"` hardwood. This is the axis the
#'     package falls back on for any species it does not recognise.}
#'   \item{sg}{Green specific gravity, dimensionless.}
#'   \item{shade}{Shade tolerance, 0 (intolerant) to 5 (tolerant).}
#'   \item{widest}{Position of the widest point of the crown as a proportion of
#'     crown length measured from the apex. 0 is a crown widest at the apex, 1 a
#'     crown widest at the base.}
#'   \item{up_exp, lo_exp}{Profile exponents for the crown surface above and
#'     below the widest point.}
#'   \item{shape}{`"p"` paraboloid, `"g"` gothic dome, `"e"` ellipsoid.}
#'   \item{mcw_a1, mcw_a2}{Maximum crown width coefficients,
#'     MCW = a1 * DBH^a2, with MCW in m and DBH in cm.}
#'   \item{lcw_b1, lcw_b2}{Largest crown width coefficients,
#'     LCW = MCW / (b1 * DBH^b2).}
#'   \item{source_crown}{`"species"` where the crown constants are specific to
#'     that species, `"type_default"` where they are the softwood or hardwood
#'     default.}
#' }
#'
#' @param extra Optional `data.frame` of the same shape as the built-in table.
#'   Rows whose `code` matches a built-in row replace it; rows with a new
#'   `code` are appended. Columns absent from `extra` are taken from the
#'   softwood or hardwood default, so a minimal override needs only `code`,
#'   `sp_type` and the columns actually being changed.
#'
#' @return A `data.table` of species traits, keyed on `code`.
#' @references
#' Niinemets, U., Valladares, F. (2006) Tolerance to shade, drought, and
#' waterlogging of temperate northern hemisphere trees and shrubs.
#' Ecological Monographs 76: 521-547.
#'
#' Russell, M.B., Weiskittel, A.R. (2011) Maximum and largest crown width
#' equations for 15 tree species in Maine. Northern Journal of Applied
#' Forestry 28: 84-91.
#'
#' Waskiewicz, J.D. (2011) Influence of neighborhood structure on growth in
#' northern red oak and eastern white pine stands. PhD dissertation,
#' University of Maine.
#' @param source Coefficient source for maximum and largest crown width,
#'   either `"acadian"`, the default, which is Russell and Weiskittel (2011)
#'   for 15 species in Maine, or `"conus"`, the CONUS library keyed on FIA
#'   species code. Under `"conus"` the 38 species that carry an FIA code take
#'   the CONUS coefficients and the domain hold and clade ceiling that go with
#'   them; the genus level and unknown codes (AS, HI, OH, OS, 99) keep the
#'   Acadian type defaults, because no single FIA code is defensible for them.
#'   Every row records which it carries in `cw_form` and `source_crown`. See
#'   [max_crown_width()] for the two forms and for the caveats on the CONUS
#'   one. `extra` is applied after the source swap, so a user supplied
#'   coefficient still wins.
#' @export
#' @examples
#' head(species_traits())
#' # the CONUS library instead of the Acadian one
#' table(species_traits(source = "conus")$cw_form)
#' # add a western species using hardwood defaults for crown shape
#' species_traits(data.frame(code = "DF", sp_type = "SW", sg = 0.45,
#'                           shade = 2.8, common_name = "Douglas-fir"))
species_traits <- function(extra = NULL,
                           source = c("acadian", "conus")) {
  source <- match.arg(source)
  tr <- data.table::copy(forestCI::species_traits_default)
  if (identical(source, "conus")) tr <- apply_conus_crown_width(tr)
  if (is.null(extra)) {
    data.table::setkeyv(tr, "code")
    return(tr[])
  }
  extra <- data.table::as.data.table(extra)
  need_cols(extra, c("code", "sp_type"), "`extra`")
  if (!all(extra$sp_type %in% c("SW", "HW"))) {
    stop("`sp_type` in `extra` must be \"SW\" or \"HW\".", call. = FALSE)
  }
  defaults <- type_defaults(tr)
  filled <- defaults[match(extra$sp_type, defaults$sp_type), ]
  filled[, "code" := extra$code]
  for (nm in setdiff(names(extra), "code")) {
    filled[[nm]] <- extra[[nm]]
  }
  tr <- tr[!tr$code %in% filled$code, ]
  out <- data.table::rbindlist(list(tr, filled), use.names = TRUE, fill = TRUE)
  data.table::setkeyv(out, "code")
  out[]
}

# Softwood and hardwood fallback rows, derived from the built-in table
type_defaults <- function(tr) {
  d <- tr[tr$source_crown == "type_default" & tr$code %in% c("OS", "OH"), ]
  d <- data.table::copy(d)
  d[, "code" := c("HW", "SW")[match(d$sp_type, c("HW", "SW"))]]
  d[, "source_crown" := "type_default"]
  data.table::setorder(d, "sp_type")
  d[]
}

#' Attach species traits to a tree table
#'
#' Joins the trait table onto trees by species code. Codes with no match fall
#' back to the softwood or hardwood default, chosen from `sp_type` when the
#' caller supplies one and from the hardwood default otherwise, and the
#' substitution is reported once rather than silently. A caller supplied
#' `sp_type` survives the fallback: it is the caller's statement about the
#' tree, so it is not overwritten by the default row it selected.
#'
#' @param trees A `data.frame` or `data.table` with a `species` column and,
#'   optionally, an `sp_type` column.
#' @param traits Trait table, normally from [species_traits()].
#' @param quiet Suppress the message listing unmatched codes.
#' @return A `data.table` with the trait columns added.
#' @export
attach_traits <- function(trees, traits = species_traits(), quiet = FALSE) {
  trees <- data.table::as.data.table(trees)
  need_cols(trees, "species", "`trees`")
  traits <- data.table::as.data.table(traits)
  trait_cols <- setdiff(names(traits), "code")

  idx <- match(as.character(trees$species), traits$code)
  unmatched <- unique(as.character(trees$species)[is.na(idx)])

  if (length(unmatched)) {
    fallback_type <- if ("sp_type" %in% names(trees)) {
      as.character(trees$sp_type)
    } else {
      rep("HW", nrow(trees))
    }
    fallback_type[!fallback_type %in% c("SW", "HW")] <- "HW"
    defaults <- type_defaults(traits)
    fb_idx <- match(fallback_type, defaults$sp_type)
    if (!quiet) {
      message(sprintf(
        "forestCI: %d species code(s) not in the trait table, using %s defaults: %s",
        length(unmatched), "softwood or hardwood",
        paste(unmatched, collapse = ", ")))
    }
    for (cl in trait_cols) {
      v <- traits[[cl]][idx]
      v[is.na(idx)] <- defaults[[cl]][fb_idx][is.na(idx)]
      trees[, (cl) := v]
    }
    # a caller supplied sp_type is the caller's statement about the tree, so it
    # survives the fallback rather than being overwritten by the default row
    if ("sp_type" %in% trait_cols) {
      keep <- is.na(idx)
      trees[keep, "sp_type" := fallback_type[keep]]
    }
  } else {
    for (cl in trait_cols) trees[, (cl) := traits[[cl]][idx]]
  }
  trees[]
}

# Replace the Acadian crown width coefficients with the CONUS library for
# every species code that carries an FIA species code, leaving the rest on the
# Acadian type defaults. The five codes that keep them (AS, HI, OH, OS and 99)
# are genus level or explicitly unknown, so no single FIA code is defensible
# for them; `source_crown` and `cw_form` record which rows moved.
apply_conus_crown_width <- function(tr) {
  L <- data.table::as.data.table(forestCI::conus_crown_width)
  i <- match(tr$spcd, L$spcd)
  hit <- !is.na(i)
  if (!any(hit)) return(tr)
  j <- i[hit]
  tr[hit, "mcw_a1"        := L$mcw_a1[j]]
  tr[hit, "mcw_a2"        := L$mcw_a2[j]]
  tr[hit, "lcw_b1"        := NA_real_]
  tr[hit, "lcw_b2"        := NA_real_]
  tr[hit, "cw_form"       := "conus"]
  tr[hit, "mcw_dbh_max"   := L$dbh_max_fit[j]]
  tr[hit, "mcw_dbh_min"   := L$dbh_min_fit[j]]
  tr[hit, "mcw_ceiling"   := L$cw_ceiling[j]]
  tr[hit, "lcw_r0"        := L$lcw_r0[j]]
  tr[hit, "lcw_rcr"       := L$lcw_rcr[j]]
  tr[hit, "lcw_rdbh"      := L$lcw_rdbh[j]]
  tr[hit, "cw_provenance" := L$provenance[j]]
  tr[hit, "cw_confidence" := L$confidence[j]]
  tr[hit, "cw_n_obs"      := L$n_obs[j]]
  tr[hit, "source_crown"  := "conus"]
  tr[]
}
