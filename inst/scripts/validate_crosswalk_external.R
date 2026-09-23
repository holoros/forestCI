## validate_crosswalk_external_20260923.R
## External validation of the forestCI species crosswalk against FIA REF_SPECIES.
##
## Usage:
##   Rscript validate_crosswalk_external_20260923.R <REF_SPECIES.csv> [outdir] [libpath]
##
## Every file read is wrapped in tryCatch; failures are appended to
## error_log.txt in the output directory rather than stopping the run.
## No graphics. data.table only.

suppressMessages(library(data.table))

args    <- commandArgs(trailingOnly = TRUE)
# The reference must be a full FIA REF_SPECIES export, roughly 2,700 rows. A
# subset whose species code set is the crown width library's own set is not an
# external reference and will only confirm internal consistency.
REFPATH <- if (length(args) >= 1) args[1] else "REF_SPECIES.csv"
OUTDIR  <- if (length(args) >= 2) args[2] else "."
LIBPATH <- if (length(args) >= 3) args[3] else NA_character_

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
ERRLOG <- file.path(OUTDIR, "error_log.txt")

log_err <- function(what, e) {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time()), what,
              conditionMessage(e)), file = ERRLOG, append = TRUE)
  NULL
}
safely <- function(what, expr) tryCatch(expr, error = function(e) log_err(what, e))

## ---------------------------------------------------------------- inputs ---
ref <- safely(paste("read REF_SPECIES:", REFPATH),
              fread(REFPATH, colClasses = "character", na.strings = c("", "NA")))

pkg <- safely("load forestCI", {
  if (!is.na(LIBPATH)) .libPaths(c(LIBPATH, .libPaths()))
  suppressMessages(library(forestCI))
  TRUE
})

ccw <- safely("get conus_crown_width", as.data.table(get("conus_crown_width")))
st  <- safely("get species_traits_default", as.data.table(get("species_traits_default")))

if (is.null(ref) || is.null(ccw) || is.null(st)) {
  cat("FAIL: one or more inputs could not be read; see ", ERRLOG, "\n", sep = "")
  quit(save = "no", status = 1)
}

ref[, SPCD := as.integer(SPCD)]

## --------------------------------------- FVS Northeast variant convention ---
## What each 2-letter code is supposed to name. This table is supplied
## independently of the package so gate G5 is not self-referential.
fvs_ne <- data.table(
  code = c("AB","AE","BA","BC","BF","BK","BP","BS","BT","BW","EC","EH","GA",
           "GB","HH","JP","NC","NS","PB","PC","PR","QA","RB","RM","RN","RO",
           "RP","RS","SB","SC","SM","ST","TA","WA","WC","WP","WS","YB"),
  fvs_common = c("American beech","American elm","black ash","black cherry",
    "balsam fir","black locust","balsam poplar","black spruce","bigtooth aspen",
    "American basswood","eastern cottonwood","eastern hemlock","green ash",
    "gray birch","eastern hophornbeam","jack pine","northern white cedar",
    "Norway spruce","paper birch","pin cherry","pin cherry","quaking aspen",
    "river birch","red maple","red pine","northern red oak","red pine",
    "red spruce","sweet birch","Scots pine","sugar maple","striped maple",
    "tamarack","white ash","northern white cedar","eastern white pine",
    "white spruce","yellow birch"),
  fvs_genus = c("Fagus","Ulmus","Fraxinus","Prunus","Abies","Robinia","Populus",
    "Picea","Populus","Tilia","Populus","Tsuga","Fraxinus","Betula","Ostrya",
    "Pinus","Thuja","Picea","Betula","Prunus","Prunus","Populus","Betula",
    "Acer","Pinus","Quercus","Pinus","Picea","Betula","Pinus","Acer","Acer",
    "Larix","Fraxinus","Thuja","Pinus","Picea","Betula"),
  fvs_species = c("grandifolia","americana","nigra","serotina","balsamea",
    "pseudoacacia","balsamifera","mariana","grandidentata","americana",
    "deltoides","canadensis","pennsylvanica","populifolia","virginiana",
    "banksiana","occidentalis","abies","papyrifera","pensylvanica",
    "pensylvanica","tremuloides","nigra","rubrum","resinosa","rubra",
    "resinosa","rubens","lenta","sylvestris","saccharum","pensylvanicum",
    "laricina","americana","occidentalis","strobus","glauca","alleghaniensis"),
  fvs_type = c("HW","HW","HW","HW","SW","HW","HW","SW","HW","HW","HW","SW",
    "HW","HW","HW","SW","SW","SW","HW","HW","HW","HW","HW","HW","SW","HW",
    "SW","SW","HW","SW","HW","HW","SW","HW","SW","SW","SW","HW")
)

## ------------------------------------------------------------- normalisers --
norm <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("\\(native\\)|\\(introduced\\)", "", x)
  x <- gsub("[-_/]", " ", x)
  x <- gsub("[^a-z ]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}
## accepted common-name synonyms (FIA vintage vs FVS usage)
syn <- function(x) {
  x <- norm(x)
  x <- gsub("^scotch pine$", "scots pine", x)
  x <- gsub("^northern white cedar$", "northern whitecedar", x)
  x <- gsub("^northern whitecedar$", "northern whitecedar", x)
  x <- gsub("^hophornbeam$", "eastern hophornbeam", x)
  x <- gsub("^eastern hophornbeam$", "eastern hophornbeam", x)
  x
}
tokset <- function(x) lapply(strsplit(syn(x), " "), function(z) sort(unique(z[nzchar(z)])))
name_ok <- function(a, b) {
  A <- tokset(a); B <- tokset(b)
  mapply(function(u, v) {
    if (!length(u) || !length(v)) return(NA)
    identical(u, v) || length(intersect(u, v)) >= min(length(u), length(v))
  }, A, B)
}

## ------------------------------------------------------------------ build ---
X <- merge(st[!is.na(spcd), .(code, spcd, sp_type,
                              st_sci = scientific_name, st_common = common_name)],
           fvs_ne, by = "code", all.x = TRUE)
X <- merge(X, ref[, .(spcd = SPCD, ref_genus = GENUS, ref_species = SPECIES,
                      ref_common = COMMON_NAME, ref_sh = SFTWD_HRDWD)],
           by = "spcd", all.x = TRUE)
X <- merge(X, ccw[, .(spcd, ccw_genus = genus, ccw_clade = clade,
                      ccw_common = common_name)],
           by = "spcd", all.x = TRUE)
setorder(X, code)

## ------------------------------------------------------------------ gates ---
## G1 SPCD exists in REF_SPECIES
X[, G1_in_ref := !is.na(ref_genus)]
## G2 REF GENUS == conus_crown_width genus, case-insensitive exact
X[, G2_genus := G1_in_ref & !is.na(ccw_genus) &
      tolower(trimws(ref_genus)) == tolower(trimws(ccw_genus))]
## G3 SFTWD_HRDWD agrees with sp_type and with ccw clade
X[, ref_type := fifelse(ref_sh == "S", "SW", fifelse(ref_sh == "H", "HW", NA_character_))]
X[, G3_type := G1_in_ref & !is.na(ref_type) & ref_type == sp_type &
      !is.na(ccw_clade) & substr(ref_type, 1, 1) == ccw_clade]
## G4 REF COMMON_NAME consistent with what the 2-letter code means
X[, G4_common := G1_in_ref & name_ok(ref_common, fvs_common) %in% TRUE]
X[, G4_common := G1_in_ref & mapply(function(a, b) isTRUE(name_ok(a, b)),
                                    ref_common, fvs_common)]
## G5 binomial GENUS + SPECIES is the species the FVS NE code names
X[, ref_binom := trimws(paste(ref_genus, ref_species))]
X[, fvs_binom := trimws(paste(fvs_genus, fvs_species))]
X[, G5_binom := G1_in_ref & norm(ref_binom) == norm(fvs_binom)]
## cross-check: package scientific_name vs FVS NE expectation
X[, pkg_binom_ok := norm(st_sci) == norm(fvs_binom)]

X[, all_pass := G1_in_ref & G2_genus & G3_type & G4_common & G5_binom]

out <- X[, .(code, spcd, sp_type,
             ref_genus, ref_species, ref_common, ref_sh,
             ccw_genus, ccw_clade, ccw_common,
             fvs_common, fvs_binom, st_sci, ref_binom,
             G1_in_ref, G2_genus, G3_type, G4_common, G5_binom,
             pkg_binom_ok, all_pass)]

csv <- file.path(OUTDIR, "crosswalk_external_results_20260923.csv")
safely(paste("write", csv), fwrite(out, csv))

## ---------------------------------------------------------------- summary ---
n  <- nrow(out); np <- sum(out$all_pass, na.rm = TRUE)
cat(sprintf(
  "%s: %d/%d forestCI codes pass all 5 external gates against %s (G1 %d, G2 %d, G3 %d, G4 %d, G5 %d)\n",
  if (np == n) "PASS" else "FAIL", np, n, basename(REFPATH),
  sum(out$G1_in_ref, na.rm = TRUE), sum(out$G2_genus, na.rm = TRUE),
  sum(out$G3_type, na.rm = TRUE), sum(out$G4_common, na.rm = TRUE),
  sum(out$G5_binom, na.rm = TRUE)))
if (np < n) {
  cat("\nFailing rows:\n")
  print(out[all_pass %in% c(FALSE, NA),
            .(code, spcd, ref_binom, ref_common, ccw_genus, ccw_clade,
              G1_in_ref, G2_genus, G3_type, G4_common, G5_binom)])
}
