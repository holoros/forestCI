# forestCI

Individual tree competition indices for forest stands, distance-dependent and
distance-independent, with a pluggable crown profile that generalises across
softwoods and hardwoods.

The package is a generalisation of the competition index code developed for the
Acadian Forest Ecoregion by Waskiewicz (2011) and evaluated by Kuehne,
Weiskittel and Waskiewicz (2019). The original implementation was a set of
scripts tied to one species pool and one plot design. This package keeps the
methods, drops the global variables, and makes three things general: any species
pool through a trait table with softwood and hardwood fallback, any plot design
through an explicit plot geometry, and any crown shape through a crown profile
interface.

## Installation

```r
# install.packages("remotes")
remotes::install_github("holoros/forestCI")
```

## What it computes

Distance-independent, available for any stand with species, diameter and an
expansion factor:

| Index | Sided | Reference |
|:---|:---|:---|
| Basal area, `ba` | two | conventional |
| Basal area of larger trees, `bal`, `bal_sw`, `bal_hw` | one | conventional |
| Stand density index, `sdi`, `sdi_reineke` | two | Reineke (1933) |
| Relative density, `rd` | two | Woodall et al. (2005), Acadian mixed species |
| Crown competition factor, `ccf` | two | Krajicek et al. (1961) |
| Crown competition factor of larger trees, `ccfl`, `ccfl_sw`, `ccfl_hw` | one | conventional |

Distance-dependent, requiring stem coordinates:

| Index | Sided | Reference |
|:---|:---|:---|
| Hegyi, `hegyi`, `hegyi_ba` | two | Hegyi (1974) |
| Martin and Ek, `martin_ek` | two | Martin and Ek (1984) |
| Daniels, `daniels` | two | Daniels (1976) |
| Lorimer, `lorimer` | two | Lorimer (1983) |
| Rouvinen and Kuuluvainen, `rouvinen` | two | Rouvinen and Kuuluvainen (1997) |
| Spurr point density, `spurr` | two | Spurr (1962) |
| Neighbourhood competition index, `nci` | two | Canham et al. (2004) |
| Local basal area, `local_ba`, `local_bal` | two, one | Kuehne et al. (2019) |
| Mean neighbour distance, `mean_dist` | two | Motz et al. (2010) |
| Area potentially available, `apa` | one | Moore et al. (1973) |
| Rasterised APA, `rapa` | one | Mu (2004) |
| Exposed crown surface area, `csax`, `csax_rel` | one | Cole and Lorimer (1994) |
| Exposed crown projection area, `cpax`, `cpax_rel` | one | Cole and Lorimer (1994) |
| Open sky view, `osv`, `osv_rel` | one | Van Pelt and North (1996) |
| Clark and Evans aggregation, `clark_evans` | plot | Clark and Evans (1954) |
| Mean directional index, `mdi` | plot | Corral-Rivas et al. (2010) |

## Quick start

```r
library(forestCI)

st <- pef_stand()          # two stem-mapped Penobscot plots
st
#> <stand>
#>   2 plot(s), 74 tree(s), 7 species
#>   stem coordinates: yes
#>   crown data (height and hcb): complete

ci <- competition_indices(st, radius = 6, edge = "exclude")
head(ci[, .(tree_id, species, dbh, bal, ccfl, hegyi, apa, csax_rel)])
```

Without stem coordinates nothing breaks. The distance-independent family is
computed, the rest is reported as skipped, and the returned table has the same
shape.

```r
st_flat <- as_stand(my_inventory, plot = "plot", tree = "tree",
                    species = "spp", dbh = "dbh", plot_area = 0.04)
competition_indices(st_flat)
```

## Verification against the original implementation

The package is a rewrite, so it is checked against the scripts it replaces
rather than only against itself. `inst/scripts/regression_vs_original.R` runs
the original CI.R, SPP.R and Acadian pipeline on the Penobscot plot and compares
every tree. On R 4.5.2:

| Quantity | Agreement |
|:---|:---|
| Maximum crown width, maximum crown area, crown competition factor | identical |
| Hegyi index | identical (1.7e-14 percent) |
| Spurr point density | identical (2.8e-14 percent) |
| Angle gauge competitor count | identical |
| Rasterised APA, `align = "node"` | identical |
| Crown volume, measured radii | 6e-6 percent |
| Plot basal area | 2.3e-4 percent |
| Crown surface area, measured radii | 0.26 percent |
| Basal area and crown competition factor of larger trees | one tied diameter pair |
| Crown surface area and volume from predicted rather than measured crown width | different input, not different method |

The crown exposure family and the exact polygon area potentially available were
added to that record on 2026-09-23, using `inst/scripts/verify_crown_exposure_apa.R`,
`inst/scripts/verify_crown_exposure_matched.R` and
`inst/scripts/verify_open_sky_view.R`. The original's `OVERLAP()` and `POVii()`
are commented out in its own driver and had never been run, so this is the first
time these quantities have been compared to anything outside the package.

| Quantity | Agreement |
|:---|:---|
| Polygon APA, every tree whose polygon is closed by neighbours | identical (15 of 25 trees, below 1e-10) |
| Polygon APA, trees with an unbounded polygon | not comparable; the original returns 0 or NA |
| Analytic crown projection area, measured radii | identical (3.8e-13 percent) |
| Analytic crown surface area, measured radii | 0.264 percent |
| Exposed crown projection ratio | r = 0.996, stand means within 0.7 percent |
| Exposed crown surface ratio | r = 0.984, stand means within 3.1 percent |
| Open sky view, sky fraction | r = 0.834, forestCI 13.5 percent high |

The last three rows are characterizations rather than certifications, and the
reason is in the reference. Halving the original's raster cell side from 0.25 m
to 0.125 m moves its own exposed surface ratio by up to 27.7 percent on an
individual tree, so it is not a precise reference at the resolution it is
normally run at. Open sky view differs by construction: `forestCI` casts one ray
per crown facet along the outward normal, while the original integrates the sky
fraction over the whole upward hemisphere at each facet, and a single normal ray
escapes more readily, so the high bias is the expected direction.

The 38 species codes that carry an FIA species code are checked against a full
FIA `REF_SPECIES` export by `inst/scripts/validate_crosswalk_external.R`. All 38
pass on species code, genus, softwood or hardwood class, common name and
binomial. Supply a full export, roughly 2,700 rows; a subset whose species code
set is the crown width library's own set confirms internal consistency only.

The last two rows are deliberate. `forestCI` treats trees of equal diameter
symmetrically, so a tied pair each see only the trees strictly larger than both,
where a plain cumulative sum gives the second one of the pair the first one's
basal area as well. And the original takes crown radius from the four cardinal
field measurements while the package predicts largest crown width from diameter;
feed the package the measured radii and the geometry agrees to a quarter of a
percent.

Two changes came out of this comparison rather than out of the test suite: crown
dimensions are integrated adaptively by default, because a fixed Simpson rule
was up to 3.7 percent low where the crown profile turns vertical at the tip, and
`rapa()` can align its grid on nodes rather than cell centres.

## Scale

Timing on one plot at 0.12 trees per m², R 4.5.2, seconds:

| Trees | `ci_distance_independent` | `ci_distance_dependent` | `apa` | `rapa` | `crown_exposure` |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.006 | 0.013 | 0.11 | 0.01 | 0.99 |
| 400 | 0.007 | 0.032 | 0.44 | 0.11 | 4.24 |
| 1000 | 0.007 | 0.094 | 1.23 | 0.47 | 11.04 |

Crown exposure and area potentially available are the two expensive routines,
and both are near linear in stand size rather than quadratic, because a
neighbour that cannot reach the subject is excluded before any geometry is
evaluated. Both exclusions are exact bounds, so nothing is traded for the speed.
Open sky view roughly doubles the crown exposure cost; pass `crown_osv = FALSE`
to skip it.

## Crown profiles

Every crown based index reads its geometry from one object, so changing the
crown model changes crown surface area, exposed crown surface area, open sky
view and crown volume together.

```r
crown_profile("dual_exponent")      # species constants, Waskiewicz (2011)
crown_profile("variable_exponent")  # crown ratio and softwood/hardwood only
crown_profile("geometric", solid = "paraboloid")   # transparent baseline
crown_profile("custom", fun = function(z) (1 - z)^0.4)
```

`crown_profile_info()` reports where each constant came from. The
`variable_exponent` coefficients shipped here are a calibration to the package's
own dual exponent defaults, not a field fit, and the function says so.

## Moving outside the Acadian species pool

The built-in trait table covers 43 codes of the northeastern United States and
eastern Canada. Unknown codes fall back to the softwood or hardwood default with
a message, never silently. To add species, pass a table:

```r
tr <- species_traits(data.frame(
  code = c("DF", "WH"), sp_type = "SW",
  sg = c(0.45, 0.42), shade = c(2.8, 4.8),
  widest = 0.75, up_exp = 3, lo_exp = 3, shape = "g"))

st <- as_stand(my_inventory, ..., traits = tr)
```

Only `code` and `sp_type` are required; anything omitted takes the softwood or
hardwood default.

## Units

Metric throughout. Diameter cm, heights and distances m, basal area m² ha⁻¹,
stand density index trees ha⁻¹, crown competition factor percent, areas m².

## References

Kuehne, C., Weiskittel, A.R., Waskiewicz, J. (2019) Comparing performance of
contrasting distance-independent and distance-dependent competition metrics in
predicting individual tree diameter increment and survival within
structurally-heterogeneous, mixed-species forests of Northeastern United States.
*Forest Ecology and Management* 433: 205-216.
doi:10.1016/j.foreco.2018.11.002

Waskiewicz, J.D. (2011) *Influence of neighborhood structure on growth in
northern red oak and eastern white pine stands.* PhD dissertation, University of
Maine.

Weiskittel, A.R., Hann, D.W., Kershaw, J.A., Vanclay, J.K. (2011) *Forest Growth
and Yield Modeling.* Wiley.

## Licence

MIT.
