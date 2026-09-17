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
