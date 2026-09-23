# forestCI (development version)

No change to any exported function or any number the package returns. Four
verification scripts are added to `inst/scripts`, and the README records what
they found.

* `verify_crown_exposure_apa.R`, `verify_crown_exposure_matched.R` and
  `verify_open_sky_view.R` compare the crown exposure family and the exact
  polygon area potentially available against the original CI.R implementation.
  These were the least verified parts of the package, because the original's
  `OVERLAP()` and `POVii()` are commented out in its own driver and had never
  been exercised. The polygon APA is identical for every tree whose polygon is
  closed by neighbours, the analytic crown projection area is exact and the
  analytic crown surface area agrees to 0.264 percent once the crown radius
  input is matched, and the exposure ratios correlate at 0.984 and above.

* `validate_crosswalk_external.R` checks the 38 species codes that carry an FIA
  species code against a full FIA `REF_SPECIES` export, on species code, genus,
  softwood or hardwood class, common name and binomial. All 38 pass. This
  replaces a check that used the crown width library's own columns, which could
  not have caught a congeneric same clade mismapping.

* `?conus_crown_width` is unchanged; the library's SPCD 376, *Betula
  neoalaskana* Sarg., is an Alaska and boreal species in a nominally
  conterminous library. It is not one of the 38 mapped codes and affects no
  prediction.

# forestCI 0.2.0

A second crown width coefficient source, and one real defect fixed. The Acadian
source remains the default, so every number the package returned in 0.1.2 it
still returns, and the regression against the original scripts is unchanged.

* **`species_traits(source = "conus")` swaps in the CONUS maximum and largest
  crown width library.** Thirty eight of the 43 species codes carry an FIA
  species code and take the CONUS coefficients; the genus level and unknown
  codes (`AS`, `HI`, `OH`, `OS`, `99`) keep the Acadian type defaults, because
  no single FIA code is defensible for them. Every row records which it carries
  in `cw_form` and `source_crown`. Nothing else in the package needs an
  argument: `as_stand(traits = species_traits(source = "conus"))` carries the
  source through to crown competition factor, the crown exposure family, the
  influence zone neighbourhood and LCW weighted APA.

* The CONUS form is a different model of a different quantity, not a refit of
  Russell and Weiskittel (2011). Maximum crown width is the tau 0.95 quantile
  of crown width on diameter at crown ratio 1, evaluated through a three step
  recipe whose order matters: the trait corrected curve, a hold at
  `mcw_dbh_max`, then a cap at the clade crown ceiling, 22.86 m hardwood and
  19.87 m softwood. Largest crown width is a bounded fraction of that envelope
  rather than a second allometry, so it lies below maximum crown width by
  construction. `?max_crown_width` documents the three caveats that come with
  it: the quantile estimand against the mean that crown competition factor was
  calibrated on, the domain hold going flat inside a normal inventory, and a
  largest crown width ratio whose crown ratio and diameter slopes are global
  constants rather than species specific.

* This does **not** settle the Russell and Weiskittel Eq. 3 crown ratio
  coefficients. That remains open and still needs a reader with the rendered
  paper. The documentation now says so in both places.

* **An unknown species code no longer always takes hardwood coefficients.**
  `lookup_traits()` selected the hardwood default row unconditionally, so an
  unrecognised softwood got hardwood maximum crown width, about 37 percent high
  at 20 cm DBH and about 1.9 times high on crown area, crown competition factor
  and CCFL. The 0.1.2 `sp_type` fix repaired `attach_traits()` but `as_stand()`
  discarded that resolution one line later. `max_crown_width()`,
  `largest_crown_width()`, `max_crown_area()` and `crown_dimensions()` now take
  an optional `sp_type`, `as_stand()` passes the stand's own, and with no
  `sp_type` the hardwood default is still used, so no existing call changes.

* `largest_crown_width()` gains `cr`, used by the CONUS form and ignored by the
  Acadian one. `as_stand()` now computes crown ratio before crown width and
  passes it through.

* `as_stand()` no longer inlines the crown area arithmetic that
  `max_crown_area()` defines; both call one internal function, so they cannot
  drift apart.

* New data object `conus_crown_width`, 466 FIA species codes, metric, carrying
  the fitted domain, the clade ceiling, provenance and confidence class for
  every code. The English unit library was converted analytically and checked
  against the source project's own independent metric export: agreement to
  5.0e-15 relative in `mcw_a1`, exact in `mcw_a2`, 1.4e-13 cm in
  `dbh_max_fit`. Adds about 15 KB to the installed package. Seventeen of
  the 466 rows carry a single crown measured tree, so `dbh_min_fit` and
  `dbh_max_fit` are equal and the fitted range is a point. All seventeen
  are donor carried or borrowed and none of them is an FIA code that the
  forestCI trait table maps to, so no shipped prediction is held at a
  point by them; `?conus_crown_width` records the count.

# forestCI 0.1.2

Stress and scaling release. An 84 probe edge case suite was run over degenerate
stand sizes, rejected inputs, graceful degradation, ties, neighbourhood rules,
multi plot stands, weighting options, determinism, scaling and 40 random stands.
It found one defect and one gap, both fixed. Two performance problems it
measured are also fixed, and the regression against the original scripts is
unchanged by any of it.

* **An empty neighbourhood now means zero competition everywhere.** A tree with
  no competitor in a stand where other trees have competitors already returned
  zero. A stand in which no tree has a competitor returned `NA` for every
  distance-dependent index. Those two cases now agree, and only `mean_dist`,
  which is genuinely undefined without a competitor, stays `NA`.
* **`as_stand()` gains `sp_type`.** A caller who knows an unrecognised species
  is a softwood had no way to say so, and every unknown code fell back to the
  hardwood defaults. A supplied `sp_type` now steers that fallback and survives
  it. Recognised codes still take their type from the trait table.
* **Crown exposure is about twenty times faster on a large stand.** A neighbour
  can only overtop a facet when the two crowns overlap horizontally, and a ray
  can only be blocked within its own horizontal travel. Restricting the
  neighbour set by those two exact bounds takes 1000 trees from 224 s to 10.5 s
  and turns the cost from quadratic into roughly linear. No result changes.
* **Area potentially available is about eight times faster.** Neighbours are now
  clipped in order of distance and the loop stops once the nearest remaining
  bisector cannot reach the polygon, which is an exact early exit. 1000 trees
  goes from 9.0 s to 1.1 s. No result changes.

Timing on a stand at 0.12 trees per m2, seconds:

| Trees | independent | dependent | apa | rapa | crown exposure |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.006 | 0.013 | 0.11 | 0.01 | 0.99 |
| 400 | 0.007 | 0.032 | 0.44 | 0.11 | 4.24 |
| 1000 | 0.007 | 0.094 | 1.23 | 0.47 | 11.04 |

99 tests passing, R CMD check clean, all 84 stress probes passing, and the
regression against the original scripts identical to 0.1.1 in every row.

# forestCI 0.1.1

Verification release. The package was run against the original Acadian
competition index scripts (CI.R, SPP.R, AcadianGY, Run.R) on the Penobscot
plot, tree by tree. Seven quantities now agree to machine precision, three to
rounding, and the two that differ do so by design. Closing the last two gaps
required two real changes to the package.

* `crown_dimensions()` gains `method`. The new default `"adaptive"` integrates
  each crown section with [stats::integrate()], which is what the original
  scripts do and which handles the vertical tangent the profile develops at the
  crown tip. The previous fixed composite Simpson rule is still available as
  `method = "simpson"` and is several times faster, but it was up to 3.7 percent
  low on the most sharply pointed crowns. With the adaptive rule crown surface
  area agrees with the original to 0.26 percent and crown volume to 6e-6
  percent.
* `rapa()` gains `align`. The original rasterises on grid nodes from `-extent`
  to `+extent` inclusive rather than on cell centres. `align = "node"`
  reproduces that grid exactly, and with it rasterised APA is bit-identical to
  the original. The default stays `"center"`, where the sampled points tile the
  domain and the areas sum to it.
* `rapa()` and `apa()` gain `boundary = "square"` with `extent`, so a square
  domain can be specified directly rather than inferred from the stems.
* Species table corrections and confirmations, all against external sources
  rather than the code being tested: paper birch (`PB`) is a hardwood and the
  `"SW"` entry in the Acadian source is an error; `BK` is black locust
  (*Robinia pseudoacacia* L.) and `SC` is Scots pine (*Pinus sylvestris* L.),
  both confirmed twice in the Forest Vegetation Simulator species crosswalk and
  the Northeast variant overview.
* `largest_crown_width()` documentation now states which published equation it
  implements. It is Eq. 2 of Russell and Weiskittel (2011), the form without a
  crown ratio term. Eq. 3 of the same paper adds crown ratio in the denominator
  and fits better for 13 of the 15 species; it is not shipped, because the
  coefficient signs in the Acadian implementation and in the published table as
  machine-extracted disagree and the discrepancy is unresolved.
* `inst/scripts/regression_vs_original.R` runs the whole comparison and is the
  script to re-run after any change to the crown or neighbourhood code.

# forestCI 0.1.0

First release.

* `as_stand()` standardises an inventory table once, and every other function
  reads from it. Stem coordinates are optional: without them the
  distance-independent family still computes and the rest is reported as
  skipped rather than failing.
* Distance-independent indices: `ci_distance_independent()` covers basal area,
  basal area of larger trees with softwood and hardwood partitions, stand
  density index in both summation and quadratic mean diameter forms, relative
  density against either the Acadian mixed species or the Woodall maximum, and
  crown competition factor with its one-sided counterpart.
* Distance-dependent indices: `ci_distance_dependent()` covers Hegyi in both
  diameter and basal area forms, Martin and Ek, Daniels, Lorimer, Rouvinen and
  Kuuluvainen, Spurr point density, the Canham neighbourhood competition index,
  local basal area and local basal area of larger trees, and mean neighbour
  distance.
* Competitor selection is separated from competitor filtering in
  `neighbors()`: fixed radius, k nearest, angle gauge and crown influence zone
  selection, crossed with no filter, larger diameter, taller, or the
  height-fraction rule Waskiewicz (2011) found best.
* Growing space: `apa()` builds exact weighted Voronoi polygons by half-plane
  clipping, with the Moore et al. (1973) bisector displacement, and `rapa()`
  builds the rasterised multiplicatively weighted form with an optional
  directional modification.
* Crown geometry: `crown_profile()` supplies four families and
  `crown_exposure()` derives crown surface area, exposed crown surface area,
  exposed crown projection area and open sky view by sampling whichever
  profile was given.
* `spatial_pattern()` returns the Clark and Evans aggregation index with
  Donnelly edge correction, and the mean directional index.
* `edge_flag()` implements both the buffer rule and the exterior sector rule.
* Species traits: 43 codes shipped, softwood and hardwood fallback for anything
  else, and user tables merge in through `species_traits()`.
* Example data: two stem-mapped Penobscot Experimental Forest plots, `pef`.
