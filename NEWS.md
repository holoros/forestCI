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
