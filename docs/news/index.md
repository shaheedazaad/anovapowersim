# Changelog

## anovapowersim (development version)

- Added the experimental, development-version-only
  [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md),
  [`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md),
  and
  [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
  functions for simulation-only power analysis of a fixed unbalanced
  allocation with user-defined cell means, marginal standard deviations,
  and within-subject correlations. Results include simulated power and
  partial eta-squared diagnostics, but deliberately omit calculated
  power.
- [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
  now derives a worst-case population Greenhouse–Geisser epsilon across
  between-subject cells for the tested term, reports it as `$epsilon`,
  and bases `power_sim` on the Greenhouse–Geisser-corrected simulated
  p-value whenever that epsilon is below 1 (requires `ss_type` `"III"`
  or `"II"`; a warning is issued if `ss_type = "I"` is combined with a
  non-spherical design).
- **Breaking (experimental):**
  [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)
  now takes a `within` argument (character vector of within-subject
  factor names, or `NULL`) and stores it on the returned design;
  [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
  no longer accepts `within` and reads it from the design instead. Move
  `within = ...` from
  [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
  into
  [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md).
- [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)
  gained `default_n`, `default_m`, and `default_sd`. Supply all three to
  auto-fill any missing cells in the complete factorial design;
  supplying only some of the three is an error, and supplying none
  requires every cell to be defined explicitly (as before).
- [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)
  now reports the exact missing factor-level combinations when a design
  is incomplete, instead of only a count, and errors clearly when a
  factor has fewer than two observed levels (previously this only
  surfaced later, inside
  [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md),
  with an unhelpful low-level contrast-fitting error).
- The within-subject `n`-consistency check (that `n` is identical across
  all within-subject rows of the same between-subject cell) now runs in
  [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)
  at construction time; it previously only surfaced inside
  [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md).
- Added the experimental, development-version-only
  [`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md)
  function for simulation-based achieved-power estimation at a fixed
  sample size and partial eta squared.
- Added the experimental, development-version-only
  [`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
  function for simulation-based minimum-detectable partial eta-squared
  searches at a fixed sample size and target power.
- Added experimental, development-version-only
  [`power_achieved_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved_calc.md)
  and
  [`power_sensitivity_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity_calc.md)
  functions for equivalent fixed-sample analyses using calculated
  noncentral-F power without simulations.
- Added
  [`power_n_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n_calc.md)
  for calculated-power, simulation-free sample-size searches in balanced
  ANOVA designs.
- Added an `epsilon` argument to
  [`power_n_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n_calc.md)
  for calculated-power nonsphericity corrections on terms containing
  within-subject factors.
- Added
  [`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
  and a `covariance` argument for
  [`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
  and
  [`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md)
  so simulations can use custom within-subject covariance structures.
  These functions now derive a term-specific population
  Greenhouse–Geisser epsilon from that covariance and apply it to their
  calculated power.
- [`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md),
  [`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md),
  [`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md),
  [`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md),
  [`power_n_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n_calc.md),
  [`power_achieved_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved_calc.md),
  [`power_sensitivity_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity_calc.md),
  and
  [`design_term_means()`](https://shaheedazaad.github.io/anovapowersim/reference/design_term_means.md)
  now warn when `gpower = TRUE` is combined with a term whose
  within-subject component has more than one degree of freedom (i.e. a
  within factor with more than two levels). In that case `target_pes`
  under `gpower = TRUE` does not equal the partial eta squared actually
  achieved – this mirrors a property of G*Power’s own “as in Cohen
  (1988)” repeated-measures convention, which does not adjust for the
  number of measurements, rather than a bug in this package
  (`gpower = TRUE` remains an exact replica of G*Power’s own
  noncentrality formula). Use the default `gpower = FALSE` when
  `target_pes` should match your reported or expected partial eta
  squared exactly.
- When a supplied covariance yields a population Greenhouse–Geisser
  epsilon below 1,
  [`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md),
  [`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md),
  [`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md),
  and
  [`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
  now base `power_sim` on each simulated dataset’s
  Greenhouse–Geisser-corrected p-value instead of the uncorrected
  univariate test, so `power_sim` and `power_calc` estimate the same
  corrected test rather than diverging under non-sphericity. This
  correction requires `ss_type` `"III"` or `"II"`; under `"I"`,
  simulated p-values remain uncorrected, and these functions now warn
  when `ss_type = "I"` is combined with a covariance whose derived
  epsilon is below 1.

## anovapowersim 1.1.0

CRAN release: 2026-05-31

- Added a tolerance argument to
  [`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
  for more precise control over the adaptive search.

## anovapowersim 1.0.0

CRAN release: 2026-05-28

- First official release
- Fixed a bug where adaptive search for purely between-subjects designs
  would fail if the starting N was too small

## anovapowersim 0.2.0

- Added parallel processing for simulation runs in
  [`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md)
  and
  [`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md).
  Use `parallel = TRUE` to enable parallel simulations and `cores` to
  control the number of cores.
