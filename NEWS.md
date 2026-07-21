# anovapowersim (development version)

* Added the experimental, development-version-only `cell_design()`,
  `unbalanced_covariance()`, and `power_unbalanced()` functions for
  simulation-only power analysis of a fixed unbalanced allocation with
  user-defined cell means, marginal standard deviations, and within-subject
  correlations. Results include simulated power and partial eta-squared
  diagnostics, but deliberately omit calculated power.
* `power_unbalanced()` now derives a worst-case population
  Greenhouse--Geisser epsilon across between-subject cells for the tested
  term, reports it as `$epsilon`, and bases `power_sim` on the
  Greenhouse--Geisser-corrected simulated p-value whenever that epsilon is
  below 1 (requires `ss_type` `"III"` or `"II"`; a warning is issued if
  `ss_type = "I"` is combined with a non-spherical design).
* Added the experimental, development-version-only `power_achieved()` function
  for simulation-based achieved-power estimation at a fixed sample size and
  partial eta squared.
* Added the experimental, development-version-only `power_sensitivity()`
  function for simulation-based minimum-detectable partial eta-squared
  searches at a fixed sample size and target power.
* Added experimental, development-version-only `power_achieved_calc()` and
  `power_sensitivity_calc()` functions for equivalent fixed-sample analyses
  using calculated noncentral-F power without simulations.
* Added `power_n_calc()` for calculated-power, simulation-free sample-size searches in
  balanced ANOVA designs.
* Added an `epsilon` argument to `power_n_calc()` for calculated-power nonsphericity
  corrections on terms containing within-subject factors.
* Added `within_covariance()` and a `covariance` argument for `power_n()` and
  `power_curve()` so simulations can use custom within-subject covariance
  structures. These functions now derive a term-specific population
  Greenhouse--Geisser epsilon from that covariance and apply it to their
  calculated power.
* When a supplied covariance yields a population Greenhouse--Geisser epsilon
  below 1, `power_curve()`, `power_n()`, `power_achieved()`, and
  `power_sensitivity()` now base `power_sim` on each simulated dataset's
  Greenhouse--Geisser-corrected p-value instead of the uncorrected univariate
  test, so `power_sim` and `power_calc` estimate the same corrected test
  rather than diverging under non-sphericity. This correction requires
  `ss_type` `"III"` or `"II"`; under `"I"`, simulated p-values remain
  uncorrected, and these functions now warn when `ss_type = "I"` is combined
  with a covariance whose derived epsilon is below 1.

# anovapowersim 1.1.0

* Added a tolerance argument to `power_n()` for more precise control over the adaptive search.

# anovapowersim 1.0.0

* First official release
* Fixed a bug where adaptive search for purely between-subjects designs would fail if the starting N was too small

# anovapowersim 0.2.0

* Added parallel processing for simulation runs in `power_curve()` and
  `power_n()`. Use `parallel = TRUE` to enable
  parallel simulations and `cores` to control the number of cores.
