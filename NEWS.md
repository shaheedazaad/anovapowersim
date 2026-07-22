# anovapowersim (development version)

* Added `sim_correction = c("auto", "GG", "none")` to all simulation power
  functions. The default `"auto"` preserves existing behavior, while users can
  now prespecify corrected or uncorrected simulated tests. Uncorrected tests
  under nonsphericity warn that excess rejection reflects alpha inflation.
* Fixed `power_n()` and `power_n_calc()` treating `n_start` as an implicit
  lower bound when power at that value already met the target. Both searches
  now probe the smallest valid sample size and refine the resulting lower
  bracket before reporting `n_needed`.
* Balanced simulation functions now issue a once-per-session message when a
  custom `means_pattern` is resolved, clarifying that its values are projected,
  normalized, and rescaled to `target_pes`, unlike the literal means supplied
  through `cell_design()`. Both documentation pages now cross-reference this
  semantic distinction.
* `cell_design()` now messages the count and exact factor-level combinations
  of cells created by `default_n` and `default_m`, making accidental factor
  levels visible instead of silently expanding the design.
* `power_unbalanced()` now warns when the deterministic reference data imply
  essentially zero partial eta squared for the tested term, pointing users to
  possible mean typos or a mismatched `term`.
* Unbalanced within-subject designs now reject `:` in level values and reject
  duplicate cell names produced by joining multi-factor levels with `_`, with
  errors that identify the problematic levels or colliding cells before any
  correlations are assigned.
* Unbalanced power print and summary output now explain that simulated sample
  partial eta squared is upward-biased and that its mean, median, and interval
  are diagnostics rather than population/reference effects.
* `power_unbalanced()` now warns when `ss_type = "I"` is used with unequal
  sample sizes, explains that sequential sums of squares are order-dependent,
  and reports the factor order inherited from `cell_design()`.
* Balanced simulation power functions now require custom covariance inputs to
  be created by `within_covariance()` and reject raw matrices, eliminating
  ambiguous assumptions about within-cell row and column order.
* Added `means_pattern()` and an optional `means_pattern` argument to
  `power_curve()`, `power_n()`, `power_achieved()`, `power_sensitivity()`, and
  `design_term_means()`. Sparse relative cell values accept one-based indices
  or exact generated balanced-design level names, are broadcast over omitted
  factors, and are projected onto the requested ANOVA term before uniform
  calibration to `target_pes`.
* Balanced simulations now use a normalized centered-linear/Kronecker
  direction when no explicit pattern is supplied. Results record and print
  whether this documented default or a custom pattern was used.
* Balanced simulations now warn when an implicit default direction is
  consequential: the tested within-subject component has more than one degree
  of freedom and its population Greenhouse--Geisser epsilon is below one.
  Under nonsphericity, `power_sim` can depend on mean direction even when
  `target_pes` and covariance are fixed. Calculation-only functions retain
  the conventional direction-insensitive noncentral-F approximation and do
  not warn; `power_unbalanced()` already receives literal means.
* The `gpower = TRUE` warning is now issued whenever `gpower = TRUE` is used,
  not only for within-subject terms with more than one degree of freedom.
  G*Power's estimates can differ from `target_pes` more broadly than that;
  the default `gpower = FALSE` remains recommended.
* `power_n()` now rejects `n_start` values above `n_max` instead of running
  the first simulation outside the requested search range.
* Balanced simulation results now retain the full-precision simulated power
  used by adaptive searches and report valid/failed fit counts. Printed power
  values remain formatted to three decimals.
* **Breaking:** simulation APIs now require one common marginal variance.
  `within_covariance()` replaces `default_sd` with `sd` and removes
  measurement-specific `standard_deviations`; direct covariance matrices must
  have equal diagonal variances. For unbalanced designs, remove cell-level
  `sd` and `default_sd` from `cell_design()` and supply the common SD through
  `unbalanced_covariance(sd = ...)`. Unequal correlations and
  Greenhouse--Geisser correction remain supported.
* Simulation functions now warn when an omitted covariance causes the common
  `sd = 1` or within-subject correlation `0.5` defaults to be used. The
  covariance constructors warn when `sd` is omitted, and resolved covariance
  specifications warn when `default_correlation` fills unnamed pairs while
  preserving every explicitly supplied correlation.
* Added the experimental, development-version-only `cell_design()`,
  `unbalanced_covariance()`, and `power_unbalanced()` functions for
  simulation-only power analysis of a fixed unbalanced allocation with
  user-defined cell means and sample sizes under a common standard deviation
  and optional within-subject correlations. Results include simulated power
  and partial eta-squared diagnostics, but deliberately omit calculated power.
* `power_unbalanced()` derives the population Greenhouse--Geisser epsilon from
  the covariance matrix shared across between-subject cells, reports it as
  `$epsilon`, and bases `power_sim` on the
  Greenhouse--Geisser-corrected simulated p-value whenever that epsilon is
  below 1 (requires `ss_type` `"III"` or `"II"`; a warning is issued if
  `ss_type = "I"` is combined with a non-spherical design).
* **Breaking (experimental):** `cell_design()` now takes a `within` argument
  (character vector of within-subject factor names, or `NULL`) and stores it
  on the returned design; `power_unbalanced()` no longer accepts `within` and
  reads it from the design instead. Move `within = ...` from
  `power_unbalanced()` into `cell_design()`.
* `cell_design()` gained `default_n` and `default_m`. Supply both to auto-fill
  any missing cells in the complete factorial design; supplying only one is
  an error, and supplying neither requires every cell to be defined explicitly
  (as before).
* `cell_design()` now reports the exact missing factor-level combinations
  when a design is incomplete, instead of only a count, and errors clearly
  when a factor has fewer than two observed levels (previously this only
  surfaced later, inside `power_unbalanced()`, with an unhelpful low-level
  contrast-fitting error).
* The within-subject `n`-consistency check (that `n` is identical across all
  within-subject rows of the same between-subject cell) now runs in
  `cell_design()` at construction time; it previously only surfaced inside
  `power_unbalanced()`.
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
  `power_curve()` so simulations can use a custom common SD and within-subject
  correlation structure. These functions now derive a term-specific population
  Greenhouse--Geisser epsilon from that covariance and apply it to their
  calculated power.
* `power_curve()`, `power_n()`, `power_achieved()`, `power_sensitivity()`,
  `power_n_calc()`, `power_achieved_calc()`, `power_sensitivity_calc()`, and
  `design_term_means()` now warn when `gpower = TRUE` is combined with a term
  whose within-subject component has more than one degree of freedom (i.e. a
  within factor with more than two levels). In that case `target_pes` under
  `gpower = TRUE` does not equal the partial eta squared actually achieved --
  this mirrors a property of G*Power's own "as in Cohen (1988)"
  repeated-measures convention, which does not adjust for the number of
  measurements, rather than a bug in this package (`gpower = TRUE` remains an
  exact replica of G*Power's own noncentrality formula). Use the default
  `gpower = FALSE` when `target_pes` should match your reported or expected
  partial eta squared exactly.
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
