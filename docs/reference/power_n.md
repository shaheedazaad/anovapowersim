# Search for the sample size needed for target ANOVA power

Adaptive simulation search for the per-between-cell sample size needed
to reach a requested power for a balanced factorial ANOVA design. The
search searches upward from `n_start` until it brackets the target or
reaches `n_max`. If `n_start` already reaches the target, the search
probes the smallest sample size supported by the design to establish a
lower bracket. It then refines the bracket using interpolation with
midpoint bisection as a fallback.

## Usage

``` r
power_n(
  between = NULL,
  within = NULL,
  term,
  target_pes,
  power = 0.9,
  n_sims = 10000,
  alpha = 0.05,
  ss_type = "III",
  n_start = NULL,
  n_max = 5000,
  tol = 0.03,
  gpower = FALSE,
  progress = interactive(),
  parallel = FALSE,
  cores = NULL,
  seed = NULL,
  covariance = NULL,
  means_pattern = NULL,
  sim_correction = c("auto", "GG", "none")
)
```

## Arguments

- between:

  Named integer vector of between-subject factor level counts, e.g.
  `c(group = 2)`. Use `NULL` for no between-subject factors.

- within:

  Named integer vector of within-subject factor level counts, e.g.
  `c(time = 3, condition = 4)`. Use `NULL` for no within-subject
  factors.

- term:

  Character scalar naming the ANOVA term to test, e.g. `"group:time"`.
  Interaction terms are order-insensitive; `"time:group"` resolves to
  `"group:time"` when that is the design's factor order.

- target_pes:

  Target partial eta squared for `term`.

- power:

  Desired target power.

- n_sims:

  Number of simulated datasets per sample size.

- alpha:

  Significance threshold.

- ss_type:

  Sums-of-squares type for the tested ANOVA term. `"III"` is the default
  for order-invariant tests in unbalanced designs. Use `"I"` to
  reproduce sequential
  [`stats::aov()`](https://rdrr.io/r/stats/aov.html) tests.
  Greenhouse–Geisser-corrected simulated p-values are available only for
  `"III"` and `"II"`.

- n_start:

  Starting sample size per between-subject cell, not a lower bound for
  the search. If `NULL`, an initial value is estimated from calculated
  power and constrained to values that support empirical calibration for
  the requested design.

- n_max:

  Maximum sample size per between-subject cell.

- tol:

  Acceptable precision above target power. If no simulated value at or
  above `power` is also no more than `power + tol`, `power_n()` warns
  that the requested precision band was not reached.

- gpower:

  Logical; if `TRUE`, calibrate means to the G*Power-style noncentrality
  convention `lambda = total_n * f^2`. The default `FALSE` calibrates
  the empirical reference dataset to `target_pes`, equivalent to
  `lambda = den_df * f^2` for the fitted ANOVA. G*Power's estimates can
  differ from `target_pes`, especially for small samples or terms with
  more degrees of freedom; a warning is issued when `gpower = TRUE`. The
  default `gpower = FALSE` is recommended.

- progress:

  Logical; if `TRUE`, show a text progress bar.

- parallel:

  Logical; if `TRUE`, run simulations for each sample size via the
  `future` ecosystem.

- cores:

  Optional positive integer number of cores to use when
  `parallel = TRUE`. If `NULL`, uses one fewer than the number of
  available cores, with a minimum of one.

- seed:

  Optional integer seed for reproducibility.

- covariance:

  Optional within-subject covariance specification created by
  [`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md).
  Raw covariance matrices are not accepted, which avoids silently
  assuming a within-cell order. The default `NULL` uses standard
  deviations of `1` and a compound-symmetric correlation of `0.5` and
  issues a warning stating those defaults. A
  [`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
  specification issues a warning when correlation pairs are omitted: its
  `default_correlation` applies only to those undefined pairs, while
  explicitly defined correlations are unchanged. All measurements use
  the specification's common marginal variance, while unequal
  correlations remain supported. For terms containing within-subject
  factors, the resolved covariance matrix is also used to derive a
  term-specific population Greenhouse–Geisser epsilon for `power_calc`.
  With the default `sim_correction = "auto"`, a population epsilon below
  `1 - 1e-8` also selects Greenhouse–Geisser-corrected simulated
  p-values for `ss_type` `"II"` or `"III"`.

- means_pattern:

  Optional relative cell-mean shape created by
  [`means_pattern()`](https://shaheedazaad.github.io/anovapowersim/reference/means_pattern.md).
  The sparse values are projected onto `term`, normalized, and uniformly
  rescaled to reach `target_pes`. If `NULL`, simulations use the
  package's deterministic linear/Kronecker pattern. For multi-df
  nonspherical within-subject terms, simulated power is conditional on
  this direction, so an explicit pattern is recommended when the
  expected shape is known.

- sim_correction:

  Sphericity correction for simulated p-values: `"auto"` (the default)
  uses Greenhouse–Geisser correction when the term-specific population
  epsilon is below `1 - 1e-8` and `ss_type` is `"II"` or `"III"`; `"GG"`
  requests correction for every simulated dataset; and `"none"` always
  uses the uncorrected univariate test. `"GG"` is an error with
  `ss_type = "I"`. For a between-only term or a within component with
  one degree of freedom, `"GG"` silently resolves to `"none"` because no
  sphericity correction applies.

## Value

An `anovapowersim_curve` object with `n_needed` and `total_n_needed`.
For `power_n()`, `n_needed` is always an explicitly simulated
`n_per_cell` value, never an interpolated sample size. If the search
reaches target power but no simulated value lands inside
`[power, power + tol]`, `power_n()` reports the smallest explicitly
simulated value at or above target power and warns that the requested
precision band was not reached.

## Examples

    power_n(
      between = c(cond = 2),
      within = c(stim = 4),
      term = "cond:stim",
      target_pes = 0.14,
      alpha = 0.05,
      power = 0.90,
      n_sims = 1000, # use 5000+ for a more precise estimate
      seed = 123 # for reproducibility
    )
