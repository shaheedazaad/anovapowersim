# Simulate ANOVA power from a balanced factorial design

Simulation-based power estimation for balanced factorial designs. Users
specify the between- and within-subject factors, the ANOVA term to test,
a target partial eta squared, and explicit sample sizes. The function
projects an explicit relative means pattern (or uses the documented
linear/Kronecker default), scales it to the requested partial eta
squared, simulates datasets, refits the ANOVA, and estimates power by
counting `p < alpha`.

## Usage

``` r
power_curve(
  between = NULL,
  within = NULL,
  term,
  target_pes,
  n_range,
  n_sims = 10000,
  alpha = 0.05,
  ss_type = "III",
  gpower = FALSE,
  progress = interactive(),
  parallel = FALSE,
  cores = NULL,
  seed = NULL,
  covariance = NULL,
  means_pattern = NULL
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

- n_range:

  Integer vector of sample sizes per between-subject cell. For pure
  within-subject designs, this is the total sample size.

- n_sims:

  Number of simulated datasets per sample size.

- alpha:

  Significance threshold.

- ss_type:

  Sums-of-squares type for the tested ANOVA term. `"III"` is the default
  for order-invariant tests in unbalanced designs. Use `"I"` to
  reproduce sequential
  [`stats::aov()`](https://rdrr.io/r/stats/aov.html) tests.
  Greenhouse–Geisser-corrected simulated p-values (see `covariance`) are
  only available for `"III"` and `"II"`; under `"I"`, simulated p-values
  always use the uncorrected univariate test, and a warning is issued if
  the supplied covariance yields a population epsilon below `1`.

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
  [`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
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
  If that population
  epsilon is below `1`, `power_sim` is also based on each simulated
  dataset's Greenhouse–Geisser-corrected p-value (from
  [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html)) rather than
  the uncorrected univariate test, so `power_sim` and `power_calc`
  estimate the same corrected test.

- means_pattern:

  Optional relative cell-mean shape created by
  [`means_pattern()`](https://shaheedazaad.github.io/anovapowersim/reference/means_pattern.md).
  The sparse values are projected onto `term`, normalized, and uniformly
  rescaled to reach `target_pes`. If `NULL`, simulations use the
  package's deterministic linear/Kronecker pattern. For multi-df
  nonspherical within-subject terms, simulated power is conditional on
  this direction, so an explicit pattern is recommended when the
  expected shape is known.

## Value

An `anovapowersim_curve` object. The `$results` tibble contains
`n_per_cell`, `total_n`, `n_sims`, successful and failed simulation
counts (`valid_sims`, `failed_sims`), the population nonsphericity
correction (`epsilon`), numerator and denominator degrees of freedom
(`num_df`, `den_df`), the noncentrality parameter (`ncp`), calculated
power (`power_calc`), and simulated power (`power_sim`). The
full-precision `power_sim` value, not its printed three-decimal
representation, is used by adaptive searches. When `epsilon < 1`, the
reported degrees of freedom and noncentrality are the corrected values
used for `power_calc`, and `power_sim` (for `ss_type` `"III"`/`"II"`) is
based on the Greenhouse–Geisser-corrected simulated p-value rather than
the uncorrected univariate test. Balanced simulation result objects also
include `custom_means_pattern`, indicating whether the relative
direction was supplied explicitly.

## Examples

    power_curve(
      between = c(cond = 2),
      within = c(stim = 2),
      term = "cond:stim",
      target_pes = 0.14,
      n_range = c(16, 20, 23, 28), # n per between-subject cell
      n_sims = 1000,
      seed = 123
    )

    power_curve(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.14,
      n_range = c(12, 16, 20),
      n_sims = 5000,
      parallel = TRUE,
      cores = 4,
      seed = 123
    )

    power_curve(
      within = c(time = 4),
      term = "time",
      target_pes = 0.15,
      n_range = 30,
      means_pattern = means_pattern(
        time = 1, value = 0,
        time = 2, value = 0.3,
        time = 3, value = 0.5,
        time = 4, value = 0.6
      ),
      n_sims = 1000,
      seed = 123
    )
