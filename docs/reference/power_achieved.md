# Estimate achieved ANOVA power at a fixed sample size

Simulates one balanced factorial ANOVA design point for a fixed partial
eta squared and sample size. The simulation estimate is the primary
achieved power result; noncentral-F calculated power is retained as a
diagnostic.

## Usage

``` r
power_achieved(
  between = NULL,
  within = NULL,
  term,
  target_pes,
  n,
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

- n:

  Number of subjects per between-subject cell. For a purely
  within-subject design, this is the total sample size.

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

  Optional within-subject covariance specification from
  [`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
  or a numeric covariance matrix. The default `NULL` uses standard
  deviations of `1` and a compound-symmetric correlation of `0.5` and
  issues a warning stating those defaults. A
  [`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
  specification issues a warning when correlation pairs are omitted: its
  `default_correlation` applies only to those undefined pairs, while
  explicitly defined correlations are unchanged. All measurements must
  have one common marginal variance. A supplied matrix must therefore
  have equal diagonal entries and one row and column per within-subject
  cell; named matrices are reordered to the design's cell order. Unequal
  correlations remain supported. For terms containing within-subject
  factors, the matrix is also used to derive a term-specific population
  Greenhouse–Geisser epsilon for `power_calc`. If that population
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

An `anovapowersim_achieved_power` object. `$results` contains the
standard one-row power diagnostics. `$achieved_power` is the simulated
power estimate and `$calculated_power` is the calculated-power
diagnostic.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`power_achieved()` is experimental and is available only in the
development version of `anovapowersim`. Its API and reporting format may
change.

## Examples

``` r
# \donttest{
power_achieved(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  target_pes = 0.14,
  n = 20,
  n_sims = 100,
  seed = 123
)
#> Warning: No `covariance` was supplied; using one common `sd = 1` and `correlation = 0.5` for every within-subject pair.
#> <anovapowersim_achieved_power>
#>   term:             'group:time'
#>   fixed n per cell: 20
#>   fixed total N:    40
#>   target pes:       0.1400
#>   alpha:            0.05
#>   simulations:      100
#>   means pattern:    default linear/Kronecker
#>   achieved power:   0.640
#>   calculated power: 0.679
#>   SS type:          III
#> 
#>  n_per_cell total_n n_sims valid_sims failed_sims epsilon num_df den_df   ncp
#>          20      40    100        100           0       1      1     38 6.186
#>  power_calc power_sim
#>       0.679     0.640
# }
```
