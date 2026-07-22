# Estimate ANOVA effect-size sensitivity at a fixed sample size

Searches for the minimum detectable partial eta squared at a fixed
sample size and target power. Calculated power supplies an efficient
initial estimate. Explicit simulations then bracket the target and
refine the effect-size bracket using interpolation with midpoint
fallback.

## Usage

``` r
power_sensitivity(
  between = NULL,
  within = NULL,
  term,
  n,
  power = 0.9,
  n_sims = 10000,
  alpha = 0.05,
  ss_type = "III",
  pes_min = 1e-06,
  pes_max = 0.99,
  pes_tol = 0.001,
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

- n:

  Number of subjects per between-subject cell. For a purely
  within-subject design, this is the total sample size.

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
  Greenhouse–Geisser-corrected simulated p-values (see `covariance`) are
  only available for `"III"` and `"II"`; under `"I"`, simulated p-values
  always use the uncorrected univariate test, and a warning is issued if
  the supplied covariance yields a population epsilon below `1`.

- pes_min:

  Lower bound of the partial eta-squared search interval.

- pes_max:

  Upper bound of the partial eta-squared search interval.

- pes_tol:

  Maximum width of the final simulated partial eta-squared bracket.

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

An `anovapowersim_sensitivity` object. `$pes_needed` is the explicitly
simulated upper effect-size bracket, or `NA` when `pes_max` does not
achieve target power. `$results` contains every explicitly simulated
effect size and its standard power diagnostics.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`power_sensitivity()` is experimental and is available only in the
development version of `anovapowersim`. Its API and reporting format may
change.

## Examples

``` r
# \donttest{
power_sensitivity(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  n = 20,
  power = 0.90,
  n_sims = 100,
  pes_tol = 0.01,
  seed = 123
)
#> Warning: No `covariance` was supplied; using one common `sd = 1` and `correlation = 0.5` for every within-subject pair.
#> Warning: `power_sim` and `power_calc` differ by more than 5 percentage points for target_pes = 0.204402088273331 (largest difference = 0.061). Try increasing `n_sims` for a more stable simulation estimate.
#> <anovapowersim_sensitivity>
#>   term:             'group:time'
#>   fixed n per cell: 20
#>   fixed total N:    40
#>   target power:     0.900
#>   alpha:            0.05
#>   detectable pes:   0.211450
#>   search interval:  [1e-06, 0.99]
#>   requested width:  0.01
#>   simulated points: 7
#>   final width:      0.007048313
#>   converged:        yes
#>   simulations/point: 100 
#>   means pattern:    default linear/Kronecker
#>   SS type:          III
#> 
#>  target_pes n_per_cell total_n n_sims valid_sims failed_sims epsilon num_df
#>    0.000001         20      40    100        100           0       1      1
#>    0.112774         20      40    100        100           0       1      1
#>    0.169161         20      40    100        100           0       1      1
#>    0.197354         20      40    100        100           0       1      1
#>    0.204402         20      40    100        100           0       1      1
#>    0.211450         20      40    100        100           0       1      1
#>    0.225547         20      40    100        100           0       1      1
#>  den_df    ncp power_calc power_sim
#>      38  0.000      0.050     0.050
#>      38  4.830      0.572     0.580
#>      38  7.737      0.774     0.780
#>      38  9.343      0.846     0.880
#>      38  9.763      0.861     0.800
#>      38 10.190      0.875     0.920
#>      38 11.067      0.900     0.900
# }
```
