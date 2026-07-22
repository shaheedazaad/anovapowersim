# Simulate power for a fixed unbalanced ANOVA design

Estimates achieved power for exact, potentially unequal cell sizes and
user-supplied population means and standard deviations. This function is
simulation-only: it does not calculate power from a noncentral F
distribution and does not scale the supplied sample sizes.

## Usage

``` r
power_unbalanced(
  design,
  term,
  covariance = NULL,
  n_sims = 10000,
  alpha = 0.05,
  ss_type = "III",
  progress = interactive(),
  parallel = FALSE,
  cores = NULL,
  seed = NULL
)
```

## Arguments

- design:

  A complete design table created by
  [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md).
  Within-subject factors are read from the design's `within` attribute
  (set via
  [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)'s
  `within` argument), not supplied here.

- term:

  Character scalar naming the ANOVA term to test.

- covariance:

  Optional correlation specification created by
  [`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md).
  It can only be supplied when `design` has within-subject factors. The
  default `NULL` uses a correlation of `0.5` between within-subject
  measurements. Standard deviations always come from `design` and may
  differ between between-subject cells. If the term's population
  covariance is non-spherical for any cell, `power_sim` is based on each
  simulated dataset's Greenhouse–Geisser-corrected p-value rather than
  the uncorrected univariate test. This correction requires `ss_type`
  `"III"` or `"II"`; under `"I"`, simulated p-values remain uncorrected,
  and a warning is issued.

- n_sims:

  Number of simulated datasets.

- alpha:

  Significance threshold.

- ss_type:

  Sums-of-squares type: `"III"`, `"II"`, or `"I"`.

- progress:

  Logical; if `TRUE`, show a text progress bar.

- parallel:

  Logical; if `TRUE`, run simulations in parallel.

- cores:

  Optional positive integer number of parallel workers.

- seed:

  Optional integer seed for reproducibility.

## Value

An `anovapowersim_unbalanced_power` object. `$power` and
`$achieved_power` contain simulated power. `$partial_eta_squared` is the
term effect size in a deterministic reference dataset. `$epsilon` is the
worst-case (smallest) population Greenhouse–Geisser epsilon across
between-subject cells for the tested term. `$results` also reports the
simulated partial eta-squared distribution and failed fits.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`power_unbalanced()` is experimental and is available only in the
development version of `anovapowersim`. Its API and reporting format may
change.

## Examples

``` r
# \donttest{
design <- cell_design(
  group = "control", time = "pre",  n = 12, m = 10, sd = 2.0,
  group = "control", time = "post", n = 12, m = 11, sd = 2.2,
  group = "treatment", time = "pre",  n = 18, m = 10, sd = 2.4,
  group = "treatment", time = "post", n = 18, m = 13, sd = 2.8,
  within = "time"
)
power_unbalanced(
  design = design,
  term = "group:time",
  covariance = unbalanced_covariance(
    correlations = c("pre:post" = 0.7)
  ),
  n_sims = 100,
  seed = 123
)
#> <anovapowersim_unbalanced_power>
#>   term:                  'group:time'
#>   fixed total N:         30
#>   between-cell n range:  12 to 18
#>   simulations:           100
#>   simulated power:       0.820
#>   reference pes:         0.2223
#>   mean simulated pes:    0.2556
#>   median simulated pes:  0.2497
#>   simulated pes 95%:     [0.0469, 0.5434]
#>   within correlations:   custom
#>   SS type:               III
#> 
#> Cell design:
#> # A tibble: 4 × 5
#>   group     time      n     m    sd
#>   <chr>     <chr> <int> <dbl> <dbl>
#> 1 control   pre      12    10   2  
#> 2 control   post     12    11   2.2
#> 3 treatment pre      18    10   2.4
#> 4 treatment post     18    13   2.8
#> 
#> Power and effect-size diagnostics:
#>  total_n min_cell_n max_cell_n n_sims valid_sims failed_sims epsilon num_df
#>       30         12         18    100        100           0       1      1
#>  den_df partial_eta_squared mean_pes_sim median_pes_sim pes_sim_lower
#>      28           0.2223183    0.2555821      0.2497243    0.04690494
#>  pes_sim_upper power_sim
#>      0.5434147     0.820
# }
```
