# Simulate power for a fixed unbalanced ANOVA design

Estimates achieved power for exact, potentially unequal cell sizes and
user-supplied population means under one common standard deviation. This
function is simulation-only: it does not calculate power from a
noncentral F distribution and does not scale the supplied sample sizes.

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

  Optional common covariance specification created by
  [`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md).
  The default `NULL` uses `sd = 1` and a correlation of `0.5` between
  within-subject measurements and issues a warning stating those
  defaults. For purely between-subject designs, use this argument to
  change the common SD; correlation settings are not applicable. When an
  [`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md)
  specification omits some correlation pairs, a warning states that
  `default_correlation` is used only for those undefined pairs. If the
  term's population covariance is non-spherical, `power_sim` is based on
  each simulated dataset's Greenhouse–Geisser-corrected p-value rather
  than the uncorrected univariate test. This correction requires
  `ss_type` `"III"` or `"II"`; under `"I"`, simulated p-values remain
  uncorrected, and a warning is issued.

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
population Greenhouse–Geisser epsilon for the tested term. `$results`
also reports the common SD, simulated partial eta-squared distribution,
and failed fits.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`power_unbalanced()` is experimental and is available only in the
development version of `anovapowersim`. Its API and reporting format may
change.

## Examples

``` r
# \donttest{
design <- cell_design(
  group = "control", time = "pre",  n = 12, m = 10,
  group = "control", time = "post", n = 12, m = 11,
  group = "treatment", time = "pre",  n = 18, m = 10,
  group = "treatment", time = "post", n = 18, m = 13,
  within = "time"
)
power_unbalanced(
  design = design,
  term = "group:time",
  covariance = unbalanced_covariance(
    sd = 2,
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
#>   common SD:             2
#>   simulated power:       0.970
#>   reference pes:         0.3000
#>   mean simulated pes:    0.3278
#>   median simulated pes:  0.3222
#>   simulated pes 95%:     [0.1132, 0.5841]
#>   within correlations:   custom
#>   SS type:               III
#> 
#> Cell design:
#> # A tibble: 4 × 4
#>   group     time      n     m
#>   <chr>     <chr> <int> <dbl>
#> 1 control   pre      12    10
#> 2 control   post     12    11
#> 3 treatment pre      18    10
#> 4 treatment post     18    13
#> 
#> Power and effect-size diagnostics:
#>  total_n min_cell_n max_cell_n n_sims valid_sims failed_sims sd epsilon num_df
#>       30         12         18    100        100           0  2       1      1
#>  den_df partial_eta_squared mean_pes_sim median_pes_sim pes_sim_lower
#>      28                 0.3    0.3278059      0.3222229      0.113171
#>  pes_sim_upper power_sim
#>      0.5841011     0.970
# }
```
