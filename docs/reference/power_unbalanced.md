# Simulate power for a fixed unbalanced ANOVA design

Estimates achieved power for exact, potentially unequal cell sizes and
user-supplied population means under one common standard deviation. This
function is simulation-only: it does not calculate power from a
noncentral F distribution and does not scale the supplied sample sizes.
A warning is issued when the deterministic reference data imply
essentially zero effect for the tested term, which often indicates a
typo, a wrong `term`, or means that contain only other effects.

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
  seed = NULL,
  sim_correction = c("auto", "GG", "none")
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
  `default_correlation` is used only for those undefined pairs. The
  resolved covariance determines the term-specific population
  Greenhouse–Geisser epsilon used by `sim_correction = "auto"`.

- n_sims:

  Number of simulated datasets.

- alpha:

  Significance threshold.

- ss_type:

  Sums-of-squares type: `"III"`, `"II"`, or `"I"`. For unequal-N
  designs, `"I"` uses sequential, order-dependent hypotheses; a warning
  reports the factor order inherited from
  [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md).

- progress:

  Logical; if `TRUE`, show a text progress bar.

- parallel:

  Logical; if `TRUE`, run simulations in parallel.

- cores:

  Optional positive integer number of parallel workers.

- seed:

  Optional integer seed for reproducibility.

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

An `anovapowersim_unbalanced_power` object. `$power` and
`$achieved_power` contain simulated power. `$partial_eta_squared` is the
term effect size in a deterministic reference dataset. `$epsilon` is the
population Greenhouse–Geisser epsilon for the tested term. `$results`
also reports the common SD, simulated partial eta-squared distribution,
and failed fits. Sample partial eta squared is upward-biased in finite
samples, so `mean_pes_sim`, `median_pes_sim`, and the simulated interval
are sampling diagnostics rather than estimates of the supplied
population effect; use `$partial_eta_squared` as the reference effect.
The object stores the requested `sim_correction` and applied
`sim_correction_resolved` values.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`power_unbalanced()` is experimental and is available only in the
development version of `anovapowersim`. Its API and reporting format may
change.

## Simulated sphericity correction

`sim_correction` governs only the simulated test. With `"GG"`, each
dataset uses its own sample-estimated Greenhouse–Geisser epsilon from
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html). Under a truly
spherical population, that sample correction is mildly conservative.
Power is estimated for the prespecified corrected or uncorrected test;
conditional Mauchly-then-correct procedures are not simulated. Unlike
the balanced simulation functions, `power_unbalanced()` does not report
a `power_calc` diagnostic.

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
#>   simulated test:        uncorrected (auto)
#>   common SD:             2
#>   simulated power:       0.970
#>   reference pes:         0.3000
#>   mean simulated pes:    0.3278
#>   median simulated pes:  0.3222
#>   simulated pes 95%:     [0.1132, 0.5841]
#>   pes note:              sample pes is upward-biased; simulated pes summaries are diagnostics, not the population/reference effect.
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
