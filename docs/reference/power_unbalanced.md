# Simulate power for a fixed unbalanced ANOVA design

Estimates achieved power for exact, potentially unequal cell sizes and
user-supplied population means and standard deviations. This function is
simulation-only: it does not calculate power from a noncentral F
distribution and does not scale the supplied sample sizes.

## Usage

``` r
power_unbalanced(
  design,
  within = NULL,
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

- within:

  Character vector naming factors in `design` that are measured within
  subjects. Use `NULL` for a purely between-subject design.

- term:

  Character scalar naming the ANOVA term to test.

- covariance:

  Optional correlation specification created by
  [`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md).
  It can only be supplied when `within` is used. The default `NULL` uses
  a correlation of `0.5` between within-subject measurements. Standard
  deviations always come from `design`.

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
term effect size in a deterministic reference dataset. `$results` also
reports the simulated partial eta-squared distribution and failed fits.

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
  group = "treatment", time = "post", n = 18, m = 13, sd = 2.8
)
power_unbalanced(
  design = design,
  within = "time",
  term = "group:time",
  covariance = unbalanced_covariance(
    correlations = c("pre:post" = 0.7)
  ),
  n_sims = 100,
  seed = 123
)
#> Error: All simulated ANOVA fits failed. This usually indicates an internal model-fitting error rather than zero power.
# }
```
