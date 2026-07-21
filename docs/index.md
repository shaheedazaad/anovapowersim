# anovapowersim

`anovapowersim` is designed to make determining a priori power for
ANOVAs as easy as possible. You can add as many within/between factors
with as many levels as you would like. There’s no need to estimate
condition means, SDs, or repeated-measures correlations; just enter the
target partial eta squared.

The package simulates data and estimates power based on the specified
design. It also provides direct power calculations for comparison.

Getting a priori power for a 2 × 2 × 3 mixed interaction effect is as
simple as running the following:

``` r

install.packages("anovapowersim") # if not already installed

library(anovapowersim)

power_n(
  between = c(group = 2), # group has 2 levels
  within = c(stim = 2, cond = 3), # stim has 2 levels, cond has 3
  term = "group:stim:cond", # three-way interaction term
  target_pes = 0.08, # target effect size
  n_sims = 5000, # increase to 10000+ for more precise estimates
  power = .90,
  alpha = .05,
  parallel = TRUE, # simulations will be run in parallel for speed
  seed = 123 # for reproducibility
)
```

In the development version, the experimental
[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md)
function estimates power for a chosen effect size at a fixed sample
size, while the experimental
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
function finds the minimum detectable partial eta squared. These
functions are not yet available in the CRAN release.

``` r

power_achieved(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  target_pes = 0.08,
  n = 30,
  seed = 123
)

power_sensitivity(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  n = 30,
  power = 0.90,
  seed = 123
)
```

To skip simulations, use the experimental calculation-only counterparts.
They use calculated noncentral-F power and support `gpower` and planned
nonsphericity through `epsilon`:

``` r

power_achieved_calc(
  between = c(group = 2),
  within = c(time = 3),
  term = "group:time",
  target_pes = 0.08,
  n = 30,
  gpower = TRUE,
  epsilon = 0.80
)

power_sensitivity_calc(
  between = c(group = 2),
  within = c(time = 3),
  term = "group:time",
  n = 30,
  power = 0.90,
  gpower = TRUE,
  epsilon = 0.80
)
```

## Power for unbalanced designs

The experimental, development-version-only
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
function simulates one exact unbalanced allocation from population cell
means and standard deviations. For repeated measures,
[`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md)
supplies correlations only; each cell’s standard deviation remains
defined in
[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md).

``` r

design <- cell_design(
  group = "control",   time = "pre",  n = 22, m = 10.0, sd = 2.0,
  group = "control",   time = "post", n = 22, m = 11.0, sd = 2.2,
  group = "treatment", time = "pre",  n = 31, m = 10.1, sd = 2.4,
  group = "treatment", time = "post", n = 31, m = 12.4, sd = 2.8
)

power_unbalanced(
  design = design,
  within = "time",
  term = "group:time",
  covariance = unbalanced_covariance(
    correlations = c("pre:post" = 0.7)
  ),
  n_sims = 5000,
  parallel = TRUE,
  seed = 123
)
```

This function is simulation-only. It reports simulated power and
reference and simulation-based partial eta-squared summaries for the
exact design; it does not return calculated power or extrapolate how the
cell sizes should scale.

``` text
#><anovapowersim_curve>
#>  term:          'group:stim:cond'
#>  target power:  0.900
#>  alpha:         0.05
#>  effect size:   pes = 0.08
#>  n values:      8 per-cell sample sizes visited
#>  sims per cell size: 5000
#>  SS type:       III
#>  n needed for between-subjects cell: 38
#>  total N needed: 76
#>
#> n_per_cell total_n n_sims num_df den_df    ncp power_calc power_sim
#>         31      62   5000      2    120 10.435      0.823     0.825
#>         37      74   5000      2    144 12.522      0.890     0.885
#>         38      76   5000      2    148 12.870      0.899     0.903
#>         39      78   5000      2    152 13.217      0.907     0.901
#>         40      80   5000      2    156 13.565      0.915     0.918
#>         41      82   5000      2    160 13.913      0.922     0.916
#>         46      92   5000      2    180 15.652      0.949     0.947
#>         62     124   5000      2    244 21.217      0.989     0.988
```

## Installation

`anovapowersim` can be installed from CRAN:

``` r

install.packages("anovapowersim")
```

You can install the development version from GitHub:

``` r

install.packages("pak")
pak::pak("shaheedazaad/anovapowersim")
```

Or, with `remotes`:

``` r

install.packages("remotes")
remotes::install_github("shaheedazaad/anovapowersim")
```

## Citation

Azaad, S. (2026). A priori power analysis for ANOVA interaction effects
with the anovapowersim R package: a short introduction.
<https://doi.org/10.31234/osf.io/86rsy_v1>.

## Limitations

`anovapowersim` is designed to be simple and easy to use first, which
means it has some limitations for now. It does not support:

- Covariates (ANCOVAs)
- Sample-size searches or power curves for unbalanced designs
- Sample-estimated Greenhouse-Geisser or Huynh-Feldt corrections in
  power simulations. Custom nonspherical covariance structures and their
  derived population Greenhouse-Geisser corrections are supported by
  [`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md),
  [`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md),
  [`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md),
  and
  [`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md),
  and planned epsilon corrections are supported by the calculation-only
  functions.
- Specific interaction shapes (based on means)
- Simple main effects/pairwise comparisons

## Other packages

I recommend checking out
[`Superpower`](https://aaroncaldwell.us/Superpower/), which handles some
of the limitations above.
