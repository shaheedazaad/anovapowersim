# anovapowersim

`anovapowersim` is designed to make determining a priori power for ANOVAs as easy as possible. You can add as many within/between factors with as many levels as you would like. There's no need to estimate condition means, SDs, or repeated-measures correlations; just enter the target partial eta squared.

The package simulates data and estimates power based on the specified design. It also provides direct power calculations for comparison.

Getting a priori power for a 2 × 2 × 3 mixed interaction effect is as simple as running the following:

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

```text
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

## Additional power analyses

The development version includes several experimental power-analysis options.
Their full examples and guidance are kept in the dedicated guides linked below.

### Achieved power and sensitivity

At a fixed sample size, `power_achieved()` estimates power for a chosen partial
eta squared, while `power_sensitivity()` estimates the minimum detectable
partial eta squared. See the
[fixed-sample tutorial](https://shaheedazaad.github.io/anovapowersim/articles/fixed-sample-power.html).

### Calculation-only functions

The `_calc()` functions skip simulations and use calculated noncentral-F power.
They also support planned nonsphericity through `epsilon`. See the
[calculated-power tutorial](https://shaheedazaad.github.io/anovapowersim/articles/calculated-power.html).

### Power for unbalanced designs

`power_unbalanced()` simulates one exact allocation from user-defined cell
means and sample sizes under a common standard deviation and optional
within-subject correlations. It
is simulation-only and does not extrapolate how unequal cell sizes should
scale. See the
[unbalanced-design tutorial](https://shaheedazaad.github.io/anovapowersim/articles/unbalanced-designs.html).

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

Azaad, S. (2026). A priori power analysis for ANOVA interaction effects with the anovapowersim R package: a short introduction. https://doi.org/10.31234/osf.io/86rsy_v1.

## Limitations

`anovapowersim` is designed to be simple and easy to use first, which means it has some limitations for now. It does not support:

- Covariates (ANCOVAs)
- Sample-size searches or power curves for unbalanced designs
- Huynh-Feldt corrections in power simulations. Greenhouse-Geisser-corrected
  simulated tests are supported for sums-of-squares type II or III when a
  custom covariance implies `epsilon < 1`; type I tests remain uncorrected.
- Heteroskedastic ANOVA. Simulation functions require one common marginal
  variance; unequal correlations and Greenhouse--Geisser corrections remain
  supported for repeated-measures designs.

Simulation functions warn when their default common `sd = 1` or default
within-subject correlation of `0.5` is used. Covariance specifications also
warn when only some correlations are defined; the default correlation fills
only the undefined pairs.
- Simple main effects/pairwise comparisons

## Other packages

I recommend checking out [`Superpower`](https://aaroncaldwell.us/Superpower/), which handles some of the limitations above.
