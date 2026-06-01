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

## Citation

Azaad, S. (2026). A priori power analysis for ANOVA interaction effects
with the anovapowersim R package: a short introduction.
<https://doi.org/10.31234/osf.io/86rsy_v1>.

## Limitations

`anovapowersim` is designed to be simple and easy to use first, which
means it has some limitations for now. It does not support:

- Covariates (ANCOVAs)
- Means-based simulations for unbalanced designs
- Nonsphericity corrections (though this might change)
- Specific interaction shapes (based on means)
- Simple main effects/pairwise comparisons

## Other packages

I recommend checking out
[`Superpower`](https://aaroncaldwell.us/Superpower/), which handles some
of the limitations above.
