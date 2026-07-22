# Getting started with anovapowersim

`anovapowersim` simulates power for balanced factorial ANOVA designs.
Specify the factors and levels, the term of interest, and a target
partial eta squared. The package generates default term-specific cell
means, simulates datasets, refits the ANOVA, and estimates power.

``` r

library(anovapowersim)
```

## Search for the required sample size

Use
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
to search for the sample size needed to reach the requested power. This
example is a mixed design with one two-level between-subject factor
(`cond`) and one four-level within-subject factor (`stim`). It tests the
`cond:stim` interaction with 90% power to detect a partial eta squared
of 0.14.

``` r

power_n(
  between = c(cond = 2), # cond has 2 levels
  within = c(stim = 4), # stim has 4 levels
  term = "cond:stim",
  target_pes = 0.14,
  alpha = 0.05,
  power = 0.90,
  n_sims = 1000, # use 5000+ for a more precise estimate
  seed = 123 # for reproducibility
)
```

    #> <anovapowersim_curve>
    #>   term:          'cond:stim'
    #>   target power:  0.900
    #>   alpha:         0.05
    #>   effect size:   pes = 0.14
    #>   n values:      6 per-cell sample sizes visited
    #>   sims per cell size: 1000
    #>   SS type:       III
    #>   n needed for between-subjects cell: 17
    #>   total N needed: 34
    #> 
    #>  n_per_cell total_n n_sims num_df den_df    ncp power_calc power_sim
    #>          13      26   1000      3     72 11.721      0.808     0.795
    #>          16      32   1000      3     90 14.651      0.897     0.885
    #>          17      34   1000      3     96 15.628      0.917     0.912
    #>          18      36   1000      3    102 16.605      0.934     0.936
    #>          20      40   1000      3    114 18.558      0.958     0.946
    #>          26      52   1000      3    150 24.419      0.991     0.991

[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
reports the required number of subjects per between-subject cell and the
corresponding total N. For a purely within-subject design, the reported
`n_per_cell` is the total sample size.

The output table uses compact column names: `num_df` and `den_df` are
the ANOVA degrees of freedom, `ncp` is the noncentrality parameter,
`power_calc` is the calculated noncentral-F power, and `power_sim` is
the simulation estimate.

The example uses 1000 simulations for speed. Use at least 5000
simulations for a more stable estimate; the package default is 10000.

### Adding factors and levels

Add as many factors and levels as required and name the term to test.
For example, this design includes a three-level between-subject factor
and tests a three-way interaction:

``` r

power_n(
  between = c(cond = 2, age = 3),
  within = c(stim = 4),
  term = "cond:stim:age",
  target_pes = 0.14,
  power = 0.90,
  n_sims = 5000,
  seed = 123
)
```

## Simulate a power curve

Use
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md)
to estimate power across explicitly chosen sample sizes. The result is a
tidy table that can be passed to
[`plot_power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/plot_power_curve.md).

``` r

pc <- power_curve(
  between = c(cond = 2),
  within = c(stim = 2),
  term = "cond:stim",
  target_pes = 0.14,
  n_range = c(16, 20, 23, 28),
  n_sims = 1000,
  seed = 123
)
pc
```

    #> <anovapowersim_curve>
    #>   term:          'cond:stim'
    #>   target power:  <not specified>
    #>   alpha:         0.05
    #>   effect size:   pes = 0.14
    #>   n values:      4 per-cell sample sizes visited
    #>   sims per cell size: 1000
    #>   SS type:       III
    #>   n needed for between-subjects cell: <not reached>
    #>   total N needed: <not reached>
    #> 
    #>  n_per_cell total_n n_sims num_df den_df   ncp power_calc power_sim
    #>          16      32   1000      1     30 4.884      0.571     0.557
    #>          20      40   1000      1     38 6.186      0.679     0.666
    #>          23      46   1000      1     44 7.163      0.745     0.735
    #>          28      56   1000      1     54 8.791      0.829     0.831

``` r

plot_power_curve(
  pc,
  power_lines = c(0.80, 0.90)
)
```

![](anovapowersim_files/figure-html/plot-fixed-1.png)

## Run simulations in parallel

For larger simulation runs, set `parallel = TRUE`. If `cores` is
omitted, `anovapowersim` uses one fewer than the available cores and
reports the number selected. Set `cores` explicitly when a fixed number
of workers is required.

``` r

power_n(
  between = c(cond = 2),
  within = c(stim = 2),
  term = "cond:stim",
  target_pes = 0.14,
  power = 0.90,
  n_sims = 5000,
  parallel = TRUE,
  cores = 4,
  seed = 123
)
```

## Next steps

- [Fixed-sample power and
  sensitivity](https://shaheedazaad.github.io/anovapowersim/articles/fixed-sample-power.md)
- [Calculated power without
  simulations](https://shaheedazaad.github.io/anovapowersim/articles/calculated-power.md)
- [Power for unbalanced
  designs](https://shaheedazaad.github.io/anovapowersim/articles/unbalanced-designs.md)
- [Covariance and
  nonsphericity](https://shaheedazaad.github.io/anovapowersim/articles/covariance.md)
- [Comparison with
  G\*Power](https://shaheedazaad.github.io/anovapowersim/articles/comparison-with-gpower.md)
