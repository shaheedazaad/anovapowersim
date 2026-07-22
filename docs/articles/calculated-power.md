# Calculated power without simulations

``` r

library(anovapowersim)
```

The calculation-only functions provide fast power calculations for
balanced factorial ANOVA designs. They use the noncentral F distribution
directly and do not simulate datasets or fit ANOVA models.

These functions are experimental and available only in the development
version of `anovapowersim`:

- [`power_n_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n_calc.md)
  finds the sample size required for a target effect size and power.
- [`power_achieved_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved_calc.md)
  calculates power for a fixed sample size and effect size.
- [`power_sensitivity_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity_calc.md)
  finds the minimum detectable partial eta squared for a fixed sample
  size and target power.

The result objects use `power_calc` for calculated power. Their `n_sims`
and `power_sim` fields are `NA` because no simulations are run.

## Required sample size

Use
[`power_n_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n_calc.md)
when the sample size is unknown. As with
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md),
the reported `n_needed` is the number of subjects per between-subject
cell. For a purely within-subject design, it is the total sample size.

``` r

power_n_calc(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  target_pes = 0.08,
  power = 0.90,
  alpha = 0.05
)
#> <anovapowersim_curve>
#>   term:          'group:time'
#>   target power:  0.900
#>   alpha:         0.05
#>   effect size:   pes = 0.08
#>   n values:      11 per-cell sample sizes visited
#>   calculation:   calculated power only
#>   n needed for between-subjects cell: 63
#>   total N needed: 126
#> 
#>  n_per_cell total_n n_sims epsilon num_df den_df    ncp power_calc power_sim
#>           2       4     NA       1      1      2  0.174      0.058      <NA>
#>           4       8     NA       1      1      6  0.522      0.094      <NA>
#>           8      16     NA       1      1     14  1.217      0.177      <NA>
#>          16      32     NA       1      1     30  2.609      0.346      <NA>
#>          32      64     NA       1      1     62  5.391      0.628      <NA>
#>          48      96     NA       1      1     94  8.174      0.808      <NA>
#>          56     112     NA       1      1    110  9.565      0.866      <NA>
#>          60     120     NA       1      1    118 10.261      0.888      <NA>
#>          62     124     NA       1      1    122 10.609      0.898      <NA>
#>          63     126     NA       1      1    124 10.783      0.903      <NA>
#>          64     128     NA       1      1    126 10.957      0.907      <NA>
```

## Achieved power

Use
[`power_achieved_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved_calc.md)
when both the sample size and partial eta squared are fixed.

``` r

power_achieved_calc(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  target_pes = 0.08,
  n = 30,
  alpha = 0.05
)
#> <anovapowersim_achieved_power>
#>   term:             'group:time'
#>   fixed n per cell: 30
#>   fixed total N:    60
#>   target pes:       0.0800
#>   alpha:            0.05
#>   calculation:      calculated power only
#>   achieved power:   0.598 (calculated)
#>   calculated power: 0.598
#> 
#>  n_per_cell total_n n_sims epsilon num_df den_df   ncp power_calc power_sim
#>          30      60     NA       1      1     58 5.043      0.598      <NA>
```

## Sensitivity

Use
[`power_sensitivity_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity_calc.md)
when the sample size is fixed but the minimum detectable effect size is
unknown. The result’s `pes_needed` is the calculated upper bracket that
reaches the requested power.

``` r

power_sensitivity_calc(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  n = 30,
  power = 0.90,
  pes_tol = 0.001,
  alpha = 0.05
)
#> <anovapowersim_sensitivity>
#>   term:             'group:time'
#>   fixed n per cell: 30
#>   fixed total N:    60
#>   target power:     0.900
#>   alpha:            0.05
#>   detectable pes:   0.158556
#>   search interval:  [1e-06, 0.99]
#>   requested width:  0.001
#>   calculated points: 12
#>   final width:      0.0009667959
#>   converged:        yes
#>   calculation:      calculated power only
#> 
#>  target_pes n_per_cell total_n n_sims epsilon num_df den_df      ncp power_calc
#>    0.000001         30      60     NA       1      1     58    0.000      0.050
#>    0.123751         30      60     NA       1      1     58    8.191      0.804
#>    0.154688         30      60     NA       1      1     58   10.614      0.893
#>    0.156622         30      60     NA       1      1     58   10.771      0.897
#>    0.157589         30      60     NA       1      1     58   10.850      0.900
#>    0.158556         30      60     NA       1      1     58   10.929      0.902
#>    0.162423         30      60     NA       1      1     58   11.247      0.910
#>    0.170157         30      60     NA       1      1     58   11.893      0.924
#>    0.185626         30      60     NA       1      1     58   13.220      0.947
#>    0.247501         30      60     NA       1      1     58   19.076      0.990
#>    0.495000         30      60     NA       1      1     58   56.852      1.000
#>    0.990000         30      60     NA       1      1     58 5742.000      1.000
#>  power_sim
#>       <NA>
#>       <NA>
#>       <NA>
#>       <NA>
#>       <NA>
#>       <NA>
#>       <NA>
#>       <NA>
#>       <NA>
#>       <NA>
#>       <NA>
#>       <NA>
```

## Planned nonsphericity

For a term containing a within-subject factor, `epsilon` supplies a
planned population nonsphericity correction. It defaults to `1`, which
assumes sphericity. Values below `1` reduce the numerator degrees of
freedom, denominator degrees of freedom, and noncentrality parameter
used in the power calculation.

The value must be at least the theoretical lower bound for the tested
term, `1 / within_term_df`. It must remain `1` for a purely
between-subject term.

``` r

power_n_calc(
  between = c(group = 2),
  within = c(time = 3),
  term = "group:time",
  target_pes = 0.08,
  power = 0.90,
  epsilon = 0.70
)
#> <anovapowersim_curve>
#>   term:          'group:time'
#>   target power:  0.900
#>   alpha:         0.05
#>   effect size:   pes = 0.08
#>   n values:      11 per-cell sample sizes visited
#>   calculation:   calculated power only
#>   epsilon:       0.7
#>   n needed for between-subjects cell: 50
#>   total N needed: 100
#> 
#>  n_per_cell total_n n_sims epsilon num_df den_df    ncp power_calc power_sim
#>           2       4     NA     0.7    1.4    2.8  0.243      0.061      <NA>
#>           4       8     NA     0.7    1.4    8.4  0.730      0.104      <NA>
#>           8      16     NA     0.7    1.4   19.6  1.704      0.205      <NA>
#>          16      32     NA     0.7    1.4   42.0  3.652      0.412      <NA>
#>          32      64     NA     0.7    1.4   86.8  7.548      0.730      <NA>
#>          48      96     NA     0.7    1.4  131.6 11.443      0.893      <NA>
#>          49      98     NA     0.7    1.4  134.4 11.687      0.900      <NA>
#>          50     100     NA     0.7    1.4  137.2 11.930      0.906      <NA>
#>          52     104     NA     0.7    1.4  142.8 12.417      0.917      <NA>
#>          56     112     NA     0.7    1.4  154.0 13.391      0.936      <NA>
#>          64     128     NA     0.7    1.4  176.4 15.339      0.962      <NA>
```

Calculation-only functions accept `epsilon`, not a
[`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
object. Use the simulation-based functions when you want to specify a
complete within-subject covariance structure and generate data from it.

## Comparison with G\*Power

The `gpower` option is shared by the balanced simulation and
calculation-only functions. See the short [Comparison with
G\*Power](https://shaheedazaad.github.io/anovapowersim/articles/comparison-with-gpower.md)
guide for matching G\*Power and the case where the default should be
preferred.

## When to use simulations instead

Calculated power is useful for quick planning, sensitivity checks, and
comparison with other software. Prefer the simulation-based functions
when you need:

- a complete within-subject covariance structure;
- simulated ANOVA fitting and Greenhouse–Geisser-adjusted tests;
- direct comparison between calculated and simulated power; or
- exact unequal sample sizes, cell means, and standard deviations
  through
  [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md).
