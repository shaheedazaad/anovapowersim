# Specify a within-subject covariance structure

Creates a readable covariance specification for
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
and
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md).
Supply one common standard deviation and default correlation, then
override only the individual correlation pairs that differ. Measurement
names are generated from the `within` design; for example,
`within = c(time = 2, condition = 2)` creates `time1_condition1`,
`time1_condition2`, `time2_condition1`, and `time2_condition2`.

## Usage

``` r
within_covariance(sd = 1, default_correlation = 0.5, correlations = NULL)
```

## Arguments

- sd:

  Common positive finite standard deviation for every measurement. If
  omitted, `1` is used and a warning is issued.

- default_correlation:

  Finite correlation in `(-1, 1)` used for pairs not listed in
  `correlations`. When the specification is resolved for a design, a
  warning identifies how many pairs were not defined and makes clear
  that this default applies only to those pairs.

- correlations:

  Optional named numeric vector of pair-specific correlations. Name each
  pair as `"cell1:cell2"`, for example `"time1:time2" = 0.7`. Pair order
  does not matter.

## Value

An `anovapowersim_covariance_spec` object.
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
and
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md)
resolve it against their `within` design and construct the covariance
matrix as `sd^2 * R`, where `R` is the correlation matrix. They also
derive the tested term's population Greenhouse–Geisser epsilon from the
resolved matrix and apply it to `power_calc`.

## Examples

``` r
covariance <- within_covariance(
  sd = 1,
  default_correlation = 0.5,
  correlations = c(
    "time1_condition1:time1_condition2" = 0.6,
    "time2_condition1:time2_condition2" = 0.7
  )
)

# \donttest{
power_n(
  within = c(time = 3, condition = 2),
  term = "time:condition",
  target_pes = 0.14,
  power = 0.90,
  n_sims = 100,
  covariance = covariance,
  seed = 123
)
#> Warning: 13 of 15 within-subject correlations were not explicitly defined; `default_correlation = 0.5` is used only for those undefined pairs. Explicitly defined correlations are unchanged.
#> Warning: This multi-df within-subject term is nonspherical, so simulated power depends on the relative cell-mean pattern. `power_sim` currently uses the package's default linear/Kronecker pattern; a different pattern with the same `target_pes` and covariance can produce different power. Supply `means_pattern` to describe the expected relative mean shape.
#> Warning: Requested precision band was not reached for target power 0.900 with tolerance 0.030; reporting the closest explicitly simulated n_per_cell at or above target power: 47 (power_sim = 0.940).
#> <anovapowersim_curve>
#>   term:          'time:condition'
#>   target power:  0.900
#>   alpha:         0.05
#>   effect size:   pes = 0.14
#>   n values:      5 per-cell sample sizes visited
#>   sims per cell size: 100
#>   means pattern: default linear/Kronecker
#>   epsilon:       0.9795918
#>   covariance:    custom 6 x 6 within-subject matrix
#>   SS type:       III
#>   n needed for between-subjects cell: 47
#>   total N needed: 47
#> 
#>  n_per_cell total_n n_sims valid_sims failed_sims   epsilon   num_df    den_df
#>          33      33    100        100           0 0.9795918 1.959184  62.69388
#>          46      46    100        100           0 0.9795918 1.959184  88.16327
#>          47      47    100        100           0 0.9795918 1.959184  90.12245
#>          53      53    100        100           0 0.9795918 1.959184 101.87755
#>          66      66    100        100           0 0.9795918 1.959184 127.34694
#>     ncp power_calc power_sim
#>  10.206      0.807     0.770
#>  14.352      0.927     0.890
#>  14.671      0.933     0.940
#>  16.585      0.959     0.980
#>  20.731      0.987     0.990
# }
```
