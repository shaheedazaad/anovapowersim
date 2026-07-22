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
#> Warning: `power_sim` and `power_calc` differ by more than 5 percentage points for n_per_cell = 45 (largest difference = 0.061). Try increasing `n_sims` for a more stable simulation estimate.
#> <anovapowersim_curve>
#>   term:          'time:condition'
#>   target power:  0.900
#>   alpha:         0.05
#>   effect size:   pes = 0.14
#>   n values:      6 per-cell sample sizes visited
#>   sims per cell size: 100
#>   epsilon:       0.9795918
#>   covariance:    custom 6 x 6 within-subject matrix
#>   SS type:       III
#>   n needed for between-subjects cell: 46
#>   total N needed: 46
#> 
#>  n_per_cell total_n n_sims valid_sims failed_sims   epsilon   num_df    den_df
#>          33      33    100        100           0 0.9795918 1.959184  62.69388
#>          45      45    100        100           0 0.9795918 1.959184  86.20408
#>          46      46    100        100           0 0.9795918 1.959184  88.16327
#>          50      50    100        100           0 0.9795918 1.959184  96.00000
#>          52      52    100        100           0 0.9795918 1.959184  99.91837
#>          66      66    100        100           0 0.9795918 1.959184 127.34694
#>     ncp power_calc power_sim
#>  10.206      0.807     0.770
#>  14.033      0.921     0.860
#>  14.352      0.927     0.910
#>  15.628      0.947     0.950
#>  16.266      0.955     0.920
#>  20.731      0.987     1.000
# }
```
