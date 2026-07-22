# Covariance and nonsphericity

``` r

library(anovapowersim)
```

The balanced simulation functions
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md),
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md),
[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md),
and
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
use a standard deviation of 1 for every within-subject cell and a
correlation of 0.5 between every pair by default. All four accept
[`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
when repeated measurements have different standard deviations or
correlations.

## Name the within-subject cells

Cell names combine within-factor levels with underscores. For
`within = c(time = 3, condition = 2)`, the expected names are
`time1_condition1`, `time1_condition2`, `time2_condition1`,
`time2_condition2`, `time3_condition1`, and `time3_condition2`.

Correlation names join a pair of cells with `:`.

``` r

covariance <- within_covariance(
  default_sd = 1,
  default_correlation = 0.5,
  standard_deviations = c(
    "time3_condition2" = 1.2
  ),
  correlations = c(
    "time1_condition1:time1_condition2" = 0.6,
    "time2_condition1:time2_condition2" = 0.7
  )
)
```

Unlisted measurements use `default_sd`, and unlisted pairs use
`default_correlation`.

## Use the covariance specification

Pass the object to any balanced simulation function:

``` r

power_n(
  within = c(time = 3, condition = 2),
  term = "time:condition",
  target_pes = 0.14,
  power = 0.90,
  n_sims = 5000,
  covariance = covariance,
  seed = 123
)
```

The custom covariance structure calibrates the target effect and
generates every simulated dataset. The package projects the covariance
into the tested term’s contrast space and calculates a population
Greenhouse–Geisser epsilon. `power_calc` applies this epsilon to its
degrees of freedom and noncentrality.

When the population epsilon is below 1, `power_sim` uses each simulated
dataset’s sample-estimated Greenhouse–Geisser-adjusted p-value for
sums-of- squares type II or III. Type I tests remain uncorrected and
issue a warning. Huynh–Feldt corrections are not currently implemented.

## Calculation-only and unbalanced designs

Calculation-only functions accept a planned `epsilon` rather than a
[`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
object; see [Calculated
power](https://shaheedazaad.github.io/anovapowersim/articles/calculated-power.md).

For
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md),
standard deviations belong in
[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md),
so its covariance argument accepts
[`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md)
correlations only. See [Power for unbalanced
designs](https://shaheedazaad.github.io/anovapowersim/articles/unbalanced-designs.md).
