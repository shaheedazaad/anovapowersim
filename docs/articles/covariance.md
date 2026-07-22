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
use a common standard deviation of 1 and a correlation of 0.5 between
every within-subject pair by default. All four accept
[`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
to change that common SD or specify pair-specific correlations.

## Name the within-subject cells

Cell names combine within-factor levels with underscores. For
`within = c(time = 3, condition = 2)`, the expected names are
`time1_condition1`, `time1_condition2`, `time2_condition1`,
`time2_condition2`, `time3_condition1`, and `time3_condition2`.

Correlation names join a pair of cells with `:`.

``` r

covariance <- within_covariance(
  sd = 1,
  default_correlation = 0.5,
  correlations = c(
    "time1_condition1:time1_condition2" = 0.6,
    "time2_condition1:time2_condition2" = 0.7
  )
)
```

Every measurement uses the common `sd`; unlisted pairs use
`default_correlation`. Pair-specific correlations can produce
nonsphericity even though all diagonal variances are equal.

The package warns whenever these defaults are actually needed:

- Omitting `covariance` from a balanced simulation uses `sd = 1` and a
  correlation of `0.5` for every within-subject pair.
- Calling
  [`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
  without `sd` uses the common `sd = 1`.
- If `correlations` does not name every pair, the warning reports how
  many pairs are undefined. `default_correlation` fills only those
  undefined pairs; correlations explicitly listed by the user are never
  replaced.

A raw covariance matrix defines both variances and correlations directly
and does not use these defaults. Its diagonal variances must all be
equal.

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
the common SD and correlations belong in
[`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md);
individual cells specify only sample sizes and means. See [Power for
unbalanced
designs](https://shaheedazaad.github.io/anovapowersim/articles/unbalanced-designs.md).
