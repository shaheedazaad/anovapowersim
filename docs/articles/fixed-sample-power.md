# Fixed-sample power and sensitivity

``` r

library(anovapowersim)
```

[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md)
and
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
are experimental and available only in the development version of
`anovapowersim`. They use the same balanced design construction, ANOVA
fitting, covariance handling, and simulation controls as
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md).

## Achieved power

Use
[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md)
when the sample size and partial eta squared are fixed. The `n` argument
is the number of subjects per between-subject cell, or total N for a
purely within-subject design.

``` r

power_achieved(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  target_pes = 0.14,
  n = 20,
  n_sims = 5000,
  parallel = TRUE,
  seed = 123
)
```

The result reports simulated achieved power as the primary estimate and
calculated power as a diagnostic.

## Sensitivity

Use
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
when the sample size is fixed and the minimum detectable partial eta
squared is unknown. It simulates effect sizes until it finds a bracket
around the requested power, then reports an explicitly simulated upper
bracket as `pes_needed`.

``` r

power_sensitivity(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  n = 20,
  power = 0.90,
  n_sims = 5000,
  pes_tol = 0.001,
  parallel = TRUE,
  seed = 123
)
```

Because simulated power has Monte Carlo variability, use enough
simulations for the precision you need and inspect all visited effect
sizes in `$results`.

## Covariance and calculated-power alternatives

Both functions accept a
[`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
object for repeated-measures designs. See [Covariance and
nonsphericity](https://shaheedazaad.github.io/anovapowersim/articles/covariance.md)
for cell naming, custom standard deviations and correlations, and
Greenhouse–Geisser handling.

To skip simulations, use
[`power_achieved_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved_calc.md)
or
[`power_sensitivity_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity_calc.md).
These are covered in the [Calculated
power](https://shaheedazaad.github.io/anovapowersim/articles/calculated-power.md)
guide.
