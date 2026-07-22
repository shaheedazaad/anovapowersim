# Power for unbalanced designs

``` r

library(anovapowersim)
```

[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
is experimental and available only in the development version of
`anovapowersim`. Use it when the exact allocation, population means, and
standard deviations are already known. It simulates one fixed design; it
does not search over sample sizes or return calculated power.

## Define the cells

Define every cell with
[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md).
Factor values identify the cell, `n` is the number of subjects in its
between-subject group, and `m` is its population mean. Name
repeated-measures factors in `within`; all remaining factors are
between-subject factors. The common population SD is specified
separately with
[`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md).

For repeated-measures designs, the same `n` must appear on every
within-subject row belonging to a given between-subject group.

``` r

unbalanced_design <- cell_design(
  group = "control",   time = "pre",  n = 22, m = 10.0,
  group = "control",   time = "post", n = 22, m = 11.0,
  group = "treatment", time = "pre",  n = 31, m = 10.1,
  group = "treatment", time = "post", n = 31, m = 12.4,
  within = "time"
)
```

The means define the complete effect pattern, including main effects and
interactions. A single marginal SD is shared across groups and
within-subject cells, preventing unequal cell sizes from being combined
with unequal variances in the classical ANOVA test.

Every factor must have at least two observed levels, and every
factor-level combination must be defined exactly once. If cells are
missing,
[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)
reports which ones. Supply `default_n` and `default_m` together to fill
missing cells automatically:

``` r

cell_design(
  group = "control",   time = "pre",  n = 22, m = 10.0,
  group = "control",   time = "post", n = 22, m = 11.0,
  group = "treatment", time = "pre",  n = 31, m = 10.1,
  within = "time",
  default_n = 31, default_m = 12.4
)
```

## Define within-subject correlations

Use
[`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md)
to define the common SD and, for repeated measures, the correlation
structure.
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
constructs one covariance matrix and uses it for every between-subject
group.

If `covariance` is omitted,
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
warns that it is using the common `sd = 1` and, for repeated measures,
correlation `0.5` for every pair. Calling
[`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md)
without `sd` likewise warns that `sd = 1` is being used. When only some
correlations are named, a separate warning reports the number of
undefined pairs: `default_correlation` applies only to those undefined
pairs and does not alter explicitly supplied correlations. For a purely
between-subject design, correlations do not apply and the warning
mentions only the common SD.

With one within-subject factor, pair names use its levels, such as
`"pre:post"`.

``` r

power_unbalanced(
  design = unbalanced_design,
  term = "group:time",
  covariance = unbalanced_covariance(
    sd = 2,
    default_correlation = 0.5,
    correlations = c("pre:post" = 0.7)
  ),
  n_sims = 5000,
  parallel = TRUE,
  seed = 123
)
```

## Interpret the result

The result reports simulated power, the common SD, partial eta squared
from a deterministic reference dataset matching the design assumptions,
and the mean, median, and 95% interval of partial eta squared across
successful simulations. These effect-size summaries describe the exact
allocation, means, shared variance, correlations, tested term, and
sums-of-squares type supplied by the user.

## Multiple within-subject factors

Name every repeated-measures factor in
[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)
and join their level combinations with underscores in correlation-pair
names:

``` r

multi_within_design <- cell_design(
  group = "A", time = "pre",  cond = "control", n = 10, m = 0.0,
  group = "A", time = "pre",  cond = "treat",   n = 10, m = 0.5,
  group = "A", time = "post", cond = "control", n = 10, m = 0.2,
  group = "A", time = "post", cond = "treat",   n = 10, m = 1.0,
  group = "B", time = "pre",  cond = "control", n = 15, m = 0.0,
  group = "B", time = "pre",  cond = "treat",   n = 15, m = 0.6,
  group = "B", time = "post", cond = "control", n = 15, m = 0.3,
  group = "B", time = "post", cond = "treat",   n = 15, m = 1.4,
  within = c("time", "cond")
)

power_unbalanced(
  design = multi_within_design,
  term = "group:time:cond",
  covariance = unbalanced_covariance(
    correlations = c("pre_control:post_control" = 0.6)
  ),
  n_sims = 5000,
  seed = 123
)
```
