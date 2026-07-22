# Simulate data from a balanced ANOVA design

Generates one long-format dataset from a balanced design. Supply means
from
[`design_term_means()`](https://shaheedazaad.github.io/anovapowersim/reference/design_term_means.md)
or any conformable matrix with one row per between-subject cell and one
column per within-subject cell.

## Usage

``` r
simulate_design_dataset(design, n, means, sd = 1, r = 0.5, empirical = FALSE)
```

## Arguments

- design:

  An `anovapowersim_design_spec` from
  [`balanced_anova_design()`](https://shaheedazaad.github.io/anovapowersim/reference/balanced_anova_design.md).

- n:

  Sample size per between-subject cell. For pure within designs, this is
  the total sample size.

- means:

  Numeric matrix of population cell means.

- sd:

  Common outcome standard deviation.

- r:

  Compound-symmetric correlation among within-subject cells.

- empirical:

  Logical; if `TRUE`, use `MASS::mvrnorm(empirical = TRUE)` so the
  generated sample closely matches the requested means/covariance.

## Value

A tibble ready for [`stats::aov()`](https://rdrr.io/r/stats/aov.html)
with columns `id`, factor columns, and `value`.

## Covariance limitation

This manual helper does not accept
[`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
specifications. It always simulates from the compound-symmetric covariance
defined by `sd` and `r`. It therefore cannot reproduce a balanced
power-function call that uses a custom `covariance`; use the power functions
directly for that workflow.

## Examples

``` r
d <- balanced_anova_design(between = c(group = 2), within = c(time = 2))
m <- design_term_means(d, term = "group:time", target_pes = 0.2, n = 20)
sim <- simulate_design_dataset(d, n = 20, means = m)
head(sim)
#> # A tibble: 6 × 4
#>   id    group   value time 
#>   <fct> <fct>   <dbl> <fct>
#> 1 1     group1  1.78  time1
#> 2 1     group1  1.04  time2
#> 3 2     group1  1.18  time1
#> 4 2     group1 -0.840 time2
#> 5 3     group1  2.61  time1
#> 6 3     group1  0.639 time2
```
