# Calculate the sample size needed for target ANOVA power

Calculation-only search for the per-between-cell sample size needed to
reach a requested power for a balanced factorial ANOVA design. Unlike
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md),
this function does not run simulations, fit ANOVA models, or call `car`;
numerator degrees of freedom, denominator degrees of freedom,
noncentrality, and power are computed analytically from the balanced
design.

## Usage

``` r
power_n_calc(
  between = NULL,
  within = NULL,
  term,
  target_pes,
  power = 0.9,
  alpha = 0.05,
  n_start = NULL,
  n_max = 5000,
  gpower = FALSE
)
```

## Arguments

- between:

  Named integer vector of between-subject factor level counts, e.g.
  `c(group = 2)`. Use `NULL` for no between-subject factors.

- within:

  Named integer vector of within-subject factor level counts, e.g.
  `c(time = 3, condition = 4)`. Use `NULL` for no within-subject
  factors.

- term:

  Character scalar naming the ANOVA term to test, e.g. `"group:time"`.
  Interaction terms are order-insensitive; `"time:group"` resolves to
  `"group:time"` when that is the design's factor order.

- target_pes:

  Target partial eta squared for `term`.

- power:

  Desired target power.

- alpha:

  Significance threshold.

- n_start:

  Starting sample size per between-subject cell. If `NULL`, starts from
  the smallest analytically valid value.

- n_max:

  Maximum sample size per between-subject cell.

- gpower:

  Logical; if `TRUE`, use the G\*Power-style noncentrality convention
  `lambda = total_n * f^2`. The default `FALSE` uses
  `lambda = den_df * f^2`.

## Value

An `anovapowersim_curve` object with `n_needed` and `total_n_needed`.
The `$results` tibble contains `n_per_cell`, `total_n`, `n_sims`,
numerator and denominator degrees of freedom (`num_df`, `den_df`), the
noncentrality parameter (`ncp`), calculated power (`power_calc`), and
simulated power (`power_sim`). For `power_n_calc()`, `n_sims` and
`power_sim` are always `NA`.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`power_n_calc()` is experimental while the analytic search API and
reporting format are refined.

## Examples

``` r
power_n_calc(
  between = c(cond = 2),
  within = c(stim = 4),
  term = "cond:stim",
  target_pes = 0.14,
  power = 0.90
)
#> <anovapowersim_curve>
#>   term:          'cond:stim'
#>   target power:  0.900
#>   alpha:         0.05
#>   effect size:   pes = 0.14
#>   n values:      9 per-cell sample sizes visited
#>   calculation:   analytic only
#>   n needed for between-subjects cell: 17
#>   total N needed: 34
#> 
#>  n_per_cell total_n n_sims num_df den_df    ncp power_calc power_sim
#>           2       4     NA      3      6  0.977      0.085      <NA>
#>           4       8     NA      3     18  2.930      0.223      <NA>
#>           8      16     NA      3     42  6.837      0.535      <NA>
#>          16      32     NA      3     90 14.651      0.897      <NA>
#>          17      34     NA      3     96 15.628      0.917      <NA>
#>          18      36     NA      3    102 16.605      0.934      <NA>
#>          20      40     NA      3    114 18.558      0.958      <NA>
#>          24      48     NA      3    138 22.465      0.984      <NA>
#>          32      64     NA      3    186 30.279      0.998      <NA>
```
