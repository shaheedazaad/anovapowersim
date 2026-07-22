# Calculate achieved ANOVA power at a fixed sample size

Calculation-only counterpart to
[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md).
Numerator and denominator degrees of freedom, noncentrality, and
achieved power are calculated directly without simulating data or
fitting ANOVA models.

## Usage

``` r
power_achieved_calc(
  between = NULL,
  within = NULL,
  term,
  target_pes,
  n,
  alpha = 0.05,
  gpower = FALSE,
  epsilon = 1
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

- n:

  Number of subjects per between-subject cell. For a purely
  within-subject design, this is the total sample size.

- alpha:

  Significance threshold.

- gpower:

  Logical; if `TRUE`, use the G*Power-style noncentrality convention
  `lambda = total_n * f^2`. The default `FALSE` uses
  `lambda = den_df * f^2`. For a term whose within-subject component has
  more than one degree of freedom, `gpower`'s `target_pes` does not
  equal the partial eta squared actually achieved (this mirrors a
  property of G*Power's own "as in Cohen (1988)" repeated-measures
  convention, which does not adjust for the number of measurements); a
  warning is issued in that case. Use the default if you want
  `target_pes` to match your reported or expected partial eta squared
  exactly.

- epsilon:

  Population nonsphericity correction for the within-subject component
  of `term`. Must lie between the theoretical lower bound
  `1 / within_term_df` and `1`. The default `1` assumes sphericity.
  Values below `1` multiply the numerator degrees of freedom,
  denominator degrees of freedom, and noncentrality parameter.
  Nonsphericity corrections do not apply to purely between-subject
  terms.

## Value

An `anovapowersim_achieved_power` object. `$achieved_power` and
`$calculated_power` contain the calculated-power estimate. In
`$results`, `n_sims` and `power_sim` are `NA` because no simulations are
run.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`power_achieved_calc()` is experimental and is available only in the
development version of `anovapowersim`. Its API and reporting format may
change.

## Examples

``` r
power_achieved_calc(
  between = c(group = 2),
  within = c(time = 3),
  term = "group:time",
  target_pes = 0.08,
  n = 30,
  gpower = TRUE,
  epsilon = 0.80
)
#> Warning: `gpower = TRUE` for term 'group:time' does not calibrate `target_pes` to the partial eta squared you will actually observe in a fitted ANOVA, because its within-subject component has 2 degrees of freedom (more than one). This mirrors a property of G*Power's own 'as in Cohen (1988)' repeated-measures convention, which does not adjust for the number of measurements. Use the default `gpower = FALSE` if you want `target_pes` to match your reported or expected partial eta squared exactly.
#> <anovapowersim_achieved_power>
#>   term:             'group:time'
#>   fixed n per cell: 30
#>   fixed total N:    60
#>   target pes:       0.0800
#>   alpha:            0.05
#>   calculation:      calculated power only
#>   achieved power:   0.453 (calculated)
#>   calculated power: 0.453
#>   G*Power convention: TRUE
#>   epsilon:          0.8
#> 
#>  n_per_cell total_n n_sims epsilon num_df den_df   ncp power_calc power_sim
#>          30      60     NA     0.8    1.6   92.8 4.174      0.453      <NA>
```
