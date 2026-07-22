# Calculate ANOVA effect-size sensitivity at a fixed sample size

Calculation-only counterpart to
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md).
The function searches for the minimum partial eta squared that reaches
target power using calculated noncentral-F power, without simulating
data or fitting ANOVA models.

## Usage

``` r
power_sensitivity_calc(
  between = NULL,
  within = NULL,
  term,
  n,
  power = 0.9,
  alpha = 0.05,
  pes_min = 1e-06,
  pes_max = 0.99,
  pes_tol = 0.001,
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

- n:

  Number of subjects per between-subject cell. For a purely
  within-subject design, this is the total sample size.

- power:

  Desired target power.

- alpha:

  Significance threshold.

- pes_min:

  Lower bound of the partial eta-squared search interval.

- pes_max:

  Upper bound of the partial eta-squared search interval.

- pes_tol:

  Maximum width of the final calculated partial eta-squared bracket.

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

An `anovapowersim_sensitivity` object. `$pes_needed` is the calculated
upper effect-size bracket, or `NA` when `pes_max` does not achieve
target power. `$results` contains every effect size evaluated by the
calculated-power search; `n_sims` and `power_sim` are always `NA`.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`power_sensitivity_calc()` is experimental and is available only in the
development version of `anovapowersim`. Its API and reporting format may
change.

## Examples

``` r
power_sensitivity_calc(
  between = c(group = 2),
  within = c(time = 3),
  term = "group:time",
  n = 30,
  power = 0.90,
  pes_tol = 0.001,
  gpower = TRUE,
  epsilon = 0.80
)
#> Warning: `gpower = TRUE` for term 'group:time' does not calibrate `target_pes` to the partial eta squared you will actually observe in a fitted ANOVA, because its within-subject component has 2 degrees of freedom (more than one). This mirrors a property of G*Power's own 'as in Cohen (1988)' repeated-measures convention, which does not adjust for the number of measurements. Use the default `gpower = FALSE` if you want `target_pes` to match your reported or expected partial eta squared exactly.
#> <anovapowersim_sensitivity>
#>   term:             'group:time'
#>   fixed n per cell: 30
#>   fixed total N:    60
#>   target power:     0.900
#>   alpha:            0.05
#>   detectable pes:   0.203995
#>   search interval:  [1e-06, 0.99]
#>   requested width:  0.001
#>   calculated points: 12
#>   final width:      0.0009667959
#>   converged:        yes
#>   calculation:      calculated power only
#>   G*Power convention: TRUE
#>   epsilon:          0.8
#> 
#>  target_pes n_per_cell total_n n_sims epsilon num_df den_df      ncp power_calc
#>    0.000001         30      60     NA     0.8    1.6   92.8    0.000      0.050
#>    0.123751         30      60     NA     0.8    1.6   92.8    6.779      0.662
#>    0.185626         30      60     NA     0.8    1.6   92.8   10.941      0.863
#>    0.201095         30      60     NA     0.8    1.6   92.8   12.082      0.896
#>    0.203028         30      60     NA     0.8    1.6   92.8   12.228      0.899
#>    0.203995         30      60     NA     0.8    1.6   92.8   12.301      0.901
#>    0.204962         30      60     NA     0.8    1.6   92.8   12.374      0.903
#>    0.208829         30      60     NA     0.8    1.6   92.8   12.670      0.910
#>    0.216563         30      60     NA     0.8    1.6   92.8   13.269      0.922
#>    0.247501         30      60     NA     0.8    1.6   92.8   15.787      0.959
#>    0.495000         30      60     NA     0.8    1.6   92.8   47.050      1.000
#>    0.990000         30      60     NA     0.8    1.6   92.8 4752.000      1.000
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
