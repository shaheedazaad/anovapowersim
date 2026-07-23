# Calculate the sample size needed for target ANOVA power

Calculation-only search for the per-between-cell sample size needed to
reach a requested power for a balanced factorial ANOVA design. Unlike
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md),
this function does not run simulations, fit ANOVA models, or call `car`;
numerator degrees of freedom, denominator degrees of freedom,
noncentrality, and calculated power are obtained directly from the
balanced design.

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

- power:

  Desired target power.

- alpha:

  Significance threshold.

- n_start:

  Starting sample size per between-subject cell, not a lower bound for
  the search. If `NULL`, starts from the smallest value with valid
  calculated-power degrees of freedom.

- n_max:

  Maximum sample size per between-subject cell.

- gpower:

  Logical; if `TRUE`, use the G*Power-style noncentrality convention
  `lambda = total_n * f^2`. The default `FALSE` uses
  `lambda = den_df * f^2`. G*Power's estimates can differ from
  `target_pes`, especially for small samples or terms with more degrees
  of freedom; a warning is issued when `gpower = TRUE`. The default
  `gpower = FALSE` is recommended.

- epsilon:

  Population nonsphericity correction for the within-subject component
  of `term`. Must lie between the theoretical lower bound
  `1 / within_term_df` and `1`. The default `1` assumes sphericity.
  Values below `1` multiply the numerator degrees of freedom,
  denominator degrees of freedom, and noncentrality parameter.
  Nonsphericity corrections do not apply to purely between-subject
  terms.

## Value

An `anovapowersim_curve` object with `n_needed` and `total_n_needed`.
The `$results` tibble contains `n_per_cell`, `total_n`, `n_sims`,
`valid_sims`, `failed_sims`, numerator and denominator degrees of
freedom (`num_df`, `den_df`), the nonsphericity correction (`epsilon`),
the noncentrality parameter (`ncp`), calculated power (`power_calc`),
and simulated power (`power_sim`). For `power_n_calc()`, the
simulation-specific columns are always `NA`. When `epsilon < 1`,
`num_df` and `den_df` are the corrected degrees of freedom used in the
power calculation.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`power_n_calc()` is experimental while the calculated-power search API
and reporting format are refined.

## Examples

``` r
power_n_calc(
  between = c(cond = 2),
  within = c(stim = 4),
  term = "cond:stim",
  target_pes = 0.14,
  power = 0.90,
  epsilon = 0.70
)
#> <anovapowersim_curve>
#>   term:          'cond:stim'
#>   target power:  0.900
#>   alpha:         0.05
#>   effect size:   pes = 0.14
#>   n values:      9 per-cell sample sizes visited
#>   calculation:   calculated power only
#>   epsilon:       0.7
#>   n needed for between-subjects cell: 21
#>   total N needed: 42
#> 
#>  n_per_cell total_n n_sims valid_sims failed_sims epsilon num_df den_df    ncp
#>           2       4     NA         NA          NA     0.7    2.1    4.2  0.684
#>           4       8     NA         NA          NA     0.7    2.1   12.6  2.051
#>           8      16     NA         NA          NA     0.7    2.1   29.4  4.786
#>          16      32     NA         NA          NA     0.7    2.1   63.0 10.256
#>          20      40     NA         NA          NA     0.7    2.1   79.8 12.991
#>          21      42     NA         NA          NA     0.7    2.1   84.0 13.674
#>          22      44     NA         NA          NA     0.7    2.1   88.2 14.358
#>          24      48     NA         NA          NA     0.7    2.1   96.6 15.726
#>          32      64     NA         NA          NA     0.7    2.1  130.2 21.195
#>  power_calc power_sim
#>       0.077      <NA>
#>       0.185      <NA>
#>       0.436      <NA>
#>       0.799      <NA>
#>       0.892      <NA>
#>       0.908      <NA>
#>       0.922      <NA>
#>       0.945      <NA>
#>       0.987      <NA>
```
