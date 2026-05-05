# Create a balanced factorial ANOVA design specification

Builds the design object used by
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md),
[`design_term_means()`](https://shaheedazaad.github.io/anovapowersim/reference/design_term_means.md),
and
[`simulate_design_dataset()`](https://shaheedazaad.github.io/anovapowersim/reference/simulate_design_dataset.md).
This object stores factor names, level counts, generated factor levels,
and the between/within cell grids.

## Usage

``` r
balanced_anova_design(between = NULL, within = NULL)
```

## Arguments

- between:

  Named integer vector of between-subject factor level counts, e.g.
  `c(group = 2)`. Use `NULL` for no between-subject factors.

- within:

  Named integer vector of within-subject factor level counts, e.g.
  `c(time = 3)`. Use `NULL` for no within-subject factors.

## Value

An object of class `anovapowersim_design_spec`.

## Examples

``` r
d <- balanced_anova_design(between = c(group = 2), within = c(time = 3))
d$between_cells
d$within_cells
```
