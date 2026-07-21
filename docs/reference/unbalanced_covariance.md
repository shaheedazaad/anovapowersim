# Specify correlations for a means-based unbalanced design

Defines the common within-subject correlation structure used by
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md).
Marginal standard deviations are deliberately absent: they come from
[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)
and may differ between groups.

## Usage

``` r
unbalanced_covariance(default_correlation = 0.5, correlations = NULL)
```

## Arguments

- default_correlation:

  Correlation in `(-1, 1)` used for unlisted pairs.

- correlations:

  Optional named numeric vector of pair-specific correlations. Name
  pairs as `"cell1:cell2"`; pair order does not matter.

## Value

An `anovapowersim_unbalanced_covariance_spec` object.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`unbalanced_covariance()` is experimental and is available only in the
development version of `anovapowersim`. Its API may change.

## Examples

``` r
unbalanced_covariance(
  default_correlation = 0.5,
  correlations = c("pre:post" = 0.7)
)
#> $default_correlation
#> [1] 0.5
#> 
#> $correlations
#> pre:post 
#>      0.7 
#> 
#> $correlation_pairs
#> # A tibble: 1 × 4
#>   pair_name cell1 cell2 pair_key   
#>   <chr>     <chr> <chr> <chr>      
#> 1 pre:post  pre   post  "post\rpre"
#> 
#> attr(,"class")
#> [1] "anovapowersim_unbalanced_covariance_spec"
```
