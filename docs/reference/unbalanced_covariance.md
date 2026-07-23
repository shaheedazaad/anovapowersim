# Specify covariance for a means-based unbalanced design

Defines the common marginal standard deviation and within-subject
correlation structure used by
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md).
The resulting covariance is shared by every between-subject cell.

## Usage

``` r
unbalanced_covariance(sd = 1, default_correlation = 0.5, correlations = NULL)
```

## Arguments

- sd:

  Common positive finite marginal standard deviation. If omitted, `1` is
  used and a warning is issued.

- default_correlation:

  Correlation in `(-1, 1)` used for unlisted pairs. When the
  specification is resolved for a design, a warning identifies how many
  pairs were not defined and makes clear that this default applies only
  to those pairs.

- correlations:

  Optional named numeric vector of pair-specific correlations. Name
  pairs as `"cell1:cell2"`; pair order does not matter. For multiple
  within factors, cell names join their level values with `_`.
  Constructed names must be unique, and level values must not contain
  `:` because it separates the two cells in a pair name.

## Value

An `anovapowersim_unbalanced_covariance_spec` object.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`unbalanced_covariance()` is experimental and is available only in the
development version of `anovapowersim`. Its API may change.

## Examples

``` r
unbalanced_covariance(
  sd = 2,
  default_correlation = 0.5,
  correlations = c("pre:post" = 0.7)
)
#> $sd
#> [1] 2
#> 
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
