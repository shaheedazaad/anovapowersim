# Define cells for a means-based unbalanced ANOVA design

Creates the complete cell table used by
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md).
Each cell is defined by its factor levels, sample size (`n`), and
population mean (`m`). End each cell after both reserved values have
been supplied. The common population standard deviation belongs in
[`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md).

## Usage

``` r
cell_design(..., within = NULL, default_n = NULL, default_m = NULL)
```

## Arguments

- ...:

  Repeated named cell definitions. Each cell must contain the same
  factor names in the same order, plus `n` and `m`. Every factor must
  have at least 2 observed levels, and every combination of factor
  levels must appear exactly once (or be filled automatically; see
  `default_n`).

- within:

  Character vector naming factors in `...` that are measured within
  subjects, or `NULL` for a purely between-subject design. Stored on the
  returned object and read by
  [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md).

- default_n, default_m:

  Optional scalars used to fill any missing cells in the complete
  factorial design. Supply both to auto-fill missing cells with these
  values; supply none to require every cell to be defined explicitly
  (the default). Supplying only one is an error.

## Value

An `anovapowersim_cell_design` tibble with one row per design cell.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`cell_design()` is experimental and is available only in the development
version of `anovapowersim`. Its API may change.

## Examples

``` r
design <- cell_design(
  group = "control", time = "pre",  n = 22, m = 10.0,
  group = "control", time = "post", n = 22, m = 11.0,
  group = "treatment", time = "pre",  n = 31, m = 10.1,
  group = "treatment", time = "post", n = 31, m = 12.4,
  within = "time"
)
design
#> # A tibble: 4 × 4
#>   group     time      n     m
#>   <chr>     <chr> <int> <dbl>
#> 1 control   pre      22  10  
#> 2 control   post     22  11  
#> 3 treatment pre      31  10.1
#> 4 treatment post     31  12.4
```
