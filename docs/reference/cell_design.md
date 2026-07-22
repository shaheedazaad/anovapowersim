# Define cells for a means-based unbalanced ANOVA design

Creates the complete cell table used by
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md).
Each cell is defined by its factor levels, sample size (`n`), population
mean (`m`), and population standard deviation (`sd`). End each cell
after all three reserved values have been supplied.

## Usage

``` r
cell_design(
  ...,
  within = NULL,
  default_n = NULL,
  default_m = NULL,
  default_sd = NULL
)
```

## Arguments

- ...:

  Repeated named cell definitions. Each cell must contain the same
  factor names in the same order, plus `n`, `m`, and `sd`. Every factor
  must have at least 2 observed levels, and every combination of factor
  levels must appear exactly once (or be filled automatically; see
  `default_n`).

- within:

  Character vector naming factors in `...` that are measured within
  subjects, or `NULL` for a purely between-subject design. Stored on the
  returned object and read by
  [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md).

- default_n, default_m, default_sd:

  Optional scalars used to fill any missing cells in the complete
  factorial design. Supply all three to auto-fill missing cells with
  these values; supply none to require every cell to be defined
  explicitly (the default). Supplying only some of the three is an
  error.

## Value

An `anovapowersim_cell_design` tibble with one row per design cell.

## Lifecycle

[![\[Experimental\]](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`cell_design()` is experimental and is available only in the development
version of `anovapowersim`. Its API may change.

## Examples

``` r
design <- cell_design(
  group = "control", time = "pre",  n = 22, m = 10.0, sd = 2.0,
  group = "control", time = "post", n = 22, m = 11.0, sd = 2.2,
  group = "treatment", time = "pre",  n = 31, m = 10.1, sd = 2.4,
  group = "treatment", time = "post", n = 31, m = 12.4, sd = 2.8,
  within = "time"
)
design
#> # A tibble: 4 × 5
#>   group     time      n     m    sd
#>   <chr>     <chr> <int> <dbl> <dbl>
#> 1 control   pre      22  10     2  
#> 2 control   post     22  11     2.2
#> 3 treatment pre      31  10.1   2.4
#> 4 treatment post     31  12.4   2.8
```
