# Build calibrated default means for a design term

Creates the default contrast pattern for one ANOVA term and scales it so
an exact reference dataset has the requested partial eta squared under
the supplied balanced design assumptions.

## Usage

``` r
design_term_means(design, term, target_pes, n, sd = 1, r = 0.5, gpower = FALSE)
```

## Arguments

- design:

  An `anovapowersim_design_spec` from
  [`balanced_anova_design()`](https://shaheedazaad.github.io/anovapowersim/reference/balanced_anova_design.md).

- term:

  Character scalar naming the ANOVA term to target. Interaction terms
  are order-insensitive.

- target_pes:

  Target partial eta squared.

- n:

  Sample size per between-subject cell. For pure within designs, this is
  the total sample size.

- sd:

  Common outcome standard deviation.

- r:

  Compound-symmetric correlation among within-subject cells.

- gpower:

  Logical; if `TRUE`, calibrate to the G\*Power-style noncentrality
  convention `lambda = total_n * f^2`.

## Value

A numeric matrix of cell means, with rows indexing between cells and
columns indexing within cells.

## Examples

``` r
d <- balanced_anova_design(between = c(group = 2), within = c(time = 2))
design_term_means(d, term = "group:time", target_pes = 0.2, n = 20)
#>            [,1]       [,2]
#> [1,]  0.2436699 -0.2436699
#> [2,] -0.2436699  0.2436699
```
