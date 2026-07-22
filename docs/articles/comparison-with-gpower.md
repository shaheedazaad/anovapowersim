# Comparison with G\*Power

``` r

library(anovapowersim)
```

The balanced simulation and calculation-only functions accept
`gpower = TRUE` when results should follow G\*Power’s noncentrality
convention. For within-subject designs, this corresponds to selecting
G\*Power’s **as in Cohen (1988)** option.

## A typical 2 x 2 design

In a 2 x 2 mixed design, the default and G\*Power conventions usually
give very similar results. Here they differ by one participant per
group:

``` r

default_2x2 <- power_n_calc(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  target_pes = 0.08,
  power = 0.90
)

gpower_2x2 <- power_n_calc(
  between = c(group = 2),
  within = c(time = 2),
  term = "group:time",
  target_pes = 0.08,
  power = 0.90,
  gpower = TRUE
)

c(default = default_2x2$n_needed, gpower = gpower_2x2$n_needed)
#> default  gpower 
#>      63      62
```

Set `gpower = TRUE` in
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md),
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md),
[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md),
or
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
in the same way when running simulations.

## More than three within-subject levels

The conventions can diverge substantially when a within-subject term has
more than three levels. With four within-subject levels, the same inputs
produce:

``` r

default_4 <- power_n_calc(
  between = c(group = 2),
  within = c(time = 4),
  term = "group:time",
  target_pes = 0.08,
  power = 0.90
)

gpower_4 <- power_n_calc(
  between = c(group = 2),
  within = c(time = 4),
  term = "group:time",
  target_pes = 0.08,
  power = 0.90,
  gpower = TRUE
)

c(default = default_4$n_needed, gpower = gpower_4$n_needed)
#> default  gpower 
#>      29      83
```

In this case, use the package default, `gpower = FALSE`, so `target_pes`
continues to match the partial eta squared you supplied. Use
`gpower = TRUE` only when reproducing G\*Power’s **as in Cohen (1988)**
result is specifically required.
