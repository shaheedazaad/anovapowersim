# Simulate ANOVA power from a balanced factorial design

Simulation-based power estimation for balanced factorial designs under
sphericity. Users specify the between- and within-subject factors, the
ANOVA term to test, a target partial eta squared, and explicit sample
sizes. The function creates a default contrast pattern for the target
term, scales it to the requested partial eta squared, simulates
datasets, refits [`stats::aov()`](https://rdrr.io/r/stats/aov.html), and
estimates power by counting `p < alpha`.

## Usage

``` r
power_curve(
  between = NULL,
  within = NULL,
  term,
  target_pes,
  n_range,
  n_sims = 10000,
  alpha = 0.05,
  gpower = FALSE,
  progress = interactive(),
  seed = NULL
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

- n_range:

  Integer vector of sample sizes per between-subject cell. For pure
  within-subject designs, this is the total sample size.

- n_sims:

  Number of simulated datasets per sample size.

- alpha:

  Significance threshold.

- gpower:

  Logical; if `TRUE`, calibrate means to the G\*Power-style
  noncentrality convention `lambda = total_n * f^2`. The default `FALSE`
  calibrates the empirical reference dataset to `target_pes`, equivalent
  to `lambda = den_df * f^2` for the fitted ANOVA.

- progress:

  Logical; if `TRUE`, show a text progress bar.

- seed:

  Optional integer seed for reproducibility.

## Value

An `anovapowersim_curve` object. The `$results` tibble contains
`n_per_cell`, `total_n`, `n_sims`, numerator and denominator degrees of
freedom (`num_df`, `den_df`), the noncentrality parameter (`ncp`),
calculated power (`power_calc`), and simulated power (`power_sim`).

## Examples

    power_curve(
      between = c(cond = 2),
      within = c(stim = 2),
      term = "cond:stim",
      target_pes = 0.14,
      n_range = c(16, 20, 23, 28), # n per between-subject cell
      n_sims = 1000,
      seed = 123
    )
