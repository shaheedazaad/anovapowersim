# Search for the sample size needed for target ANOVA power

Adaptive simulation search for the per-between-cell sample size needed
to reach a requested power for a balanced factorial ANOVA design. The
search doubles upward from `n_start` until it brackets the target or
reaches `n_max`, then bisects the bracket.

## Usage

``` r
power_n(
  between = NULL,
  within = NULL,
  term,
  target_pes,
  power = 0.8,
  n_sims = 10000,
  alpha = 0.05,
  sd = 1,
  r = 0.5,
  n_start = NULL,
  n_max = 1000,
  tol = 0.01,
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

- power:

  Desired target power.

- n_sims:

  Number of simulated datasets per sample size.

- alpha:

  Significance threshold.

- sd:

  Common outcome standard deviation.

- r:

  Compound-symmetric correlation among within-subject cells.

- n_start:

  Starting sample size per between-subject cell. If `NULL`, starts at
  the smallest value that can support empirical calibration for the
  requested design.

- n_max:

  Maximum sample size per between-subject cell.

- tol:

  Stop when estimated power is within `tol` of `power`.

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

An `anovapowersim_curve` object with `n_needed` and `total_n_needed`.

## Examples

``` r
if (FALSE) { # \dontrun{
power_n(
  between = c(color = 2),
  within = c(age = 2),
  term = "age:color",
  target_pes = 0.20721,
  power = 0.90,
  n_sims = 10000
)
} # }
```
