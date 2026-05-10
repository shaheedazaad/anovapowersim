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
  n_start = NULL,
  n_max = 1000,
  tol = 0.01,
  gpower = FALSE,
  progress = interactive(),
  parallel = FALSE,
  cores = NULL,
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

- parallel:

  Logical; if `TRUE`, run simulations for each sample size via the
  `future` ecosystem.

- cores:

  Optional positive integer number of cores to use when
  `parallel = TRUE`. If `NULL`, uses one fewer than the number of
  available cores, with a minimum of one.

- seed:

  Optional integer seed for reproducibility.

## Value

An `anovapowersim_curve` object with `n_needed` and `total_n_needed`.

## Examples

    power_n(
      between = c(cond = 2),
      within = c(stim = 4),
      term = "cond:stim",
      target_pes = 0.14,
      alpha = 0.05,
      power = 0.80,
      n_sims = 1000, # use 5000+ for a more precise estimate
      seed = 123 # for reproducibility
    )
