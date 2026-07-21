# Getting started with anovapowersim

`anovapowersim` primarily simulates power for balanced factorial ANOVA
designs. Specify the factors/levels, the term of interest, and a target
partial eta squared. `anovapowersim` generates default term-specific
cell means, simulates datasets, refits the ANOVA, and estimates power.
An experimental interface for one fixed unbalanced design is described
below.

``` r

library(anovapowersim)
```

## Search for the required sample size

The easiest way to get your required sample size is to use
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
to search for the sample size needed to reach the requested `power`.

This example is a 2 x 2 mixed design with one between-subjects factor
(`cond`) and one within-subject factor (`stim`).

We specify that we are interested in the `cond:stim` interaction, and
that we want to have 90% power to detect a partial eta squared of 0.14.
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
will search for the required sample size per between-subject cell, so
`n = 17` gives total `N = 34`.

``` r

power_n(
  between = c(cond = 2), # cond has 2 levels
  within = c(stim = 4), # stim has 4 levels
  term = "cond:stim",
  target_pes = 0.14,
  alpha = 0.05,
  power = 0.90,
  n_sims = 1000, # use 5000+ for a more precise estimate
  seed = 123 # for reproducibility
)
```

    #> <anovapowersim_curve>
    #>   term:          'cond:stim'
    #>   target power:  0.900
    #>   alpha:         0.05
    #>   effect size:   pes = 0.14
    #>   n values:      6 per-cell sample sizes visited
    #>   sims per cell size: 1000
    #>   SS type:       III
    #>   n needed for between-subjects cell: 17
    #>   total N needed: 34
    #> 
    #>  n_per_cell total_n n_sims num_df den_df    ncp power_calc power_sim
    #>          13      26   1000      3     72 11.721      0.808     0.795
    #>          16      32   1000      3     90 14.651      0.897     0.885
    #>          17      34   1000      3     96 15.628      0.917     0.912
    #>          18      36   1000      3    102 16.605      0.934     0.936
    #>          20      40   1000      3    114 18.558      0.958     0.946
    #>          26      52   1000      3    150 24.419      0.991     0.991

Note: here we use 1000 simulations for a quick example, but the package
defaults to 10000 simulations for more precise estimates.

The output table uses compact column names: `n_per_cell` is the sample
size per between-subject cell, `total_n` is the full sample size,
`num_df` and `den_df` are the ANOVA degrees of freedom, `ncp` is the
noncentrality parameter, `power_calc` is the noncentral F power
calculation, and `power_sim` is the simulation estimate.

## Fixed-sample achieved power and sensitivity

[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md)
and
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
are experimental and available only in the development version of
`anovapowersim`. Their APIs and reporting formats may change.

When sample size is already fixed,
[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md)
estimates simulated power for a specified partial eta squared. The `n`
argument is subjects per between-subject cell, or total N for a purely
within-subject design.

``` r

power_achieved(
  between = c(cond = 2),
  within = c(stim = 2),
  term = "cond:stim",
  target_pes = 0.14,
  n = 20,
  n_sims = 1000,
  seed = 123
)
```

Use
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
when the effect size is unknown. It searches for the smallest explicitly
simulated partial eta squared that reaches the requested power at the
fixed sample size.

``` r

power_sensitivity(
  between = c(cond = 2),
  within = c(stim = 2),
  term = "cond:stim",
  n = 20,
  power = 0.90,
  n_sims = 1000,
  pes_tol = 0.001,
  seed = 123
)
```

When you want to skip simulations, use
[`power_achieved_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved_calc.md)
and
[`power_sensitivity_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity_calc.md).
These experimental functions use calculated noncentral-F power and
accept a planned nonsphericity correction through `epsilon`.

``` r

power_achieved_calc(
  between = c(cond = 2),
  within = c(stim = 3),
  term = "cond:stim",
  target_pes = 0.14,
  n = 20,
  gpower = TRUE,
  epsilon = 0.80
)

power_sensitivity_calc(
  between = c(cond = 2),
  within = c(stim = 3),
  term = "cond:stim",
  n = 20,
  power = 0.90,
  gpower = TRUE,
  epsilon = 0.80
)
```

## Power for unbalanced designs

[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
is an experimental, development-version-only function for the less
common case where the exact allocation, population means, and standard
deviations are already known. It simulates that one fixed design; it
does not search over sample sizes and does not return calculated power.

Define every cell in the complete factorial design with
[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md).
The factor values identify the cell, `n` is the number of subjects in
its between-subject group, `m` is its population mean, and `sd` is its
marginal population standard deviation. For repeated-measures designs,
the same `n` must therefore appear on every within-subject row belonging
to a given between-subject group.

``` r

unbalanced_design <- cell_design(
  group = "control",   time = "pre",  n = 22, m = 10.0, sd = 2.0,
  group = "control",   time = "post", n = 22, m = 11.0, sd = 2.2,
  group = "treatment", time = "pre",  n = 31, m = 10.1, sd = 2.4,
  group = "treatment", time = "post", n = 31, m = 12.4, sd = 2.8
)
```

Here, the means define the complete effect pattern, including any main
effects and interactions. The standard deviations can differ across
groups and within-subject cells. Name the repeated-measures factors in
`within`; all other factors are treated as between-subject factors.

For a within-subject design, use
[`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md)
to define the correlations among repeated measurements. It deliberately
has no standard deviation argument:
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
combines the correlations with the cell-specific `sd` values above to
construct the covariance matrix separately for each between-subject
group. With one within-subject factor, correlation pair names use its
level names, such as `"pre:post"`. With multiple within- subject
factors, join the level combinations with underscores, as described in
the nonsphericity section below.

``` r

power_unbalanced(
  design = unbalanced_design,
  within = "time",
  term = "group:time",
  covariance = unbalanced_covariance(
    default_correlation = 0.5,
    correlations = c("pre:post" = 0.7)
  ),
  n_sims = 5000,
  parallel = TRUE,
  seed = 123
)
```

The result reports simulated power, the partial eta squared from a
deterministic reference dataset matching the supplied design
assumptions, and the mean, median, and 95% interval of partial eta
squared across successful simulations. Because power in an unbalanced
design depends on the allocation, means, variances, covariance
structure, tested term, and sums-of-squares type, these effect-size
summaries describe this exact design rather than a separate target
effect size.

### Adding factors and levels

You can add factors and levels as needed, and specify any term of
interest. For, example if we want to add a between condition with 3
levels, and we are interested in the 3-way interaction, we can do:

``` r

power_n(
  between = c(cond = 2, age = 3), # cond has 2 levels, age has 3 levels
  within = c(stim = 4), # stim has 4 levels
  term = "cond:stim:age",
  target_pes = 0.14,
  alpha = 0.05,
  power = 0.90,
  n_sims = 1000, # use 5000+ for a more precise estimate
  seed = 123 # for reproducibility
)
```

## Simulate a power curve

You might want to see how power changes across a range of sample sizes.
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md)
simulates power across a range of sample sizes, which you can specify
with `n_range`. The result is a tidy data frame that you can plot with
[`plot_power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/plot_power_curve.md).

``` r

pc <- power_curve(
  between = c(cond = 2),
  within = c(stim = 2),
  term = "cond:stim",
  target_pes = 0.14,
  n_range = c(16, 20, 23, 28), # n per between-subject cell
  n_sims = 1000,
  seed = 123
)
pc
```

    #> <anovapowersim_curve>
    #>   term:          'cond:stim'
    #>   target power:  <not specified>
    #>   alpha:         0.05
    #>   effect size:   pes = 0.14
    #>   n values:      4 per-cell sample sizes visited
    #>   sims per cell size: 1000
    #>   SS type:       III
    #>   n needed for between-subjects cell: <not reached>
    #>   total N needed: <not reached>
    #> 
    #>  n_per_cell total_n n_sims num_df den_df   ncp power_calc power_sim
    #>          16      32   1000      1     30 4.884      0.571     0.557
    #>          20      40   1000      1     38 6.186      0.679     0.666
    #>          23      46   1000      1     44 7.163      0.745     0.735
    #>          28      56   1000      1     54 8.791      0.829     0.831

``` r

plot_power_curve(
  pc,
  power_lines = c(.80, .90) # adds horizontal lines at 80% and 90% power
)
```

![](anovapowersim_files/figure-html/plot-fixed-1.png)

## Nonspherical designs

By default, the balanced simulation functions
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md),
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md),
[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md),
and
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
use a standard deviation of 1 for every within-subject cell and a
correlation of 0.5 between every pair. All four accept a
[`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
object when repeated measurements have different standard deviations or
correlations. The simulation-free `_calc()` functions instead accept a
planned `epsilon` value.

[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
is the exception: its standard deviations are supplied in
[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md),
so its `covariance` argument accepts an
[`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md)
object containing correlations only.

The within-subject cell names combine factor levels with underscores.
For `within = c(time = 3, condition = 2)`, the cells are
`time1_condition1`, `time1_condition2`, `time2_condition1`,
`time2_condition2`, `time3_condition1`, and `time3_condition2`.
Correlation names join a pair of cells with `:`.

Define the covariance assumptions separately so they remain easy to read
and check:

``` r

covariance <- within_covariance(
  default_sd = 1,
  default_correlation = 0.5,
  standard_deviations = c(
    "time3_condition2" = 1.2
  ),
  correlations = c(
    "time1_condition1:time1_condition2" = 0.6,
    "time2_condition1:time2_condition2" = 0.7
  )
)
```

Unlisted measurements use `default_sd`, and unlisted pairs use
`default_correlation`. Feed the resulting specification into
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md):

``` r

power_n(
  within = c(time = 3, condition = 2),
  term = "time:condition",
  target_pes = 0.14,
  power = 0.90,
  n_sims = 1000,
  covariance = covariance,
  seed = 123
)
```

The custom covariance structure is used to calibrate the target effect
and to generate every simulated dataset. For the tested term, the
package also projects the covariance matrix into the term’s
within-subject contrast space and calculates a population
Greenhouse–Geisser epsilon. `power_calc` applies that epsilon to its
numerator degrees of freedom, denominator degrees of freedom, and
noncentrality parameter. This is equivalent to supplying the derived
value explicitly as `power_n_calc(epsilon = ...)`.

`power_sim` reflects the custom data-generating covariance, but its
simulated ANOVA currently uses the ordinary univariate p-value rather
than a sample-estimated Greenhouse–Geisser or Huynh–Feldt correction.

## Advanced options

### Run simulations in parallel

For larger simulation runs, set `parallel = TRUE`. If you do not set
`cores`, `anovapowersim` uses one fewer than the number of available
cores and prints a message with the chosen count. Set `cores` explicitly
when you want a fixed number of cores.

``` r

power_curve(
  between = c(cond = 2),
  within = c(stim = 2),
  term = "cond:stim",
  target_pes = 0.14,
  n_range = c(16, 20, 23, 28),
  n_sims = 5000,
  parallel = TRUE,
  cores = 4,
  seed = 123
)
```

### Match the G\*Power convention

By default, `anovapowersim` calibrates the simulated cell means so the
empirical reference dataset has the requested partial eta squared under
the fitted [`stats::aov()`](https://rdrr.io/r/stats/aov.html) model.
This corresponds to the fitted ANOVA denominator-df noncentrality
convention.

Set `gpower = TRUE` when you want the G\*Power-style convention (when
using the ’as in Cohen (1988) option for within-subjects designs)
`lambda = total_n * f^2`.

``` r

power_n(
  between = c(cond = 2),
  within = c(stim = 4),
  term = "cond:stim",
  target_pes = 0.14,
  alpha = 0.05,
  power = 0.90,
  n_sims = 1000,
  seed = 123,
  gpower = TRUE
)
```

    #> <anovapowersim_curve>
    #>   term:          'cond:stim'
    #>   target power:  0.900
    #>   alpha:         0.05
    #>   effect size:   pes = 0.14
    #>   n values:      6 per-cell sample sizes visited
    #>   sims per cell size: 1000
    #>   G*Power convention: TRUE
    #>   SS type:       III
    #>   n needed for between-subjects cell: 45
    #>   total N needed: 90
    #> 
    #>  n_per_cell total_n n_sims num_df den_df    ncp power_calc power_sim
    #>          36      72   1000      3    210 11.721      0.823     0.829
    #>          44      88   1000      3    258 14.326      0.899     0.889
    #>          45      90   1000      3    264 14.651      0.906     0.907
    #>          47      94   1000      3    276 15.302      0.919     0.919
    #>          52     104   1000      3    306 16.930      0.944     0.936
    #>          72     144   1000      3    426 23.442      0.989     0.990
