# anovapowersim: Design-Based Power Simulation for Factorial ANOVA

Simulation-based power analysis for balanced factorial ANOVA designs.
Given between- and within-subject factor structures, a target partial
eta squared, and sample sizes, `anovapowersim` generates default
term-specific cell means, simulates datasets using default or custom
within-subject covariance structures, refits the ANOVA with stats, and
returns tidy power estimates and a `ggplot2` power curve.

## Details

The main entry points are
[`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
for required sample size,
[`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md)
for achieved power,
[`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
for a minimum detectable effect size, and
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md)
for power across explicit sample sizes.
[`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
provides experimental means-based simulation for one fixed unbalanced
design. Experimental calculation-only counterparts are provided by
[`power_n_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n_calc.md),
[`power_achieved_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved_calc.md),
and
[`power_sensitivity_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity_calc.md).

The lower-level building blocks are also exported so users can compose
custom simulations:

- [`balanced_anova_design()`](https://shaheedazaad.github.io/anovapowersim/reference/balanced_anova_design.md)

- [`design_term_means()`](https://shaheedazaad.github.io/anovapowersim/reference/design_term_means.md)

- [`simulate_design_dataset()`](https://shaheedazaad.github.io/anovapowersim/reference/simulate_design_dataset.md)

- [`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)

- [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)

- [`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md)

- [`compute_scale_factor()`](https://shaheedazaad.github.io/anovapowersim/reference/compute_scale_factor.md)

## See also

Useful links:

- <https://shaheedazaad.github.io/anovapowersim/>

- <https://github.com/shaheedazaad/anovapowersim>

- Report bugs at <https://github.com/shaheedazaad/anovapowersim/issues>

## Author

**Maintainer**: Shaheed Azaad <sazaad@uni-muenster.de>
