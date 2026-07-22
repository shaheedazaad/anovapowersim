# Package index

## Core simulation

Required-sample-size searches and power curves for balanced designs.

- [`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
  : Search for the sample size needed for target ANOVA power
- [`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md)
  : Simulate ANOVA power from a balanced factorial design
- [`plot_power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/plot_power_curve.md)
  : Plot a simulation-based power curve

## Fixed-sample simulation

Achieved power and effect-size sensitivity at a fixed sample size.

- [`power_achieved()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved.md)
  : Estimate achieved ANOVA power at a fixed sample size
- [`power_sensitivity()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity.md)
  : Estimate ANOVA effect-size sensitivity at a fixed sample size

## Calculated power

Calculation-only analyses that skip simulations.

- [`power_n_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n_calc.md)
  : Calculate the sample size needed for target ANOVA power
- [`power_achieved_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_achieved_calc.md)
  : Calculate achieved ANOVA power at a fixed sample size
- [`power_sensitivity_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_sensitivity_calc.md)
  : Calculate ANOVA effect-size sensitivity at a fixed sample size

## Unbalanced designs

Fixed-design simulation with explicit sample sizes, means, and SDs.

- [`power_unbalanced()`](https://shaheedazaad.github.io/anovapowersim/reference/power_unbalanced.md)
  : Simulate power for a fixed unbalanced ANOVA design
- [`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)
  : Define cells for a means-based unbalanced ANOVA design
- [`unbalanced_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/unbalanced_covariance.md)
  : Specify covariance for a means-based unbalanced design

## Design components

Component functions for power analyses.

- [`balanced_anova_design()`](https://shaheedazaad.github.io/anovapowersim/reference/balanced_anova_design.md)
  : Create a balanced factorial ANOVA design specification
- [`within_covariance()`](https://shaheedazaad.github.io/anovapowersim/reference/within_covariance.md)
  : Specify a within-subject covariance structure
- [`means_pattern()`](https://shaheedazaad.github.io/anovapowersim/reference/means_pattern.md)
  : Define a sparse relative cell-mean pattern
- [`design_term_means()`](https://shaheedazaad.github.io/anovapowersim/reference/design_term_means.md)
  : Build calibrated means for a design term
- [`simulate_design_dataset()`](https://shaheedazaad.github.io/anovapowersim/reference/simulate_design_dataset.md)
  : Simulate data from a balanced ANOVA design
- [`compute_scale_factor()`](https://shaheedazaad.github.io/anovapowersim/reference/compute_scale_factor.md)
  : Compute the mean-deviation scaling factor from a change in partial
  eta squared

## Output methods

Printing and summarising power-analysis results.

- [`print(`*`<anovapowersim_curve>`*`)`](https://shaheedazaad.github.io/anovapowersim/reference/print.anovapowersim_curve.md)
  : Print an anovapowersim power curve
- [`summary(`*`<anovapowersim_curve>`*`)`](https://shaheedazaad.github.io/anovapowersim/reference/summary.anovapowersim_curve.md)
  : Summarise an anovapowersim power curve
- [`print(`*`<anovapowersim_achieved_power>`*`)`](https://shaheedazaad.github.io/anovapowersim/reference/print.anovapowersim_achieved_power.md)
  : Print a fixed-sample achieved-power result
- [`summary(`*`<anovapowersim_achieved_power>`*`)`](https://shaheedazaad.github.io/anovapowersim/reference/summary.anovapowersim_achieved_power.md)
  : Summarise a fixed-sample achieved-power result
- [`print(`*`<anovapowersim_sensitivity>`*`)`](https://shaheedazaad.github.io/anovapowersim/reference/print.anovapowersim_sensitivity.md)
  : Print a fixed-sample sensitivity result
- [`summary(`*`<anovapowersim_sensitivity>`*`)`](https://shaheedazaad.github.io/anovapowersim/reference/summary.anovapowersim_sensitivity.md)
  : Summarise a fixed-sample sensitivity result
- [`print(`*`<anovapowersim_unbalanced_power>`*`)`](https://shaheedazaad.github.io/anovapowersim/reference/print.anovapowersim_unbalanced_power.md)
  : Print simulated power for an unbalanced design
- [`summary(`*`<anovapowersim_unbalanced_power>`*`)`](https://shaheedazaad.github.io/anovapowersim/reference/summary.anovapowersim_unbalanced_power.md)
  : Summarise simulated power for an unbalanced design
