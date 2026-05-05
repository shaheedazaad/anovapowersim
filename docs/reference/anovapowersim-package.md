# anovapowersim: Design-Based Power Simulation for Factorial ANOVA

Simulation-based power analysis for balanced factorial ANOVA designs.
Given between- and within-subject factor structures, a target partial
eta squared, and sample sizes, `anovapowersim` generates default
term-specific cell means, simulates datasets under sphericity, refits
the ANOVA with stats, and returns tidy power estimates and a `ggplot2`
power curve.

## Details

The primary entry point is
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md).
It accepts a design specification, a term name, a target partial eta
squared, and sample sizes. It returns an object of class
`anovapowersim_curve` that can be printed, summarised, or plotted with
[`plot_power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/plot_power_curve.md).

The lower-level building blocks are also exported so users can compose
custom simulations:

- [`balanced_anova_design()`](https://shaheedazaad.github.io/anovapowersim/reference/balanced_anova_design.md)

- [`design_term_means()`](https://shaheedazaad.github.io/anovapowersim/reference/design_term_means.md)

- [`simulate_design_dataset()`](https://shaheedazaad.github.io/anovapowersim/reference/simulate_design_dataset.md)

- [`compute_scale_factor()`](https://shaheedazaad.github.io/anovapowersim/reference/compute_scale_factor.md)

## See also

Useful links:

- <https://shaheedazaad.github.io/anovapowersim/>

- <https://github.com/shaheedazaad/anovapowersim>

- Report bugs at <https://github.com/shaheedazaad/anovapowersim/issues>

## Author

**Maintainer**: Shaheed Azaad <sazaad@uni-muenster.de>
