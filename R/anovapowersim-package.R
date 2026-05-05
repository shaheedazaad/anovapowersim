#' anovapowersim: Design-Based Power Simulation for Factorial ANOVA
#'
#' Simulation-based power analysis for balanced factorial ANOVA designs. Given
#' between- and within-subject factor structures, a target partial eta squared,
#' and sample sizes, `anovapowersim` generates default term-specific cell means,
#' simulates datasets under sphericity, refits the ANOVA with \pkg{stats}, and
#' returns tidy power estimates and a `ggplot2` power curve.
#'
#' The primary entry point is [power_curve()]. It accepts a design
#' specification, a term name, a target partial eta squared, and sample sizes.
#' It returns an object of class `anovapowersim_curve` that can be printed,
#' summarised, or plotted with [plot_power_curve()].
#'
#' The lower-level building blocks are also exported so users can compose
#' custom simulations:
#'
#' * [balanced_anova_design()]
#' * [design_term_means()]
#' * [simulate_design_dataset()]
#' * [compute_scale_factor()]
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data :=
#' @importFrom stats model.matrix qf reformulate setNames terms
## usethis namespace: end
NULL
