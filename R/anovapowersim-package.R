#' anovapowersim: Design-Based Power Simulation for Factorial ANOVA
#'
#' Simulation-based power analysis for balanced factorial ANOVA designs. Given
#' between- and within-subject factor structures, a target partial eta squared,
#' and sample sizes, `anovapowersim` generates term-specific cell means from an
#' explicit [means_pattern()] or a documented linear/Kronecker default,
#' simulates datasets using default or custom within-subject covariance
#' structures, refits the ANOVA, and returns tidy power estimates and a
#' `ggplot2` power curve.
#'
#' The main entry points are [power_n()] for required sample size,
#' [power_achieved()] for achieved power, [power_sensitivity()] for a minimum
#' detectable effect size, and [power_curve()] for power across explicit sample
#' sizes. [power_unbalanced()] provides experimental means-based simulation for
#' one fixed unbalanced design. Experimental calculation-only counterparts are
#' provided by [power_n_calc()], [power_achieved_calc()], and
#' [power_sensitivity_calc()].
#'
#' The lower-level building blocks are also exported so users can compose
#' custom simulations:
#'
#' * [balanced_anova_design()]
#' * [design_term_means()]
#' * [means_pattern()]
#' * [simulate_design_dataset()]
#' * [within_covariance()]
#' * [cell_design()]
#' * [unbalanced_covariance()]
#' * [compute_scale_factor()]
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data :=
#' @importFrom stats model.matrix qf reformulate setNames terms
## usethis namespace: end
NULL
