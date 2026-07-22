#' Print an anovapowersim power curve
#'
#' Compact one-screen summary: target, term, effective effect size,
#' estimated per-cell and total sample sizes, and the first and last rows of
#' the power curve.
#'
#' @param x An `anovapowersim_curve` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.anovapowersim_curve <- function(x, ...) {
  calc_only <- is_calculation_only_curve(x)
  cat("<anovapowersim_curve>\n")
  cat("  term:          '", x$term, "'\n", sep = "")
  cat("  target power:  ",
      if (is.null(x$power) || !is.finite(x$power)) "<not specified>"
      else sprintf("%.3f", x$power),
      "\n", sep = "")
  cat("  alpha:         ", format(x$alpha), "\n", sep = "")
  cat("  effect size:   pes = ", format(round(x$target_pes, 4)), sep = "")
  if (!is.null(x$scale_factor) && is.finite(x$scale_factor) &&
      !isTRUE(all.equal(x$scale_factor, 1))) {
    cat("  (rescaled, k = ", format(round(x$scale_factor, 3)), ")", sep = "")
  }
  cat("\n")
  if (inherits(x$design, "anovapowersim_unbalanced_design_spec")) {
    cat("  n values:      explicit unbalanced cell counts\n", sep = "")
  } else {
    cat("  n values:      ", nrow(x$results), " per-cell sample sizes visited\n", sep = "")
  }
  if (calc_only) {
    cat("  calculation:   calculated power only\n", sep = "")
  } else {
    cat("  sims per cell size: ", x$n_sims, "\n", sep = "")
  }
  if (isTRUE(x$gpower)) {
    cat("  G*Power convention: TRUE\n", sep = "")
  }
  if (!is.null(x$epsilon) && is.finite(x$epsilon) && x$epsilon < 1) {
    cat("  epsilon:       ", format(x$epsilon), "\n", sep = "")
  }
  if (isTRUE(x$custom_covariance)) {
    cat(
      "  covariance:    custom ", nrow(x$covariance), " x ",
      ncol(x$covariance), " within-subject matrix\n",
      sep = ""
    )
  }
  if (!is.null(x$ss_type)) {
    cat("  SS type:       ", x$ss_type, "\n", sep = "")
  }

  if (inherits(x$design, "anovapowersim_unbalanced_design_spec")) {
    cat("  design:        unbalanced between-subject cells\n", sep = "")
  } else {
    cat("  n needed for between-subjects cell: ",
        if (is.na(x$n_needed)) "<not reached>" else x$n_needed, "\n",
        sep = "")
    cat("  total N needed: ",
        if (is.na(x$total_n_needed)) "<not reached>" else x$total_n_needed, "\n",
        sep = "")
  }
  cat("\n")
  print(format_power_results(x$results), row.names = FALSE)
  invisible(x)
}


#' Summarise an anovapowersim power curve
#'
#' Returns the full `$results` tibble along with a small header containing
#' the target, effective effect size, and estimated `n_needed`.
#'
#' @param object An `anovapowersim_curve` object.
#' @param ... Unused.
#'
#' @return A list with elements `header` (named character) and `curve`
#'   (tibble), invisibly; printed to console as well.
#' @export
summary.anovapowersim_curve <- function(object, ...) {
  calc_only <- is_calculation_only_curve(object)
  header <- c(
    term         = object$term,
    target_power = if (is.null(object$power) || !is.finite(object$power)) {
      "<not specified>"
    } else sprintf("%.3f", object$power),
    alpha        = sprintf("%.3f", object$alpha),
    target_pes   = sprintf("%.4f", object$target_pes),
    scale_factor = if (is.null(object$scale_factor) ||
                       !is.finite(object$scale_factor)) {
      "<not recorded>"
    } else sprintf("%.3f", object$scale_factor),
    calculation  = if (calc_only) "calculated power only" else "simulation",
    n_sims       = if (calc_only) "<not run>" else as.character(object$n_sims),
    n_needed_between_subjects_cell =
      if (is.na(object$n_needed)) "<not reached>"
      else as.character(object$n_needed),
    total_n_needed = if (is.na(object$total_n_needed)) "<not reached>"
                     else as.character(object$total_n_needed)
  )
  if (isTRUE(object$gpower)) {
    header <- append(header, c(gpower_convention = "TRUE"), after = 4L)
  }
  if (!is.null(object$epsilon) && is.finite(object$epsilon) &&
      object$epsilon < 1) {
    header <- append(
      header,
      c(epsilon = sprintf("%.4f", object$epsilon)),
      after = 4L
    )
  }
  if (isTRUE(object$custom_covariance)) {
    header <- append(
      header,
      c(covariance = paste0(
        "custom ", nrow(object$covariance), " x ", ncol(object$covariance),
        " within-subject matrix"
      )),
      after = 4L
    )
  }
  if (!is.null(object$ss_type)) {
    header <- append(header, c(ss_type = object$ss_type), after = 4L)
  }
  out <- list(header = header, curve = object$results)
  cat(if (calc_only) {
    "anovapowersim calculated-power summary\n"
  } else {
    "anovapowersim power simulation summary\n"
  })
  cat("----------------------------------\n")
  for (nm in names(header)) {
    cat(sprintf("  %-13s %s\n", paste0(nm, ":"), header[[nm]]))
  }
  cat("\nPower curve:\n")
  print(format_power_results(object$results), row.names = FALSE)
  invisible(out)
}


#' Print a fixed-sample achieved-power result
#'
#' @param x An `anovapowersim_achieved_power` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.anovapowersim_achieved_power <- function(x, ...) {
  calc_only <- isTRUE(x$calculation_only)
  cat("<anovapowersim_achieved_power>\n")
  cat("  term:             '", x$term, "'\n", sep = "")
  if (!length(x$design$between)) {
    cat("  fixed total N:    ", x$total_n, "\n", sep = "")
  } else {
    cat("  fixed n per cell: ", x$n, "\n", sep = "")
    cat("  fixed total N:    ", x$total_n, "\n", sep = "")
  }
  cat("  target pes:       ", sprintf("%.4f", x$target_pes), "\n", sep = "")
  cat("  alpha:            ", format(x$alpha), "\n", sep = "")
  if (calc_only) {
    cat("  calculation:      calculated power only\n")
  } else {
    cat("  simulations:      ", x$n_sims, "\n", sep = "")
  }
  cat("  achieved power:   ", sprintf("%.3f", x$achieved_power),
      if (calc_only) " (calculated)" else "", "\n", sep = "")
  cat("  calculated power: ", sprintf("%.3f", x$calculated_power), "\n", sep = "")
  if (isTRUE(x$gpower)) cat("  G*Power convention: TRUE\n")
  if (!is.null(x$epsilon) && is.finite(x$epsilon) && x$epsilon < 1) {
    cat("  epsilon:          ", format(x$epsilon), "\n", sep = "")
  }
  if (isTRUE(x$custom_covariance)) {
    cat("  covariance:       custom ", nrow(x$covariance), " x ",
        ncol(x$covariance), " within-subject matrix\n", sep = "")
  }
  if (!is.null(x$ss_type)) {
    cat("  SS type:          ", x$ss_type, "\n", sep = "")
  }
  cat("\n")
  print(format_power_results(x$results), row.names = FALSE)
  invisible(x)
}


#' Summarise a fixed-sample achieved-power result
#'
#' @param object An `anovapowersim_achieved_power` object.
#' @param ... Unused.
#'
#' @return A list with `header` and `results`, invisibly; printed to the
#'   console as well.
#' @export
summary.anovapowersim_achieved_power <- function(object, ...) {
  calc_only <- isTRUE(object$calculation_only)
  fixed_n <- if (!length(object$design$between)) {
    c(fixed_total_n = as.character(object$total_n))
  } else {
    c(
      fixed_n_per_cell = as.character(object$n),
      fixed_total_n = as.character(object$total_n)
    )
  }
  header <- c(
    term = object$term,
    fixed_n,
    target_pes = sprintf("%.4f", object$target_pes),
    alpha = sprintf("%.3f", object$alpha),
    calculation = if (calc_only) "calculated power only" else "simulation",
    n_sims = if (calc_only) "<not run>" else as.character(object$n_sims),
    achieved_power = sprintf("%.3f", object$achieved_power),
    calculated_power = sprintf("%.3f", object$calculated_power),
    ss_type = object$ss_type
  )
  if (isTRUE(object$gpower)) {
    header <- append(header, c(gpower_convention = "TRUE"), after = 5L)
  }
  if (!is.null(object$epsilon) && is.finite(object$epsilon) &&
      object$epsilon < 1) {
    header <- append(header, c(epsilon = sprintf("%.4f", object$epsilon)),
                     after = 5L)
  }
  if (isTRUE(object$custom_covariance)) {
    header <- append(
      header,
      c(covariance = paste0(
        "custom ", nrow(object$covariance), " x ", ncol(object$covariance),
        " within-subject matrix"
      )),
      after = 5L
    )
  }
  out <- list(header = header, results = object$results)
  cat(if (calc_only) {
    "anovapowersim calculated achieved-power summary\n"
  } else {
    "anovapowersim achieved-power summary\n"
  })
  cat("-------------------------------------\n")
  for (nm in names(header)) {
    cat(sprintf("  %-18s %s\n", paste0(nm, ":"), header[[nm]]))
  }
  cat("\nPower diagnostics:\n")
  print(format_power_results(object$results), row.names = FALSE)
  invisible(out)
}


#' Print a fixed-sample sensitivity result
#'
#' @param x An `anovapowersim_sensitivity` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.anovapowersim_sensitivity <- function(x, ...) {
  calc_only <- isTRUE(x$calculation_only)
  cat("<anovapowersim_sensitivity>\n")
  cat("  term:             '", x$term, "'\n", sep = "")
  if (!length(x$design$between)) {
    cat("  fixed total N:    ", x$total_n, "\n", sep = "")
  } else {
    cat("  fixed n per cell: ", x$n, "\n", sep = "")
    cat("  fixed total N:    ", x$total_n, "\n", sep = "")
  }
  cat("  target power:     ", sprintf("%.3f", x$power), "\n", sep = "")
  cat("  alpha:            ", format(x$alpha), "\n", sep = "")
  cat(
    "  detectable pes:   ",
    if (is.na(x$pes_needed)) "<not reached>" else sprintf("%.6f", x$pes_needed),
    "\n", sep = ""
  )
  cat("  search interval:  [", format(x$pes_min), ", ",
      format(x$pes_max), "]\n", sep = "")
  cat("  requested width:  ", format(x$pes_tol), "\n", sep = "")
  cat(if (calc_only) "  calculated points: " else "  simulated points: ",
      nrow(x$results), "\n", sep = "")
  cat("  final width:      ",
      if (is.finite(x$bracket_width)) format(x$bracket_width) else "<none>",
      "\n", sep = "")
  cat("  converged:        ", if (isTRUE(x$converged)) "yes" else "no", "\n",
      sep = "")
  if (calc_only) {
    cat("  calculation:      calculated power only\n")
  } else {
    cat("  simulations/point:", x$n_sims, "\n", sep = " ")
  }
  if (isTRUE(x$gpower)) cat("  G*Power convention: TRUE\n")
  if (!is.null(x$epsilon) && is.finite(x$epsilon) && x$epsilon < 1) {
    cat("  epsilon:          ", format(x$epsilon), "\n", sep = "")
  }
  if (isTRUE(x$custom_covariance)) {
    cat("  covariance:       custom ", nrow(x$covariance), " x ",
        ncol(x$covariance), " within-subject matrix\n", sep = "")
  }
  if (!is.null(x$ss_type)) {
    cat("  SS type:          ", x$ss_type, "\n", sep = "")
  }
  cat("\n")
  print(format_power_results(x$results), row.names = FALSE)
  invisible(x)
}


#' Summarise a fixed-sample sensitivity result
#'
#' @param object An `anovapowersim_sensitivity` object.
#' @param ... Unused.
#'
#' @return A list with `header` and `results`, invisibly; printed to the
#'   console as well.
#' @export
summary.anovapowersim_sensitivity <- function(object, ...) {
  calc_only <- isTRUE(object$calculation_only)
  fixed_n <- if (!length(object$design$between)) {
    c(fixed_total_n = as.character(object$total_n))
  } else {
    c(
      fixed_n_per_cell = as.character(object$n),
      fixed_total_n = as.character(object$total_n)
    )
  }
  header <- c(
    term = object$term,
    fixed_n,
    target_power = sprintf("%.3f", object$power),
    alpha = sprintf("%.3f", object$alpha),
    detectable_pes = if (is.na(object$pes_needed)) {
      "<not reached>"
    } else sprintf("%.6f", object$pes_needed),
    search_interval = sprintf("[%.6g, %.6g]", object$pes_min, object$pes_max),
    requested_bracket_width = format(object$pes_tol),
    final_bracket_width = if (is.finite(object$bracket_width)) {
      format(object$bracket_width)
    } else "<none>",
    evaluated_points = as.character(nrow(object$results)),
    converged = if (isTRUE(object$converged)) "yes" else "no",
    calculation = if (calc_only) "calculated power only" else "simulation",
    n_sims_per_point = if (calc_only) {
      "<not run>"
    } else as.character(object$n_sims),
    ss_type = object$ss_type
  )
  if (isTRUE(object$gpower)) {
    header <- append(header, c(gpower_convention = "TRUE"), after = 5L)
  }
  if (!is.null(object$epsilon) && is.finite(object$epsilon) &&
      object$epsilon < 1) {
    header <- append(header, c(epsilon = sprintf("%.4f", object$epsilon)),
                     after = 5L)
  }
  if (isTRUE(object$custom_covariance)) {
    header <- append(
      header,
      c(covariance = paste0(
        "custom ", nrow(object$covariance), " x ", ncol(object$covariance),
        " within-subject matrix"
      )),
      after = 5L
    )
  }
  out <- list(header = header, results = object$results)
  cat(if (calc_only) {
    "anovapowersim calculated-power sensitivity summary\n"
  } else {
    "anovapowersim sensitivity summary\n"
  })
  cat("---------------------------------\n")
  for (nm in names(header)) {
    cat(sprintf("  %-25s %s\n", paste0(nm, ":"), header[[nm]]))
  }
  cat(if (calc_only) "\nCalculated-power sensitivity search:\n" else
      "\nSensitivity search:\n")
  print(format_power_results(object$results), row.names = FALSE)
  invisible(out)
}


#' @keywords internal
#' @noRd
format_power_results <- function(x) {
  x <- as.data.frame(x)
  if ("target_pes" %in% names(x)) {
    x$target_pes <- ifelse(
      is.na(x$target_pes), NA_character_, sprintf("%.6f", x$target_pes)
    )
  }
  numeric_3dp <- intersect(
    c("ncp", "power_calc", "power_sim"),
    names(x)
  )
  for (nm in numeric_3dp) {
    x[[nm]] <- ifelse(is.na(x[[nm]]), NA_character_, sprintf("%.3f", x[[nm]]))
  }
  x
}


#' @keywords internal
#' @noRd
is_calculation_only_curve <- function(x) {
  is.null(x$n_sims) || length(x$n_sims) != 1L || is.na(x$n_sims)
}
