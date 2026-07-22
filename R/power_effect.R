#' Estimate achieved ANOVA power at a fixed sample size
#'
#' Simulates one balanced factorial ANOVA design point for a fixed partial eta
#' squared and sample size. The simulation estimate is the primary achieved
#' power result; noncentral-F calculated power is retained as a diagnostic.
#'
#' @section Lifecycle:
#' \ifelse{html}{\out{<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg" alt="[Experimental]"></a>}}{\strong{Experimental}}
#'
#' `power_achieved()` is experimental and is available only in the development
#' version of `anovapowersim`. Its API and reporting format may change.
#'
#' @inheritParams power_curve
#' @param n Number of subjects per between-subject cell. For a purely
#'   within-subject design, this is the total sample size.
#'
#' @return An `anovapowersim_achieved_power` object. `$results` contains the
#'   standard one-row power diagnostics. `$achieved_power` is the simulated
#'   power estimate and `$calculated_power` is the calculated-power diagnostic.
#'
#' @examples
#' \donttest{
#' power_achieved(
#'   between = c(group = 2),
#'   within = c(time = 2),
#'   term = "group:time",
#'   target_pes = 0.14,
#'   n = 20,
#'   n_sims = 100,
#'   seed = 123
#' )
#' }
#'
#' @export
power_achieved <- function(between = NULL,
                           within = NULL,
                           term,
                           target_pes,
                           n,
                           n_sims = 10000,
                           alpha = 0.05,
                           ss_type = "III",
                           gpower = FALSE,
                           progress = interactive(),
                           parallel = FALSE,
                           cores = NULL,
                           seed = NULL,
                           covariance = NULL,
                           means_pattern = NULL) {
  sd <- 1
  r <- 0.5
  setup <- prepare_power_curve_inputs(
    between = between,
    within = within,
    term = term,
    target_pes = target_pes,
    n_sims = n_sims,
    alpha = alpha,
    ss_type = ss_type,
    sd = sd,
    r = r,
    covariance = covariance,
    gpower = gpower,
    progress = progress,
    parallel = parallel,
    cores = cores,
    means_pattern = means_pattern
  )
  n <- validate_fixed_design_n(n, setup$spec)
  message_long_serial_run(setup$n_sims, setup$parallel)
  if (!is.null(seed)) set.seed(seed)

  progress_bar <- make_progress_bar(
    enabled = setup$progress,
    total = if (setup$parallel) 1L else setup$n_sims,
    label = "Simulating achieved power"
  )
  on.exit(close_progress_bar(progress_bar), add = TRUE)

  row <- run_design_power_at_n(
    spec = setup$spec,
    term = setup$term,
    target_pes = target_pes,
    n = n,
    n_sims = setup$n_sims,
    alpha = alpha,
    ss_type = setup$ss_type,
    sd = sd,
    r = r,
    covariance = setup$covariance,
    epsilon = setup$epsilon,
    gpower = setup$gpower,
    progress_bar = if (setup$parallel) NULL else progress_bar,
    parallel = setup$parallel,
    cores = setup$cores,
    resolved_means_pattern = setup$means_pattern
  )
  if (setup$parallel) tick_progress_bar(progress_bar)
  warn_power_disagreement(row, setup$n_sims)

  structure(
    list(
      results = row,
      term = setup$term,
      alpha = alpha,
      target_pes = target_pes,
      n = n,
      total_n = as.integer(n * max(1L, setup$spec$n_between_cells)),
      n_sims = setup$n_sims,
      achieved_power = row$power_sim[[1L]],
      calculated_power = row$power_calc[[1L]],
      gpower = setup$gpower,
      epsilon = setup$epsilon,
      covariance = setup$covariance,
      custom_covariance = setup$custom_covariance,
      custom_means_pattern = setup$custom_means_pattern,
      ss_type = setup$ss_type,
      design = setup$spec,
      call = match.call()
    ),
    class = "anovapowersim_achieved_power"
  )
}


#' Estimate ANOVA effect-size sensitivity at a fixed sample size
#'
#' Searches for the minimum detectable partial eta squared at a fixed sample
#' size and target power. Calculated power supplies an efficient initial
#' estimate. Explicit simulations then bracket the target and refine the
#' effect-size bracket using interpolation with midpoint fallback.
#'
#' @section Lifecycle:
#' \ifelse{html}{\out{<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg" alt="[Experimental]"></a>}}{\strong{Experimental}}
#'
#' `power_sensitivity()` is experimental and is available only in the
#' development version of `anovapowersim`. Its API and reporting format may
#' change.
#'
#' @inheritParams power_achieved
#' @param power Desired target power.
#' @param pes_min Lower bound of the partial eta-squared search interval.
#' @param pes_max Upper bound of the partial eta-squared search interval.
#' @param pes_tol Maximum width of the final simulated partial eta-squared
#'   bracket.
#'
#' @return An `anovapowersim_sensitivity` object. `$pes_needed` is the
#'   explicitly simulated upper effect-size bracket, or `NA` when `pes_max`
#'   does not achieve target power. `$results` contains every explicitly
#'   simulated effect size and its standard power diagnostics.
#'
#' @examples
#' \donttest{
#' power_sensitivity(
#'   between = c(group = 2),
#'   within = c(time = 2),
#'   term = "group:time",
#'   n = 20,
#'   power = 0.90,
#'   n_sims = 100,
#'   pes_tol = 0.01,
#'   seed = 123
#' )
#' }
#'
#' @export
power_sensitivity <- function(between = NULL,
                              within = NULL,
                              term,
                              n,
                              power = 0.90,
                              n_sims = 10000,
                              alpha = 0.05,
                              ss_type = "III",
                              pes_min = 1e-6,
                              pes_max = 0.99,
                              pes_tol = 0.001,
                              gpower = FALSE,
                              progress = interactive(),
                              parallel = FALSE,
                              cores = NULL,
                              seed = NULL,
                              covariance = NULL,
                              means_pattern = NULL) {
  sd <- 1
  r <- 0.5
  setup <- prepare_balanced_power_inputs(
    between = between,
    within = within,
    term = term,
    n_sims = n_sims,
    alpha = alpha,
    ss_type = ss_type,
    sd = sd,
    r = r,
    covariance = covariance,
    gpower = gpower,
    progress = progress,
    parallel = parallel,
    cores = cores,
    means_pattern = means_pattern
  )
  assert_unit_interval(power, "power")
  if (power < 0.90) {
    warning("Power greater than or equal to .90 is recommended.",
            call. = FALSE, immediate. = TRUE)
  }
  n <- validate_fixed_design_n(n, setup$spec)
  bounds <- validate_pes_search_controls(pes_min, pes_max, pes_tol)
  pes_min <- bounds$pes_min
  pes_max <- bounds$pes_max
  pes_tol <- bounds$pes_tol
  message_long_serial_run(setup$n_sims, setup$parallel)

  initial_pes <- estimate_calculated_pes_needed(
    spec = setup$spec,
    term = setup$term,
    n = n,
    target_power = power,
    alpha = alpha,
    ss_type = setup$ss_type,
    sd = sd,
    r = r,
    covariance = setup$covariance,
    epsilon = setup$epsilon,
    gpower = setup$gpower,
    pes_min = pes_min,
    pes_max = pes_max
  )
  if (!is.null(seed)) set.seed(seed)

  max_iter <- 25L
  progress_bar <- make_progress_bar(
    enabled = setup$progress,
    total = max_iter + 3L,
    label = "Searching effect size"
  )
  on.exit(close_progress_bar(progress_bar), add = TRUE)

  run_one <- function(pes) {
    row <- run_design_power_at_n(
      spec = setup$spec,
      term = setup$term,
      target_pes = pes,
      n = n,
      n_sims = setup$n_sims,
      alpha = alpha,
      ss_type = setup$ss_type,
      sd = sd,
      r = r,
      covariance = setup$covariance,
      epsilon = setup$epsilon,
      gpower = setup$gpower,
      progress_bar = NULL,
      parallel = setup$parallel,
      cores = setup$cores,
      resolved_means_pattern = setup$means_pattern
    )
    tibble::add_column(row, target_pes = as.numeric(pes), .before = 1L)
  }

  search <- adaptive_effect_search(
    run_one = run_one,
    target = power,
    pes_start = initial_pes,
    pes_min = pes_min,
    pes_max = pes_max,
    pes_tol = pes_tol,
    max_iter = max_iter,
    progress_bar = progress_bar
  )
  warn_power_disagreement(search$results, setup$n_sims)

  if (!search$reached) {
    warning(
      sprintf(
        paste(
          "Target power %.3f was not reached by `pes_max = %.6g`.",
          "Increase `pes_max` and rerun the search."
        ),
        power, pes_max
      ),
      call. = FALSE
    )
  } else if (!search$converged) {
    warning(
      sprintf(
        paste(
          "The sensitivity search reached its iteration cap before",
          "`pes_tol = %.6g`; reporting the best explicitly simulated",
          "upper bracket (pes = %.6g)."
        ),
        pes_tol, search$pes_upper
      ),
      call. = FALSE
    )
  }

  structure(
    list(
      results = search$results,
      term = setup$term,
      power = power,
      alpha = alpha,
      n = n,
      total_n = as.integer(n * max(1L, setup$spec$n_between_cells)),
      n_sims = setup$n_sims,
      pes_needed = if (search$reached) search$pes_upper else NA_real_,
      pes_lower = search$pes_lower,
      pes_upper = search$pes_upper,
      pes_min = pes_min,
      pes_max = pes_max,
      pes_tol = pes_tol,
      bracket_width = search$bracket_width,
      converged = search$converged,
      iterations = search$iterations,
      initial_pes = initial_pes,
      gpower = setup$gpower,
      epsilon = setup$epsilon,
      covariance = setup$covariance,
      custom_covariance = setup$custom_covariance,
      custom_means_pattern = setup$custom_means_pattern,
      ss_type = setup$ss_type,
      design = setup$spec,
      call = match.call()
    ),
    class = "anovapowersim_sensitivity"
  )
}


#' @keywords internal
#' @noRd
validate_fixed_design_n <- function(n, spec) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != as.integer(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  validate_calibration_n(n, spec, "n")
  n
}


#' @keywords internal
#' @noRd
validate_pes_search_controls <- function(pes_min, pes_max, pes_tol) {
  assert_unit_interval(pes_min, "pes_min")
  assert_unit_interval(pes_max, "pes_max")
  if (pes_min >= pes_max) {
    stop("`pes_min` must be smaller than `pes_max`.", call. = FALSE)
  }
  if (!is.numeric(pes_tol) || length(pes_tol) != 1L ||
      !is.finite(pes_tol) || pes_tol <= 0) {
    stop("`pes_tol` must be a single positive finite number.", call. = FALSE)
  }
  list(
    pes_min = as.numeric(pes_min),
    pes_max = as.numeric(pes_max),
    pes_tol = as.numeric(pes_tol)
  )
}


#' Use calculated power to choose a sensitivity-search starting point
#'
#' @keywords internal
#' @noRd
estimate_calculated_pes_needed <- function(spec, term, n, target_power, alpha,
                                           ss_type, sd, r, covariance,
                                           epsilon, gpower, pes_min, pes_max) {
  calculated_power <- function(pes) {
    power_calc_at_n(
      spec = spec,
      term = term,
      target_pes = pes,
      n = n,
      alpha = alpha,
      ss_type = ss_type,
      sd = sd,
      r = r,
      covariance = covariance,
      epsilon = epsilon,
      gpower = gpower
    )
  }
  lower_power <- calculated_power(pes_min)
  if (is.finite(lower_power) && lower_power >= target_power) return(pes_min)
  upper_power <- calculated_power(pes_max)
  if (!is.finite(upper_power) || upper_power <= target_power) return(pes_max)

  root <- tryCatch(
    stats::uniroot(
      function(pes) calculated_power(pes) - target_power,
      interval = c(pes_min, pes_max),
      tol = min(1e-8, (pes_max - pes_min) / 1000)
    )$root,
    error = function(e) NA_real_
  )
  if (is.finite(root)) min(pes_max, max(pes_min, root)) else {
    (pes_min + pes_max) / 2
  }
}


#' Continuously search an effect-size interval using explicit simulations
#'
#' @keywords internal
#' @noRd
adaptive_effect_search <- function(run_one, target, pes_start, pes_min,
                                   pes_max, pes_tol, max_iter = 25L,
                                   progress_bar = NULL) {
  visited <- list()
  simulate <- function(pes) {
    row <- run_one(as.numeric(pes))
    tick_progress_bar(progress_bar)
    visited[[length(visited) + 1L]] <<- row
    row
  }

  pes_start <- min(pes_max, max(pes_min, as.numeric(pes_start)))
  start_row <- simulate(pes_start)
  start_power <- start_row$power_sim[[1L]]
  lo <- NA_real_
  lo_power <- NA_real_
  hi <- NA_real_
  hi_power <- NA_real_

  if (is.finite(start_power) && start_power >= target) {
    hi <- pes_start
    hi_power <- start_power
    if (pes_start > pes_min) {
      lower_row <- simulate(pes_min)
      lower_power <- lower_row$power_sim[[1L]]
      if (is.finite(lower_power) && lower_power >= target) {
        hi <- pes_min
        hi_power <- lower_power
      } else {
        lo <- pes_min
        lo_power <- lower_power
      }
    }
  } else {
    lo <- pes_start
    lo_power <- start_power
    if (pes_start < pes_max) {
      upper_row <- simulate(pes_max)
      upper_power <- upper_row$power_sim[[1L]]
      if (is.finite(upper_power) && upper_power >= target) {
        hi <- pes_max
        hi_power <- upper_power
      } else {
        lo <- pes_max
        lo_power <- upper_power
      }
    }
  }

  reached <- is.finite(hi)
  iter <- 0L
  if (reached && is.finite(lo)) {
    while ((hi - lo) > pes_tol && iter < max_iter) {
      next_pes <- next_adaptive_pes(
        lo_pes = lo,
        lo_power = lo_power,
        hi_pes = hi,
        hi_power = hi_power,
        target = target
      )
      row <- simulate(next_pes)
      observed_power <- row$power_sim[[1L]]
      if (is.finite(observed_power) && observed_power >= target) {
        hi <- next_pes
        hi_power <- observed_power
      } else {
        lo <- next_pes
        lo_power <- observed_power
      }
      iter <- iter + 1L
    }
  }

  results <- dplyr::bind_rows(visited) |>
    dplyr::arrange(.data$target_pes) |>
    dplyr::distinct(.data$target_pes, .keep_all = TRUE)
  width <- if (reached && is.finite(lo)) hi - lo else if (reached) 0 else NA_real_
  list(
    results = results,
    reached = reached,
    converged = reached && (!is.finite(lo) || width <= pes_tol),
    pes_lower = if (is.finite(lo)) lo else NA_real_,
    pes_upper = if (reached) hi else NA_real_,
    bracket_width = width,
    iterations = iter
  )
}


#' @keywords internal
#' @noRd
next_adaptive_pes <- function(lo_pes, lo_power, hi_pes, hi_power, target) {
  midpoint <- lo_pes + (hi_pes - lo_pes) / 2
  candidate <- NA_real_
  if (is.finite(lo_power) && is.finite(hi_power) && hi_power > lo_power) {
    fraction <- (target - lo_power) / (hi_power - lo_power)
    candidate <- lo_pes + fraction * (hi_pes - lo_pes)
  }
  if (!is.finite(candidate) || candidate <= lo_pes || candidate >= hi_pes ||
      candidate == lo_pes || candidate == hi_pes) {
    candidate <- midpoint
  }
  as.numeric(candidate)
}
