#' Simulate ANOVA power from a balanced factorial design
#'
#' Simulation-based power estimation for balanced factorial designs. Users
#' specify the between- and within-subject factors, the ANOVA term to test, a
#' target partial eta squared, and explicit sample sizes. The function projects
#' an explicit relative means pattern (or uses the documented
#' linear/Kronecker default), scales it to the requested partial eta squared,
#' simulates datasets, refits the ANOVA, and estimates power by counting
#' `p < alpha`.
#'
#' @param between Named integer vector of between-subject factor level counts,
#'   e.g. `c(group = 2)`. Use `NULL` for no between-subject factors.
#' @param within Named integer vector of within-subject factor level counts,
#'   e.g. `c(time = 3, condition = 4)`. Use `NULL` for no within-subject
#'   factors.
#' @param term Character scalar naming the ANOVA term to test, e.g.
#'   `"group:time"`. Interaction terms are order-insensitive; `"time:group"`
#'   resolves to `"group:time"` when that is the design's factor order.
#' @param target_pes Target partial eta squared for `term`.
#' @param n_range Integer vector of sample sizes per between-subject cell. For
#'   pure within-subject designs, this is the total sample size.
#' @param n_sims Number of simulated datasets per sample size.
#' @param alpha Significance threshold.
#' @param ss_type Sums-of-squares type for the tested ANOVA term. `"III"` is
#'   the default for order-invariant tests in unbalanced designs. Use `"I"` to
#'   reproduce sequential `stats::aov()` tests. Greenhouse--Geisser-corrected
#'   simulated p-values are available only for `"III"` and `"II"`.
#' @param sim_correction Sphericity correction for simulated p-values:
#'   `"auto"` (the default) uses Greenhouse--Geisser correction when the
#'   term-specific population epsilon is below `1 - 1e-8` and `ss_type` is
#'   `"II"` or `"III"`; `"GG"` requests correction for every simulated
#'   dataset; and `"none"` always uses the uncorrected univariate test.
#'   `"GG"` is an error with `ss_type = "I"`. For a between-only term or a
#'   within component with one degree of freedom, `"GG"` silently resolves to
#'   `"none"` because no sphericity correction applies.
#' @param gpower Logical; if `TRUE`, calibrate means to the G*Power-style
#'   noncentrality convention `lambda = total_n * f^2`. The default `FALSE`
#'   calibrates the empirical reference dataset to `target_pes`, equivalent to
#'   `lambda = den_df * f^2` for the fitted ANOVA. G*Power's estimates can
#'   differ from `target_pes`, especially for small samples or terms with more
#'   degrees of freedom; a warning is issued when `gpower = TRUE`. The default
#'   `gpower = FALSE` is recommended.
#' @param progress Logical; if `TRUE`, show a text progress bar.
#' @param parallel Logical; if `TRUE`, run simulations for each sample size via
#'   the `future` ecosystem.
#' @param cores Optional positive integer number of cores to use when
#'   `parallel = TRUE`. If `NULL`, uses one fewer than the number of available
#'   cores, with a minimum of one.
#' @param seed Optional integer seed for reproducibility.
#' @param covariance Optional within-subject covariance specification created
#'   by [within_covariance()]. Raw covariance matrices are not accepted, which
#'   avoids silently assuming a within-cell order. The default `NULL` uses
#'   standard deviations of `1` and a compound-symmetric correlation of `0.5`
#'   and issues a warning stating those defaults. A
#'   [within_covariance()] specification issues a warning when correlation
#'   pairs are omitted: its `default_correlation` applies only to those
#'   undefined pairs, while explicitly defined correlations are unchanged.
#'   All measurements use the specification's common marginal variance, while
#'   unequal correlations remain supported. For terms containing
#'   within-subject factors, the resolved covariance matrix is also used to
#'   derive a term-specific population Greenhouse--Geisser epsilon for
#'   `power_calc`. With the default `sim_correction = "auto"`, a population
#'   epsilon below `1 - 1e-8` also selects Greenhouse--Geisser-corrected
#'   simulated p-values for `ss_type` `"II"` or `"III"`.
#' @param means_pattern Optional relative cell-mean shape created by
#'   [means_pattern()]. The sparse values are projected onto `term`, normalized,
#'   and uniformly rescaled to reach `target_pes`. If `NULL`, simulations use
#'   the package's deterministic linear/Kronecker pattern. For multi-df
#'   nonspherical within-subject terms, simulated power is conditional on this
#'   direction, so an explicit pattern is recommended when the expected shape
#'   is known.
#'
#' @return An `anovapowersim_curve` object. The `$results` tibble contains
#'   `n_per_cell`, `total_n`, `n_sims`, successful and failed simulation counts
#'   (`valid_sims`, `failed_sims`), the population nonsphericity correction
#'   (`epsilon`), numerator and denominator degrees of freedom (`num_df`,
#'   `den_df`), the noncentrality parameter (`ncp`), calculated power
#'   (`power_calc`), and simulated power (`power_sim`). The full-precision
#'   `power_sim` value, not its printed three-decimal representation, is used
#'   by adaptive searches. When `epsilon < 1`, the reported degrees of freedom
#'   and noncentrality are the corrected values used for `power_calc`. With the
#'   default `sim_correction = "auto"`, `power_sim` uses the
#'   Greenhouse--Geisser-corrected simulated p-value when
#'   `epsilon < 1 - 1e-8` and `ss_type` is `"III"` or `"II"`; otherwise it
#'   uses the uncorrected univariate test. Balanced simulation result objects
#'   also include `custom_means_pattern`, indicating whether the relative
#'   direction was supplied explicitly, plus `sim_correction` and
#'   `sim_correction_resolved` for the requested and applied simulated-test
#'   correction.
#'
#' @section Simulated sphericity correction:
#' `sim_correction` changes only `power_sim`. When Greenhouse--Geisser
#' correction is selected, each simulated dataset is tested using its own
#' sample-estimated epsilon from `car::Anova()`. `power_calc` is unchanged and
#' always models the population-epsilon-adjusted test. Consequently, forcing
#' `sim_correction = "GG"` under a truly spherical population can make
#' `power_sim` slightly smaller than `power_calc`, because sample-epsilon GG
#' correction is mildly conservative under sphericity.
#'
#' Power is estimated for the prespecified corrected or uncorrected test.
#' Conditional procedures that first run Mauchly's test and then decide whether
#' to correct are not simulated.
#'
#' @section Examples:
#' ```{r, eval = FALSE}
#' power_curve(
#'   between = c(cond = 2),
#'   within = c(stim = 2),
#'   term = "cond:stim",
#'   target_pes = 0.14,
#'   n_range = c(16, 20, 23, 28), # n per between-subject cell
#'   n_sims = 1000,
#'   seed = 123
#' )
#'
#' power_curve(
#'   between = c(group = 2),
#'   within = c(time = 2),
#'   term = "group:time",
#'   target_pes = 0.14,
#'   n_range = c(12, 16, 20),
#'   n_sims = 5000,
#'   parallel = TRUE,
#'   cores = 4,
#'   seed = 123
#' )
#'
#' power_curve(
#'   within = c(time = 4),
#'   term = "time",
#'   target_pes = 0.15,
#'   n_range = 30,
#'   means_pattern = means_pattern(
#'     time = 1, value = 0,
#'     time = 2, value = 0.3,
#'     time = 3, value = 0.5,
#'     time = 4, value = 0.6
#'   ),
#'   n_sims = 1000,
#'   seed = 123
#' )
#' ```
#'
#' @export
power_curve <- function(between = NULL,
                        within = NULL,
                        term,
                        target_pes,
                        n_range,
                        n_sims = 10000,
                        alpha = 0.05,
                        ss_type = "III",
                        gpower = FALSE,
                        progress = interactive(),
                        parallel = FALSE,
                        cores = NULL,
                        seed = NULL,
                        covariance = NULL,
                        means_pattern = NULL,
                        sim_correction = c("auto", "GG", "none")) {
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
    means_pattern = means_pattern,
    sim_correction = sim_correction
  )
  message_long_serial_run(setup$n_sims, setup$parallel)

  if (!is.null(seed)) set.seed(seed)

  ns <- normalize_n_range(n_range)
  validate_calibration_n(ns, setup$spec, "n_range")
  progress_bar <- make_progress_bar(
    enabled = setup$progress,
    total = if (setup$parallel) length(ns) else length(ns) * setup$n_sims,
    label = "Simulating power"
  )
  on.exit(close_progress_bar(progress_bar), add = TRUE)

  run_one <- function(n) {
    row <- run_design_power_at_n(
      spec = setup$spec,
      term = setup$term,
      target_pes = target_pes,
      n = n,
      n_sims = setup$n_sims,
      alpha = alpha,
      ss_type = setup$ss_type,
      sim_correction_resolved = setup$sim_correction_resolved,
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
    row
  }

  curve <- purrr::map_dfr(ns, run_one)
  warn_power_disagreement(
    curve, setup$n_sims,
    sim_correction_resolved = setup$sim_correction_resolved,
    epsilon = setup$epsilon
  )

  structure(
    list(
      results = curve,
      term = setup$term,
      power = NA_real_,
      alpha = alpha,
      target_pes = target_pes,
      scale_factor = NA_real_,
      n_sims = setup$n_sims,
      n_needed = NA_integer_,
      total_n_needed = NA_integer_,
      gpower = setup$gpower,
      epsilon = setup$epsilon,
      covariance = setup$covariance,
      custom_covariance = setup$custom_covariance,
      custom_means_pattern = setup$custom_means_pattern,
      ss_type = setup$ss_type,
      sim_correction = setup$sim_correction,
      sim_correction_resolved = setup$sim_correction_resolved,
      design = setup$spec,
      call = match.call()
    ),
    class = "anovapowersim_curve"
  )
}


#' Search for the sample size needed for target ANOVA power
#'
#' Adaptive simulation search for the per-between-cell sample size needed to
#' reach a requested power for a balanced factorial ANOVA design. The search
#' searches upward from `n_start` until it brackets the target or reaches
#' `n_max`. If `n_start` already reaches the target, the search probes the
#' smallest sample size supported by the design to establish a lower bracket.
#' It then refines the bracket using interpolation with midpoint bisection as
#' a fallback.
#'
#' @inheritParams power_curve
#' @param power Desired target power.
#' @param n_start Starting sample size per between-subject cell, not a lower
#'   bound for the search. If `NULL`, an initial value is estimated from
#'   calculated power and constrained to values that support empirical
#'   calibration for the requested design.
#' @param n_max Maximum sample size per between-subject cell.
#' @param tol Acceptable precision above target power. If no simulated value at
#'   or above `power` is also no more than `power + tol`, `power_n()` warns that
#'   the requested precision band was not reached.
#'
#' @return An `anovapowersim_curve` object with `n_needed` and
#'   `total_n_needed`. For `power_n()`, `n_needed` is always an explicitly
#'   simulated `n_per_cell` value, never an interpolated sample size. If the
#'   search reaches target power but no simulated value lands inside
#'   `[power, power + tol]`, `power_n()` reports the smallest explicitly
#'   simulated value at or above target power and warns that the requested
#'   precision band was not reached.
#'
#' @section Examples:
#' ```{r, eval = FALSE}
#' power_n(
#'   between = c(cond = 2),
#'   within = c(stim = 4),
#'   term = "cond:stim",
#'   target_pes = 0.14,
#'   alpha = 0.05,
#'   power = 0.90,
#'   n_sims = 1000, # use 5000+ for a more precise estimate
#'   seed = 123 # for reproducibility
#' )
#' ```
#'
#' @export
power_n <- function(between = NULL,
                    within = NULL,
                    term,
                    target_pes,
                    power = 0.90,
                    n_sims = 10000,
                    alpha = 0.05,
                    ss_type = "III",
                    n_start = NULL,
                    n_max = 5000,
                    tol = 0.03,
                    gpower = FALSE,
                    progress = interactive(),
                    parallel = FALSE,
                    cores = NULL,
                    seed = NULL,
                    covariance = NULL,
                    means_pattern = NULL,
                    sim_correction = c("auto", "GG", "none")) {
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
    means_pattern = means_pattern,
    sim_correction = sim_correction
  )
  assert_unit_interval(power, "power")
  if (is.finite(power) && power < 0.90) {
    warning("Power greater than or equal to .90 is recommended.",
            call. = FALSE,
            immediate. = TRUE)
  }
  message_long_serial_run(setup$n_sims, setup$parallel)
  if (!is.numeric(n_max) || length(n_max) != 1L ||
      n_max < 1 || n_max != as.integer(n_max)) {
    stop("`n_max` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite number.", call. = FALSE)
  }

  min_n <- minimum_calibration_n(setup$spec)
  if (is.null(n_start)) {
    ncp_n <- estimate_ncp_n_needed(
      spec = setup$spec,
      term = setup$term,
      target_pes = target_pes,
      target_power = power,
      alpha = alpha,
      ss_type = setup$ss_type,
      sd = sd,
      r = r,
      covariance = setup$covariance,
      epsilon = setup$epsilon,
      gpower = setup$gpower,
      n_min = min_n,
      n_max = as.integer(n_max)
    )
    n_start <- if (is.na(ncp_n)) {
      min_n
    } else {
      max(min_n, as.integer(floor(ncp_n * 0.8)))
    }
  } else if (!is.numeric(n_start) || length(n_start) != 1L ||
             n_start < 1 || n_start != as.integer(n_start)) {
    stop("`n_start` must be a single positive integer.", call. = FALSE)
  }
  n_start <- as.integer(n_start)
  n_max <- as.integer(n_max)
  if (n_max < n_start) {
    stop("`n_max` must be greater than or equal to `n_start`.",
         call. = FALSE)
  }
  validate_calibration_n(n_start, setup$spec, "n_start")
  if (n_max < min_n) {
    stop(
      "`n_max` is too small for this design. This design has ",
      setup$spec$n_within_cells, " within-subject cell",
      if (setup$spec$n_within_cells == 1L) "" else "s",
      ", so the adaptive search needs `n_max >= ", min_n,
      "` per between-subject cell.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) set.seed(seed)

  progress_bar <- make_progress_bar(
    enabled = setup$progress,
    total = estimate_adaptive_progress_total(n_start, n_max, n_min = min_n),
    label = "Searching sample size"
  )
  on.exit(close_progress_bar(progress_bar), add = TRUE)

  run_one <- function(n) {
    run_design_power_at_n(
      spec = setup$spec,
      term = setup$term,
      target_pes = target_pes,
      n = n,
      n_sims = setup$n_sims,
      alpha = alpha,
      ss_type = setup$ss_type,
      sim_correction_resolved = setup$sim_correction_resolved,
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
  }

  curve <- adaptive_design_search(
    run_one = run_one,
    target = power,
    n_start = n_start,
    n_max = n_max,
    tol = tol,
    n_min = min_n,
    progress_bar = progress_bar
  )
  warn_power_disagreement(
    curve, setup$n_sims,
    sim_correction_resolved = setup$sim_correction_resolved,
    epsilon = setup$epsilon
  )

  n_needed <- estimate_design_n_needed(curve, target = power)
  warn_target_power_not_reached(
    n_needed = n_needed,
    target = power,
    n_max = n_max
  )
  warn_precision_band_not_reached(
    curve = curve,
    target = power,
    tol = tol,
    n_needed = n_needed
  )
  total_n_needed <- if (is.na(n_needed)) {
    NA_integer_
  } else {
    as.integer(n_needed * max(1L, setup$spec$n_between_cells))
  }

  structure(
    list(
      results = curve,
      term = setup$term,
      power = power,
      alpha = alpha,
      target_pes = target_pes,
      scale_factor = NA_real_,
      n_sims = setup$n_sims,
      n_needed = n_needed,
      total_n_needed = total_n_needed,
      gpower = setup$gpower,
      epsilon = setup$epsilon,
      covariance = setup$covariance,
      custom_covariance = setup$custom_covariance,
      custom_means_pattern = setup$custom_means_pattern,
      ss_type = setup$ss_type,
      sim_correction = setup$sim_correction,
      sim_correction_resolved = setup$sim_correction_resolved,
      design = setup$spec,
      call = match.call()
    ),
    class = "anovapowersim_curve"
  )
}


#' Define a sparse relative cell-mean pattern
#'
#' Creates a sparse mean-shape specification for the balanced simulation
#' functions. End each cell definition with `value`. Unlisted cells have raw
#' value zero, and factors in the tested term that are omitted from every row
#' are broadcast when the pattern is resolved against a design.
#'
#' Pattern values describe relative shape, not effect magnitude. The selected
#' power function projects the raw values onto the requested ANOVA term,
#' normalizes that component, and rescales it uniformly to reach `target_pes`.
#' Multiplying all values by one positive constant, adding an intercept or a
#' lower-order component, or reversing every sign therefore leaves the same
#' target-term direction (up to sign). Under nonsphericity, different
#' directions within a multi-df term can nevertheless produce different
#' simulated power.
#'
#' This differs from [cell_design()], where each `m` is a literal population
#' mean whose magnitude directly determines the simulated effect.
#'
#' @section Default direction:
#' When no pattern is supplied, balanced simulations use centered scores in
#' generated level order, `i - (L + 1) / 2` for levels `i = 1, ..., L`,
#' normalized to unit length. Interactions use the Kronecker product of their
#' factors' normalized score vectors, followed by one final normalization after
#' broadcasting. This is an ordered, reproducible convention rather than a
#' neutral scientific assumption; an explicit pattern is recommended whenever
#' the expected shape is known.
#'
#' @param ... Repeated named sparse-cell definitions. Each definition must use
#'   the same factor names in the same order and end in a finite numeric scalar
#'   named `value`. Factor levels may be supplied as one-based integer indices
#'   (for example, `time = 3`) or as the generated balanced-design names (for
#'   example, `time = "time3"`). The two forms are equivalent after the
#'   pattern is resolved against a design.
#'
#' @return An object of class `anovapowersim_means_pattern`, retaining the
#'   sparse definitions until a balanced simulation function resolves them
#'   against its design and tested term.
#'
#' @seealso [cell_design()] for unbalanced designs with literal cell means.
#'
#' @examples
#' trend <- means_pattern(
#'   time = 1, value = 0,
#'   time = 2, value = 0.3,
#'   time = 3, value = 0.5,
#'   time = 4, value = 0.6
#' )
#'
#' interaction_shape <- means_pattern(
#'   group = "group1", time = "time3", value = 1,
#'   group = "group2", time = "time3", value = -1
#' )
#'
#' @export
means_pattern <- function(...) {
  dots <- list(...)
  nms <- names(dots)
  if (!length(dots)) {
    stop(
      "Enter at least one sparse cell, ending each cell with `value`.",
      call. = FALSE
    )
  }
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("Every value in `...` must be named.", call. = FALSE)
  }

  rows <- list()
  current <- list()
  factor_names <- NULL
  for (i in seq_along(dots)) {
    nm <- nms[[i]]
    x <- dots[[i]]
    if (length(x) != 1L || is.null(x) || is.na(x)) {
      stop("`", nm, "` values must be single non-missing values.",
           call. = FALSE)
    }

    if (identical(nm, "value")) {
      if (!length(current)) {
        stop("Each `value` must follow one or more factor values.",
             call. = FALSE)
      }
      if (!is.numeric(x) || !is.finite(x)) {
        stop("`value` must be a finite numeric scalar.", call. = FALSE)
      }
      row_factor_names <- names(current)
      if (is.null(factor_names)) {
        factor_names <- row_factor_names
      } else if (!identical(row_factor_names, factor_names)) {
        stop("Every cell must use the same factor names in the same order.",
             call. = FALSE)
      }
      current$value <- as.numeric(x)
      rows[[length(rows) + 1L]] <- current
      current <- list()
    } else {
      if (nm %in% names(current)) {
        stop("Factor `", nm,
             "` appears more than once before the next `value`.",
             call. = FALSE)
      }
      current[[nm]] <- x
    }
  }

  if (length(current)) {
    stop("The final sparse cell is missing `value`.", call. = FALSE)
  }
  bad_names <- factor_names[make.names(factor_names) != factor_names]
  if (length(bad_names)) {
    stop(
      "Factor names must be syntactic R names. Problem name",
      if (length(bad_names) == 1L) "" else "s", ": ",
      paste(shQuote(bad_names), collapse = ", "), ".",
      call. = FALSE
    )
  }

  raw_key <- vapply(rows, function(row) {
    paste(vapply(factor_names, function(nm) {
      paste0(typeof(row[[nm]]), ":", as.character(row[[nm]]))
    }, character(1L)), collapse = "\r")
  }, character(1L))
  if (anyDuplicated(raw_key)) {
    stop("Each sparse cell may be defined only once.", call. = FALSE)
  }

  structure(
    list(definitions = rows, factor_names = factor_names),
    class = "anovapowersim_means_pattern"
  )
}


#' Create a balanced factorial ANOVA design specification
#'
#' Builds the design object used by [power_curve()], [design_term_means()], and
#' [simulate_design_dataset()]. This object stores factor names, level counts,
#' generated factor levels, and the between/within cell grids.
#'
#' @param between Named integer vector of between-subject factor level counts,
#'   e.g. `c(group = 2)`. Use `NULL` for no between-subject factors.
#' @param within Named integer vector of within-subject factor level counts,
#'   e.g. `c(time = 3)`. Use `NULL` for no within-subject factors.
#'
#' @return An object of class `anovapowersim_design_spec`.
#'
#' @examples
#' d <- balanced_anova_design(between = c(group = 2), within = c(time = 3))
#' d$between_cells
#' d$within_cells
#'
#' @export
balanced_anova_design <- function(between = NULL, within = NULL) {
  spec <- validate_design_spec(between, within)
  class(spec) <- c("anovapowersim_design_spec", class(spec))
  spec
}


#' Build calibrated means for a design term
#'
#' Projects a supplied relative pattern (or uses the documented
#' linear/Kronecker default) for one ANOVA term and scales it so an exact
#' reference dataset has the requested partial eta squared under the supplied
#' balanced design assumptions.
#'
#' @section Covariance limitation:
#' This manual helper does not accept [within_covariance()] specifications.
#' Calibration always uses the compound-symmetric covariance defined by `sd`
#' and `r`. Consequently, its calibrated means can differ from those used by
#' [power_curve()], [power_n()], or [power_achieved()] with a custom
#' `covariance`, because the covariance affects the reference residual sum of
#' squares and therefore the mean scale factor.
#'
#' @param design An `anovapowersim_design_spec` from [balanced_anova_design()].
#' @param term Character scalar naming the ANOVA term to target. Interaction
#'   terms are order-insensitive.
#' @param target_pes Target partial eta squared.
#' @param n Sample size per between-subject cell. For pure within designs, this
#'   is the total sample size.
#' @param sd Common outcome standard deviation.
#' @param r Compound-symmetric correlation among within-subject cells.
#' @param gpower Logical; if `TRUE`, calibrate to the G*Power-style
#'   noncentrality convention `lambda = total_n * f^2` (as in the "Cohen
#'   (1988)" option for within-subjects designs in G*Power). G*Power's
#'   estimates can differ from `target_pes`, especially for small samples or
#'   terms with more degrees of freedom; a warning is issued when
#'   `gpower = TRUE`. The default `gpower = FALSE` is recommended.
#' @param ss_type Sums-of-squares type for the tested ANOVA term. `"III"` is
#'   the default for order-invariant tests in unbalanced designs. Use `"I"` to
#'   reproduce sequential `stats::aov()` tests.
#' @param means_pattern Optional relative mean shape from [means_pattern()].
#'   Sparse values are projected onto `term` before calibration. If `NULL`, the
#'   deterministic linear/Kronecker default is used.
#'
#' @return A numeric matrix of cell means, with rows indexing between cells and
#'   columns indexing within cells.
#'
#' @examples
#' d <- balanced_anova_design(between = c(group = 2), within = c(time = 2))
#' design_term_means(d, term = "group:time", target_pes = 0.2, n = 20)
#'
#' @export
design_term_means <- function(design, term, target_pes, n, sd = 1, r = 0.5,
                              gpower = FALSE, ss_type = "III",
                              means_pattern = NULL) {
  assert_design_spec(design)
  term <- resolve_design_term(term, design)
  ss_type <- validate_ss_type(ss_type)
  assert_unit_interval(target_pes, "target_pes")
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != as.integer(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  validate_calibration_n(as.integer(n), design, "n")
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r <= -1 || r >= 1) {
    stop("`r` must be a single finite correlation in (-1, 1).",
         call. = FALSE)
  }
  if (!is.logical(gpower) || length(gpower) != 1L || is.na(gpower)) {
    stop("`gpower` must be TRUE or FALSE.", call. = FALSE)
  }
  warn_gpower_within_term_df(gpower = gpower)
  resolved_means_pattern <- resolve_means_pattern(
    means_pattern = means_pattern,
    spec = design,
    term = term
  )
  message_custom_means_pattern(
    custom_means_pattern = !is.null(means_pattern),
    term = term,
    target_pes = target_pes
  )
  calibrate_design_means(
    spec = design,
    term = term,
    target_pes = target_pes,
    n = as.integer(n),
    sd = sd,
    r = r,
    gpower = gpower,
    ss_type = ss_type,
    resolved_means_pattern = resolved_means_pattern
  )
}


#' Simulate data from a balanced ANOVA design
#'
#' Generates one long-format dataset from a balanced design. Supply means from
#' [design_term_means()] or any conformable matrix with one row per
#' between-subject cell and one column per within-subject cell.
#'
#' @section Covariance limitation:
#' This manual helper does not accept [within_covariance()] specifications. It
#' always simulates from the compound-symmetric covariance defined by `sd` and
#' `r`. It therefore cannot reproduce a balanced power-function call that uses
#' a custom `covariance`; use the power functions directly for that workflow.
#'
#' @param design An `anovapowersim_design_spec` from [balanced_anova_design()].
#' @param n Sample size per between-subject cell. For pure within designs, this
#'   is the total sample size.
#' @param means Numeric matrix of population cell means.
#' @param sd Common outcome standard deviation.
#' @param r Compound-symmetric correlation among within-subject cells.
#' @param empirical Logical; if `TRUE`, use `MASS::mvrnorm(empirical = TRUE)`
#'   so the generated sample closely matches the requested means/covariance.
#'
#' @return A tibble ready for `stats::aov()` with columns `id`, factor
#'   columns, and `value`.
#'
#' @examples
#' d <- balanced_anova_design(between = c(group = 2), within = c(time = 2))
#' m <- design_term_means(d, term = "group:time", target_pes = 0.2, n = 20)
#' sim <- simulate_design_dataset(d, n = 20, means = m)
#' head(sim)
#'
#' @export
simulate_design_dataset <- function(design, n, means, sd = 1, r = 0.5,
                                    empirical = FALSE) {
  assert_design_spec(design)
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != as.integer(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  if (missing(means)) {
    stop("`means` is required. Use design_term_means() to generate defaults.",
         call. = FALSE)
  }
  if (!is.matrix(means)) means <- as.matrix(means)
  expected_dim <- c(design$n_between_cells, design$n_within_cells)
  if (!identical(dim(means), expected_dim)) {
    stop("`means` must be a ", expected_dim[[1]], " x ", expected_dim[[2]],
         " matrix for this design.", call. = FALSE)
  }
  if (any(!is.finite(means))) {
    stop("`means` must contain only finite values.", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r <= -1 || r >= 1) {
    stop("`r` must be a single finite correlation in (-1, 1).",
         call. = FALSE)
  }
  if (!is.logical(empirical) || length(empirical) != 1L || is.na(empirical)) {
    stop("`empirical` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(empirical)) {
    validate_calibration_n(as.integer(n), design, "n")
  }

  simulate_balanced_design_data(
    spec = design,
    n = as.integer(n),
    means = means,
    sd = sd,
    r = r,
    empirical = empirical
  )
}


#' @keywords internal
#' @noRd
validate_design_spec <- function(between, within) {
  parse_counts <- function(x, arg) {
    if (is.null(x)) return(stats::setNames(integer(0), character(0)))
    if (!is.numeric(x) || is.null(names(x)) || any(names(x) == "")) {
      stop("`", arg, "` must be a named integer vector of level counts.",
           call. = FALSE)
    }
    bad_names <- names(x)[make.names(names(x)) != names(x)]
    if (length(bad_names)) {
      stop(
        "Factor names in `", arg, "` must be syntactic R names. ",
        "Problem name", if (length(bad_names) == 1L) "" else "s", ": ",
        paste(shQuote(bad_names), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    if (any(x < 2) || any(x != as.integer(x))) {
      stop("Every entry in `", arg, "` must be an integer >= 2.",
           call. = FALSE)
    }
    stats::setNames(as.integer(x), names(x))
  }

  between <- parse_counts(between, "between")
  within <- parse_counts(within, "within")
  all_names <- c(names(between), names(within))
  if (!length(all_names)) {
    stop("At least one between- or within-subject factor is required.",
         call. = FALSE)
  }
  if (anyDuplicated(all_names)) {
    stop("Factor names must be unique across `between` and `within`.",
         call. = FALSE)
  }

  levels <- stats::setNames(
    lapply(all_names, function(nm) {
      k <- c(between, within)[[nm]]
      factor(paste0(nm, seq_len(k)), levels = paste0(nm, seq_len(k)))
    }),
    all_names
  )

  between_cells <- if (length(between)) {
    tidyr::expand_grid(!!!levels[names(between)])
  } else {
    tibble::tibble(.dummy_between = factor("all"))
  }
  within_cells <- if (length(within)) {
    tidyr::expand_grid(!!!levels[names(within)])
  } else {
    tibble::tibble(.dummy_within = factor("dv"))
  }

  list(
    between = names(between),
    within = names(within),
    factor_names = all_names,
    level_counts = c(between, within),
    levels = levels,
    between_cells = between_cells,
    within_cells = within_cells,
    n_between_cells = if (length(between)) nrow(between_cells) else 1L,
    n_within_cells = if (length(within)) nrow(within_cells) else 1L
  )
}


#' @keywords internal
#' @noRd
validate_cell_count_values <- function(x, arg) {
  if (!is.numeric(x) || any(is.na(x)) || any(x < 1) ||
      any(x != as.integer(x))) {
    stop("`", arg, "` must contain positive integer counts.", call. = FALSE)
  }
  as.integer(x)
}


#' @keywords internal
#' @noRd
interaction_key <- function(data, columns) {
  if (!length(columns)) return(rep(".all", nrow(data)))
  do.call(paste, c(lapply(data[columns], as.character), sep = "\r"))
}


#' @keywords internal
#' @noRd
assert_design_spec <- function(design) {
  if (!inherits(design, "anovapowersim_design_spec")) {
    stop("`design` must be created by balanced_anova_design().",
         call. = FALSE)
  }
  invisible(design)
}


#' @keywords internal
#' @noRd
assert_term_name <- function(term) {
  if (!is.character(term) || length(term) != 1L || is.na(term) || !nzchar(term)) {
    stop("`term` must be a single non-empty character string.", call. = FALSE)
  }
  invisible(term)
}


#' @keywords internal
#' @noRd
resolve_design_term <- function(term, spec) {
  assert_term_name(term)

  requested <- strsplit(term, ":", fixed = TRUE)[[1L]]
  if (any(!nzchar(requested))) {
    stop("`term` must not contain empty ':' components.", call. = FALSE)
  }
  if (anyDuplicated(requested)) {
    stop("`term` must not repeat factor names.", call. = FALSE)
  }

  unknown <- setdiff(requested, spec$factor_names)
  if (length(unknown)) {
    stop("Unknown factor", if (length(unknown) > 1L) "s" else "", " in `term`: ",
         paste(shQuote(unknown), collapse = ", "), ". Available factors: ",
         paste(shQuote(spec$factor_names), collapse = ", "), call. = FALSE)
  }

  paste(spec$factor_names[spec$factor_names %in% requested], collapse = ":")
}


#' @keywords internal
#' @noRd
prepare_power_curve_inputs <- function(between, within, term, target_pes,
                                       n_sims, alpha, ss_type, sd, r,
                                       covariance = NULL, gpower, progress,
                                       parallel, cores,
                                       means_pattern = NULL,
                                       sim_correction = c("auto", "GG", "none")) {
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
    means_pattern = means_pattern,
    sim_correction = sim_correction
  )
  validate_target_pes(target_pes)
  message_custom_means_pattern(
    custom_means_pattern = setup$custom_means_pattern,
    term = setup$term,
    target_pes = target_pes
  )
  setup
}


#' Prepare inputs shared by balanced simulation-based power functions
#'
#' @keywords internal
#' @noRd
prepare_balanced_power_inputs <- function(between, within, term, n_sims,
                                          alpha, ss_type, sd, r,
                                          covariance = NULL, gpower, progress,
                                          parallel, cores,
                                          means_pattern = NULL,
                                          sim_correction = c("auto", "GG", "none")) {
  spec <- balanced_anova_design(between = between, within = within)
  term <- resolve_design_term(term, spec)
  if (!is.null(covariance) &&
      !inherits(covariance, "anovapowersim_covariance_spec")) {
    stop(
      "`covariance` must be created by within_covariance(); raw covariance ",
      "matrices are not accepted by balanced power functions because their ",
      "within-cell order can be ambiguous.",
      call. = FALSE
    )
  }
  resolved_means_pattern <- resolve_means_pattern(
    means_pattern = means_pattern,
    spec = spec,
    term = term
  )
  custom_means_pattern <- !is.null(means_pattern)
  resolved_covariance <- resolve_within_covariance(
    covariance = covariance,
    spec = spec,
    sd = sd,
    r = r
  )
  epsilon <- covariance_term_epsilon(
    covariance = resolved_covariance,
    spec = spec,
    term = term
  )
  ss_type <- validate_ss_type(ss_type)
  correction <- resolve_sim_correction(
    sim_correction = sim_correction,
    ss_type = ss_type,
    spec = spec,
    term = term,
    epsilon = epsilon
  )
  warn_uncorrected_nonsphericity(
    sim_correction = correction$requested,
    sim_correction_resolved = correction$resolved,
    ss_type = ss_type,
    epsilon = epsilon
  )
  warn_direction_sensitivity(
    spec = spec,
    term = term,
    epsilon = epsilon,
    custom_means_pattern = custom_means_pattern
  )
  assert_unit_interval(alpha, "alpha")
  if (!is.numeric(n_sims) || length(n_sims) != 1L || !is.finite(n_sims) ||
      n_sims < 1 || n_sims != as.integer(n_sims)) {
    stop("`n_sims` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r <= -1 || r >= 1) {
    stop("`r` must be a single finite correlation in (-1, 1).",
         call. = FALSE)
  }
  if (!is.logical(gpower) || length(gpower) != 1L || is.na(gpower)) {
    stop("`gpower` must be TRUE or FALSE.", call. = FALSE)
  }
  warn_gpower_within_term_df(gpower = gpower)
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  cores <- validate_parallel_cores(cores = cores, parallel = parallel)

  list(
    spec = spec,
    term = term,
    n_sims = as.integer(n_sims),
    ss_type = ss_type,
    covariance = resolved_covariance,
    custom_covariance = !is.null(covariance),
    means_pattern = resolved_means_pattern,
    custom_means_pattern = custom_means_pattern,
    epsilon = epsilon,
    sim_correction = correction$requested,
    sim_correction_resolved = correction$resolved,
    gpower = gpower,
    progress = progress,
    parallel = parallel,
    cores = cores
  )
}


#' Validate a target partial eta squared and issue the established advisory
#'
#' @keywords internal
#' @noRd
validate_target_pes <- function(target_pes) {
  assert_unit_interval(target_pes, "target_pes")
  if (target_pes == 0.06) {
    warning(
      paste(
        "It looks like you are using a rule-of-thumb \"medium\" effect size.",
        "This might overestimate the true effect size, rendering your study",
        "underpowered. Consider basing your power calculations on previous",
        "research or empirically-derived guidelines."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  invisible(target_pes)
}


#' @keywords internal
#' @noRd
normalize_n_range <- function(n_range) {
  if (!is.numeric(n_range) || any(n_range < 1) ||
      any(n_range != as.integer(n_range))) {
    stop("`n_range` must contain positive integers.", call. = FALSE)
  }
  sort(unique(as.integer(n_range)))
}


#' @keywords internal
#' @noRd
validate_ss_type <- function(ss_type) {
  if (!is.character(ss_type) || length(ss_type) != 1L || is.na(ss_type)) {
    stop("`ss_type` must be one of 'III', 'II', or 'I'.", call. = FALSE)
  }
  ss_type <- toupper(ss_type)
  if (!ss_type %in% c("III", "II", "I")) {
    stop("`ss_type` must be one of 'III', 'II', or 'I'.", call. = FALSE)
  }
  ss_type
}


#' Resolve the sphericity correction used for simulated tests
#'
#' @keywords internal
#' @noRd
resolve_sim_correction <- function(sim_correction, ss_type, spec, term,
                                   epsilon) {
  requested <- match.arg(sim_correction, c("auto", "GG", "none"))
  if (identical(requested, "GG") && identical(ss_type, "I")) {
    stop(
      "`sim_correction = \"GG\"` is not available with `ss_type = \"I\"` ",
      "because Type I tests do not provide Greenhouse-Geisser-corrected ",
      "p-values. Use `ss_type = \"II\"` or `\"III\"`, or set ",
      "`sim_correction = \"none\"`.",
      call. = FALSE
    )
  }

  has_multi_df_within <- within_term_df(spec = spec, term = term) > 1L
  resolved <- if (identical(requested, "GG")) {
    if (has_multi_df_within) "GG" else "none"
  } else if (identical(requested, "auto") &&
             isTRUE(epsilon < 1 - 1e-8) &&
             ss_type %in% c("II", "III")) {
    "GG"
  } else {
    "none"
  }

  list(requested = requested, resolved = resolved)
}


#' Warn when an uncorrected simulated test is used under nonsphericity
#'
#' @keywords internal
#' @noRd
warn_uncorrected_nonsphericity <- function(sim_correction,
                                           sim_correction_resolved,
                                           ss_type, epsilon) {
  if (!identical(sim_correction_resolved, "none") ||
      !isTRUE(epsilon < 1 - 1e-8)) {
    return(invisible(NULL))
  }

  type_i_note <- if (identical(ss_type, "I")) {
    paste0(
      " `ss_type = \"I\"` cannot provide GG-corrected p-values; use ",
      "`ss_type = \"II\"` or `\"III\"` with `sim_correction = \"GG\"` ",
      "to simulate the corrected test."
    )
  } else {
    " Set `sim_correction = \"GG\"` to simulate the corrected test."
  }

  warning(
    "The population covariance is nonspherical (Greenhouse-Geisser epsilon = ",
    signif(epsilon, 4), "), but the simulated test is uncorrected. Its true ",
    "Type I error exceeds the nominal alpha under this covariance, so ",
    "`power_sim` will typically exceed corrected power; that excess is alpha ",
    "inflation, not real power.",
    type_i_note,
    call. = FALSE,
    immediate. = TRUE
  )
  invisible(NULL)
}


#' Select the p-value used for one simulated ANOVA test
#'
#' @keywords internal
#' @noRd
select_simulated_p <- function(stats, use_gg_correction) {
  if (isTRUE(use_gg_correction)) {
    if (length(stats$p_value_gg) != 1L ||
        !is.finite(stats$p_value_gg)) {
      return(NA_real_)
    }
    return(as.numeric(stats$p_value_gg))
  }

  if (length(stats$p_value) != 1L || !is.finite(stats$p_value)) {
    return(NA_real_)
  }
  as.numeric(stats$p_value)
}


#' Assert that a reference fit supplied the required GG p-value
#'
#' @keywords internal
#' @noRd
assert_reference_gg_p <- function(stats, use_gg_correction, term) {
  if (!isTRUE(use_gg_correction)) return(invisible(NULL))
  if (length(stats$p_value_gg) == 1L && is.finite(stats$p_value_gg)) {
    return(invisible(NULL))
  }

  stop(
    "Internal error: Greenhouse-Geisser correction is required for term '",
    term, "', but the empirical reference fit did not return a finite ",
    "GG-corrected p-value. Please file a bug report at ",
    "https://github.com/shaheedazaad/anovapowersim/issues and include the ",
    "function call.",
    call. = FALSE
  )
}


#' Warn when the default direction is consequential for simulated power
#'
#' @keywords internal
#' @noRd
warn_direction_sensitivity <- function(spec, term, epsilon,
                                       custom_means_pattern) {
  if (isTRUE(custom_means_pattern) ||
      within_term_df(spec = spec, term = term) <= 1L ||
      !isTRUE(epsilon < 1 - 1e-8)) {
    return(invisible(NULL))
  }
  warning(
    "This multi-df within-subject term is nonspherical, so simulated power ",
    "depends on the relative cell-mean pattern. `power_sim` currently uses ",
    "the package's default linear/Kronecker pattern; a different pattern ",
    "with the same `target_pes` and covariance can produce different power. ",
    "Supply `means_pattern` to describe the expected relative mean shape.",
    call. = FALSE,
    immediate. = TRUE
  )
  invisible(NULL)
}


#' Explain the shape-only semantics of custom balanced means patterns
#'
#' @keywords internal
#' @noRd
message_custom_means_pattern <- function(custom_means_pattern, term,
                                         target_pes = NULL) {
  if (!isTRUE(custom_means_pattern)) return(invisible(NULL))

  rescaling <- if (is.null(target_pes)) {
    "It will be rescaled to each candidate `target_pes` during the search."
  } else {
    paste0(
      "It will be rescaled to `target_pes = ",
      format(target_pes, trim = TRUE), "`."
    )
  }
  rlang::inform(
    paste0(
      "Custom `means_pattern` values are shape-only. The pattern was ",
      "projected onto term '", term, "' and normalized; its supplied ",
      "magnitude and any components outside term '", term,
      "' are discarded. ",
      rescaling, " By contrast, `cell_design()` uses literal `m` values."
    ),
    .frequency = "once",
    .frequency_id = "anovapowersim_custom_means_pattern_semantics"
  )
  invisible(NULL)
}


#' Warn that gpower = TRUE may not calibrate target_pes exactly
#'
#' @keywords internal
#' @noRd
warn_gpower_within_term_df <- function(gpower) {
  if (!isTRUE(gpower)) return(invisible(NULL))
  warning(
    "`gpower = TRUE` calibrates means to G*Power's noncentrality ",
    "convention, so the partial eta squared actually achieved can differ ",
    "from `target_pes` -- this is more pronounced for small samples and ",
    "terms with more degrees of freedom. The default `gpower = FALSE` is ",
    "recommended if you want `target_pes` to match your reported or ",
    "expected partial eta squared exactly.",
    call. = FALSE,
    immediate. = TRUE
  )
  invisible(NULL)
}


#' @keywords internal
#' @noRd
parallel_worker_helpers <- function(names) {
  env <- new.env(parent = globalenv())
  for (nm in names) {
    value <- get(nm, mode = "function")
    if (is.function(value)) environment(value) <- env
    assign(nm, value, envir = env)
  }
  stats::setNames(lapply(names, get, envir = env), names)
}


#' @keywords internal
#' @noRd
validate_parallel_cores <- function(cores, parallel) {
  available <- as.integer(future::availableCores()[[1L]])
  if (!is.finite(available) || available < 1L) available <- 1L

  if (!is.null(cores) &&
      (!is.numeric(cores) || length(cores) != 1L ||
       !is.finite(cores) || cores < 1 || cores != as.integer(cores))) {
    stop("`cores` must be a single positive integer.", call. = FALSE)
  }

  if (!is.null(cores)) {
    cores <- as.integer(cores)
    if (cores > available) {
      stop(
        "`cores` must not exceed the number of available cores (",
        available, ").",
        call. = FALSE
      )
    }
    return(cores)
  }

  if (!isTRUE(parallel)) return(NULL)

  cores <- max(1L, available - 1L)
  message(
    "`parallel = TRUE` and `cores` was not set; using ",
    cores,
    " of ",
    available,
    " available core",
    if (available == 1L) "" else "s",
    "."
  )
  cores
}


#' @keywords internal
#' @noRd
message_long_serial_run <- function(n_sims, parallel) {
  if (isTRUE(parallel) || n_sims < 5000L) return(invisible(NULL))
  message(
    "`n_sims` is 5000 or more and `parallel = FALSE`; this may take a ",
    "while. Consider setting `parallel = TRUE` to run simulations in parallel."
  )
  invisible(NULL)
}


#' @keywords internal
#' @noRd
run_design_power_at_n <- function(spec, term, target_pes, n, n_sims,
                                  alpha, ss_type, sim_correction_resolved,
                                  sd, r, gpower,
                                  covariance = NULL,
                                  epsilon = 1,
                                  progress_bar = NULL,
                                  parallel = FALSE, cores = NULL,
                                  resolved_means_pattern = NULL) {
  means <- calibrate_design_means(
    spec = spec,
    term = term,
    target_pes = target_pes,
    n = n,
    sd = sd,
    r = r,
    covariance = covariance,
    gpower = gpower,
    ss_type = ss_type,
    resolved_means_pattern = resolved_means_pattern
  )
  sanity <- sanity_check_term_effect(
    spec = spec,
    term = term,
    target_pes = target_pes,
    n = n,
    means = means,
    sd = sd,
    r = r,
    covariance = covariance,
    alpha = alpha,
    gpower = gpower,
    epsilon = epsilon,
    ss_type = ss_type
  )
  use_gg_correction <- identical(sim_correction_resolved, "GG")
  assert_reference_gg_p(
    stats = sanity,
    use_gg_correction = use_gg_correction,
    term = term
  )

  helpers <- parallel_worker_helpers(c(
    "simulate_balanced_design_data",
    "fit_design_term_stats",
    "validate_ss_type",
    "fit_design_model",
    "fit_car_term_stats",
    "car_gg_p_values",
    "extract_term_stats",
    "set_sum_contrasts",
    "design_aov_formula",
    "extract_aov_rows",
    "car_between_rows",
    "car_repeated_rows",
    "repeated_measures_wide_data",
    "make_cell_labels",
    "interaction_key",
    "validate_calibration_n",
    "minimum_calibration_n",
    "compound_symmetric_sigma",
    "select_simulated_p",
    "tick_progress_bar"
  ))
  fit_one_stats <- helpers$fit_design_term_stats
  simulate_one_dataset <- helpers$simulate_balanced_design_data
  select_p <- helpers$select_simulated_p
  tick_progress <- helpers$tick_progress_bar
  simulate_one <- function(i) {
    sim <- simulate_one_dataset(
      spec = spec,
      n = n,
      means = means,
      sd = sd,
      r = r,
      covariance = covariance,
      empirical = FALSE
    )
    tick_progress(progress_bar)
    stats <- tryCatch(
      fit_one_stats(sim, spec, term, ss_type = ss_type),
      error = function(e) NULL
    )
    if (is.null(stats)) return(NA)
    p_value <- select_p(stats, use_gg_correction)
    if (!is.finite(p_value)) return(NA)
    isTRUE(p_value < alpha)
  }

  successes <- if (parallel) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
    unlist(
      future.apply::future_lapply(
        seq_len(n_sims),
        simulate_one,
        future.seed = TRUE
      ),
      use.names = FALSE
    )
  } else {
    purrr::map_lgl(seq_len(n_sims), simulate_one)
  }

  failed_count <- sum(is.na(successes))
  valid_count <- n_sims - failed_count
  if (valid_count == 0L) {
    stop(
      "All simulated ANOVA fits failed. This usually indicates an internal ",
      "model-fitting error rather than zero power.",
      call. = FALSE
    )
  }
  if (failed_count > 0L) {
    warning(
      failed_count, " of ", n_sims,
      " simulated ANOVA fits failed and were excluded from `power_sim`.",
      call. = FALSE
    )
  }
  success_count <- sum(successes, na.rm = TRUE)
  power <- if (valid_count > 0L) success_count / valid_count else NA_real_

  tibble::tibble(
    n_per_cell = as.integer(n),
    total_n = as.integer(n * max(1L, spec$n_between_cells)),
    n_sims = n_sims,
    valid_sims = as.integer(valid_count),
    failed_sims = as.integer(failed_count),
    epsilon = epsilon,
    num_df = sanity$num_df,
    den_df = sanity$den_df,
    ncp = round(sanity$ncp, 3),
    power_calc = round(sanity$power_calc, 3),
    power_sim = power
  )
}


#' @keywords internal
#' @noRd
adaptive_design_search <- function(run_one, target, n_start, n_max, tol,
                                   max_iter = 25L, progress_bar = NULL,
                                   n_min = 1L) {
  visited <- list()
  n <- max(1L, as.integer(n_start))
  n_min <- max(1L, min(n, as.integer(n_min)))
  lo <- NULL
  lo_p <- NA_real_
  hi <- NULL
  hi_p <- NA_real_

  repeat {
    row <- run_one(n)
    tick_progress_bar(progress_bar)
    visited[[length(visited) + 1L]] <- row
    p <- row$power_sim
    if (!is.na(p) && p >= target) {
      hi <- n
      hi_p <- p
      if (is.null(lo) && n > n_min) {
        lower_row <- run_one(n_min)
        tick_progress_bar(progress_bar)
        visited[[length(visited) + 1L]] <- lower_row
        lower_p <- lower_row$power_sim
        if (!is.na(lower_p) && lower_p >= target) {
          hi <- n_min
          hi_p <- lower_p
        } else {
          lo <- n_min
          lo_p <- lower_p
        }
      }
      break
    }
    if (n >= n_max) {
      break
    }
    if (!is.na(p)) {
      lo <- n
      lo_p <- p
    }
    n <- min(n_max, max(n + 1L, n * 2L))
  }

  if (!is.null(lo) && !is.null(hi)) {
    lo_n <- lo
    hi_n <- hi
    iter <- 0L
    while (hi_n > lo_n + 1L && iter < max_iter) {
      visited_n <- vapply(
        visited,
        function(x) as.integer(x$n_per_cell[[1L]]),
        integer(1L)
      )
      next_n <- next_adaptive_n(
        lo_n = lo_n,
        lo_p = lo_p,
        hi_n = hi_n,
        hi_p = hi_p,
        target = target,
        visited_n = visited_n
      )
      if (is.na(next_n)) break

      row <- run_one(next_n)
      tick_progress_bar(progress_bar)
      visited[[length(visited) + 1L]] <- row
      p <- row$power_sim
      if (is.na(p)) break
      if (p >= target) {
        hi_n <- next_n
        hi_p <- p
      } else {
        lo_n <- next_n
        lo_p <- p
      }
      iter <- iter + 1L
    }
  }

  dplyr::bind_rows(visited) |>
    dplyr::arrange(.data$n_per_cell) |>
    dplyr::distinct(.data$n_per_cell, .keep_all = TRUE)
}


#' @keywords internal
#' @noRd
estimate_design_n_needed <- function(curve, target) {
  curve <- dplyr::arrange(curve, .data$n_per_cell)
  above <- which(curve$power_sim >= target)
  if (length(above) == 0L) return(NA_integer_)
  as.integer(curve$n_per_cell[[above[1L]]])
}


#' @keywords internal
#' @noRd
power_is_in_precision_band <- function(power_sim, target, tol) {
  is.finite(power_sim) && power_sim >= target && power_sim <= target + tol
}


#' @keywords internal
#' @noRd
next_adaptive_n <- function(lo_n, lo_p, hi_n, hi_p, target, visited_n) {
  candidates <- unsimulated_bracket_n(lo_n, hi_n, visited_n)
  if (!length(candidates)) return(NA_integer_)

  interp <- NA_integer_
  if (is.finite(lo_p) && is.finite(hi_p) && hi_p > lo_p) {
    frac <- (target - lo_p) / (hi_p - lo_p)
    raw_interp <- lo_n + frac * (hi_n - lo_n)
    interp <- as.integer(ceiling(raw_interp - sqrt(.Machine$double.eps)))
    interp <- min(max(interp, lo_n + 1L), hi_n - 1L)
  }
  if (!is.na(interp)) {
    return(as.integer(candidates[which.min(abs(candidates - interp))]))
  }

  mid <- as.integer(floor((lo_n + hi_n) / 2L))
  as.integer(candidates[which.min(abs(candidates - mid))])
}


#' @keywords internal
#' @noRd
unsimulated_bracket_n <- function(lo_n, hi_n, visited_n) {
  if (hi_n <= lo_n + 1L) return(integer())
  candidates <- seq.int(lo_n + 1L, hi_n - 1L)
  candidates[!(candidates %in% visited_n)]
}


#' @keywords internal
#' @noRd
warn_precision_band_not_reached <- function(curve, target, tol, n_needed) {
  if (is.na(n_needed)) return(invisible(NULL))
  if (any(vapply(
    curve$power_sim,
    power_is_in_precision_band,
    logical(1L),
    target = target,
    tol = tol
  ))) {
    return(invisible(NULL))
  }

  reported <- curve[curve$n_per_cell == n_needed, , drop = FALSE]
  if (!nrow(reported)) return(invisible(NULL))
  warning(
    sprintf(
      paste(
        "Requested precision band was not reached for target power %.3f",
        "with tolerance %.3f; reporting the closest explicitly simulated",
        "n_per_cell at or above target power: %d (power_sim = %.3f)."
      ),
      target,
      tol,
      as.integer(reported$n_per_cell[[1L]]),
      reported$power_sim[[1L]]
    ),
    call. = FALSE
  )
  invisible(NULL)
}


#' @keywords internal
#' @noRd
warn_power_disagreement <- function(results, n_sims, threshold = 0.05,
                                    sim_correction_resolved = NULL,
                                    epsilon = 1) {
  if (!all(c("power_sim", "power_calc") %in% names(results))) return(invisible(NULL))
  diff <- abs(results$power_sim - results$power_calc)
  bad <- which(is.finite(diff) & diff > threshold)
  if (!length(bad)) return(invisible(NULL))
  if (identical(sim_correction_resolved, "none") &&
      isTRUE(epsilon < 1 - 1e-8)) {
    return(invisible(NULL))
  }

  point_label <- if ("target_pes" %in% names(results)) {
    "target_pes"
  } else if ("n_per_cell" %in% names(results)) {
    "n_per_cell"
  } else {
    "total_n"
  }
  points <- paste(results[[point_label]][bad], collapse = ", ")
  max_diff <- max(diff[bad], na.rm = TRUE)
  msg <- paste0(
    "`power_sim` and `power_calc` differ by more than ",
    threshold * 100,
    " percentage points for ", point_label, " = ",
    points,
    " (largest difference = ",
    sprintf("%.3f", max_diff),
    "). "
  )
  if (n_sims < 5000) {
    msg <- paste0(msg, "Try increasing `n_sims` for a more stable simulation estimate.")
  } else {
    msg <- paste0(
      msg,
      "You are already using 5000+ simulations; consider raising a GitHub issue ",
      "and treat this power estimate as potentially unstable."
    )
  }
  warning(msg, call. = FALSE)
  invisible(NULL)
}


#' Build the complete design grid in the package's canonical cell order
#'
#' @keywords internal
#' @noRd
complete_design_grid <- function(spec) {
  grid <- tidyr::expand_grid(!!!spec$levels[spec$factor_names])
  for (nm in spec$factor_names) {
    grid[[nm]] <- factor(grid[[nm]], levels = levels(spec$levels[[nm]]))
  }
  grid
}


#' Select the sum-contrast model-matrix columns for one design term
#'
#' @keywords internal
#' @noRd
term_model_matrix <- function(spec, term, grid = complete_design_grid(spec)) {
  for (nm in spec$factor_names) {
    stats::contrasts(grid[[nm]]) <- stats::contr.sum(nlevels(grid[[nm]]))
  }
  form <- stats::reformulate(paste(spec$factor_names, collapse = "*"))
  mm <- stats::model.matrix(form, data = grid)
  term_labels <- attr(stats::terms(form), "term.labels")
  idx <- match(term, term_labels)
  if (is.na(idx)) {
    stop("Term '", term, "' is not available. Available terms: ",
         paste(shQuote(term_labels), collapse = ", "), call. = FALSE)
  }
  cols <- which(attr(mm, "assign") == idx)
  mm[, cols, drop = FALSE]
}


#' Align a complete-grid vector to the package's means-matrix layout
#'
#' @keywords internal
#' @noRd
align_grid_pattern <- function(pattern, grid, spec, label = "mean") {
  if (length(pattern) != nrow(grid)) {
    stop("Internal ", label, " pattern length does not match the design grid.",
         call. = FALSE)
  }

  between_key <- if (length(spec$between)) {
    do.call(paste, c(grid[spec$between], sep = "\r"))
  } else {
    rep(".all_between", nrow(grid))
  }
  within_key <- if (length(spec$within)) {
    do.call(paste, c(grid[spec$within], sep = "\r"))
  } else {
    rep(".all_within", nrow(grid))
  }
  between_levels <- if (length(spec$between)) {
    do.call(paste, c(spec$between_cells[spec$between], sep = "\r"))
  } else {
    ".all_between"
  }
  within_levels <- if (length(spec$within)) {
    do.call(paste, c(spec$within_cells[spec$within], sep = "\r"))
  } else {
    ".all_within"
  }

  out <- matrix(
    NA_real_,
    nrow = spec$n_between_cells,
    ncol = spec$n_within_cells
  )
  for (i in seq_along(pattern)) {
    row <- match(between_key[[i]], between_levels)
    col <- match(within_key[[i]], within_levels)
    out[row, col] <- pattern[[i]]
  }
  if (anyNA(out)) {
    stop("Could not align the ", label, " pattern to the design cells.",
         call. = FALSE)
  }
  out
}


#' Resolve one balanced factor level to its canonical one-based index
#'
#' @keywords internal
#' @noRd
resolve_pattern_level <- function(value, factor_name, spec) {
  n_levels <- spec$level_counts[[factor_name]]
  generated <- levels(spec$levels[[factor_name]])
  valid_message <- paste0(
    "Valid levels for factor `", factor_name, "` are indices 1 through ",
    n_levels, " or generated names ",
    paste(shQuote(generated), collapse = ", "), "."
  )

  if (is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value == floor(value) && value >= 1 && value <= n_levels) {
    return(as.integer(value))
  }
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    idx <- match(value, generated)
    if (!is.na(idx)) return(as.integer(idx))
  }
  stop(
    "Invalid level ", shQuote(as.character(value)), " for factor `",
    factor_name, "`. ", valid_message,
    call. = FALSE
  )
}


#' Resolve and project a sparse relative means pattern
#'
#' @keywords internal
#' @noRd
resolve_means_pattern <- function(means_pattern, spec, term) {
  if (is.null(means_pattern)) return(NULL)
  if (!inherits(means_pattern, "anovapowersim_means_pattern")) {
    stop("`means_pattern` must be created by means_pattern().", call. = FALSE)
  }
  definitions <- means_pattern$definitions
  factor_names <- means_pattern$factor_names
  if (!is.list(definitions) || !length(definitions) ||
      !is.character(factor_names) || !length(factor_names)) {
    stop("`means_pattern` contains invalid sparse definitions.", call. = FALSE)
  }

  term_factors <- strsplit(term, ":", fixed = TRUE)[[1L]]
  outside_term <- setdiff(factor_names, term_factors)
  if (length(outside_term)) {
    stop(
      "Factor", if (length(outside_term) == 1L) "" else "s", " ",
      paste(shQuote(outside_term), collapse = ", "),
      " in `means_pattern` ",
      if (length(outside_term) == 1L) "is" else "are",
      " not part of term '", term, "'. Pattern factors must belong to the ",
      "tested term.",
      call. = FALSE
    )
  }

  canonical <- vector("list", length(definitions))
  values <- numeric(length(definitions))
  for (i in seq_along(definitions)) {
    row <- definitions[[i]]
    if (!is.list(row) || !identical(names(row), c(factor_names, "value"))) {
      stop("`means_pattern` contains invalid sparse definitions.", call. = FALSE)
    }
    value <- row$value
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("`means_pattern` values must be finite numeric scalars.",
           call. = FALSE)
    }
    canonical[[i]] <- stats::setNames(
      lapply(factor_names, function(nm) {
        resolve_pattern_level(row[[nm]], factor_name = nm, spec = spec)
      }),
      factor_names
    )
    values[[i]] <- as.numeric(value)
  }
  canonical_key <- vapply(canonical, function(row) {
    paste(unlist(row, use.names = FALSE), collapse = "\r")
  }, character(1L))
  if (anyDuplicated(canonical_key)) {
    stop(
      "Each sparse cell may be defined only once; indices and generated ",
      "level names that identify the same cell are duplicates.",
      call. = FALSE
    )
  }

  grid <- complete_design_grid(spec)
  raw <- numeric(nrow(grid))
  for (i in seq_along(canonical)) {
    selected <- rep(TRUE, nrow(grid))
    for (nm in factor_names) {
      selected <- selected & as.integer(grid[[nm]]) == canonical[[i]][[nm]]
    }
    raw[selected] <- values[[i]]
  }

  term_mm <- term_model_matrix(spec = spec, term = term, grid = grid)
  term_qr <- qr(term_mm)
  term_basis <- qr.Q(term_qr, complete = FALSE)
  projected <- as.numeric(term_basis %*% crossprod(term_basis, raw))
  projected_norm <- sqrt(sum(projected^2))
  zero_tolerance <- sqrt(.Machine$double.eps) * max(1, sqrt(sum(raw^2)))
  if (!is.finite(projected_norm) || projected_norm <= zero_tolerance) {
    stop(
      "The supplied `means_pattern` has a zero projection onto term '",
      term, "'. Supply values with a nonzero component for the tested term.",
      call. = FALSE
    )
  }
  projected <- projected / projected_norm
  align_grid_pattern(projected, grid = grid, spec = spec, label = "resolved")
}


#' Normalized centered linear scores for one factor
#'
#' @keywords internal
#' @noRd
centered_linear_scores <- function(n_levels) {
  scores <- seq_len(n_levels) - (n_levels + 1) / 2
  scores / sqrt(sum(scores^2))
}


#' Deterministic normalized linear/Kronecker pattern for a design term
#'
#' @keywords internal
#' @noRd
default_term_pattern <- function(spec, term) {
  grid <- complete_design_grid(spec)
  term_factors <- strsplit(term, ":", fixed = TRUE)[[1L]]
  pattern <- rep(1, nrow(grid))
  for (nm in term_factors) {
    scores <- centered_linear_scores(spec$level_counts[[nm]])
    pattern <- pattern * scores[as.integer(grid[[nm]])]
  }
  pattern_norm <- sqrt(sum(pattern^2))
  if (!is.finite(pattern_norm) || pattern_norm <= sqrt(.Machine$double.eps)) {
    stop("The default linear/Kronecker pattern for term '", term,
         "' is zero.", call. = FALSE)
  }
  pattern <- pattern / pattern_norm
  align_grid_pattern(pattern, grid = grid, spec = spec, label = "default")
}


#' @keywords internal
#' @noRd
calibrate_design_means <- function(spec, term, target_pes, n, sd, r,
                                   covariance = NULL, gpower = FALSE,
                                   ss_type = "III",
                                   resolved_means_pattern = NULL) {
  base <- if (is.null(resolved_means_pattern)) {
    default_term_pattern(spec, term)
  } else {
    resolved_means_pattern
  }
  exact <- simulate_balanced_design_data(
    spec = spec,
    n = n,
    means = base,
    sd = sd,
    r = r,
    covariance = covariance,
    empirical = TRUE
  )
  stats <- fit_design_term_stats(exact, spec, term, ss_type = ss_type)
  old_pes <- stats$pes
  if (is.na(old_pes) || old_pes <= 0 || old_pes >= 1) {
    stop("Could not calibrate the means pattern for term '", term, "'.",
         call. = FALSE)
  }
  calibration_pes <- calibration_pes_for_ncp(
    target_pes = target_pes,
    total_n = n * max(1L, spec$n_between_cells),
    num_df = stats$num_df,
    den_df = stats$den_df,
    gpower = gpower
  )
  k <- compute_scale_factor(old_pes, calibration_pes)
  base * k
}


#' @keywords internal
#' @noRd
calibration_pes_for_ncp <- function(target_pes, total_n, num_df, den_df,
                                    gpower) {
  if (!isTRUE(gpower)) return(target_pes)
  f2 <- target_pes / (1 - target_pes)
  target_ncp <- total_n * f2
  calibration_f2 <- target_ncp / den_df
  calibration_pes <- calibration_f2 / (1 + calibration_f2)
  if (!is.finite(calibration_pes) || calibration_pes <= 0 ||
      calibration_pes >= 1) {
    stop("Could not calibrate G*Power-style means for this design.",
         call. = FALSE)
  }
  calibration_pes
}


#' @keywords internal
#' @noRd
sanity_check_term_effect <- function(spec, term, target_pes, n, means, sd, r,
                                     alpha, gpower, covariance = NULL,
                                     epsilon = 1,
                                     ss_type = "III") {
  exact <- simulate_balanced_design_data(
    spec = spec,
    n = n,
    means = means,
    sd = sd,
    r = r,
    covariance = covariance,
    empirical = TRUE
  )
  stats <- fit_design_term_stats(exact, spec, term, ss_type = ss_type)
  total_n <- n * max(1L, spec$n_between_cells)
  expected_pes <- calibration_pes_for_ncp(
    target_pes = target_pes,
    total_n = total_n,
    num_df = stats$num_df,
    den_df = stats$den_df,
    gpower = gpower
  )
  tolerance <- max(1e-6, sqrt(.Machine$double.eps) * 10)
  if (!is.finite(stats$pes) || abs(stats$pes - expected_pes) > tolerance) {
    stop(
      "Sanity check failed: empirical reference data produced partial eta ",
      "squared ", signif(stats$pes, 6), " for term '", term,
      "', expected ", signif(expected_pes, 6), ".",
      call. = FALSE
    )
  }
  uncorrected_ncp <- ncp_from_pes(
    pes = target_pes,
    total_n = total_n,
    den_df = stats$den_df,
    gpower = gpower
  )
  stats$num_df <- epsilon * stats$num_df
  stats$den_df <- epsilon * stats$den_df
  stats$ncp <- epsilon * uncorrected_ncp
  stats$power_calc <- stats::pf(
    stats::qf(1 - alpha, stats$num_df, stats$den_df),
    stats$num_df,
    stats$den_df,
    ncp = stats$ncp,
    lower.tail = FALSE
  )
  stats
}


#' @keywords internal
#' @noRd
ncp_from_pes <- function(pes, total_n, den_df, gpower) {
  f2 <- pes / (1 - pes)
  if (isTRUE(gpower)) total_n * f2 else den_df * f2
}


#' @keywords internal
#' @noRd
estimate_ncp_n_needed <- function(spec, term, target_pes, target_power, alpha,
                                  ss_type, sd, r, gpower, n_min, n_max,
                                  covariance = NULL, epsilon = 1) {
  power_at <- function(n) {
    power_calc_at_n(
      spec = spec,
      term = term,
      target_pes = target_pes,
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

  p_min <- power_at(n_min)
  if (!is.finite(p_min)) return(NA_integer_)
  if (p_min >= target_power) return(as.integer(n_min))

  lo <- n_min
  hi <- min(n_max, max(n_min + 1L, n_min * 2L))
  repeat {
    p_hi <- power_at(hi)
    if (!is.finite(p_hi)) return(NA_integer_)
    if (p_hi >= target_power || hi >= n_max) break
    lo <- hi
    hi <- min(n_max, max(hi + 1L, hi * 2L))
  }
  if (hi >= n_max && power_at(hi) < target_power) return(NA_integer_)

  while (hi > lo + 1L) {
    mid <- as.integer(floor((lo + hi) / 2L))
    p_mid <- power_at(mid)
    if (!is.finite(p_mid)) return(NA_integer_)
    if (p_mid < target_power) lo <- mid else hi <- mid
  }
  as.integer(hi)
}


#' @keywords internal
#' @noRd
power_calc_at_n <- function(spec, term, target_pes, n, alpha, ss_type, sd, r,
                            gpower, covariance = NULL, epsilon = 1) {
  base <- default_term_pattern(spec, term)
  exact <- simulate_balanced_design_data(
    spec = spec,
    n = n,
    means = base,
    sd = sd,
    r = r,
    covariance = covariance,
    empirical = TRUE
  )
  term_stats <- fit_design_term_stats(exact, spec, term, ss_type = ss_type)
  uncorrected_ncp <- ncp_from_pes(
    pes = target_pes,
    total_n = n * max(1L, spec$n_between_cells),
    den_df = term_stats$den_df,
    gpower = gpower
  )
  num_df <- epsilon * term_stats$num_df
  den_df <- epsilon * term_stats$den_df
  ncp <- epsilon * uncorrected_ncp
  stats::pf(
    stats::qf(1 - alpha, num_df, den_df),
    num_df,
    den_df,
    ncp = ncp,
    lower.tail = FALSE
  )
}


#' @keywords internal
#' @noRd
simulate_balanced_design_data <- function(spec, n, means, sd, r,
                                          empirical = FALSE,
                                          covariance = NULL) {
  if (isTRUE(empirical)) {
    validate_calibration_n(n, spec, "n")
  }

  sigma <- if (is.null(covariance)) {
    compound_symmetric_sigma(spec$n_within_cells, sd = sd, r = r)
  } else {
    covariance
  }
  between_labels <- make_cell_labels("b", spec$n_between_cells)
  within_labels <- make_cell_labels("w", spec$n_within_cells)
  rownames(means) <- between_labels
  colnames(means) <- within_labels

  subject_rows <- vector("list", spec$n_between_cells)
  y_rows <- vector("list", spec$n_between_cells)

  for (i in seq_len(spec$n_between_cells)) {
    y <- MASS::mvrnorm(
      n = n,
      mu = means[i, ],
      Sigma = sigma,
      empirical = empirical
    )
    if (spec$n_within_cells == 1L) {
      y <- matrix(y, nrow = n, ncol = 1L)
    } else if (n == 1L) {
      y <- matrix(y, nrow = 1L)
    }
    colnames(y) <- within_labels

    b <- spec$between_cells[i, spec$between, drop = FALSE]
    subject_rows[[i]] <- dplyr::bind_cols(
      tibble::tibble(id = seq_len(n) + (i - 1L) * n),
      b[rep(1L, n), , drop = FALSE]
    )
    y_rows[[i]] <- y
  }

  subjects <- dplyr::bind_rows(subject_rows)
  y_mat <- do.call(rbind, y_rows)
  wide <- dplyr::bind_cols(subjects, tibble::as_tibble(y_mat))

  if (!length(spec$within)) {
    out <- wide |>
      dplyr::rename(value = dplyr::all_of(within_labels[1L]))
  } else {
    within_map <- dplyr::bind_cols(
      tibble::tibble(.within_cell = within_labels),
      spec$within_cells[, spec$within, drop = FALSE]
    )
    out <- wide |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(within_labels),
        names_to = ".within_cell",
        values_to = "value"
      ) |>
      dplyr::left_join(within_map, by = ".within_cell") |>
      dplyr::select(-".within_cell")
  }

  out$id <- factor(out$id)
  for (nm in spec$factor_names) {
    out[[nm]] <- factor(out[[nm]], levels = levels(spec$levels[[nm]]))
  }
  tibble::as_tibble(out)
}


#' @keywords internal
#' @noRd
minimum_calibration_n <- function(spec) {
  if (spec$n_within_cells > 1L) spec$n_within_cells + 1L else 2L
}


#' @keywords internal
#' @noRd
validate_calibration_n <- function(n, spec, arg) {
  min_n <- minimum_calibration_n(spec)
  too_small <- n < min_n
  if (!any(too_small)) return(invisible(n))

  bad <- paste(sort(unique(n[too_small])), collapse = ", ")
  stop(
    "`", arg, "` is too small for this design. This design has ",
    spec$n_within_cells, " within-subject cell",
    if (spec$n_within_cells == 1L) "" else "s",
    ", so empirical calibration needs at least ", min_n,
    " participant", if (min_n == 1L) "" else "s",
    " per between-subject cell. ",
    "Increase `", arg, "`",
    if (length(n) > 1L) " values" else "",
    " (too small: ", bad, ").",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
compound_symmetric_sigma <- function(n, sd, r) {
  sigma <- matrix(r * sd^2, nrow = n, ncol = n)
  diag(sigma) <- sd^2
  sigma
}


#' @keywords internal
#' @noRd
fit_design_model <- function(data, spec) {
  data <- set_sum_contrasts(data, spec)
  formula <- design_aov_formula(spec)
  suppressWarnings(stats::aov(formula, data = data))
}


#' @keywords internal
#' @noRd
fit_design_term <- function(data, spec, term, alpha, ss_type = "III") {
  stats <- tryCatch(
    fit_design_term_stats(data, spec, term, ss_type = ss_type),
    error = function(e) NULL
  )
  if (is.null(stats) || !is.finite(stats$p_value)) return(NA)
  isTRUE(stats$p_value < alpha)
}


#' @keywords internal
#' @noRd
fit_design_term_stats <- function(data, spec, term, ss_type = "III") {
  ss_type <- validate_ss_type(ss_type)
  if (identical(ss_type, "I")) {
    fit <- fit_design_model(data, spec)
    return(extract_term_stats(fit, term))
  }
  fit_car_term_stats(data, spec, term, ss_type = ss_type)
}


#' @keywords internal
#' @noRd
fit_car_term_stats <- function(data, spec, term, ss_type) {
  data <- set_sum_contrasts(data, spec)
  type <- switch(ss_type, II = 2, III = 3)
  gg_p_values <- NULL
  if (!length(spec$within)) {
    fixed <- paste(spec$factor_names, collapse = " * ")
    formula <- stats::as.formula(paste("value ~", fixed))
    fit <- stats::lm(formula, data = data)
    tab <- as.data.frame(car::Anova(fit, type = type))
    rows <- car_between_rows(tab)
  } else {
    wide_setup <- repeated_measures_wide_data(data, spec)
    fit <- stats::lm(wide_setup$formula, data = wide_setup$wide)
    idesign <- stats::reformulate(paste(spec$within, collapse = " * "))
    av <- car::Anova(
      fit,
      idata = wide_setup$idata,
      idesign = idesign,
      type = type,
      icontrasts = c("contr.sum", "contr.poly")
    )
    av_summary <- suppressWarnings(summary(av, multivariate = FALSE))
    tab <- as.data.frame.matrix(av_summary$univariate.tests)
    rows <- car_repeated_rows(tab)
    gg_p_values <- car_gg_p_values(av_summary)
  }

  idx <- which(rows$term == term)
  if (!length(idx)) {
    stop("Term '", term, "' was not found in the fitted ANOVA table.",
         call. = FALSE)
  }
  row <- rows[idx[[1L]], , drop = FALSE]
  if (!is.finite(row$f_value) || !is.finite(row$den_df)) {
    stop("Term '", term, "' does not have a finite F test in the fitted ANOVA.",
         call. = FALSE)
  }
  pes <- row$f_value * row$num_df / (row$f_value * row$num_df + row$den_df)
  p_value_gg <- if (!is.null(gg_p_values) && term %in% names(gg_p_values)) {
    as.numeric(gg_p_values[[term]])
  } else {
    NA_real_
  }
  list(
    num_df = as.numeric(row$num_df),
    den_df = as.numeric(row$den_df),
    pes = as.numeric(pes),
    p_value = as.numeric(row$p_value),
    p_value_gg = p_value_gg
  )
}


#' Extract Greenhouse--Geisser-adjusted p-values from a car::Anova summary
#'
#' @keywords internal
#' @noRd
car_gg_p_values <- function(av_summary) {
  adj <- av_summary$pval.adjustments
  if (is.null(adj) || !nrow(adj)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  terms <- trimws(rownames(adj))
  values <- as.numeric(adj[, "Pr(>F[GG])"])
  stats::setNames(values, terms)
}


#' @keywords internal
#' @noRd
car_between_rows <- function(tab) {
  row_terms <- trimws(rownames(tab))
  residual_idx <- which(row_terms == "Residuals")
  den_df <- if (length(residual_idx)) {
    as.numeric(tab[residual_idx[[1L]], "Df"])
  } else {
    NA_real_
  }
  keep <- row_terms != "Residuals" & row_terms != "(Intercept)" & nzchar(row_terms)
  tibble::tibble(
    term = row_terms[keep],
    num_df = as.numeric(tab[keep, "Df"]),
    den_df = den_df,
    f_value = as.numeric(tab[keep, "F value"]),
    p_value = as.numeric(tab[keep, "Pr(>F)"])
  )
}


#' @keywords internal
#' @noRd
car_repeated_rows <- function(tab) {
  row_terms <- trimws(rownames(tab))
  keep <- row_terms != "(Intercept)" & nzchar(row_terms)
  tibble::tibble(
    term = row_terms[keep],
    num_df = as.numeric(tab[keep, "num Df"]),
    den_df = as.numeric(tab[keep, "den Df"]),
    f_value = as.numeric(tab[keep, "F value"]),
    p_value = as.numeric(tab[keep, "Pr(>F)"])
  )
}


#' @keywords internal
#' @noRd
repeated_measures_wide_data <- function(data, spec) {
  within_labels <- make_cell_labels("w", spec$n_within_cells)
  within_keys <- interaction_key(spec$within_cells, spec$within)
  data$.within_key <- interaction_key(data, spec$within)
  data$.within_cell <- within_labels[match(data$.within_key, within_keys)]
  data$.within_key <- NULL

  id_cols <- c("id", spec$between)
  wide <- data |>
    dplyr::select(dplyr::all_of(c(id_cols, ".within_cell", "value"))) |>
    tidyr::pivot_wider(
      names_from = ".within_cell",
      values_from = "value"
    )
  for (nm in spec$between) {
    wide[[nm]] <- factor(wide[[nm]], levels = levels(spec$levels[[nm]]))
    stats::contrasts(wide[[nm]]) <- stats::contr.sum(nlevels(wide[[nm]]))
  }

  response <- paste0("cbind(", paste(within_labels, collapse = ", "), ")")
  rhs <- if (length(spec$between)) {
    paste(spec$between, collapse = " * ")
  } else {
    "1"
  }
  idata <- spec$within_cells[, spec$within, drop = FALSE]
  for (nm in spec$within) {
    idata[[nm]] <- factor(idata[[nm]], levels = levels(spec$levels[[nm]]))
    stats::contrasts(idata[[nm]]) <- stats::contr.sum(nlevels(idata[[nm]]))
  }

  list(
    wide = wide,
    idata = as.data.frame(idata),
    formula = stats::as.formula(paste(response, "~", rhs))
  )
}


#' @keywords internal
#' @noRd
extract_term_stats <- function(fit, term) {
  rows <- extract_aov_rows(fit)
  idx <- which(rows$term == term)
  if (!length(idx)) {
    stop("Term '", term, "' was not found in the fitted ANOVA table.",
         call. = FALSE)
  }
  row <- rows[idx[[1L]], , drop = FALSE]
  if (!is.finite(row$f_value) || !is.finite(row$den_df)) {
    stop("Term '", term, "' does not have a finite F test in the fitted ANOVA.",
         call. = FALSE)
  }
  pes <- row$f_value * row$num_df / (row$f_value * row$num_df + row$den_df)
  list(
    num_df = as.numeric(row$num_df),
    den_df = as.numeric(row$den_df),
    pes = as.numeric(pes),
    p_value = as.numeric(row$p_value),
    p_value_gg = NA_real_
  )
}


#' @keywords internal
#' @noRd
set_sum_contrasts <- function(data, spec) {
  for (nm in spec$factor_names) {
    if (nm %in% names(data)) {
      data[[nm]] <- factor(data[[nm]], levels = levels(spec$levels[[nm]]))
      stats::contrasts(data[[nm]]) <- stats::contr.sum(nlevels(data[[nm]]))
    }
  }
  data$id <- factor(data$id)
  data
}


#' @keywords internal
#' @noRd
design_aov_formula <- function(spec) {
  fixed <- paste(spec$factor_names, collapse = " * ")
  if (!length(spec$within)) {
    return(stats::as.formula(paste("value ~", fixed)))
  }
  within <- paste(spec$within, collapse = " * ")
  stats::as.formula(
    paste("value ~", fixed, "+ Error(id / (", within, "))")
  )
}


#' @keywords internal
#' @noRd
extract_aov_rows <- function(fit) {
  s <- summary(fit)
  if (inherits(fit, "aovlist")) {
    tables <- lapply(s, `[[`, 1L)
  } else {
    tables <- list(s[[1L]])
  }
  rows <- vector("list", length(tables))
  for (i in seq_along(tables)) {
    tab <- as.data.frame(tables[[i]])
    row_terms <- trimws(rownames(tab))
    residual_idx <- which(row_terms == "Residuals")
    den_df <- if (length(residual_idx)) {
      as.numeric(tab[residual_idx[[1L]], "Df"])
    } else {
      NA_real_
    }
    keep <- row_terms != "Residuals" & nzchar(row_terms)
    rows[[i]] <- tibble::tibble(
      term = row_terms[keep],
      num_df = as.numeric(tab[keep, "Df"]),
      den_df = den_df,
      f_value = as.numeric(tab[keep, "F value"]),
      p_value = as.numeric(tab[keep, "Pr(>F)"])
    )
  }
  dplyr::bind_rows(rows)
}


#' @keywords internal
#' @noRd
make_progress_bar <- function(enabled, total, label) {
  if (!isTRUE(enabled)) return(NULL)
  total <- max(1L, as.integer(total))
  bar <- utils::txtProgressBar(min = 0, max = total, style = 3)
  env <- new.env(parent = emptyenv())
  env$bar <- bar
  env$total <- total
  env$value <- 0L
  env$label <- label
  env
}


#' @keywords internal
#' @noRd
tick_progress_bar <- function(progress_bar, by = 1L) {
  if (is.null(progress_bar)) return(invisible(NULL))
  progress_bar$value <- min(progress_bar$total, progress_bar$value + by)
  utils::setTxtProgressBar(progress_bar$bar, progress_bar$value)
  invisible(NULL)
}


#' @keywords internal
#' @noRd
close_progress_bar <- function(progress_bar) {
  if (is.null(progress_bar)) return(invisible(NULL))
  utils::setTxtProgressBar(progress_bar$bar, progress_bar$total)
  close(progress_bar$bar)
  invisible(NULL)
}


#' @keywords internal
#' @noRd
estimate_adaptive_progress_total <- function(n_start, n_max, max_iter = 25L,
                                             n_min = 1L) {
  n <- max(1L, as.integer(n_start))
  probes_lower_bound <- n > as.integer(n_min)
  n_max <- as.integer(n_max)
  steps <- 1L
  while (n < n_max) {
    n_next <- min(n_max, max(n + 1L, n * 2L))
    steps <- steps + 1L
    n <- n_next
  }
  steps + as.integer(probes_lower_bound) + max_iter
}
