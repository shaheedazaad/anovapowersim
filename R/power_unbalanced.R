#' Define cells for a means-based unbalanced ANOVA design
#'
#' Creates the complete cell table used by [power_unbalanced()]. Each cell is
#' defined by its factor levels, sample size (`n`), and population mean (`m`).
#' End each cell after both reserved values have been supplied. The common
#' population standard deviation belongs in [unbalanced_covariance()].
#'
#' The `m` values are literal population cell means: their magnitudes and all
#' effects they contain are used as supplied. This differs from
#' [means_pattern()], whose values specify only a relative shape that balanced
#' simulation functions project onto the tested term, normalize, and rescale
#' to `target_pes`.
#'
#' @section Lifecycle:
#' \ifelse{html}{\out{<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg" alt="[Experimental]"></a>}}{\strong{Experimental}}
#'
#' `cell_design()` is experimental and is available only in the development
#' version of `anovapowersim`. Its API may change.
#'
#' @param ... Repeated named cell definitions. Each cell must contain the same
#'   factor names in the same order, plus `n` and `m`. Every factor must
#'   have at least 2 observed levels, and every combination of factor levels
#'   must appear exactly once (or be filled automatically; see `default_n`).
#' @param within Character vector naming factors in `...` that are measured
#'   within subjects, or `NULL` for a purely between-subject design. Stored on
#'   the returned object and read by [power_unbalanced()]. Within-cell names
#'   used by [unbalanced_covariance()] join level values with `_`; these names
#'   must be unique, and within-factor levels must not contain `:`.
#' @param default_n,default_m Optional scalars used to fill any missing cells
#'   in the complete factorial design. Supply both to
#'   auto-fill missing cells with these values; supply none to require every
#'   cell to be defined explicitly (the default). Supplying only one is an
#'   error. When cells are auto-filled, a message reports their count and exact
#'   factor-level combinations so that unintended levels can be spotted.
#'
#' @return An `anovapowersim_cell_design` tibble with one row per design cell.
#'
#' @seealso [means_pattern()] for shape-only patterns used by balanced
#'   simulation functions.
#'
#' @examples
#' design <- cell_design(
#'   group = "control", time = "pre",  n = 22, m = 10.0,
#'   group = "control", time = "post", n = 22, m = 11.0,
#'   group = "treatment", time = "pre",  n = 31, m = 10.1,
#'   group = "treatment", time = "post", n = 31, m = 12.4,
#'   within = "time"
#' )
#' design
#'
#' @export
cell_design <- function(...,
                        within = NULL,
                        default_n = NULL,
                        default_m = NULL) {
  dots <- list(...)
  nms <- names(dots)
  reserved <- c("n", "m")
  if (!length(dots)) {
    stop("Enter at least one cell with factor values, `n`, and `m`.",
         call. = FALSE)
  }
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("Every value in `...` must be named.", call. = FALSE)
  }
  legacy_sd <- intersect(nms, c("sd", "default_sd"))
  if (length(legacy_sd)) {
    stop(
      "Standard deviations no longer belong in `cell_design()`. Remove ",
      paste(paste0("`", unique(legacy_sd), "`"), collapse = " and "),
      " and supply one common SD with `unbalanced_covariance(sd = ...)`.",
      call. = FALSE
    )
  }

  rows <- list()
  current <- list()
  factor_names <- NULL
  for (i in seq_along(dots)) {
    nm <- nms[[i]]
    value <- dots[[i]]
    if (length(value) != 1L || is.null(value) || is.na(value)) {
      stop("`", nm, "` values must be single non-missing values.",
           call. = FALSE)
    }
    if (nm %in% names(current)) {
      stop("`", nm, "` appears more than once in the same cell.",
           call. = FALSE)
    }
    current[[nm]] <- value

    if (all(reserved %in% names(current))) {
      row_factor_names <- names(current)[!names(current) %in% reserved]
      if (!length(row_factor_names)) {
        stop("Each cell must define at least one factor.", call. = FALSE)
      }
      if (is.null(factor_names)) {
        factor_names <- row_factor_names
      } else if (!identical(row_factor_names, factor_names)) {
        stop("Every cell must use the same factor names in the same order.",
             call. = FALSE)
      }
      current <- current[c(factor_names, reserved)]
      rows[[length(rows) + 1L]] <- current
      current <- list()
    }
  }
  if (length(current)) {
    missing <- setdiff(reserved, names(current))
    stop(
      "The final cell is incomplete; missing ",
      paste(paste0("`", missing, "`"), collapse = ", "), ".",
      call. = FALSE
    )
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
  rows <- lapply(rows, function(row) {
    for (nm in factor_names) row[[nm]] <- as.character(row[[nm]])
    row
  })
  out <- dplyr::bind_rows(rows)
  out$n <- validate_cell_count_values(out$n, "n")
  if (!is.numeric(out$m) || any(!is.finite(out$m))) {
    stop("`m` must contain finite numeric means.", call. = FALSE)
  }
  key <- interaction_key(out, factor_names)
  if (anyDuplicated(key)) {
    stop("Each design cell may be defined only once.", call. = FALSE)
  }

  within <- validate_unbalanced_within(within, factor_names)

  default_n <- validate_optional_default_n(default_n)
  default_m <- validate_optional_default_m(default_m)
  defaults_supplied <- c(!is.null(default_n), !is.null(default_m))
  if (any(defaults_supplied) && !all(defaults_supplied)) {
    stop(
      "Supply both `default_n` and `default_m` to fill ",
      "missing cells automatically, or none of them.",
      call. = FALSE
    )
  }
  use_defaults <- all(defaults_supplied)

  factor_levels <- stats::setNames(
    lapply(factor_names, function(nm) unique(out[[nm]])),
    factor_names
  )
  for (nm in factor_names) {
    lv <- factor_levels[[nm]]
    if (length(lv) < 2L) {
      stop(
        "Factor `", nm, "` must have at least 2 levels; only ",
        paste(shQuote(lv), collapse = ", "),
        " was supplied. Add cells for every level of `", nm,
        "`, or remove it from the design if it is not meant to vary.",
        call. = FALSE
      )
    }
  }

  expected <- tidyr::expand_grid(!!!factor_levels)
  expected_key <- interaction_key(expected, factor_names)
  missing_idx <- which(!expected_key %in% key)

  if (length(missing_idx)) {
    missing_cells <- expected[missing_idx, , drop = FALSE]
    labels <- vapply(seq_len(nrow(missing_cells)), function(i) {
      paste(
        paste0(factor_names, ' = "',
               unlist(missing_cells[i, factor_names]), '"'),
        collapse = ", "
      )
    }, character(1))
    if (use_defaults) {
      message(
        "Auto-filled ", length(labels), " missing cell",
        if (length(labels) == 1L) "" else "s",
        " using `default_n = ", default_n, "` and `default_m = ",
        format(default_m, trim = TRUE), "`:\n",
        paste0("  ", labels, collapse = "\n")
      )
      fill <- missing_cells
      fill$n <- default_n
      fill$m <- default_m
      out <- dplyr::bind_rows(
        out, fill[, c(factor_names, reserved), drop = FALSE]
      )
    } else {
      stop(
        "Every combination of factor levels must be defined exactly once. ",
        "Missing cell", if (length(labels) == 1L) "" else "s", ":\n",
        paste0("  ", labels, collapse = "\n"),
        "\nSupply the missing cell", if (length(labels) == 1L) "" else "s",
        ", or provide `default_n` and `default_m` to fill ",
        "missing cells automatically.",
        call. = FALSE
      )
    }
  }

  grid <- resolve_unbalanced_grid(out, factor_names, within)
  within_cell_names(list(
    within = grid$within,
    within_cells = grid$within_cells
  ))
  n_matrix <- matrix(
    grid$ordered$n, nrow = grid$n_between, ncol = grid$n_within, byrow = TRUE
  )
  inconsistent <- apply(n_matrix, 1L, function(x) length(unique(x)) != 1L)
  if (any(inconsistent)) {
    bad_rows <- which(inconsistent)
    bad_labels <- if (length(grid$between)) {
      vapply(bad_rows, function(r) {
        paste(
          paste0(grid$between, ' = "',
                 as.character(unlist(grid$between_cells[r, ])), '"'),
          collapse = ", "
        )
      }, character(1))
    } else {
      "the design"
    }
    stop(
      "`n` must be identical across all within-subject rows belonging to ",
      "the same between-subject cell. Inconsistent group",
      if (length(bad_labels) == 1L) "" else "s", ": ",
      paste(bad_labels, collapse = "; "), ".",
      call. = FALSE
    )
  }

  out <- out[match(expected_key, interaction_key(out, factor_names)), ,
             drop = FALSE]
  out <- tibble::as_tibble(out[, c(factor_names, reserved), drop = FALSE])
  class(out) <- c("anovapowersim_cell_design", class(out))
  attr(out, "within") <- grid$within
  out
}


#' @keywords internal
#' @noRd
validate_unbalanced_within <- function(within, factor_names) {
  if (is.null(within)) within <- character(0)
  if (!is.character(within) || anyNA(within) || any(!nzchar(within)) ||
      anyDuplicated(within)) {
    stop("`within` must be NULL or a character vector of unique factor names.",
         call. = FALSE)
  }
  unknown <- setdiff(within, factor_names)
  if (length(unknown)) {
    stop("Unknown factor", if (length(unknown) == 1L) "" else "s",
         " in `within`: ", paste(shQuote(unknown), collapse = ", "), ".",
         call. = FALSE)
  }
  within
}


#' @keywords internal
#' @noRd
validate_optional_default_n <- function(default_n) {
  if (is.null(default_n)) return(NULL)
  if (!is.numeric(default_n) || length(default_n) != 1L ||
      !is.finite(default_n) || default_n < 1 ||
      default_n != as.integer(default_n)) {
    stop("`default_n` must be a single positive integer.", call. = FALSE)
  }
  as.integer(default_n)
}


#' @keywords internal
#' @noRd
validate_optional_default_m <- function(default_m) {
  if (is.null(default_m)) return(NULL)
  if (!is.numeric(default_m) || length(default_m) != 1L ||
      !is.finite(default_m)) {
    stop("`default_m` must be a single finite number.", call. = FALSE)
  }
  as.numeric(default_m)
}


#' Partition factors into a between x within grid for an unbalanced design
#'
#' Shared by [cell_design()] (to check `n`-consistency at construction time)
#' and `prepare_unbalanced_means_design()` (to build the design spec), so the
#' between/within partitioning and row ordering are always derived identically.
#'
#' @keywords internal
#' @noRd
resolve_unbalanced_grid <- function(data, factor_names, within) {
  between <- factor_names[!factor_names %in% within]
  within <- factor_names[factor_names %in% within]
  level_values <- stats::setNames(
    lapply(factor_names, function(nm) unique(as.character(data[[nm]]))),
    factor_names
  )
  levels <- stats::setNames(
    lapply(factor_names, function(nm) {
      factor(level_values[[nm]], levels = level_values[[nm]])
    }),
    factor_names
  )
  between_cells <- if (length(between)) {
    tidyr::expand_grid(!!!levels[between])
  } else {
    tibble::tibble(.dummy_between = factor("all"))
  }
  within_cells <- if (length(within)) {
    tidyr::expand_grid(!!!levels[within])
  } else {
    tibble::tibble(.dummy_within = factor("dv"))
  }
  ordered_grid <- if (length(between) && length(within)) {
    tidyr::expand_grid(!!!levels[c(between, within)])
  } else if (length(between)) {
    tidyr::expand_grid(!!!levels[between])
  } else {
    tidyr::expand_grid(!!!levels[within])
  }
  design_key <- interaction_key(data, c(between, within))
  grid_key <- interaction_key(ordered_grid, c(between, within))
  ordered <- data[match(grid_key, design_key), , drop = FALSE]
  n_between <- if (length(between)) nrow(between_cells) else 1L
  n_within <- if (length(within)) nrow(within_cells) else 1L

  list(
    between = between,
    within = within,
    levels = levels,
    between_cells = between_cells,
    within_cells = within_cells,
    n_between = n_between,
    n_within = n_within,
    ordered = ordered
  )
}


#' Specify covariance for a means-based unbalanced design
#'
#' Defines the common marginal standard deviation and within-subject
#' correlation structure used by [power_unbalanced()]. The resulting
#' covariance is shared by every between-subject cell.
#'
#' @section Lifecycle:
#' \ifelse{html}{\out{<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg" alt="[Experimental]"></a>}}{\strong{Experimental}}
#'
#' `unbalanced_covariance()` is experimental and is available only in the
#' development version of `anovapowersim`. Its API may change.
#'
#' @param sd Common positive finite marginal standard deviation. If omitted,
#'   `1` is used and a warning is issued.
#' @param default_correlation Correlation in `(-1, 1)` used for unlisted pairs.
#'   When the specification is resolved for a design, a warning identifies how
#'   many pairs were not defined and makes clear that this default applies only
#'   to those pairs.
#' @param correlations Optional named numeric vector of pair-specific
#'   correlations. Name pairs as `"cell1:cell2"`; pair order does not matter.
#'   For multiple within factors, cell names join their level values with `_`.
#'   Constructed names must be unique, and level values must not contain `:`
#'   because it separates the two cells in a pair name.
#'
#' @return An `anovapowersim_unbalanced_covariance_spec` object.
#'
#' @examples
#' unbalanced_covariance(
#'   sd = 2,
#'   default_correlation = 0.5,
#'   correlations = c("pre:post" = 0.7)
#' )
#'
#' @export
unbalanced_covariance <- function(sd = 1,
                                  default_correlation = 0.5,
                                  correlations = NULL) {
  sd_missing <- missing(sd)
  if (!is.numeric(sd) || length(sd) != 1L ||
      !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a single positive finite number.", call. = FALSE)
  }
  if (!is.numeric(default_correlation) ||
      length(default_correlation) != 1L ||
      !is.finite(default_correlation) ||
      default_correlation <= -1 || default_correlation >= 1) {
    stop("`default_correlation` must be a single finite number in (-1, 1).",
         call. = FALSE)
  }
  correlations <- validate_named_covariance_values(
    correlations,
    arg = "correlations",
    predicate = function(x) x > -1 & x < 1,
    requirement = "finite numbers in (-1, 1)"
  )
  pairs <- parse_correlation_pairs(names(correlations))
  if (sd_missing) {
    warning(
      "`sd` was not supplied to `unbalanced_covariance()`; using one common ",
      "`sd = 1` for every design cell.",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  structure(
    list(
      sd = as.numeric(sd),
      default_correlation = as.numeric(default_correlation),
      correlations = correlations,
      correlation_pairs = pairs
    ),
    class = "anovapowersim_unbalanced_covariance_spec"
  )
}


#' Simulate power for a fixed unbalanced ANOVA design
#'
#' Estimates achieved power for exact, potentially unequal cell sizes and
#' user-supplied population means under one common standard deviation. This
#' function is simulation-only: it does not calculate power from a noncentral
#' F distribution and does not scale the supplied sample sizes. A warning is
#' issued when the deterministic reference data imply essentially zero effect
#' for the tested term, which often indicates a typo, a wrong `term`, or means
#' that contain only other effects.
#'
#' @section Lifecycle:
#' \ifelse{html}{\out{<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg" alt="[Experimental]"></a>}}{\strong{Experimental}}
#'
#' `power_unbalanced()` is experimental and is available only in the
#' development version of `anovapowersim`. Its API and reporting format may
#' change.
#'
#' @param design A complete design table created by [cell_design()].
#'   Within-subject factors are read from the design's `within` attribute
#'   (set via [cell_design()]'s `within` argument), not supplied here.
#' @param term Character scalar naming the ANOVA term to test.
#' @param covariance Optional common covariance specification created by
#'   [unbalanced_covariance()]. The default `NULL` uses `sd = 1` and a
#'   correlation of `0.5` between within-subject measurements and issues a
#'   warning stating those defaults. For purely between-subject designs, use
#'   this argument to change the common SD; correlation settings are not
#'   applicable. When an [unbalanced_covariance()] specification omits some
#'   correlation pairs, a warning states that `default_correlation` is used
#'   only for those undefined pairs. The resolved covariance determines the
#'   term-specific population Greenhouse--Geisser epsilon used by
#'   `sim_correction = "auto"`.
#' @param n_sims Number of simulated datasets.
#' @param alpha Significance threshold.
#' @param ss_type Sums-of-squares type: `"III"`, `"II"`, or `"I"`.
#'   For unequal-N designs, `"I"` uses sequential, order-dependent hypotheses;
#'   a warning reports the factor order inherited from [cell_design()].
#' @param sim_correction Sphericity correction for simulated p-values:
#'   `"auto"` (the default) uses Greenhouse--Geisser correction when the
#'   term-specific population epsilon is below `1 - 1e-8` and `ss_type` is
#'   `"II"` or `"III"`; `"GG"` requests correction for every simulated
#'   dataset; and `"none"` always uses the uncorrected univariate test.
#'   `"GG"` is an error with `ss_type = "I"`. For a between-only term or a
#'   within component with one degree of freedom, `"GG"` silently resolves to
#'   `"none"` because no sphericity correction applies.
#' @param progress Logical; if `TRUE`, show a text progress bar.
#' @param parallel Logical; if `TRUE`, run simulations in parallel.
#' @param cores Optional positive integer number of parallel workers.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An `anovapowersim_unbalanced_power` object. `$power` and
#'   `$achieved_power` contain simulated power. `$partial_eta_squared` is the
#'   term effect size in a deterministic reference dataset. `$epsilon` is the
#'   population Greenhouse--Geisser epsilon for the tested term. `$results`
#'   also reports the common SD, simulated partial eta-squared distribution,
#'   and failed fits. Sample partial eta squared is upward-biased in finite
#'   samples, so `mean_pes_sim`, `median_pes_sim`, and the simulated interval
#'   are sampling diagnostics rather than estimates of the supplied population
#'   effect; use `$partial_eta_squared` as the reference effect. The object
#'   stores the requested `sim_correction` and applied
#'   `sim_correction_resolved` values.
#'
#' @section Simulated sphericity correction:
#' `sim_correction` governs only the simulated test. With `"GG"`, each dataset
#' uses its own sample-estimated Greenhouse--Geisser epsilon from
#' `car::Anova()`. Under a truly spherical population, that sample correction
#' is mildly conservative. Power is estimated for the prespecified corrected
#' or uncorrected test; conditional Mauchly-then-correct procedures are not
#' simulated. Unlike the balanced simulation functions, `power_unbalanced()`
#' does not report a `power_calc` diagnostic.
#'
#' @examples
#' \donttest{
#' design <- cell_design(
#'   group = "control", time = "pre",  n = 12, m = 10,
#'   group = "control", time = "post", n = 12, m = 11,
#'   group = "treatment", time = "pre",  n = 18, m = 10,
#'   group = "treatment", time = "post", n = 18, m = 13,
#'   within = "time"
#' )
#' power_unbalanced(
#'   design = design,
#'   term = "group:time",
#'   covariance = unbalanced_covariance(
#'     sd = 2,
#'     correlations = c("pre:post" = 0.7)
#'   ),
#'   n_sims = 100,
#'   seed = 123
#' )
#' }
#'
#' @export
power_unbalanced <- function(design,
                             term,
                             covariance = NULL,
                             n_sims = 10000,
                             alpha = 0.05,
                             ss_type = "III",
                             progress = interactive(),
                             parallel = FALSE,
                             cores = NULL,
                             seed = NULL,
                             sim_correction = c("auto", "GG", "none")) {
  spec <- prepare_unbalanced_means_design(design)
  term <- resolve_design_term(term, spec)
  ss_type <- validate_ss_type(ss_type)
  warn_unbalanced_type_i_order(ss_type = ss_type, spec = spec)
  assert_unit_interval(alpha, "alpha")
  if (!is.numeric(n_sims) || length(n_sims) != 1L || !is.finite(n_sims) ||
      n_sims < 1 || n_sims != as.integer(n_sims)) {
    stop("`n_sims` must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  cores <- validate_parallel_cores(cores, parallel)
  custom_covariance <- !is.null(covariance)
  if (is.null(covariance)) {
    if (length(spec$within)) {
      warning(
        "No `covariance` was supplied; using one common `sd = 1` and ",
        "`correlation = 0.5` for every within-subject pair.",
        call. = FALSE,
        immediate. = TRUE
      )
    } else {
      warning(
        "No SD or `covariance` was supplied; using one common `sd = 1`. ",
        "Correlations do not apply to this purely between-subject design.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    covariance <- unbalanced_covariance(sd = 1)
  }
  if (!inherits(covariance, "anovapowersim_unbalanced_covariance_spec")) {
    stop(
      "`covariance` must be created by unbalanced_covariance(); ",
      "within_covariance() objects and raw matrices are not accepted.",
      call. = FALSE
    )
  }
  spec$sd <- covariance$sd
  spec$within_correlation <- resolve_unbalanced_correlation(
    covariance, spec, warn_undefined = custom_covariance
  )
  validate_calibration_n(spec$cell_n, spec, "design$n")
  epsilon <- unbalanced_term_epsilon(spec, term)
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
  message_long_serial_run(as.integer(n_sims), parallel)

  reference <- simulate_unbalanced_means_data(spec, empirical = TRUE)
  reference_stats <- fit_design_term_stats(
    reference, spec, term, ss_type = ss_type
  )
  use_gg_correction <- identical(correction$resolved, "GG")
  assert_reference_gg_p(
    stats = reference_stats,
    use_gg_correction = use_gg_correction,
    term = term
  )
  warn_negligible_reference_effect(reference_stats$pes, term = term)
  if (!is.null(seed)) set.seed(seed)

  progress_bar <- make_progress_bar(
    enabled = progress,
    total = if (parallel) 1L else as.integer(n_sims),
    label = "Simulating unbalanced power"
  )
  on.exit(close_progress_bar(progress_bar), add = TRUE)

  simulation <- run_unbalanced_means_simulations(
    spec = spec,
    term = term,
    n_sims = as.integer(n_sims),
    alpha = alpha,
    ss_type = ss_type,
    sim_correction_resolved = correction$resolved,
    progress_bar = if (parallel) NULL else progress_bar,
    parallel = parallel,
    cores = cores
  )
  if (parallel) tick_progress_bar(progress_bar)

  total_n <- as.integer(sum(spec$cell_n))
  pes_interval <- stats::quantile(
    simulation$pes,
    probs = c(0.025, 0.975),
    names = FALSE,
    na.rm = TRUE
  )
  row <- tibble::tibble(
    total_n = total_n,
    min_cell_n = as.integer(min(spec$cell_n)),
    max_cell_n = as.integer(max(spec$cell_n)),
    n_sims = as.integer(n_sims),
    valid_sims = simulation$valid_count,
    failed_sims = simulation$failed_count,
    sd = spec$sd,
    epsilon = epsilon,
    num_df = reference_stats$num_df,
    den_df = reference_stats$den_df,
    partial_eta_squared = reference_stats$pes,
    mean_pes_sim = mean(simulation$pes),
    median_pes_sim = stats::median(simulation$pes),
    pes_sim_lower = pes_interval[[1L]],
    pes_sim_upper = pes_interval[[2L]],
    power_sim = simulation$power
  )

  structure(
    list(
      results = row,
      term = term,
      alpha = alpha,
      n_sims = as.integer(n_sims),
      power = simulation$power,
      achieved_power = simulation$power,
      partial_eta_squared = reference_stats$pes,
      mean_pes_sim = row$mean_pes_sim[[1L]],
      median_pes_sim = row$median_pes_sim[[1L]],
      pes_sim_interval = unname(pes_interval),
      total_n = total_n,
      min_cell_n = as.integer(min(spec$cell_n)),
      max_cell_n = as.integer(max(spec$cell_n)),
      valid_sims = simulation$valid_count,
      failed_sims = simulation$failed_count,
      sd = spec$sd,
      epsilon = epsilon,
      design = design,
      design_spec = spec,
      covariance = covariance,
      correlation = spec$within_correlation,
      custom_covariance = custom_covariance,
      ss_type = ss_type,
      sim_correction = correction$requested,
      sim_correction_resolved = correction$resolved,
      call = match.call()
    ),
    class = "anovapowersim_unbalanced_power"
  )
}


#' Warn that sequential sums of squares depend on factor order under unequal N
#'
#' @keywords internal
#' @noRd
warn_unbalanced_type_i_order <- function(ss_type, spec) {
  if (!identical(ss_type, "I") || length(unique(spec$cell_n)) <= 1L) {
    return(invisible(NULL))
  }

  warning(
    "`ss_type = \"I\"` uses sequential sums of squares in this unequal-N ",
    "design, so the tested hypothesis depends on factor order. The factor ",
    "order inherited from `cell_design()` is: ",
    paste(spec$factor_names, collapse = ", "),
    ". Use `ss_type = \"II\"` or `\"III\"` unless this sequential, ",
    "order-dependent hypothesis is intentional.",
    call. = FALSE,
    immediate. = TRUE
  )
  invisible(NULL)
}


#' Warn when literal cell means contain essentially no tested-term effect
#'
#' @keywords internal
#' @noRd
warn_negligible_reference_effect <- function(pes, term, threshold = 1e-8) {
  if (!is.finite(pes) || pes > threshold) return(invisible(NULL))

  warning(
    "The supplied `m` values imply essentially no effect for tested term '",
    term, "' (reference partial eta squared = ", format(pes, digits = 3),
    ", at or below ", format(threshold, scientific = TRUE), "). Check for ",
    "typos, confirm `term`, and verify that the cell means contain the ",
    "tested effect rather than only other main effects or interactions.",
    call. = FALSE,
    immediate. = TRUE
  )
  invisible(NULL)
}


#' @keywords internal
#' @noRd
prepare_unbalanced_means_design <- function(design) {
  if (!inherits(design, "anovapowersim_cell_design")) {
    stop("`design` must be created by cell_design().", call. = FALSE)
  }
  factor_names <- setdiff(names(design), c("n", "m"))
  within <- attr(design, "within")
  if (is.null(within)) within <- character(0)

  grid <- resolve_unbalanced_grid(design, factor_names, within)
  n_matrix <- matrix(
    grid$ordered$n, nrow = grid$n_between, ncol = grid$n_within, byrow = TRUE
  )

  spec <- list(
    between = grid$between,
    within = grid$within,
    factor_names = c(grid$between, grid$within),
    level_counts = stats::setNames(
      vapply(grid$levels[c(grid$between, grid$within)], nlevels, integer(1)),
      c(grid$between, grid$within)
    ),
    levels = grid$levels[c(grid$between, grid$within)],
    between_cells = grid$between_cells,
    within_cells = grid$within_cells,
    n_between_cells = grid$n_between,
    n_within_cells = grid$n_within,
    cell_n = as.integer(n_matrix[, 1L]),
    cell_means = matrix(
      grid$ordered$m, nrow = grid$n_between, ncol = grid$n_within, byrow = TRUE
    ),
    cell_design = grid$ordered
  )
  class(spec) <- c(
    "anovapowersim_unbalanced_means_design_spec",
    "anovapowersim_design_spec",
    class(spec)
  )
  spec
}


#' @keywords internal
#' @noRd
resolve_unbalanced_correlation <- function(covariance, spec,
                                           warn_undefined = TRUE) {
  if (!length(spec$within)) {
    if (length(covariance$correlations)) {
      stop(
        "`correlations` cannot be supplied for a purely between-subject ",
        "design.",
        call. = FALSE
      )
    }
    return(matrix(1, nrow = 1L, ncol = 1L,
                  dimnames = list("outcome", "outcome")))
  }
  cell_names <- within_cell_names(spec)
  named_cells <- unique(c(
    covariance$correlation_pairs$cell1,
    covariance$correlation_pairs$cell2
  ))
  unknown <- setdiff(named_cells, cell_names)
  if (length(unknown)) {
    stop(
      "Unknown within-subject cell", if (length(unknown) == 1L) "" else "s",
      " in `covariance`: ", paste(shQuote(unknown), collapse = ", "),
      ". Expected cells: ", paste(shQuote(cell_names), collapse = ", "), ".",
      call. = FALSE
    )
  }
  correlation <- matrix(
    covariance$default_correlation,
    nrow = length(cell_names),
    ncol = length(cell_names),
    dimnames = list(cell_names, cell_names)
  )
  diag(correlation) <- 1
  if (nrow(covariance$correlation_pairs)) {
    for (i in seq_len(nrow(covariance$correlation_pairs))) {
      cell1 <- covariance$correlation_pairs$cell1[[i]]
      cell2 <- covariance$correlation_pairs$cell2[[i]]
      correlation[cell1, cell2] <- covariance$correlations[[i]]
      correlation[cell2, cell1] <- covariance$correlations[[i]]
    }
  }
  positive_definite <- tryCatch({ chol(correlation); TRUE },
                                error = function(e) FALSE)
  if (!positive_definite) {
    stop(
      "The resolved within-subject correlation matrix must be positive ",
      "definite. Check the supplied correlations.",
      call. = FALSE
    )
  }
  if (warn_undefined) {
    warn_undefined_correlations(
      defined_pairs = nrow(covariance$correlation_pairs),
      n_cells = length(cell_names),
      default_correlation = covariance$default_correlation
    )
  }
  correlation
}


#' Population Greenhouse-Geisser epsilon for an unbalanced means design
#'
#' Every between-subject cell shares one covariance matrix, constructed from
#' the common SD and within-subject correlation matrix.
#'
#' @keywords internal
#' @noRd
unbalanced_term_epsilon <- function(spec, term) {
  if (!length(spec$within)) return(1)
  covariance <- spec$sd^2 * spec$within_correlation
  covariance_term_epsilon(covariance = covariance, spec = spec, term = term)
}


#' @keywords internal
#' @noRd
simulate_unbalanced_means_data <- function(spec, empirical = FALSE) {
  between_labels <- make_cell_labels("b", spec$n_between_cells)
  within_labels <- make_cell_labels("w", spec$n_within_cells)
  subject_rows <- vector("list", spec$n_between_cells)
  y_rows <- vector("list", spec$n_between_cells)
  id_offset <- 0L
  sigma <- spec$sd^2 * spec$within_correlation

  for (i in seq_len(spec$n_between_cells)) {
    n_i <- spec$cell_n[[i]]
    y <- MASS::mvrnorm(
      n = n_i,
      mu = spec$cell_means[i, ],
      Sigma = sigma,
      empirical = empirical
    )
    if (spec$n_within_cells == 1L) {
      y <- matrix(y, nrow = n_i, ncol = 1L)
    } else if (n_i == 1L) {
      y <- matrix(y, nrow = 1L)
    }
    colnames(y) <- within_labels
    b <- spec$between_cells[i, spec$between, drop = FALSE]
    subject_rows[[i]] <- dplyr::bind_cols(
      tibble::tibble(id = seq_len(n_i) + id_offset),
      b[rep(1L, n_i), , drop = FALSE]
    )
    y_rows[[i]] <- y
    id_offset <- id_offset + n_i
  }
  subjects <- dplyr::bind_rows(subject_rows)
  y_mat <- do.call(rbind, y_rows)
  wide <- dplyr::bind_cols(subjects, tibble::as_tibble(y_mat))
  if (!length(spec$within)) {
    out <- wide |>
      dplyr::rename(value = dplyr::all_of(within_labels[[1L]]))
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
run_unbalanced_means_simulations <- function(spec, term, n_sims, alpha,
                                             ss_type,
                                             sim_correction_resolved,
                                             progress_bar = NULL,
                                             parallel = FALSE, cores = NULL) {
  helpers <- parallel_worker_helpers(c(
    "simulate_unbalanced_means_data",
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
    "select_simulated_p",
    "tick_progress_bar"
  ))
  simulate_one_dataset <- helpers$simulate_unbalanced_means_data
  fit_one_stats <- helpers$fit_design_term_stats
  select_p <- helpers$select_simulated_p
  tick_progress <- helpers$tick_progress_bar
  use_gg_correction <- identical(sim_correction_resolved, "GG")
  simulate_one <- function(i) {
    sim <- simulate_one_dataset(spec, empirical = FALSE)
    tick_progress(progress_bar)
    stats <- tryCatch(
      fit_one_stats(sim, spec, term, ss_type = ss_type),
      error = function(e) NULL
    )
    if (is.null(stats) || !is.finite(stats$pes)) {
      return(list(reject = NA, pes = NA_real_))
    }
    p_value <- select_p(stats, use_gg_correction)
    if (!is.finite(p_value)) return(list(reject = NA, pes = NA_real_))
    list(reject = isTRUE(p_value < alpha), pes = stats$pes)
  }
  simulations <- if (parallel) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
    future.apply::future_lapply(
      seq_len(n_sims), simulate_one, future.seed = TRUE
    )
  } else {
    purrr::map(seq_len(n_sims), simulate_one)
  }
  rejected <- vapply(simulations, `[[`, logical(1L), "reject")
  pes <- vapply(simulations, `[[`, numeric(1L), "pes")
  failed_count <- sum(is.na(rejected) | is.na(pes))
  valid <- !is.na(rejected) & !is.na(pes)
  valid_count <- sum(valid)
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
      " simulated ANOVA fits failed and were excluded from the results.",
      call. = FALSE
    )
  }
  list(
    power = mean(rejected[valid]),
    pes = pes[valid],
    valid_count = as.integer(valid_count),
    failed_count = as.integer(failed_count)
  )
}


#' Print simulated power for an unbalanced design
#'
#' @param x An `anovapowersim_unbalanced_power` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.anovapowersim_unbalanced_power <- function(x, ...) {
  cat("<anovapowersim_unbalanced_power>\n")
  cat("  term:                  '", x$term, "'\n", sep = "")
  cat("  fixed total N:         ", x$total_n, "\n", sep = "")
  cat("  between-cell n range:  ", x$min_cell_n, " to ",
      x$max_cell_n, "\n", sep = "")
  cat("  simulations:           ", x$n_sims, "\n", sep = "")
  cat("  simulated test:        ", simulated_test_label(
    x$sim_correction, x$sim_correction_resolved
  ), "\n", sep = "")
  cat("  common SD:             ", format(x$sd), "\n", sep = "")
  cat("  simulated power:       ", sprintf("%.3f", x$power), "\n",
      sep = "")
  if (!is.null(x$epsilon) && is.finite(x$epsilon) && x$epsilon < 1) {
    cat("  epsilon:               ", format(x$epsilon), "\n", sep = "")
  }
  cat("  reference pes:         ",
      sprintf("%.4f", x$partial_eta_squared), "\n", sep = "")
  cat("  mean simulated pes:    ", sprintf("%.4f", x$mean_pes_sim),
      "\n", sep = "")
  cat("  median simulated pes:  ", sprintf("%.4f", x$median_pes_sim),
      "\n", sep = "")
  cat("  simulated pes 95%:     [",
      sprintf("%.4f", x$pes_sim_interval[[1L]]), ", ",
      sprintf("%.4f", x$pes_sim_interval[[2L]]), "]\n", sep = "")
  cat("  pes note:              sample pes is upward-biased; simulated pes ",
      "summaries are diagnostics, not the population/reference effect.\n",
      sep = "")
  correlation_label <- if (!length(x$design_spec$within)) {
    "not applicable"
  } else if (isTRUE(x$custom_covariance)) {
    "custom"
  } else {
    "default (0.5)"
  }
  cat("  within correlations:   ", correlation_label, "\n", sep = "")
  cat("  SS type:               ", x$ss_type, "\n", sep = "")
  if (x$failed_sims > 0L) {
    cat("  failed fits:           ", x$failed_sims, "\n", sep = "")
  }
  cat("\nCell design:\n")
  print(x$design, row.names = FALSE)
  cat("\nPower and effect-size diagnostics:\n")
  print(format_power_results(x$results), row.names = FALSE)
  invisible(x)
}


#' Summarise simulated power for an unbalanced design
#'
#' @param object An `anovapowersim_unbalanced_power` object.
#' @param ... Unused.
#'
#' @return A list with `header`, `design`, and `results`, invisibly; printed to
#'   the console as well.
#' @export
summary.anovapowersim_unbalanced_power <- function(object, ...) {
  header <- c(
    term = object$term,
    fixed_total_n = as.character(object$total_n),
    between_cell_n_range = paste(object$min_cell_n, object$max_cell_n,
                                 sep = " to "),
    alpha = sprintf("%.3f", object$alpha),
    simulations = as.character(object$n_sims),
    valid_simulations = as.character(object$valid_sims),
    failed_simulations = as.character(object$failed_sims),
    common_sd = format(object$sd),
    simulated_power = sprintf("%.3f", object$power),
    simulated_test = simulated_test_label(
      object$sim_correction, object$sim_correction_resolved
    ),
    reference_partial_eta_squared = sprintf(
      "%.4f", object$partial_eta_squared
    ),
    mean_simulated_partial_eta_squared = sprintf("%.4f", object$mean_pes_sim),
    median_simulated_partial_eta_squared = sprintf(
      "%.4f", object$median_pes_sim
    ),
    simulated_pes_95_interval = sprintf(
      "[%.4f, %.4f]", object$pes_sim_interval[[1L]],
      object$pes_sim_interval[[2L]]
    ),
    within_correlations = if (!length(object$design_spec$within)) {
      "not applicable"
    } else if (isTRUE(object$custom_covariance)) {
      "custom"
    } else {
      "default (0.5)"
    },
    ss_type = object$ss_type
  )
  if (!is.null(object$epsilon) && is.finite(object$epsilon) &&
      object$epsilon < 1) {
    header <- append(header, c(epsilon = sprintf("%.4f", object$epsilon)),
                     after = 8L)
  }
  out <- list(header = header, design = object$design, results = object$results)
  cat("anovapowersim unbalanced-design power summary\n")
  cat("---------------------------------------------\n")
  for (nm in names(header)) {
    cat(sprintf("  %-38s %s\n", paste0(nm, ":"), header[[nm]]))
  }
  cat("  note:                                  sample pes is upward-biased; ",
      "simulated pes summaries are diagnostics, not the population/",
      "reference effect.\n", sep = "")
  cat("\nCell design:\n")
  print(object$design, row.names = FALSE)
  cat("\nPower and effect-size diagnostics:\n")
  print(format_power_results(object$results), row.names = FALSE)
  invisible(out)
}
