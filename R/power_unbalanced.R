#' Define cells for a means-based unbalanced ANOVA design
#'
#' Creates the complete cell table used by [power_unbalanced()]. Each cell is
#' defined by its factor levels, sample size (`n`), population mean (`m`), and
#' population standard deviation (`sd`). End each cell after all three reserved
#' values have been supplied.
#'
#' @section Lifecycle:
#' \ifelse{html}{\out{<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg" alt="[Experimental]"></a>}}{\strong{Experimental}}
#'
#' `cell_design()` is experimental and is available only in the development
#' version of `anovapowersim`. Its API may change.
#'
#' @param ... Repeated named cell definitions. Each cell must contain the same
#'   factor names in the same order, plus `n`, `m`, and `sd`. Every factor must
#'   have at least 2 observed levels, and every combination of factor levels
#'   must appear exactly once (or be filled automatically; see `default_n`).
#' @param within Character vector naming factors in `...` that are measured
#'   within subjects, or `NULL` for a purely between-subject design. Stored on
#'   the returned object and read by [power_unbalanced()].
#' @param default_n,default_m,default_sd Optional scalars used to fill any
#'   missing cells in the complete factorial design. Supply all three to
#'   auto-fill missing cells with these values; supply none to require every
#'   cell to be defined explicitly (the default). Supplying only some of the
#'   three is an error.
#'
#' @return An `anovapowersim_cell_design` tibble with one row per design cell.
#'
#' @examples
#' design <- cell_design(
#'   group = "control", time = "pre",  n = 22, m = 10.0, sd = 2.0,
#'   group = "control", time = "post", n = 22, m = 11.0, sd = 2.2,
#'   group = "treatment", time = "pre",  n = 31, m = 10.1, sd = 2.4,
#'   group = "treatment", time = "post", n = 31, m = 12.4, sd = 2.8,
#'   within = "time"
#' )
#' design
#'
#' @export
cell_design <- function(...,
                        within = NULL,
                        default_n = NULL,
                        default_m = NULL,
                        default_sd = NULL) {
  dots <- list(...)
  nms <- names(dots)
  reserved <- c("n", "m", "sd")
  if (!length(dots)) {
    stop("Enter at least one cell with factor values, `n`, `m`, and `sd`.",
         call. = FALSE)
  }
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("Every value in `...` must be named.", call. = FALSE)
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
  if (!is.numeric(out$sd) || any(!is.finite(out$sd)) || any(out$sd <= 0)) {
    stop("`sd` must contain positive finite standard deviations.",
         call. = FALSE)
  }

  key <- interaction_key(out, factor_names)
  if (anyDuplicated(key)) {
    stop("Each design cell may be defined only once.", call. = FALSE)
  }

  within <- validate_unbalanced_within(within, factor_names)

  default_n <- validate_optional_default_n(default_n)
  default_m <- validate_optional_default_m(default_m)
  default_sd <- validate_optional_default_sd(default_sd)
  defaults_supplied <- c(
    !is.null(default_n), !is.null(default_m), !is.null(default_sd)
  )
  if (any(defaults_supplied) && !all(defaults_supplied)) {
    stop(
      "Supply all of `default_n`, `default_m`, and `default_sd` to fill ",
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
    if (use_defaults) {
      fill <- expected[missing_idx, , drop = FALSE]
      fill$n <- default_n
      fill$m <- default_m
      fill$sd <- default_sd
      out <- dplyr::bind_rows(
        out, fill[, c(factor_names, reserved), drop = FALSE]
      )
    } else {
      missing_cells <- expected[missing_idx, , drop = FALSE]
      labels <- vapply(seq_len(nrow(missing_cells)), function(i) {
        paste(
          paste0(factor_names, ' = "',
                 unlist(missing_cells[i, factor_names]), '"'),
          collapse = ", "
        )
      }, character(1))
      stop(
        "Every combination of factor levels must be defined exactly once. ",
        "Missing cell", if (length(labels) == 1L) "" else "s", ":\n",
        paste0("  ", labels, collapse = "\n"),
        "\nSupply the missing cell", if (length(labels) == 1L) "" else "s",
        ", or provide `default_n`, `default_m`, and `default_sd` to fill ",
        "missing cells automatically.",
        call. = FALSE
      )
    }
  }

  grid <- resolve_unbalanced_grid(out, factor_names, within)
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


#' @keywords internal
#' @noRd
validate_optional_default_sd <- function(default_sd) {
  if (is.null(default_sd)) return(NULL)
  if (!is.numeric(default_sd) || length(default_sd) != 1L ||
      !is.finite(default_sd) || default_sd <= 0) {
    stop("`default_sd` must be a single positive finite number.", call. = FALSE)
  }
  as.numeric(default_sd)
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


#' Specify correlations for a means-based unbalanced design
#'
#' Defines the common within-subject correlation structure used by
#' [power_unbalanced()]. Marginal standard deviations are deliberately absent:
#' they come from [cell_design()] and may differ between groups.
#'
#' @section Lifecycle:
#' \ifelse{html}{\out{<a href="https://lifecycle.r-lib.org/articles/stages.html#experimental"><img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg" alt="[Experimental]"></a>}}{\strong{Experimental}}
#'
#' `unbalanced_covariance()` is experimental and is available only in the
#' development version of `anovapowersim`. Its API may change.
#'
#' @param default_correlation Correlation in `(-1, 1)` used for unlisted pairs.
#' @param correlations Optional named numeric vector of pair-specific
#'   correlations. Name pairs as `"cell1:cell2"`; pair order does not matter.
#'
#' @return An `anovapowersim_unbalanced_covariance_spec` object.
#'
#' @examples
#' unbalanced_covariance(
#'   default_correlation = 0.5,
#'   correlations = c("pre:post" = 0.7)
#' )
#'
#' @export
unbalanced_covariance <- function(default_correlation = 0.5,
                                  correlations = NULL) {
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
  structure(
    list(
      default_correlation = as.numeric(default_correlation),
      correlations = correlations,
      correlation_pairs = parse_correlation_pairs(names(correlations))
    ),
    class = "anovapowersim_unbalanced_covariance_spec"
  )
}


#' Simulate power for a fixed unbalanced ANOVA design
#'
#' Estimates achieved power for exact, potentially unequal cell sizes and
#' user-supplied population means and standard deviations. This function is
#' simulation-only: it does not calculate power from a noncentral F
#' distribution and does not scale the supplied sample sizes.
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
#' @param covariance Optional correlation specification created by
#'   [unbalanced_covariance()]. It can only be supplied when `design` has
#'   within-subject factors. The default `NULL` uses a correlation of `0.5`
#'   between within-subject
#'   measurements. Standard deviations always come from `design` and may
#'   differ between between-subject cells. If the term's population
#'   covariance is non-spherical for any cell, `power_sim` is based on each
#'   simulated dataset's Greenhouse--Geisser-corrected p-value rather than the
#'   uncorrected univariate test. This correction requires `ss_type` `"III"`
#'   or `"II"`; under `"I"`, simulated p-values remain uncorrected, and a
#'   warning is issued.
#' @param n_sims Number of simulated datasets.
#' @param alpha Significance threshold.
#' @param ss_type Sums-of-squares type: `"III"`, `"II"`, or `"I"`.
#' @param progress Logical; if `TRUE`, show a text progress bar.
#' @param parallel Logical; if `TRUE`, run simulations in parallel.
#' @param cores Optional positive integer number of parallel workers.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An `anovapowersim_unbalanced_power` object. `$power` and
#'   `$achieved_power` contain simulated power. `$partial_eta_squared` is the
#'   term effect size in a deterministic reference dataset. `$epsilon` is the
#'   worst-case (smallest) population Greenhouse--Geisser epsilon across
#'   between-subject cells for the tested term. `$results` also reports the
#'   simulated partial eta-squared distribution and failed fits.
#'
#' @examples
#' \donttest{
#' design <- cell_design(
#'   group = "control", time = "pre",  n = 12, m = 10, sd = 2.0,
#'   group = "control", time = "post", n = 12, m = 11, sd = 2.2,
#'   group = "treatment", time = "pre",  n = 18, m = 10, sd = 2.4,
#'   group = "treatment", time = "post", n = 18, m = 13, sd = 2.8,
#'   within = "time"
#' )
#' power_unbalanced(
#'   design = design,
#'   term = "group:time",
#'   covariance = unbalanced_covariance(
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
                             seed = NULL) {
  spec <- prepare_unbalanced_means_design(design)
  term <- resolve_design_term(term, spec)
  ss_type <- validate_ss_type(ss_type)
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
  spec$within_correlation <- resolve_unbalanced_correlation(covariance, spec)
  validate_calibration_n(spec$cell_n, spec, "design$n")
  epsilon <- unbalanced_term_epsilon(spec, term)
  warn_ss_type_i_uncorrected_gg(ss_type = ss_type, epsilon = epsilon)
  message_long_serial_run(as.integer(n_sims), parallel)

  reference <- simulate_unbalanced_means_data(spec, empirical = TRUE)
  reference_stats <- fit_design_term_stats(
    reference, spec, term, ss_type = ss_type
  )
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
    epsilon = epsilon,
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
      epsilon = epsilon,
      design = design,
      design_spec = spec,
      covariance = covariance,
      correlation = spec$within_correlation,
      custom_covariance = !is.null(covariance),
      ss_type = ss_type,
      call = match.call()
    ),
    class = "anovapowersim_unbalanced_power"
  )
}


#' @keywords internal
#' @noRd
prepare_unbalanced_means_design <- function(design) {
  if (!inherits(design, "anovapowersim_cell_design")) {
    stop("`design` must be created by cell_design().", call. = FALSE)
  }
  factor_names <- setdiff(names(design), c("n", "m", "sd"))
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
    cell_sds = matrix(
      grid$ordered$sd, nrow = grid$n_between, ncol = grid$n_within,
      byrow = TRUE
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
resolve_unbalanced_correlation <- function(covariance, spec) {
  if (!length(spec$within)) {
    if (!is.null(covariance)) {
      stop(
        "`covariance` can only be supplied when `within` identifies at ",
        "least one within-subject factor.",
        call. = FALSE
      )
    }
    return(matrix(1, nrow = 1L, ncol = 1L,
                  dimnames = list("outcome", "outcome")))
  }
  if (is.null(covariance)) covariance <- unbalanced_covariance()
  if (!inherits(covariance, "anovapowersim_unbalanced_covariance_spec")) {
    stop(
      "`covariance` must be created by unbalanced_covariance(); ",
      "within_covariance() objects and raw matrices are not accepted.",
      call. = FALSE
    )
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
  correlation
}


#' Population Greenhouse-Geisser epsilon for an unbalanced means design
#'
#' Between-subject cells may have different standard deviations, so there is
#' no single shared covariance matrix for the term. Each cell's own implied
#' covariance (its standard deviations combined with the shared correlation)
#' is checked, and the worst-case (smallest) epsilon across cells is used, so
#' that a Greenhouse-Geisser-corrected simulated p-value is used whenever any
#' cell's population covariance is non-spherical for this term.
#'
#' @keywords internal
#' @noRd
unbalanced_term_epsilon <- function(spec, term) {
  if (!length(spec$within)) return(1)
  epsilons <- vapply(seq_len(spec$n_between_cells), function(i) {
    covariance <- outer(spec$cell_sds[i, ], spec$cell_sds[i, ]) *
      spec$within_correlation
    covariance_term_epsilon(covariance = covariance, spec = spec, term = term)
  }, numeric(1))
  min(epsilons)
}


#' @keywords internal
#' @noRd
simulate_unbalanced_means_data <- function(spec, empirical = FALSE) {
  between_labels <- make_cell_labels("b", spec$n_between_cells)
  within_labels <- make_cell_labels("w", spec$n_within_cells)
  subject_rows <- vector("list", spec$n_between_cells)
  y_rows <- vector("list", spec$n_between_cells)
  id_offset <- 0L

  for (i in seq_len(spec$n_between_cells)) {
    n_i <- spec$cell_n[[i]]
    standard_deviations <- spec$cell_sds[i, ]
    sigma <- outer(standard_deviations, standard_deviations) *
      spec$within_correlation
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
                                             ss_type, epsilon = 1,
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
    "tick_progress_bar"
  ))
  simulate_one_dataset <- helpers$simulate_unbalanced_means_data
  fit_one_stats <- helpers$fit_design_term_stats
  tick_progress <- helpers$tick_progress_bar
  use_gg_correction <- isTRUE(epsilon < 1 - 1e-8)
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
    p_value <- if (use_gg_correction && is.finite(stats$p_value_gg)) {
      stats$p_value_gg
    } else {
      stats$p_value
    }
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
    simulated_power = sprintf("%.3f", object$power),
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
  cat("\nCell design:\n")
  print(object$design, row.names = FALSE)
  cat("\nPower and effect-size diagnostics:\n")
  print(format_power_results(object$results), row.names = FALSE)
  invisible(out)
}
