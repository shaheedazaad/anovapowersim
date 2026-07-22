#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Assert that a number is in the open interval (0, 1)
#'
#' @param x Numeric.
#' @param arg Argument name used in the error message.
#' @keywords internal
#' @noRd
assert_unit_interval <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    stop("`", arg, "` must be a single number in (0, 1). Got: ",
         paste(x, collapse = ", "), call. = FALSE)
  }
  invisible(x)
}


#' Coerce and assert that a value is in the open interval (0, 1)
#'
#' Accepts numeric scalars and numeric-looking character scalars. This keeps
#' low-level helpers usable with formatted effect-size output such as ".310"
#' while preserving strict interval validation.
#'
#' @param x Numeric or character scalar.
#' @param arg Argument name used in the error message.
#' @keywords internal
#' @noRd
as_unit_interval <- function(x, arg) {
  x_num <- if (is.character(x) && length(x) == 1L) {
    suppressWarnings(as.numeric(trimws(x)))
  } else {
    x
  }
  assert_unit_interval(x_num, arg)
  x_num
}


#' Warn when a sample-size search does not reach target power
#'
#' @keywords internal
#' @noRd
warn_target_power_not_reached <- function(n_needed, target, n_max) {
  if (!is.na(n_needed)) return(invisible(NULL))
  warning(
    sprintf(
      paste(
        "Target power %.3f was not reached by `n_max = %d`.",
        "Increase `n_max` and rerun the search."
      ),
      target,
      as.integer(n_max)
    ),
    call. = FALSE
  )
  invisible(NULL)
}


#' Wilson confidence interval for a binomial proportion
#'
#' @keywords internal
#' @noRd
binomial_wilson_interval <- function(successes, trials, conf_level = 0.95) {
  if (trials <= 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  proportion <- successes / trials
  denominator <- 1 + z^2 / trials
  centre <- (proportion + z^2 / (2 * trials)) / denominator
  half_width <- z * sqrt(
    proportion * (1 - proportion) / trials + z^2 / (4 * trials^2)
  ) / denominator

  c(
    lower = max(0, centre - half_width),
    upper = min(1, centre + half_width)
  )
}


#' Create internal, collision-proof cell labels
#'
#' @keywords internal
#' @noRd
make_cell_labels <- function(prefix, n) {
  sprintf(".anovapowersim_%s_%03d", prefix, seq_len(n))
}
