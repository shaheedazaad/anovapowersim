#' Specify a within-subject covariance structure
#'
#' Creates a readable covariance specification for [power_n()] and
#' [power_curve()]. Supply one common standard deviation and default
#' correlation, then override only the individual correlation pairs that
#' differ. Measurement names are generated from the `within` design; for example,
#' `within = c(time = 2, condition = 2)` creates `time1_condition1`,
#' `time1_condition2`, `time2_condition1`, and `time2_condition2`.
#'
#' @param sd Common positive finite standard deviation for every measurement.
#'   If omitted, `1` is used and a warning is issued.
#' @param default_correlation Finite correlation in `(-1, 1)` used for pairs
#'   not listed in `correlations`. When the specification is resolved for a
#'   design, a warning identifies how many pairs were not defined and makes
#'   clear that this default applies only to those pairs.
#' @param correlations Optional named numeric vector of pair-specific
#'   correlations. Name each pair as `"cell1:cell2"`, for example
#'   `"time1:time2" = 0.7`. Pair order does not matter.
#'
#' @return An `anovapowersim_covariance_spec` object. [power_n()] and
#'   [power_curve()] resolve it against their `within` design and construct the
#'   covariance matrix as `sd^2 * R`, where `R` is the correlation matrix.
#'   They also derive the tested term's
#'   population Greenhouse--Geisser epsilon from the resolved matrix and apply
#'   it to `power_calc`.
#'
#' @examples
#' covariance <- within_covariance(
#'   sd = 1,
#'   default_correlation = 0.5,
#'   correlations = c(
#'     "time1_condition1:time1_condition2" = 0.6,
#'     "time2_condition1:time2_condition2" = 0.7
#'   )
#' )
#'
#' \donttest{
#' power_n(
#'   within = c(time = 3, condition = 2),
#'   term = "time:condition",
#'   target_pes = 0.14,
#'   power = 0.90,
#'   n_sims = 100,
#'   covariance = covariance,
#'   seed = 123
#' )
#' }
#'
#' @export
within_covariance <- function(sd = 1,
                              default_correlation = 0.5,
                              correlations = NULL) {
  sd_missing <- missing(sd)
  if (!is.numeric(sd) || length(sd) != 1L ||
      !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a single positive finite number.",
         call. = FALSE)
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
      "`sd` was not supplied to `within_covariance()`; using one common ",
      "`sd = 1` for every measurement.",
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
    class = "anovapowersim_covariance_spec"
  )
}


#' @keywords internal
#' @noRd
validate_named_covariance_values <- function(x, arg, predicate, requirement) {
  if (is.null(x)) return(stats::setNames(numeric(0), character(0)))
  if (!is.numeric(x) || any(!is.finite(x)) || any(!predicate(x))) {
    stop("`", arg, "` must contain only ", requirement, ".",
         call. = FALSE)
  }
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    stop("Every value in `", arg, "` must be named.", call. = FALSE)
  }
  if (anyDuplicated(names(x))) {
    stop("Names in `", arg, "` must be unique.", call. = FALSE)
  }
  x
}


#' @keywords internal
#' @noRd
parse_correlation_pairs <- function(pair_names) {
  if (!length(pair_names)) {
    return(tibble::tibble(
      pair_name = character(0),
      cell1 = character(0),
      cell2 = character(0),
      pair_key = character(0)
    ))
  }

  pieces <- strsplit(pair_names, ":", fixed = TRUE)
  valid <- lengths(pieces) == 2L & vapply(
    pieces,
    function(x) all(nzchar(x)),
    logical(1)
  )
  if (!all(valid)) {
    stop(
      "Names in `correlations` must identify exactly two within-subject ",
      "cells as 'cell1:cell2'. Problem name",
      if (sum(!valid) == 1L) "" else "s", ": ",
      paste(shQuote(pair_names[!valid]), collapse = ", "), ".",
      call. = FALSE
    )
  }

  cell1 <- vapply(pieces, `[[`, character(1), 1L)
  cell2 <- vapply(pieces, `[[`, character(1), 2L)
  self <- cell1 == cell2
  if (any(self)) {
    stop(
      "`correlations` must not define a cell's correlation with itself. ",
      "Problem name", if (sum(self) == 1L) "" else "s", ": ",
      paste(shQuote(pair_names[self]), collapse = ", "), ".",
      call. = FALSE
    )
  }

  pair_key <- vapply(
    seq_along(cell1),
    function(i) paste(sort(c(cell1[[i]], cell2[[i]])), collapse = "\r"),
    character(1)
  )
  duplicated_pair <- duplicated(pair_key) | duplicated(pair_key, fromLast = TRUE)
  if (any(duplicated_pair)) {
    stop(
      "Each pair in `correlations` may be defined only once, regardless of ",
      "order. Duplicated definitions: ",
      paste(shQuote(pair_names[duplicated_pair]), collapse = ", "), ".",
      call. = FALSE
    )
  }

  tibble::tibble(
    pair_name = pair_names,
    cell1 = cell1,
    cell2 = cell2,
    pair_key = pair_key
  )
}


#' @keywords internal
#' @noRd
within_cell_names <- function(spec) {
  if (!length(spec$within)) return("outcome")
  values <- lapply(spec$within_cells[spec$within], as.character)
  do.call(paste, c(values, sep = "_"))
}


#' @keywords internal
#' @noRd
resolve_within_covariance <- function(covariance, spec, sd = 1, r = 0.5) {
  cell_names <- within_cell_names(spec)
  n_cells <- length(cell_names)

  if (is.null(covariance)) {
    if (length(spec$within)) {
      warning(
        "No `covariance` was supplied; using one common `sd = ", sd,
        "` and `correlation = ", r,
        "` for every within-subject pair.",
        call. = FALSE,
        immediate. = TRUE
      )
    } else {
      warning(
        "No SD or `covariance` was supplied; using one common `sd = ", sd,
        "`. Correlations do not apply to this purely between-subject design.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    sigma <- compound_symmetric_sigma(n_cells, sd = sd, r = r)
    dimnames(sigma) <- list(cell_names, cell_names)
    return(sigma)
  }
  if (!length(spec$within)) {
    stop("`covariance` can only be supplied for a design with a `within` factor.",
         call. = FALSE)
  }

  if (inherits(covariance, "anovapowersim_covariance_spec")) {
    return(resolve_within_covariance_spec(covariance, cell_names))
  }
  validate_direct_covariance_matrix(covariance, cell_names)
}


#' @keywords internal
#' @noRd
resolve_within_covariance_spec <- function(covariance, cell_names) {
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
      value <- covariance$correlations[[i]]
      correlation[cell1, cell2] <- value
      correlation[cell2, cell1] <- value
    }
  }

  sigma <- covariance$sd^2 * correlation
  dimnames(sigma) <- list(cell_names, cell_names)
  sigma <- validate_positive_definite_covariance(sigma)
  warn_undefined_correlations(
    defined_pairs = nrow(covariance$correlation_pairs),
    n_cells = length(cell_names),
    default_correlation = covariance$default_correlation
  )
  sigma
}


#' Warn when the default correlation fills unnamed within-subject pairs
#'
#' @keywords internal
#' @noRd
warn_undefined_correlations <- function(defined_pairs, n_cells,
                                        default_correlation) {
  total_pairs <- choose(n_cells, 2L)
  undefined_pairs <- total_pairs - defined_pairs
  if (undefined_pairs <= 0L) return(invisible(NULL))

  warning(
    undefined_pairs, " of ", total_pairs,
    " within-subject correlation",
    if (total_pairs == 1L) " was" else
      if (undefined_pairs == 1L) "s was" else "s were",
    " not explicitly defined; `default_correlation = ",
    default_correlation, "` is used only for ",
    if (undefined_pairs == 1L) "that undefined pair" else
      "those undefined pairs",
    ". Explicitly defined correlations are unchanged.",
    call. = FALSE,
    immediate. = TRUE
  )
  invisible(NULL)
}


#' @keywords internal
#' @noRd
validate_direct_covariance_matrix <- function(covariance, cell_names) {
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      any(!is.finite(covariance))) {
    stop(
      "`covariance` must be created by within_covariance() or be a finite ",
      "numeric matrix.",
      call. = FALSE
    )
  }
  expected_dim <- c(length(cell_names), length(cell_names))
  if (!identical(dim(covariance), expected_dim)) {
    stop(
      "`covariance` must be a ", expected_dim[[1]], " x ", expected_dim[[2]],
      " matrix for this design.",
      call. = FALSE
    )
  }
  row_names <- rownames(covariance)
  col_names <- colnames(covariance)
  has_row_names <- !is.null(row_names)
  has_col_names <- !is.null(col_names)
  if (xor(has_row_names, has_col_names)) {
    stop("`covariance` must have both row and column names, or neither.",
         call. = FALSE)
  }
  if (has_row_names) {
    if (anyDuplicated(row_names) || anyDuplicated(col_names) ||
        !setequal(row_names, cell_names) || !setequal(col_names, cell_names)) {
      stop(
        "Row and column names in `covariance` must match the expected cells: ",
        paste(shQuote(cell_names), collapse = ", "), ".",
        call. = FALSE
      )
    }
    covariance <- covariance[cell_names, cell_names, drop = FALSE]
  } else {
    dimnames(covariance) <- list(cell_names, cell_names)
  }
  if (!isTRUE(all.equal(covariance, t(covariance), tolerance = 1e-10,
                        check.attributes = FALSE))) {
    stop("`covariance` must be symmetric.", call. = FALSE)
  }

  validate_equal_covariance_diagonal(covariance)
  validate_positive_definite_covariance(covariance)
}


#' Validate that a covariance matrix has one common marginal variance
#'
#' Equal marginal variances do not imply sphericity: unequal correlations can
#' still yield unequal variances of pairwise differences.
#'
#' @keywords internal
#' @noRd
validate_equal_covariance_diagonal <- function(covariance,
                                               tolerance = 1e-10) {
  variances <- diag(covariance)
  common <- rep(variances[[1L]], length(variances))
  if (!isTRUE(all.equal(
    variances, common, tolerance = tolerance, check.attributes = FALSE
  ))) {
    stop(
      "`covariance` must have equal diagonal variances. anovapowersim ",
      "requires one common marginal standard deviation; unequal ",
      "correlations are still allowed.",
      call. = FALSE
    )
  }
  invisible(covariance)
}


#' @keywords internal
#' @noRd
validate_positive_definite_covariance <- function(covariance) {
  if (any(diag(covariance) <= 0)) {
    stop("`covariance` must have positive variances on its diagonal.",
         call. = FALSE)
  }
  positive_definite <- tryCatch(
    {
      chol(covariance)
      TRUE
    },
    error = function(e) FALSE
  )
  if (!positive_definite) {
    stop(
      "The resolved `covariance` matrix must be positive definite. Check ",
      "the supplied standard deviations and correlations.",
      call. = FALSE
    )
  }
  covariance
}


#' Degrees of freedom for the within-subject component of a term
#'
#' Purely between-subject terms have a within component of one for the purpose
#' of shared epsilon and direction-sensitivity checks.
#'
#' @keywords internal
#' @noRd
within_term_df <- function(spec, term) {
  term_factors <- strsplit(term, ":", fixed = TRUE)[[1L]]
  within_factors <- intersect(term_factors, spec$within)
  if (!length(within_factors)) return(1L)
  as.integer(prod(spec$level_counts[within_factors] - 1L))
}


#' Calculate population Greenhouse--Geisser epsilon for a design term
#'
#' @keywords internal
#' @noRd
covariance_term_epsilon <- function(covariance, spec, term) {
  term_factors <- strsplit(term, ":", fixed = TRUE)[[1L]]
  within_term_factors <- intersect(term_factors, spec$within)
  if (!length(within_term_factors)) return(1)

  term_df <- within_term_df(spec = spec, term = term)
  if (term_df == 1L) return(1)

  contrast_matrix <- within_term_contrast_matrix(
    spec = spec,
    within_term_factors = within_term_factors
  )
  projected <- contrast_matrix %*% covariance %*% t(contrast_matrix)
  trace_projected <- sum(diag(projected))
  epsilon <- trace_projected^2 /
    (term_df * sum(projected * projected))

  # Floating-point error can put a theoretically valid epsilon just outside
  # its bounds. Clamp only after the population value has been calculated.
  lower_bound <- 1 / term_df
  epsilon <- min(1, max(lower_bound, epsilon))
  if (abs(epsilon - 1) < 1e-12) epsilon <- 1
  as.numeric(epsilon)
}


#' Build orthonormal contrasts for the within component of a term
#'
#' @keywords internal
#' @noRd
within_term_contrast_matrix <- function(spec, within_term_factors) {
  factor_matrices <- lapply(spec$within, function(factor_name) {
    n_levels <- spec$level_counts[[factor_name]]
    if (!factor_name %in% within_term_factors) {
      return(matrix(rep(1 / sqrt(n_levels), n_levels), nrow = 1L))
    }

    contrasts <- stats::contr.helmert(n_levels)
    contrasts <- sweep(
      contrasts,
      MARGIN = 2L,
      STATS = sqrt(colSums(contrasts^2)),
      FUN = "/"
    )
    t(contrasts)
  })

  Reduce(function(left, right) kronecker(left, right), factor_matrices)
}
