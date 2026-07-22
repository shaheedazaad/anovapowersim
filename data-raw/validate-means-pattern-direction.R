# Independent development-time audit of mean-direction sensitivity.
#
# This deliberately computes the one-factor repeated-measures F statistic
# directly from multivariate-normal draws instead of calling the package's
# ANOVA-fitting helpers. It is not part of the automated test suite because the
# high simulation count is intended for release audits rather than every check.

library(MASS)

direct_rm_power <- function(direction, covariance, n, target_pes,
                            n_sims = 20000, alpha = 0.05,
                            gg_corrected = TRUE, seed = 1) {
  n_levels <- length(direction)
  contrast_projection <- diag(n_levels) - 1 / n_levels
  contrast_covariance <- contrast_projection %*% covariance %*%
    contrast_projection
  term_df <- n_levels - 1L
  epsilon <- sum(diag(contrast_covariance))^2 /
    (term_df * sum(contrast_covariance^2))

  direction <- as.numeric(contrast_projection %*% direction)
  direction <- direction / sqrt(sum(direction^2))
  target_norm_squared <- target_pes / (1 - target_pes) *
    (n - 1) / n * sum(diag(contrast_covariance))
  population_mean <- direction * sqrt(target_norm_squared)

  numerator_df <- if (gg_corrected) epsilon * term_df else term_df
  denominator_df <- if (gg_corrected) {
    epsilon * term_df * (n - 1)
  } else {
    term_df * (n - 1)
  }
  critical_f <- stats::qf(
    1 - alpha, df1 = numerator_df, df2 = denominator_df
  )

  set.seed(seed)
  rejected <- replicate(n_sims, {
    outcomes <- MASS::mvrnorm(n, mu = population_mean, Sigma = covariance)
    sample_mean <- colMeans(outcomes)
    residuals <- sweep(outcomes, 2L, sample_mean)
    effect_ss <- n * sum((contrast_projection %*% sample_mean)^2)
    error_ss <- sum((residuals %*% contrast_projection)^2)
    f_statistic <- (effect_ss / term_df) /
      (error_ss / (term_df * (n - 1)))
    f_statistic > critical_f
  })

  c(power = mean(rejected), epsilon = epsilon)
}

ar_covariance <- outer(1:4, 1:4, function(i, j) 0.8^abs(i - j))
directions <- list(
  default_linear = c(-3, -1, 1, 3),
  curved = c(1, -1, -1, 1),
  level_dominated = c(0, 0, 0, 1)
)

corrected <- vapply(
  directions,
  direct_rm_power,
  numeric(2L),
  covariance = ar_covariance,
  n = 25,
  target_pes = 0.20,
  gg_corrected = TRUE
)
uncorrected <- vapply(
  directions,
  direct_rm_power,
  numeric(2L),
  covariance = ar_covariance,
  n = 25,
  target_pes = 0.20,
  gg_corrected = FALSE
)

print(corrected)
print(uncorrected)
stopifnot(
  max(corrected["power", ]) - min(corrected["power", ]) > 0.03,
  max(uncorrected["power", ]) - min(uncorrected["power", ]) > 0.03,
  max(abs(corrected["epsilon", ] - corrected["epsilon", 1L])) < 1e-12
)
