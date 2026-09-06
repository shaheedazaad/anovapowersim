# Independent one-df reference: square a shifted normal, then integrate over
# the independent chi-squared denominator. The critical value comes from t.
tail_power_reference <- function(n, pes, alpha) {
  df <- 2 * n - 2
  critical <- stats::qt(alpha / 2, df, lower.tail = FALSE)^2
  shift <- sqrt(df * pes / (1 - pes))
  stats::integrate(function(z) {
    stats::dnorm(z) * stats::pchisq(df * (z + shift)^2 / critical, df)
  }, -Inf, Inf, rel.tol = 1e-10)$value
}

test_that("calculated achieved power preserves very small upper tails", {
  for (alpha in c(0.05, 1e-20)) {
    result <- power_achieved_calc(
      between = c(group = 2), term = "group", n = 100,
      target_pes = 0.5, alpha = alpha
    )
    expect_equal(result$achieved_power,
      tail_power_reference(100, 0.5, alpha),
      tolerance = 1e-8
    )
  }
})

test_that("calculated searches reach targets with very small alpha", {
  result <- suppressWarnings(power_n_calc(
    between = c(group = 2), term = "group", target_pes = 0.5,
    alpha = 1e-20, power = 0.90, n_max = 100
  ))
  expect_equal(result$n_needed, 78L)
  expect_lt(tail_power_reference(77, 0.5, 1e-20), 0.90)
  expect_gte(tail_power_reference(78, 0.5, 1e-20), 0.90)

  sensitivity <- suppressWarnings(power_sensitivity_calc(
    between = c(group = 2), term = "group", n = 100,
    alpha = 1e-20, power = 0.90, pes_tol = 1e-6
  ))
  expect_true(sensitivity$converged)
  if (is.finite(sensitivity$pes_needed)) {
    expect_gte(tail_power_reference(100, sensitivity$pes_needed, 1e-20), 0.90)
    expect_lt(tail_power_reference(100, sensitivity$pes_lower, 1e-20), 0.90)
  }
})

test_that("simulation diagnostics and starting estimates preserve upper tails", {
  spec <- balanced_anova_design(between = c(group = 2))
  expected <- tail_power_reference(100, 0.5, 1e-20)
  starting_power <- anovapowersim:::power_calc_at_n(
    spec = spec, term = "group", target_pes = 0.5, n = 100,
    alpha = 1e-20, ss_type = "III", sd = 1, r = 0.5, gpower = FALSE
  )
  expect_equal(starting_power, expected, tolerance = 1e-8)

  achieved <- suppressWarnings(power_achieved(
    between = c(group = 2), term = "group", target_pes = 0.5,
    n = 100, alpha = 1e-20, n_sims = 1, seed = 123, progress = FALSE
  ))
  expect_equal(achieved$calculated_power, round(expected, 3))
})
