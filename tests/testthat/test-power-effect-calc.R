quiet_power_achieved_calc <- function(...) {
  suppressWarnings(power_achieved_calc(...))
}

quiet_power_sensitivity_calc <- function(...) {
  suppressWarnings(power_sensitivity_calc(...))
}

test_that("fixed-sample calculated-power functions are exported", {
  exports <- getNamespaceExports("anovapowersim")
  expect_true("power_achieved_calc" %in% exports)
  expect_true("power_sensitivity_calc" %in% exports)
  expect_equal(formals(power_sensitivity_calc)$power, 0.90)
  expect_equal(formals(power_sensitivity_calc)$pes_tol, 0.001)
})

test_that("power_achieved_calc returns the direct calculated-power row", {
  designs <- list(
    list(between = c(group = 2), within = NULL, term = "group"),
    list(between = NULL, within = c(time = 3), term = "time"),
    list(
      between = c(group = 2), within = c(time = 3), term = "group:time"
    )
  )

  for (design in designs) {
    result <- do.call(
      quiet_power_achieved_calc,
      c(design, list(target_pes = 0.12, n = 8))
    )
    expected <- anovapowersim:::analytic_power_row(
      n = 8,
      spec = result$design,
      term = result$term,
      target_pes = 0.12,
      alpha = 0.05,
      gpower = FALSE,
      epsilon = 1
    )

    expect_s3_class(result, "anovapowersim_achieved_power")
    expect_true(result$calculation_only)
    expect_equal(result$results, expected)
    expect_equal(result$achieved_power, expected$power_calc[[1L]])
    expect_equal(result$calculated_power, expected$power_calc[[1L]])
    expect_true(is.na(result$results$n_sims[[1L]]))
    expect_true(is.na(result$results$power_sim[[1L]]))
  }
})

test_that("power_achieved_calc matches power_n_calc at a fixed n", {
  achieved <- quiet_power_achieved_calc(
    between = c(group = 2),
    within = c(time = 3),
    term = "group:time",
    target_pes = 0.14,
    n = 12,
    gpower = TRUE,
    epsilon = 0.8
  )
  sample_size <- suppressWarnings(power_n_calc(
    between = c(group = 2),
    within = c(time = 3),
    term = "group:time",
    target_pes = 0.14,
    power = 0.99,
    n_start = 12,
    n_max = 12,
    gpower = TRUE,
    epsilon = 0.8
  ))

  expect_equal(achieved$results, sample_size$results)
})

test_that("power_sensitivity_calc returns a calculated upper bracket", {
  result <- quiet_power_sensitivity_calc(
    between = c(group = 2),
    within = c(time = 3),
    term = "group:time",
    n = 30,
    power = 0.90,
    pes_tol = 1e-5,
    gpower = TRUE,
    epsilon = 0.8
  )

  expect_s3_class(result, "anovapowersim_sensitivity")
  expect_true(result$calculation_only)
  expect_true(result$pes_needed %in% result$results$target_pes)
  expect_lte(result$bracket_width, result$pes_tol)
  upper <- result$results[result$results$target_pes == result$pes_needed, ]
  lower <- result$results[result$results$target_pes == result$pes_lower, ]
  expect_gte(upper$power_calc[[1L]], result$power)
  expect_lt(lower$power_calc[[1L]], result$power)
  expect_true(all(is.na(result$results$n_sims)))
  expect_true(all(is.na(result$results$power_sim)))
})

test_that("calculated sensitivity agrees with calculated achieved power", {
  sensitivity <- quiet_power_sensitivity_calc(
    within = c(time = 3),
    term = "time",
    n = 20,
    power = 0.90,
    pes_tol = 1e-6,
    epsilon = 0.8
  )
  achieved <- quiet_power_achieved_calc(
    within = c(time = 3),
    term = "time",
    target_pes = sensitivity$pes_needed,
    n = 20,
    epsilon = 0.8
  )

  expect_gte(achieved$calculated_power, sensitivity$power)
})

test_that("power_sensitivity_calc reports unreachable targets", {
  expect_warning(
    result <- power_sensitivity_calc(
      between = c(group = 2),
      term = "group",
      n = 5,
      power = 0.90,
      pes_min = 1e-6,
      pes_max = 2e-6
    ),
    "pes_max"
  )
  expect_true(is.na(result$pes_needed))
})

test_that("fixed-sample calculated power validates n, bounds, and epsilon", {
  expect_error(
    quiet_power_achieved_calc(
      between = c(group = 2), term = "group", target_pes = 0.1, n = 1
    ),
    "too small for calculated power"
  )
  expect_error(
    quiet_power_sensitivity_calc(
      between = c(group = 2), term = "group", n = 5,
      pes_min = 0.2, pes_max = 0.1
    ),
    "smaller"
  )
  expect_error(
    quiet_power_achieved_calc(
      between = c(group = 2), term = "group", target_pes = 0.1,
      n = 5, epsilon = 0.8
    ),
    "purely between-subject"
  )
})

test_that("calculated fixed-N output distinguishes total and per-cell n", {
  within_result <- quiet_power_achieved_calc(
    within = c(time = 2), term = "time", target_pes = 0.1, n = 10
  )
  mixed_result <- quiet_power_achieved_calc(
    between = c(group = 2), within = c(time = 2), term = "group:time",
    target_pes = 0.1, n = 10
  )

  expect_equal(within_result$total_n, 10L)
  expect_equal(mixed_result$total_n, 20L)
  expect_output(print(within_result), "fixed total N")
  expect_output(print(mixed_result), "fixed n per cell")
})

test_that("calculated fixed-N methods identify calculated power", {
  achieved <- quiet_power_achieved_calc(
    between = c(group = 2), term = "group", target_pes = 0.1, n = 10
  )
  sensitivity <- quiet_power_sensitivity_calc(
    between = c(group = 2), term = "group", n = 10, pes_tol = 0.01
  )

  expect_output(print(achieved), "calculated power only")
  expect_output(summary(achieved), "calculated achieved-power summary")
  expect_output(print(sensitivity), "calculated power only")
  expect_output(summary(sensitivity), "calculated-power sensitivity summary")
})
