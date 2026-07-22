quiet_power_achieved <- function(...) suppressWarnings(power_achieved(...))
quiet_power_sensitivity <- function(...) suppressWarnings(power_sensitivity(...))

test_that("power_achieved exactly matches a one-point power_curve", {
  designs <- list(
    list(between = c(group = 2), within = NULL, term = "group", n = 5),
    list(between = NULL, within = c(time = 2), term = "time", n = 5),
    list(
      between = c(group = 2), within = c(time = 2),
      term = "group:time", n = 5
    )
  )

  for (design in designs) {
    common <- c(
      design[c("between", "within", "term")],
      list(target_pes = 0.14, n_sims = 8, progress = FALSE, seed = 42)
    )
    achieved <- do.call(
      quiet_power_achieved,
      c(common, list(n = design$n))
    )
    curve <- do.call(
      function(...) suppressWarnings(power_curve(...)),
      c(common, list(n_range = design$n))
    )

    expect_equal(achieved$results, curve$results)
    expect_equal(achieved$achieved_power, curve$results$power_sim[[1L]])
    expect_equal(achieved$calculated_power, curve$results$power_calc[[1L]])
  }
})

test_that("power_achieved returns its dedicated class and fixed-N fields", {
  result <- quiet_power_achieved(
    between = c(group = 2),
    term = "group",
    target_pes = 0.1,
    n = 5,
    n_sims = 4,
    progress = FALSE,
    seed = 3
  )

  expect_s3_class(result, "anovapowersim_achieved_power")
  expect_equal(result$n, 5L)
  expect_equal(result$total_n, 10L)
  expect_equal(nrow(result$results), 1L)
  expect_true(all(c("achieved_power", "calculated_power") %in% names(result)))
})

test_that("power_sensitivity returns an explicitly simulated upper bracket", {
  result <- quiet_power_sensitivity(
    between = c(group = 2),
    term = "group",
    n = 8,
    power = 0.80,
    n_sims = 30,
    pes_min = 0.001,
    pes_max = 0.9,
    pes_tol = 0.02,
    progress = FALSE,
    seed = 12
  )

  expect_s3_class(result, "anovapowersim_sensitivity")
  expect_true(result$pes_needed %in% result$results$target_pes)
  upper <- result$results[result$results$target_pes == result$pes_needed, ]
  expect_gte(upper$power_sim[[1L]], result$power)
  expect_lte(result$bracket_width, result$pes_tol)
  if (is.finite(result$pes_lower)) {
    lower <- result$results[result$results$target_pes == result$pes_lower, ]
    expect_lt(lower$power_sim[[1L]], result$power)
  }
})

test_that("power_sensitivity is reproducible for a fixed seed", {
  args <- list(
    within = c(time = 2),
    term = "time",
    n = 7,
    power = 0.80,
    n_sims = 12,
    pes_tol = 0.05,
    progress = FALSE,
    seed = 99
  )
  first <- do.call(quiet_power_sensitivity, args)
  second <- do.call(quiet_power_sensitivity, args)

  expect_equal(first$results, second$results)
  expect_equal(first$pes_needed, second$pes_needed)
})

test_that("power_sensitivity reports an unreachable upper bound", {
  warnings <- testthat::capture_warnings(
    result <- power_sensitivity(
      between = c(group = 2),
      term = "group",
      n = 5,
      power = 0.99,
      n_sims = 4,
      alpha = 1e-8,
      pes_min = 1e-6,
      pes_max = 2e-6,
      progress = FALSE,
      seed = 7
    )
  )
  expect_true(any(grepl("pes_max", warnings, fixed = TRUE)))
  expect_true(any(grepl(
    "No SD or `covariance` was supplied", warnings, fixed = TRUE
  )))
  expect_true(is.na(result$pes_needed))
  expect_true(2e-6 %in% result$results$target_pes)
})

test_that("fixed-N APIs validate n and sensitivity controls", {
  common <- list(between = c(group = 2), term = "group", progress = FALSE)
  expect_error(
    suppressWarnings(do.call(
      power_achieved, c(common, list(target_pes = 0.1, n = 1.5))
    )),
    "single positive integer"
  )
  expect_error(
    suppressWarnings(do.call(
      power_sensitivity,
      c(common, list(n = 5, pes_min = 0.2, pes_max = 0.1))
    )),
    "smaller"
  )
  expect_error(
    suppressWarnings(do.call(
      power_sensitivity, c(common, list(n = 5, pes_tol = 0))
    )),
    "positive finite"
  )
})

test_that("fixed-N APIs retain covariance and analysis controls", {
  sigma <- matrix(c(1, 0.8, 0, 0.8, 1, 0, 0, 0, 1), nrow = 3)
  dimnames(sigma) <- list(paste0("time", 1:3), paste0("time", 1:3))
  result <- quiet_power_achieved(
    within = c(time = 3),
    term = "time",
    target_pes = 0.15,
    n = 5,
    n_sims = 3,
    covariance = sigma,
    ss_type = "I",
    gpower = TRUE,
    progress = FALSE,
    seed = 4
  )

  expect_equal(result$epsilon, 0.6540541, tolerance = 1e-7)
  expect_true(result$custom_covariance)
  expect_identical(result$ss_type, "I")
  expect_true(result$gpower)
})

test_that("new result classes have dedicated print and summary methods", {
  achieved <- quiet_power_achieved(
    between = c(group = 2), term = "group", target_pes = 0.1,
    n = 5, n_sims = 2, progress = FALSE, seed = 1
  )
  sensitivity <- quiet_power_sensitivity(
    between = c(group = 2), term = "group", n = 5, power = 0.8,
    n_sims = 2, pes_tol = 0.2, progress = FALSE, seed = 1
  )

  expect_output(print(achieved), "achieved power")
  expect_output(summary(achieved), "achieved-power summary")
  expect_output(print(sensitivity), "detectable pes")
  expect_output(summary(sensitivity), "sensitivity summary")
})

test_that("continuous effect search warns through non-convergence metadata", {
  run_one <- function(pes) {
    tibble::tibble(target_pes = pes, power_sim = pes, power_calc = pes)
  }
  search <- anovapowersim:::adaptive_effect_search(
    run_one = run_one,
    target = 0.8,
    pes_start = 0.8,
    pes_min = 0.1,
    pes_max = 0.9,
    pes_tol = 1e-8,
    max_iter = 1L
  )

  expect_true(search$reached)
  expect_false(search$converged)
  expect_true(search$pes_upper %in% search$results$target_pes)
})
