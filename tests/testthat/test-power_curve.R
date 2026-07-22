quiet_power_curve <- function(...) suppressWarnings(power_curve(...))
quiet_power_n <- function(...) suppressWarnings(power_n(...))
quiet_power_n_calc <- function(...) suppressWarnings(power_n_calc(...))

capture_warning_messages <- function(expr) {
  warnings <- character()
  withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  warnings
}

test_that("power_n defaults to 90 percent power", {
  expect_equal(formals(power_n)$power, 0.90)
  expect_equal(formals(power_n)$n_max, 5000)
})

test_that("simulation APIs expose auto correction by default", {
  simulation_functions <- list(
    power_curve, power_n, power_achieved, power_sensitivity, power_unbalanced
  )
  for (fn in simulation_functions) {
    expect_identical(
      formals(fn)$sim_correction,
      quote(c("auto", "GG", "none"))
    )
  }
  expect_null(formals(power_n_calc)$sim_correction)
})

test_that("power_n_calc is exported and defaults to 90 percent power", {
  expect_true("power_n_calc" %in% getNamespaceExports("anovapowersim"))
  expect_equal(formals(power_n_calc)$power, 0.90)
  expect_equal(formals(power_n_calc)$n_max, 5000)
  expect_equal(formals(power_n_calc)$epsilon, 1)
})

test_that("power_n warns when requested power is below 90 percent", {
  low_power_warnings <- capture_warning_messages(
    expect_error(
      power_n(
        between = c(group = 2),
        within = c(time = 2),
        term = "group:time",
        target_pes = 0.1,
        power = 0.80,
        n_sims = 1,
        n_max = 1,
        progress = FALSE
      )
    )
  )
  expect_true(any(grepl(
    "Power greater than or equal to .90 is recommended.",
    low_power_warnings,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "No `covariance` was supplied",
    low_power_warnings,
    fixed = TRUE
  )))

  warnings <- capture_warning_messages(
    expect_error(
      power_n(
        between = c(group = 2),
        within = c(time = 2),
        term = "group:time",
        target_pes = 0.1,
        power = 0.90,
        n_sims = 1,
        n_max = 1,
        progress = FALSE
      )
    )
  )
  expect_false(any(grepl("Power greater than or equal to .90", warnings,
                         fixed = TRUE)))
})

test_that("power_n_calc warns when requested power is below 90 percent", {
  expect_warning(
    power_n_calc(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.1,
      power = 0.80,
      n_max = 5000
    ),
    "Power greater than or equal to .90 is recommended.",
    fixed = TRUE
  )
})

test_that("rule-of-thumb medium partial eta squared warns at call time", {
  medium_message <- paste(
    "It looks like you are using a rule-of-thumb \"medium\" effect size.",
    "This might overestimate the true effect size, rendering your study",
    "underpowered. Consider basing your power calculations on previous",
    "research or empirically-derived guidelines."
  )

  curve_warnings <- capture_warning_messages(
    expect_error(
      power_curve(
        between = c(group = 2),
        within = c(time = 2),
        term = "group:time",
        target_pes = 0.06,
        n_range = 0,
        n_sims = 1,
        progress = FALSE
      )
    )
  )
  expect_true(any(grepl(medium_message, curve_warnings, fixed = TRUE)))
  expect_true(any(grepl(
    "No `covariance` was supplied", curve_warnings, fixed = TRUE
  )))

  n_warnings <- capture_warning_messages(
    expect_error(
      power_n(
        between = c(group = 2),
        within = c(time = 2),
        term = "group:time",
        target_pes = 0.06,
        n_sims = 1,
        n_max = 1,
        progress = FALSE
      )
    )
  )
  expect_true(any(grepl(medium_message, n_warnings, fixed = TRUE)))
  expect_true(any(grepl(
    "No `covariance` was supplied", n_warnings, fixed = TRUE
  )))
  expect_warning(
    expect_error(
      power_n_calc(
        between = c(group = 2),
        within = c(time = 2),
        term = "group:time",
        target_pes = 0.06,
        n_start = 1,
        n_max = 1
      )
    ),
    medium_message,
    fixed = TRUE
  )

  for (nearby_pes in c(0.059, 0.061)) {
    curve_warnings <- capture_warning_messages(
      expect_error(
        power_curve(
          between = c(group = 2),
          within = c(time = 2),
          term = "group:time",
          target_pes = nearby_pes,
          n_range = 0,
          n_sims = 1,
          progress = FALSE
        )
      )
    )
    n_warnings <- capture_warning_messages(
      expect_error(
        power_n(
          between = c(group = 2),
          within = c(time = 2),
          term = "group:time",
          target_pes = nearby_pes,
          n_sims = 1,
          n_max = 1,
          progress = FALSE
        )
      )
    )
    n_calc_warnings <- capture_warning_messages(
      expect_error(
        power_n_calc(
          between = c(group = 2),
          within = c(time = 2),
          term = "group:time",
          target_pes = nearby_pes,
          n_start = 1,
          n_max = 1
        )
      )
    )

    expect_false(any(grepl("rule-of-thumb \"medium\" effect size",
                           curve_warnings, fixed = TRUE)))
    expect_false(any(grepl("rule-of-thumb \"medium\" effect size",
                           n_warnings, fixed = TRUE)))
    expect_false(any(grepl("rule-of-thumb \"medium\" effect size",
                           n_calc_warnings, fixed = TRUE)))
  }
})

test_that("power_curve simulates a balanced mixed design", {
  set.seed(1)
  pc <- quiet_power_curve(
    between = c(color = 2),
    within = c(age = 2),
    term = "color:age",
    target_pes = 0.20721,
    n_range = c(8L, 10L),
    n_sims = 20,
    seed = 1
  )

  expect_s3_class(pc, "anovapowersim_curve")
  expect_equal(pc$results$n_per_cell, c(8L, 10L))
  expect_equal(pc$results$total_n, c(16L, 20L))
  expect_true(all(c(
    "valid_sims", "failed_sims", "num_df", "den_df", "ncp",
    "power_calc", "power_sim"
  ) %in% names(pc$results)))
  expect_false(any(c("power_sim_lower", "power_sim_upper") %in%
                   names(pc$results)))
  expect_equal(pc$results$valid_sims + pc$results$failed_sims,
               pc$results$n_sims)
  expect_true(all(pc$results$power_sim >= 0 & pc$results$power_sim <= 1))
  expect_equal(pc$target_pes, 0.20721)
  expect_equal(pc$ss_type, "III")
  expect_identical(pc$sim_correction, "auto")
  expect_identical(pc$sim_correction_resolved, "none")
})

test_that("ss_type is validated and stored", {
  pc <- quiet_power_curve(
    between = c(color = 2),
    within = c(age = 2),
    term = "color:age",
    target_pes = 0.1,
    n_range = 8,
    n_sims = 5,
    ss_type = "I",
    seed = 11
  )
  expect_equal(pc$ss_type, "I")

  expect_error(
    quiet_power_curve(
      between = c(color = 2),
      term = "color",
      target_pes = 0.1,
      n_range = 8,
      n_sims = 5,
      ss_type = "IV"
    ),
    "`ss_type`"
  )
})

test_that("sim_correction is validated by public balanced APIs", {
  expect_error(
    quiet_power_curve(
      within = c(time = 3), term = "time", target_pes = 0.1,
      n_range = 5, n_sims = 1, sim_correction = "HF"
    ),
    "one of"
  )
  expect_error(
    quiet_power_curve(
      within = c(time = 3), term = "time", target_pes = 0.1,
      n_range = 5, n_sims = 1, ss_type = "I", sim_correction = "GG"
    ),
    'sim_correction = "GG".*ss_type = "I"'
  )
})

test_that("printed output uses concise public labels", {
  pc <- quiet_power_curve(
    between = c(color = 2),
    within = c(age = 2),
    term = "color:age",
    target_pes = 0.20721,
    n_range = 8,
    n_sims = 10,
    seed = 13
  )
  printed <- capture.output(print(pc))

  expect_true(any(grepl("n needed for between-subjects cell", printed,
                       fixed = TRUE)))
  expect_false(any(grepl("G\\*Power NCP", printed)))
  expect_false(any(grepl("G\\*Power convention", printed)))
  expect_true(any(grepl("\\b[0-9]\\.[0-9]{3}\\b", printed)))
  expect_snapshot_output(print(pc))

  pc_gpower <- quiet_power_curve(
    between = c(color = 2),
    within = c(age = 2),
    term = "color:age",
    target_pes = 0.20721,
    n_range = 8,
    n_sims = 10,
    gpower = TRUE,
    seed = 14
  )
  printed_gpower <- capture.output(print(pc_gpower))
  expect_true(any(grepl("G\\*Power convention: TRUE", printed_gpower)))
})

test_that("power_n_calc returns calculated-power-only result columns", {
  pc <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    power = 0.90,
    n_max = 40
  )

  expect_s3_class(pc, "anovapowersim_curve")
  expect_true(all(c(
    "n_per_cell", "total_n", "n_sims", "num_df", "den_df", "ncp",
    "power_calc", "power_sim"
  ) %in% names(pc$results)))
  expect_true(all(is.na(pc$results$n_sims)))
  expect_true(all(is.na(pc$results$power_sim)))
  expect_true(all(is.finite(pc$results$power_calc)))
  expect_true(is.na(pc$n_sims))
  expect_null(pc$ss_type)
})

test_that("power_n_calc print output does not report simulations", {
  pc <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    n_max = 20
  )
  printed <- capture.output(print(pc))
  summarized <- capture.output(summary(pc))

  expect_false(any(grepl("sims per cell size", printed, fixed = TRUE)))
  expect_false(any(grepl("simulated test", printed, fixed = TRUE)))
  expect_true(any(grepl("calculation:   calculated power only", printed,
                       fixed = TRUE)))
  expect_true(any(grepl("calculated-power summary", summarized,
                       fixed = TRUE)))
})

test_that("power_n_calc reports the smallest calculated n meeting target", {
  pc <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    power = 0.90,
    n_max = 80
  )

  expect_false(is.na(pc$n_needed))
  expect_equal(pc$n_needed, min(pc$results$n_per_cell[
    pc$results$power_calc >= pc$power
  ]))
  previous <- pc$results[pc$results$n_per_cell < pc$n_needed, , drop = FALSE]
  expect_true(nrow(previous) == 0L || all(previous$power_calc < pc$power))
  expect_equal(pc$total_n_needed, pc$n_needed * 2L)
})

test_that("power_n_calc searches below an adequate n_start", {
  baseline <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    power = 0.90,
    n_max = 100
  )
  pc <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    power = 0.90,
    n_start = 100,
    n_max = 100
  )

  expect_equal(pc$n_needed, baseline$n_needed)
  expect_lt(pc$n_needed, 100L)
  expect_true(all(c(pc$n_needed - 1L, pc$n_needed, 100L) %in%
                  pc$results$n_per_cell))
  previous_power <- pc$results$power_calc[
    pc$results$n_per_cell == pc$n_needed - 1L
  ]
  expect_lt(previous_power, pc$power)
})

test_that("power_n_calc does not exhaustively visit n_max after finding target", {
  pc <- quiet_power_n_calc(
    between = c(cond = 2),
    within = c(stim = 4),
    term = "cond:stim",
    target_pes = 0.14,
    power = 0.90
  )
  previous <- anovapowersim:::analytic_power_row(
    n = pc$n_needed - 1L,
    spec = pc$design,
    term = pc$term,
    target_pes = pc$target_pes,
    alpha = pc$alpha,
    gpower = pc$gpower
  )

  expect_equal(pc$n_needed, 17L)
  expect_equal(pc$total_n_needed, 34L)
  expect_lt(nrow(pc$results), 20L)
  expect_false(1000L %in% pc$results$n_per_cell)
  expect_lt(previous$power_calc, pc$power)
})

test_that("power_n_calc reports unreached target when n_max is too small", {
  expect_warning(
    pc <- power_n_calc(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.01,
      power = 0.99,
      n_max = 4
    ),
    "Target power 0.990 was not reached by `n_max = 4`. Increase `n_max`",
    fixed = TRUE
  )

  expect_true(is.na(pc$n_needed))
  expect_true(is.na(pc$total_n_needed))
  expect_true(all(pc$results$power_calc < pc$power))
})

test_that("power_n warns when target power is not reached by n_max", {
  warnings <- capture_warning_messages(
    pc <- power_n(
      between = c(group = 2),
      term = "group",
      target_pes = 1e-6,
      power = 0.99,
      n_sims = 10,
      n_start = 2,
      n_max = 2,
      progress = FALSE,
      seed = 202
    )
  )

  expect_true(any(grepl(
    "Target power 0.990 was not reached by `n_max = 2`. Increase `n_max`",
    warnings,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "No SD or `covariance` was supplied", warnings, fixed = TRUE
  )))
  expect_true(is.na(pc$n_needed))
  expect_true(is.na(pc$total_n_needed))
})

test_that("power_n_calc computes calculated-power dfs for common balanced designs", {
  between_pc <- quiet_power_n_calc(
    between = c(group = 3),
    term = "group",
    target_pes = 0.2,
    n_start = 8,
    n_max = 8
  )
  expect_equal(between_pc$results$num_df, 2)
  expect_equal(between_pc$results$den_df, 21)

  within_pc <- quiet_power_n_calc(
    within = c(time = 4),
    term = "time",
    target_pes = 0.2,
    n_start = 8,
    n_max = 8
  )
  expect_equal(within_pc$results$total_n, 8L)
  expect_equal(within_pc$results$num_df, 3)
  expect_equal(within_pc$results$den_df, 21)

  mixed_pc <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 3),
    term = "group:time",
    target_pes = 0.2,
    n_start = 8,
    n_max = 8
  )
  expect_equal(mixed_pc$results$total_n, 16L)
  expect_equal(mixed_pc$results$num_df, 2)
  expect_equal(mixed_pc$results$den_df, 28)
})

test_that("power_n_calc applies epsilon to dfs, ncp, and power", {
  epsilon <- 0.5
  target_pes <- 0.15
  f2 <- target_pes / (1 - target_pes)

  pc <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 3),
    term = "group:time",
    target_pes = target_pes,
    n_start = 8,
    n_max = 8,
    epsilon = epsilon
  )
  expected_ncp <- epsilon * 28 * f2
  expected_power <- stats::pf(
    stats::qf(0.95, df1 = 1, df2 = 14),
    df1 = 1,
    df2 = 14,
    ncp = expected_ncp,
    lower.tail = FALSE
  )

  expect_equal(pc$epsilon, epsilon)
  expect_equal(pc$results$epsilon, epsilon)
  expect_equal(pc$results$num_df, 1)
  expect_equal(pc$results$den_df, 14)
  expect_equal(pc$results$ncp, expected_ncp, tolerance = 1e-12)
  expect_equal(pc$results$power_calc, expected_power, tolerance = 1e-12)

  gp <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 3),
    term = "group:time",
    target_pes = target_pes,
    n_start = 8,
    n_max = 8,
    gpower = TRUE,
    epsilon = epsilon
  )
  expect_equal(gp$results$ncp, epsilon * 16 * f2, tolerance = 1e-12)
})

test_that("power_n_calc validates epsilon for the tested term", {
  common_args <- list(
    within = c(time = 4),
    term = "time",
    target_pes = 0.2,
    n_start = 8,
    n_max = 8
  )

  for (bad_epsilon in list(0, 1.1, NA_real_, Inf, c(0.5, 0.6), "0.5")) {
    expect_error(
      do.call(quiet_power_n_calc, c(common_args, list(epsilon = bad_epsilon))),
      "`epsilon` must be a single finite number"
    )
  }

  expect_error(
    quiet_power_n_calc(
      between = c(group = 2),
      term = "group",
      target_pes = 0.2,
      n_start = 8,
      n_max = 8,
      epsilon = 0.8
    ),
    "purely between-subject term"
  )
  expect_error(
    quiet_power_n_calc(
      within = c(time = 2),
      term = "time",
      target_pes = 0.2,
      n_start = 8,
      n_max = 8,
      epsilon = 0.9
    ),
    "must be at least 1"
  )
  expect_error(
    do.call(quiet_power_n_calc, c(common_args, list(epsilon = 0.3))),
    "must be at least"
  )

  at_lower_bound <- do.call(
    quiet_power_n_calc,
    c(common_args, list(epsilon = 1 / 3))
  )
  expect_equal(at_lower_bound$epsilon, 1 / 3)
})

test_that("power_n_calc reports a non-unit epsilon", {
  pc <- quiet_power_n_calc(
    within = c(time = 4),
    term = "time",
    target_pes = 0.2,
    n_start = 8,
    n_max = 8,
    epsilon = 0.7
  )
  printed <- capture.output(print(pc))
  summarized <- capture.output(summary(pc))

  expect_true(any(grepl("epsilon:       0.7", printed, fixed = TRUE)))
  expect_true(any(grepl("epsilon:      0.7000", summarized, fixed = TRUE)))
})

test_that("power_n_calc G*Power convention changes ncp and remains finite", {
  base <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 2, condition = 3),
    term = "group:time",
    target_pes = 0.15,
    n_start = 16,
    n_max = 16
  )
  gp <- quiet_power_n_calc(
    between = c(group = 2),
    within = c(time = 2, condition = 3),
    term = "group:time",
    target_pes = 0.15,
    n_start = 16,
    n_max = 16,
    gpower = TRUE
  )

  expect_true(gp$gpower)
  expect_false(isTRUE(all.equal(base$results$ncp, gp$results$ncp)))
  expect_equal(gp$results$ncp, 32 * 0.15 / (1 - 0.15), tolerance = 1e-12)
  expect_true(is.finite(gp$results$power_calc))
})

test_that("power_n_calc matches car-backed calculated power for balanced designs", {
  compare_at_n <- function(between = NULL, within = NULL, term, n) {
    calc <- quiet_power_n_calc(
      between = between,
      within = within,
      term = term,
      target_pes = 0.15,
      power = 0.90,
      n_start = n,
      n_max = n
    )
    car_backed <- quiet_power_n(
      between = between,
      within = within,
      term = term,
      target_pes = 0.15,
      power = 0.90,
      n_sims = 1,
      n_start = n,
      n_max = n,
      seed = 123
    )

    car_at_n <- car_backed$results[
      car_backed$results$n_per_cell == n, , drop = FALSE
    ]
    expect_equal(calc$results$num_df, car_at_n$num_df)
    expect_equal(calc$results$den_df, car_at_n$den_df)
    expect_equal(calc$results$ncp, car_at_n$ncp, tolerance = 0.001)
    expect_equal(round(calc$results$power_calc, 3),
                 car_at_n$power_calc,
                 tolerance = 1e-12)
  }

  compare_at_n(between = c(group = 3), term = "group", n = 8)
  compare_at_n(within = c(time = 2), term = "time", n = 8)
  compare_at_n(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    n = 8
  )
})

test_that("power disagreement warning suggests the right next step", {
  unstable <- tibble::tibble(
    n_per_cell = c(10L, 20L),
    power_calc = c(0.70, 0.80),
    power_sim = c(0.62, 0.79)
  )

  expect_warning(
    anovapowersim:::warn_power_disagreement(unstable, n_sims = 1000),
    "Try increasing `n_sims`"
  )
  expect_warning(
    anovapowersim:::warn_power_disagreement(unstable, n_sims = 5000),
    "raising a GitHub issue"
  )
})

make_power_search_runner <- function(powers) {
  env <- new.env(parent = emptyenv())
  env$visited <- integer()
  env$run_one <- function(n) {
    env$visited <- c(env$visited, as.integer(n))
    key <- as.character(n)
    if (!key %in% names(powers)) stop("Unexpected n: ", n, call. = FALSE)
    tibble::tibble(
      n_per_cell = as.integer(n),
      total_n = as.integer(n),
      n_sims = 1L,
      num_df = 1,
      den_df = 1,
      ncp = NA_real_,
      power_calc = NA_real_,
      power_sim = powers[[key]]
    )
  }
  env
}

make_formula_power_search_runner <- function(power_fn) {
  env <- new.env(parent = emptyenv())
  env$visited <- integer()
  env$run_one <- function(n) {
    env$visited <- c(env$visited, as.integer(n))
    tibble::tibble(
      n_per_cell = as.integer(n),
      total_n = as.integer(n),
      n_sims = 1L,
      num_df = 1,
      den_df = 1,
      ncp = NA_real_,
      power_calc = NA_real_,
      power_sim = power_fn(n)
    )
  }
  env
}

test_that("adaptive search uses one-sided precision above target", {
  below_close <- make_power_search_runner(c(
    `10` = 0.79,
    `12` = 0.79,
    `13` = 0.79,
    `14` = 0.81,
    `20` = 0.82
  ))
  below_curve <- anovapowersim:::adaptive_design_search(
    run_one = below_close$run_one,
    target = 0.80,
    n_start = 10,
    n_max = 20,
    tol = 0.03
  )

  expect_equal(below_close$visited, c(10L, 20L, 14L, 12L, 13L))
  expect_equal(below_curve$n_per_cell, c(10L, 12L, 13L, 14L, 20L))
  expect_equal(anovapowersim:::estimate_design_n_needed(below_curve, 0.80),
               14L)

  in_band <- make_power_search_runner(c(
    `2` = 0.70,
    `8` = 0.79,
    `9` = 0.81,
    `10` = 0.82
  ))
  in_band_curve <- anovapowersim:::adaptive_design_search(
    run_one = in_band$run_one,
    target = 0.80,
    n_start = 10,
    n_max = 20,
    tol = 0.03,
    n_min = 2
  )

  expect_equal(in_band$visited, c(10L, 2L, 9L, 8L))
  expect_equal(in_band_curve$n_per_cell, c(2L, 8L, 9L, 10L))
  expect_equal(anovapowersim:::estimate_design_n_needed(in_band_curve, 0.80),
               9L)
})

test_that("adaptive search searches below n_start when it already reaches target", {
  runner <- make_formula_power_search_runner(function(n) n / 100)
  curve <- anovapowersim:::adaptive_design_search(
    run_one = runner$run_one,
    target = 0.40,
    n_start = 100,
    n_max = 100,
    tol = 0.03,
    n_min = 2
  )

  expect_equal(runner$visited, c(100L, 2L, 40L, 39L))
  expect_equal(anovapowersim:::estimate_design_n_needed(curve, 0.40), 40L)
})

test_that("adaptive search simulates interpolated candidates after overshoot", {
  runner <- make_power_search_runner(c(
    `10` = 0.70,
    `13` = 0.79,
    `14` = 0.82,
    `20` = 0.95
  ))
  curve <- anovapowersim:::adaptive_design_search(
    run_one = runner$run_one,
    target = 0.80,
    n_start = 10,
    n_max = 20,
    tol = 0.03
  )

  expect_equal(runner$visited, c(10L, 20L, 14L, 13L))
  expect_equal(curve$n_per_cell, c(10L, 13L, 14L, 20L))
  expect_equal(anovapowersim:::estimate_design_n_needed(curve, 0.80), 14L)
})

test_that("adaptive search narrows the bracket without exhaustive scanning", {
  runner <- make_formula_power_search_runner(function(n) n / 1000)
  curve <- anovapowersim:::adaptive_design_search(
    run_one = runner$run_one,
    target = 0.90,
    n_start = 100,
    n_max = 1000,
    tol = 0.03
  )

  expect_equal(runner$visited, c(100L, 200L, 400L, 800L, 1000L, 900L, 899L))
  expect_equal(anovapowersim:::estimate_design_n_needed(curve, 0.90), 900L)
  expect_true(all(c(899L, 900L) %in% curve$n_per_cell))
  expect_lt(length(runner$visited), length(800:1000))
})

test_that("adaptive search warns when no simulated value reaches precision band", {
  runner <- make_power_search_runner(c(`10` = 0.70, `11` = 0.86, `12` = 0.95))
  curve <- anovapowersim:::adaptive_design_search(
    run_one = runner$run_one,
    target = 0.80,
    n_start = 10,
    n_max = 12,
    tol = 0.03
  )
  n_needed <- anovapowersim:::estimate_design_n_needed(curve, 0.80)

  expect_equal(runner$visited, c(10L, 12L, 11L))
  expect_equal(n_needed, 11L)
  expect_warning(
    anovapowersim:::warn_precision_band_not_reached(
      curve = curve,
      target = 0.80,
      tol = 0.03,
      n_needed = n_needed
    ),
    "Requested precision band was not reached.*0.800.*0.030.*11.*0.860"
  )
})

test_that("estimate_design_n_needed reports only simulated sample sizes", {
  interpolated_curve <- tibble::tibble(
    n_per_cell = c(10L, 20L),
    power_sim = c(0.70, 0.90)
  )
  precision_curve <- tibble::tibble(
    n_per_cell = c(8L, 12L, 16L),
    power_sim = c(0.78, 0.81, 0.84)
  )
  overshoot_curve <- tibble::tibble(
    n_per_cell = c(8L, 12L),
    power_sim = c(0.78, 0.94)
  )
  unreached_curve <- tibble::tibble(
    n_per_cell = c(8L, 12L),
    power_sim = c(0.70, 0.79)
  )

  expect_equal(
    anovapowersim:::estimate_design_n_needed(interpolated_curve, 0.80),
    20L
  )
  expect_equal(
    anovapowersim:::estimate_design_n_needed(precision_curve, 0.80),
    12L
  )
  expect_equal(
    anovapowersim:::estimate_design_n_needed(overshoot_curve, 0.80),
    12L
  )
  expect_true(is.na(
    anovapowersim:::estimate_design_n_needed(unreached_curve, 0.80)
  ))
})

test_that("power_n default tolerance is one-sided above target", {
  expect_equal(formals(power_n)$tol, 0.03)
  expect_false(anovapowersim:::power_is_in_precision_band(0.79, 0.80, 0.03))
  expect_true(anovapowersim:::power_is_in_precision_band(0.82, 0.80, 0.03))
  expect_false(anovapowersim:::power_is_in_precision_band(0.84, 0.80, 0.03))
})

test_that("design components are available for direct simulation", {
  d <- balanced_anova_design(
    between = c(group = 2),
    within = c(time = 2)
  )
  means <- design_term_means(
    d,
    term = "group:time",
    target_pes = 0.2,
    n = 12
  )
  sim <- simulate_design_dataset(d, n = 12, means = means)

  expect_s3_class(d, "anovapowersim_design_spec")
  expect_equal(dim(means), c(2L, 2L))
  expect_equal(dplyr::n_distinct(sim$id), 24L)
  expect_true(all(c("id", "group", "time", "value") %in% names(sim)))
})

test_that("default term means do not create non-target condition interactions", {
  d <- balanced_anova_design(
    between = c(group = 2),
    within = c(time = 2, condition = 3)
  )
  means <- design_term_means(
    d,
    term = "group:time",
    target_pes = 0.15,
    n = 16
  )
  sim <- simulate_design_dataset(d, n = 16, means = means, empirical = TRUE)
  fit <- anovapowersim:::fit_design_model(sim, d)
  group_time <- anovapowersim:::extract_term_stats(fit, "group:time")
  time_condition <- anovapowersim:::extract_term_stats(fit, "time:condition")
  three_way <- anovapowersim:::extract_term_stats(fit, "group:time:condition")

  expect_equal(group_time$pes, 0.15, tolerance = 1e-8)
  expect_equal(time_condition$pes, 0, tolerance = 1e-8)
  expect_equal(three_way$pes, 0, tolerance = 1e-8)
})

test_that("default term means calibrate common design classes exactly", {
  check_calibration <- function(between = NULL, within = NULL, term, n,
                                target_pes = 0.08) {
    d <- balanced_anova_design(between = between, within = within)
    resolved_term <- anovapowersim:::resolve_design_term(term, d)
    means <- design_term_means(
      d,
      term = term,
      target_pes = target_pes,
      n = n
    )
    sim <- simulate_design_dataset(d, n = n, means = means, empirical = TRUE)
    fit <- anovapowersim:::fit_design_model(sim, d)
    stats <- anovapowersim:::extract_term_stats(fit, resolved_term)

    expect_equal(stats$pes, target_pes, tolerance = 1e-8)
  }

  check_calibration(between = c(group = 3), term = "group", n = 4)
  check_calibration(
    between = c(group = 2, condition = 3),
    term = "group:condition",
    n = 4
  )
  check_calibration(within = c(time = 4), term = "time", n = 5)
  check_calibration(
    within = c(time = 2, condition = 3),
    term = "time:condition",
    n = 7
  )
  check_calibration(
    between = c(group = 2),
    within = c(time = 3),
    term = "group",
    n = 5
  )
  check_calibration(
    between = c(group = 2),
    within = c(time = 3),
    term = "time",
    n = 5
  )
  check_calibration(
    between = c(group = 2),
    within = c(time = 3),
    term = "group:time",
    n = 5
  )
  check_calibration(
    between = c(group = 2),
    within = c(time = 2, condition = 3),
    term = "group:time:condition",
    n = 7
  )
})

test_that("G*Power option uses total N noncentrality convention", {
  pc <- quiet_power_curve(
    between = c(group = 2),
    within = c(time = 2, condition = 3),
    term = "group:time",
    target_pes = 0.15,
    n_range = 16,
    n_sims = 200,
    gpower = TRUE,
    seed = 10
  )
  expected <- stats::pf(
    stats::qf(0.95, df1 = 1, df2 = 30),
    df1 = 1,
    df2 = 30,
    ncp = 32 * 0.15 / (1 - 0.15),
    lower.tail = FALSE
  )

  expect_true(pc$gpower)
  expect_equal(pc$results$power_calc, round(expected, 3), tolerance = 1e-8)
})

test_that("design_term_means can calibrate to G*Power convention", {
  d <- balanced_anova_design(
    between = c(group = 2),
    within = c(time = 2, condition = 3)
  )
  expect_warning(
    means <- design_term_means(
      d,
      term = "group:time",
      target_pes = 0.15,
      n = 16,
      gpower = TRUE
    ),
    "gpower = TRUE"
  )
  sim <- simulate_design_dataset(d, n = 16, means = means, empirical = TRUE)
  fit <- anovapowersim:::fit_design_model(sim, d)
  group_time <- anovapowersim:::extract_term_stats(fit, "group:time")
  expected_f2 <- (32 * 0.15 / (1 - 0.15)) / 30
  expected_pes <- expected_f2 / (1 + expected_f2)

  expect_equal(group_time$pes, expected_pes, tolerance = 1e-8)
})

test_that("power_curve aligns with noncentral F using denominator df", {
  pc <- quiet_power_curve(
    between = c(group = 2),
    within = c(time = 2, condition = 3),
    term = "group:time",
    target_pes = 0.15,
    n_range = 16,
    n_sims = 1000,
    seed = 9
  )
  expected <- stats::pf(
    stats::qf(0.95, df1 = 1, df2 = 30),
    df1 = 1,
    df2 = 30,
    ncp = 30 * 0.15 / (1 - 0.15),
    lower.tail = FALSE
  )

  expect_equal(pc$results$power_sim, expected, tolerance = 0.05)
})

test_that("design terms are order-insensitive for interactions", {
  d <- balanced_anova_design(
    between = c(group = 2),
    within = c(time = 2, condition = 2)
  )

  means_a <- design_term_means(
    d,
    term = "group:time:condition",
    target_pes = 0.1,
    n = 12
  )
  means_b <- design_term_means(
    d,
    term = "condition:group:time",
    target_pes = 0.1,
    n = 12
  )

  expect_equal(means_a, means_b)

  pc <- quiet_power_curve(
    between = c(group = 2),
    within = c(time = 2),
    term = "time:group",
    target_pes = 0.1,
    n_range = 8,
    n_sims = 10,
    seed = 7
  )
  expect_equal(pc$term, "group:time")
})

test_that("design terms reject unknown and repeated factors", {
  d <- balanced_anova_design(
    between = c(group = 2),
    within = c(time = 2)
  )

  expect_error(
    design_term_means(d, term = "time:missing", target_pes = 0.1, n = 8),
    "Unknown factor"
  )
  expect_error(
    design_term_means(d, term = "time:time", target_pes = 0.1, n = 8),
    "repeat"
  )
})

test_that("power_curve handles pure within and pure between designs", {
  within_pc <- quiet_power_curve(
    within = c(time = 2),
    term = "time",
    target_pes = 0.2,
    n_range = 8,
    n_sims = 10,
    seed = 2
  )
  expect_equal(within_pc$results$total_n, 8L)

  between_pc <- quiet_power_curve(
    between = c(group = 3),
    term = "group",
    target_pes = 0.2,
    n_range = 8,
    n_sims = 10,
    seed = 3
  )
  expect_equal(between_pc$results$total_n, 24L)
})

test_that("power_curve handles common factorial design classes", {
  between_interaction <- quiet_power_curve(
    between = c(group = 2, condition = 2),
    term = "group:condition",
    target_pes = 0.08,
    n_range = 2,
    n_sims = 5,
    seed = 21
  )
  expect_equal(between_interaction$results$total_n, 8L)
  expect_equal(between_interaction$results$num_df, 1)

  within_interaction <- quiet_power_curve(
    within = c(time = 2, condition = 3),
    term = "time:condition",
    target_pes = 0.08,
    n_range = 7,
    n_sims = 5,
    seed = 22
  )
  expect_equal(within_interaction$results$total_n, 7L)
  expect_equal(within_interaction$results$num_df, 2)

  mixed_between_main <- quiet_power_curve(
    between = c(group = 2),
    within = c(time = 3),
    term = "group",
    target_pes = 0.08,
    n_range = 4,
    n_sims = 5,
    seed = 23
  )
  expect_equal(mixed_between_main$results$total_n, 8L)
  expect_equal(mixed_between_main$results$num_df, 1)

  mixed_within_main <- quiet_power_curve(
    between = c(group = 2),
    within = c(time = 3),
    term = "time",
    target_pes = 0.08,
    n_range = 4,
    n_sims = 5,
    seed = 24
  )
  expect_equal(mixed_within_main$results$total_n, 8L)
  expect_equal(mixed_within_main$results$num_df, 2)
})

test_that("power_curve validates design inputs", {
  expect_error(
    quiet_power_curve(term = "x", target_pes = 0.2, n_range = 10),
    "At least one"
  )
  expect_error(
    quiet_power_curve(between = c(2), term = "x", target_pes = 0.2,
                       n_range = 10),
    "named integer"
  )
  expect_error(
    quiet_power_curve(between = c(group = 2), within = c(group = 2),
                       term = "group", target_pes = 0.2, n_range = 10),
    "unique"
  )
  expect_error(
    quiet_power_curve(between = c("group condition" = 2), term = "group condition",
                       target_pes = 0.2, n_range = 10),
    "syntactic R names"
  )
  expect_error(
    quiet_power_curve(within = c("time-point" = 2), term = "time-point",
                       target_pes = 0.2, n_range = 10),
    "syntactic R names"
  )
  expect_error(
    quiet_power_curve(between = c(group = 2), term = "missing",
                       target_pes = 0.2, n_range = 10),
    "Unknown factor"
  )
})

test_that("power_curve is reproducible with a seed", {
  args <- list(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.1,
    n_range = 8,
    n_sims = 15,
    seed = 42
  )

  pc1 <- do.call(quiet_power_curve, args)
  pc2 <- do.call(quiet_power_curve, args)

  expect_equal(pc1$results, pc2$results)
})

test_that("power_curve supports parallel simulations", {
  pc <- quiet_power_curve(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.1,
    n_range = c(8, 10),
    n_sims = 10,
    parallel = TRUE,
    cores = min(2L, as.integer(future::availableCores()[[1L]])),
    seed = 99
  )

  expect_s3_class(pc, "anovapowersim_curve")
  expect_equal(pc$results$n_per_cell, c(8L, 10L))
  expect_true(all(c("n_per_cell", "total_n", "n_sims", "power_sim") %in%
                    names(pc$results)))
  expect_true(all(pc$results$power_sim >= 0 & pc$results$power_sim <= 1))
})

test_that("power_curve validates parallel controls", {
  available <- as.integer(future::availableCores()[[1L]])

  expect_error(
    quiet_power_curve(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.1,
      n_range = 8,
      parallel = NA
    ),
    "`parallel` must be TRUE or FALSE"
  )
  expect_error(
    quiet_power_curve(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.1,
      n_range = 8,
      cores = 0
    ),
    "`cores` must be a single positive integer"
  )
  expect_error(
    quiet_power_curve(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.1,
      n_range = 8,
      cores = available + 1L
    ),
    "`cores` must not exceed"
  )
})

test_that("power_curve reports the default core count in parallel mode", {
  expect_message(
    quiet_power_curve(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.1,
      n_range = 8,
      n_sims = 4,
      parallel = TRUE,
      seed = 42
    ),
    "`parallel = TRUE` and `cores` was not set"
  )
})

test_that("power_curve advises parallel mode for long serial runs", {
  expect_message(
    expect_error(
      quiet_power_curve(
        between = c(group = 2),
        within = c(time = 2),
        term = "group:time",
        target_pes = 0.1,
        n_range = 1,
        n_sims = 5000,
        parallel = FALSE
      ),
      "n_range"
    ),
    "Consider setting `parallel = TRUE`"
  )
})

test_that("power_curve is reproducible with a seed in parallel mode", {
  args <- list(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.1,
    n_range = c(8, 10),
    n_sims = 15,
    parallel = TRUE,
    cores = min(2L, as.integer(future::availableCores()[[1L]])),
    seed = 42
  )

  pc1 <- do.call(quiet_power_curve, args)
  pc2 <- do.call(quiet_power_curve, args)

  expect_equal(pc1$results, pc2$results)
})

test_that("power_curve requires explicit n_range", {
  expect_error(
    quiet_power_curve(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.2,
      n_sims = 10
    ),
    "n_range"
  )
})

test_that("sample sizes below the calibration minimum are rejected clearly", {
  expect_error(
    quiet_power_curve(
      between = c(group = 2, stim = 2),
      term = "group:stim",
      target_pes = 0.08,
      n_range = 1,
      n_sims = 5
    ),
    "n_range.*at least 2"
  )

  expect_error(
    quiet_power_n(
      between = c(group = 2, stim = 2),
      term = "group:stim",
      target_pes = 0.08,
      n_start = 1,
      n_sims = 5
    ),
    "n_start.*at least 2"
  )

  expect_error(
    quiet_power_curve(
      between = c(group = 2),
      within = c(stimulus = 2, condition = 3),
      term = "group:stimulus:condition",
      target_pes = 0.08,
      n_range = 4,
      n_sims = 5
    ),
    "at least 7"
  )

  expect_error(
    quiet_power_n(
      between = c(group = 2),
      within = c(stimulus = 2, condition = 3),
      term = "group:stimulus:condition",
      target_pes = 0.08,
      n_start = 4,
      n_sims = 5
    ),
    "n_start.*at least 7"
  )

  d <- balanced_anova_design(
    between = c(group = 2),
    within = c(stimulus = 2, condition = 3)
  )
  expect_error(
    design_term_means(
      d,
      term = "group:stimulus:condition",
      target_pes = 0.08,
      n = 4
    ),
    "at least 7"
  )
})

test_that("power_n defaults to the calibration minimum", {
  pc <- quiet_power_n(
    between = c(group = 2),
    within = c(stimulus = 2, condition = 3),
    term = "group:stimulus:condition",
    target_pes = 0.08,
    n_sims = 5,
    n_max = 8,
    seed = 11
  )

  expect_true(min(pc$results$n_per_cell) >= 7L)
})

test_that("power_n default n_start avoids saturated between-subject models", {
  pc <- quiet_power_n(
    between = c(group = 2, stim = 2),
    term = "group:stim",
    target_pes = 0.08,
    n_sims = 5,
    n_max = 4,
    seed = 15
  )

  expect_s3_class(pc, "anovapowersim_curve")
  expect_true(min(pc$results$n_per_cell) >= 2L)
})

test_that("power_n default n_start uses the NCP estimate when available", {
  pc_auto <- quiet_power_n(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    power = 0.8,
    n_sims = 10,
    n_max = 80,
    seed = 12
  )
  pc_min <- quiet_power_n(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    power = 0.8,
    n_sims = 10,
    n_start = 3,
    n_max = 80,
    seed = 12
  )

  expect_gt(min(pc_auto$results$n_per_cell), 3L)
  expect_lt(min(pc_auto$results$n_per_cell), max(pc_auto$results$n_per_cell))
  expect_equal(pc_min$results$n_per_cell[[1L]], 3L)
})

test_that("power_n adaptively searches for required n", {
  pc <- quiet_power_n(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    power = 0.8,
    n_sims = 20,
    n_start = 4,
    n_max = 40,
    seed = 8
  )

  expect_s3_class(pc, "anovapowersim_curve")
  expect_true(nrow(pc$results) >= 1L)
  expect_equal(pc$power, 0.8)
  expect_true(is.na(pc$n_needed) || pc$n_needed <= 40L)
  expect_true(is.na(pc$n_needed) || pc$n_needed %in% pc$results$n_per_cell)
  expect_identical(pc$sim_correction, "auto")
  expect_identical(pc$sim_correction_resolved, "none")
})

test_that("power_n is reproducible with a seed", {
  args <- list(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.1,
    power = 0.8,
    n_sims = 15,
    n_start = 4,
    n_max = 30,
    seed = 42
  )

  pc1 <- do.call(quiet_power_n, args)
  pc2 <- do.call(quiet_power_n, args)

  expect_equal(pc1$results, pc2$results)
  expect_equal(pc1$n_needed, pc2$n_needed)
})

test_that("power_n supports parallel simulations within each searched n", {
  pc <- quiet_power_n(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.1,
    power = 0.8,
    n_sims = 6,
    n_start = 4,
    n_max = 8,
    parallel = TRUE,
    cores = min(2L, as.integer(future::availableCores()[[1L]])),
    seed = 101
  )

  expect_s3_class(pc, "anovapowersim_curve")
  expect_true(nrow(pc$results) >= 1L)
  expect_true(all(pc$results$power_sim >= 0 & pc$results$power_sim <= 1))
})

test_that("parallel power_n returns finite simulation power for mixed interactions", {
  pc <- quiet_power_n(
    between = c(group = 2),
    within = c(stim = 2, cond = 3),
    term = "group:stim:cond",
    target_pes = 0.08,
    power = 0.80,
    n_sims = 10,
    n_max = 30,
    parallel = TRUE,
    cores = min(2L, as.integer(future::availableCores()[[1L]])),
    seed = 123
  )

  expect_true(all(!is.na(pc$results$power_sim)))
  expect_true(all(pc$results$power_sim >= 0 & pc$results$power_sim <= 1))
})

test_that("power_n advises parallel mode for long serial runs", {
  expect_message(
    expect_error(
      quiet_power_n(
        between = c(group = 2),
        within = c(time = 2),
        term = "group:time",
        target_pes = 0.1,
        n_sims = 5000,
        n_start = 3,
        n_max = 1,
        parallel = FALSE
      ),
      "n_max"
    ),
    "Consider setting `parallel = TRUE`"
  )
})

test_that("power_n validates parallel controls", {
  available <- as.integer(future::availableCores()[[1L]])

  expect_error(
    quiet_power_n(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.1,
      parallel = c(TRUE, FALSE)
    ),
    "`parallel` must be TRUE or FALSE"
  )
  expect_error(
    quiet_power_n(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.1,
      cores = 1.5
    ),
    "`cores` must be a single positive integer"
  )
  expect_error(
    quiet_power_n(
      between = c(group = 2),
      within = c(time = 2),
      term = "group:time",
      target_pes = 0.1,
      cores = available + 1L
    ),
    "`cores` must not exceed"
  )
})

test_that("power_n rejects a starting sample size above n_max", {
  expect_error(
    quiet_power_n(
      between = c(group = 2),
      term = "group",
      target_pes = 0.1,
      n_sims = 1,
      n_start = 10,
      n_max = 9,
      progress = FALSE
    ),
    "`n_max` must be greater than or equal to `n_start`",
    fixed = TRUE
  )
})

test_that("balanced simulation power retains full precision", {
  result <- quiet_power_curve(
    between = c(group = 2),
    term = "group",
    target_pes = 0.1,
    n_range = 8,
    n_sims = 7,
    seed = 808,
    progress = FALSE
  )$results

  expect_equal(result$valid_sims + result$failed_sims, result$n_sims)
  expect_equal(result$power_sim * result$valid_sims,
               round(result$power_sim * result$valid_sims))
  expect_false(isTRUE(all.equal(result$power_sim, round(result$power_sim, 3))))
})

test_that("Wilson intervals cover the observed simulation proportion", {
  interval <- anovapowersim:::binomial_wilson_interval(5, 10)

  expect_equal(unname(interval), c(0.2365931, 0.7634069), tolerance = 1e-7)
  expect_lte(interval[["lower"]], 0.5)
  expect_gte(interval[["upper"]], 0.5)
})

test_that("warn_gpower_within_term_df warns whenever gpower = TRUE", {
  expect_warning(
    anovapowersim:::warn_gpower_within_term_df(TRUE),
    "gpower = TRUE"
  )
  expect_silent(
    anovapowersim:::warn_gpower_within_term_df(FALSE)
  )
})

test_that("gpower = TRUE warns across all gpower-accepting functions", {
  # each incidental warning (target power not reached with a tiny n_max,
  # power < .90 recommended, sim/calc disagreement from n_sims = 1 noise) is
  # an expected side effect of keeping these calls fast; capture_warnings()
  # checks only for the gpower-specific message among whatever fires.
  expect_gpower_warning <- function(expr) {
    warnings <- testthat::capture_warnings(expr)
    expect_true(any(grepl("gpower = TRUE", warnings, fixed = TRUE)))
  }
  expect_no_gpower_warning <- function(expr) {
    warnings <- testthat::capture_warnings(expr)
    expect_false(any(grepl("gpower = TRUE", warnings, fixed = TRUE)))
  }

  # calculated-power functions (fast, no simulation)
  expect_gpower_warning(
    power_n_calc(
      within = c(time = 2), term = "time", target_pes = 0.1, power = 0.5,
      n_start = 4, n_max = 4, gpower = TRUE
    )
  )
  expect_no_gpower_warning(
    power_n_calc(
      within = c(time = 2), term = "time", target_pes = 0.1, power = 0.5,
      n_start = 4, n_max = 4, gpower = FALSE
    )
  )

  expect_gpower_warning(
    power_achieved_calc(
      within = c(time = 3), term = "time", target_pes = 0.1, n = 4,
      gpower = TRUE
    )
  )
  expect_gpower_warning(
    power_sensitivity_calc(
      within = c(time = 3), term = "time", n = 4, power = 0.9,
      pes_tol = 0.3, gpower = TRUE
    )
  )

  expect_gpower_warning(
    design_term_means(
      balanced_anova_design(within = c(time = 2)), term = "time",
      target_pes = 0.1, n = 4, gpower = TRUE
    )
  )
  expect_no_gpower_warning(
    design_term_means(
      balanced_anova_design(within = c(time = 2)), term = "time",
      target_pes = 0.1, n = 4, gpower = FALSE
    )
  )

  # simulation-based functions
  expect_gpower_warning(
    power_curve(
      within = c(time = 2), term = "time", target_pes = 0.1, n_range = 5,
      n_sims = 1, gpower = TRUE, progress = FALSE
    )
  )
  expect_no_gpower_warning(
    power_curve(
      within = c(time = 2), term = "time", target_pes = 0.1, n_range = 5,
      n_sims = 1, gpower = FALSE, progress = FALSE
    )
  )

  expect_gpower_warning(
    power_n(
      within = c(time = 3), term = "time", target_pes = 0.1, power = 0.5,
      n_sims = 1, n_start = 4, n_max = 4, gpower = TRUE, progress = FALSE
    )
  )

  expect_gpower_warning(
    power_achieved(
      within = c(time = 3), term = "time", target_pes = 0.1, n = 4,
      n_sims = 1, gpower = TRUE, progress = FALSE
    )
  )

  expect_gpower_warning(
    power_sensitivity(
      within = c(time = 3), term = "time", n = 4, power = 0.5,
      n_sims = 1, pes_tol = 0.5, gpower = TRUE, progress = FALSE
    )
  )
})
