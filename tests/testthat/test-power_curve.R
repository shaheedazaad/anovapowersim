quiet_power_curve <- function(...) suppressWarnings(power_curve(...))
quiet_power_n <- function(...) suppressWarnings(power_n(...))

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
  expect_true(all(c("num_df", "den_df", "ncp", "power_calc", "power_sim") %in%
                    names(pc$results)))
  hidden_columns <- paste0(c("ci", "ci", "n_success", "n_fail"),
                           c("_low", "_high", "es", "ed"))
  expect_false(any(hidden_columns %in% names(pc$results)))
  expect_true(all(pc$results$power_sim >= 0 & pc$results$power_sim <= 1))
  expect_equal(pc$target_pes, 0.20721)
  expect_equal(pc$ss_type, "III")
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
  means <- design_term_means(
    d,
    term = "group:time",
    target_pes = 0.15,
    n = 16,
    gpower = TRUE
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
