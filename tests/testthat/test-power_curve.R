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
