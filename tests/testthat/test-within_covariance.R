quiet_covariance_power_curve <- function(...) suppressWarnings(power_curve(...))
quiet_covariance_power_n <- function(...) suppressWarnings(power_n(...))

test_that("within_covariance records flat SD and correlation overrides", {
  covariance <- within_covariance(
    default_sd = 1,
    default_correlation = 0.5,
    standard_deviations = c("time2_condition2" = 1.2),
    correlations = c(
      "time1_condition1:time1_condition2" = 0.6,
      "time2_condition1:time2_condition2" = 0.7
    )
  )

  expect_s3_class(covariance, "anovapowersim_covariance_spec")
  expect_equal(covariance$default_sd, 1)
  expect_equal(covariance$default_correlation, 0.5)
  expect_equal(covariance$standard_deviations[["time2_condition2"]], 1.2)
  expect_equal(
    covariance$correlation_pairs$pair_name,
    c(
      "time1_condition1:time1_condition2",
      "time2_condition1:time2_condition2"
    )
  )
})

test_that("covariance specifications resolve defaults and overrides", {
  design <- balanced_anova_design(within = c(time = 2, condition = 2))
  covariance <- within_covariance(
    default_sd = 1,
    default_correlation = 0.5,
    standard_deviations = c("time2_condition2" = 1.2),
    correlations = c(
      "time1_condition1:time1_condition2" = 0.6,
      "time2_condition1:time2_condition2" = 0.7
    )
  )

  sigma <- anovapowersim:::resolve_within_covariance(covariance, design)
  expected_names <- c(
    "time1_condition1", "time1_condition2",
    "time2_condition1", "time2_condition2"
  )

  expect_equal(rownames(sigma), expected_names)
  expect_equal(colnames(sigma), expected_names)
  expect_equal(unname(diag(sigma)), c(1, 1, 1, 1.44), tolerance = 1e-12)
  expect_equal(sigma["time1_condition1", "time1_condition2"], 0.6)
  expect_equal(sigma["time1_condition1", "time2_condition1"], 0.5)
  expect_equal(sigma["time1_condition1", "time2_condition2"], 0.6)
  expect_equal(sigma["time2_condition1", "time2_condition2"], 0.84)
  expect_equal(sigma, t(sigma))
})

test_that("correlation pair order is ignored and duplicate pairs are rejected", {
  design <- balanced_anova_design(within = c(time = 3))
  forward <- within_covariance(
    correlations = c("time1:time2" = 0.7)
  )
  reverse <- within_covariance(
    correlations = c("time2:time1" = 0.7)
  )

  expect_equal(
    anovapowersim:::resolve_within_covariance(forward, design),
    anovapowersim:::resolve_within_covariance(reverse, design)
  )
  expect_error(
    within_covariance(
      correlations = c(
        "time1:time2" = 0.7,
        "time2:time1" = 0.6
      )
    ),
    "defined only once"
  )
})

test_that("within_covariance validates its public inputs", {
  expect_error(within_covariance(default_sd = 0), "`default_sd`")
  expect_error(within_covariance(default_correlation = 1),
               "`default_correlation`")
  expect_error(
    within_covariance(standard_deviations = c(time1 = -1)),
    "positive finite numbers"
  )
  expect_error(
    within_covariance(standard_deviations = c(1)),
    "must be named"
  )
  expect_error(
    within_covariance(correlations = c("time1:time2" = 1)),
    "finite numbers"
  )
  expect_error(
    within_covariance(correlations = c("time1-time2" = 0.5)),
    "cell1:cell2"
  )
  expect_error(
    within_covariance(correlations = c("time1:time1" = 0.5)),
    "itself"
  )
})

test_that("covariance resolution validates cells and positive definiteness", {
  design <- balanced_anova_design(within = c(time = 3))

  expect_error(
    anovapowersim:::resolve_within_covariance(
      within_covariance(standard_deviations = c("missing" = 1)),
      design
    ),
    "Unknown within-subject cell"
  )
  expect_error(
    anovapowersim:::resolve_within_covariance(
      within_covariance(
        default_correlation = 0,
        correlations = c(
          "time1:time2" = 0.9,
          "time1:time3" = 0.9,
          "time2:time3" = -0.9
        )
      ),
      design
    ),
    "positive definite"
  )
  between_design <- balanced_anova_design(between = c(group = 2))
  expect_error(
    anovapowersim:::resolve_within_covariance(within_covariance(), between_design),
    "design with a `within` factor"
  )
})

test_that("named covariance matrices are reordered to the design", {
  design <- balanced_anova_design(within = c(time = 3))
  sigma <- matrix(
    c(1, 0.3, 0.2, 0.3, 1, 0.4, 0.2, 0.4, 1),
    nrow = 3,
    dimnames = list(
      c("time3", "time1", "time2"),
      c("time3", "time1", "time2")
    )
  )
  resolved <- anovapowersim:::resolve_within_covariance(sigma, design)

  expect_equal(rownames(resolved), c("time1", "time2", "time3"))
  expect_equal(resolved["time1", "time2"], 0.4)
  expect_equal(resolved["time1", "time3"], 0.3)
})

test_that("power_curve and power_n accept and retain custom covariance", {
  covariance <- within_covariance(
    default_sd = 1,
    default_correlation = 0.4,
    correlations = c("time1:time2" = 0.7)
  )
  curve <- quiet_covariance_power_curve(
    within = c(time = 3),
    term = "time",
    target_pes = 0.14,
    n_range = 5,
    n_sims = 3,
    covariance = covariance,
    progress = FALSE,
    seed = 101
  )
  searched <- quiet_covariance_power_n(
    within = c(time = 3),
    term = "time",
    target_pes = 0.14,
    power = 0.99,
    n_sims = 3,
    n_start = 5,
    n_max = 5,
    covariance = covariance,
    progress = FALSE,
    seed = 101
  )

  for (result in list(curve, searched)) {
    expect_s3_class(result, "anovapowersim_curve")
    expect_true(result$custom_covariance)
    expect_equal(result$covariance["time1", "time2"], 0.7)
    expect_equal(result$covariance["time1", "time3"], 0.4)
    expect_true(is.finite(result$results$power_sim))
  }
})

test_that("the balanced simulator uses the resolved covariance matrix", {
  design <- balanced_anova_design(within = c(time = 3))
  sigma <- matrix(
    c(1, 0.6, 0.3, 0.6, 1.44, 0.48, 0.3, 0.48, 1),
    nrow = 3,
    dimnames = list(
      c("time1", "time2", "time3"),
      c("time1", "time2", "time3")
    )
  )
  simulated <- anovapowersim:::simulate_balanced_design_data(
    spec = design,
    n = 5,
    means = matrix(0, nrow = 1, ncol = 3),
    sd = 1,
    r = 0.5,
    covariance = sigma,
    empirical = TRUE
  )
  wide <- tidyr::pivot_wider(
    simulated,
    id_cols = "id",
    names_from = "time",
    values_from = "value"
  )
  observed <- stats::cov(as.data.frame(wide[c("time1", "time2", "time3")]))

  expect_equal(observed, sigma, tolerance = 1e-10)
})

test_that("omitting covariance preserves the current default matrix", {
  curve <- quiet_covariance_power_curve(
    within = c(time = 3),
    term = "time",
    target_pes = 0.14,
    n_range = 5,
    n_sims = 1,
    progress = FALSE,
    seed = 9
  )

  expect_false(curve$custom_covariance)
  expect_equal(unname(diag(curve$covariance)), rep(1, 3))
  expect_equal(curve$covariance[lower.tri(curve$covariance)], rep(0.5, 3))
  expect_equal(curve$epsilon, 1)
  expect_equal(curve$results$epsilon, 1)
})

test_that("population epsilon is derived for the tested within term", {
  design <- balanced_anova_design(
    between = c(group = 2),
    within = c(time = 3, condition = 2)
  )
  time_covariance <- diag(c(1, 1, 4))
  condition_covariance <- diag(2)
  sigma <- kronecker(time_covariance, condition_covariance)
  dimnames(sigma) <- list(
    anovapowersim:::within_cell_names(design),
    anovapowersim:::within_cell_names(design)
  )

  expect_equal(
    anovapowersim:::covariance_term_epsilon(sigma, design, "time"),
    0.8,
    tolerance = 1e-12
  )
  expect_equal(
    anovapowersim:::covariance_term_epsilon(
      sigma, design, "time:condition"
    ),
    0.8,
    tolerance = 1e-12
  )
  expect_equal(
    anovapowersim:::covariance_term_epsilon(sigma, design, "group:time"),
    0.8,
    tolerance = 1e-12
  )
  expect_equal(
    anovapowersim:::covariance_term_epsilon(sigma, design, "condition"),
    1
  )
  expect_equal(
    anovapowersim:::covariance_term_epsilon(sigma, design, "group"),
    1
  )
})

test_that("covariance-corrected power_calc matches power_n_calc", {
  sigma <- diag(c(1, 1, 4))
  dimnames(sigma) <- list(paste0("time", 1:3), paste0("time", 1:3))

  curve <- quiet_covariance_power_curve(
    within = c(time = 3),
    term = "time",
    target_pes = 0.15,
    n_range = 8,
    n_sims = 2,
    covariance = sigma,
    progress = FALSE,
    seed = 42
  )
  calculated <- suppressWarnings(power_n_calc(
    within = c(time = 3),
    term = "time",
    target_pes = 0.15,
    n_start = 8,
    n_max = 8,
    epsilon = curve$epsilon
  ))

  expect_equal(curve$epsilon, 0.8, tolerance = 1e-12)
  expect_equal(curve$results$epsilon, calculated$results$epsilon)
  expect_equal(curve$results$num_df, calculated$results$num_df)
  expect_equal(curve$results$den_df, calculated$results$den_df)
  expect_equal(curve$results$ncp, round(calculated$results$ncp, 3))
  expect_equal(
    curve$results$power_calc,
    round(calculated$results$power_calc, 3)
  )
})

test_that("power_n retains and applies covariance-derived epsilon", {
  sigma <- diag(c(1, 1, 4))
  dimnames(sigma) <- list(paste0("time", 1:3), paste0("time", 1:3))

  searched <- quiet_covariance_power_n(
    within = c(time = 3),
    term = "time",
    target_pes = 0.15,
    power = 0.99,
    n_sims = 2,
    n_start = 8,
    n_max = 8,
    covariance = sigma,
    progress = FALSE,
    seed = 42
  )

  expect_equal(searched$epsilon, 0.8, tolerance = 1e-12)
  expect_equal(searched$results$epsilon, 0.8, tolerance = 1e-12)
  expect_equal(searched$results$num_df, 1.6, tolerance = 1e-12)
  expect_equal(searched$results$den_df, 11.2, tolerance = 1e-12)
})

test_that("fit_design_term_stats reports a Greenhouse-Geisser p-value for ss_type III/II that matches car directly, and NA for ss_type I", {
  design <- balanced_anova_design(within = c(time = 3))
  sigma <- diag(c(1, 1, 4))
  dimnames(sigma) <- list(paste0("time", 1:3), paste0("time", 1:3))

  set.seed(321)
  sim <- anovapowersim:::simulate_balanced_design_data(
    spec = design,
    n = 15,
    means = matrix(0, nrow = 1, ncol = 3),
    sd = 1,
    r = 0.5,
    covariance = sigma,
    empirical = FALSE
  )

  stats_iii <- anovapowersim:::fit_design_term_stats(sim, design, "time", ss_type = "III")
  stats_i <- anovapowersim:::fit_design_term_stats(sim, design, "time", ss_type = "I")

  wide <- tidyr::pivot_wider(
    sim,
    id_cols = "id", names_from = "time", values_from = "value"
  )
  fit <- stats::lm(cbind(time1, time2, time3) ~ 1, data = wide)
  idata <- data.frame(time = factor(paste0("t", 1:3)))
  stats::contrasts(idata$time) <- stats::contr.sum(3)
  av <- car::Anova(
    fit, idata = idata, idesign = ~time, type = 3,
    icontrasts = c("contr.sum", "contr.poly")
  )
  expected_gg <- suppressWarnings(
    summary(av, multivariate = FALSE)
  )$pval.adjustments["time", "Pr(>F[GG])"]

  expect_equal(stats_iii$p_value_gg, unname(expected_gg), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(stats_iii$p_value_gg, stats_iii$p_value)))
  expect_true(is.na(stats_i$p_value_gg))
})

test_that("power_curve GG-corrects power_sim to match power_calc under severe non-sphericity", {
  sigma <- diag(c(1, 1, 1, 16))
  dimnames(sigma) <- list(paste0("time", 1:4), paste0("time", 1:4))

  corrected <- quiet_covariance_power_curve(
    within = c(time = 4),
    term = "time",
    target_pes = 0.15,
    n_range = 25,
    n_sims = 1500,
    covariance = sigma,
    progress = FALSE,
    seed = 55
  )
  uncorrected <- quiet_covariance_power_curve(
    within = c(time = 4),
    term = "time",
    target_pes = 0.15,
    n_range = 25,
    n_sims = 1500,
    covariance = sigma,
    ss_type = "I",
    progress = FALSE,
    seed = 55
  )

  expect_equal(corrected$epsilon, 0.4451295, tolerance = 1e-6)
  expect_equal(
    corrected$results$power_sim,
    corrected$results$power_calc,
    tolerance = 0.08
  )
  # ss_type = "I" cannot supply a Greenhouse-Geisser p-value, so power_sim
  # stays uncorrected and diverges noticeably from power_calc.
  expect_true(
    abs(uncorrected$results$power_sim - uncorrected$results$power_calc) > 0.05
  )
})

test_that("warn_ss_type_i_uncorrected_gg only warns for ss_type = 'I' with epsilon < 1", {
  expect_warning(
    anovapowersim:::warn_ss_type_i_uncorrected_gg(ss_type = "I", epsilon = 0.8),
    "ss_type = \"I\""
  )
  expect_silent(
    anovapowersim:::warn_ss_type_i_uncorrected_gg(ss_type = "I", epsilon = 1)
  )
  expect_silent(
    anovapowersim:::warn_ss_type_i_uncorrected_gg(ss_type = "III", epsilon = 0.5)
  )
})

test_that("power_curve warns when ss_type = 'I' is combined with a non-spherical covariance", {
  sigma <- diag(c(1, 1, 4))
  dimnames(sigma) <- list(paste0("time", 1:3), paste0("time", 1:3))

  # ss_type = "I" also diverges from power_calc in a real way here (it stays
  # uncorrected), so the pre-existing power_sim/power_calc disagreement
  # warning legitimately co-fires alongside the new one; capture all
  # warnings rather than asserting on a single one.
  warnings <- testthat::capture_warnings(
    power_curve(
      within = c(time = 3),
      term = "time",
      target_pes = 0.14,
      n_range = 5,
      n_sims = 1,
      covariance = sigma,
      ss_type = "I",
      progress = FALSE,
      seed = 1
    )
  )
  expect_true(any(grepl("ss_type = \"I\"", warnings, fixed = TRUE)))
})
