test_that("means_pattern parses and retains sparse cell definitions", {
  pattern <- means_pattern(
    time = 1, value = 0,
    time = 3, value = 1
  )

  expect_s3_class(pattern, "anovapowersim_means_pattern")
  expect_identical(pattern$factor_names, "time")
  expect_length(pattern$definitions, 2L)
  expect_equal(pattern$definitions[[2L]]$value, 1)
})

test_that("means_pattern validates its sparse syntax", {
  expect_error(means_pattern(), "at least one sparse cell")
  expect_error(means_pattern(1, value = 2), "must be named")
  expect_error(means_pattern(value = 1), "must follow")
  expect_error(means_pattern(time = 1), "missing `value`")
  expect_error(means_pattern(time = 1, time = 2, value = 1),
               "more than once")
  expect_error(
    means_pattern(time = 1, value = 0, group = 1, value = 1),
    "same factor names"
  )
  expect_error(means_pattern(time = 1, value = Inf), "finite numeric")
  expect_error(means_pattern(time = 1, value = "high"), "finite numeric")
  expect_error(
    means_pattern(time = 1, value = 0, time = 1, value = 1),
    "only once"
  )
})

test_that("custom pattern semantics are explained once per session", {
  design <- balanced_anova_design(within = c(time = 2))
  pattern <- means_pattern(time = 1, value = -10, time = 2, value = 10)

  expect_message(
    design_term_means(
      design, term = "time", target_pes = 0.15, n = 4,
      means_pattern = pattern
    ),
    "shape-only.*term 'time'.*target_pes = 0.15.*literal `m` values"
  )
  expect_silent(design_term_means(
    design, term = "time", target_pes = 0.20, n = 4,
    means_pattern = pattern
  ))
})

test_that("sparse patterns fill zero cells, broadcast, and preserve ordering", {
  within_design <- balanced_anova_design(within = c(time = 4))
  sparse <- means_pattern(time = 3, value = 1)
  resolved <- anovapowersim:::resolve_means_pattern(
    sparse, within_design, "time"
  )
  expected <- c(-0.25, -0.25, 0.75, -0.25)
  expected <- expected / sqrt(sum(expected^2))
  expect_equal(as.numeric(resolved), expected, tolerance = 1e-12)

  mixed_design <- balanced_anova_design(
    between = c(group = 2), within = c(time = 4)
  )
  broadcast <- anovapowersim:::resolve_means_pattern(
    sparse, mixed_design, "time"
  )
  expect_equal(broadcast[1, ], broadcast[2, ], tolerance = 1e-12)
  expect_equal(sqrt(sum(broadcast^2)), 1, tolerance = 1e-12)

  interaction <- means_pattern(
    group = 1, time = 3, value = 1,
    group = 2, time = 3, value = -1
  )
  interaction_resolved <- anovapowersim:::resolve_means_pattern(
    interaction, mixed_design, "group:time"
  )
  expect_equal(dim(interaction_resolved), c(2L, 4L))
  expect_equal(colSums(interaction_resolved), rep(0, 4), tolerance = 1e-12)
  expect_equal(rowSums(interaction_resolved), c(0, 0), tolerance = 1e-12)
  expect_gt(abs(interaction_resolved[1, 3]),
            abs(interaction_resolved[1, 1]))
})

test_that("pattern levels accept indices and generated names canonically", {
  design <- balanced_anova_design(within = c(time = 4))
  by_index <- anovapowersim:::resolve_means_pattern(
    means_pattern(time = 3, value = 1), design, "time"
  )
  by_name <- anovapowersim:::resolve_means_pattern(
    means_pattern(time = "time3", value = 1), design, "time"
  )
  expect_equal(by_index, by_name)

  mixed_duplicate <- means_pattern(
    time = 3, value = 1,
    time = "time3", value = 2
  )
  expect_error(
    anovapowersim:::resolve_means_pattern(mixed_duplicate, design, "time"),
    "generated level names.*duplicates"
  )
})

test_that("pattern resolution rejects invalid factors and levels clearly", {
  design <- balanced_anova_design(
    between = c(group = 2), within = c(time = 4)
  )
  expect_error(
    anovapowersim:::resolve_means_pattern(
      means_pattern(group = 1, value = 1), design, "time"
    ),
    "not part of term"
  )
  expect_error(
    anovapowersim:::resolve_means_pattern(
      means_pattern(missing = 1, value = 1), design, "time"
    ),
    "not part of term"
  )
  for (bad in list(0, 5, 1.5, "third", TRUE)) {
    pattern <- means_pattern(time = bad, value = 1)
    expect_error(
      anovapowersim:::resolve_means_pattern(pattern, design, "time"),
      "indices 1 through 4.*'time1'.*'time4'"
    )
  }
  expect_error(
    anovapowersim:::resolve_means_pattern(list(time = 1), design, "time"),
    "created by means_pattern"
  )
})

test_that("projection rejects zero target components and removes nuisances", {
  design <- balanced_anova_design(
    between = c(group = 2), within = c(time = 3)
  )
  no_interaction <- means_pattern(
    group = 1, value = -1,
    group = 2, value = 1
  )
  expect_error(
    anovapowersim:::resolve_means_pattern(
      no_interaction, design, "group:time"
    ),
    "zero projection"
  )

  sparse_interaction <- means_pattern(
    group = 1, time = 3, value = 1,
    group = 2, time = 3, value = -1
  )
  resolved <- anovapowersim:::resolve_means_pattern(
    sparse_interaction, design, "group:time"
  )
  expect_equal(rowSums(resolved), c(0, 0), tolerance = 1e-12)
  expect_equal(colSums(resolved), rep(0, 3), tolerance = 1e-12)
  expect_equal(sqrt(sum(resolved^2)), 1, tolerance = 1e-12)
})

test_that("default patterns use normalized centered linear Kronecker scores", {
  expected_scores <- list(
    c(-1, 1),
    c(-1, 0, 1),
    c(-3, -1, 1, 3),
    c(-2, -1, 0, 1, 2)
  )
  for (n_levels in 2:5) {
    scores <- anovapowersim:::centered_linear_scores(n_levels)
    expected <- expected_scores[[n_levels - 1L]]
    expect_equal(scores / scores[[length(scores)]],
                 expected / expected[[length(expected)]])
    expect_equal(sum(scores), 0, tolerance = 1e-12)
    expect_equal(sqrt(sum(scores^2)), 1, tolerance = 1e-12)
  }

  two_by_two <- balanced_anova_design(
    between = c(group = 2), within = c(time = 2)
  )
  default_22 <- anovapowersim:::default_term_pattern(
    two_by_two, "group:time"
  )
  expected_22 <- matrix(c(1, -1, -1, 1), 2, 2, byrow = TRUE)
  expect_equal(default_22 / default_22[1, 1], expected_22)

  two_by_four <- balanced_anova_design(
    between = c(group = 2), within = c(time = 4)
  )
  default_24 <- anovapowersim:::default_term_pattern(
    two_by_four, "group:time"
  )
  expected_24 <- rbind(c(3, 1, -1, -3), c(-3, -1, 1, 3))
  expect_equal(default_24 / default_24[1, 2], expected_24)
  expect_equal(sqrt(sum(default_24^2)), 1, tolerance = 1e-12)

  higher <- balanced_anova_design(
    between = c(group = 2), within = c(time = 3, condition = 2)
  )
  higher_default <- anovapowersim:::default_term_pattern(
    higher, "group:time:condition"
  )
  expected_vector <- as.vector(outer(
    anovapowersim:::centered_linear_scores(2),
    as.vector(outer(
      anovapowersim:::centered_linear_scores(3),
      anovapowersim:::centered_linear_scores(2)
    ))
  ))
  expect_equal(as.numeric(t(higher_default)), expected_vector,
               tolerance = 1e-12)
})

test_that("projection equivalences and calibration scale are invariant", {
  design <- balanced_anova_design(within = c(time = 4))
  base <- means_pattern(
    time = 1, value = -1,
    time = 2, value = 0,
    time = 3, value = 0.2,
    time = 4, value = 0.8
  )
  scaled <- means_pattern(
    time = 1, value = -3,
    time = 2, value = 0,
    time = 3, value = 0.6,
    time = 4, value = 2.4
  )
  shifted <- means_pattern(
    time = 1, value = 4,
    time = 2, value = 5,
    time = 3, value = 5.2,
    time = 4, value = 5.8
  )
  reversed <- means_pattern(
    time = 1, value = 1,
    time = 2, value = 0,
    time = 3, value = -0.2,
    time = 4, value = -0.8
  )

  resolved <- lapply(list(base, scaled, shifted, reversed), function(x) {
    anovapowersim:::resolve_means_pattern(x, design, "time")
  })
  expect_equal(resolved[[1]], resolved[[2]], tolerance = 1e-12)
  expect_equal(resolved[[1]], resolved[[3]], tolerance = 1e-12)
  expect_equal(resolved[[1]], -resolved[[4]], tolerance = 1e-12)

  direction_a <- resolved[[1]]
  direction_b <- anovapowersim:::resolve_means_pattern(
    means_pattern(time = 2, value = 1), design, "time"
  )
  sigma <- matrix(c(
    1, .7, .2, 0,
    .7, 1, .1, 0,
    .2, .1, 1, .4,
    0, 0, .4, 1
  ), 4, 4, byrow = TRUE)
  calibrated_a <- anovapowersim:::calibrate_design_means(
    design, "time", .15, 8, 1, .5, covariance = sigma,
    resolved_means_pattern = direction_a
  )
  calibrated_b <- anovapowersim:::calibrate_design_means(
    design, "time", .15, 8, 1, .5, covariance = sigma,
    resolved_means_pattern = direction_b
  )
  expect_equal(sqrt(sum(calibrated_a^2)), sqrt(sum(calibrated_b^2)),
               tolerance = 1e-10)
})

test_that("design_term_means calibrates supplied patterns", {
  design <- balanced_anova_design(within = c(time = 4))
  pattern <- means_pattern(time = 2, value = 1)
  means <- design_term_means(
    design, "time", target_pes = .15, n = 8,
    means_pattern = pattern
  )
  resolved <- anovapowersim:::resolve_means_pattern(pattern, design, "time")
  expect_equal(
    as.numeric(means) / sqrt(sum(means^2)),
    as.numeric(resolved),
    tolerance = 1e-12
  )
  exact <- simulate_design_dataset(design, 8, means, empirical = TRUE)
  achieved <- anovapowersim:::fit_design_term_stats(
    exact, design, "time"
  )$pes
  expect_equal(achieved, .15, tolerance = 1e-8)
})

test_that("direction warning is limited to implicit nonspherical multi-df simulations", {
  sigma3 <- matrix(c(1, .8, 0, .8, 1, 0, 0, 0, 1), 3, 3)
  dimnames(sigma3) <- list(paste0("time", 1:3), paste0("time", 1:3))
  design3 <- balanced_anova_design(within = c(time = 3))
  epsilon3 <- anovapowersim:::covariance_term_epsilon(sigma3, design3, "time")

  expect_warning(
    anovapowersim:::warn_direction_sensitivity(
      design3, "time", epsilon3, FALSE
    ),
    "power_sim.*default linear/Kronecker.*Supply `means_pattern`"
  )
  expect_silent(anovapowersim:::warn_direction_sensitivity(
    design3, "time", epsilon3, TRUE
  ))
  expect_silent(anovapowersim:::warn_direction_sensitivity(
    balanced_anova_design(within = c(time = 2)), "time", .75, FALSE
  ))
  expect_silent(anovapowersim:::warn_direction_sensitivity(
    balanced_anova_design(between = c(group = 3)), "group", .75, FALSE
  ))
  expect_silent(anovapowersim:::warn_direction_sensitivity(
    design3, "time", 1, FALSE
  ))

  custom <- means_pattern(time = 1, value = -1, time = 3, value = 1)
  warnings <- testthat::capture_warnings(
    result <- power_curve(
      within = c(time = 3), term = "time", target_pes = .1,
      n_range = 4, n_sims = 1,
      covariance = test_covariance_spec_from_matrix(sigma3),
      means_pattern = custom, progress = FALSE, seed = 1
    )
  )
  expect_false(any(grepl("relative cell-mean pattern", warnings, fixed = TRUE)))
  expect_true(result$custom_means_pattern)

  type_i_warnings <- testthat::capture_warnings(power_curve(
    within = c(time = 3), term = "time", target_pes = .1,
    n_range = 4, n_sims = 1,
    covariance = test_covariance_spec_from_matrix(sigma3),
    ss_type = "I", progress = FALSE, seed = 1
  ))
  expect_true(any(grepl("cannot provide GG-corrected", type_i_warnings,
                        fixed = TRUE)))
  expect_true(any(grepl("relative cell-mean pattern", type_i_warnings,
                        fixed = TRUE)))
})

test_that("balanced simulation APIs retain custom-pattern metadata", {
  pattern <- means_pattern(time = 1, value = -1, time = 2, value = 1)
  curve <- suppressWarnings(power_curve(
    within = c(time = 2), term = "time", target_pes = .1,
    n_range = 4, n_sims = 1, means_pattern = pattern,
    progress = FALSE, seed = 1
  ))
  achieved <- suppressWarnings(power_achieved(
    within = c(time = 2), term = "time", target_pes = .1,
    n = 4, n_sims = 1, means_pattern = pattern,
    progress = FALSE, seed = 1
  ))
  searched_n <- suppressWarnings(power_n(
    within = c(time = 2), term = "time", target_pes = .1,
    power = .5, n_start = 4, n_max = 4, n_sims = 1,
    means_pattern = pattern, progress = FALSE, seed = 1
  ))
  sensitivity <- suppressWarnings(power_sensitivity(
    within = c(time = 2), term = "time", n = 4, power = .5,
    pes_min = .01, pes_max = .5, pes_tol = .5, n_sims = 1,
    means_pattern = pattern, progress = FALSE, seed = 1
  ))

  expect_true(curve$custom_means_pattern)
  expect_true(achieved$custom_means_pattern)
  expect_true(searched_n$custom_means_pattern)
  expect_true(sensitivity$custom_means_pattern)
  expect_output(print(curve), "means pattern: custom", fixed = TRUE)
  summary_result <- NULL
  expect_output(
    summary_result <- summary(curve),
    "means pattern: custom",
    fixed = TRUE
  )
  expect_match(summary_result$header[["means pattern"]], "custom")
})
