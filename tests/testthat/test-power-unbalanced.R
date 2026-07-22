quiet_power_unbalanced <- function(...) suppressWarnings(power_unbalanced(...))


test_that("cell_design creates and validates a complete factorial cell table", {
  design <- cell_design(
    group = "a", time = "pre",  n = 6, m = 1.0, sd = 1.1,
    group = "a", time = "post", n = 6, m = 1.5, sd = 1.2,
    group = "b", time = "pre",  n = 9, m = 1.1, sd = 1.3,
    group = "b", time = "post", n = 9, m = 2.0, sd = 1.4
  )

  expect_s3_class(design, "anovapowersim_cell_design")
  expect_named(design, c("group", "time", "n", "m", "sd"))
  expect_equal(design$n, c(6L, 6L, 9L, 9L))
  expect_error(
    cell_design(group = "a", n = 4, m = 0, sd = 1,
                group = "a", n = 4, m = 0, sd = 1),
    "only once"
  )
  expect_error(
    cell_design(group = "a", time = "pre", n = 4, m = 0, sd = 1,
                group = "b", time = "post", n = 4, m = 0, sd = 1),
    "missing"
  )
  expect_error(
    cell_design(group = "a", n = 4, m = 0, sd = 0),
    "positive finite"
  )
})


test_that("unbalanced_covariance stores correlations without SD inputs", {
  covariance <- unbalanced_covariance(
    default_correlation = 0.4,
    correlations = c("pre:post" = 0.7)
  )

  expect_s3_class(
    covariance, "anovapowersim_unbalanced_covariance_spec"
  )
  expect_equal(covariance$default_correlation, 0.4)
  expect_equal(covariance$correlations[[1L]], 0.7)
  expect_false("default_sd" %in% names(covariance))
  expect_false("standard_deviations" %in% names(covariance))
  expect_error(unbalanced_covariance(default_correlation = 1), "in \\(-1, 1\\)")
  expect_error(
    unbalanced_covariance(
      correlations = c("pre:post" = 0.6, "post:pre" = 0.5)
    ),
    "defined only once"
  )
})


test_that("power_unbalanced handles fixed pure-between designs", {
  design <- cell_design(
    group = "control", n = 7,  m = 0.0, sd = 1.0,
    group = "treated", n = 11, m = 0.8, sd = 1.4
  )
  result <- quiet_power_unbalanced(
    design = design,
    term = "group",
    n_sims = 6,
    progress = FALSE,
    seed = 10
  )

  expect_s3_class(result, "anovapowersim_unbalanced_power")
  expect_equal(result$total_n, 18L)
  expect_equal(result$min_cell_n, 7L)
  expect_equal(result$max_cell_n, 11L)
  expect_equal(result$power, result$results$power_sim[[1L]])
  expect_equal(result$achieved_power, result$power)
  expect_true(is.finite(result$partial_eta_squared))
  expect_true(result$partial_eta_squared >= 0 &&
                result$partial_eta_squared <= 1)
  expect_false(any(c("power_calc", "calculated_power", "ncp") %in%
                     names(result$results)))
  expect_false(any(c("power_calc", "calculated_power", "ncp") %in%
                     names(result)))
})


test_that("power_unbalanced combines cell SDs with within correlations", {
  design <- cell_design(
    group = "control", time = "pre",  n = 7, m = 10, sd = 2.0,
    group = "control", time = "post", n = 7, m = 11, sd = 2.5,
    group = "treated", time = "pre",  n = 9, m = 10, sd = 3.0,
    group = "treated", time = "post", n = 9, m = 13, sd = 4.0,
    within = "time"
  )
  covariance <- unbalanced_covariance(
    correlations = c("pre:post" = 0.6)
  )
  spec <- anovapowersim:::prepare_unbalanced_means_design(design)
  spec$within_correlation <-
    anovapowersim:::resolve_unbalanced_correlation(covariance, spec)
  simulated <- anovapowersim:::simulate_unbalanced_means_data(
    spec, empirical = TRUE
  )
  wide <- tidyr::pivot_wider(
    simulated,
    id_cols = c("id", "group"),
    names_from = "time",
    values_from = "value"
  )
  control <- wide[wide$group == "control", ]
  treated <- wide[wide$group == "treated", ]

  expect_equal(
    unname(stats::cov(as.data.frame(control[c("pre", "post")]))),
    matrix(c(4, 3, 3, 6.25), 2),
    tolerance = 1e-10
  )
  expect_equal(
    unname(stats::cov(as.data.frame(treated[c("pre", "post")]))),
    matrix(c(9, 7.2, 7.2, 16), 2),
    tolerance = 1e-10
  )

  result <- quiet_power_unbalanced(
    design = design,
    term = "group:time",
    covariance = covariance,
    n_sims = 4,
    progress = FALSE,
    seed = 11
  )
  expect_equal(result$total_n, 16L)
  expect_equal(result$correlation["pre", "post"], 0.6)
  expect_true(result$custom_covariance)
  expect_true(all(is.finite(unlist(result$results[c(
    "partial_eta_squared", "mean_pes_sim", "median_pes_sim",
    "pes_sim_lower", "pes_sim_upper", "power_sim"
  )]))))
})


test_that("power_unbalanced supports pure within designs and fixed seeds", {
  design <- cell_design(
    time = "pre",  n = 7, m = 0.0, sd = 1.0,
    time = "post", n = 7, m = 0.7, sd = 1.2,
    within = "time"
  )
  args <- list(
    design = design,
    term = "time",
    n_sims = 5,
    progress = FALSE,
    seed = 99
  )
  first <- do.call(quiet_power_unbalanced, args)
  second <- do.call(quiet_power_unbalanced, args)

  expect_equal(first$results, second$results)
  expect_equal(first$total_n, 7L)
})


test_that("power_unbalanced supports SS type and parallel controls", {
  design <- cell_design(
    group = "a", n = 5, m = 0, sd = 1,
    group = "b", n = 7, m = 1, sd = 1
  )
  result <- quiet_power_unbalanced(
    design,
    term = "group",
    ss_type = "II",
    n_sims = 2,
    parallel = TRUE,
    cores = 1,
    progress = FALSE,
    seed = 17
  )

  expect_identical(result$ss_type, "II")
  expect_equal(result$valid_sims + result$failed_sims, 2L)
  expect_true(is.finite(result$power))
})


test_that("power_unbalanced validates design, covariance, and controls", {
  expect_error(
    cell_design(
      group = "a", time = "pre",  n = 5, m = 0, sd = 1,
      group = "a", time = "post", n = 6, m = 1, sd = 1,
      group = "b", time = "pre",  n = 7, m = 0, sd = 1,
      group = "b", time = "post", n = 7, m = 2, sd = 1,
      within = "time"
    ),
    "identical across"
  )
  between <- cell_design(
    group = "a", n = 5, m = 0, sd = 1,
    group = "b", n = 7, m = 1, sd = 1
  )
  expect_error(
    power_unbalanced(between, term = "group",
                     covariance = unbalanced_covariance(),
                     n_sims = 1, progress = FALSE),
    "can only be supplied"
  )
  within_design <- cell_design(
    time = "pre", n = 5, m = 0, sd = 1,
    time = "post", n = 5, m = 1, sd = 1,
    within = "time"
  )
  expect_error(
    power_unbalanced(within_design, term = "time",
                     covariance = within_covariance(),
                     n_sims = 1, progress = FALSE),
    "unbalanced_covariance"
  )
  expect_error(
    power_unbalanced(within_design, term = "time",
                     covariance = diag(2), n_sims = 1, progress = FALSE),
    "raw matrices are not accepted"
  )
  expect_error(
    power_unbalanced(between, term = "group", n_sims = 1.5,
                     progress = FALSE),
    "positive integer"
  )
})


test_that("unbalanced power has dedicated print and summary methods", {
  design <- cell_design(
    group = "a", n = 5, m = 0, sd = 1,
    group = "b", n = 7, m = 1, sd = 1
  )
  result <- quiet_power_unbalanced(
    design, term = "group", n_sims = 2, progress = FALSE, seed = 1
  )

  expect_output(print(result), "simulated power")
  expect_output(print(result), "reference pes")
  expect_output(summary(result), "unbalanced-design power summary")
  summary_value <- suppressMessages(capture.output(
    summary_object <- summary(result)
  ))
  expect_named(summary_object, c("header", "design", "results"))
})


test_that("new unbalanced-design functions are exported", {
  exports <- getNamespaceExports("anovapowersim")
  expect_true(all(c(
    "cell_design", "unbalanced_covariance", "power_unbalanced"
  ) %in% exports))
})


test_that("unbalanced_term_epsilon uses the worst-case epsilon across cells", {
  design <- cell_design(
    group = "A", time = "t1", n = 10, m = 0, sd = 1,
    group = "A", time = "t2", n = 10, m = 0, sd = 1,
    group = "A", time = "t3", n = 10, m = 0, sd = 1,
    group = "B", time = "t1", n = 10, m = 0, sd = 1,
    group = "B", time = "t2", n = 10, m = 0, sd = 1,
    group = "B", time = "t3", n = 10, m = 0, sd = 4,
    within = "time"
  )
  spec <- anovapowersim:::prepare_unbalanced_means_design(design)
  spec$within_correlation <-
    anovapowersim:::resolve_unbalanced_correlation(NULL, spec)
  eps <- anovapowersim:::unbalanced_term_epsilon(spec, "time")
  eps_group_b <- anovapowersim:::covariance_term_epsilon(
    covariance = outer(spec$cell_sds[2, ], spec$cell_sds[2, ]) *
      spec$within_correlation,
    spec = spec,
    term = "time"
  )

  expect_true(eps < 1)
  expect_equal(eps, eps_group_b)

  between_design <- cell_design(
    group = "A", n = 5, m = 0, sd = 1,
    group = "B", n = 5, m = 1, sd = 1
  )
  between_spec <- anovapowersim:::prepare_unbalanced_means_design(
    between_design
  )
  between_spec$within_correlation <-
    anovapowersim:::resolve_unbalanced_correlation(NULL, between_spec)
  expect_equal(
    anovapowersim:::unbalanced_term_epsilon(between_spec, "group"), 1
  )
})


test_that("power_unbalanced GG-corrects power_sim under severe non-sphericity", {
  design <- cell_design(
    group = "A", time = "t1", n = 18, m = 0.0, sd = 1,
    group = "A", time = "t2", n = 18, m = 0.5, sd = 1,
    group = "A", time = "t3", n = 18, m = 1.0, sd = 1,
    group = "A", time = "t4", n = 18, m = 1.5, sd = 4,
    group = "B", time = "t1", n = 25, m = 0.0, sd = 1,
    group = "B", time = "t2", n = 25, m = 0.5, sd = 1,
    group = "B", time = "t3", n = 25, m = 1.0, sd = 1,
    group = "B", time = "t4", n = 25, m = 1.5, sd = 4,
    within = "time"
  )
  covariance <- unbalanced_covariance(default_correlation = 0)

  corrected <- quiet_power_unbalanced(
    design = design,
    term = "time",
    covariance = covariance,
    n_sims = 800,
    ss_type = "III",
    progress = FALSE,
    seed = 123
  )
  uncorrected <- quiet_power_unbalanced(
    design = design,
    term = "time",
    covariance = covariance,
    n_sims = 800,
    ss_type = "I",
    progress = FALSE,
    seed = 123
  )

  expect_equal(corrected$epsilon, 0.4451295, tolerance = 1e-6)
  expect_equal(corrected$results$epsilon, corrected$epsilon)
  # ss_type = "I" cannot supply a Greenhouse-Geisser p-value, so its power
  # stays on the uncorrected, Type-I-inflated univariate test and is
  # noticeably higher than the corrected estimate.
  expect_true(corrected$power < uncorrected$power - 0.05)
})


test_that("power_unbalanced warns when ss_type = 'I' meets a non-spherical design", {
  design <- cell_design(
    group = "A", time = "t1", n = 10, m = 0, sd = 1,
    group = "A", time = "t2", n = 10, m = 0, sd = 1,
    group = "A", time = "t3", n = 10, m = 0, sd = 4,
    group = "B", time = "t1", n = 12, m = 0, sd = 1,
    group = "B", time = "t2", n = 12, m = 0, sd = 1,
    group = "B", time = "t3", n = 12, m = 0, sd = 4,
    within = "time"
  )
  expect_warning(
    power_unbalanced(
      design = design,
      term = "time",
      n_sims = 1,
      ss_type = "I",
      progress = FALSE,
      seed = 1
    ),
    "ss_type = \"I\""
  )
})


test_that("cell_design lists exact missing cells instead of a bare count", {
  expect_error(
    cell_design(
      group = "a", time = "pre",  n = 10, m = 0, sd = 1,
      group = "a", time = "post", n = 10, m = 1, sd = 1,
      group = "b", time = "pre",  n = 10, m = 0, sd = 1,
      within = "time"
    ),
    'group = "b", time = "post"',
    fixed = TRUE
  )
})


test_that("cell_design fills missing cells when all three defaults are supplied", {
  filled <- cell_design(
    group = "a", time = "pre",  n = 10, m = 0, sd = 1,
    group = "a", time = "post", n = 10, m = 1, sd = 1,
    group = "b", time = "pre",  n = 10, m = 0, sd = 1,
    within = "time",
    default_n = 10, default_m = 99, default_sd = 5
  )
  filled_row <- filled[filled$group == "b" & filled$time == "post", ]

  expect_equal(nrow(filled), 4L)
  expect_equal(filled_row$n, 10L)
  expect_equal(filled_row$m, 99)
  expect_equal(filled_row$sd, 5)
})


test_that("cell_design requires all three defaults together", {
  common <- list(
    group = "a", time = "pre",  n = 10, m = 0, sd = 1,
    group = "a", time = "post", n = 10, m = 1, sd = 1,
    group = "b", time = "pre",  n = 10, m = 0, sd = 1,
    within = "time"
  )
  expect_error(
    do.call(cell_design, c(common, list(default_n = 10))),
    "Supply all of"
  )
  expect_error(
    do.call(cell_design, c(common, list(default_n = 10, default_m = 0))),
    "Supply all of"
  )
})


test_that("default-filled cells still participate in the n-consistency check", {
  expect_error(
    cell_design(
      group = "A", time = "pre",  n = 10, m = 0, sd = 1,
      group = "A", time = "post", n = 10, m = 1, sd = 1,
      group = "B", time = "pre",  n = 15, m = 0, sd = 1,
      within = "time",
      default_n = 99, default_m = 0, default_sd = 1
    ),
    "identical across"
  )
})


test_that("cell_design errors early when a factor has fewer than 2 levels", {
  expect_error(
    cell_design(
      group = "A", time = "pre",  n = 10, m = 0, sd = 1,
      group = "A", time = "post", n = 10, m = 1, sd = 1,
      within = "time"
    ),
    "group.*must have at least 2 levels"
  )
})


test_that("cell_design stores within as a retrievable attribute", {
  within_design <- cell_design(
    group = "A", time = "pre",  n = 10, m = 0, sd = 1,
    group = "A", time = "post", n = 10, m = 1, sd = 1,
    group = "B", time = "pre",  n = 10, m = 0, sd = 1,
    group = "B", time = "post", n = 10, m = 1, sd = 1,
    within = "time"
  )
  between_design <- cell_design(
    group = "A", n = 10, m = 0, sd = 1,
    group = "B", n = 10, m = 1, sd = 1
  )

  expect_identical(attr(within_design, "within"), "time")
  expect_identical(attr(between_design, "within"), character(0))
})


test_that("cell_design supports multiple within-subject factors end to end", {
  design <- cell_design(
    group = "A", time = "pre",  cond = "control", n = 10, m = 0.0, sd = 1,
    group = "A", time = "pre",  cond = "treat",   n = 10, m = 0.5, sd = 1,
    group = "A", time = "post", cond = "control", n = 10, m = 0.2, sd = 1,
    group = "A", time = "post", cond = "treat",   n = 10, m = 1.0, sd = 1,
    group = "B", time = "pre",  cond = "control", n = 15, m = 0.0, sd = 1,
    group = "B", time = "pre",  cond = "treat",   n = 15, m = 0.6, sd = 1,
    group = "B", time = "post", cond = "control", n = 15, m = 0.3, sd = 1,
    group = "B", time = "post", cond = "treat",   n = 15, m = 1.4, sd = 1,
    within = c("time", "cond")
  )

  expect_identical(attr(design, "within"), c("time", "cond"))

  result <- quiet_power_unbalanced(
    design = design,
    term = "group:time:cond",
    covariance = unbalanced_covariance(
      correlations = c("pre_control:post_control" = 0.6)
    ),
    n_sims = 3,
    progress = FALSE,
    seed = 1
  )
  expect_s3_class(result, "anovapowersim_unbalanced_power")
  expect_true(is.finite(result$power))
})
