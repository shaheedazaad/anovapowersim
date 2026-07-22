quiet_power_unbalanced <- function(...) suppressWarnings(power_unbalanced(...))


test_that("cell_design creates a complete n-and-means cell table", {
  design <- cell_design(
    group = "a", time = "pre",  n = 6, m = 1.0,
    group = "a", time = "post", n = 6, m = 1.5,
    group = "b", time = "pre",  n = 9, m = 1.1,
    group = "b", time = "post", n = 9, m = 2.0
  )

  expect_s3_class(design, "anovapowersim_cell_design")
  expect_named(design, c("group", "time", "n", "m"))
  expect_equal(design$n, c(6L, 6L, 9L, 9L))
  expect_error(
    cell_design(group = "a", n = 4, m = 0,
                group = "a", n = 4, m = 0),
    "only once"
  )
  expect_error(
    cell_design(group = "a", time = "pre", n = 4, m = 0,
                group = "b", time = "post", n = 4, m = 0),
    "missing"
  )
})


test_that("cell_design rejects legacy cell-level SD inputs clearly", {
  expect_error(
    cell_design(
      group = "a", n = 4, m = 0, sd = 1,
      group = "b", n = 4, m = 1, sd = 1
    ),
    "unbalanced_covariance\\(sd = \\.\\.\\.\\)"
  )
  expect_error(
    cell_design(
      group = "a", n = 4, m = 0,
      group = "b", n = 4, m = 1,
      default_sd = 1
    ),
    "Standard deviations no longer belong"
  )
})


test_that("unbalanced_covariance stores one common SD and correlations", {
  covariance <- unbalanced_covariance(
    sd = 2,
    default_correlation = 0.4,
    correlations = c("pre:post" = 0.7)
  )

  expect_s3_class(
    covariance, "anovapowersim_unbalanced_covariance_spec"
  )
  expect_equal(covariance$sd, 2)
  expect_equal(covariance$default_correlation, 0.4)
  expect_equal(covariance$correlations[[1L]], 0.7)
  expect_false(any(c("default_sd", "standard_deviations") %in%
                     names(covariance)))
  expect_warning(
    unbalanced_covariance(),
    "using one common `sd = 1`"
  )
  expect_error(unbalanced_covariance(sd = 0), "`sd`")
  expect_error(unbalanced_covariance(default_correlation = 1), "in \\(-1, 1\\)")
  expect_error(
    unbalanced_covariance(
      correlations = c("pre:post" = 0.6, "post:pre" = 0.5)
    ),
    "defined only once"
  )
})


test_that("power_unbalanced handles unequal-N pure-between designs", {
  design <- cell_design(
    group = "control", n = 7,  m = 0.0,
    group = "treated", n = 11, m = 0.8
  )
  result <- quiet_power_unbalanced(
    design = design,
    term = "group",
    covariance = unbalanced_covariance(sd = 2),
    n_sims = 6,
    progress = FALSE,
    seed = 10
  )

  expect_s3_class(result, "anovapowersim_unbalanced_power")
  expect_equal(result$total_n, 18L)
  expect_equal(result$min_cell_n, 7L)
  expect_equal(result$max_cell_n, 11L)
  expect_equal(result$sd, 2)
  expect_equal(result$results$sd, 2)
  expect_equal(result$power, result$results$power_sim[[1L]])
  expect_equal(result$achieved_power, result$power)
  expect_true(is.finite(result$partial_eta_squared))
  expect_true(result$partial_eta_squared >= 0 &&
                result$partial_eta_squared <= 1)
  expect_false(any(c("power_calc", "calculated_power", "ncp") %in%
                     names(result$results)))
})


test_that("power_unbalanced warns when covariance defaults are used", {
  design <- cell_design(
    group = "small", n = 5, m = 0,
    group = "large", n = 12, m = 0
  )
  warnings <- testthat::capture_warnings(
    power_unbalanced(
      design,
      term = "group",
      n_sims = 1,
      progress = FALSE,
      seed = 14
    )
  )

  expect_true(any(grepl(
    "No SD or `covariance` was supplied; using one common `sd = 1`",
    warnings,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Correlations do not apply to this purely between-subject design",
    warnings,
    fixed = TRUE
  )))
})


test_that("power_unbalanced warns when supplied means omit the tested effect", {
  design <- cell_design(
    group = "control", time = "pre",  n = 6, m = 0,
    group = "control", time = "post", n = 6, m = 1,
    group = "treated", time = "pre",  n = 8, m = 0,
    group = "treated", time = "post", n = 8, m = 1,
    within = "time"
  )

  expect_warning(
    result <- power_unbalanced(
      design,
      term = "group:time",
      covariance = unbalanced_covariance(
        sd = 1,
        correlations = c("pre:post" = 0.5)
      ),
      n_sims = 1,
      progress = FALSE,
      seed = 15
    ),
    "essentially no effect.*group:time.*partial eta squared"
  )
  expect_lte(result$partial_eta_squared, 1e-8)
})


test_that("unbalanced designs share one covariance across groups", {
  design <- cell_design(
    group = "control", time = "pre",  n = 7, m = 10,
    group = "control", time = "post", n = 7, m = 11,
    group = "treated", time = "pre",  n = 9, m = 10,
    group = "treated", time = "post", n = 9, m = 13,
    within = "time"
  )
  covariance <- unbalanced_covariance(
    sd = 2,
    correlations = c("pre:post" = 0.6)
  )
  spec <- anovapowersim:::prepare_unbalanced_means_design(design)
  spec$sd <- covariance$sd
  expect_silent(
    spec$within_correlation <-
      anovapowersim:::resolve_unbalanced_correlation(covariance, spec)
  )
  simulated <- anovapowersim:::simulate_unbalanced_means_data(
    spec, empirical = TRUE
  )
  wide <- tidyr::pivot_wider(
    simulated,
    id_cols = c("id", "group"),
    names_from = "time",
    values_from = "value"
  )
  expected <- matrix(c(4, 2.4, 2.4, 4), 2)

  for (group in c("control", "treated")) {
    rows <- wide[wide$group == group, ]
    expect_equal(
      unname(stats::cov(as.data.frame(rows[c("pre", "post")]))),
      expected,
      tolerance = 1e-10
    )
  }

  result <- quiet_power_unbalanced(
    design = design,
    term = "group:time",
    covariance = covariance,
    n_sims = 4,
    progress = FALSE,
    seed = 11
  )
  expect_equal(result$correlation["pre", "post"], 0.6)
  expect_equal(result$sd, 2)
  expect_true(result$custom_covariance)
})


test_that("power_unbalanced supports pure within designs and fixed seeds", {
  design <- cell_design(
    time = "pre",  n = 7, m = 0.0,
    time = "post", n = 7, m = 0.7,
    within = "time"
  )
  args <- list(
    design = design,
    term = "time",
    covariance = unbalanced_covariance(sd = 1.5),
    n_sims = 5,
    progress = FALSE,
    seed = 99
  )
  first <- do.call(quiet_power_unbalanced, args)
  second <- do.call(quiet_power_unbalanced, args)

  expect_equal(first$results, second$results)
  expect_equal(first$total_n, 7L)
  expect_equal(first$sd, 1.5)
})


test_that("power_unbalanced supports SS type and parallel controls", {
  design <- cell_design(
    group = "a", n = 5, m = 0,
    group = "b", n = 7, m = 1
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
      group = "a", time = "pre",  n = 5, m = 0,
      group = "a", time = "post", n = 6, m = 1,
      group = "b", time = "pre",  n = 7, m = 0,
      group = "b", time = "post", n = 7, m = 2,
      within = "time"
    ),
    "identical across"
  )
  between <- cell_design(
    group = "a", n = 5, m = 0,
    group = "b", n = 7, m = 1
  )
  expect_silent(
    quiet_power_unbalanced(
      between,
      term = "group",
      covariance = unbalanced_covariance(sd = 2),
      n_sims = 1,
      progress = FALSE
    )
  )
  expect_error(
    power_unbalanced(
      between,
      term = "group",
      covariance = unbalanced_covariance(
        sd = 1,
        correlations = c("pre:post" = 0.5)
      ),
      n_sims = 1,
      progress = FALSE
    ),
    "purely between-subject"
  )
  within_design <- cell_design(
    time = "pre", n = 5, m = 0,
    time = "post", n = 5, m = 1,
    within = "time"
  )
  expect_error(
    power_unbalanced(
      within_design,
      term = "time",
      covariance = within_covariance(sd = 1),
      n_sims = 1,
      progress = FALSE
    ),
    "unbalanced_covariance"
  )
  expect_error(
    power_unbalanced(
      within_design,
      term = "time",
      covariance = diag(2),
      n_sims = 1,
      progress = FALSE
    ),
    "raw matrices are not accepted"
  )
  expect_error(
    power_unbalanced(between, term = "group", n_sims = 1.5,
                     progress = FALSE),
    "positive integer"
  )
})


test_that("unbalanced power prints and summarises the common SD", {
  design <- cell_design(
    group = "a", n = 5, m = 0,
    group = "b", n = 7, m = 1
  )
  result <- quiet_power_unbalanced(
    design,
    term = "group",
    covariance = unbalanced_covariance(sd = 2),
    n_sims = 2,
    progress = FALSE,
    seed = 1
  )

  expect_output(print(result), "common SD")
  expect_output(print(result), "simulated power")
  expect_output(
    print(result),
    "sample pes is upward-biased.*diagnostics, not the population/reference"
  )
  expect_output(summary(result), "common_sd")
  expect_output(
    summary(result),
    "sample pes is upward-biased.*diagnostics, not the population/reference"
  )
  suppressMessages(capture.output(summary_object <- summary(result)))
  expect_named(summary_object, c("header", "design", "results"))
  expect_equal(unname(summary_object$header[["common_sd"]]), "2")
})


test_that("Type I SS warns about factor order for unequal-N designs", {
  design <- cell_design(
    treatment = "A", site = "north", n = 5,  m = 0.0,
    treatment = "A", site = "south", n = 7,  m = 0.2,
    treatment = "B", site = "north", n = 9,  m = 1.0,
    treatment = "B", site = "south", n = 11, m = 1.2
  )

  warnings <- testthat::capture_warnings(power_unbalanced(
    design,
    term = "treatment",
    covariance = unbalanced_covariance(sd = 1),
    n_sims = 1,
    ss_type = "I",
    progress = FALSE,
    seed = 18
  ))

  expect_true(any(grepl(
    "sequential sums of squares.*factor order.*treatment, site",
    warnings
  )))
})


test_that("new unbalanced-design functions are exported", {
  exports <- getNamespaceExports("anovapowersim")
  expect_true(all(c(
    "cell_design", "unbalanced_covariance", "power_unbalanced"
  ) %in% exports))
})


test_that("unbalanced_term_epsilon uses the shared covariance", {
  design <- cell_design(
    group = "A", time = "t1", n = 10, m = 0,
    group = "A", time = "t2", n = 10, m = 0,
    group = "A", time = "t3", n = 10, m = 0,
    group = "B", time = "t1", n = 12, m = 0,
    group = "B", time = "t2", n = 12, m = 0,
    group = "B", time = "t3", n = 12, m = 0,
    within = "time"
  )
  covariance <- unbalanced_covariance(
    sd = 2,
    default_correlation = 0,
    correlations = c("t1:t2" = 0.8)
  )
  spec <- anovapowersim:::prepare_unbalanced_means_design(design)
  spec$sd <- covariance$sd
  expect_warning(
    spec$within_correlation <-
      anovapowersim:::resolve_unbalanced_correlation(covariance, spec),
    "used only for those undefined pairs"
  )
  epsilon <- anovapowersim:::unbalanced_term_epsilon(spec, "time")
  expected <- anovapowersim:::covariance_term_epsilon(
    covariance = covariance$sd^2 * spec$within_correlation,
    spec = spec,
    term = "time"
  )

  expect_lt(epsilon, 1)
  expect_equal(epsilon, expected)

  between_design <- cell_design(
    group = "A", n = 5, m = 0,
    group = "B", n = 5, m = 1
  )
  between_spec <- anovapowersim:::prepare_unbalanced_means_design(
    between_design
  )
  between_spec$sd <- 1
  between_spec$within_correlation <- matrix(1)
  expect_equal(
    anovapowersim:::unbalanced_term_epsilon(between_spec, "group"), 1
  )
})


test_that("power_unbalanced GG-corrects equal-variance nonsphericity", {
  design <- cell_design(
    group = "A", time = "t1", n = 18, m = 0.0,
    group = "A", time = "t2", n = 18, m = 0.5,
    group = "A", time = "t3", n = 18, m = 1.0,
    group = "B", time = "t1", n = 25, m = 0.0,
    group = "B", time = "t2", n = 25, m = 0.5,
    group = "B", time = "t3", n = 25, m = 1.0,
    within = "time"
  )
  covariance <- unbalanced_covariance(
    sd = 1,
    default_correlation = 0,
    correlations = c("t1:t2" = 0.8)
  )
  corrected <- quiet_power_unbalanced(
    design = design,
    term = "time",
    covariance = covariance,
    n_sims = 10,
    ss_type = "III",
    progress = FALSE,
    seed = 123
  )

  expect_lt(corrected$epsilon, 1)
  expect_true(is.finite(corrected$power))
  expect_identical(corrected$failed_sims, 0L)
  warnings <- testthat::capture_warnings(
    power_unbalanced(
      design = design,
      term = "time",
      covariance = covariance,
      n_sims = 1,
      ss_type = "I",
      progress = FALSE,
      seed = 1
    )
  )
  expect_true(any(grepl("ss_type = \"I\"", warnings, fixed = TRUE)))
  expect_true(any(grepl("undefined pairs", warnings, fixed = TRUE)))
})


test_that("cell_design lists, fills, and validates missing cells", {
  expect_error(
    cell_design(
      group = "a", time = "pre",  n = 10, m = 0,
      group = "a", time = "post", n = 10, m = 1,
      group = "b", time = "pre",  n = 10, m = 0,
      within = "time"
    ),
    'group = "b", time = "post"',
    fixed = TRUE
  )

  filled <- suppressMessages(cell_design(
    group = "a", time = "pre",  n = 10, m = 0,
    group = "a", time = "post", n = 10, m = 1,
    group = "b", time = "pre",  n = 10, m = 0,
    within = "time",
    default_n = 10,
    default_m = 99
  ))
  filled_row <- filled[filled$group == "b" & filled$time == "post", ]
  expect_equal(nrow(filled), 4L)
  expect_equal(filled_row$n, 10L)
  expect_equal(filled_row$m, 99)
  expect_false("sd" %in% names(filled_row))

  common <- list(
    group = "a", time = "pre",  n = 10, m = 0,
    group = "a", time = "post", n = 10, m = 1,
    group = "b", time = "pre",  n = 10, m = 0,
    within = "time"
  )
  expect_error(
    do.call(cell_design, c(common, list(default_n = 10))),
    "Supply both"
  )
  expect_error(
    do.call(cell_design, c(common, list(default_m = 0))),
    "Supply both"
  )
})


test_that("cell_design reports exact cells created by defaults", {
  messages <- testthat::capture_messages(
    design <- cell_design(
      group = "control", time = "pre",  n = 10, m = 0,
      group = "contrl", time = "post",  n = 10, m = 1,
      group = "treatment", time = "pre",  n = 10, m = 0,
      group = "treatment", time = "post", n = 10, m = 1,
      within = "time",
      default_n = 10,
      default_m = 99
    )
  )

  expect_true(any(grepl("Auto-filled 2 missing cells", messages,
                        fixed = TRUE)))
  expect_true(any(grepl('group = "control", time = "post"', messages,
                        fixed = TRUE)))
  expect_true(any(grepl('group = "contrl", time = "pre"', messages,
                        fixed = TRUE)))
  expect_equal(nrow(design), 6L)
  expect_true("contrl" %in% design$group)
})


test_that("default-filled cells participate in n-consistency checks", {
  expect_error(
    suppressMessages(cell_design(
      group = "A", time = "pre",  n = 10, m = 0,
      group = "A", time = "post", n = 10, m = 1,
      group = "B", time = "pre",  n = 15, m = 0,
      within = "time",
      default_n = 99,
      default_m = 0
    )),
    "identical across"
  )
})


test_that("cell_design validates factor levels and stores within factors", {
  expect_error(
    cell_design(
      group = "A", time = "pre",  n = 10, m = 0,
      group = "A", time = "post", n = 10, m = 1,
      within = "time"
    ),
    "group.*must have at least 2 levels"
  )

  within_design <- cell_design(
    group = "A", time = "pre",  n = 10, m = 0,
    group = "A", time = "post", n = 10, m = 1,
    group = "B", time = "pre",  n = 10, m = 0,
    group = "B", time = "post", n = 10, m = 1,
    within = "time"
  )
  between_design <- cell_design(
    group = "A", n = 10, m = 0,
    group = "B", n = 10, m = 1
  )
  expect_identical(attr(within_design, "within"), "time")
  expect_identical(attr(between_design, "within"), character(0))
})


test_that("cell_design rejects unsafe constructed within-cell names", {
  expect_error(
    cell_design(
      time = "pre:baseline", n = 5, m = 0,
      time = "post", n = 5, m = 1,
      within = "time"
    ),
    "must not contain ':'.*`time` = 'pre:baseline'"
  )

  expect_error(
    cell_design(
      first = "a_b", second = "c",   n = 5, m = 0,
      first = "a_b", second = "b_c", n = 5, m = 1,
      first = "a",   second = "c",   n = 5, m = 2,
      first = "a",   second = "b_c", n = 5, m = 3,
      within = c("first", "second")
    ),
    "not unique.*'a_b_c'.*first = 'a_b'.*second = 'c'.*first = 'a'.*second = 'b_c'"
  )
})


test_that("cell_design supports multiple within-subject factors end to end", {
  design <- cell_design(
    group = "A", time = "pre",  cond = "control", n = 10, m = 0.0,
    group = "A", time = "pre",  cond = "treat",   n = 10, m = 0.5,
    group = "A", time = "post", cond = "control", n = 10, m = 0.2,
    group = "A", time = "post", cond = "treat",   n = 10, m = 1.0,
    group = "B", time = "pre",  cond = "control", n = 15, m = 0.0,
    group = "B", time = "pre",  cond = "treat",   n = 15, m = 0.6,
    group = "B", time = "post", cond = "control", n = 15, m = 0.3,
    group = "B", time = "post", cond = "treat",   n = 15, m = 1.4,
    within = c("time", "cond")
  )

  expect_identical(attr(design, "within"), c("time", "cond"))
  result <- quiet_power_unbalanced(
    design = design,
    term = "group:time:cond",
    covariance = unbalanced_covariance(
      sd = 1,
      correlations = c("pre_control:post_control" = 0.6)
    ),
    n_sims = 3,
    progress = FALSE,
    seed = 1
  )
  expect_s3_class(result, "anovapowersim_unbalanced_power")
  expect_true(is.finite(result$power))
})


test_that("common-SD simulations are scale invariant", {
  design_one <- cell_design(
    group = "A", n = 10, m = 0,
    group = "B", n = 25, m = 0.5
  )
  design_two <- cell_design(
    group = "A", n = 10, m = 0,
    group = "B", n = 25, m = 1
  )
  first <- quiet_power_unbalanced(
    design_one,
    term = "group",
    covariance = unbalanced_covariance(sd = 1),
    n_sims = 20,
    progress = FALSE,
    seed = 144
  )
  second <- quiet_power_unbalanced(
    design_two,
    term = "group",
    covariance = unbalanced_covariance(sd = 2),
    n_sims = 20,
    progress = FALSE,
    seed = 144
  )

  expect_equal(first$power, second$power)
  expect_equal(first$partial_eta_squared, second$partial_eta_squared,
               tolerance = 1e-12)
  expect_equal(first$results$mean_pes_sim, second$results$mean_pes_sim,
               tolerance = 1e-12)
})


test_that("unequal N with a common variance controls null rejection", {
  design <- cell_design(
    group = "small", n = 10, m = 0,
    group = "large", n = 40, m = 0
  )
  result <- quiet_power_unbalanced(
    design,
    term = "group",
    covariance = unbalanced_covariance(sd = 2),
    n_sims = 2000,
    progress = FALSE,
    seed = 20260722
  )

  expect_lte(abs(result$power - 0.05), 0.02)
})
