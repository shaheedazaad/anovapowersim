test_that("calibration preserves measurement scale and target noncentrality", {
  designs <- list(
    list(design = balanced_anova_design(between = c(A = 3)), term = "A"),
    list(design = balanced_anova_design(within = c(A = 3)), term = "A"),
    list(design = balanced_anova_design(between = c(A = 2),
                                      within = c(B = 3)), term = "A:B")
  )
  for (s in designs) for (ss in c("I", "II", "III")) {
    for (g in c(FALSE, TRUE)) {
      set.seed(300602)
      reference <- suppressWarnings(design_term_means(
        s$design, s$term, target_pes = f_to_pes(.25), n = 8,
        ss_type = ss, gpower = g
      ))
      # car::Anova.lm separately rejects absolute residual SS < sqrt(eps).
      scales <- if (!length(s$design$within) && ss != "I") {
        c(1e-3, 1e6)
      } else c(1e-6, 1e-8, 1e6)
      for (scale in scales) {
        set.seed(300602)
        means <- suppressWarnings(design_term_means(
          s$design, s$term, target_pes = f_to_pes(.25), n = 8,
          sd = scale, ss_type = ss, gpower = g
        ))
        expect_equal(means / scale, reference, tolerance = 1e-6)
        exact <- simulate_design_dataset(s$design, 8, means,
                                        sd = scale, empirical = TRUE)
        fit <- anovapowersim:::fit_design_term_stats(
          exact, s$design, s$term, ss_type = ss
        )
        expected <- if (g) 8 * max(1L, s$design$n_between_cells) * .25^2 else {
          fit$den_df * .25^2
        }
        observed <- fit$pes / (1 - fit$pes) * fit$den_df
        expect_equal(observed, expected, tolerance = 1e-6)
      }
    }
  }
})

test_that("small-SD public simulation reproducer succeeds with the same power", {
  run <- function(sd) suppressWarnings(power_achieved(
    within = c(A = 2), term = "A", n = 3, target_pes = .1,
    covariance = within_covariance(sd = sd),
    n_sims = 20, seed = 300602, progress = FALSE
  ))
  reference <- run(1)$results
  for (scale in c(1e-6, 1e-8)) {
    actual <- run(scale)$results
    expect_equal(actual$power_sim, reference$power_sim)
    expect_equal(actual$power_calc, reference$power_calc)
  }
})
