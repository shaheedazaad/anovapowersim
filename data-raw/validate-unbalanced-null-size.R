# Long-running validation of nominal alpha for common-variance unbalanced ANOVA
# designs. Run from the package root with:
# Rscript data-raw/validate-unbalanced-null-size.R

devtools::load_all(quiet = TRUE)

n_sims <- 10000L
alpha <- 0.05

run_null_case <- function(label, design, term, covariance = NULL,
                          ss_type = "III", sim_correction = "auto",
                          gg_corrected = FALSE, expect_liberal = FALSE, seed) {
  result <- suppressMessages(suppressWarnings(power_unbalanced(
    design = design,
    term = term,
    covariance = covariance,
    n_sims = n_sims,
    alpha = alpha,
    ss_type = ss_type,
    sim_correction = sim_correction,
    progress = FALSE,
    seed = seed
  )))
  interval <- anovapowersim:::binomial_wilson_interval(
    successes = round(result$power * result$valid_sims),
    trials = result$valid_sims
  )
  passes <- if (expect_liberal) {
    result$power > 0.06
  } else if (gg_corrected) {
    result$power <= 0.06
  } else {
    abs(result$power - alpha) <= 0.01
  }
  data.frame(
    scenario = label,
    ss_type = ss_type,
    sim_correction = sim_correction,
    gg_corrected = gg_corrected,
    rejection_rate = result$power,
    lower_95 = interval[["lower"]],
    upper_95 = interval[["upper"]],
    valid_sims = result$valid_sims,
    failed_sims = result$failed_sims,
    passes = passes
  )
}

two_group_10_40 <- cell_design(
  group = "small", n = 10, m = 0,
  group = "large", n = 40, m = 0
)
two_group_40_10 <- cell_design(
  group = "small", n = 40, m = 0,
  group = "large", n = 10, m = 0
)
factorial <- cell_design(
  a = "a1", b = "b1", n = 8,  m = 0,
  a = "a1", b = "b2", n = 15, m = 0,
  a = "a2", b = "b1", n = 30, m = 0,
  a = "a2", b = "b2", n = 45, m = 0
)
mixed <- cell_design(
  group = "small", time = "t1", n = 12, m = 0,
  group = "small", time = "t2", n = 12, m = 0,
  group = "small", time = "t3", n = 12, m = 0,
  group = "large", time = "t1", n = 35, m = 0,
  group = "large", time = "t2", n = 35, m = 0,
  group = "large", time = "t3", n = 35, m = 0,
  within = "time"
)

spherical <- unbalanced_covariance(sd = 2, default_correlation = 0.5)
nonspherical <- unbalanced_covariance(
  sd = 2,
  default_correlation = 0,
  correlations = c("t1:t2" = 0.8)
)

results <- rbind(
  run_null_case(
    "two groups n=10/40", two_group_10_40, "group", seed = 101
  ),
  run_null_case(
    "two groups n=40/10", two_group_40_10, "group", seed = 102
  ),
  run_null_case(
    "2x2 unequal cells", factorial, "a:b", ss_type = "I", seed = 201
  ),
  run_null_case(
    "2x2 unequal cells", factorial, "a:b", ss_type = "II", seed = 202
  ),
  run_null_case(
    "2x2 unequal cells", factorial, "a:b", ss_type = "III", seed = 203
  ),
  run_null_case(
    "mixed spherical", mixed, "group:time", covariance = spherical,
    seed = 301
  ),
  run_null_case(
    "mixed nonspherical GG", mixed, "group:time",
    covariance = nonspherical, sim_correction = "GG",
    gg_corrected = TRUE, seed = 302
  ),
  run_null_case(
    "mixed nonspherical uncorrected", mixed, "group:time",
    covariance = nonspherical, sim_correction = "none",
    expect_liberal = TRUE, seed = 303
  )
)

print(results, row.names = FALSE)
stopifnot(all(results$passes), all(results$failed_sims == 0L))
