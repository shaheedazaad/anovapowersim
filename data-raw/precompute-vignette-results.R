if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(anovapowersim)
}

adaptive <- power_n(
  between = c(cond = 2),
  within = c(stim = 4),
  term = "cond:stim",
  target_pes = 0.14,
  alpha = 0.05,
  power = 0.90,
  n_sims = 1000,
  seed = 123,
  progress = FALSE
)

curve <- power_curve(
  between = c(cond = 2),
  within = c(stim = 2),
  term = "cond:stim",
  target_pes = 0.14,
  n_range = c(16, 20, 23, 28),
  n_sims = 1000,
  seed = 123,
  progress = FALSE
)

gpower_adaptive <- power_n(
  between = c(cond = 2),
  within = c(stim = 4),
  term = "cond:stim",
  target_pes = 0.14,
  alpha = 0.05,
  power = 0.90,
  n_sims = 1000,
  seed = 123,
  gpower = TRUE,
  progress = FALSE
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    adaptive = adaptive,
    curve = curve,
    gpower_adaptive = gpower_adaptive
  ),
  "inst/extdata/anovapowersim-vignette-results.rds"
)
