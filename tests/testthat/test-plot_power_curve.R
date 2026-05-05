quiet_power_curve <- function(...) suppressWarnings(power_curve(...))

test_that("plot_power_curve returns a ggplot with the expected basic layers", {
  pc <- quiet_power_curve(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    n_range = c(6L, 12L, 24L),
    n_sims = 10,
    seed = 5
  )
  g <- plot_power_curve(pc)
  expect_s3_class(g, "ggplot")
  # Layer classes: at least one GeomLine and one GeomPoint.
  geoms <- vapply(g$layers, function(layer) class(layer$geom)[1L], character(1))
  expect_true("GeomLine"   %in% geoms)
  expect_true("GeomPoint"  %in% geoms)
  expect_false("GeomRibbon" %in% geoms)
})

test_that("plot_power_curve honours show_target / show_n_needed flags", {
  pc <- quiet_power_curve(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    n_range = c(6L, 12L),
    n_sims = 10,
    seed = 6
  )

  g <- plot_power_curve(pc, show_target = FALSE, show_n_needed = FALSE)
  geoms <- vapply(g$layers, function(layer) class(layer$geom)[1L], character(1))
  expect_false("GeomRibbon" %in% geoms)
  expect_false("GeomHline"  %in% geoms)
  expect_false("GeomVline"  %in% geoms)
})

test_that("plot_power_curve can add extra power reference lines", {
  pc <- quiet_power_curve(
    between = c(group = 2),
    within = c(time = 2),
    term = "group:time",
    target_pes = 0.2,
    n_range = c(6L, 12L),
    n_sims = 10,
    seed = 8
  )

  g <- plot_power_curve(pc, power_lines = c(.80, .90))
  geoms <- vapply(g$layers, function(layer) class(layer$geom)[1L], character(1))
  expect_true("GeomHline" %in% geoms)
})

test_that("plot_power_curve errors on the wrong input", {
  expect_error(plot_power_curve(iris), "anovapowersim_curve")
})
