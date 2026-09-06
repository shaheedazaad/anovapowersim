test_that("f_to_pes converts known effect sizes and handles numeric extremes", {
  expect_equal(f_to_pes(0.25), 1 / 17)
  expect_equal(f_to_pes(1), 0.5)
  expect_equal(f_to_pes(2), 0.8)
  expect_identical(f_to_pes(0), 0)
  expect_identical(f_to_pes(.Machine$double.xmax), 1)
  expect_true("f_to_pes" %in% getNamespaceExports("anovapowersim"))
})

test_that("f_to_pes requires one finite nonnegative numeric value", {
  invalid <- list(-0.25, NA_real_, NaN, Inf, -Inf, numeric(), c(0.1, 0.2),
                  "0.25", TRUE, NULL, list(0.25), 0.25 + 1i)
  for (f in invalid) {
    expect_error(f_to_pes(f),
                 "`f` must be a single finite, nonnegative numeric value.",
                 fixed = TRUE)
  }
})

test_that("f_to_pes works inline for between, within, and interaction effects", {
  designs <- list(
    list(between = c(group = 2), term = "group"),
    list(within = c(time = 3), term = "time"),
    list(between = c(group = 2), within = c(time = 3), term = "group:time")
  )
  for (design in designs) {
    converted <- do.call(power_achieved_calc,
                         c(design, list(n = 30, target_pes = f_to_pes(0.25))))
    direct <- do.call(power_achieved_calc,
                      c(design, list(n = 30, target_pes = 1 / 17)))
    expect_equal(converted$calculated_power, direct$calculated_power)
    expect_equal(converted$target_pes, 1 / 17)
  }
})
