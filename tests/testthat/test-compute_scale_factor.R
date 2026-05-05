test_that("compute_scale_factor matches the closed-form derivation", {
  k <- compute_scale_factor(old_pes = 0.10, new_pes = 0.05)
  expected <- sqrt((0.05 / 0.95) / (0.10 / 0.90))
  expect_equal(k, expected)
})

test_that("compute_scale_factor returns 1 when old and new match", {
  expect_equal(compute_scale_factor(0.25, 0.25), 1)
})

test_that("compute_scale_factor accepts numeric-looking formatted values", {
  expect_equal(
    compute_scale_factor(".310", ".200"),
    compute_scale_factor(0.310, 0.200)
  )
  expect_equal(
    compute_scale_factor(" 0.310 ", " 0.200 "),
    compute_scale_factor(0.310, 0.200)
  )
})

test_that("compute_scale_factor scales monotonically", {
  k_shrink    <- compute_scale_factor(0.20, 0.10)
  k_identity  <- compute_scale_factor(0.20, 0.20)
  k_amplify   <- compute_scale_factor(0.20, 0.40)

  expect_lt(k_shrink, k_identity)
  expect_lt(k_identity, k_amplify)
})

test_that("compute_scale_factor validates its inputs", {
  expect_error(compute_scale_factor(0,   0.1), "old_pes")
  expect_error(compute_scale_factor(1.1, 0.1), "old_pes")
  expect_error(compute_scale_factor(0.1, -0.1), "new_pes")
  expect_error(compute_scale_factor(0.1, 1),    "new_pes")
  expect_error(compute_scale_factor("foo", 0.1), "old_pes")
})
