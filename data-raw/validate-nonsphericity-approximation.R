# Run from the repository root: Rscript data-raw/validate-nonsphericity-approximation.R
# Quantifies the approximation; Monte Carlo estimates are not regression thresholds.
devtools::load_all(quiet = TRUE)
set.seed(71021)
s <- balanced_anova_design(within = c(time = 4))
W <- t(contr.poly(4))
Sigma <- .8^abs(outer(1:4, 1:4, "-"))
Psi <- W %*% Sigma %*% t(W)
d <- nrow(W)
epsilon <- sum(diag(Psi))^2 / (d * sum(Psi^2))
n <- 12L
reps <- 100000L
# Wishart residual SSP and independent sample means are exact sufficient
# statistics for the normal model. No package fitting helper is used here.
E <- rWishart(reps, n - 1, Psi)
traces <- apply(E, 3, function(x) sum(diag(x)))
eps_hat <- traces^2 / (d * apply(E, 3, function(x) sum(x^2)))
noise <- matrix(rnorm(reps * d), reps, d) %*% chol(Psi / n)
vectors <- eigen(Psi)$vectors
scale <- sqrt((n - 1) / n * .2 / (1 - .2) * sum(diag(Psi)))
for (direction in c(1L, 3L)) {
  theta <- scale * vectors[, direction]
  Z <- sweep(noise, 2, theta, "+")
  f <- n * (n - 1) * rowSums(Z^2) / traces
  p <- pf(f, eps_hat * d, eps_hat * (n - 1) * d, lower.tail = FALSE)
  mc_power <- mean(p < .05)
  calc <- power_achieved_calc(within = c(time = 4), term = "time", target_pes = .2, n = n, epsilon = epsilon)$calculated_power
  cat("Direction", direction, "epsilon", epsilon, "calculated power", calc, "GG power", mc_power, "SE", sqrt(mc_power * (1 - mc_power) / reps), "\n")
  # Cross-check the actual fitting/extraction on a subset of full datasets.
  mu <- as.numeric(t(W) %*% theta)
  observed <- matrix(NA_real_, 50, 2)
  for (i in 1:50) {
    Y <- matrix(rnorm(n * 4), n, 4) %*% chol(Sigma)
    Y <- sweep(Y, 2, mu, "+")
    Z <- Y %*% t(W)
    centred <- sweep(Z, 2, colMeans(Z), "-")
    R <- crossprod(centred)
    F <- n * (n - 1) * sum(colMeans(Z)^2) / sum(diag(R))
    eps <- sum(diag(R))^2 / (d * sum(R^2))
    expected <- pf(F, eps * d, eps * (n - 1) * d, lower.tail = FALSE)
    dat <- data.frame(id = factor(rep(1:n, each = 4)), time = factor(rep(paste0("time", 1:4), n)), value = as.numeric(t(Y)))
    actual <- fit_design_term_stats(dat, s, "time", "III")$p_value_gg
    observed[i, ] <- c(expected, actual)
  }
  stopifnot(max(abs(observed[, 1] - observed[, 2])) < 1e-10)
  cat("Full-dataset GG oracle checks: 50; max p error", max(abs(observed[, 1] - observed[, 2])), "\n")
}
