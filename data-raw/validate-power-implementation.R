# Run from the repository root: Rscript data-raw/validate-power-implementation.R
# Independent statistical implementation audit, September 2026.
devtools::load_all(quiet = TRUE)
set.seed(64291)
# Orthonormal factorial contrasts, independently constructed from polynomial
# scores rather than the package's Helmert covariance and sum-model helpers.
contrast <- function(counts, active) {
  if (!length(counts)) {
    return(matrix(1, 1, 1))
  }
  Reduce(kronecker, lapply(names(counts), function(nm) {
    k <- counts[[nm]]
    if (!(nm %in% active)) {
      return(matrix(1 / sqrt(k), 1, k))
    }
    t(contr.poly(k))
  }))
}
term_names <- function(factors) {
  if (!length(factors)) {
    return(character())
  }
  unlist(lapply(seq_along(factors), function(k) apply(combn(factors, k), 2, paste, collapse = ":")))
}
rss <- function(X, Y) {
  if (!ncol(X)) {
    return(sum(Y^2))
  }
  sum(qr.resid(qr(X), Y)^2)
}
cases <- list(
  list(c(a = 2, b = 3, c = 2), NULL), list(NULL, c(t = 2, u = 3, v = 2)),
  list(c(a = 2, b = 3), c(t = 3, u = 2)), list(c(a = 3), c(t = 3)),
  list(c(a = 3), NULL)
)
checks <- list()
for (case_i in seq_along(cases)) {
  case <- cases[[case_i]]
  for (unequal in c(FALSE, TRUE)) {
    for (ar in c(FALSE, TRUE)) {
      s <- balanced_anova_design(case[[1]], case[[2]])
      b <- s$n_between_cells
      w <- s$n_within_cells
      ns <- if (unequal) (seq_len(b) * 3 + w + 2) else rep(w + 8, b)
      s$cell_n <- ns
      s$cell_means <- matrix(rnorm(b * w, sd = .3), b, w)
      s$sd <- 1
      sigma <- if (ar) .7^abs(outer(seq_len(w), seq_len(w), "-")) else .4 + diag(.6, w)
      s$within_correlation <- sigma
      dat <- simulate_unbalanced_means_data(s)
      Y <- matrix(dat$value, ncol = w, byrow = TRUE)
      group <- rep(seq_len(b), ns)
      bt <- c("", term_names(s$between))
      bmat <- lapply(bt, function(term) t(contrast(case[[1]], strsplit(term, ":", fixed = TRUE)[[1]])))
      Xcell <- do.call(cbind, bmat)
      X <- Xcell[group, , drop = FALSE]
      assignment <- rep(seq_along(bt), vapply(bmat, ncol, integer(1)))
      for (term in term_names(s$factor_names)) {
        active <- strsplit(term, ":", fixed = TRUE)[[1]]
        bterm <- paste(intersect(s$between, active), collapse = ":")
        index <- match(bterm, bt)
        W <- contrast(case[[2]], active)
        Z <- Y %*% t(W)
        residual <- qr.resid(qr(X), Z)
        error_ss <- sum(residual^2)
        df1 <- ncol(bmat[[index]]) * nrow(W)
        df2 <- (sum(ns) - b) * nrow(W)
        S <- crossprod(residual)
        eps <- sum(diag(S))^2 / (nrow(W) * sum(S^2))
        popS <- W %*% sigma %*% t(W)
        pop_eps <- sum(diag(popS))^2 / (nrow(W) * sum(popS^2))
        eps_code <- covariance_term_epsilon(sigma, s, term)
        for (type in c("I", "II", "III")) {
          included <- switch(type,
            I = seq_len(index),
            II = which(vapply(bt, function(other) {
              other_factors <- strsplit(other, ":", fixed = TRUE)[[1]]
              target_factors <- strsplit(bterm, ":", fixed = TRUE)[[1]]
              !(all(target_factors %in% other_factors) && length(other_factors) > length(target_factors))
            }, logical(1))),
            III = seq_along(bt)
          )
          full_cols <- assignment %in% included
          reduced_cols <- full_cols & assignment != index
          numerator <- rss(X[, reduced_cols, drop = FALSE], Z) - rss(X[, full_cols, drop = FALSE], Z)
          f <- (numerator / df1) / (error_ss / df2)
          p <- pf(f, df1, df2, lower.tail = FALSE)
          pes <- numerator / (numerator + error_ss)
          expected_gg <- if (nrow(W) > 1) pf(f, eps * df1, eps * df2, lower.tail = FALSE) else NA_real_
          observed <- suppressMessages(fit_design_term_stats(dat, s, term, type))
          checks[[length(checks) + 1]] <- data.frame(
            case = case_i, unequal, ar, term, type,
            p_error = abs(p - observed$p_value), pes_error = abs(pes - observed$pes),
            df_error = abs(df1 - observed$num_df) + abs(df2 - observed$den_df),
            epsilon_error = abs(pop_eps - eps_code),
            gg_error = if (type != "I" && nrow(W) > 1) abs(expected_gg - observed$p_value_gg) else NA_real_
          )
        }
      }
    }
  }
}
checks <- do.call(rbind, checks)
print(aggregate(checks[c("p_error", "pes_error", "df_error", "epsilon_error", "gg_error")], by = checks["type"], FUN = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)))
cat("Total test-statistic cases:", nrow(checks), "\n")
stopifnot(max(checks$p_error) < 1e-8, max(checks$pes_error) < 1e-8, max(checks$df_error) == 0, max(checks$epsilon_error) < 1e-8, max(checks$gg_error, na.rm = TRUE) < 1e-8)
# Check both NCP conventions from population means and covariance, not F tables.
calibration <- list()
for (case_i in seq_along(cases)) {
  case <- cases[[case_i]]
  s <- balanced_anova_design(case[[1]], case[[2]])
  n <- s$n_within_cells + 7
  for (term in term_names(s$factor_names)) {
    for (gpower in c(FALSE, TRUE)) {
      for (type in c("I", "II", "III")) {
        active <- strsplit(term, ":", fixed = TRUE)[[1]]
        B <- contrast(case[[1]], active)
        W <- contrast(case[[2]], active)
        mu <- suppressMessages(calibrate_design_means(s, term, .17, n, 1.4, -.02, gpower = gpower, ss_type = type))
        sigma <- 1.4^2 * (-.02 + diag(1.02, s$n_within_cells))
        variance <- sum(diag(W %*% sigma %*% t(W))) / nrow(W)
        lambda <- n * sum((B %*% mu %*% t(W))^2) / variance
        expected <- analytic_power_row(n, s, term, .17, .05, gpower)$ncp
        calibration[[length(calibration) + 1]] <- data.frame(case = case_i, term, type, gpower, error = abs(lambda - expected))
      }
    }
  }
}
calibration <- do.call(rbind, calibration)
cat("Calibrated NCP checks:", nrow(calibration), "maximum error:", max(calibration$error), "\n")
stopifnot(max(calibration$error) < 1e-8)
# Literal population means versus the balanced empirical-reference convention.
s <- balanced_anova_design(c(group = 2))
for (n in c(4, 40)) {
  mu <- design_term_means(s, "group", .2, n)
  pop_f2 <- mean((mu - mean(mu))^2)
  lambda <- n * sum((mu - mean(mu))^2)
  df <- 2 * n - 2
  expected_pes <- integrate(function(f) (f / (f + df)) * stats::df(f, 1, df, ncp = lambda), 0, Inf)$value
  cat("Reference convention: n", n, "mean difference", diff(as.numeric(mu)), "population eta squared", pop_f2 / (1 + pop_f2), "mean sample PES", expected_pes, "\n")
}
