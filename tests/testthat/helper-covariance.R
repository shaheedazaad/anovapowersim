test_covariance_spec_from_matrix <- function(sigma) {
  stopifnot(
    is.matrix(sigma),
    is.numeric(sigma),
    !is.null(rownames(sigma)),
    identical(rownames(sigma), colnames(sigma)),
    length(unique(diag(sigma))) == 1L
  )
  pairs <- which(upper.tri(sigma), arr.ind = TRUE)
  correlations <- sigma[pairs] / diag(sigma)[[1L]]
  names(correlations) <- paste(
    rownames(sigma)[pairs[, "row"]],
    colnames(sigma)[pairs[, "col"]],
    sep = ":"
  )
  within_covariance(
    sd = sqrt(diag(sigma)[[1L]]),
    default_correlation = 0,
    correlations = correlations
  )
}
