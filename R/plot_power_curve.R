#' Plot a simulation-based power curve
#'
#' Renders an `anovapowersim_curve` as a `ggplot2` line + ribbon with a
#' horizontal reference at requested power values and, when auto-search was
#' used, a vertical marker at the estimated required total sample size.
#'
#' @param x An `anovapowersim_curve` object from [power_curve()].
#' @param show_target Logical; draw the horizontal target power line
#'   (default `TRUE`).
#' @param power_lines Optional numeric vector of additional power reference
#'   lines, e.g. `c(.80, .90)`.
#' @param show_n_needed Logical; draw the vertical line at `n_needed`
#'   (default `TRUE`).
#' @param ... Unused, for S3 consistency.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' pc <- power_curve(
#'   between = c(group = 2),
#'   within = c(time = 2),
#'   term = "group:time",
#'   target_pes = 0.2,
#'   n_range = c(10, 20, 30),
#'   n_sims = 200
#' )
#' plot_power_curve(pc)
#' }
#'
#' @seealso [power_curve()]
#' @export
plot_power_curve <- function(x,
                             show_target = TRUE,
                             power_lines = NULL,
                             show_n_needed = TRUE,
                             ...) {
  if (!inherits(x, "anovapowersim_curve")) {
    stop("`x` must be an `anovapowersim_curve` object.",
         call. = FALSE)
  }

  df <- x$results
  n_between_cells <- nrow(x$design$between_cells)
  refs <- numeric(0)
  if (isTRUE(show_target) && is.finite(x$power)) refs <- c(refs, x$power)
  if (!is.null(power_lines)) {
    if (!is.numeric(power_lines) || any(!is.finite(power_lines)) ||
        any(power_lines <= 0 | power_lines >= 1)) {
      stop("`power_lines` must be a numeric vector with values in (0, 1).",
           call. = FALSE)
    }
    refs <- c(refs, power_lines)
  }
  refs <- sort(unique(refs))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$total_n, y = .data$power_sim))

  p <- p +
    ggplot2::geom_line(colour = "#4477AA", linewidth = 0.8) +
    ggplot2::geom_point(colour = "#4477AA", size = 2.2)

  if (length(refs)) {
    p <- p +
      ggplot2::geom_hline(
        data = tibble::tibble(power_ref = refs),
        ggplot2::aes(yintercept = .data$power_ref),
        inherit.aes = FALSE,
        linetype = "dashed",
        colour = "grey35",
        linewidth = 0.4
      ) +
      ggplot2::geom_text(
        data = tibble::tibble(power_ref = refs),
        ggplot2::aes(
          x = max(df$total_n, na.rm = TRUE),
          y = .data$power_ref,
          label = paste0(round(.data$power_ref * 100), "%")
        ),
        inherit.aes = FALSE,
        hjust = 1,
        vjust = -0.4,
        size = 3,
        colour = "grey25"
      )
  }

  if (isTRUE(show_n_needed) && !is.null(x$n_needed) && !is.na(x$n_needed)) {
    p <- p + ggplot2::geom_vline(
      xintercept = x$total_n_needed,
      linetype   = "dotted",
      colour     = "grey25",
      linewidth = 0.5
    )
  }

  subtitle <- sprintf(
    "term: %s | target pes: %.3f | alpha: %.2f | sims/n: %d",
    x$term, x$target_pes, x$alpha, x$n_sims
  )

  p +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(
      name = "Total sample size (N)",
      sec.axis = ggplot2::sec_axis(
        ~ . / n_between_cells,
        name = "Subjects per between-cell"
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      labels = function(z) paste0(round(z * 100), "%")
    ) +
    ggplot2::labs(
      y        = "Simulated power",
      title    = "anovapowersim: simulated power curve",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title.x.top = ggplot2::element_text(colour = "grey35")
    )
}
