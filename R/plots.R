#' Plot power curves by scenario
#'
#' @description
#' Bar chart of power by effect size, faceted by monitoring horizon.  One
#' panel per horizon, bars grouped by scenario.
#'
#' For hybrid scenarios the output contains both `"combined"` (legacy + new
#' sites) and `"legacy_only"` rows.  Use `sim_type_filter` to select which
#' to display, or call [plot_power_gain()] to show both traces together.
#'
#' @param results A `monpwr_results` data frame from [run_power_sim()].
#' @param scale_filter Character scalar.  Reporting scale to display.
#'   Default `"Overall"`.
#' @param group_filter Character scalar or `NULL`.  Group within the scale.
#'   `NULL` includes all groups.
#' @param sim_type_filter Character scalar.  Which simulation type to display:
#'   `"combined"` (default) or `"legacy_only"`.
#' @param power_target Numeric scalar.  Reference line for target power.
#'   Default 0.80.
#' @param alpha Numeric scalar.  Passed to the caption.  Default 0.10.
#' @param title Character scalar.  Plot title.  Default auto-generated.
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_power_gain()], [plot_mdc()], [plot_cv()], [run_power_sim()]
#' @export
plot_power <- function(results,
                       scale_filter    = "Overall",
                       group_filter    = NULL,
                       sim_type_filter = "combined",
                       power_target    = 0.80,
                       alpha           = 0.10,
                       title           = NULL) {

  # Back-compat: older results without sim_type column
  if (!"sim_type" %in% names(results)) {
    results <- results |> mutate(sim_type = "combined")
  }

  dat <- results |>
    filter(.data$scale    == scale_filter,
           .data$sim_type == sim_type_filter,
           !is.na(.data$power))

  if (!is.null(group_filter)) dat <- dat |> filter(.data$group == group_filter)

  if (nrow(dat) == 0) {
    warn("No data to plot after filtering. Check `scale_filter`, `group_filter`, and `sim_type_filter`.")
    return(invisible(NULL))
  }

  n_iter_approx <- if (all(c("n_converged", "conv_rate") %in% names(dat))) {
    round(mean(dat$n_converged / pmax(dat$conv_rate, 1e-6), na.rm = TRUE))
  } else NA

  ggplot2::ggplot(
    dat,
    ggplot2::aes(x    = factor(.data$effect_pct),
                 y    = .data$power,
                 fill = .data$label)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.7) +
    ggplot2::geom_hline(yintercept = power_target, linetype = "dashed",
                        colour = "firebrick", linewidth = 0.7) +
    ggplot2::facet_wrap(~ paste0("Horizon: ", .data$horizon, " yr")) +
    ggplot2::scale_fill_viridis_d(end = 0.85, name = "Scenario") +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(
      x       = "Effect size (% change per visit)",
      y       = "Power",
      title   = title %||% paste0("Power by scenario - ", scale_filter,
                                  " [", sim_type_filter, "]"),
      caption = paste0("Alpha = ", alpha,
                       if (!is.na(n_iter_approx))
                         paste0("  |  ~", n_iter_approx, " replicates per cell")
                       else "")
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}


#' Plot hybrid power gain - combined vs legacy-only traces
#'
#' @description
#' For hybrid scenarios (those with `n_new_sites > 0`), shows both the
#' combined (legacy + new sites) and the legacy-only power curves on the same
#' panel.  The vertical gap between the solid and dashed lines at any effect
#' size is the power gain attributable solely to the new sites.
#'
#' Non-hybrid scenarios are silently dropped unless `include_legacy_only` is
#' `TRUE`, in which case they appear as a single solid line (combined =
#' legacy-only).
#'
#' @param results A `monpwr_results` data frame from [run_power_sim()].
#' @param scale_filter Character scalar.  Reporting scale to display.
#'   Default `"Overall"`.
#' @param group_filter Character scalar or `NULL`.  Group within the scale.
#' @param power_target Numeric scalar.  Reference line.  Default 0.80.
#' @param title Character scalar.  Plot title.  Default auto-generated.
#'
#' @return A `ggplot` object, or `NULL` (invisibly) if no hybrid scenarios
#'   are found after filtering.
#'
#' @seealso [plot_power()], [run_power_sim()]
#' @export
plot_power_gain <- function(results,
                            scale_filter  = "Overall",
                            group_filter  = NULL,
                            power_target  = 0.80,
                            title         = NULL) {

  if (!"sim_type" %in% names(results)) {
    warn("No `sim_type` column found - this plot requires hybrid scenario results.")
    return(invisible(NULL))
  }

  dat <- results |>
    filter(.data$scale == scale_filter, !is.na(.data$power))
  if (!is.null(group_filter)) dat <- dat |> filter(.data$group == group_filter)

  # Keep only scenarios that have both sim_types (i.e. true hybrid scenarios)
  hybrid_labels <- dat |>
    group_by(.data$label) |>
    summarise(has_both = n_distinct(.data$sim_type) > 1, .groups = "drop") |>
    filter(.data$has_both) |>
    pull(.data$label)

  dat <- dat |> filter(.data$label %in% hybrid_labels)

  if (nrow(dat) == 0) {
    warn("No hybrid scenarios found after filtering.")
    return(invisible(NULL))
  }

  dat <- dat |>
    mutate(
      trace_label = if_else(.data$sim_type == "combined",
                            paste0(.data$label, " (combined)"),
                            paste0(.data$label, " (legacy only)"))
    )

  ggplot2::ggplot(
    dat,
    ggplot2::aes(x        = factor(.data$effect_pct),
                 y        = .data$power,
                 colour   = .data$label,
                 linetype = .data$sim_type,
                 group    = interaction(.data$label, .data$sim_type))
  ) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(.data$label, .data$sim_type))) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_hline(yintercept = power_target, linetype = "dashed",
                        colour = "firebrick", linewidth = 0.7) +
    ggplot2::facet_wrap(~ paste0("Horizon: ", .data$horizon, " yr")) +
    ggplot2::scale_colour_viridis_d(end = 0.85, name = "Scenario") +
    ggplot2::scale_linetype_manual(
      values = c(combined = "solid", legacy_only = "dashed"),
      labels = c(combined = "Combined (legacy + new)", legacy_only = "Legacy only"),
      name   = NULL
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(
      x     = "Effect size (% change per visit)",
      y     = "Power",
      title = title %||% paste0("Power gain from new sites - ", scale_filter)
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}


#' MDC heatmap by scenario and horizon
#'
#' @description
#' Tile heatmap of minimum detectable change (MDC) for each scenario x
#' monitoring horizon.  Darker tiles = lower MDC = more sensitive design.
#' Red tiles indicate that the power target was not achieved within the
#' tested effect sizes; the maximum power achieved is shown instead.
#'
#' @param mdc A data frame from [compute_mdc()].
#' @param scale_filter Character scalar.  Reporting scale to display.
#'   Default `"Overall"`.
#' @param group_filter Character scalar or `NULL`.  Group to display.
#' @param sim_type_filter Character scalar.  `"combined"` (default) or
#'   `"legacy_only"`.
#' @param power_target Numeric scalar.  Used in the subtitle.  Default 0.80.
#' @param alpha Numeric scalar.  Used in the caption.  Default 0.10.
#' @param title Character scalar.  Plot title.  Default auto-generated.
#'
#' @return A `ggplot` object.
#'
#' @seealso [compute_mdc()], [plot_power()], [plot_cv()]
#' @export
plot_mdc <- function(mdc,
                     scale_filter    = "Overall",
                     group_filter    = NULL,
                     sim_type_filter = "combined",
                     power_target    = 0.80,
                     alpha           = 0.10,
                     title           = NULL) {

  if (!"sim_type" %in% names(mdc)) mdc <- mdc |> mutate(sim_type = "combined")

  dat <- mdc |>
    filter(.data$scale    == scale_filter,
           .data$sim_type == sim_type_filter)
  if (!is.null(group_filter)) dat <- dat |> filter(.data$group == group_filter)

  if (nrow(dat) == 0) {
    warn("No data to plot after filtering.")
    return(invisible(NULL))
  }

  dat <- dat |>
    mutate(
      label_txt = if_else(
        is.na(.data$mdc_pct),
        paste0(round(.data$max_power * 100), "%\nmax"),
        paste0(.data$mdc_pct, "%")
      )
    )

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x     = factor(.data$horizon),
      y     = reorder(.data$label, -.data$mdc_pct, na.rm = TRUE),
      fill  = .data$mdc_pct,
      label = .data$label_txt
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(size = 3.5, colour = "white", fontface = "bold") +
    ggplot2::scale_fill_viridis_c(
      name = "MDC\n(%/visit)", option = "D",
      direction = -1, na.value = "#d73027"
    ) +
    ggplot2::labs(
      x        = "Monitoring horizon (yr)",
      y        = NULL,
      title    = title %||% paste0("MDC by scenario and horizon - ", scale_filter,
                                   " [", sim_type_filter, "]"),
      subtitle = paste0("Darker = lower MDC = more sensitive  |  ",
                        "Red = ", round(power_target * 100), "% power not achieved"),
      caption  = paste0("Alpha = ", alpha)
    ) +
    ggplot2::theme_bw(base_size = 11)
}


#' CV bar chart by scenario
#'
#' @description
#' Horizontal bar chart of bootstrap CV of mean response by scenario.  A
#' dashed reference line at `cv_target` marks the adequacy threshold.
#'
#' @param precision A data frame from [run_precision()].
#' @param scale_filter Character scalar.  Reporting scale to display.
#'   Default `"Overall"`.
#' @param group_filter Character scalar or `NULL`.  Group to display.
#' @param cv_target Numeric scalar.  CV adequacy threshold reference line.
#'   Default 0.20.
#' @param title Character scalar.  Plot title.  Default auto-generated.
#'
#' @return A `ggplot` object.
#'
#' @seealso [run_precision()], [bootstrap_cv()], [plot_power()]
#' @export
plot_cv <- function(precision,
                    scale_filter = "Overall",
                    group_filter = NULL,
                    cv_target    = 0.20,
                    title        = NULL) {

  dat <- precision |>
    filter(!is.na(.data$cv), .data$scale == scale_filter)
  if (!is.null(group_filter)) dat <- dat |> filter(.data$group == group_filter)

  if (nrow(dat) == 0) {
    warn("No data to plot after filtering.")
    return(invisible(NULL))
  }

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x    = reorder(.data$label, .data$cv),
      y    = .data$cv,
      fill = .data$label
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = cv_target, colour = "firebrick",
                        linetype = "dashed", linewidth = 0.7) +
    ggplot2::scale_fill_viridis_d(end = 0.85, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(
      x       = NULL,
      y       = "CV of mean response",
      title   = title %||% paste0("Precision of state estimates - ", scale_filter),
      caption = paste0("Red dashed = CV ", cv_target, " adequacy threshold")
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11)
}
