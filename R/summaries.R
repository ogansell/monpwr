#' Compute minimum detectable change (MDC) from power results
#'
#' @description
#' Derives the minimum detectable change (MDC) for each scenario × sim_type ×
#' group × horizon combination: the smallest effect size (% change per visit)
#' at which power reaches `power_target`.
#'
#' For hybrid scenarios, the `sim_type` column distinguishes between the
#' combined (legacy + new sites) and the legacy-only traces, allowing the
#' MDC improvement attributable to new sites to be read off directly.
#'
#' @param results A `monpwr_results` data frame from [run_power_sim()].
#' @param power_target Numeric scalar.  Target power level.  Default 0.80.
#'
#' @return A data frame with one row per scenario × sim_type × group × horizon,
#'   and columns: `scenario`, `label`, `sim_type`, `group`, `horizon`,
#'   `n_legacy`, `n_new`, `n_total`, `n_future`,
#'   `mdc_pct` (`NA` if target not achieved within tested effect sizes),
#'   `max_power`.
#'
#' @seealso [run_power_sim()], [plot_mdc()]
#' @export
compute_mdc <- function(results, power_target = 0.80) {
  stopifnot(inherits(results, "monpwr_results") || is.data.frame(results))

  # sim_type may not exist in results from older versions — default to "combined"
  if (!"sim_type" %in% names(results)) {
    results <- results |> mutate(sim_type = "combined")
  }
  if (!"n_legacy" %in% names(results)) {
    results <- results |> mutate(n_legacy = .data$n_plots, n_new = 0L,
                                 n_total = .data$n_plots)
  }

  results |>
    group_by(
      .data$scenario, .data$label, .data$sim_type,
      .data$group, .data$horizon
    ) |>
    summarise(
      mdc_pct   = {
        candidates <- .data$effect_pct[.data$power >= power_target]
        if (length(candidates) == 0) NA_real_ else min(candidates)
      },
      max_power = max(.data$power, na.rm = TRUE),
      n_legacy  = first(.data$n_legacy),
      n_new     = first(.data$n_new),
      n_total   = first(.data$n_total),
      n_future  = first(.data$n_future),
      .groups   = "drop"
    )
}


#' Bootstrap coefficient of variation of mean response
#'
#' @description
#' Estimates the precision of a state estimate for a set of sites using a
#' nonparametric bootstrap over sites.  The CV (coefficient of variation of
#' the bootstrap distribution of the mean) is a common design-adequacy
#' criterion (typical threshold: CV < 0.20).
#'
#' Only **legacy** sites with an observed historical record contribute to the
#' CV.  New sites in hybrid scenarios have no observed data and should not be
#' included in `site_ids`.
#'
#' @param data A data frame — the modelling dataset, one row per site × visit.
#' @param site_ids Character vector.  Site IDs to include.
#' @param n_boot Integer scalar.  Bootstrap replicates.  Default 1000.
#' @param place_var Character scalar.  Site ID column in `data`.
#'   Default `"Place"`.
#' @param season_var Character scalar.  Season/year column used to select
#'   the most recent monitoring round.  Default `"Season"`.  If not found
#'   in `data`, all rows for the retained sites are used.
#' @param response_fn A function applied to each site's rows (a data frame
#'   subset) that returns a single numeric value representing that site's
#'   state estimate.  Default: proportion of rows with `count > 0`
#'   (`function(d) mean(d$count > 0)`).  Supply a custom function for
#'   indices defined differently in your monitoring programme.
#'
#' @return A one-row data frame with columns: `n_sites`, `mean_resp`, `cv`.
#'   Returns `NA` for `mean_resp` and `cv` if fewer than 2 sites have data.
#'
#' @seealso [run_precision()], [plot_cv()]
#' @export
bootstrap_cv <- function(data, site_ids, n_boot = 1000L,
                         place_var   = "Place",
                         season_var  = "Season",
                         response_fn = function(d) mean(d$count > 0)) {

  if (!season_var %in% names(data)) {
    warn(paste0("`season_var` column '", season_var,
                "' not found in `data`. Using all rows for retained sites."))
    recent_dat <- data |> filter(.data[[place_var]] %in% site_ids)
  } else {
    recent_season <- max(data[[season_var]], na.rm = TRUE)
    recent_dat    <- data |>
      filter(.data[[place_var]] %in% site_ids,
             .data[[season_var]] == recent_season)
  }

  site_dat <- recent_dat |>
    group_by(.data[[place_var]]) |>
    summarise(resp = response_fn(pick(everything())), .groups = "drop")

  if (nrow(site_dat) < 2) {
    return(data.frame(n_sites = nrow(site_dat), mean_resp = NA_real_, cv = NA_real_))
  }

  boot_means <- map_dbl(seq_len(n_boot), function(i) {
    mean(slice_sample(site_dat, n = nrow(site_dat), replace = TRUE)$resp,
         na.rm = TRUE)
  })

  data.frame(
    n_sites   = nrow(site_dat),
    mean_resp = mean(site_dat$resp, na.rm = TRUE),
    cv        = sd(boot_means) / mean(boot_means)
  )
}


#' Run precision estimation across scenarios
#'
#' @description
#' Convenience wrapper that calls [bootstrap_cv()] for each scenario ×
#' reporting scale combination and returns a combined data frame suitable for
#' [plot_cv()].
#'
#' Only **legacy** sites contribute to CV estimates.  New sites in hybrid
#' scenarios are excluded because they have no observed historical record.
#'
#' @param scenarios Named list of `monpwr_scenario` objects from [scenario()].
#' @param data A data frame — the full modelling dataset.
#' @param plot_metadata A data frame with one row per legacy site.
#' @param n_boot Integer scalar.  Bootstrap replicates.  Default 1000.
#' @param reporting_groups Named list mapping display labels to column names
#'   in `plot_metadata` — the same list passed to [run_power_sim()].
#'   Default `list()` (overall scale only).
#' @param n_min_group Integer scalar.  Minimum retained sites for a group-level
#'   CV.  Default 5.
#' @param place_var Character scalar.  Site ID column.  Default `"Place"`.
#' @param ... Passed to [bootstrap_cv()].
#'
#' @return A data frame with columns `scenario`, `label`, `scale`, `group`,
#'   `n_sites`, `mean_resp`, `cv`.
#'
#' @seealso [bootstrap_cv()], [plot_cv()]
#' @export
run_precision <- function(scenarios,
                          data,
                          plot_metadata,
                          n_boot           = 1000L,
                          reporting_groups = list(),
                          n_min_group      = 5L,
                          place_var        = "Place",
                          ...) {

  map_dfr(names(scenarios), function(sc_name) {
    sc    <- scenarios[[sc_name]]
    sites <- sc$site_selector(plot_metadata)  # legacy sites only

    rows <- list()

    # Overall
    rows <- c(rows, list(
      bootstrap_cv(data, sites, n_boot, place_var = place_var, ...) |>
        mutate(scale = "Overall", group = "All")
    ))

    # Reporting group scales
    for (grp_label in names(reporting_groups)) {
      grp_col <- reporting_groups[[grp_label]]
      if (!grp_col %in% names(plot_metadata)) next

      levels_df <- plot_metadata |>
        filter(.data[[place_var]] %in% sites) |>
        count(.data[[grp_col]], name = "n") |>
        filter(.data$n >= n_min_group)

      for (i in seq_len(nrow(levels_df))) {
        lvl       <- levels_df[[grp_col]][i]
        grp_sites <- plot_metadata[[place_var]][
          plot_metadata[[grp_col]] == lvl &
            plot_metadata[[place_var]] %in% sites
        ]
        rows <- c(rows, list(
          bootstrap_cv(data, grp_sites, n_boot, place_var = place_var, ...) |>
            mutate(scale = grp_label, group = as.character(lvl))
        ))
      }
    }

    bind_rows(rows) |> mutate(scenario = sc_name, label = sc$label)
  })
}
