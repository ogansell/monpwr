#' Re-test power at a different alpha level
#'
#' @description
#' Recomputes power, CIs, and convergence stats from the stored raw p-values
#' without re-running any simulations.  Requires results produced by
#' [run_power_sim()] (which stores p-values as a list-column).
#'
#' @param results A `monpwr_results` data frame from [run_power_sim()].
#' @param alpha Numeric scalar.  New significance threshold.
#'
#' @return A `monpwr_results` data frame with updated `power`, `power_lower`,
#'   `power_upper` columns.  The `p_values` list-column is preserved.
#'
#' @seealso [run_power_sim()], [extend()]
#' @export
retest <- function(results, alpha) {
  if (!"p_values" %in% names(results)) {
    abort(c(
      "Cannot retest: `p_values` column not found.",
      i = "Raw p-values are only stored by `run_power_sim()` >= 0.4.0."
    ))
  }
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1)

  out <- results |>
    mutate(
      n_converged = vapply(.data$p_values, function(pv) sum(!is.na(pv)), integer(1)),
      power       = vapply(.data$p_values, function(pv) {
        nc <- sum(!is.na(pv))
        if (nc == 0) return(NA_real_)
        sum(pv < alpha, na.rm = TRUE) / nc
      }, double(1)),
      power_lower = vapply(.data$p_values, function(pv) {
        nc <- sum(!is.na(pv))
        if (nc == 0) return(NA_real_)
        binom.test(sum(pv < alpha, na.rm = TRUE), nc)$conf.int[1]
      }, double(1)),
      power_upper = vapply(.data$p_values, function(pv) {
        nc <- sum(!is.na(pv))
        if (nc == 0) return(NA_real_)
        binom.test(sum(pv < alpha, na.rm = TRUE), nc)$conf.int[2]
      }, double(1)),
      conv_rate   = round(.data$n_converged /
        vapply(.data$p_values, length, integer(1)), 3)
    )

  class(out) <- union("monpwr_results", class(out))
  attr(out, "alpha") <- alpha
  out
}


#' Extend power results with additional iterations
#'
#' @description
#' Combines two `monpwr_results` data frames by concatenating their raw
#' p-values for matching cells (same scenario × group × effect × horizon),
#' then recomputes power and CIs from the pooled replicates.  This lets you
#' accumulate iterations incrementally rather than re-running from scratch.
#'
#' Both inputs must have a `p_values` list-column.  Cells present in only
#' one input are kept as-is.
#'
#' @param results A `monpwr_results` data frame.
#' @param additional A `monpwr_results` data frame with the same structure
#'   (typically from a second [run_power_sim()] call with different seeds).
#' @param alpha Numeric scalar.  Significance threshold for the re-summary.
#'   Default uses `attr(results, "alpha")`, falling back to 0.10.
#'
#' @return A `monpwr_results` data frame with pooled p-values and updated
#'   power estimates.
#'
#' @seealso [run_power_sim()], [retest()]
#' @export
extend <- function(results, additional, alpha = NULL) {
  for (x in list(results, additional)) {
    if (!"p_values" %in% names(x)) {
      abort(c(
        "Cannot extend: `p_values` column not found.",
        i = "Raw p-values are only stored by `run_power_sim()` >= 0.4.0."
      ))
    }
  }

  alpha <- alpha %||% attr(results, "alpha") %||% 0.10
  key_cols <- c("scenario", "label", "group", "effect_pct", "horizon")

  merged <- dplyr::full_join(
    results  |> select(all_of(key_cols), p_values_a = "p_values",
                       "n_plots", "n_future"),
    additional |> select(all_of(key_cols), p_values_b = "p_values",
                         "n_plots", "n_future"),
    by = c(key_cols, "n_plots", "n_future"),
    suffix = c("", ".b")
  ) |>
    mutate(
      p_values = mapply(function(a, b) {
        c(if (!is.null(a)) a, if (!is.null(b)) b)
      }, .data$p_values_a, .data$p_values_b, SIMPLIFY = FALSE)
    ) |>
    select(-"p_values_a", -"p_values_b")

  merged <- merged |>
    mutate(
      n_converged = vapply(.data$p_values, function(pv) sum(!is.na(pv)), integer(1)),
      power       = vapply(.data$p_values, function(pv) {
        nc <- sum(!is.na(pv))
        if (nc == 0) return(NA_real_)
        sum(pv < alpha, na.rm = TRUE) / nc
      }, double(1)),
      power_lower = vapply(.data$p_values, function(pv) {
        nc <- sum(!is.na(pv))
        if (nc == 0) return(NA_real_)
        binom.test(sum(pv < alpha, na.rm = TRUE), nc)$conf.int[1]
      }, double(1)),
      power_upper = vapply(.data$p_values, function(pv) {
        nc <- sum(!is.na(pv))
        if (nc == 0) return(NA_real_)
        binom.test(sum(pv < alpha, na.rm = TRUE), nc)$conf.int[2]
      }, double(1)),
      conv_rate   = round(.data$n_converged /
        vapply(.data$p_values, length, integer(1)), 3)
    )

  class(merged) <- union("monpwr_results", class(merged))
  attr(merged, "alpha") <- alpha
  merged
}


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


#' Estimate parametric bias in power estimates
#'
#' @description
#' Compares monpwr's prospective power estimate against a brute-force
#' ground truth at a single design point.  The brute-force estimator
#' generates fresh data from `ref_params` each replicate (re-drawing
#' random effects and observation noise), fits the test model, and
#' collects p-values — identical to monpwr's DGP but without conditioning
#' on a fixed variance estimate.
#'
#' The difference between monpwr and brute-force power is the parametric
#' bias: the optimism introduced by conditioning on estimated variance
#' components from finite data.  This bias is approximately constant
#' across the design grid for a given `ref_params`, so calibrating at
#' one point is sufficient.
#'
#' @param ref_params A `monpwr_params` object from [extract_params()].
#' @param n_plots Integer scalar.  Number of plots in the design.
#' @param n_visits Integer scalar.  Number of future visits (time points).
#' @param effect_pct Numeric scalar.  Effect size as percent change per
#'   visit (e.g. `3` for 3%).
#' @param n_cal Integer scalar.  Number of calibration replicates.
#'   Default 200.  Higher values give a tighter bias estimate but take
#'   longer.
#' @param alpha Numeric scalar.  Significance threshold.  Default 0.05.
#' @param test Character scalar.  `"wald"` (default) or `"lrt"`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`monpwr_power`}{Power estimated by monpwr prospective simulation.}
#'     \item{`truth_power`}{Brute-force ground-truth power.}
#'     \item{`bias`}{`monpwr_power - truth_power`.  Positive = optimistic.}
#'     \item{`truth_ci`}{95% binomial CI on the ground-truth estimate.}
#'     \item{`n_cal`}{Number of calibration replicates used.}
#'     \item{`sigma_cond`}{Extracted sigma (for reference).}
#'   }
#'
#' @seealso [run_power_sim()], [extract_params()]
#' @export
calibrate_bias <- function(ref_params, n_plots, n_visits, effect_pct,
                           n_cal = 200L, alpha = 0.05, test = c("wald", "lrt")) {
  stopifnot(inherits(ref_params, "monpwr_params"))
  test <- match.arg(test)
  n_cal <- as.integer(n_cal)

  eff_log <- log(1 + effect_pct / 100)

  # --- monpwr prospective ---
  plot_state <- init_prospective_marginal(ref_params, n_plots)

  cli::cli_alert_info("Calibration: monpwr ({n_cal} reps)...")
  monpwr_pvals <- vapply(seq_len(n_cal), function(i) {
    future_dat <- simulate_visits(plot_state, n_visits, eff_log, ref_params,
                                  draw_re = TRUE)
    fit_and_test(future_dat, ref_params, test = test)
  }, double(1))

  monpwr_conv  <- sum(!is.na(monpwr_pvals))
  monpwr_power <- if (monpwr_conv > 0) sum(monpwr_pvals < alpha, na.rm = TRUE) / monpwr_conv else NA_real_

  # --- brute-force ground truth ---
  cli::cli_alert_info("Calibration: brute-force ({n_cal} reps)...")

  ps <- ref_params$plot_state
  marginal_int <- mean(ps$eta_last_cond - ps$blup_cond -
                         ref_params$beta_visit * ps$visit_num)
  marginal_zi  <- if (ref_params$sigma_zi > 0) {
    mean(ps$eta_last_zi - ps$blup_zi)
  } else {
    0
  }

  bf_pvals <- vapply(seq_len(n_cal), function(i) {
    re_cond <- rnorm(n_plots, 0, ref_params$sigma_cond)
    re_zi   <- if (ref_params$sigma_zi > 0) {
      rnorm(n_plots, 0, ref_params$sigma_zi)
    } else {
      rep(0, n_plots)
    }

    rows <- vector("list", n_plots)
    for (p in seq_len(n_plots)) {
      vis <- seq_len(n_visits)
      eta_c <- marginal_int + eff_log * vis + re_cond[p]
      eta_z <- marginal_zi + re_zi[p]
      mu  <- exp(eta_c)
      pzi <- plogis(eta_z)
      counts <- .draw_counts(ref_params$family, mu, pzi, n_visits,
                              ref_params$disp_par)
      rows[[p]] <- data.frame(
        plotid     = paste0("plot_", p),
        visit_num  = vis,
        log_effort = ref_params$log_effort_future,
        count      = counts,
        source     = "future",
        stringsAsFactors = FALSE
      )
    }
    bf_dat <- do.call(rbind, rows)
    fit_and_test(bf_dat, ref_params, test = test)
  }, double(1))

  bf_conv  <- sum(!is.na(bf_pvals))
  bf_power <- if (bf_conv > 0) sum(bf_pvals < alpha, na.rm = TRUE) / bf_conv else NA_real_
  bf_ci    <- if (bf_conv > 0) {
    stats::binom.test(sum(bf_pvals < alpha, na.rm = TRUE), bf_conv)$conf.int
  } else {
    c(NA_real_, NA_real_)
  }

  bias <- monpwr_power - bf_power

  cli::cli_alert_info(
    "Bias: {round(bias * 100, 1)} percentage points (monpwr {round(monpwr_power * 100, 1)}% vs truth {round(bf_power * 100, 1)}%)"
  )

  list(
    monpwr_power = monpwr_power,
    truth_power  = bf_power,
    bias         = bias,
    truth_ci     = as.numeric(bf_ci),
    n_cal        = n_cal,
    sigma_cond   = ref_params$sigma_cond
  )
}
