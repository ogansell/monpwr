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
      power_all   = vapply(.data$p_values, function(pv) {
        sum(pv < alpha, na.rm = TRUE) / length(pv)
      }, double(1)),
      power_all_lower = vapply(.data$p_values, function(pv) {
        binom.test(sum(pv < alpha, na.rm = TRUE), length(pv))$conf.int[1]
      }, double(1)),
      power_all_upper = vapply(.data$p_values, function(pv) {
        binom.test(sum(pv < alpha, na.rm = TRUE), length(pv))$conf.int[2]
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
      power_all   = vapply(.data$p_values, function(pv) {
        sum(pv < alpha, na.rm = TRUE) / length(pv)
      }, double(1)),
      power_all_lower = vapply(.data$p_values, function(pv) {
        binom.test(sum(pv < alpha, na.rm = TRUE), length(pv))$conf.int[1]
      }, double(1)),
      power_all_upper = vapply(.data$p_values, function(pv) {
        binom.test(sum(pv < alpha, na.rm = TRUE), length(pv))$conf.int[2]
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


#' Bracket variance-estimate uncertainty by scaling the random-effect SD
#'
#' @description
#' Variance components (`sigma_cond`, `sigma_zi`) estimated from a small pilot
#' are uncertain, and power is sensitive to them.  This helper rescales the
#' random-effect SD of a `monpwr_params` object by one or more multipliers so
#' the user can bracket how much that uncertainty moves the power estimate —
#' operationalising the sensitivity-analysis guidance in the package
#' documentation.
#'
#' Only the random-effect SD is scaled.  The dispersion parameter `disp_par`
#' is left unchanged, because it represents a different (observation-level)
#' source of uncertainty.
#'
#' @param ref_params A `monpwr_params` object from [extract_params()].
#' @param scales Numeric vector of strictly positive multipliers applied to
#'   `sigma_cond` (and `sigma_zi` when it is `> 0`).  Default `c(0.8, 1, 1.2)`.
#' @param run_fn Optional function of a single `monpwr_params` argument that
#'   returns a `monpwr_results` data frame (e.g. a partially-applied
#'   [run_power_sim()]).  If supplied, the helper runs it for each scale and
#'   row-binds the results with a numeric `sigma_scale` column.  If `NULL`
#'   (default), the helper returns a named list of scaled `monpwr_params`
#'   objects instead.
#'
#' @return If `run_fn` is `NULL`: a named list of `monpwr_params` objects, one
#'   per scale (names are the scale values).  Otherwise: a single data frame —
#'   the row-bound `run_fn` outputs with an added `sigma_scale` column.
#'
#' @seealso [calibrate_bias()], [run_power_sim()], [extract_params()]
#' @export
with_sigma_scaling <- function(ref_params,
                               scales = c(0.8, 1, 1.2),
                               run_fn = NULL) {
  stopifnot(inherits(ref_params, "monpwr_params"))
  if (!is.numeric(scales) || any(scales <= 0)) {
    abort(c(
      "`scales` must be a numeric vector of strictly positive multipliers.",
      i = "A scale of 0 would zero out the random-effect SD and break simulation."
    ))
  }

  scale_one <- function(s) {
    out <- ref_params
    out$sigma_cond <- ref_params$sigma_cond * s
    if (ref_params$sigma_zi > 0) {
      out$sigma_zi <- ref_params$sigma_zi * s
    }
    class(out) <- "monpwr_params"
    out
  }

  scaled <- lapply(scales, scale_one)
  names(scaled) <- as.character(scales)

  if (is.null(run_fn)) {
    return(scaled)
  }

  if (!is.function(run_fn)) {
    abort("`run_fn` must be a function of one argument (a monpwr_params object), or NULL.")
  }

  map_dfr(seq_along(scaled), function(i) {
    res <- run_fn(scaled[[i]])
    res |> mutate(sigma_scale = scales[i])
  })
}


#' Estimate parametric bias in power estimates
#'
#' @description
#' Measures the optimism introduced by conditioning on variance components
#' estimated from a finite pilot.  The truth arm uses the variance components
#' in `ref_params` as if they were the true population values.  The monpwr arm
#' refits a fresh pilot of `n_pilot` plots each replicate, re-extracts the
#' variance components, and conditions on those noisy estimates.  The
#' difference is the optimism from conditioning on variance components
#' estimated from a finite pilot, and it shrinks toward zero as `n_pilot`
#' grows.
#'
#' This is the in-package equivalent of Experiment 6 in the simr-comparison
#' study.
#'
#' @param ref_params A `monpwr_params` object from [extract_params()].
#' @param n_plots Integer scalar.  Number of plots in the design.
#' @param n_visits Integer scalar.
#'   Number of future visits (time points).
#' @param effect_pct Numeric scalar.  Effect size as percent change per
#'   visit (e.g. `3` for 3%).
#' @param n_cal Integer scalar.  Number of calibration replicates.
#'   Default 200.  Higher values give a tighter bias estimate but take
#'   longer.
#' @param alpha Numeric scalar.  Significance threshold.  Default 0.05.
#' @param test Character scalar.  `"wald"` (default) or `"lrt"`.
#' @param n_pilot Integer scalar or `NULL`.  Number of plots in each refit
#'   pilot.  Default `NULL` uses `nrow(ref_params$plot_state)` (the size
#'   of the pilot the params actually came from).
#' @param pilot_visits Integer scalar or `NULL`.  Visits per pilot plot.
#'   Default `NULL` uses `n_visits`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`monpwr_power`}{Power estimated by monpwr conditioning on
#'       re-estimated pilot parameters.}
#'     \item{`truth_power`}{Power using the true (known) variance components.}
#'     \item{`bias`}{`monpwr_power - truth_power`.  Positive = optimistic.}
#'     \item{`truth_ci`}{95% binomial CI on the ground-truth estimate.}
#'     \item{`n_cal`}{Number of calibration replicates used.}
#'     \item{`n_pilot`}{Pilot size used.}
#'     \item{`sigma_cond`}{Extracted sigma (for reference).}
#'   }
#'
#' @note This is the in-package equivalent of Experiment 6 in the
#'   simr-comparison study (`simr_extend_experiment.R`).
#'
#' @seealso [run_power_sim()], [extract_params()], [with_sigma_scaling()]
#' @export
calibrate_bias <- function(ref_params, n_plots, n_visits, effect_pct,
                           n_cal = 200L, alpha = 0.05,
                           test = c("wald", "lrt"),
                           n_pilot = NULL, pilot_visits = NULL) {
  stopifnot(inherits(ref_params, "monpwr_params"))
  test  <- match.arg(test)
  n_cal <- as.integer(n_cal)

  eff_log <- log(1 + effect_pct / 100)

  n_pilot      <- if (is.null(n_pilot)) nrow(ref_params$plot_state) else as.integer(n_pilot)
  pilot_visits <- if (is.null(pilot_visits)) n_visits else as.integer(pilot_visits)
  if (n_pilot < 2L) abort("`n_pilot` must be >= 2.")

  ps <- ref_params$plot_state
  marginal_int <- mean(ps$eta_last_cond - ps$blup_cond -
                         ref_params$beta_visit * ps$visit_num)

  # ---- TRUTH ARM: power using the TRUE variance components ----
  cli::cli_alert_info("Calibration: ground truth ({n_cal} reps, known variance)...")
  truth_state <- data.frame(
    plotid        = paste0("plot_", seq_len(n_plots)),
    visit_num     = 0L,
    eta_last_cond = marginal_int,
    eta_last_zi   = if (ref_params$sigma_zi > 0) {
      mean(ps$eta_last_zi - ps$blup_zi)
    } else 0,
    stringsAsFactors = FALSE
  )
  truth_pvals <- vapply(seq_len(n_cal), function(i) {
    dat <- simulate_visits(truth_state, n_visits, eff_log, ref_params,
                           draw_re = TRUE)
    fit_and_test(dat, ref_params, test = test)
  }, double(1))
  truth_conv  <- sum(!is.na(truth_pvals))
  truth_power <- if (truth_conv > 0) {
    sum(truth_pvals < alpha, na.rm = TRUE) / truth_conv
  } else NA_real_
  truth_ci <- if (truth_conv > 0) {
    stats::binom.test(sum(truth_pvals < alpha, na.rm = TRUE), truth_conv)$conf.int
  } else c(NA_real_, NA_real_)

  # ---- monpwr-AS-USED ARM: condition on sigma RE-ESTIMATED from a finite pilot ----
  cli::cli_alert_info("Calibration: monpwr conditioning on pilot estimates ({n_cal} reps)...")
  monpwr_pvals <- vapply(seq_len(n_cal), function(i) {
    pilot <- .simulate_pilot(marginal_int, eff_log, ref_params,
                             n_pilot, pilot_visits)
    ref_hat <- tryCatch(
      .refit_pilot_params(pilot, ref_params),
      error = function(e) NULL
    )
    if (is.null(ref_hat)) return(NA_real_)
    state_hat <- init_prospective_marginal(ref_hat, n_plots)
    dat <- simulate_visits(state_hat, n_visits, eff_log, ref_hat, draw_re = TRUE)
    fit_and_test(dat, ref_hat, test = test)
  }, double(1))
  monpwr_conv  <- sum(!is.na(monpwr_pvals))
  monpwr_power <- if (monpwr_conv > 0) {
    sum(monpwr_pvals < alpha, na.rm = TRUE) / monpwr_conv
  } else NA_real_

  bias <- monpwr_power - truth_power

  cli::cli_alert_info(
    "Bias: {round(bias * 100, 1)} pp (monpwr {round(monpwr_power * 100, 1)}% vs truth {round(truth_power * 100, 1)}%) | pilot n_plots = {n_pilot}"
  )

  list(
    monpwr_power = monpwr_power,
    truth_power  = truth_power,
    bias         = bias,
    truth_ci     = as.numeric(truth_ci),
    n_cal        = n_cal,
    n_pilot      = n_pilot,
    sigma_cond   = ref_params$sigma_cond
  )
}


.simulate_pilot <- function(marginal_int, eff_log, ref_params,
                            n_pilot, pilot_visits) {
  vnum  <- ref_params$visit_num_var
  pid   <- ref_params$plotid_var
  place <- ref_params$place_var
  resp  <- ref_params$count_var
  off   <- ref_params$offset_var

  re_cond <- rnorm(n_pilot, 0, ref_params$sigma_cond)
  re_zi   <- if (ref_params$sigma_zi > 0) {
    rnorm(n_pilot, 0, ref_params$sigma_zi)
  } else rep(0, n_pilot)

  marginal_zi <- if (ref_params$sigma_zi > 0) {
    ps <- ref_params$plot_state
    mean(ps$eta_last_zi - ps$blup_zi)
  } else 0

  rows <- vector("list", n_pilot)
  for (p in seq_len(n_pilot)) {
    vis   <- seq_len(pilot_visits)
    eta_c <- marginal_int + eff_log * vis + re_cond[p]
    eta_z <- marginal_zi + re_zi[p]
    counts <- .draw_counts(ref_params$family, exp(eta_c), plogis(eta_z),
                           pilot_visits, ref_params$disp_par)
    df <- data.frame(
      plot_id_tmp = paste0("pilot_", p),
      vnum_tmp    = vis,
      resp_tmp    = counts,
      stringsAsFactors = FALSE
    )
    if (!is.null(off)) df$off_tmp <- ref_params$log_effort_future
    rows[[p]] <- df
  }
  out <- do.call(rbind, rows)

  names(out)[names(out) == "plot_id_tmp"] <- pid
  names(out)[names(out) == "vnum_tmp"]    <- vnum
  names(out)[names(out) == "resp_tmp"]    <- resp
  out[[place]] <- out[[pid]]
  if (!is.null(off)) names(out)[names(out) == "off_tmp"] <- off
  out
}


.refit_pilot_params <- function(pilot, ref_params) {
  vnum  <- ref_params$visit_num_var
  pid   <- ref_params$plotid_var
  resp  <- ref_params$count_var
  off   <- ref_params$offset_var
  fam   <- ref_params$family

  off_term <- if (!is.null(off)) paste0(" + offset(", off, ")") else ""
  form_txt <- paste0(resp, " ~ ", vnum, off_term, " + (1 | ", pid, ")")

  suppress_w <- function(expr) {
    withCallingHandlers(expr,
      warning = function(w) invokeRestart("muffleWarning"))
  }

  fit <- switch(fam,
    poisson = suppress_w(lme4::glmer(stats::as.formula(form_txt),
                                     family = poisson, data = pilot)),
    binomial = suppress_w(lme4::glmer(
      stats::as.formula(paste0(resp, " ~ ", vnum, " + (1 | ", pid, ")")),
      family = binomial, data = pilot)),
    gaussian = suppress_w(lme4::lmer(
      stats::as.formula(paste0(resp, " ~ ", vnum, " + (1 | ", pid, ")")),
      data = pilot)),
    nbinom2 = suppress_w(glmmTMB::glmmTMB(stats::as.formula(form_txt),
                                          family = glmmTMB::nbinom2, data = pilot)),
    hurdle_nbinom2 = suppress_w(glmmTMB::glmmTMB(
      stats::as.formula(form_txt),
      ziformula = stats::as.formula(paste0("~ (1 | ", pid, ")")),
      family    = glmmTMB::truncated_nbinom2, data = pilot)),
    abort(paste0("calibrate_bias: no pilot refit implemented for family '", fam, "'."))
  )

  extract_params(
    fit,
    data          = pilot,
    visit_num_var = vnum,
    plotid_var    = pid,
    place_var     = ref_params$place_var,
    offset_var    = off
  )
}
