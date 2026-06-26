#' Simulate future visits for a set of plots
#'
#' @description
#' The core simulation engine.  Given a plot-state data frame (from either
#' [init_conditional()] or [init_prospective_marginal()] / a custom `init_fn`),
#' simulates `n_future` additional visits per plot under a specified effect size,
#' using the parameters in `ref_params`.
#'
#' Counts are drawn from the appropriate distribution for the model family
#' (`ref_params$family`):
#' * `"hurdle_nbinom2"` — truncated NB2 for occupied plots, Bernoulli for occupancy
#' * `"nbinom2"` — NB2 (including zeros)
#' * `"poisson"` — Poisson
#' * `"binomial"` — Bernoulli (occupancy only; count = 0 or 1)
#' * `"gaussian"` — Normal with SD = `ref_params$disp_par`
#'
#' @param plot_state Data frame, one row per plot.  Must contain:
#'   `plotid` (character), `visit_num` (integer), `eta_last_cond` (numeric),
#'   `eta_last_zi` (numeric).  Produced by [init_conditional()],
#'   [init_prospective_marginal()], or a custom `init_fn`.
#' @param n_future Integer scalar.  Number of additional visits to simulate
#'   beyond each plot's current `visit_num`.
#' @param eff_log Numeric scalar.  The hypothetical effect size on the
#'   log scale: `log(1 + effect_pct / 100)`.  This is the trend coefficient
#'   applied to future visits.
#' @param ref_params A `monpwr_params` object from [extract_params()].
#' @param draw_re Logical.  If `TRUE` (default for prospective mode), fresh
#'   random effects are drawn from `N(0, sigma_cond)` and `N(0, sigma_zi)`
#'   for each plot.  If `FALSE` (conditional mode), random effects are already
#'   absorbed into `eta_last_cond` / `eta_last_zi` and are not re-drawn.
#'
#' @return A data frame with one row per plot × future visit, containing:
#'   `plotid`, `visit_num`, `log_effort` (`ref_params$log_effort_future` —
#'   median observed effort, or `0` if no offset), `count`, `source` (`"future"`).
#'
#' @seealso [build_historical()], [fit_and_test()], [run_power_sim()]
#' @export
simulate_visits <- function(plot_state, n_future, eff_log, ref_params,
                            draw_re = FALSE) {
  stopifnot(
    is.data.frame(plot_state),
    n_future >= 1L,
    inherits(ref_params, "monpwr_params")
  )

  n_plots  <- nrow(plot_state)
  fam      <- ref_params$family
  draw_re  <- isTRUE(draw_re)

  # Draw plot-level random effects if prospective mode
  re_cond <- if (draw_re) rnorm(n_plots, 0, ref_params$sigma_cond) else rep(0, n_plots)
  re_zi   <- if (draw_re && ref_params$sigma_zi > 0) {
    rnorm(n_plots, 0, ref_params$sigma_zi)
  } else {
    rep(0, n_plots)
  }

  rows <- vector("list", n_plots)

  for (i in seq_len(n_plots)) {
    ps           <- plot_state[i, ]
    future_vis   <- seq(ps$visit_num + 1L, ps$visit_num + n_future)

    # Linear predictor for future visits:
    eta_c <- ps$eta_last_cond +
      eff_log * (future_vis - ps$visit_num) +
      ref_params$beta_gap_cond * ref_params$visit_gap_med +
      re_cond[i]

    eta_z <- ps$eta_last_zi +
      ref_params$beta_gap_zi * ref_params$visit_gap_med +
      re_zi[i]

    mu  <- exp(eta_c)
    pzi <- plogis(eta_z)   # P(structural zero)

    counts <- .draw_counts(fam, mu, pzi, n_future, ref_params$disp_par)

    rows[[i]] <- data.frame(
      plotid     = as.character(ps$plotid),
      visit_num  = future_vis,
      log_effort = ref_params$log_effort_future,
      count      = counts,
      source     = "future",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}


# ------------------------------------------------------------------------------
# Internal count sampler — dispatches on family
# ------------------------------------------------------------------------------

.draw_counts <- function(family, mu, pzi, n, disp) {
  switch(family,

    hurdle_nbinom2 = {
      occupied <- rbinom(n, 1, 1 - pzi)
      counts   <- integer(n)
      for (j in seq_len(n)) {
        if (occupied[j] == 0L) next
        repeat {
          x <- rnbinom(1, size = disp, mu = mu[j])
          if (x > 0L) { counts[j] <- x; break }
        }
      }
      counts
    },

    nbinom2 = rnbinom(n, size = disp, mu = mu),

    poisson = rpois(n, lambda = mu),

    binomial = rbinom(n, 1, plogis(log(mu))),  # mu = exp(eta); log(mu) = eta recovers the logit -> probability

    gaussian = rnorm(n, mean = mu, sd = disp),

    # Unknown family — warn once and treat as Poisson
    {
      warn(paste0("Unknown family '", family, "'; simulating as Poisson."))
      rpois(n, lambda = mu)
    }
  )
}

# NULL-coalescing operator (avoid dependency on rlang::`%||%` at call sites)
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Build the historical data stub for retained plots
#'
#' @description
#' Extracts the observed monitoring record for a set of retained plots from
#' the full modelling dataset.  The historical stub is stacked with simulated
#' future data in [run_power_sim()] (conditional mode) so that the trend model
#' sees the full temporal record — existing visits + future visits.
#'
#' In prospective mode this function is not called; the trend model sees only
#' simulated future data.
#'
#' @param data A data frame — the full modelling dataset, one row per plot ×
#'   visit.  Must contain `Place`, `visit_num`, `log_effort`, and `count`
#'   columns (or the equivalents specified via `place_var` etc. in
#'   [extract_params()]; the defaults are used here).
#' @param site_ids Character vector.  Plot `Place` IDs to retain.
#' @param place_var Character scalar.  Name of the plot ID column in `data`.
#'   Default `"Place"`.
#' @param plotid_var Character scalar.  Name of the model grouping column.
#'   Default `"plotid_model"`.  This becomes the `plotid` column in the
#'   output (used by [fit_and_test()]).
#'
#' @return A data frame with columns `plotid`, `visit_num`, `log_effort`,
#'   `count`, `source` (`"observed"`).
#'
#' @seealso [simulate_visits()], [run_power_sim()]
#' @export
build_historical <- function(data, site_ids,
                             place_var     = "Place",
                             plotid_var    = "plotid_model",
                             visit_num_var = "visit_num",
                             count_var     = "count",
                             offset_var    = NULL) {
  log_eff_vals <- if (!is.null(offset_var) && offset_var %in% names(data)) {
    data[[offset_var]]
  } else {
    rep(0, nrow(data))
  }

  data |>
    mutate(.log_effort = log_eff_vals) |>
    filter(.data[[place_var]] %in% site_ids) |>
    transmute(
      plotid     = as.character(.data[[plotid_var]]),
      visit_num  = .data[[visit_num_var]],
      log_effort = .data$.log_effort,
      count      = .data[[count_var]],
      source     = "observed"
    )
}


#' Fit a trend model and return the p-value for the visit coefficient
#'
#' @description
#' Fits a mixed model to the combined observed + simulated data and returns
#' the p-value for the visit-sequence coefficient (`visit_num`).
#' This is the test statistic used in power estimation.
#'
#' Two test types are supported:
#'
#' * **`"wald"`** (default) — extracts the z/t-statistic from model
#'   coefficients.  Fast.
#' * **`"lrt"`** — likelihood-ratio test comparing the full model to a
#'   reduced model without the visit coefficient.  More reliable at small
#'   sample sizes but approximately twice as slow (two model fits per
#'   replicate).
#'
#' The model family used for testing is chosen based on `ref_params$family`:
#' * `"hurdle_nbinom2"` — `glmmTMB` hurdle (truncated NB2 + ZI random intercept)
#' * `"nbinom2"` — `glmmTMB` NB2
#' * `"poisson"` — `glmer` Poisson
#' * `"binomial"` — `glmer` binomial
#' * `"gaussian"` — `lmer` Gaussian
#'
#' Convergence failures return `NA_real_` silently; the proportion of
#' `NA` values is tracked as `conv_rate` in [run_power_sim()] output.
#'
#' @param data A data frame — the combined observed + simulated dataset for
#'   one simulation replicate.  Must contain `plotid`, `visit_num`,
#'   `log_effort`, `count`.
#' @param ref_params A `monpwr_params` object.  Used to select the
#'   appropriate model family for the test.
#' @param test Character scalar.  `"wald"` (default) or `"lrt"`.
#'
#' @return Numeric scalar — the p-value for `visit_num`, or `NA_real_`
#'   on convergence failure.
#'
#' @details
#' The test model is intentionally simpler than the data-generating model
#' (no fixed effects for park, forest, spatial terms) — these are absorbed
#' into the plot-level random intercept.  This gives conservative power
#' estimates and is consistent with typical FPI reporting practice.
#'
#' @seealso [simulate_visits()], [run_power_sim()]
#' @export
fit_and_test <- function(data, ref_params, test = c("wald", "lrt")) {
  stopifnot(inherits(ref_params, "monpwr_params"))
  test <- match.arg(test)

  withCallingHandlers(
    tryCatch({
      fit <- .fit_test_model(data, ref_params$family)
      if (test == "lrt") {
        fit_null <- .fit_null_model(data, ref_params$family)
        .extract_lrt_pval(fit, fit_null)
      } else {
        .extract_pval(fit, ref_params$family)
      }
    },
    error = function(e) NA_real_),
    warning = function(w) invokeRestart("muffleWarning")
  )
}


# ------------------------------------------------------------------------------
# Internal: fit the test model
# ------------------------------------------------------------------------------

.fit_test_model <- function(data, family) {
  suppress_w <- function(expr) {
    withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
  }

  switch(family,

    hurdle_nbinom2 = suppress_w(
      glmmTMB::glmmTMB(
        formula   = count ~ visit_num + offset(log_effort) + (1 | plotid),
        ziformula =       ~ (1 | plotid),
        family    = glmmTMB::truncated_nbinom2,
        data      = data,
        control   = glmmTMB::glmmTMBControl(
          optCtrl = list(iter.max = 200, eval.max = 300))
      )
    ),

    nbinom2 = suppress_w(
      glmmTMB::glmmTMB(
        count ~ visit_num + offset(log_effort) + (1 | plotid),
        family  = glmmTMB::nbinom2,
        data    = data,
        control = glmmTMB::glmmTMBControl(
          optCtrl = list(iter.max = 200, eval.max = 300))
      )
    ),

    poisson = suppress_w(
      lme4::glmer(
        count ~ visit_num + offset(log_effort) + (1 | plotid),
        family = poisson, data = data
      )
    ),

    binomial = suppress_w(
      lme4::glmer(
        count ~ visit_num + (1 | plotid),
        family = binomial, data = data
      )
    ),

    gaussian = suppress_w(
      lme4::lmer(
        count ~ visit_num + (1 | plotid),
        data = data
      )
    ),

    abort(paste0("No test model implemented for family '", family, "'."))
  )
}


# ------------------------------------------------------------------------------
# Internal: extract p-value from test model
# ------------------------------------------------------------------------------

.extract_pval <- function(fit, family) {
  if (inherits(fit, "glmmTMB")) {
    ct <- summary(fit)$coefficients$cond
    if (!"visit_num" %in% rownames(ct)) return(NA_real_)
    return(ct["visit_num", "Pr(>|z|)"])
  }

  if (inherits(fit, c("glmerMod", "lmerMod"))) {
    ct <- summary(fit)$coefficients
    if (!"visit_num" %in% rownames(ct)) return(NA_real_)
    col <- if ("Pr(>|z|)" %in% colnames(ct)) "Pr(>|z|)" else "Pr(>|t|)"
    return(ct["visit_num", col])
  }

  NA_real_
}


# ------------------------------------------------------------------------------
# Internal: fit null model (no visit_num) for LRT
# ------------------------------------------------------------------------------

.fit_null_model <- function(data, family) {
  suppress_w <- function(expr) {
    withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
  }

  switch(family,

    hurdle_nbinom2 = suppress_w(
      glmmTMB::glmmTMB(
        formula   = count ~ 1 + offset(log_effort) + (1 | plotid),
        ziformula =       ~ (1 | plotid),
        family    = glmmTMB::truncated_nbinom2,
        data      = data,
        control   = glmmTMB::glmmTMBControl(
          optCtrl = list(iter.max = 200, eval.max = 300))
      )
    ),

    nbinom2 = suppress_w(
      glmmTMB::glmmTMB(
        count ~ 1 + offset(log_effort) + (1 | plotid),
        family  = glmmTMB::nbinom2,
        data    = data,
        control = glmmTMB::glmmTMBControl(
          optCtrl = list(iter.max = 200, eval.max = 300))
      )
    ),

    poisson = suppress_w(
      lme4::glmer(
        count ~ 1 + offset(log_effort) + (1 | plotid),
        family = poisson, data = data
      )
    ),

    binomial = suppress_w(
      lme4::glmer(
        count ~ 1 + (1 | plotid),
        family = binomial, data = data
      )
    ),

    gaussian = suppress_w(
      lme4::lmer(
        count ~ 1 + (1 | plotid),
        data = data
      )
    ),

    abort(paste0("No null model implemented for family '", family, "'."))
  )
}


# ------------------------------------------------------------------------------
# Internal: extract LRT p-value from full vs null model comparison
# ------------------------------------------------------------------------------

.extract_lrt_pval <- function(fit_full, fit_null) {
  if (inherits(fit_full, "glmmTMB")) {
    ll_full <- as.numeric(stats::logLik(fit_full))
    ll_null <- as.numeric(stats::logLik(fit_null))
  } else {
    ll_full <- as.numeric(stats::logLik(fit_full, REML = FALSE))
    ll_null <- as.numeric(stats::logLik(fit_null, REML = FALSE))
  }
  chisq <- 2 * (ll_full - ll_null)
  stats::pchisq(chisq, df = 1, lower.tail = FALSE)
}
