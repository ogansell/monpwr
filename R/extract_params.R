#' Extract standardised parameters from a fitted model
#'
#' @description
#' Extracts all parameters needed to drive the `monpwr` simulation engine
#' from a fitted model object.  Dispatches on `class(fit)`, with methods
#' provided for `glmmTMB` and `lme4` model objects.
#'
#' The returned `monpwr_params` list has a **fixed structure regardless of
#' model family** — this is the contract that allows the simulation engine,
#' initialisers, and plotting functions to be completely model-agnostic.
#'
#' @param fit A fitted model object.  Supported classes:
#'   * `"glmmTMB"` — via [extract_params.glmmTMB()]
#'   * `"glmerMod"` or `"lmerMod"` — via [extract_params.glmerMod()]
#'   * Anything else — via [extract_params.default()], which stops with
#'     an informative error explaining how to write a custom extractor.
#' @param data A **plain `data.frame`** — the dataset used to fit `fit`,
#'   one row per plot × visit.  Must contain at minimum the columns named
#'   by `visit_num_var`, `plotid_var`, and `place_var`.
#' @param visit_num_var Character scalar.  Name of the visit sequence column.
#'   Default `"visit_num"`.
#' @param plotid_var Character scalar.  Name of the random-effect grouping
#'   column.  Default `"plotid_model"`.
#' @param place_var Character scalar.  Name of the site identifier column.
#'   Default `"Place"`.
#' @param visit_gap_var Character scalar.  Name of the visit-gap column used
#'   as a nuisance covariate.  Default `"visit_gap"`.  Set to `NULL` or a
#'   column not present in `data` to treat the gap effect as zero.
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A named list of class `"monpwr_params"` with elements:
#' \describe{
#'   \item{`beta_visit`}{Numeric scalar. Fixed-effect coefficient for
#'     `visit_num_var` in the conditional/count component — the trend
#'     parameter under test.}
#'   \item{`beta_gap_cond`}{Numeric scalar. Coefficient for the visit-gap
#'     variable in the conditional component; `0` if not in model.}
#'   \item{`beta_gap_zi`}{Numeric scalar. Coefficient for the visit-gap
#'     variable in the ZI component; `0` if not in model or model has no ZI.}
#'   \item{`disp_par`}{Numeric scalar. Dispersion parameter. NB2 `phi`;
#'     `1` for Poisson; residual SD for Gaussian; `1` for binomial.}
#'   \item{`sigma_cond`}{Numeric scalar. SD of the plot-level random
#'     intercept in the conditional/count component.}
#'   \item{`sigma_zi`}{Numeric scalar. SD of the plot-level random intercept
#'     in the ZI component; `0` if model has no ZI component.}
#'   \item{`visit_gap_med`}{Numeric scalar. Median visit-gap across `data`,
#'     used as a fixed nuisance value when simulating future visits.}
#'   \item{`marginal_int_cond`}{Numeric scalar. Population-level intercept
#'     for the conditional component at `visit_num = 0`.  Used to initialise
#'     new sites in hybrid scenarios via [init_new_sites()].}
#'   \item{`marginal_int_zi`}{Numeric scalar. Population-level intercept for
#'     the ZI component.  Used to initialise new sites in hybrid scenarios.}
#'   \item{`family`}{Character scalar. One of `"hurdle_nbinom2"`,
#'     `"nbinom2"`, `"poisson"`, `"binomial"`, `"gaussian"`.}
#'   \item{`visit_num_var`}{Character scalar. Passed through from argument.}
#'   \item{`plotid_var`}{Character scalar. Passed through from argument.}
#'   \item{`place_var`}{Character scalar. Passed through from argument.}
#'   \item{`visit_gap_var`}{Character scalar. Passed through from argument.}
#'   \item{`plot_state`}{Data frame, one row per unique plot. Columns:
#'     `Place`, `plotid_model`, `visit_num` (last observed), `eta_last_cond`,
#'     `eta_last_zi`, `blup_cond`, `blup_zi`. Used by [init_conditional()]
#'     and as a reference for prospective/hybrid initialisers.}
#' }
#'
#' @seealso [init_conditional()], [init_new_sites()],
#'   [init_prospective_marginal()], [run_power_sim()]
#'
#' @examples
#' \dontrun{
#' # glmmTMB hurdle model
#' fit <- glmmTMB::glmmTMB(
#'   count ~ visit_num + offset(log_effort) + (1 | plotid_model),
#'   ziformula = ~ (1 | plotid_model),
#'   family = glmmTMB::truncated_nbinom2,
#'   data = long_model
#' )
#' ref <- extract_params(fit, data = long_model)
#'
#' # lme4 binomial model
#' fit2 <- lme4::glmer(
#'   cbind(occupied, n_subplots - occupied) ~ visit_num + (1 | plotid_model),
#'   family = binomial, data = long_model
#' )
#' ref2 <- extract_params(fit2, data = long_model)
#' }
#'
#' @export
extract_params <- function(fit, data,
                           visit_num_var = "visit_num",
                           plotid_var    = "plotid_model",
                           place_var     = "Place",
                           visit_gap_var = "visit_gap",
                           ...) {
  UseMethod("extract_params")
}


# ------------------------------------------------------------------------------
# Shared validation helpers (not exported)
# ------------------------------------------------------------------------------

.validate_extract_inputs <- function(data, visit_num_var, plotid_var, place_var) {
  if (!is.data.frame(data) || inherits(data, "tbl_df")) {
    abort(c(
      "`data` must be a plain data.frame.",
      i = "Call `as.data.frame(data)` before passing to `extract_params()`.",
      i = "Tibbles and data.tables can cause silent failures in BLUP extraction."
    ))
  }
  required <- c(visit_num_var, plotid_var, place_var)
  missing  <- setdiff(required, names(data))
  if (length(missing) > 0) {
    abort(c(
      paste0("Required columns missing from `data`: ",
             paste(missing, collapse = ", ")),
      i = "Check `visit_num_var`, `plotid_var`, and `place_var` arguments."
    ))
  }
  invisible(NULL)
}


# Derive visit_gap_med and gap slopes from data and fixed-effect vectors.
# Returns list(beta_gap_cond, beta_gap_zi, visit_gap_med).
.extract_gap_params <- function(fe_cond, fe_zi = NULL,
                                data, visit_gap_var = "visit_gap") {
  has_gap <- !is.null(visit_gap_var) && visit_gap_var %in% names(data)
  list(
    beta_gap_cond = if (!is.null(fe_cond) && visit_gap_var %in% names(fe_cond))
      fe_cond[[visit_gap_var]] else 0,
    beta_gap_zi   = if (!is.null(fe_zi) && visit_gap_var %in% names(fe_zi))
      fe_zi[[visit_gap_var]] else 0,
    visit_gap_med = if (has_gap) median(data[[visit_gap_var]], na.rm = TRUE) else 0
  )
}


# Extract plot_state from fitted predictions + BLUPs.
.extract_plot_state <- function(fit, data, visit_num_var, plotid_var,
                                place_var, blups_cond, blups_zi,
                                predict_cond_fn, predict_zi_fn) {
  last_visit <- data |>
    group_by(.data[[place_var]], .data[[plotid_var]]) |>
    slice_max(.data[[visit_num_var]], n = 1, with_ties = FALSE) |>
    ungroup()

  eta_cond <- predict_cond_fn(last_visit)
  eta_zi   <- predict_zi_fn(last_visit)

  plot_state <- last_visit |>
    select(
      Place        = .data[[place_var]],
      plotid_model = .data[[plotid_var]],
      visit_num    = .data[[visit_num_var]]
    ) |>
    mutate(
      eta_last_cond = eta_cond,
      eta_last_zi   = eta_zi,
      blup_cond     = blups_cond[as.character(.data$plotid_model)],
      blup_zi       = blups_zi[as.character(.data$plotid_model)]
    )

  n_na <- sum(is.na(plot_state$blup_cond))
  if (n_na > 0) {
    abort(c(
      paste0("BLUP extraction failed for ", n_na, " plot(s)."),
      i = "Check that `plotid_var` values in `data` match the model's grouping factor levels.",
      i = "Ensure `data` is a plain data.frame (not a tibble)."
    ))
  }

  plot_state
}


# Derive marginal intercepts for new-site initialisation in hybrid scenarios.
# Recovers the population-average LP at visit_num = 0 by back-computing
# from plots at their earliest observed visit.  This is more robust than
# directly reading the "(Intercept)" fixed effect, which may be confounded
# with other baseline covariates in complex models.
.extract_marginal_intercepts <- function(plot_state, beta_visit) {
  ps <- plot_state
  if (nrow(ps) == 0) return(list(marginal_int_cond = 0, marginal_int_zi = 0))

  min_vis <- min(ps$visit_num, na.rm = TRUE)
  early   <- ps[ps$visit_num == min_vis, , drop = FALSE]

  list(
    marginal_int_cond = mean(early$eta_last_cond - beta_visit * (min_vis - 1L),
                             na.rm = TRUE),
    marginal_int_zi   = mean(early$eta_last_zi, na.rm = TRUE)
  )
}


# ------------------------------------------------------------------------------
# Method: glmmTMB
# ------------------------------------------------------------------------------

#' @rdname extract_params
#' @export
#'
#' @details
#' ## `glmmTMB` method
#'
#' Supports `truncated_nbinom2` (hurdle), `nbinom2`, `poisson`, and
#' `binomial` families.  The family string returned in `ref_params$family`
#' is derived from `family(fit)$family` and the presence of a non-trivial
#' ZI formula.
#'
#' `sigma_cond` and `sigma_zi` are extracted via [lme4::VarCorr()].  The
#' function stops with an informative error on VarCorr failure — it does
#' **not** fall back to positional TMB theta parameters, which would be
#' silently wrong.
#'
#' `marginal_int_cond` and `marginal_int_zi` are recovered by back-computing
#' from the earliest observed visit rather than reading the raw fixed-effect
#' intercept, which may be confounded with other baseline covariates.
#'
#' If `visit_num_var` is involved in an interaction term, `beta_visit` is the
#' main effect only.  Interaction contributions are absorbed into the per-plot
#' `eta_last_cond` via `predict(..., re.form = NULL)`.
extract_params.glmmTMB <- function(fit, data,
                                   visit_num_var = "visit_num",
                                   plotid_var    = "plotid_model",
                                   place_var     = "Place",
                                   visit_gap_var = "visit_gap",
                                   ...) {
  .validate_extract_inputs(data, visit_num_var, plotid_var, place_var)

  fe_cond <- glmmTMB::fixef(fit)$cond
  fe_zi   <- glmmTMB::fixef(fit)$zi

  if (!visit_num_var %in% names(fe_cond)) {
    abort(c(
      paste0("`", visit_num_var, "` not found in conditional fixed effects."),
      i = "Check `visit_num_var` argument or model formula."
    ))
  }
  beta_visit <- fe_cond[[visit_num_var]]

  gap <- .extract_gap_params(fe_cond, fe_zi, data, visit_gap_var)

  disp_par <- sigma(fit)

  vc <- tryCatch(
    lme4::VarCorr(fit),
    error = function(e) abort(c(
      "Failed to extract random-effect variances via `VarCorr()`.",
      i = "Inspect the model fit for convergence warnings.",
      x = conditionMessage(e)
    ))
  )

  sigma_cond <- tryCatch(
    sqrt(as.numeric(vc$cond[[plotid_var]])),
    error = function(e) abort(c(
      paste0("Could not find RE variance for `", plotid_var,
             "` in conditional component."),
      i = "Check that `plotid_var` matches the grouping factor name in the model."
    ))
  )

  sigma_zi <- tryCatch(
    sqrt(as.numeric(vc$zi[[plotid_var]])),
    error = function(e) 0
  )

  re_cond <- tryCatch(
    glmmTMB::ranef(fit)$cond[[plotid_var]],
    error = function(e) abort("Failed to extract conditional BLUPs via `ranef()`.")
  )
  re_zi <- tryCatch(glmmTMB::ranef(fit)$zi[[plotid_var]], error = function(e) NULL)

  blups_cond <- setNames(re_cond[, "(Intercept)"], rownames(re_cond))
  blups_zi   <- if (!is.null(re_zi)) {
    setNames(re_zi[, "(Intercept)"], rownames(re_zi))
  } else {
    setNames(rep(0, length(blups_cond)), names(blups_cond))
  }

  fam_raw <- family(fit)$family
  has_zi  <- !identical(fit$modelInfo$allForm$ziformula, ~0)
  fam_str <- if (grepl("truncated", fam_raw) || (grepl("nbinom2", fam_raw) && has_zi)) {
    "hurdle_nbinom2"
  } else if (grepl("nbinom2", fam_raw)) {
    "nbinom2"
  } else if (grepl("poisson", fam_raw, ignore.case = TRUE)) {
    "poisson"
  } else if (grepl("binomial", fam_raw, ignore.case = TRUE)) {
    "binomial"
  } else {
    fam_raw
  }

  plot_state <- .extract_plot_state(
    fit             = fit,
    data            = data,
    visit_num_var   = visit_num_var,
    plotid_var      = plotid_var,
    place_var       = place_var,
    blups_cond      = blups_cond,
    blups_zi        = blups_zi,
    predict_cond_fn = function(nd) predict(fit, newdata = nd,
                                                     type = "link", re.form = NULL),
    predict_zi_fn   = function(nd) tryCatch(
      predict(fit, newdata = nd, type = "zlink", re.form = NULL),
      error = function(e) rep(0, nrow(nd))
    )
  )

  marg <- .extract_marginal_intercepts(plot_state, beta_visit)

  structure(
    list(
      beta_visit      = beta_visit,
      beta_gap_cond   = gap$beta_gap_cond,
      beta_gap_zi     = gap$beta_gap_zi,
      disp_par        = disp_par,
      sigma_cond      = sigma_cond,
      sigma_zi        = sigma_zi,
      visit_gap_med   = gap$visit_gap_med,
      marginal_int_cond = marg$marginal_int_cond,
      marginal_int_zi   = marg$marginal_int_zi,
      family          = fam_str,
      visit_num_var   = visit_num_var,
      plotid_var      = plotid_var,
      place_var       = place_var,
      visit_gap_var   = visit_gap_var,
      plot_state      = plot_state
    ),
    class = "monpwr_params"
  )
}


# ------------------------------------------------------------------------------
# Method: lme4 (glmerMod / lmerMod)
# ------------------------------------------------------------------------------

#' @rdname extract_params
#' @export
#'
#' @details
#' ## `lme4` method (`glmerMod` / `lmerMod`)
#'
#' Handles `binomial`, `poisson`, and negative-binomial (`glmer.nb`) families.
#' For `lmerMod` (Gaussian), `disp_par` is set to `sigma(fit)` (residual SD)
#' and `family` is `"gaussian"`.
#'
#' There is no ZI component in `lme4` models; `beta_gap_zi`, `sigma_zi`,
#' `marginal_int_zi`, and `blup_zi` are all `0`.
extract_params.glmerMod <- function(fit, data,
                                    visit_num_var = "visit_num",
                                    plotid_var    = "plotid_model",
                                    place_var     = "Place",
                                    visit_gap_var = "visit_gap",
                                    ...) {
  .validate_extract_inputs(data, visit_num_var, plotid_var, place_var)

  fe <- lme4::fixef(fit)

  if (!visit_num_var %in% names(fe)) {
    abort(c(
      paste0("`", visit_num_var, "` not found in fixed effects."),
      i = "Check `visit_num_var` argument or model formula."
    ))
  }
  beta_visit <- fe[[visit_num_var]]
  gap        <- .extract_gap_params(fe, NULL, data, visit_gap_var)

  fam_obj  <- family(fit)
  fam_name <- fam_obj$family
  disp_par <- if (grepl("Negative Binomial", fam_name, ignore.case = TRUE)) {
    lme4::getME(fit, "glmer.nb.theta")
  } else 1

  fam_str <- if (grepl("Negative Binomial", fam_name, ignore.case = TRUE)) {
    "nbinom2"
  } else if (grepl("poisson", fam_name, ignore.case = TRUE)) {
    "poisson"
  } else if (grepl("binomial", fam_name, ignore.case = TRUE)) {
    "binomial"
  } else {
    fam_name
  }

  vc      <- lme4::VarCorr(fit)
  re_name <- grep(plotid_var, names(vc), value = TRUE)[1]
  if (is.na(re_name)) {
    abort(c(
      paste0("Could not find RE variance for `", plotid_var, "` in model."),
      i = "Check that `plotid_var` matches the grouping factor name in the model."
    ))
  }
  sigma_cond <- sqrt(as.numeric(vc[[re_name]]))

  re_df      <- lme4::ranef(fit)[[plotid_var]]
  blups_cond <- setNames(re_df[, "(Intercept)"], rownames(re_df))
  blups_zi   <- setNames(rep(0, length(blups_cond)), names(blups_cond))

  plot_state <- .extract_plot_state(
    fit             = fit,
    data            = data,
    visit_num_var   = visit_num_var,
    plotid_var      = plotid_var,
    place_var       = place_var,
    blups_cond      = blups_cond,
    blups_zi        = blups_zi,
    predict_cond_fn = function(nd) predict(fit, newdata = nd,
                                                  type = "link", re.form = NULL),
    predict_zi_fn   = function(nd) rep(0, nrow(nd))
  )

  marg <- .extract_marginal_intercepts(plot_state, beta_visit)

  structure(
    list(
      beta_visit        = beta_visit,
      beta_gap_cond     = gap$beta_gap_cond,
      beta_gap_zi       = 0,
      disp_par          = disp_par,
      sigma_cond        = sigma_cond,
      sigma_zi          = 0,
      visit_gap_med     = gap$visit_gap_med,
      marginal_int_cond = marg$marginal_int_cond,
      marginal_int_zi   = 0,
      family            = fam_str,
      visit_num_var     = visit_num_var,
      plotid_var        = plotid_var,
      place_var         = place_var,
      visit_gap_var     = visit_gap_var,
      plot_state        = plot_state
    ),
    class = "monpwr_params"
  )
}

#' @rdname extract_params
#' @export
extract_params.lmerMod <- function(fit, data,
                                   visit_num_var = "visit_num",
                                   plotid_var    = "plotid_model",
                                   place_var     = "Place",
                                   visit_gap_var = "visit_gap",
                                   ...) {
  .validate_extract_inputs(data, visit_num_var, plotid_var, place_var)

  fe <- lme4::fixef(fit)
  if (!visit_num_var %in% names(fe)) {
    abort(c(paste0("`", visit_num_var, "` not found in fixed effects.")))
  }

  beta_visit <- fe[[visit_num_var]]
  gap        <- .extract_gap_params(fe, NULL, data, visit_gap_var)
  disp_par   <- sigma(fit)

  vc      <- lme4::VarCorr(fit)
  re_name <- grep(plotid_var, names(vc), value = TRUE)[1]
  sigma_cond <- sqrt(as.numeric(vc[[re_name]]))

  re_df      <- lme4::ranef(fit)[[plotid_var]]
  blups_cond <- setNames(re_df[, "(Intercept)"], rownames(re_df))
  blups_zi   <- setNames(rep(0, length(blups_cond)), names(blups_cond))

  plot_state <- .extract_plot_state(
    fit             = fit,
    data            = data,
    visit_num_var   = visit_num_var,
    plotid_var      = plotid_var,
    place_var       = place_var,
    blups_cond      = blups_cond,
    blups_zi        = blups_zi,
    predict_cond_fn = function(nd) predict(fit, newdata = nd,
                                                  type = "link", re.form = NULL),
    predict_zi_fn   = function(nd) rep(0, nrow(nd))
  )

  marg <- .extract_marginal_intercepts(plot_state, beta_visit)

  structure(
    list(
      beta_visit        = beta_visit,
      beta_gap_cond     = gap$beta_gap_cond,
      beta_gap_zi       = 0,
      disp_par          = disp_par,
      sigma_cond        = sigma_cond,
      sigma_zi          = 0,
      visit_gap_med     = gap$visit_gap_med,
      marginal_int_cond = marg$marginal_int_cond,
      marginal_int_zi   = 0,
      family            = "gaussian",
      visit_num_var     = visit_num_var,
      plotid_var        = plotid_var,
      place_var         = place_var,
      visit_gap_var     = visit_gap_var,
      plot_state        = plot_state
    ),
    class = "monpwr_params"
  )
}


# ------------------------------------------------------------------------------
# Default method — informative error
# ------------------------------------------------------------------------------

#' @rdname extract_params
#' @export
extract_params.default <- function(fit, data,
                                   visit_num_var = "visit_num",
                                   plotid_var    = "plotid_model",
                                   place_var     = "Place",
                                   visit_gap_var = "visit_gap",
                                   ...) {
  abort(c(
    paste0("No `extract_params` method for class: `",
           paste(class(fit), collapse = "`, `"), "`."),
    i = "Supported classes: `glmmTMB`, `glmerMod`, `lmerMod`.",
    i = paste0(
      "To use a custom model, implement `extract_params.myclass()` returning ",
      "a list with the `monpwr_params` structure (including `marginal_int_cond` ",
      "and `marginal_int_zi`). See `vignette('custom-extractor')`."
    )
  ))
}


# ------------------------------------------------------------------------------
# print method
# ------------------------------------------------------------------------------

#' @export
print.monpwr_params <- function(x, ...) {
  cli::cli_h2("monpwr_params")
  cli::cli_bullets(c(
    "*" = paste0("Family:                ", x$family),
    "*" = paste0("visit_num slope:       ", round(x$beta_visit,        4)),
    "*" = paste0("RE SD (cond):          ", round(x$sigma_cond,        4)),
    "*" = paste0("RE SD (zi):            ", round(x$sigma_zi,          4)),
    "*" = paste0("Dispersion:            ", round(x$disp_par,          4)),
    "*" = paste0("Median visit gap:      ", round(x$visit_gap_med,     1)),
    "*" = paste0("Marginal int (cond):   ", round(x$marginal_int_cond, 4)),
    "*" = paste0("Marginal int (zi):     ", round(x$marginal_int_zi,   4)),
    "*" = paste0("Plots in state:        ", nrow(x$plot_state))
  ))
  invisible(x)
}
