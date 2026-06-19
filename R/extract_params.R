#' Extract standardised parameters from a fitted model
#'
#' @description
#' Extracts all parameters needed to drive the `monpwr` simulation engine
#' from a fitted model object.  Dispatches on `class(fit)`, with methods
#' provided for `glmmTMB` and `lme4` model objects.
#'
#' The returned `ref_params` list has a **fixed structure regardless of
#' model family** — this is the contract that allows the simulation engine,
#' initialisers, and plotting functions to be completely model-agnostic.
#'
#' @param fit A fitted model object.  Supported classes:
#'   * `"glmmTMB"` — via [extract_params.glmmTMB()]
#'   * `"glmerMod"` or `"lmerMod"` — via [extract_params.glmerMod()]
#'   * Anything else — via [extract_params.default()], which stops with
#'     an informative error explaining how to write a custom extractor.
#' @param data A **plain `data.frame`** — the dataset used to fit `fit`,
#'   one row per plot × visit.  Must contain at minimum:
#'   * `Place` — plot identifier (character or factor)
#'   * `plotid_model` — plot ID as used in the model random effect (character
#'     or factor; often the same as `Place`)
#'   * `visit_num` — visit sequence number (integer, 1 = first visit)
#'   * `visit_gap` — years since previous visit (numeric; used for nuisance
#'     fixing in simulation; can be `NA` or a constant if not in model)
#' @param visit_num_var Character scalar.  Name of the column in `data`
#'   representing the visit sequence number.  Default `"visit_num"`.
#' @param plotid_var Character scalar.  Name of the random-effect grouping
#'   column in `data`.  Default `"plotid_model"`.
#' @param place_var Character scalar.  Name of the plot identifier column
#'   in `data`.  Default `"Place"`.
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A named list of class `"monpwr_params"` with elements:
#' \describe{
#'   \item{`beta_visit`}{Numeric scalar. Fixed-effect coefficient for the
#'     visit sequence variable in the conditional/count component. This is
#'     the trend parameter under test.}
#'   \item{`beta_gap_cond`}{Numeric scalar. Coefficient for `visit_gap` in
#'     the conditional component; `0` if not in model.}
#'   \item{`beta_gap_zi`}{Numeric scalar. Coefficient for `visit_gap` in
#'     the ZI component; `0` if not in model or model has no ZI.}
#'   \item{`disp_par`}{Numeric scalar. Dispersion parameter. NB2 `phi`;
#'     `1` for Poisson; `Inf` for Gaussian (unused); binomial `1`.}
#'   \item{`sigma_cond`}{Numeric scalar. SD of the plot-level random
#'     intercept in the conditional/count component.}
#'   \item{`sigma_zi`}{Numeric scalar. SD of the plot-level random intercept
#'     in the ZI component; `0` if model has no ZI component.}
#'   \item{`visit_gap_med`}{Numeric scalar. Median `visit_gap` across `data`,
#'     used as a fixed nuisance value when simulating future visits.}
#'   \item{`family`}{Character scalar. One of `"hurdle_nbinom2"`,
#'     `"nbinom2"`, `"poisson"`, `"binomial"`, `"gaussian"`.}
#'   \item{`visit_num_var`}{Character scalar. Passed through from argument.}
#'   \item{`plotid_var`}{Character scalar. Passed through from argument.}
#'   \item{`place_var`}{Character scalar. Passed through from argument.}
#'   \item{`plot_state`}{Data frame, one row per unique plot. Columns:
#'     `Place`, `plotid_model`, `visit_num` (last observed), `eta_last_cond`,
#'     `eta_last_zi`, `blup_cond`, `blup_zi`. Used by [init_conditional()]
#'     and as a template for prospective initialisers.}
#' }
#'
#' @seealso [init_conditional()], [init_prospective_marginal()],
#'   [run_power_sim()]
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
                           visit_num_var     = "visit_num",
                           plotid_var        = "plotid_model",
                           place_var         = "Place",
                           offset_var        = NULL,
                           offset_transform  = NULL,
                           log_effort_future = NULL,
                           ...) {
  UseMethod("extract_params")
}


# ------------------------------------------------------------------------------
# Shared validation helpers (not exported)
# ------------------------------------------------------------------------------

# Extract the response variable name from a model formula.
# Returns a simple identifier if the LHS is one; falls back to "count".
.response_var <- function(fit) {
  lhs <- tryCatch(deparse(formula(fit)[[2]]), error = function(e) "count")
  if (grepl("^[a-zA-Z.][a-zA-Z0-9._]*$", lhs)) lhs else "count"
}

.validate_extract_inputs <- function(data, visit_num_var, plotid_var, place_var) {
  if (!identical(class(data), "data.frame")) {
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

# Materialise a log-effort vector from data given the three offset strategies.
# Returns a numeric vector of length nrow(data) and a scalar log_effort_future.
.resolve_offset <- function(data, offset_var, offset_transform, log_effort_future) {
  log_eff <- if (!is.null(offset_var)) {
    if (!offset_var %in% names(data)) {
      abort(c(
        paste0("`offset_var` column '", offset_var, "' not found in `data`."),
        i = "Supply the correct column name or use `offset_transform` for inline offsets."
      ))
    }
    data[[offset_var]]
  } else if (!is.null(offset_transform)) {
    if (!is.function(offset_transform)) {
      abort("`offset_transform` must be a function of the form function(data) -> numeric vector.")
    }
    result <- offset_transform(data)
    if (!is.numeric(result) || length(result) != nrow(data)) {
      abort("`offset_transform` must return a numeric vector with one value per row of `data`.")
    }
    result
  } else {
    rep(0, nrow(data))
  }

  future_val <- log_effort_future %||% median(log_eff, na.rm = TRUE)

  list(log_eff = log_eff, log_effort_future = future_val)
}


.extract_plot_state <- function(fit, data, visit_num_var, plotid_var,
                                place_var, blups_cond, blups_zi,
                                predict_cond_fn, predict_zi_fn) {
  # Last observed visit per plot
  # unique() handles the case where place_var == plotid_var
  group_cols <- unique(c(place_var, plotid_var))
  last_visit <- data |>
    group_by(across(all_of(group_cols))) |>
    slice_max(.data[[visit_num_var]], n = 1, with_ties = FALSE) |>
    ungroup()

  eta_cond <- predict_cond_fn(last_visit)
  eta_zi   <- predict_zi_fn(last_visit)

  # transmute avoids duplicate-column issues when place_var == plotid_var
  plot_state <- last_visit |>
    transmute(
      place_id      = as.character(.data[[place_var]]),
      plotid        = as.character(.data[[plotid_var]]),
      visit_num     = .data[[visit_num_var]],
      eta_last_cond = eta_cond,
      eta_last_zi   = eta_zi,
      blup_cond     = blups_cond[as.character(.data[[plotid_var]])],
      blup_zi       = blups_zi[as.character(.data[[plotid_var]])]
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
#' `sigma_cond` and `sigma_zi` are extracted via [lme4::VarCorr()].  If
#' extraction fails (e.g. non-positive-definite Hessian), the function stops
#' with an error — it does **not** fall back to positional TMB theta
#' parameters, which would be silently wrong.
#'
#' If `visit_num` is involved in an interaction term, `beta_visit` is the
#' main effect only.  Interaction contributions are absorbed into the per-plot
#' `eta_last_cond` via `predict(..., re.form = NULL)`, so the starting
#' linear predictor is correct for each plot's observed covariates.
extract_params.glmmTMB <- function(fit, data,
                                   visit_num_var     = "visit_num",
                                   plotid_var        = "plotid_model",
                                   place_var         = "Place",
                                   visit_gap_var     = "visit_gap",
                                   offset_var        = NULL,
                                   offset_transform  = NULL,
                                   log_effort_future = NULL,
                                   ...) {
  .validate_extract_inputs(data, visit_num_var, plotid_var, place_var)
  off <- .resolve_offset(data, offset_var, offset_transform, log_effort_future)

  fe_cond <- glmmTMB::fixef(fit)$cond
  fe_zi   <- glmmTMB::fixef(fit)$zi

  # visit_num slope — main effect only
  if (!visit_num_var %in% names(fe_cond)) {
    abort(c(
      paste0("`", visit_num_var, "` not found in conditional fixed effects."),
      i = "Check `visit_num_var` argument or model formula."
    ))
  }
  beta_visit <- fe_cond[[visit_num_var]]

  # visit_gap slopes
  beta_gap_cond <- if (visit_gap_var %in% names(fe_cond)) fe_cond[[visit_gap_var]] else 0
  beta_gap_zi   <- if (visit_gap_var %in% names(fe_zi))   fe_zi[[visit_gap_var]]   else 0

  # Dispersion
  disp_par <- glmmTMB::sigma(fit)

  # RE SDs — hard fail if VarCorr doesn't work
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
    error = function(e) 0  # no ZI RE is fine
  )

  # BLUPs
  re_cond <- tryCatch(
    glmmTMB::ranef(fit)$cond[[plotid_var]],
    error = function(e) abort("Failed to extract conditional BLUPs via `ranef()`.")
  )
  re_zi <- tryCatch(
    glmmTMB::ranef(fit)$zi[[plotid_var]],
    error = function(e) NULL
  )

  blups_cond <- setNames(re_cond[, "(Intercept)"], rownames(re_cond))
  blups_zi   <- if (!is.null(re_zi)) {
    setNames(re_zi[, "(Intercept)"], rownames(re_zi))
  } else {
    setNames(rep(0, length(blups_cond)), names(blups_cond))
  }

  # visit_gap median
  visit_gap_med <- if (visit_gap_var %in% names(data)) {
    median(data[[visit_gap_var]], na.rm = TRUE)
  } else 0

  count_var <- .response_var(fit)

  # Family string
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
    fam_raw  # pass through unknown families
  }

  # Per-plot state
  plot_state <- .extract_plot_state(
    fit          = fit,
    data         = data,
    visit_num_var = visit_num_var,
    plotid_var   = plotid_var,
    place_var    = place_var,
    blups_cond   = blups_cond,
    blups_zi     = blups_zi,
    predict_cond_fn = function(nd) predict(fit, newdata = nd,
                                           type = "link",  re.form = NULL),
    predict_zi_fn   = function(nd) tryCatch(
      predict(fit, newdata = nd, type = "zlink", re.form = NULL),
      error = function(e) rep(0, nrow(nd))
    )
  )

  structure(
    list(
      beta_visit        = beta_visit,
      beta_gap_cond     = beta_gap_cond,
      beta_gap_zi       = beta_gap_zi,
      disp_par          = disp_par,
      sigma_cond        = sigma_cond,
      sigma_zi          = sigma_zi,
      visit_gap_med     = visit_gap_med,
      family            = fam_str,
      visit_num_var     = visit_num_var,
      plotid_var        = plotid_var,
      place_var         = place_var,
      visit_gap_var     = visit_gap_var,
      count_var         = count_var,
      offset_var        = offset_var,
      log_effort_future = off$log_effort_future,
      plot_state        = plot_state
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
#' and `blup_zi` are all `0`.
extract_params.glmerMod <- function(fit, data,
                                    visit_num_var     = "visit_num",
                                    plotid_var        = "plotid_model",
                                    place_var         = "Place",
                                    visit_gap_var     = "visit_gap",
                                    offset_var        = NULL,
                                    offset_transform  = NULL,
                                    log_effort_future = NULL,
                                    ...) {
  .validate_extract_inputs(data, visit_num_var, plotid_var, place_var)
  off <- .resolve_offset(data, offset_var, offset_transform, log_effort_future)

  fe <- lme4::fixef(fit)

  if (!visit_num_var %in% names(fe)) {
    abort(c(
      paste0("`", visit_num_var, "` not found in fixed effects."),
      i = "Check `visit_num_var` argument or model formula."
    ))
  }
  beta_visit    <- fe[[visit_num_var]]
  beta_gap_cond <- if (visit_gap_var %in% names(fe)) fe[[visit_gap_var]] else 0

  # Dispersion / family
  fam_obj  <- family(fit)
  fam_name <- fam_obj$family
  disp_par <- if (grepl("Negative Binomial", fam_name, ignore.case = TRUE)) {
    lme4::getME(fit, "glmer.nb.theta")
  } else {
    1  # Poisson / binomial
  }

  fam_str <- if (grepl("Negative Binomial", fam_name, ignore.case = TRUE)) {
    "nbinom2"
  } else if (grepl("poisson", fam_name, ignore.case = TRUE)) {
    "poisson"
  } else if (grepl("binomial", fam_name, ignore.case = TRUE)) {
    "binomial"
  } else {
    fam_name
  }

  # RE SD
  vc <- lme4::VarCorr(fit)
  re_name <- grep(plotid_var, names(vc), value = TRUE)[1]
  if (is.na(re_name)) {
    abort(c(
      paste0("Could not find RE variance for `", plotid_var, "` in model."),
      i = "Check that `plotid_var` matches the grouping factor name in the model."
    ))
  }
  sigma_cond <- sqrt(as.numeric(vc[[re_name]]))

  # BLUPs
  re_df <- lme4::ranef(fit)[[plotid_var]]
  blups_cond <- setNames(re_df[, "(Intercept)"], rownames(re_df))
  blups_zi   <- setNames(rep(0, length(blups_cond)), names(blups_cond))

  visit_gap_med <- if (visit_gap_var %in% names(data)) {
    median(data[[visit_gap_var]], na.rm = TRUE)
  } else 0

  count_var <- .response_var(fit)

  plot_state <- .extract_plot_state(
    fit           = fit,
    data          = data,
    visit_num_var = visit_num_var,
    plotid_var    = plotid_var,
    place_var     = place_var,
    blups_cond    = blups_cond,
    blups_zi      = blups_zi,
    predict_cond_fn = function(nd) lme4::predict(fit, newdata = nd,
                                                  type = "link", re.form = NULL),
    predict_zi_fn   = function(nd) rep(0, nrow(nd))
  )

  structure(
    list(
      beta_visit        = beta_visit,
      beta_gap_cond     = beta_gap_cond,
      beta_gap_zi       = 0,
      disp_par          = disp_par,
      sigma_cond        = sigma_cond,
      sigma_zi          = 0,
      visit_gap_med     = visit_gap_med,
      family            = fam_str,
      visit_num_var     = visit_num_var,
      plotid_var        = plotid_var,
      place_var         = place_var,
      visit_gap_var     = visit_gap_var,
      count_var         = count_var,
      offset_var        = offset_var,
      log_effort_future = off$log_effort_future,
      plot_state        = plot_state
    ),
    class = "monpwr_params"
  )
}

#' @rdname extract_params
#' @export
extract_params.lmerMod <- function(fit, data,
                                   visit_num_var     = "visit_num",
                                   plotid_var        = "plotid_model",
                                   place_var         = "Place",
                                   visit_gap_var     = "visit_gap",
                                   offset_var        = NULL,
                                   offset_transform  = NULL,
                                   log_effort_future = NULL,
                                   ...) {
  .validate_extract_inputs(data, visit_num_var, plotid_var, place_var)
  off <- .resolve_offset(data, offset_var, offset_transform, log_effort_future)

  fe <- lme4::fixef(fit)
  if (!visit_num_var %in% names(fe)) {
    abort(c(paste0("`", visit_num_var, "` not found in fixed effects.")))
  }

  beta_visit    <- fe[[visit_num_var]]
  beta_gap_cond <- if (visit_gap_var %in% names(fe)) fe[[visit_gap_var]] else 0
  disp_par      <- lme4::sigma(fit)  # residual SD for Gaussian

  vc      <- lme4::VarCorr(fit)
  re_name <- grep(plotid_var, names(vc), value = TRUE)[1]
  sigma_cond <- sqrt(as.numeric(vc[[re_name]]))

  re_df      <- lme4::ranef(fit)[[plotid_var]]
  blups_cond <- setNames(re_df[, "(Intercept)"], rownames(re_df))
  blups_zi   <- setNames(rep(0, length(blups_cond)), names(blups_cond))

  visit_gap_med <- if (visit_gap_var %in% names(data)) {
    median(data[[visit_gap_var]], na.rm = TRUE)
  } else 0

  count_var <- .response_var(fit)

  plot_state <- .extract_plot_state(
    fit           = fit,
    data          = data,
    visit_num_var = visit_num_var,
    plotid_var    = plotid_var,
    place_var     = place_var,
    blups_cond    = blups_cond,
    blups_zi      = blups_zi,
    predict_cond_fn = function(nd) lme4::predict(fit, newdata = nd,
                                                  type = "link", re.form = NULL),
    predict_zi_fn   = function(nd) rep(0, nrow(nd))
  )

  structure(
    list(
      beta_visit        = beta_visit,
      beta_gap_cond     = beta_gap_cond,
      beta_gap_zi       = 0,
      disp_par          = disp_par,
      sigma_cond        = sigma_cond,
      sigma_zi          = 0,
      visit_gap_med     = visit_gap_med,
      family            = "gaussian",
      visit_num_var     = visit_num_var,
      plotid_var        = plotid_var,
      place_var         = place_var,
      visit_gap_var     = visit_gap_var,
      count_var         = count_var,
      offset_var        = offset_var,
      log_effort_future = off$log_effort_future,
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
                                   ...) {
  abort(c(
    paste0("No `extract_params` method for class: `",
           paste(class(fit), collapse = "`, `"), "`."),
    i = "Supported classes: `glmmTMB`, `glmerMod`, `lmerMod`.",
    i = paste0(
      "To use a custom model, implement `extract_params.myclass()` that ",
      "returns a list with the `monpwr_params` structure. ",
      "See `vignette('custom-extractor')` for the required fields and a worked example."
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
    "*" = paste0("Family:            ", x$family),
    "*" = paste0("visit_num slope:   ", round(x$beta_visit, 4)),
    "*" = paste0("RE SD (cond):      ", round(x$sigma_cond, 4)),
    "*" = paste0("RE SD (zi):        ", round(x$sigma_zi,   4)),
    "*" = paste0("Dispersion:        ", round(x$disp_par,   4)),
    "*" = paste0("Median visit gap:  ", round(x$visit_gap_med, 1)),
    "*" = paste0("Offset var:        ", x$offset_var %||% "(none)"),
    "*" = paste0("log_effort_future: ", round(x$log_effort_future, 4)),
    "*" = paste0("Plots in state:    ", nrow(x$plot_state))
  ))
  invisible(x)
}
