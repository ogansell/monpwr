#' monpwr: Power Analysis for Monitoring Programme Design
#'
#' @description
#' `monpwr` evaluates statistical power and precision for monitoring programme
#' design scenarios under three complementary modes:
#'
#' * **Conditional** — each plot is initialised at its observed visit number
#'   and BLUP from a fitted model; future visits are projected forward.
#'   The existing temporal record contributes to the trend estimate.
#'   Dropping a plot loses that record — a cost invisible to from-scratch
#'   analysis.
#'
#' * **Prospective** — from-scratch simulation using model parameters but
#'   no existing record.  Equivalent to `simr`-style power analysis.
#'   Useful for evaluating entirely new designs or comparing what power
#'   would have been without the existing investment.
#'
#' * **Hybrid** — a mix of legacy plots (initialised conditionally from
#'   their existing record) and new plots (initialised from the population
#'   marginal intercept with freshly sampled BLUPs each replicate).
#'   This is the natural mode when evaluating designs that expand an
#'   existing network with new sites.  Both the combined (legacy + new)
#'   and the legacy-only power are returned so the marginal contribution
#'   of new sites can be quantified explicitly.
#'
#' ## Key design principles
#'
#' * **Fully generalised** — no assumption about response variable, covariate
#'   structure, or grouping taxonomy.  Reporting scales are user-defined via
#'   the `reporting_groups` argument to [run_power_sim()].
#' * **Model-agnostic** — [extract_params()] dispatches on model class and
#'   returns a fixed-structure `monpwr_params` object.  All downstream
#'   functions consume only this object, not the original model.
#' * **Separable** — every public function can be called independently.
#'   The outer loop [run_power_sim()] is a convenience wrapper; advanced
#'   users can call [simulate_visits()], [build_historical()], and
#'   [fit_and_test()] directly.
#'
#' ## Typical workflow
#'
#' ```r
#' # 1. Fit your model (glmmTMB or lme4)
#' fit <- glmmTMB::glmmTMB(count ~ visit_num + (1|plot), ...)
#'
#' # 2. Extract standardised parameters — dispatches on class(fit)
#' ref <- extract_params(fit, data = long_model)
#'
#' # 3. Define design scenarios
#' scenarios <- list(
#'   baseline = scenario(
#'     label         = "Full grid, 5-yr",
#'     remeasure_yrs = 5,
#'     n_new_sites   = 0L,
#'     site_selector = function(sp) sp$site_id
#'   ),
#'   hybrid = scenario(
#'     label         = "Thinned grid + 30 new sites, 5-yr",
#'     remeasure_yrs = 5,
#'     n_new_sites   = 30L,
#'     site_selector = function(sp) sp$site_id[sp$in_coarse_grid]
#'   )
#' )
#'
#' # 4. Run simulation — reporting scales are fully user-defined
#' results <- run_power_sim(
#'   ref_params       = ref,
#'   scenarios        = scenarios,
#'   plot_metadata    = site_meta,
#'   data             = long_model,
#'   effect_sizes_pct = c(10, 20, 30),
#'   horizons         = c(10, 20),
#'   reporting_groups = list("Region" = "region_col",
#'                           "Subregion" = "subregion_col"),
#'   place_var        = "site_id"
#' )
#'
#' # 5. Summarise
#' mdc <- compute_mdc(results)
#' plot_power(results)
#' plot_mdc(mdc)
#' ```
#'
#' ## Model support
#'
#' Out of the box:
#' * `glmmTMB` — hurdle NB2, NB2, Poisson, binomial
#' * `lme4` — `glmerMod` (binomial, Poisson, negative binomial via
#'   `glmer.nb`), `lmerMod` (Gaussian, for log-transformed indices)
#'
#' For other families, supply a custom extractor; see
#' `vignette("custom-extractor")`.
#'
#' ## Key assumption
#'
#' `data` passed to [extract_params()] **must be a plain `data.frame`**,
#' not a tibble or `data.table`.  BLUP extraction via `ranef()` can
#' silently fail on non-data-frame inputs.  Call `as.data.frame(data)`
#' first if needed.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr filter mutate select group_by summarise ungroup arrange left_join full_join bind_rows distinct pull slice_max slice_sample n n_distinct count first if_else transmute rename pick everything all_of
#' @importFrom tidyr expand_grid pivot_wider
#' @importFrom purrr map map_dfr map_dbl pmap pmap_dfr walk
#' @importFrom furrr future_pmap_dfr furrr_options
#' @importFrom future plan multisession sequential
#' @importFrom parallelly availableCores
#' @importFrom rlang abort warn inform .data `%||%`
#' @importFrom cli cli_alert_info cli_alert_warning cli_alert_success cli_h2 cli_bullets cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom stats binom.test binomial family logLik median pchisq plogis poisson rbinom reorder rnbinom rnorm rpois sd setNames predict sigma
#' @importFrom tibble tibble
## usethis namespace: end
NULL
