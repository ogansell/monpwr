#' Initialise legacy plot state for conditional simulation
#'
#' @description
#' Builds the plot-state data frame for **conditional** simulation by
#' filtering `ref_params$plot_state` to the retained site IDs.  Each plot
#' starts at its last observed visit number and fitted linear predictor,
#' so future simulated visits build on the existing monitoring record.
#'
#' This is the default initialiser used by [run_power_sim()] for legacy
#' plots when `mode = "conditional"` or in the legacy component of a hybrid
#' scenario.  You will not normally need to call it directly.
#'
#' @param ref_params A `monpwr_params` object from [extract_params()].
#' @param site_ids Character vector of site IDs (matching the `Place` column
#'   in `ref_params$plot_state`) for the retained plots under the current
#'   scenario.
#' @param ... Ignored; present for interface consistency.
#'
#' @return A data frame with one row per retained plot, containing:
#'   `Place`, `plotid_model`, `visit_num`, `eta_last_cond`, `eta_last_zi`,
#'   `blup_cond`, `blup_zi`.
#'
#' @seealso [init_new_sites()], [init_prospective_marginal()], [run_power_sim()]
#' @export
init_conditional <- function(ref_params, site_ids, ...) {
  stopifnot(inherits(ref_params, "monpwr_params"))

  state <- ref_params$plot_state |>
    filter(.data$place_id %in% site_ids)

  missing_ids <- setdiff(site_ids, state$place_id)
  if (length(missing_ids) > 0) {
    warn(c(
      paste0(length(missing_ids), " site ID(s) in `site_ids` not found in ",
             "`ref_params$plot_state` and will be dropped."),
      i = "Check that `site_ids` values match the `place_id` column in `plot_state`."
    ))
  }

  if (nrow(state) < 2) {
    abort("Fewer than 2 plots retained after filtering. Cannot run simulation.")
  }

  state
}


#' Initialise new sites for hybrid scenarios
#'
#' @description
#' Creates a synthetic plot-state data frame for **new sites** that have no
#' existing monitoring record.  Each new site is initialised at `visit_num = 0`
#' with a linear predictor drawn from the population marginal distribution:
#' the marginal intercept (at `visit_num = 0`) plus a freshly sampled BLUP
#' from \eqn{N(0, \hat\sigma)}.
#'
#' This function is called **inside every simulation replicate** by
#' [run_power_sim()] so that BLUPs are independently resampled each time.
#' Fixing BLUPs across replicates would understate new-site variability and
#' overstate power.
#'
#' @param ref_params A `monpwr_params` object from [extract_params()].
#' @param n_new Integer scalar.  Number of new sites to initialise.
#' @param eta_offset_cond Numeric scalar.  Optional shift on the marginal
#'   conditional intercept (log scale).  Use a positive value if new sites
#'   are targeted at areas expected to be above the population average, or
#'   negative for below-average areas.  Default `0`.
#' @param eta_offset_zi Numeric scalar.  Same shift for the ZI component.
#'   Defaults to `eta_offset_cond`.
#' @param id_prefix Character scalar.  Prefix for synthetic site ID strings.
#'   Default `".new_"`.  Should not clash with any real site IDs.
#' @param ... Ignored; present for interface consistency.
#'
#' @return A data frame with one row per new site, with the same column
#'   structure as `ref_params$plot_state` (`Place`, `plotid_model`,
#'   `visit_num`, `eta_last_cond`, `eta_last_zi`, `blup_cond`, `blup_zi`),
#'   suitable for binding with `init_conditional()` output before passing to
#'   [simulate_visits()].
#'
#' @details
#' ## When to use `eta_offset_cond`
#'
#' By default new sites are assumed to be drawn from the same population as
#' existing sites.  If your new sites are targeted (e.g. gap-filling in areas
#' known to have lower occupancy), supply a negative offset to shift the
#' starting intercept accordingly.  The offset is on the linear-predictor
#' (log or logit) scale.
#'
#' ## Writing a custom new-site initialiser
#'
#' If you need more control — for example, initialising new sites from a
#' specific covariate profile rather than the population marginal — write a
#' function with the signature:
#'
#' ```r
#' my_new_init <- function(ref_params, n_new, ...) {
#'   # return a data.frame with columns:
#'   #   Place, plotid_model, visit_num (0L), eta_last_cond, eta_last_zi,
#'   #   blup_cond, blup_zi
#' }
#' ```
#'
#' Pass it as `new_site_init_fn` to [scenario()].
#'
#' @seealso [init_conditional()], [scenario()], [run_power_sim()]
#' @export
init_new_sites <- function(ref_params, n_new,
                           eta_offset_cond = 0,
                           eta_offset_zi   = eta_offset_cond,
                           id_prefix       = ".new_",
                           ...) {
  stopifnot(inherits(ref_params, "monpwr_params"))
  n_new <- as.integer(n_new)
  if (n_new == 0L) return(NULL)
  if (n_new < 0L)  abort("`n_new` must be a non-negative integer.")

  blup_cond <- rnorm(n_new, mean = 0, sd = ref_params$sigma_cond)
  blup_zi   <- if (ref_params$sigma_zi > 0) {
    rnorm(n_new, mean = 0, sd = ref_params$sigma_zi)
  } else {
    rep(0, n_new)
  }

  ps <- ref_params$plot_state
  marginal_int_cond <- mean(ps$eta_last_cond - ps$blup_cond -
                              ref_params$beta_visit * ps$visit_num)
  marginal_int_zi   <- if (ref_params$sigma_zi > 0) {
    mean(ps$eta_last_zi - ps$blup_zi)
  } else {
    0
  }

  ids <- paste0(id_prefix, seq_len(n_new))

  data.frame(
    place_id  = ids,
    plotid    = ids,
    visit_num = 0L,
    eta_last_cond = marginal_int_cond + eta_offset_cond + blup_cond,
    eta_last_zi   = marginal_int_zi + eta_offset_zi + blup_zi,
    blup_cond     = blup_cond,
    blup_zi       = blup_zi,
    stringsAsFactors = FALSE
  )
}


#' Initialise plot state for prospective (from-scratch) simulation
#'
#' @description
#' Builds a synthetic plot-state data frame for **prospective** simulation.
#' Each plot is initialised at `visit_num = 0` with linear predictors set to
#' the **marginal (population-average) prediction**.  Random effects are drawn
#' fresh each simulation replicate inside [simulate_visits()].
#'
#' This represents "what power would we have if we were starting today with
#' no existing monitoring record" — equivalent to `simr`-style power analysis.
#'
#' This is the default initialiser used by [run_power_sim()] when
#' `mode = "prospective"`.
#'
#' @param ref_params A `monpwr_params` object from [extract_params()].
#' @param n_plots Integer scalar.  Number of plots to simulate.
#' @param ... Ignored; present for interface consistency with custom
#'   `init_fn` functions.
#'
#' @details
#' ## Writing a custom `init_fn` for prospective mode
#'
#' If you want to initialise prospective plots with specific covariate
#' profiles (e.g. sampled from your empirical data, or set to a design
#' target), supply your own function with the signature:
#'
#' ```r
#' my_init_fn <- function(ref_params, n_plots, ...) {
#'   # return a data.frame with columns:
#'   #   plotid        <character>  unique plot identifier
#'   #   visit_num     <integer>    starting visit number (typically 0L)
#'   #   eta_last_cond <numeric>    starting linear predictor, conditional
#'   #   eta_last_zi   <numeric>    starting linear predictor, ZI (0 if no ZI)
#' }
#'
#' run_power_sim(..., mode = "prospective", init_fn = my_init_fn)
#' ```
#'
#' @return A data frame with one row per plot, containing:
#'   `plotid` (character), `visit_num` (integer, `0L`),
#'   `eta_last_cond` (numeric, marginal intercept),
#'   `eta_last_zi` (numeric, marginal ZI intercept or `0`).
#'
#' @seealso [init_conditional()], [init_new_sites()], [run_power_sim()]
#' @export
init_prospective_marginal <- function(ref_params, n_plots, ...) {
  stopifnot(inherits(ref_params, "monpwr_params"))
  n_plots <- as.integer(n_plots)
  if (n_plots < 2L) abort("`n_plots` must be >= 2.")

  ps <- ref_params$plot_state
  marginal_int_cond <- mean(ps$eta_last_cond - ps$blup_cond -
                              ref_params$beta_visit * ps$visit_num)
  marginal_int_zi   <- if (ref_params$sigma_zi > 0) {
    mean(ps$eta_last_zi - ps$blup_zi)
  } else {
    0
  }

  data.frame(
    plotid        = paste0("plot_", seq_len(n_plots)),
    visit_num     = 0L,
    eta_last_cond = marginal_int_cond,
    eta_last_zi   = marginal_int_zi,
    stringsAsFactors = FALSE
  )
}
