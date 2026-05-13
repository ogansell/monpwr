#' Define a monitoring design scenario
#'
#' @description
#' Constructs a scenario object — a named list describing a single monitoring
#' design option.  Scenarios are passed as a named list to [run_power_sim()].
#'
#' Supports three design types depending on `n_new_sites`:
#' * **Legacy-only** (`n_new_sites = 0`) — only existing plots with an
#'   observed temporal record.  Uses conditional simulation throughout.
#' * **Prospective** — set `mode = "prospective"` in [run_power_sim()];
#'   `n_new_sites` is then the total number of plots (all from scratch).
#' * **Hybrid** (`n_new_sites > 0`, `mode = "conditional"`) — a mix of
#'   existing legacy plots and new plots with no prior record.  Legacy plots
#'   use conditional initialisation; new plots use [init_new_sites()].  Both
#'   a combined (legacy + new) and a legacy-only power trace are returned.
#'
#' @param label Character scalar.  Human-readable label for plots and tables.
#' @param remeasure_yrs Positive numeric scalar.  Remeasurement interval in
#'   years.  Determines how many future visits fall within each monitoring
#'   horizon: `n_future = max(1L, floor(horizon / remeasure_yrs))`.
#' @param site_selector A function with signature
#'   `function(plot_metadata) -> character vector` that returns the IDs of
#'   **legacy** sites to retain under this scenario.  Receives the
#'   `plot_metadata` data frame passed to [run_power_sim()].  May use any
#'   column in that frame.  For a scenario retaining all plots, use
#'   `function(sp) unique(sp[[place_var]])`.
#' @param n_new_sites Non-negative integer scalar.  Number of **new** sites
#'   to add (no prior record).  Default `0L` (legacy-only scenario).
#' @param eta_offset_cond Numeric scalar.  Shift on the marginal conditional
#'   intercept for new sites (log scale).  Default `0`.  Positive values mean
#'   new sites are expected to have higher baseline counts than the population
#'   average; negative for lower.  Passed to [init_new_sites()].
#' @param eta_offset_zi Numeric scalar.  Same shift for the ZI component.
#'   Defaults to `eta_offset_cond`.
#' @param new_site_init_fn Function or `NULL`.  Custom initialiser for new
#'   sites with signature `function(ref_params, n_new, ...)`.  If `NULL`,
#'   [init_new_sites()] is used.  Ignored when `n_new_sites = 0`.
#'
#' @return A list of class `"monpwr_scenario"`.
#'
#' @examples
#' \dontrun{
#' scenarios <- list(
#'   baseline = scenario(
#'     label         = "Full grid, 5-yr",
#'     remeasure_yrs = 5,
#'     site_selector = function(sp) unique(sp$site_id)
#'   ),
#'   hybrid = scenario(
#'     label         = "Coarse grid + 30 new sites, 5-yr",
#'     remeasure_yrs = 5,
#'     site_selector = function(sp) sp$site_id[sp$in_coarse_grid],
#'     n_new_sites   = 30L
#'   ),
#'   extended = scenario(
#'     label         = "Full grid, 10-yr",
#'     remeasure_yrs = 10,
#'     site_selector = function(sp) unique(sp$site_id)
#'   )
#' )
#' }
#'
#' @seealso [run_power_sim()], [init_new_sites()]
#' @export
scenario <- function(label,
                     remeasure_yrs,
                     site_selector,
                     n_new_sites       = 0L,
                     eta_offset_cond   = 0,
                     eta_offset_zi     = eta_offset_cond,
                     new_site_init_fn  = NULL) {

  if (!is.function(site_selector)) {
    abort("`site_selector` must be a function of the form function(plot_metadata) -> character vector.")
  }
  if (!is.numeric(remeasure_yrs) || length(remeasure_yrs) != 1 || remeasure_yrs <= 0) {
    abort("`remeasure_yrs` must be a single positive number.")
  }
  n_new_sites <- as.integer(n_new_sites)
  if (n_new_sites < 0L) {
    abort("`n_new_sites` must be a non-negative integer.")
  }
  if (!is.null(new_site_init_fn) && !is.function(new_site_init_fn)) {
    abort("`new_site_init_fn` must be a function or NULL.")
  }

  structure(
    list(
      label            = label,
      remeasure_yrs    = remeasure_yrs,
      site_selector    = site_selector,
      n_new_sites      = n_new_sites,
      eta_offset_cond  = eta_offset_cond,
      eta_offset_zi    = eta_offset_zi,
      new_site_init_fn = new_site_init_fn %||% init_new_sites
    ),
    class = "monpwr_scenario"
  )
}

#' @export
print.monpwr_scenario <- function(x, ...) {
  cli::cli_bullets(c(
    "*" = paste0("Label:          ", x$label),
    "*" = paste0("Remeasure (yr): ", x$remeasure_yrs),
    "*" = paste0("New sites:      ", x$n_new_sites),
    "*" = paste0("Hybrid:         ", if (x$n_new_sites > 0L) "yes" else "no")
  ))
  invisible(x)
}


# ==============================================================================

#' Run power simulations across design scenarios
#'
#' @description
#' The main outer loop of `monpwr`.  Iterates over scenarios, user-defined
#' reporting scales, effect sizes, and monitoring horizons, running `n_iter`
#' simulation replicates per cell in parallel via `furrr`.
#'
#' Supports three simulation modes:
#'
#' * **`"conditional"`** — each legacy plot is initialised at its observed
#'   visit number and BLUP (via [init_conditional()]); historical data are
#'   stacked with simulated future data so the existing record contributes
#'   to the trend estimate.  New sites (if `n_new_sites > 0` in the scenario)
#'   are initialised via [init_new_sites()] with freshly sampled BLUPs each
#'   replicate — this is the **hybrid** mode.
#'
#' * **`"prospective"`** — all plots initialised from scratch (via
#'   [init_prospective_marginal()] or a custom `init_fn`); only simulated
#'   future data are used.  Equivalent to `simr`-style analysis.
#'
#' Reporting scales are **fully user-defined** via `reporting_groups`, a
#' named list mapping display labels to column names in `plot_metadata`.
#' No grouping taxonomy (park type, region, land cover, etc.) is assumed.
#' An "Overall" scale covering all retained sites is always included.
#'
#' For hybrid scenarios (scenarios with `n_new_sites > 0`), both a combined
#' trace (legacy + new sites) and a legacy-only trace are returned in the
#' `sim_type` column, allowing the power gain from new sites to be
#' quantified explicitly.
#'
#' @param ref_params A `monpwr_params` object from [extract_params()].
#' @param scenarios A named list of `monpwr_scenario` objects from
#'   [scenario()].
#' @param plot_metadata A data frame with one row per unique **legacy** site.
#'   Must contain a column named `place_var` and any columns used by
#'   `site_selector` functions or listed in `reporting_groups`.
#' @param mode Character scalar.  `"conditional"` (default) or
#'   `"prospective"`.
#' @param data A data frame — the full modelling dataset (required for
#'   conditional/hybrid mode to build the historical stub).  Pass `NULL` for
#'   prospective mode.
#' @param effect_sizes_pct Numeric vector.  Effect sizes as % change in counts
#'   per visit period.  Converted internally to log-scale via
#'   `log(1 + pct / 100)`.
#' @param horizons Numeric vector.  Monitoring horizons in years.  Future
#'   visits per horizon: `max(1L, floor(horizon / remeasure_yrs))`.
#' @param n_iter Integer scalar.  Simulation replicates per cell.  Default
#'   200.
#' @param alpha Numeric scalar.  Significance threshold.  Default 0.10.
#' @param reporting_groups Named list mapping display labels to column names
#'   in `plot_metadata`.  Each unique level of each listed column becomes a
#'   separate reporting group.  Groups with fewer than `n_min_group` retained
#'   sites are skipped.  Default `list()` (overall scale only).
#' @param n_min_group Integer scalar.  Minimum retained sites for a
#'   group-level analysis.  Default 5.
#' @param place_var Character scalar.  Name of the site ID column in
#'   `plot_metadata` and `data`.  Must match `ref_params$place_var`.
#'   Default `"Place"`.
#' @param init_fn Function or `NULL`.  Custom initialiser for prospective
#'   mode, with signature `function(ref_params, n_plots, ...)`.  Ignored in
#'   conditional/hybrid mode.
#' @param workers Integer scalar.  Parallel workers for `furrr`.  Default
#'   `max(1, parallelly::availableCores() - 1)`.  Set to 1 to disable
#'   parallelism.
#' @param ... Additional arguments forwarded to `init_fn` (prospective) or
#'   `new_site_init_fn` (hybrid).
#'
#' @return A data frame of class `"monpwr_results"` with one row per
#'   scenario × sim_type × scale × group × effect size × horizon, and columns:
#'   `scenario`, `label`, `sim_type` (`"combined"` or `"legacy_only"`),
#'   `scale`, `group`, `effect_pct`, `horizon`, `n_legacy`, `n_new`,
#'   `n_total`, `n_future`, `power`, `n_converged`, `conv_rate`.
#'   `sim_type = "legacy_only"` rows only appear for hybrid scenarios.
#'
#' @seealso [scenario()], [extract_params()], [compute_mdc()],
#'   [run_precision()], [plot_power()], [plot_mdc()], [plot_cv()]
#' @export
run_power_sim <- function(ref_params,
                          scenarios,
                          plot_metadata,
                          mode             = c("conditional", "prospective"),
                          data             = NULL,
                          effect_sizes_pct = c(10, 20, 30),
                          horizons         = c(10, 20),
                          n_iter           = 200L,
                          alpha            = 0.10,
                          reporting_groups = list(),
                          n_min_group      = 5L,
                          place_var        = "Place",
                          init_fn          = NULL,
                          workers          = NULL,
                          ...) {

  mode <- match.arg(mode)

  # --- Validation ------------------------------------------------------------
  stopifnot(inherits(ref_params, "monpwr_params"))
  if (!is.list(scenarios) || length(scenarios) == 0) {
    abort("`scenarios` must be a non-empty named list of `monpwr_scenario` objects.")
  }
  if (mode == "conditional" && is.null(data)) {
    abort("`data` is required for `mode = 'conditional'`.")
  }
  if (!place_var %in% names(plot_metadata)) {
    abort(paste0("`place_var` column '", place_var,
                 "' not found in `plot_metadata`."))
  }
  if (length(reporting_groups) > 0) {
    missing_cols <- setdiff(unlist(reporting_groups), names(plot_metadata))
    if (length(missing_cols) > 0) {
      abort(c(
        paste0("Column(s) listed in `reporting_groups` not found in `plot_metadata`: ",
               paste(missing_cols, collapse = ", ")),
        i = "Check column names or remove the offending entry from `reporting_groups`."
      ))
    }
  }

  # --- Parallel setup --------------------------------------------------------
  n_workers <- workers %||% max(1L, parallelly::availableCores() - 1L)
  if (n_workers > 1L) {
    future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # --- Default prospective init_fn -------------------------------------------
  if (is.null(init_fn)) {
    init_fn <- init_prospective_marginal
  }

  # --- Outer scenario loop ---------------------------------------------------
  all_results <- vector("list", length(scenarios))

  for (sc_idx in seq_along(scenarios)) {
    sc_name  <- names(scenarios)[sc_idx]
    sc       <- scenarios[[sc_idx]]
    sites    <- sc$site_selector(plot_metadata)
    n_new    <- sc$n_new_sites
    is_hybrid <- (mode == "conditional") && (n_new > 0L)

    cli::cli_alert_info(
      "Scenario {sc_idx}/{length(scenarios)}: {sc$label} | legacy = {length(sites)} | new = {n_new}"
    )

    # Filter to legacy sites retained in this scenario
    plot_state_leg <- init_conditional(ref_params, sites)
    long_dat_leg   <- if (mode == "conditional") {
      build_historical(
        data,
        sites,
        place_var      = place_var,
        plotid_var     = ref_params$plotid_var,
        visit_num_var  = ref_params$visit_num_var,
        log_effort_var = if ("log_effort" %in% names(data)) "log_effort" else ref_params$visit_num_var,
        count_var      = "count"
      )
    } else NULL

    # Build reporting scale list
    scale_list <- .build_scale_list(
      sites            = sites,
      plot_metadata    = plot_metadata,
      ref_params       = ref_params,
      data             = data,
      mode             = mode,
      init_fn          = init_fn,
      reporting_groups = reporting_groups,
      n_min_group      = n_min_group,
      place_var        = place_var,
      ...
    )

    # Simulate power for each scale × effect × horizon
    sc_results <- map_dfr(scale_list, function(sl) {
      if (nrow(sl$plot_state) < 2) return(NULL)

      cli::cli_alert_info(
        "  Scale: {sl$scale} | Group: {sl$group} | n_legacy = {nrow(sl$plot_state)} | n_new = {n_new}"
      )

      expand_grid(
        effect_pct = effect_sizes_pct,
        horizon    = horizons
      ) |>
        future_pmap_dfr(function(effect_pct, horizon) {
          n_future <- max(1L, floor(horizon / sc$remeasure_yrs))

          # Combined trace (legacy + new sites)
          res_combined <- .run_one_cell(
            plot_state_leg   = sl$plot_state,
            hist_dat         = sl$hist_dat,
            n_new            = n_new,
            n_future         = n_future,
            effect_pct       = effect_pct,
            ref_params       = ref_params,
            n_iter           = n_iter,
            alpha            = alpha,
            mode             = mode,
            draw_re          = (mode == "prospective"),
            new_site_init_fn = sc$new_site_init_fn,
            eta_offset_cond  = sc$eta_offset_cond,
            eta_offset_zi    = sc$eta_offset_zi,
            ...
          ) |>
            mutate(
              sim_type = "combined",
              scale    = sl$scale,
              group    = sl$group,
              horizon  = horizon
            )

          # Legacy-only trace for hybrid scenarios — quantifies power gain
          res_legacy <- if (is_hybrid) {
            .run_one_cell(
              plot_state_leg   = sl$plot_state,
              hist_dat         = sl$hist_dat,
              n_new            = 0L,
              n_future         = n_future,
              effect_pct       = effect_pct,
              ref_params       = ref_params,
              n_iter           = n_iter,
              alpha            = alpha,
              mode             = mode,
              draw_re          = (mode == "prospective"),
              new_site_init_fn = sc$new_site_init_fn,
              eta_offset_cond  = 0,
              eta_offset_zi    = 0,
              ...
            ) |>
              mutate(
                sim_type = "legacy_only",
                scale    = sl$scale,
                group    = sl$group,
                horizon  = horizon
              )
          } else NULL

          bind_rows(res_combined, res_legacy)

        }, .options = furrr::furrr_options(seed = TRUE))
    }) |>
      mutate(scenario = sc_name, label = sc$label)

    all_results[[sc_idx]] <- sc_results
  }

  out <- bind_rows(all_results)
  class(out) <- c("monpwr_results", class(out))
  out
}


# ------------------------------------------------------------------------------
# Internal: build the list of reporting scales for one scenario
# ------------------------------------------------------------------------------

.build_scale_list <- function(sites, plot_metadata, ref_params, data, mode,
                              init_fn, reporting_groups, n_min_group,
                              place_var, ...) {

  make_entry <- function(scale, group, these_sites) {
    ps <- if (mode == "conditional") {
      init_conditional(ref_params, these_sites)
    } else {
      init_fn(ref_params, n_plots = length(these_sites), ...)
    }

    hd <- if (mode == "conditional" && !is.null(data)) {
      build_historical(
        data, these_sites,
        place_var      = place_var,
        plotid_var     = ref_params$plotid_var,
        visit_num_var  = ref_params$visit_num_var,
        log_effort_var = if ("log_effort" %in% names(data)) "log_effort"
                         else ref_params$visit_num_var,
        count_var      = "count"
      )
    } else NULL

    list(scale = scale, group = group, plot_state = ps, hist_dat = hd)
  }

  # Always include Overall scale
  sl <- list(make_entry("Overall", "All", sites))

  # User-defined reporting group scales
  for (grp_label in names(reporting_groups)) {
    grp_col <- reporting_groups[[grp_label]]
    if (!grp_col %in% names(plot_metadata)) next

    levels_in_sc <- plot_metadata |>
      filter(.data[[place_var]] %in% sites) |>
      count(.data[[grp_col]], name = "n") |>
      filter(.data$n >= n_min_group)

    for (i in seq_len(nrow(levels_in_sc))) {
      lvl        <- levels_in_sc[[grp_col]][i]
      grp_sites  <- plot_metadata[[place_var]][
        plot_metadata[[grp_col]] == lvl &
          plot_metadata[[place_var]] %in% sites
      ]
      if (length(grp_sites) >= 2) {
        sl <- c(sl, list(make_entry(grp_label, as.character(lvl), grp_sites)))
      }
    }
  }

  sl
}


# ------------------------------------------------------------------------------
# Internal: run one scenario × scale × effect × horizon cell
# New sites are re-initialised every replicate (BLUPs resampled).
# ------------------------------------------------------------------------------

.run_one_cell <- function(plot_state_leg, hist_dat,
                          n_new, n_future, effect_pct,
                          ref_params, n_iter, alpha,
                          mode, draw_re,
                          new_site_init_fn,
                          eta_offset_cond, eta_offset_zi,
                          ...) {
  eff_log  <- log(1 + effect_pct / 100)
  n_legacy <- nrow(plot_state_leg)

  p_vals <- map_dbl(seq_len(n_iter), function(i) {
    # Re-sample new-site BLUPs every replicate
    new_state <- if (n_new > 0L) {
      new_site_init_fn(
        ref_params,
        n_new           = n_new,
        eta_offset_cond = eta_offset_cond,
        eta_offset_zi   = eta_offset_zi,
        ...
      )
    } else NULL

    all_state <- bind_rows(plot_state_leg, new_state)

    future_dat <- simulate_visits(
      all_state, n_future, eff_log, ref_params,
      draw_re = draw_re
    )
    combined <- if (!is.null(hist_dat)) {
      bind_rows(hist_dat, future_dat)
    } else {
      future_dat
    }
    fit_and_test(combined, ref_params)
  })

  n_conv <- sum(!is.na(p_vals))

  tibble(
    n_legacy    = n_legacy,
    n_new       = n_new,
    n_total     = n_legacy + n_new,
    n_future    = n_future,
    effect_pct  = effect_pct,
    power       = mean(p_vals < alpha, na.rm = TRUE),
    n_converged = n_conv,
    conv_rate   = round(n_conv / n_iter, 3)
  )
}


# ------------------------------------------------------------------------------
# print method
# ------------------------------------------------------------------------------

#' @export
print.monpwr_results <- function(x, ...) {
  cli::cli_h2("monpwr_results")
  cli::cli_bullets(c(
    "*" = paste0("Scenarios:    ", paste(unique(x$label), collapse = ", ")),
    "*" = paste0("Sim types:    ", paste(unique(x$sim_type), collapse = ", ")),
    "*" = paste0("Scales:       ", paste(unique(x$scale), collapse = ", ")),
    "*" = paste0("Effect sizes: ", paste(sort(unique(x$effect_pct)), collapse = ", "), "%"),
    "*" = paste0("Horizons:     ", paste(sort(unique(x$horizon)), collapse = ", "), " yr"),
    "*" = paste0("Rows:         ", nrow(x))
  ))
  invisible(x)
}
