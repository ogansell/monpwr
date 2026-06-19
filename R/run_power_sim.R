#' Define a monitoring design scenario
#'
#' @description
#' Constructs a scenario object — a named list describing a single monitoring
#' design option.  Scenarios are passed as a named list to [run_power_sim()].
#'
#' @param label Character scalar.  Human-readable label used in plots and
#'   tables (e.g. `"8x8km grid, 5-yr interval"`).
#' @param remeasure_yrs Numeric scalar.  Remeasurement interval in years.
#'   Determines how many future visits fall within each monitoring horizon:
#'   `n_future = floor(horizon / remeasure_yrs)`.
#' @param site_selector A function with signature `function(plot_metadata)`
#'   that takes the `plot_metadata` data frame passed to [run_power_sim()]
#'   and returns a character vector of `Place` IDs to retain under this
#'   scenario.  For a scenario retaining all plots, use
#'   `function(sp) unique(sp$Place)`.
#'
#' @return A list of class `"monpwr_scenario"` with elements `label`,
#'   `remeasure_yrs`, and `site_selector`.
#'
#' @examples
#' \dontrun{
#' scenarios <- list(
#'   baseline = scenario(
#'     label         = "Full grid, 5-yr",
#'     remeasure_yrs = 5,
#'     site_selector = function(sp) unique(sp$Place)
#'   ),
#'   coarse = scenario(
#'     label         = "Coarse grid, 5-yr",
#'     remeasure_yrs = 5,
#'     site_selector = function(sp) sp$Place[sp$in_coarse_grid]
#'   ),
#'   extended = scenario(
#'     label         = "Full grid, 10-yr",
#'     remeasure_yrs = 10,
#'     site_selector = function(sp) unique(sp$Place)
#'   )
#' )
#' }
#'
#' @seealso [run_power_sim()]
#' @export
scenario <- function(label, remeasure_yrs, site_selector) {
  if (!is.function(site_selector)) {
    abort("`site_selector` must be a function of the form function(plot_metadata) -> character vector.")
  }
  if (!is.numeric(remeasure_yrs) || remeasure_yrs <= 0) {
    abort("`remeasure_yrs` must be a positive number.")
  }
  structure(
    list(label = label, remeasure_yrs = remeasure_yrs, site_selector = site_selector),
    class = "monpwr_scenario"
  )
}

#' @export
print.monpwr_scenario <- function(x, ...) {
  cli::cli_bullets(c(
    "*" = paste0("Label:         ", x$label),
    "*" = paste0("Remeasure (yr):", x$remeasure_yrs)
  ))
  invisible(x)
}


# ------------------------------------------------------------------------------

#' Run power simulations across design scenarios
#'
#' @description
#' The main outer loop of `monpwr`.  Iterates over scenarios, reporting
#' scales, effect sizes, and monitoring horizons, running `n_iter` simulation
#' replicates per cell in parallel via `furrr`.
#'
#' Supports two simulation modes:
#'
#' * **`"conditional"`** — each plot is initialised at its observed visit
#'   number and BLUP (via [init_conditional()]); the historical data stub is
#'   stacked with simulated future data so the existing record contributes
#'   to the trend estimate.
#'
#' * **`"prospective"`** — plots are initialised from scratch (via
#'   [init_prospective_marginal()] or a custom `init_fn`); only simulated
#'   future data are used in the trend test.  Equivalent to `simr`-style
#'   analysis.
#'
#' @param ref_params A `monpwr_params` object from [extract_params()].
#' @param scenarios A named list of `monpwr_scenario` objects from
#'   [scenario()].
#' @param plot_metadata A data frame with one row per unique plot, containing
#'   at minimum a `Place` column (or the value of `place_var`) and any
#'   columns used by `site_selector` functions (e.g. park type, spatial flags).
#' @param mode Character scalar.  `"conditional"` (default) or
#'   `"prospective"`.  See Details.
#' @param data A data frame — the full modelling dataset (required for
#'   conditional mode; used to build the historical stub).  Not needed for
#'   prospective mode; pass `NULL`.
#' @param effect_sizes_pct Numeric vector.  Effect sizes as percentage change
#'   in counts per visit period (e.g. `c(10, 20, 30)`).  Converted
#'   internally to log-scale coefficients via `log(1 + pct / 100)`.
#' @param horizons Numeric vector.  Monitoring horizons in years
#'   (e.g. `c(10, 20)`).  The number of future visits per horizon is
#'   `max(1L, floor(horizon / remeasure_yrs))`.
#' @param n_iter Integer scalar.  Number of simulation replicates per cell.
#'   Default 200.  Use a smaller value (e.g. 30) for quick smoke-tests.
#' @param alpha Numeric scalar.  Significance threshold.  Default 0.10
#'   (consistent with standard Tier 1 practice).
#' @param scales Character vector.  Reporting scales to include.  Any
#'   combination of `"National"`, `"Park type"`, `"Park"`.
#'   Default `c("National", "Park type")`.
#' @param park_type_var Character scalar.  Column in `plot_metadata`
#'   distinguishing park types (e.g. `"Park_type"`).  Required if
#'   `"Park type"` or `"Park"` is in `scales`.
#' @param park_name_var Character scalar.  Column in `plot_metadata` for
#'   individual park names (e.g. `"PCL_Name"`).  Required if `"Park"` is
#'   in `scales`.
#' @param n_min_park Integer scalar.  Minimum number of retained plots for
#'   a park-level analysis to be run.  Default 5.
#' @param place_var Character scalar.  Name of the plot ID column in
#'   `plot_metadata`.  Default `"Place"`.
#' @param init_fn Function or `NULL`.  Custom initialiser with signature
#'   `function(ref_params, n_plots, ...)`.  If `NULL` (default), uses
#'   [init_conditional()] for `mode = "conditional"` and
#'   [init_prospective_marginal()] for `mode = "prospective"`.
#' @param workers Integer scalar.  Number of parallel workers for `furrr`.
#'   Default `max(1, parallelly::availableCores() - 1)`.  Set to 1 to
#'   disable parallelism.
#' @param ... Additional arguments passed to `init_fn`.
#'
#' @return A data frame (class `"monpwr_results"`) with one row per
#'   scenario × scale × group × effect size × horizon, and columns:
#'   `scenario`, `label`, `scale`, `group`, `effect_pct`, `horizon`,
#'   `n_plots`, `n_future`, `power`, `n_converged`, `conv_rate`.
#'
#' @examples
#' \dontrun{
#' ref <- extract_params(fit, data = long_model)
#'
#' scenarios <- list(
#'   baseline = scenario("Full grid, 5-yr", 5, function(sp) sp$Place),
#'   coarse   = scenario("Half grid, 5-yr", 5, function(sp) sp$Place[sp$in_coarse])
#' )
#'
#' # Conditional (default)
#' results <- run_power_sim(
#'   ref, scenarios,
#'   plot_metadata    = fpi_spatial,
#'   data             = long_model,
#'   effect_sizes_pct = c(10, 20, 30),
#'   horizons         = c(10, 20)
#' )
#'
#' # Prospective
#' results_p <- run_power_sim(
#'   ref, scenarios,
#'   plot_metadata    = fpi_spatial,
#'   mode             = "prospective",
#'   effect_sizes_pct = c(10, 20, 30),
#'   horizons         = c(10, 20)
#' )
#' }
#'
#' @seealso [extract_params()], [scenario()], [compute_mdc()],
#'   [bootstrap_cv()], [plot_power()], [plot_mdc()]
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
                          scales           = c("National", "Park type"),
                          park_type_var    = "Park_type",
                          park_name_var    = "PCL_Name",
                          n_min_park       = 5L,
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
    abort(paste0("`place_var` column '", place_var, "' not found in `plot_metadata`."))
  }

  # Park-scale checks
  want_park_type <- "Park type" %in% scales
  want_park      <- "Park" %in% scales
  if ((want_park_type || want_park) && !park_type_var %in% names(plot_metadata)) {
    abort(paste0("`park_type_var` column '", park_type_var,
                 "' not found in `plot_metadata`. ",
                 "Set `scales` to exclude 'Park type' and 'Park', or supply the column."))
  }
  if (want_park && !park_name_var %in% names(plot_metadata)) {
    abort(paste0("`park_name_var` column '", park_name_var,
                 "' not found in `plot_metadata`. ",
                 "Set `scales` to exclude 'Park', or supply the column."))
  }

  # --- Parallel setup --------------------------------------------------------
  n_workers <- workers %||% max(1L, parallelly::availableCores() - 1L)
  if (n_workers > 1L) {
    future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # --- Default init_fn -------------------------------------------------------
  if (is.null(init_fn)) {
    init_fn <- if (mode == "conditional") init_conditional else init_prospective_marginal
  }

  # --- Outer loop ------------------------------------------------------------
  all_results <- vector("list", length(scenarios))

  for (sc_idx in seq_along(scenarios)) {
    sc_name <- names(scenarios)[sc_idx]
    sc      <- scenarios[[sc_idx]]
    sites   <- sc$site_selector(plot_metadata)
    n_sc    <- length(sites)

    cli::cli_alert_info(
      "Scenario {sc_idx}/{length(scenarios)}: {sc$label} | {n_sc} plots"
    )

    # Build scale list for this scenario
    scale_list <- .build_scale_list(
      sites         = sites,
      plot_metadata = plot_metadata,
      ref_params    = ref_params,
      data          = data,
      mode          = mode,
      init_fn       = init_fn,
      scales        = scales,
      park_type_var = park_type_var,
      park_name_var = park_name_var,
      n_min_park    = n_min_park,
      place_var     = place_var,
      ...
    )

    # Simulate each scale x effect x horizon cell
    sc_results <- map_dfr(scale_list, function(sl) {
      if (nrow(sl$plot_state) < 2) return(NULL)

      cli::cli_alert_info(
        "  Scale: {sl$scale} | Group: {sl$group} | n={nrow(sl$plot_state)}"
      )

      expand_grid(
        effect_pct = effect_sizes_pct,
        horizon    = horizons
      ) |>
        future_pmap_dfr(function(effect_pct, horizon) {
          n_future <- max(1L, floor(horizon / sc$remeasure_yrs))
          .run_one_cell(
            plot_state  = sl$plot_state,
            hist_dat    = sl$hist_dat,
            n_future    = n_future,
            effect_pct  = effect_pct,
            ref_params  = ref_params,
            n_iter      = n_iter,
            alpha       = alpha,
            mode        = mode,
            draw_re     = (mode == "prospective")
          ) |>
            mutate(scale = sl$scale, group = sl$group, horizon = horizon)
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
                              init_fn, scales, park_type_var, park_name_var,
                              n_min_park, place_var, ...) {

  sl <- list()

  # Helper: build one scale entry
  make_entry <- function(scale, group, these_sites) {
    ps <- if (mode == "conditional") {
      init_conditional(ref_params, these_sites)
    } else {
      init_fn(ref_params, n_plots = length(these_sites), ...)
    }

    hd <- if (mode == "conditional" && !is.null(data)) {
      build_historical(data, these_sites, place_var = place_var,
                       plotid_var = ref_params$plotid_var,
                       offset_var = ref_params$offset_var)
    } else {
      NULL
    }

    list(scale = scale, group = group, plot_state = ps, hist_dat = hd)
  }

  # National
  if ("National" %in% scales) {
    sl <- c(sl, list(make_entry("National", "All", sites)))
  }

  # Park type
  if ("Park type" %in% scales) {
    park_types <- unique(plot_metadata[[park_type_var]])
    for (pt in park_types) {
      type_sites <- plot_metadata[[place_var]][
        plot_metadata[[park_type_var]] == pt &
          plot_metadata[[place_var]] %in% sites
      ]
      if (length(type_sites) >= 2) {
        sl <- c(sl, list(make_entry("Park type", pt, type_sites)))
      }
    }
  }

  # Individual parks
  if ("Park" %in% scales && !is.null(park_name_var)) {
    park_counts <- plot_metadata |>
      filter(.data[[place_var]] %in% sites) |>
      count(.data[[park_name_var]], name = "n_plots") |>
      filter(.data$n_plots >= n_min_park)

    for (i in seq_len(nrow(park_counts))) {
      pcl        <- park_counts[[park_name_var]][i]
      pcl_sites  <- plot_metadata[[place_var]][
        plot_metadata[[park_name_var]] == pcl &
          plot_metadata[[place_var]] %in% sites
      ]
      sl <- c(sl, list(make_entry("Park", pcl, pcl_sites)))
    }
  }

  sl
}


# ------------------------------------------------------------------------------
# Internal: run one scenario x scale x effect x horizon cell
# ------------------------------------------------------------------------------

.run_one_cell <- function(plot_state, hist_dat, n_future, effect_pct,
                          ref_params, n_iter, alpha, mode, draw_re) {
  eff_log <- log(1 + effect_pct / 100)

  p_vals <- map_dbl(seq_len(n_iter), function(i) {
    future_dat <- simulate_visits(plot_state, n_future, eff_log, ref_params,
                                  draw_re = draw_re)
    combined <- if (!is.null(hist_dat)) {
      bind_rows(hist_dat, future_dat)
    } else {
      future_dat
    }
    fit_and_test(combined, ref_params)
  })

  n_conv <- sum(!is.na(p_vals))
  power  <- mean(p_vals < alpha, na.rm = TRUE)

  tibble(
    n_plots     = nrow(plot_state),
    n_future    = n_future,
    effect_pct  = effect_pct,
    power       = power,
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
    "*" = paste0("Scales:       ", paste(unique(x$scale), collapse = ", ")),
    "*" = paste0("Effect sizes: ", paste(sort(unique(x$effect_pct)), collapse = ", "), "%"),
    "*" = paste0("Horizons:     ", paste(sort(unique(x$horizon)), collapse = ", "), " yr"),
    "*" = paste0("Rows:         ", nrow(x))
  ))
  invisible(x)
}
