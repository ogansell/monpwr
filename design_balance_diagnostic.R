# =============================================================================
# design_balance_diagnostic() — generalisable design-balance diagnostic
# =============================================================================
#
# Motivation
# ----------
# When comparing a simulated/extended design against an intended one, the
# total row count (or ESS ratio = nrow(actual) / nrow(intended)) can look
# healthy while the PER-UNIT structure is badly distorted. For example,
# simr::extend(along = <time>) on an unbalanced pilot recycles the pilot's
# per-unit visit pattern — including its gaps — so some unit x time cells get
# heavy pseudo-replication and others get zero observations, yet the totals
# still come out near the intended count.
#
# The ESS ratio cannot see this; it is an aggregate. This diagnostic looks at
# the distribution of rows per (unit x time) cell, which is where the
# distortion actually lives.
#
# What "balanced" means here
# --------------------------
# The intended design is assumed to be one row per unit per time point
# (a complete, balanced panel). Any departure — cells with 0 rows
# (missingness) or cells with >1 row (pseudo-replication) — is a distortion.
# If your intended design legitimately has a different target (e.g. k>1 rows
# per cell), pass `target_per_cell = k`.
#
# This file is self-contained (base + dplyr + tidyr + ggplot2). It is written
# to be dropped into any extend-style comparison, not just the kea analysis:
# you pass the data frame and the names of the unit and time columns.
# =============================================================================

#' Diagnose per-cell balance of a (possibly recycled/extended) design
#'
#' @param data A data frame: the realised design to inspect (e.g. the output of
#'   `simr::getData(extend(...))`, or any simulated panel).
#' @param unit_var Character. Column identifying the sampling unit (plot, site,
#'   tile, subject).
#' @param time_var Character. Column identifying the time point (visit, year,
#'   wave).
#' @param target_per_cell Numeric scalar. Expected rows per unit x time cell in
#'   the *intended* design. Default 1 (a balanced one-row-per-cell panel).
#' @param label Character or NULL. Optional label for this design (used when
#'   binding several diagnostics together, e.g. "simr extend" vs "intended").
#'
#' @return A one-row data frame summarising balance, with columns:
#'   `label`, `n_units`, `n_times`, `n_cells_possible` (units x times),
#'   `n_cells_filled` (cells with >=1 row), `n_cells_empty`,
#'   `prop_empty`, `n_rows`, `mean_per_cell`, `sd_per_cell`, `cv_per_cell`,
#'   `min_per_cell`, `max_per_cell`, `ess_ratio`
#'   (n_rows / (n_units * n_times * target_per_cell)),
#'   `imbalance_index` (see Details).
#'
#' @details
#' `cv_per_cell` (SD/mean of per-cell counts over ALL possible cells, including
#' empty ones) is the headline number: it is 0 for a perfectly balanced design
#' and grows with both pseudo-replication and missingness. `prop_empty` isolates
#' the missingness component; `max_per_cell` isolates the worst pseudo-
#' replication. `imbalance_index` is `mean(abs(count - target_per_cell)) /
#' target_per_cell` over all possible cells — a single 0-is-perfect summary of
#' total departure from the intended structure.
#'
#' Note that `ess_ratio` can sit near 1 while `cv_per_cell` and `prop_empty`
#' are large; that divergence is the entire point of the diagnostic.
design_balance_diagnostic <- function(data, unit_var, time_var,
                                      target_per_cell = 1, label = NULL) {
  stopifnot(is.data.frame(data),
            unit_var %in% names(data),
            time_var %in% names(data),
            is.numeric(target_per_cell), target_per_cell > 0)
  
  units <- unique(data[[unit_var]])
  times <- unique(data[[time_var]])
  n_units <- length(units)
  n_times <- length(times)
  n_cells_possible <- n_units * n_times
  
  # Count rows per (unit x time) cell, then pad to the full grid so empty
  # cells are explicitly counted as zero (this is what the ESS ratio misses).
  filled <- as.data.frame(
    table(unit = data[[unit_var]], time = data[[time_var]]),
    stringsAsFactors = FALSE
  )
  counts <- filled$Freq  # length == n_units * n_times already (table fills grid)
  
  n_rows         <- sum(counts)
  n_cells_filled <- sum(counts > 0)
  n_cells_empty  <- sum(counts == 0)
  
  mean_pc <- mean(counts)
  sd_pc   <- stats::sd(counts)
  cv_pc   <- if (mean_pc > 0) sd_pc / mean_pc else NA_real_
  
  imbalance <- mean(abs(counts - target_per_cell)) / target_per_cell
  
  data.frame(
    label            = label %||% NA_character_,
    n_units          = n_units,
    n_times          = n_times,
    n_cells_possible = n_cells_possible,
    n_cells_filled   = n_cells_filled,
    n_cells_empty    = n_cells_empty,
    prop_empty       = round(n_cells_empty / n_cells_possible, 3),
    n_rows           = n_rows,
    mean_per_cell    = round(mean_pc, 3),
    sd_per_cell      = round(sd_pc, 3),
    cv_per_cell      = round(cv_pc, 3),
    min_per_cell     = min(counts),
    max_per_cell     = max(counts),
    ess_ratio        = round(n_rows / (n_cells_possible * target_per_cell), 3),
    imbalance_index  = round(imbalance, 3),
    stringsAsFactors = FALSE
  )
}

# Small NULL-coalescing helper if not already in scope (rlang exports `%||%`,
# but keep the script self-contained).
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a


#' Per-cell heatmap of a design's row counts
#'
#' Visualises the distortion that the scalar diagnostic summarises. Units on one
#' axis, time on the other, fill = rows in that cell. A balanced design is a
#' uniform field at `target_per_cell`; recycling/missingness shows as bright
#' (over-replicated) and empty (missing) cells.
#'
#' @param data,unit_var,time_var As in [design_balance_diagnostic()].
#' @param target_per_cell Numeric. Midpoint of the diverging fill scale.
#' @param title,subtitle Character or NULL.
#' @param max_units Integer or NULL. If set and there are more units than this,
#'   a random sample of `max_units` units is shown (keeps the heatmap legible
#'   for large designs). Sampling is seeded by `seed`.
#' @param seed Integer. For reproducible unit sampling.
#'
#' @return A ggplot object.
design_balance_heatmap <- function(data, unit_var, time_var,
                                   target_per_cell = 1,
                                   title = NULL, subtitle = NULL,
                                   max_units = 60, seed = 1) {
  stopifnot(is.data.frame(data),
            unit_var %in% names(data), time_var %in% names(data))
  
  cell <- as.data.frame(
    table(unit = data[[unit_var]], time = data[[time_var]]),
    stringsAsFactors = FALSE
  )
  names(cell)[names(cell) == "Freq"] <- "n_rows"
  # `time` came back as character from table(); restore numeric order if possible.
  suppressWarnings({
    t_num <- as.numeric(cell$time)
    if (!any(is.na(t_num))) cell$time <- t_num
  })
  
  if (!is.null(max_units)) {
    u <- unique(cell$unit)
    if (length(u) > max_units) {
      set.seed(seed)
      keep <- sample(u, max_units)
      cell <- cell[cell$unit %in% keep, , drop = FALSE]
    }
  }
  
  ggplot2::ggplot(cell, ggplot2::aes(x = factor(.data$time),
                                     y = .data$unit,
                                     fill = .data$n_rows)) +
    ggplot2::geom_tile(colour = "grey90", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      midpoint = target_per_cell, name = "Rows in cell"
    ) +
    ggplot2::labs(
      x = time_var, y = unit_var,
      title = title %||% "Per-cell row counts",
      subtitle = subtitle %||%
        paste0("Uniform = balanced (", target_per_cell,
               "/cell). Blue = missing, red = over-replicated.")
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 5),
      panel.grid  = ggplot2::element_blank()
    )
}


# =============================================================================
# Drop-in usage for the kea (or any) extend comparison
# =============================================================================
# Replace the ESS-only block with this. It works for any number of designs and
# any unit/time column names — just change `unit_col` / `time_col`.
#
# Assumes in scope: n_sites, ed_df, horizon, take_sites(), simr, lme4, ggplot2,
# purrr, dplyr. Adapt the pilot-fit + extend call to your model; everything
# after `ext <- getData(...)` is generic.

run_balance_diagnostics <- function(n_sites, ed_df, horizon, take_sites,
                                    unit_col = "tile", time_col = "year_seq",
                                    out_dir  = "power_analysis_outputs") {
  
  diag_tbl <- purrr::map_dfr(n_sites, function(n) {
    keep <- take_sites(n)
    dat_sub <- ed_df |> dplyr::filter(.data[[unit_col]] %in% keep)
    
    fit_sub <- lme4::glmer(
      hours_with_kea ~ year_seq + offset(log_n_hours) + (1 | tile),
      family = poisson, data = dat_sub
    )
    fixef(fit_sub)["year_seq"] <- log(1.05)
    ext <- simr::getData(simr::extend(fit_sub, along = time_col, n = horizon))
    
    # Diagnose BOTH designs on the same footing.
    d_simr <- design_balance_diagnostic(
      ext, unit_var = unit_col, time_var = time_col,
      target_per_cell = 1, label = "simr_extend"
    )
    # Intended design: a complete unit x horizon panel, one row per cell.
    intended <- expand.grid(
      unit = keep, time = seq_len(horizon),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    names(intended) <- c(unit_col, time_col)
    d_intended <- design_balance_diagnostic(
      intended, unit_var = unit_col, time_var = time_col,
      target_per_cell = 1, label = "intended"
    )
    
    dplyr::bind_rows(d_simr, d_intended) |>
      dplyr::mutate(n_requested = n, .before = 1)
  })
  
  print(diag_tbl)
  
  # Heatmap for the largest design (most informative; subsample units to stay legible).
  n_big   <- max(n_sites)
  keepbig <- take_sites(n_big)
  fitbig  <- lme4::glmer(
    hours_with_kea ~ year_seq + offset(log_n_hours) + (1 | tile),
    family = poisson, data = ed_df |> dplyr::filter(.data[[unit_col]] %in% keepbig)
  )
  fixef(fitbig)["year_seq"] <- log(1.05)
  extbig <- simr::getData(simr::extend(fitbig, along = time_col, n = horizon))
  
  hm <- design_balance_heatmap(
    extbig, unit_var = unit_col, time_var = time_col, target_per_cell = 1,
    title = paste0("simr extend() per-cell structure (", n_big, " units)"),
    subtitle = "Each row = a unit; each column = a time point. Recycling propagates the pilot's gaps and pile-ups."
  )
  ggplot2::ggsave(file.path(out_dir, "design_balance_heatmap.png"),
                  hm, width = 9, height = 11)
  
  cat("\nKey reading:\n")
  cat("  ESS ratio near 1 while cv_per_cell / prop_empty are large means the\n")
  cat("  TOTAL effort matches but its DISTRIBUTION does not — the distortion the\n")
  cat("  ESS ratio alone cannot see.\n")
  
  invisible(diag_tbl)
}

# Example call (uncomment in the comparison script):
# balance_tbl <- run_balance_diagnostics(n_sites, ed_df, horizon, take_sites)
# saveRDS(balance_tbl, "power_analysis_outputs/design_balance.rds")