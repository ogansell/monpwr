## monpwr — core unit tests
## Run with: testthat::test_file("tests/testthat/test-core.R")
## or devtools::test()

library(testthat)

# ==============================================================================
# Helpers — shared fake monpwr_params used across tests
# ==============================================================================

.fake_params <- function(family = "hurdle_nbinom2") {
  structure(
    list(
      beta_visit        = 0.05,
      beta_gap_cond     = 0,
      beta_gap_zi       = 0,
      disp_par          = 2,
      sigma_cond        = 0.5,
      sigma_zi          = 0.3,
      visit_gap_med     = 5,
      marginal_int_cond = 1.2,
      marginal_int_zi   = -0.5,
      family            = family,
      visit_num_var     = "visit_num",
      plotid_var        = "plotid_model",
      place_var         = "Place",
      visit_gap_var     = "visit_gap",
      plot_state        = data.frame(
        Place         = c("A", "B", "C"),
        plotid_model  = c("A", "B", "C"),
        visit_num     = c(2L, 3L, 1L),
        eta_last_cond = c(1.2, 0.9, 1.5),
        eta_last_zi   = c(-0.5, -0.3, -0.7),
        blup_cond     = c(0.1, -0.1, 0.0),
        blup_zi       = c(0.0, 0.0, 0.0),
        stringsAsFactors = FALSE
      )
    ),
    class = "monpwr_params"
  )
}


# ==============================================================================
# extract_params — validation layer
# ==============================================================================

test_that("extract_params rejects non-data.frame", {
  expect_error(
    monpwr:::.validate_extract_inputs(tibble::tibble(x = 1), "visit_num", "plotid_model", "Place"),
    "plain data.frame"
  )
})

test_that("extract_params.default gives informative error mentioning custom-extractor", {
  fake_fit <- structure(list(), class = "unsupported_class")
  expect_error(extract_params(fake_fit, data = data.frame()), "No `extract_params` method")
  expect_error(extract_params(fake_fit, data = data.frame()), "custom-extractor",
               info = "error message should mention custom-extractor vignette")
})

test_that(".extract_marginal_intercepts recovers baseline correctly", {
  ps <- data.frame(
    visit_num     = c(1L, 1L, 2L),
    eta_last_cond = c(1.25, 1.15, 1.35),
    eta_last_zi   = c(-0.5, -0.4, -0.6)
  )
  beta_visit <- 0.10
  # Expected: mean(1.25, 1.15) - 0.10 * (1 - 1) = mean(1.25, 1.15) = 1.20
  out <- .extract_marginal_intercepts(ps, beta_visit)
  expect_equal(out$marginal_int_cond, 1.20, tolerance = 1e-6)
  expect_equal(out$marginal_int_zi,   mean(c(-0.5, -0.4)), tolerance = 1e-6)
})


# ==============================================================================
# scenario() constructor
# ==============================================================================

test_that("scenario() validates site_selector must be a function", {
  expect_error(scenario("label", 5, "not_a_function"), "function")
})

test_that("scenario() validates remeasure_yrs must be positive", {
  expect_error(scenario("label", -1, function(x) x), "positive")
  expect_error(scenario("label",  0, function(x) x), "positive")
})

test_that("scenario() validates n_new_sites must be non-negative", {
  expect_error(scenario("label", 5, function(x) x, n_new_sites = -1L), "non-negative")
})

test_that("scenario() returns monpwr_scenario with correct fields", {
  sc <- scenario("test", 5, function(sp) sp$Place, n_new_sites = 10L,
                 eta_offset_cond = 0.2)
  expect_s3_class(sc, "monpwr_scenario")
  expect_equal(sc$label, "test")
  expect_equal(sc$remeasure_yrs, 5)
  expect_equal(sc$n_new_sites, 10L)
  expect_equal(sc$eta_offset_cond, 0.2)
  expect_true(is.function(sc$site_selector))
  expect_true(is.function(sc$new_site_init_fn))
})

test_that("scenario() defaults n_new_sites to 0", {
  sc <- scenario("test", 5, function(sp) sp$Place)
  expect_equal(sc$n_new_sites, 0L)
})


# ==============================================================================
# init_conditional()
# ==============================================================================

test_that("init_conditional filters to requested site_ids", {
  ref <- .fake_params()
  out <- init_conditional(ref, c("A", "C"))
  expect_equal(nrow(out), 2)
  expect_true(all(out$Place %in% c("A", "C")))
})

test_that("init_conditional warns on missing site_ids", {
  ref <- .fake_params()
  expect_warning(init_conditional(ref, c("A","B","MISSING")), regexp = "not found")
})

test_that("init_conditional errors with fewer than 2 plots", {
  ref <- .fake_params()
  # Suppress warning from missing IDs first
  expect_error(suppressWarnings(init_conditional(ref, c("A", "MISSING2", "MISSING3"))),
               "Fewer than 2")
})


# ==============================================================================
# init_new_sites()
# ==============================================================================

test_that("init_new_sites returns NULL for n_new = 0", {
  ref <- .fake_params()
  expect_null(init_new_sites(ref, n_new = 0L))
})

test_that("init_new_sites returns correct structure for n_new > 0", {
  set.seed(1)
  ref <- .fake_params()
  out <- init_new_sites(ref, n_new = 5L)
  expect_equal(nrow(out), 5)
  expect_true(all(out$visit_num == 0L))
  expect_named(out, c("Place", "plotid_model", "visit_num",
                       "eta_last_cond", "eta_last_zi", "blup_cond", "blup_zi"))
})

test_that("init_new_sites BLUPs are centred near 0 on average", {
  set.seed(42)
  ref <- .fake_params()
  out <- init_new_sites(ref, n_new = 500L)
  expect_lt(abs(mean(out$blup_cond)), 0.1)
})

test_that("init_new_sites eta_offset shifts marginal intercept", {
  set.seed(1)
  ref  <- .fake_params()
  off  <- 0.5
  out1 <- init_new_sites(ref, n_new = 200L, eta_offset_cond = 0)
  out2 <- init_new_sites(ref, n_new = 200L, eta_offset_cond = off)
  # Mean eta_last_cond should differ by approximately `off`
  expect_equal(mean(out2$eta_last_cond) - mean(out1$eta_last_cond), off,
               tolerance = 0.1)
})

test_that("init_new_sites errors for negative n_new", {
  ref <- .fake_params()
  expect_error(init_new_sites(ref, n_new = -1L), "non-negative")
})

test_that("init_new_sites IDs have expected prefix", {
  ref <- .fake_params()
  out <- init_new_sites(ref, n_new = 3L, id_prefix = "NEWSITE_")
  expect_true(all(startsWith(out$Place, "NEWSITE_")))
})


# ==============================================================================
# init_prospective_marginal()
# ==============================================================================

test_that("init_prospective_marginal returns correct structure", {
  ref <- .fake_params()
  out <- init_prospective_marginal(ref, n_plots = 10L)
  expect_equal(nrow(out), 10)
  expect_true(all(out$visit_num == 0L))
  expect_named(out, c("plotid", "visit_num", "eta_last_cond", "eta_last_zi"))
})

test_that("init_prospective_marginal uses marginal_int_cond from ref_params", {
  ref <- .fake_params()
  out <- init_prospective_marginal(ref, n_plots = 5L)
  expect_equal(unique(out$eta_last_cond), ref$marginal_int_cond)
})

test_that("init_prospective_marginal errors for n_plots < 2", {
  ref <- .fake_params()
  expect_error(init_prospective_marginal(ref, n_plots = 1L), ">= 2")
})


# ==============================================================================
# simulate_visits()
# ==============================================================================

test_that("simulate_visits returns correct row count and columns", {
  ref <- .fake_params("poisson")
  ps  <- data.frame(
    plotid        = c("A", "B", "C"),
    Place         = c("A", "B", "C"),
    visit_num     = c(2L, 3L, 1L),
    eta_last_cond = c(1.0, 0.8, 1.2),
    eta_last_zi   = c(0.0, 0.0, 0.0)
  )
  out <- simulate_visits(ps, n_future = 2L, eff_log = log(1.2),
                         ref_params = ref, draw_re = FALSE)
  expect_equal(nrow(out), 6)   # 3 plots × 2 visits
  expect_named(out, c("plotid", "visit_num", "log_effort", "count", "source"))
  expect_true(all(out$source == "future"))
  expect_true(all(out$count >= 0))
})

test_that("simulate_visits increments visit_num correctly", {
  ref <- .fake_params("poisson")
  ps  <- data.frame(
    plotid        = "A",
    Place         = "A",
    visit_num     = 3L,
    eta_last_cond = 1.0,
    eta_last_zi   = 0.0
  )
  out <- simulate_visits(ps, n_future = 3L, eff_log = 0,
                         ref_params = ref, draw_re = FALSE)
  expect_equal(out$visit_num, c(4L, 5L, 6L))
})

test_that("simulate_visits draw_re = TRUE samples different BLUPs each call", {
  set.seed(1)
  ref <- .fake_params("poisson")
  ps  <- data.frame(
    plotid = "A", Place = "A", visit_num = 0L,
    eta_last_cond = 0.0, eta_last_zi = 0.0
  )
  counts1 <- simulate_visits(ps, 5L, 0, ref, draw_re = TRUE)$count
  counts2 <- simulate_visits(ps, 5L, 0, ref, draw_re = TRUE)$count
  # With RE sampling, trajectories will differ almost surely
  expect_false(identical(counts1, counts2))
})


# ==============================================================================
# .draw_counts() — internal family dispatch
# ==============================================================================

test_that(".draw_counts handles all supported families", {
  set.seed(1)
  expect_true(all(monpwr:::.draw_counts("hurdle_nbinom2", rep(2, 5), rep(0.2, 5), 5, 1.5) >= 0))
  expect_true(all(monpwr:::.draw_counts("nbinom2",        rep(2, 5), rep(0.2, 5), 5, 1.5) >= 0))
  expect_true(all(monpwr:::.draw_counts("poisson",        rep(2, 5), rep(0.0, 5), 5, 1)   >= 0))
  expect_true(all(monpwr:::.draw_counts("binomial",       rep(0.5,5),rep(0.0, 5), 5, 1)   %in% c(0, 1)))
  expect_length(monpwr:::.draw_counts("gaussian",         rep(1, 5), rep(0, 5),   5, 0.3), 5)
})

test_that(".draw_counts hurdle_nbinom2 never produces negative counts", {
  set.seed(99)
  counts <- monpwr:::.draw_counts("hurdle_nbinom2", rep(3, 100), rep(0.3, 100), 100, 2)
  expect_true(all(counts >= 0))
})


# ==============================================================================
# build_historical()
# ==============================================================================

test_that("build_historical returns correct columns and filters to site_ids", {
  dat <- data.frame(
    Place        = c("A", "A", "B", "C"),
    plotid_model = c("A", "A", "B", "C"),
    visit_num    = c(1L, 2L, 1L, 1L),
    log_effort   = c(0, 0, 0, 0),
    count        = c(3L, 5L, 0L, 2L)
  )
  out <- build_historical(dat, c("A", "C"))
  expect_equal(nrow(out), 3)
  expect_true(all(out$plotid %in% c("A", "C")))
  expect_named(out, c("plotid", "visit_num", "log_effort", "count", "source"))
  expect_true(all(out$source == "observed"))
})


# ==============================================================================
# compute_mdc()
# ==============================================================================

test_that("compute_mdc returns NA when power target not reached", {
  fake_results <- structure(
    data.frame(
      scenario    = "sc1", label = "Test", sim_type = "combined",
      scale = "Overall", group = "All", horizon = 10,
      effect_pct = c(10, 20, 30), power = c(0.4, 0.6, 0.7),
      n_legacy = 50, n_new = 0L, n_total = 50, n_future = 2,
      n_converged = 200, conv_rate = 1.0
    ),
    class = c("monpwr_results", "data.frame")
  )
  mdc <- compute_mdc(fake_results, power_target = 0.80)
  expect_true(is.na(mdc$mdc_pct))
  expect_equal(round(mdc$max_power, 1), 0.7)
})

test_that("compute_mdc returns correct MDC when power target reached", {
  fake_results <- structure(
    data.frame(
      scenario = "sc1", label = "Test", sim_type = "combined",
      scale = "Overall", group = "All", horizon = 10,
      effect_pct = c(10, 20, 30), power = c(0.6, 0.82, 0.95),
      n_legacy = 50, n_new = 0L, n_total = 50, n_future = 2,
      n_converged = 200, conv_rate = 1.0
    ),
    class = c("monpwr_results", "data.frame")
  )
  mdc <- compute_mdc(fake_results, power_target = 0.80)
  expect_equal(mdc$mdc_pct, 20)
})

test_that("compute_mdc handles results without sim_type column (backward compat)", {
  fake_results <- structure(
    data.frame(
      scenario = "sc1", label = "Test",
      scale = "Overall", group = "All", horizon = 10,
      effect_pct = c(10, 20, 30), power = c(0.6, 0.82, 0.95),
      n_plots = 50, n_future = 2, n_converged = 200, conv_rate = 1.0
    ),
    class = c("monpwr_results", "data.frame")
  )
  expect_no_error(compute_mdc(fake_results))
})

test_that("compute_mdc groups by sim_type for hybrid results", {
  fake_results <- structure(
    data.frame(
      scenario   = "sc1",
      label      = "Hybrid",
      sim_type   = rep(c("combined", "legacy_only"), each = 3),
      scale      = "Overall",
      group      = "All",
      horizon    = 10,
      effect_pct = rep(c(10, 20, 30), 2),
      power      = c(0.75, 0.88, 0.95, 0.5, 0.65, 0.75),
      n_legacy   = 30, n_new = 20L, n_total = 50, n_future = 2,
      n_converged = 200, conv_rate = 1.0
    ),
    class = c("monpwr_results", "data.frame")
  )
  mdc <- compute_mdc(fake_results, power_target = 0.80)
  expect_equal(nrow(mdc), 2)
  expect_equal(mdc$mdc_pct[mdc$sim_type == "combined"],    20)
  expect_true(is.na(mdc$mdc_pct[mdc$sim_type == "legacy_only"]))
})


# ==============================================================================
# bootstrap_cv()
# ==============================================================================

test_that("bootstrap_cv returns correct structure", {
  set.seed(1)
  dat <- data.frame(
    Place  = rep(c("A", "B", "C", "D"), each = 3),
    Season = 2024,
    count  = c(1,0,1, 0,0,0, 1,1,0, 1,0,1)
  )
  out <- bootstrap_cv(dat, c("A","B","C","D"), n_boot = 100)
  expect_named(out, c("n_sites", "mean_resp", "cv"))
  expect_equal(out$n_sites, 4)
  expect_true(out$cv > 0)
})

test_that("bootstrap_cv returns NA for fewer than 2 sites", {
  dat <- data.frame(Place = "A", Season = 2024, count = c(1, 0))
  out <- bootstrap_cv(dat, "A", n_boot = 50)
  expect_true(is.na(out$mean_resp))
  expect_true(is.na(out$cv))
})

test_that("bootstrap_cv accepts a custom response_fn", {
  dat <- data.frame(
    Place  = rep(c("A", "B"), each = 4),
    Season = 2024,
    count  = c(2, 4, 6, 8, 1, 3, 5, 7)
  )
  out <- bootstrap_cv(dat, c("A", "B"), n_boot = 100,
                      response_fn = function(d) mean(d$count))
  expect_equal(out$mean_resp, mean(c(mean(c(2,4,6,8)), mean(c(1,3,5,7)))))
})


# ==============================================================================
# reporting_groups — scenario summary helper
# ==============================================================================

test_that("reporting_groups columns are validated in run_power_sim", {
  ref  <- .fake_params()
  scns <- list(sc1 = scenario("sc1", 5, function(sp) sp$Place))
  meta <- data.frame(Place = c("A", "B", "C"))

  expect_error(
    run_power_sim(
      ref_params       = ref,
      scenarios        = scns,
      plot_metadata    = meta,
      data             = data.frame(Place = "A", plotid_model = "A",
                                     visit_num = 1L, log_effort = 0, count = 1L),
      reporting_groups = list("Region" = "missing_col"),
      workers          = 1L
    ),
    "missing_col"
  )
})
