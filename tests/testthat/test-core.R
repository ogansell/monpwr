test_that("extract_params rejects non-data.frame", {
  # Create a minimal fake glmmTMB-like object just to test dispatch
  # Real model tests require the model to be fitted; these test the validation layer
  expect_error(
    .validate_extract_inputs(tibble::tibble(x = 1), "visit_num", "plotid_model", "Place"),
    "plain data.frame"
  )
})

test_that("extract_params.default gives informative error", {
  fake_fit <- structure(list(), class = "unsupported_class")
  expect_error(
    extract_params(fake_fit, data = data.frame()),
    "No `extract_params` method"
  )
  expect_error(
    extract_params(fake_fit, data = data.frame()),
    "custom-extractor"
  )
})

test_that("scenario() validates inputs", {
  expect_error(scenario("label", 5, "not_a_function"), "function")
  expect_error(scenario("label", -1, function(x) x), "positive")
  expect_s3_class(
    scenario("test", 5, function(sp) sp$Place),
    "monpwr_scenario"
  )
})

test_that("init_prospective_marginal returns correct structure", {
  # Build a minimal monpwr_params without a real model
  fake_params <- structure(
    list(
      beta_visit        = 0.05,
      beta_gap_cond     = 0,
      beta_gap_zi       = 0,
      disp_par          = 2,
      sigma_cond        = 0.5,
      sigma_zi          = 0.3,
      visit_gap_med     = 5,
      family            = "hurdle_nbinom2",
      visit_num_var     = "visit_num",
      plotid_var        = "plotid_model",
      place_var         = "Place",
      offset_var        = NULL,
      log_effort_future = 0,
      plot_state    = data.frame(
        Place        = c("A", "B"),
        plotid_model = c("A", "B"),
        visit_num    = c(2L, 3L),
        eta_last_cond = c(1.2, 0.9),
        eta_last_zi   = c(-0.5, -0.3),
        blup_cond    = c(0.1, -0.1),
        blup_zi      = c(0.0, 0.0)
      )
    ),
    class = "monpwr_params"
  )

  out <- init_prospective_marginal(fake_params, n_plots = 10)
  expect_equal(nrow(out), 10)
  expect_true(all(out$visit_num == 0L))
  expect_named(out, c("plotid", "visit_num", "eta_last_cond", "eta_last_zi"))
})

test_that("simulate_visits returns correct structure", {
  fake_params <- structure(
    list(
      beta_visit        = 0.05,
      beta_gap_cond     = 0,
      beta_gap_zi       = 0,
      disp_par          = 2,
      sigma_cond        = 0.5,
      sigma_zi          = 0.3,
      visit_gap_med     = 5,
      family            = "poisson",
      visit_num_var     = "visit_num",
      plotid_var        = "plotid_model",
      place_var         = "Place",
      offset_var        = NULL,
      log_effort_future = 0,
      plot_state    = data.frame()
    ),
    class = "monpwr_params"
  )

  plot_state <- data.frame(
    plotid        = c("A", "B", "C"),
    Place         = c("A", "B", "C"),
    visit_num     = c(2L, 3L, 1L),
    eta_last_cond = c(1.0, 0.8, 1.2),
    eta_last_zi   = c(0.0, 0.0, 0.0)
  )

  out <- simulate_visits(plot_state, n_future = 2, eff_log = log(1.2),
                         ref_params = fake_params, draw_re = FALSE)

  expect_equal(nrow(out), 6)  # 3 plots x 2 future visits
  expect_named(out, c("plotid", "visit_num", "log_effort", "count", "source"))
  expect_true(all(out$source == "future"))
  expect_true(all(out$count >= 0))
  expect_true(all(out$log_effort == 0))  # no offset case
})

test_that("simulate_visits uses log_effort_future from ref_params", {
  fake_params <- structure(
    list(
      beta_visit        = 0.05,
      beta_gap_cond     = 0,
      beta_gap_zi       = 0,
      disp_par          = 2,
      sigma_cond        = 0.5,
      sigma_zi          = 0,
      visit_gap_med     = 5,
      family            = "poisson",
      visit_num_var     = "visit_num",
      plotid_var        = "plotid_model",
      place_var         = "Place",
      offset_var        = "log_n_hours",
      log_effort_future = log(8),
      plot_state        = data.frame()
    ),
    class = "monpwr_params"
  )

  plot_state <- data.frame(
    plotid        = c("A", "B"),
    Place         = c("A", "B"),
    visit_num     = c(1L, 2L),
    eta_last_cond = c(1.0, 0.8),
    eta_last_zi   = c(0.0, 0.0)
  )

  out <- simulate_visits(plot_state, n_future = 2, eff_log = log(1.1),
                         ref_params = fake_params, draw_re = FALSE)

  expect_true(all(out$log_effort == log(8)))
})

test_that(".draw_counts handles all families", {
  set.seed(1)
  expect_true(all(.draw_counts("hurdle_nbinom2", rep(2, 5), rep(0.2, 5), 5, 1.5) >= 0))
  expect_true(all(.draw_counts("nbinom2",        rep(2, 5), rep(0.2, 5), 5, 1.5) >= 0))
  expect_true(all(.draw_counts("poisson",        rep(2, 5), rep(0.0, 5), 5, 1)   >= 0))
  expect_true(all(.draw_counts("binomial",       rep(0.5,5),rep(0.0, 5), 5, 1)   %in% c(0,1)))
  # Gaussian can be negative (log-scale index)
  expect_length(.draw_counts("gaussian", rep(1, 5), rep(0, 5), 5, 0.3), 5)
})

test_that("compute_mdc returns NA when power target not reached", {
  fake_results <- structure(
    data.frame(
      scenario   = "sc1",
      label      = "Test",
      scale      = "National",
      group      = "All",
      horizon    = 10,
      effect_pct = c(10, 20, 30),
      power      = c(0.4, 0.6, 0.7),
      n_plots    = 50,
      n_future   = 2,
      n_converged = 200,
      conv_rate  = 1.0
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
      scenario   = "sc1",
      label      = "Test",
      scale      = "National",
      group      = "All",
      horizon    = 10,
      effect_pct = c(10, 20, 30),
      power      = c(0.6, 0.82, 0.95),
      n_plots    = 50,
      n_future   = 2,
      n_converged = 200,
      conv_rate  = 1.0
    ),
    class = c("monpwr_results", "data.frame")
  )

  mdc <- compute_mdc(fake_results, power_target = 0.80)
  expect_equal(mdc$mdc_pct, 20)
})
