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
      visit_gap_var     = "visit_gap",
      count_var         = "count",
      offset_var        = NULL,
      log_effort_future = 0,
      plot_state    = data.frame(
        place_id      = c("A", "B"),
        plotid        = c("A", "B"),
        visit_num     = c(2L, 3L),
        eta_last_cond = c(1.2, 0.9),
        eta_last_zi   = c(-0.5, -0.3),
        blup_cond     = c(0.1, -0.1),
        blup_zi       = c(0.0, 0.0)
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
      visit_gap_var     = "visit_gap",
      count_var         = "count",
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
      visit_gap_var     = "visit_gap",
      count_var         = "count",
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
      scenario    = "sc1",
      label       = "Test",
      group       = "All",
      horizon     = 10,
      effect_pct  = c(10, 20, 30),
      power       = c(0.4, 0.6, 0.7),
      power_lower = c(0.33, 0.53, 0.63),
      power_upper = c(0.47, 0.67, 0.77),
      n_plots     = 50,
      n_future    = 2,
      n_converged = 200,
      conv_rate   = 1.0
    ),
    class = c("monpwr_results", "data.frame")
  )

  mdc <- compute_mdc(fake_results, power_target = 0.80)
  expect_true(is.na(mdc$mdc_pct))
  expect_equal(round(mdc$max_power, 1), 0.7)
})

test_that("retest recalculates power at new alpha", {
  fake_results <- structure(
    data.frame(
      scenario    = "sc1",
      label       = "Test",
      group       = "All",
      horizon     = 10,
      effect_pct  = 20,
      power       = 0.5,
      power_lower = 0.4,
      power_upper = 0.6,
      n_plots     = 50,
      n_future    = 2,
      n_converged = 10,
      conv_rate   = 1.0
    ),
    class = c("monpwr_results", "data.frame")
  )
  fake_results$p_values <- list(c(0.01, 0.03, 0.05, 0.08, 0.12,
                                  0.15, 0.20, 0.50, 0.80, 0.95))

  retested <- retest(fake_results, alpha = 0.05)
  expect_equal(retested$power, 0.2)

  retested2 <- retest(fake_results, alpha = 0.20)
  expect_equal(retested2$power, 0.6)
})

test_that("retest aborts when p_values missing", {
  fake <- structure(
    data.frame(scenario = "s", label = "L", group = "All",
               horizon = 10, effect_pct = 10, power = 0.5,
               n_plots = 10, n_future = 2, n_converged = 10, conv_rate = 1),
    class = c("monpwr_results", "data.frame")
  )
  expect_error(retest(fake, alpha = 0.05), "p_values")
})

test_that("extend concatenates p-values and recomputes power", {
  make_result <- function(pv) {
    out <- structure(
      data.frame(
        scenario    = "sc1",
        label       = "Test",
        group       = "All",
        horizon     = 10,
        effect_pct  = 20,
        power       = NA_real_,
        power_lower = NA_real_,
        power_upper = NA_real_,
        n_plots     = 50,
        n_future    = 2,
        n_converged = length(pv),
        conv_rate   = 1.0
      ),
      class = c("monpwr_results", "data.frame")
    )
    out$p_values <- list(pv)
    out
  }

  r1 <- make_result(c(0.01, 0.05, 0.20, 0.50, 0.90))
  r2 <- make_result(c(0.02, 0.04, 0.06, 0.30, 0.80))

  combined <- extend(r1, r2, alpha = 0.10)
  expect_equal(length(combined$p_values[[1]]), 10)
  expect_equal(combined$power, 0.5)
})

test_that("fit_and_test returns a p-value despite benign warnings", {
  skip_if_not_installed("lme4")
  set.seed(2)
  d <- data.frame(
    plotid     = factor(rep(paste0("p", 1:6), each = 3)),
    visit_num  = rep(1:3, times = 6),
    log_effort = 0,
    count      = rpois(18, lambda = 2)
  )
  ref <- structure(list(
    beta_visit = 0.05, beta_gap_cond = 0, beta_gap_zi = 0, disp_par = 1,
    sigma_cond = 0.01, sigma_zi = 0, visit_gap_med = 0, family = "poisson",
    visit_num_var = "visit_num", plotid_var = "plotid", place_var = "plotid",
    visit_gap_var = "visit_gap", count_var = "count", offset_var = NULL,
    log_effort_future = 0, plot_state = data.frame()
  ), class = "monpwr_params")
  p <- fit_and_test(d, ref, test = "wald")
  expect_true(is.na(p) || (p >= 0 && p <= 1))
})

test_that(".run_one_cell reports power_all alongside convergence-conditioned power", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  set.seed(3)
  d <- expand.grid(plot = factor(paste0("p", 1:10)), visit = 1:4)
  re <- rnorm(10, 0, 0.5)
  d$count <- rpois(nrow(d), exp(0.5 + log(1.05) * d$visit + re[as.integer(d$plot)]))
  d <- as.data.frame(d)
  fit <- lme4::glmer(count ~ visit + (1 | plot), family = poisson, data = d)
  ref <- extract_params(fit, data = d, visit_num_var = "visit",
                        plotid_var = "plot", place_var = "plot")
  scenarios <- list(
    test = scenario("Test", 5, function(sp) unique(sp$plot))
  )
  plot_meta <- data.frame(plot = unique(d$plot), stringsAsFactors = FALSE)
  res <- run_power_sim(ref, scenarios, plot_metadata = plot_meta,
                       mode = "prospective", effect_sizes_pct = 10,
                       horizons = 10, n_iter = 20, alpha = 0.10,
                       place_var = "plot", workers = 1)
  expect_true(all(c("power", "power_all", "conv_rate", "n_converged",
                     "power_all_lower", "power_all_upper") %in% names(res)))
  expect_true(all(res$power_all >= 0 & res$power_all <= 1))
})

test_that("calibrate_bias returns the expected fields", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  set.seed(1)
  d <- expand.grid(plot = factor(paste0("p", 1:12)), visit = 1:5)
  re <- rnorm(12, 0, 0.6)
  d$count <- rpois(nrow(d), exp(0.4 + log(1.05) * d$visit + re[as.integer(d$plot)]))
  d <- as.data.frame(d)
  fit <- lme4::glmer(count ~ visit + (1 | plot), family = poisson, data = d)
  ref <- extract_params(fit, data = d, visit_num_var = "visit",
                        plotid_var = "plot", place_var = "plot")
  cal <- calibrate_bias(ref, n_plots = 20, n_visits = 6, effect_pct = 5,
                        n_cal = 20, n_pilot = 12)
  expect_true(all(c("monpwr_power", "truth_power", "bias",
                    "truth_ci", "n_cal", "n_pilot", "sigma_cond") %in% names(cal)))
  expect_true(is.na(cal$bias) || (cal$bias >= -1 && cal$bias <= 1))
})

test_that("compute_mdc returns correct MDC when power target reached", {
  fake_results <- structure(
    data.frame(
      scenario    = "sc1",
      label       = "Test",
      group       = "All",
      horizon     = 10,
      effect_pct  = c(10, 20, 30),
      power       = c(0.6, 0.82, 0.95),
      power_lower = c(0.53, 0.76, 0.91),
      power_upper = c(0.67, 0.87, 0.98),
      n_plots     = 50,
      n_future    = 2,
      n_converged = 200,
      conv_rate   = 1.0
    ),
    class = c("monpwr_results", "data.frame")
  )

  mdc <- compute_mdc(fake_results, power_target = 0.80)
  expect_equal(mdc$mdc_pct, 20)
})
