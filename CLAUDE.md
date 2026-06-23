# CLAUDE.md — monpwr

## What this package does

`monpwr` estimates statistical power and minimum detectable change (MDC) for
biodiversity monitoring programme design scenarios. It supports two simulation
modes:

- **Conditional** — plots initialised from their observed visit history and
  empirical BLUP. Future visits project forward from that starting point. The
  existing temporal record contributes to the trend estimate. This is the
  correct mode when evaluating which plots in an *existing* programme to
  continue visiting.
- **Prospective** — plots initialised from scratch using model parameters only
  (no existing record). Equivalent to `simr`-style analysis. Used for new
  programme design or as a comparison baseline.

The **conditional vs prospective gap** is the core intellectual contribution.
`simr` is structurally prospective (initialises every plot from visit zero).
`monpwr`'s conditional mode initialises each plot from its empirically derived
BLUP, so accumulated monitoring investment is correctly reflected in power
estimates. An empirical comparison of the two modes on real data is the
planned paper anchor. Do not make architectural changes that conflate the
two modes or make the comparison harder to run.

Primary application: NZ DOC Tier 1 Forest & Bird Index (FPI). The package
must remain fully generic — no FPI-specific or application-specific assumptions
anywhere in the package internals.


## File map

```
R/
  extract_params.R   — S3 dispatch layer; glmmTMB + lme4 methods; monpwr_params contract
  initialisers.R     — init_conditional(), init_prospective_marginal(), init_new_sites()
  simulation.R       — simulate_visits(), build_historical(), fit_and_test()
  run_power_sim.R    — scenario(), run_power_sim() outer loop, .build_scale_list(), .run_one_cell()
  summaries.R        — compute_mdc(), bootstrap_cv()
  plots.R            — plot_power(), plot_mdc(), plot_cv(), run_precision()
  monpwr-package.R   — package-level doc and @importFrom declarations

tests/
  testthat/test-core.R   — stub-based unit tests (no real model fits)

vignettes/
  getting-started.Rmd    — user-facing workflow vignette
```

Analysis documents (outside the package, in the project root or a sibling directory):
- `tier1_design_scenario_power_conditional.Rmd` — live FPI design scenario analysis


## The monpwr_params contract

`extract_params()` returns a list of class `"monpwr_params"`. Every method
must return **all** of these fields. Downstream functions key off this
structure and are model-agnostic.

```r
list(
  beta_visit        = <numeric>   # fixed-effect slope for visit sequence variable
  beta_gap_cond     = <numeric>   # visit_gap slope, conditional component (0 if absent)
  beta_gap_zi       = <numeric>   # visit_gap slope, ZI component (0 if absent/no ZI)
  disp_par          = <numeric>   # NB2 phi; 1 for Poisson/binomial; residual SD for Gaussian
  sigma_cond        = <numeric>   # SD of plot-level RE, conditional component
  sigma_zi          = <numeric>   # SD of plot-level RE, ZI component (0 if no ZI)
  visit_gap_med     = <numeric>   # median visit gap across data (0 if not in model)
  family            = <character> # one of: "hurdle_nbinom2", "nbinom2", "poisson",
                                  #         "binomial", "gaussian"
  visit_num_var     = <character> # column name of visit sequence in data
  plotid_var        = <character> # column name of RE grouping factor in data
  place_var         = <character> # column name of plot identifier in data
  visit_gap_var     = <character> # column name of visit gap variable (default "visit_gap")
  count_var         = <character> # column name of response variable; auto-detected from
                                  # model formula LHS; falls back to "count" if LHS is complex
  offset_var        = <character> # column name of pre-computed log-offset, or NULL
  log_effort_future = <numeric>   # scalar assumed log-offset for simulated future visits
                                  # (0 if no offset; median observed if offset present)
  plot_state        = <data.frame> # one row per plot; see plot_state contract below
)
```

### plot_state contract

One row per unique sampling unit. Required columns:

| Column         | Type      | Description                                      |
|----------------|-----------|--------------------------------------------------|
| `place_id`     | character | Plot identifier (renamed from user's `place_var` column) |
| `plotid`       | character | RE grouping ID (renamed from user's `plotid_var` column) |
| `visit_num`    | integer   | Last observed visit number                       |
| `eta_last_cond`| numeric   | Fitted linear predictor at last visit (cond)     |
| `eta_last_zi`  | numeric   | Fitted linear predictor at last visit (ZI)       |
| `blup_cond`    | numeric   | Plot-level BLUP, conditional component           |
| `blup_zi`      | numeric   | Plot-level BLUP, ZI component (0 if no ZI)      |

`eta_last_cond` and `eta_last_zi` are obtained via `predict(..., re.form = NULL)`
on the last observed visit row. They absorb all fixed covariate effects (park,
forest, spatial) and the BLUP, so no further covariate marginalisation is needed
in simulation.


## Offset handling

The package supports three offset situations:

1. **No offset** — `offset_var = NULL`, `offset_transform = NULL` (default).
   `log_effort_future = 0` throughout. Test model's `offset(log_effort)` term
   evaluates to 1 everywhere. Equivalent to no offset.

2. **Pre-computed log-offset column** — `offset_var = "my_log_col"`. The named
   column is read from `data` in `extract_params()` and `build_historical()`.

3. **Inline offset** (e.g. `offset(log(n_hours))`) — `offset_transform = function(d) log(d$n_hours)`.
   The function is applied to `data` in `extract_params()` to materialise values.

In all cases, `log_effort_future` (stored in `ref_params`) is the scalar used
for all simulated future visits in `simulate_visits()`. Default is the median
of observed log-effort values. Override by passing `log_effort_future = log(8)`
(or any scalar) to `extract_params()`.

`build_historical()` receives `offset_var` from `ref_params$offset_var` via
`.build_scale_list()`. It materialises the offset column before filtering to
avoid row-index misalignment.


## Simulation mechanics

### Conditional mode

1. `init_conditional(ref_params, site_ids)` — filters `ref_params$plot_state`
   to retained sites. Each plot starts at its last observed `visit_num` and
   `eta_last_cond`/`eta_last_zi`.

2. `build_historical(data, site_ids, ...)` — extracts the observed record for
   retained plots. Stacked *above* simulated future data so the trend model
   sees the full temporal record.

3. `simulate_visits(plot_state, n_future, eff_log, ref_params, draw_re = FALSE)`
   — increments `visit_num` from each plot's current position. The LP increment
   is `eff_log * steps` — `eff_log` is the full forward slope applied from
   `eta_last_cond`. For prospective mode the intercept is already de-trended
   (beta_visit stripped by `init_prospective_marginal`), so `eff_log` is the
   complete trend. For conditional mode, `eta_last_cond` is the fitted value
   at the last observed visit; `eff_log` is the hypothetical future slope.
   `draw_re = FALSE` — random effects already absorbed into `eta_last_cond`.

4. `fit_and_test(combined, ref_params, test)` — fits a simplified test model
   (random intercept only, no fixed covariates) to observed + simulated data.
   Returns p-value for `visit_num`; `NA` on convergence failure.
   `test = "wald"` (default) uses the z/t-statistic from model coefficients.
   `test = "lrt"` fits a null model (no `visit_num`) and computes a
   likelihood-ratio chi-squared test — more reliable at small sample sizes
   but ~2× slower (two model fits per replicate).

### Raw p-values, retest, and extend

`.run_one_cell()` stores every replicate's raw p-value in a list-column
`p_values` on the results tibble. This enables:

- `retest(results, alpha)` — re-summarise power at a different significance
  threshold without re-simulating.
- `extend(results, additional)` — concatenate p-values from two runs of the
  same design grid and recompute power from the pooled replicates.
  Analogous to `simr::extend()`.

Both functions recompute `power`, `power_lower`, `power_upper`, `n_converged`,
and `conv_rate` from the stored/pooled p-values.

### Prospective mode

Steps 1 and 2 differ:
1. `init_prospective_marginal(ref_params, n_plots)` — all plots start at
   `visit_num = 0` with LP set to the marginal baseline intercept, computed
   on-the-fly from `plot_state` by stripping each plot's BLUP and accumulated
   trend and averaging across plots:
   ```r
   marginal_int_cond = mean(eta_last_cond - blup_cond - beta_visit * visit_num)
   marginal_int_zi   = mean(eta_last_zi   - blup_zi)   # 0 if sigma_zi == 0
   ```
   These are **not** stored in the `monpwr_params` contract — they are derived
   inside the initialiser each call. Do not add them as contract fields.
   `init_new_sites()` uses the same derivation before adding fresh BLUPs.
2. No historical stub — `hist_dat = NULL`. Test model sees only simulated data.
3. `simulate_visits(..., draw_re = TRUE)` — fresh REs drawn per replicate from
   `N(0, sigma_cond)` and `N(0, sigma_zi)`.

### Why conditional is methodologically superior

Prospective simulation initialises every plot from visit zero, discarding
observed histories. This systematically undervalues existing monitoring
investment. Plots with 2–3 banked visits start closer to adequate power in
conditional mode — that gap is the "temporal capital" the prospective approach
ignores. Using BLUPs as point estimates (not sampling from their posterior
variance) is deliberate: sampling posterior variance subtly shifts conditional
toward prospective, undermining the framing.


## Coding conventions

### R style
- `rlang::abort()` / `rlang::warn()` with named `i`/`x` bullet vectors — never
  `stop()` / `warning()`
- `.data$col` pronoun inside all `dplyr` verbs — never bare `col`
- `cli::cli_alert_info()` for user-facing progress
- Internal helpers prefixed with `.` (e.g. `.draw_counts`, `.fit_test_model`,
  `.extract_plot_state`, `.resolve_offset`)
- `map_dfr()` / `future_pmap_dfr()` over explicit `bind_rows(map(...))`
- `stringsAsFactors = FALSE` in every `data.frame()` constructor
- No `<<-`, no `library()` inside package files, no hardcoded column names

### Parallelism
- `furrr::future_pmap_dfr()` with `furrr_options(seed = TRUE)` for the inner
  effect × horizon loop
- `pmap_dfr()` (sequential) for the outer scenario loop — `future_pmap_dfr()`
  at the outer level causes silent segfaults
- Always restore: `on.exit(future::plan(future::sequential), add = TRUE)`
- `future::multicore` on Linux; `future::multisession` on Windows/Mac

### `data` argument rules
- `data` passed to `extract_params()` **must be a plain `data.frame`** — not a
  tibble or `data.table`. BLUP extraction via `ranef()` can silently fail
  otherwise. Always enforce with `.validate_extract_inputs()`.
- The guard uses `identical(class(data), "data.frame")`, not `is.data.frame()`.
  Tibbles inherit from `data.frame` so `is.data.frame(tibble(...))` returns
  `TRUE` and would bypass the check silently.
- Call `as.data.frame(data)` in analysis scripts before passing to the package.

### Test conventions (`tests/testthat/test-core.R`)
- `testthat` 3e style
- Build minimal `monpwr_params` stubs (structure + class) — no real model fits
- Stub must include **all** fields in the `monpwr_params` contract, including
  `offset_var` and `log_effort_future`
- Test validation/`abort` paths explicitly
- `simulate_visits(..., draw_re = FALSE)` for deterministic output checks
- Currently: 0 errors, 0 warnings, 2 harmless notes from R CMD CHECK


## Key design decisions — do not reverse without discussion

1. **S3 dispatch via `extract_params()`** — all downstream code is model-agnostic.
   The `monpwr_params` object is the sole contract between extraction and simulation.
   Adding a new model class means adding a new `extract_params.myclass()` method
   that returns all required fields — no changes to simulation or plotting code.

2. **No hardcoded column names anywhere in package internals.** All user-data
   column names (`visit_num_var`, `plotid_var`, `place_var`, `visit_gap_var`,
   `count_var`, `offset_var`) are passed as arguments and stored in `ref_params`.
   The literal strings `"visit_num"`, `"plotid_model"`, `"Place"`, `"visit_gap"`,
   `"count"`, `"log_effort"` appear only as argument *defaults*, not inside logic.
   The internal `plot_state` column names (`place_id`, `plotid`, `visit_num`, etc.)
   are fixed pipeline labels — not defaults, not user-facing — and may appear as
   literals only in code that reads/writes `plot_state`.
   `plot_state` always has fixed internal column names (`place_id`, `plotid`,
   `visit_num`, `eta_last_cond`, `eta_last_zi`, `blup_cond`, `blup_zi`) — code
   that reads from `plot_state` may use these literals directly.

3. **BLUPs as point estimates.** Using point-estimate BLUPs for conditional
   simulation is correct. Sampling from the BLUP posterior variance shifts
   conditional toward prospective and undermines the core framing. Shrinkage
   is a caveat to document, not to correct architecturally.

4. **`simr` cannot support conditional mode.** `simr`'s architecture assumes a
   fixed data frame with only the response varying. `monpwr`'s conditional mode
   changes the data frame each replicate. Do not attempt to bolt onto `simr`.

5. **Retrospective power analysis excluded by design.** MDC and adequacy
   auditing are defensible when decoupled from observed test statistics. A
   `mode = "retrospective"` was considered and rejected.

6. **Test model is intentionally simpler than the data-generating model.**
   Fixed covariates (park, forest, spatial) are absorbed into the plot-level
   random intercept in the test model. This is conservative and model-agnostic.
   It is not a bug.

7. **`offset_var = NULL` default.** No offset is assumed unless explicitly
   supplied. This is a breaking change from any code that relied on a
   hardcoded `log_effort` column — callers must now pass `offset_var = "log_effort"`
   explicitly. This is intentional: it forces the user to declare their offset.


## Scope discipline

The paper anchor is an **empirical comparison of conditional vs prospective power**
on real monitoring data. That comparison has not yet been run. Do not prioritise
feature work over it.

Before adding any new feature, ask: does this serve the paper or just the package?
Be sceptical of:
- Additional model families beyond what a stress-test dataset requires
- New plotting functions before existing ones are validated on real output
- Hybrid mode implementation before conditional vs prospective is documented
- Vignette infrastructure before the core comparison is complete

Planned backlog (in priority order):
1. Conditional vs prospective comparison write-up on FPI data — **this first**
2. ~~Binomial CIs on power estimates via `binom.test`~~ — **done**
3. ~~`simulate_visits` trend formula bug~~ — **fixed** (was `(eff_log - beta_visit) * steps`,
   now `eff_log * steps`). Previous analyses underestimated power proportional to
   `beta_visit`. All kea and FPI results should be re-run.
4. `retest()` and `extend()` should preserve extra columns added by the user
   (e.g. `n_target`, `effort_hrs` from outer loops) — currently drops them
5. Investigate `simr::extend()` behaviour with unbalanced data — see note below
6. Hybrid mode — plots absent from `plot_state` initialised from marginal intercept
7. Stress test on a second dataset (bird counts, different family)
8. `trend_fn` interface replacing scalar `eff_log` with a function over visit steps


## simr::extend() investigation

During validation testing (kea data, synthetic experiments), `simr::extend()`
was found to construct the extended data by **cyclically recycling rows** from
the original data. With balanced data (equal n per time point) this produces
a clean balanced design and power estimates agree with monpwr and brute-force
ground truth. However, with unbalanced data (varying n per time point — common
in ecological monitoring), the recycled structure does not match any realistic
design:

- Kea data: original visits per year = 70, 174, 248, 245, 181. After
  `extend(along = "year_seq", n = 10)`, the extended data cycles these counts
  across 10 time points: 174, 248, 181, 70, 245, 174, 248, 181, 70, 245.
- The resulting design has unequal replication per time point, inherited from
  the original data's imbalance, which no real monitoring programme would
  implement.

**Status**: needs systematic investigation across data types:
1. Balanced (equal n per visit) — validated, simr agrees with ground truth
2. Missing visits (some plots have gaps) — not yet tested
3. Unequal plots per time (programme grew over time) — not yet tested
4. Multiple observations per plot-visit (repeat surveys) — not yet tested

**Framing for the paper**: do not claim simr is wrong. The finding is that
`extend()` — the standard convenience function — creates a data structure
that may not match the user's intended design when the original data is
unbalanced. simr produces correct results when given a properly constructed
balanced design matrix. The issue is usability, not correctness.


## Parameter estimation uncertainty

Power analysis — whether via monpwr, simr, or brute-force from estimated
parameters — is only as reliable as the parameters it's conditioned on.
Variance components (`sigma_cond`, `sigma_zi`, `disp_par`) estimated from
small pilot datasets can be substantially biased:

- 15 plots × 5 visits: `sigma_cond` estimated at 0.68 vs true 0.80 (15%
  underestimate → power overestimate)
- The bias is stochastic — sometimes over, sometimes under — but with
  small pilots, the RE variance is consistently underestimated due to
  shrinkage in mixed models

**Implications**: monpwr extracts parameters once and uses them for all
replicates. If the extracted `sigma_cond` is too small, every replicate
simulates less noisy data than reality, inflating power. simr has the
same underlying issue but it's partially masked because simr re-estimates
variance during `powerSim()`.

**This is not fixable in code** — it's a fundamental limitation of
conditioning on estimated parameters from finite data. It belongs in the
methods discussion of the paper, not as a code change.

**Practical guidance**:
- Fit the model to as much data as possible before extracting parameters.
  A model from 144 tiles is far more reliable than one from 15.
- Consider sensitivity analysis: run power with `sigma_cond * 1.2` and
  `sigma_cond * 0.8` to bracket the uncertainty.
- Report that power estimates are conditional on the estimated variance
  structure and may be optimistic if the pilot is small.


## Infrastructure notes

- **AWS EC2** (ap-southeast-2, Melbourne) for production runs
- **S3 bucket** — requires `AWS_DEFAULT_REGION` env var; use `system()` CLI
  calls rather than the `aws.s3` R package
- **Docker** `rocker/r-ver:4.5.1` for prospective runs — R 4.5.2 has a
  TMB incompatibility on bare EC2
- **GitHub token** was previously exposed; must be revoked and regenerated
  before any repository push