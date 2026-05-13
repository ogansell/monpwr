# monpwr

**Power analysis for monitoring programme design**

`monpwr` evaluates statistical power and precision for biodiversity monitoring
programme design scenarios. It supports three simulation modes and is fully
generalised — no assumption is made about response variable, covariate
structure, or grouping taxonomy.

---

## Simulation modes

| Mode | Description | Use when |
|---|---|---|
| **Conditional** | Plots initialised from observed visit history and empirical BLUP | You have existing data and are deciding which plots to continue |
| **Prospective** | Plots initialised from scratch using model parameters only | Evaluating a brand new design with no prior record |
| **Hybrid** | Mix of legacy plots (conditional) and new plots (from scratch) | Expanding an existing network with new sites |

The conditional and hybrid modes account for temporal capital already
accumulated in legacy plots. Dropping a plot with several existing visits
loses that record entirely — a cost invisible to prospective-only analysis.

---

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install monpwr from source
devtools::install("path/to/monpwr")
```

Required system dependencies: a C++ compiler (Rtools on Windows,
Xcode command line tools on Mac) is needed to install `glmmTMB`.

---

## Quick start

```r
library(monpwr)

# 1. Fit your model (glmmTMB or lme4)
fit <- glmmTMB::glmmTMB(
  count ~ visit_num + offset(log_effort) + (1 | plotid_model),
  ziformula = ~ (1 | plotid_model),
  family    = glmmTMB::truncated_nbinom2,
  data      = long_model   # must be a plain data.frame
)

# 2. Extract standardised parameters
ref <- extract_params(
  fit,
  data          = long_model,
  visit_num_var = "visit_num",
  plotid_var    = "plotid_model",
  place_var     = "site_id",
  visit_gap_var = "visit_gap"
)

# 3. Define design scenarios
scenarios <- list(
  baseline = scenario(
    label         = "Full grid, 5-yr",
    remeasure_yrs = 5,
    site_selector = function(sp) unique(sp$site_id)
  ),
  coarse = scenario(
    label         = "Coarse grid, 5-yr",
    remeasure_yrs = 5,
    site_selector = function(sp) sp$site_id[sp$in_coarse_grid]
  ),
  hybrid = scenario(
    label         = "Coarse grid + 30 new sites, 5-yr",
    remeasure_yrs = 5,
    n_new_sites   = 30L,
    site_selector = function(sp) sp$site_id[sp$in_coarse_grid]
  )
)

# 4. Run power simulation
# Reporting scales are fully user-defined — no taxonomy assumed
results <- run_power_sim(
  ref_params       = ref,
  scenarios        = scenarios,
  plot_metadata    = site_meta,
  mode             = "conditional",
  data             = long_model,
  effect_sizes_pct = c(10, 20, 30),
  horizons         = c(10, 20),
  n_iter           = 200,
  alpha            = 0.10,
  reporting_groups = list("Region" = "region_col",
                          "Subregion" = "subregion_col"),
  place_var        = "site_id"
)

# 5. Summarise and plot
mdc <- compute_mdc(results, power_target = 0.80)

plot_power(results, scale_filter = "Overall")
plot_mdc(mdc,    scale_filter = "Overall")
plot_power_gain(results, scale_filter = "Overall")  # hybrid scenarios only

# Precision (CV of state estimates)
precision <- run_precision(
  scenarios        = scenarios,
  data             = long_model,
  plot_metadata    = site_meta,
  reporting_groups = list("Region" = "region_col"),
  place_var        = "site_id"
)
plot_cv(precision, scale_filter = "Overall")
```

---

## Key functions

| Function | Description |
|---|---|
| `extract_params()` | Extract model parameters — dispatches on model class |
| `scenario()` | Define a design scenario |
| `run_power_sim()` | Run power simulations across scenarios and reporting scales |
| `compute_mdc()` | Derive minimum detectable change from power results |
| `run_precision()` | Bootstrap CV of state estimates across scenarios |
| `plot_power()` | Power curve bar chart by scenario |
| `plot_power_gain()` | Hybrid scenario power gain (combined vs legacy-only) |
| `plot_mdc()` | MDC heatmap by scenario and horizon |
| `plot_cv()` | CV bar chart by scenario |
| `init_conditional()` | Initialise legacy plots from observed record |
| `init_new_sites()` | Initialise new plots for hybrid scenarios |
| `init_prospective_marginal()` | Initialise plots from scratch (prospective mode) |
| `simulate_visits()` | Simulate future visits for a set of plots |
| `build_historical()` | Build historical data stub for legacy plots |
| `fit_and_test()` | Fit trend model and return p-value |
| `bootstrap_cv()` | Bootstrap CV for a set of sites |

---

## Supported model families

`extract_params()` dispatches on model class and returns a standardised
`monpwr_params` object. All downstream functions consume only this object.

| Model class | Families |
|---|---|
| `glmmTMB` | `truncated_nbinom2` (hurdle), `nbinom2`, `poisson`, `binomial` |
| `glmerMod` | `binomial`, `poisson`, negative binomial via `glmer.nb` |
| `lmerMod` | Gaussian (for log-transformed indices) |

For any other model class, implement an `extract_params.myclass()` S3 method
returning a list with the `monpwr_params` structure. See the custom extractor
vignette for the required fields.

---

## Important notes

- `data` passed to `extract_params()` must be a **plain `data.frame`** — not
  a tibble or data.table. Call `as.data.frame(data)` first if needed. BLUP
  extraction via `ranef()` can fail silently on non-data-frame inputs.
- New sites in hybrid scenarios have BLUPs **resampled every replicate**. This
  correctly propagates uncertainty about new-site baselines into the power
  estimate. Fixing BLUPs across replicates would overstate power.
- The trend test model is intentionally simpler than the data-generating model
  (no fixed covariates beyond `visit_num`). This gives conservative power
  estimates consistent with typical monitoring reporting practice.
- CV estimates from `run_precision()` reflect legacy sites only. New sites in
  hybrid scenarios have no observed record and do not contribute.

---

## Package structure

```
R/
├── monpwr-package.R    # package documentation and namespace imports
├── extract_params.R    # parameter extraction — S3 dispatch on model class
├── initialisers.R      # plot-state initialisation functions
├── simulation.R        # simulate_visits(), build_historical(), fit_and_test()
├── run_power_sim.R     # scenario() constructor and main outer loop
├── summaries.R         # compute_mdc(), bootstrap_cv(), run_precision()
└── plots.R             # plot_power(), plot_mdc(), plot_cv(), plot_power_gain()
```

---

## License

MIT