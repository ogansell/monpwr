# =============================================================================
# EXPERIMENT 8: Family misspecification — NB2 truth, Poisson test
# =============================================================================
# Closes the "family-matching" caveat: all prior validation cells share the DGP
# family with the test model. Here the truth is overdispersed (NB2) but the
# test model is Poisson — the mistake a user makes when they don't notice
# overdispersion.
#
# Hypothesis: the Poisson test under-estimates SE(visit), over-rejects, and
# inflates power; inflation grows as overdispersion rises.
#
# NB2 (glmmTMB nbinom2): var = mu + mu^2/size.
# Large size (phi) -> near-Poisson; small size -> strong overdispersion.
# =============================================================================

library(lme4)
library(glmmTMB)
library(dplyr)
library(purrr)
library(ggplot2)

theme_set(theme_minimal(base_size = 12))

# --- Shared constants (match other experiments where applicable) --------------
n_plots    <- 30L
n_design   <- 6L          # visits per plot
intercept  <- 0.8         # log-scale intercept
true_trend <- log(1.05)   # 5% change per visit
sigma_plot <- 0.6         # plot-level RE SD
alpha_val  <- 0.05
n_sim      <- 100L

cat("\n\n=== EXPERIMENT 8: Family misspecification (NB2 truth, Poisson test) ===\n")
cat(sprintf("  n_plots = %d, n_visits = %d, intercept = %.2f, trend = %.4f, sigma = %.2f\n",
            n_plots, n_design, intercept, true_trend, sigma_plot))
cat(sprintf("  alpha = %.2f, n_sim = %d\n\n", alpha_val, n_sim))


# --- NB2 data generator -------------------------------------------------------
make_balanced_nb2 <- function(n_plots, n_visits, intercept, trend, sigma, phi) {
  dat <- expand.grid(
    plot  = factor(paste0("p", seq_len(n_plots))),
    visit = seq_len(n_visits)
  )
  re <- rnorm(n_plots, 0, sigma)
  dat$re <- re[as.integer(dat$plot)]
  dat$y  <- rnbinom(nrow(dat), size = phi,
                    mu = exp(intercept + trend * dat$visit + dat$re))
  dat
}


# --- Phi grid: small = overdispersed, large = near-Poisson -------------------
phi_grid <- c(0.5, 1, 2, 5, 20)

set.seed(42)

results_misspec <- map_dfr(phi_grid, function(phi) {
  cat(sprintf("--- size (phi) = %.1f ---\n", phi))

  # Compute the variance-to-mean ratio at the mean count for reference
  mean_mu <- exp(intercept + true_trend * mean(1:n_design))
  vmr <- 1 + mean_mu / phi   # var/mu = 1 + mu/size for NB2

  pvals <- replicate(n_sim, {
    dat <- make_balanced_nb2(n_plots, n_design, intercept, true_trend,
                             sigma_plot, phi)

    # NB2-correct arm
    p_nb2 <- tryCatch(
      suppressWarnings({
        f <- glmmTMB::glmmTMB(y ~ visit + (1 | plot),
                              family = glmmTMB::nbinom2, data = dat)
        summary(f)$coefficients$cond["visit", "Pr(>|z|)"]
      }), error = function(e) NA_real_)

    # Poisson-misspecified arm
    p_pois <- tryCatch(
      suppressWarnings({
        f <- lme4::glmer(y ~ visit + (1 | plot), family = poisson, data = dat)
        summary(f)$coefficients["visit", "Pr(>|z|)"]
      }), error = function(e) NA_real_)

    c(nb2 = p_nb2, pois = p_pois)
  })

  conv_nb2  <- sum(!is.na(pvals["nb2", ]))
  conv_pois <- sum(!is.na(pvals["pois", ]))
  power_nb2  <- mean(pvals["nb2", ]  < alpha_val, na.rm = TRUE)
  power_pois <- mean(pvals["pois", ] < alpha_val, na.rm = TRUE)

  cat(sprintf("  NB2-correct power: %.3f (%d/%d conv) | Poisson-misspec power: %.3f (%d/%d conv) | inflation: %+.3f\n",
              power_nb2, conv_nb2, n_sim, power_pois, conv_pois, n_sim,
              power_pois - power_nb2))

  data.frame(
    phi                   = phi,
    vmr_at_mean           = round(vmr, 2),
    power_nb2_correct     = round(power_nb2, 3),
    power_poisson_misspec = round(power_pois, 3),
    inflation             = round(power_pois - power_nb2, 3),
    conv_nb2              = conv_nb2,
    conv_pois             = conv_pois,
    stringsAsFactors      = FALSE
  )
})

cat("\n--- Results table ---\n")
print(results_misspec, row.names = FALSE)
cat(sprintf("\n  Max power inflation from Poisson misspecification: %+.3f\n",
            max(results_misspec$inflation, na.rm = TRUE)))


# --- Plot: power vs phi, NB2-correct vs Poisson-misspecified ---
plot_dat <- results_misspec |>
  tidyr::pivot_longer(
    cols = c(power_nb2_correct, power_poisson_misspec),
    names_to = "method",
    values_to = "power"
  ) |>
  mutate(
    method = dplyr::recode(method,
      "power_nb2_correct"     = "NB2 (correct)",
      "power_poisson_misspec" = "Poisson (misspecified)"
    )
  )

p <- ggplot(plot_dat, aes(x = phi, y = power, colour = method, shape = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_log10(breaks = phi_grid) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title    = "Experiment 8: Family misspecification (NB2 truth, Poisson test)",
    subtitle = sprintf("n=%d plots, %d visits, 5%% trend, sigma=%.1f, alpha=%.2f, %d reps",
                        n_plots, n_design, sigma_plot, alpha_val, n_sim),
    x        = "NB2 size (phi) — larger = less overdispersion",
    y        = "Power",
    colour   = "Test model",
    shape    = "Test model"
  ) +
  theme(legend.position = "bottom")

print(p)

ggsave("power_analysis_outputs/exp8_family_misspecification.png", p,
       width = 8, height = 5, dpi = 150)
cat("\n  Plot saved to power_analysis_outputs/exp8_family_misspecification.png\n")
