# =============================================================================
# 01b_independence_simulation.R
# Module 1. Foundations of longitudinal data
#
# Purpose : Monte Carlo demonstration that treating repeated measurements as
#           independent (ordinary least squares, OLS) gives standard errors
#           that are too small for between-subject comparisons and too large
#           for within-subject comparisons. Each simulated trial is analysed
#           with OLS and with a linear mixed model, and both model-based SEs
#           are compared with the true sampling SD of the estimates across
#           replicates.
# Inputs  : none (data simulated with sim_two_arm() from R/sim_helpers.R)
# Outputs : outputs/tables/m01_sim_summary.csv
#           outputs/models/m01_sim_replicates.rds
#           outputs/figures/m01_sim_sampling.png
# Notes   : The same functions, seeds and numbers of replicates are used in
#           notes/01-foundations.qmd, so the results agree exactly.
#           Run time is a few minutes. ALL DATA ARE SIMULATED.
# =============================================================================

source(here::here("R", "setup.R"))

# ---- simulation functions -----------------------------------------------------
# Scenario A: random intercept only (sd_slope = 0), no true arm difference and
# a common linear time trend. The arm effect is a between-subject comparison;
# the time slope is a within-subject comparison.
sim_once_ri <- function() {
  d <- sim_two_arm(n_per_arm = 50, times = c(0, 3, 6, 9, 12),
                   b0 = 150, b_time = -0.5, delta = 0,
                   sd_int = 11, sd_slope = 0, cor_int_slope = 0, sd_resid = 6)
  ols <- coef(summary(lm(y ~ time + trt, data = d)))
  lmm <- coef(summary(lmer(y ~ time + trt + (1 | id), data = d,
                           control = lmer_ctrl)))
  c(arm_est = ols["trt", 1], arm_se_ols = ols["trt", 2],
    arm_se_lmm = lmm["trt", 2], arm_est_lmm = lmm["trt", 1],
    time_est = ols["time", 1], time_se_ols = ols["time", 2],
    time_se_lmm = lmm["time", 2], time_est_lmm = lmm["time", 1])
}

# Scenario B: the sim_two_arm() defaults (random intercept and random slope,
# SD of slopes 0.4 per month) with a true difference in mean slope of -0.5.
# The target is the arm-by-time interaction, a comparison of within-subject
# change between groups.
sim_once_rs <- function(sd_slope = 0.4) {
  d <- sim_two_arm(n_per_arm = 50, sd_slope = sd_slope)
  ols <- coef(summary(lm(y ~ time * trt, data = d)))
  lmm <- coef(summary(suppressMessages(
    lmer(y ~ time * trt + (time | id), data = d, control = lmer_ctrl))))
  c(int_est = ols["time:trt", 1], int_se_ols = ols["time:trt", 2],
    int_se_lmm = lmm["time:trt", 2], int_est_lmm = lmm["time:trt", 1])
}

# Summarise one estimand: true sampling SD of the estimates across replicates
# with its Monte Carlo SE (Morris et al. 2019: SD / sqrt(2(R - 1))), the mean
# model-based SEs, and the coverage of nominal 95% Wald intervals.
summarise_se <- function(est, se_ols, est_lmm, se_lmm, truth, label) {
  R <- length(est)
  cov_ols <- mean(abs(est - truth) <= 1.96 * se_ols)
  cov_lmm <- mean(abs(est_lmm - truth) <= 1.96 * se_lmm)
  data.frame(
    estimand      = label,
    replicates    = R,
    true_value    = truth,
    mean_estimate = mean(est),
    empirical_sd  = sd(est),
    mcse_sd       = sd(est) / sqrt(2 * (R - 1)),
    mean_se_ols   = mean(se_ols),
    mean_se_lmm   = mean(se_lmm),
    ratio_ols     = mean(se_ols) / sd(est),
    ratio_lmm     = mean(se_lmm) / sd(est_lmm),
    coverage_ols  = cov_ols,
    coverage_lmm  = cov_lmm,
    mcse_coverage_ols = sqrt(cov_ols * (1 - cov_ols) / R),
    mcse_coverage_lmm = sqrt(cov_lmm * (1 - cov_lmm) / R)
  )
}

# ---- run scenario A (1000 replicates) ------------------------------------------
set.seed(20261)
res_a <- as.data.frame(t(replicate(1000, sim_once_ri())))

# ---- run scenario B (500 replicates; random-slope fits are slower) ---------------
set.seed(20262)
res_b <- as.data.frame(t(replicate(500, sim_once_rs())))

# ---- summarise -------------------------------------------------------------------
sim_summary <- rbind(
  summarise_se(res_a$arm_est, res_a$arm_se_ols, res_a$arm_est_lmm,
               res_a$arm_se_lmm, 0, "A: arm difference (between-subject)"),
  summarise_se(res_a$time_est, res_a$time_se_ols, res_a$time_est_lmm,
               res_a$time_se_lmm, -0.5, "A: time slope (within-subject)"),
  summarise_se(res_b$int_est, res_b$int_se_ols, res_b$int_est_lmm,
               res_b$int_se_lmm, -0.5, "B: arm-by-time interaction")
)
print(sim_summary, digits = 3)
save_table(sim_summary, "m01_sim_summary")
save_model(list(scenario_a = res_a, scenario_b = res_b), "m01_sim_replicates")

# ---- figure: sampling distributions with naive and mixed-model SEs -----------------
plot_df <- rbind(
  data.frame(estimand = "Arm difference (between-subject)",
             est = res_a$arm_est, truth = 0,
             se_ols = mean(res_a$arm_se_ols), se_lmm = mean(res_a$arm_se_lmm)),
  data.frame(estimand = "Time slope (within-subject)",
             est = res_a$time_est, truth = -0.5,
             se_ols = mean(res_a$time_se_ols), se_lmm = mean(res_a$time_se_lmm))
)
curves <- plot_df |>
  group_by(estimand) |>
  summarise(truth = first(truth), se_ols = first(se_ols),
            se_lmm = first(se_lmm), lo = min(est), hi = max(est)) |>
  rowwise() |>
  reframe(estimand = estimand,
          x = seq(truth - 4.5 * max(se_ols, se_lmm),
                  truth + 4.5 * max(se_ols, se_lmm), length.out = 200),
          `Naive OLS SE` = dnorm(x, truth, se_ols),
          `Mixed-model SE` = dnorm(x, truth, se_lmm)) |>
  pivot_longer(c(`Naive OLS SE`, `Mixed-model SE`), names_to = "curve",
               values_to = "density")

p_sim <- ggplot(plot_df, aes(est)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "#D9D9D9", colour = "white") +
  geom_line(data = curves, aes(x, density, colour = curve), linewidth = 0.8) +
  facet_wrap(~ estimand, scales = "free") +
  scale_colour_manual(values = c(`Naive OLS SE` = "#D55E00",
                                 `Mixed-model SE` = "#0072B2"), name = NULL) +
  labs(x = "Estimate across 1,000 simulated trials", y = "Density")
save_fig(p_sim, "m01_sim_sampling", width = 9, height = 4)
