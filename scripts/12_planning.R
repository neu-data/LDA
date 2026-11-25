# =============================================================================
# 12_planning.R
# Module 12. Study planning, reporting and reproducibility
#
# Purpose
#   Design calculations for longitudinal studies: follow-up schedules, the
#   variance of a slope estimate in a balanced design, the effect of
#   correlation and dropout on precision, and simulation-based power with
#   Monte Carlo standard errors. Also records the computational environment.
#
# Inputs
#   data/walk_truth.rds (true WALK mean trajectories, for the schedule figure)
#   R/sim_helpers.R      sim_two_arm(), power_two_arm()
#
# Outputs
#   outputs/figures/m12_*.png, outputs/tables/m12_*.csv
#   outputs/models/m12_power_grid.rds      cached power simulations
#   outputs/models/m12_slope_sims.rds      cached estimates for the formula check
#   outputs/m12_session_info.txt           sessionInfo() at the time of running
#
# Run time: the power grid needs 3,000 model fits and takes roughly 8 to 12
# minutes on a laptop. Results are cached in outputs/models/, so the notes,
# workbook and answer key render quickly. Delete the .rds files to recompute
# them (the seeds are fixed, so the results are identical).
#
# All data are simulated. Results illustrate methods, not empirical findings.
# =============================================================================

source(here::here("R", "setup.R"))

cache_or_run <- function(name, expr) {
  f <- here("outputs", "models", paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  x <- expr
  save_model(x, name)
  x
}

# Design parameters: the defaults of sim_two_arm(), used throughout
pars <- list(delta = -0.5, sd_int = 11, sd_slope = 0.4, cor_int_slope = -0.2,
             sd_resid = 6, times = c(0, 3, 6, 9, 12))

# ---- 1. follow-up schedules and the shape of change ------------------------
truth <- readRDS(here("data", "walk_truth.rds"))
curve <- expand.grid(month = seq(0, 12, by = 0.1), arm = c("UC", "SP", "EP")) |>
  mutate(change = -truth$A[as.character(arm)] * (1 - exp(-month / truth$tau)),
         arm = factor(arm, levels = c("UC", "SP", "EP")))
schedules <- list(
  "WALK: 0, 1, 3, 6, 9, 12" = c(0, 1, 3, 6, 9, 12),
  "Equally spaced: 0, 3, 6, 9, 12" = c(0, 3, 6, 9, 12),
  "Sparse: 0, 6, 12" = c(0, 6, 12))
sched_df <- bind_rows(lapply(names(schedules), function(s)
  data.frame(schedule = s, month = schedules[[s]])))
sched_df$schedule <- factor(sched_df$schedule, levels = names(schedules))
sched_pts <- sched_df |>
  inner_join(curve |> mutate(month = round(month, 1)), by = "month",
             relationship = "many-to-many")
p_sched <- ggplot(curve, aes(month, change, colour = arm)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = sched_pts, size = 2) +
  geom_line(data = sched_pts, linetype = 2, linewidth = 0.4) +
  facet_wrap(~schedule, ncol = 3) +
  scale_colour_arm(name = NULL, labels = walk_arm_labels) +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12)) +
  labs(x = "Months since randomisation",
       y = "True mean change in SBP (mmHg)")
save_fig(p_sched, "m12_schedules", width = 9, height = 3.8)

# ---- 2. variance of a slope estimate (balanced design) ---------------------
# Per-participant OLS slope variance: tau1^2 + sigma^2 / S_tt, where
# S_tt = sum_j (t_j - tbar)^2. The difference in mean slopes between two arms
# of n participants each has variance 2 * (tau1^2 + sigma^2 / S_tt) / n.
s_tt <- function(times) sum((times - mean(times))^2)
var_slope_diff <- function(n_per_arm, times, sd_slope, sd_resid) {
  2 * (sd_slope^2 + sd_resid^2 / s_tt(times)) / n_per_arm
}
power_formula <- function(n_per_arm, times = pars$times, delta = pars$delta,
                          sd_slope = pars$sd_slope, sd_resid = pars$sd_resid,
                          alpha = 0.05) {
  se <- sqrt(var_slope_diff(n_per_arm, times, sd_slope, sd_resid))
  pnorm(abs(delta) / se - qnorm(1 - alpha / 2))
}
n_formula <- function(power = 0.9, times = pars$times, delta = pars$delta,
                      sd_slope = pars$sd_slope, sd_resid = pars$sd_resid,
                      alpha = 0.05) {
  ceiling(2 * (qnorm(1 - alpha / 2) + qnorm(power))^2 *
            (sd_slope^2 + sd_resid^2 / s_tt(times)) / delta^2)
}

design_tab <- data.frame(
  schedule = c("0, 12", "0, 6, 12", "0, 3, 6, 9, 12", "0, 1, 3, 6, 9, 12",
               "0, 2, 4, 6, 8, 10, 12", "0, 0.5, 1, ..., 12 (25 visits)",
               "0, 1, 2, 10, 11, 12"),
  stringsAsFactors = FALSE)
design_times <- list(c(0, 12), c(0, 6, 12), c(0, 3, 6, 9, 12),
                     c(0, 1, 3, 6, 9, 12), seq(0, 12, 2), seq(0, 12, 0.5),
                     c(0, 1, 2, 10, 11, 12))
design_tab$n_visits <- lengths(design_times)
design_tab$S_tt <- vapply(design_times, s_tt, 0)
design_tab$var_subject_slope <- pars$sd_slope^2 + pars$sd_resid^2 /
  design_tab$S_tt
design_tab$se_diff_n50 <- sqrt(2 * design_tab$var_subject_slope / 50)
design_tab$n_per_arm_90 <- vapply(design_times, function(t)
  n_formula(0.9, times = t), 0)
save_table(design_tab, "m12_design_table")

# Check the formula against simulation: 500 trials, 50 per arm, no dropout
slope_sims <- cache_or_run("m12_slope_sims", {
  set.seed(1201)
  bind_rows(lapply(1:500, function(i) {
    d <- sim_two_arm(n_per_arm = 50)
    m <- suppressMessages(lmer(y ~ time * trt + (time | id), data = d,
                               control = lmerControl(calc.derivs = FALSE)))
    s <- coef(summary(m))["time:trt", ]
    data.frame(sim = i, est = s[["Estimate"]], se = s[["Std. Error"]],
               p = s[["Pr(>|t|)"]])
  }))
})
slope_check <- data.frame(
  quantity = c("Formula SE", "Empirical SD of estimates",
               "Mean model-based SE", "Mean estimate"),
  value = c(sqrt(var_slope_diff(50, pars$times, pars$sd_slope,
                                pars$sd_resid)),
            sd(slope_sims$est), mean(slope_sims$se), mean(slope_sims$est)))
save_table(slope_check, "m12_slope_formula_check")

# SE against the number of equally spaced visits over 12 months
visits_curve <- expand.grid(n_visits = 2:13, sd_slope = c(0.2, 0.4, 0.8)) |>
  mutate(S = vapply(n_visits, function(k) s_tt(seq(0, 12, length.out = k)), 0),
         se = sqrt(2 * (sd_slope^2 + pars$sd_resid^2 / S) / 50),
         sd_slope = factor(sd_slope))
p_visits <- ggplot(visits_curve, aes(n_visits, se, colour = sd_slope)) +
  geom_line() + geom_point(size = 1.5) +
  scale_colour_manual(values = c("0.2" = "#009E73", "0.4" = "#0072B2",
                                 "0.8" = "#D55E00"),
                      name = expression(tau[1] ~ "(mmHg/month)")) +
  scale_x_continuous(breaks = 2:13) +
  labs(x = "Equally spaced visits over 12 months",
       y = "SE of difference in slopes\n(mmHg/month, 50 per arm)")
save_fig(p_visits, "m12_se_by_visits")

# Correlation and a single follow-up: change score vs ANCOVA
rho_curve <- data.frame(rho = seq(0, 0.95, by = 0.01)) |>
  mutate(`Follow-up only` = 1, `Change score` = 2 * (1 - rho),
         ANCOVA = 1 - rho^2) |>
  pivot_longer(-rho, names_to = "analysis", values_to = "rel_var")
p_rho <- ggplot(rho_curve, aes(rho, rel_var, colour = analysis)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("Follow-up only" = "#6B6B6B",
                                 "Change score" = "#0072B2",
                                 ANCOVA = "#D55E00"), name = NULL) +
  labs(x = "Correlation between baseline and follow-up",
       y = "Variance relative to follow-up only")
save_fig(p_rho, "m12_correlation")

# Approximate variance with monotone dropout as generated by sim_two_arm():
# a participant is seen at the first k visits with probability
# p (1 - p)^(k - 1) for k < J and (1 - p)^(J - 1) for k = J. Each participant
# with k >= 2 visits contributes slope information 1 / (tau1^2 + sigma^2/S_k).
# This ignores the information that flows through the intercept-slope
# correlation, so it is an approximation.
power_formula_dropout <- function(n_per_arm, dropout, times = pars$times,
                                  delta = pars$delta, sd_slope = pars$sd_slope,
                                  sd_resid = pars$sd_resid, alpha = 0.05) {
  J <- length(times)
  pk <- if (dropout == 0) c(rep(0, J - 1), 1) else
    c(dropout * (1 - dropout)^(0:(J - 2)), (1 - dropout)^(J - 1))
  info <- vapply(seq_len(J), function(k) if (k < 2) 0 else
    1 / (sd_slope^2 + sd_resid^2 / s_tt(times[1:k])), 0)
  se <- sqrt(2 / (n_per_arm * sum(pk * info)))
  pnorm(abs(delta) / se - qnorm(1 - alpha / 2))
}

# ---- 3. simulation-based power ----------------------------------------------
grid <- expand.grid(n_per_arm = c(20, 35, 50, 65, 80),
                    dropout = c(0, 0.1, 0.2))
power_grid <- cache_or_run("m12_power_grid", {
  bind_rows(lapply(seq_len(nrow(grid)), function(k) {
    pw <- power_two_arm(n_sims = 200, seed = 1200 + k,
                        n_per_arm = grid$n_per_arm[k],
                        dropout = grid$dropout[k])
    data.frame(grid[k, ], power = pw[["power"]], mcse = pw[["mcse"]],
               n_sims = pw[["n_sims"]])
  }))
})
power_grid <- power_grid |>
  mutate(formula = mapply(power_formula_dropout, n_per_arm, dropout),
         lower = pmax(0, power - 1.96 * mcse),
         upper = pmin(1, power + 1.96 * mcse))
save_table(power_grid, "m12_power_grid")

formula_lines <- expand.grid(n_per_arm = 15:85, dropout = c(0, 0.1, 0.2)) |>
  mutate(formula = mapply(power_formula_dropout, n_per_arm, dropout))
p_power <- ggplot(power_grid,
                  aes(n_per_arm, power, colour = factor(dropout))) +
  geom_line(data = formula_lines, aes(y = formula), linetype = 2,
            linewidth = 0.5) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 2.5), size = 0.3) +
  geom_hline(yintercept = c(0.8, 0.9), colour = "grey70") +
  scale_colour_manual(values = c("0" = "#222222", "0.1" = "#0072B2",
                                 "0.2" = "#D55E00"),
                      name = "Dropout probability per visit") +
  labs(x = "Participants per arm",
       y = "Power (simulated, with 95% MC interval)")
save_fig(p_power, "m12_power_curve")

# Monte Carlo error as the number of simulated trials grows, using the 500
# stored trials (50 per arm, no dropout)
mc_running <- slope_sims |>
  mutate(reject = p < 0.05,
         n_sims = seq_along(reject),
         power = cumsum(reject) / n_sims,
         mcse = sqrt(power * (1 - power) / n_sims))
save_table(mc_running, "m12_mc_running")
p_mc <- ggplot(mc_running, aes(n_sims, power)) +
  geom_ribbon(aes(ymin = power - 1.96 * mcse, ymax = power + 1.96 * mcse),
              fill = "grey85") +
  geom_line() +
  geom_hline(yintercept = power_formula(50), linetype = 2,
             colour = "#D55E00") +
  coord_cartesian(ylim = c(0.6, 1)) +
  labs(x = "Number of simulated trials", y = "Running power estimate")
save_fig(p_mc, "m12_mc_running")

# ---- 4. computational environment -------------------------------------------
writeLines(capture.output(sessionInfo()),
           here("outputs", "m12_session_info.txt"))

message("Module 12 script finished.")
