# =============================================================================
# scripts/13e_msm.R
# Module 13 (ADVANCED), section 5: time-varying confounding and a marginal
# structural model (MSM) fitted by stabilised inverse probability weighting.
#
# Purpose : a small, self-contained simulation with two treatment times.
#           L1 is a time-varying confounder that is itself affected by the
#           earlier treatment A0. Standard regression adjustment for L1 is
#           biased; a stabilised-IPW MSM recovers the true effects.
# Inputs  : none (simulated here, fixed seed)
# Outputs : outputs/tables/m13_msm_*.csv, outputs/figures/m13_msm_*.png
# Module  : 13 Advanced extensions
# =============================================================================

source(here::here("R", "setup.R"))

# ---- data-generating mechanism ---------------------------------------------
# L0: baseline severity. U: unmeasured frailty (affects L1 and Y only).
# A0: treatment at time 0 (more likely if L0 is high).
# L1: severity at time 1, LOWERED by A0 and raised by U.
# A1: treatment at time 1 (more likely if L1 is high or A0 = 1).
# Y : final outcome (higher = worse), lowered by A0 and A1 and raised by
#     L0, L1 and U.
sim_tvc <- function(n, seed) {
  set.seed(seed)
  L0 <- rnorm(n)
  U  <- rnorm(n)
  A0 <- rbinom(n, 1, plogis(-0.3 + 0.8 * L0))
  L1 <- 0.5 * L0 - 0.8 * A0 + 0.8 * U + rnorm(n, 0, 0.6)
  A1 <- rbinom(n, 1, plogis(-0.5 + 1.0 * L1 + 1.0 * A0))
  Y  <- 50 + 3 * L0 + 1.5 * L1 - 2 * A0 - 2 * A1 + 4 * U + rnorm(n, 0, 3)
  data.frame(id = seq_len(n), L0, A0, L1, A1, Y)
}

# True MSM: E[Y^(a0, a1)] = beta0 + beta1 a0 + beta2 a1.
# a0 acts directly (-2) and through L1 (-0.8 x 1.5 = -1.2): beta1 = -3.2.
# a1 acts directly only: beta2 = -2.
msm_truth <- c(A0 = -2 - 0.8 * 1.5, A1 = -2)

dat <- sim_tvc(n = 5000, seed = 2026)

# ---- standard regression approaches -----------------------------------------
fit_unadj <- lm(Y ~ A0 + A1, data = dat)
fit_L0    <- lm(Y ~ A0 + A1 + L0, data = dat)
fit_L01   <- lm(Y ~ A0 + A1 + L0 + L1, data = dat)

# ---- stabilised inverse probability weights, by hand -------------------------
# Denominator: probability of the treatment actually received given the
# past (confounder and treatment history). Numerator: given treatment
# history only. Their ratio, multiplied over time, is the stabilised weight.
ps0_den <- glm(A0 ~ L0,            family = binomial, data = dat)
ps0_num <- glm(A0 ~ 1,             family = binomial, data = dat)
ps1_den <- glm(A1 ~ A0 + L0 + L1,  family = binomial, data = dat)
ps1_num <- glm(A1 ~ A0,            family = binomial, data = dat)

p_received <- function(model, a) {
  p <- fitted(model)
  ifelse(a == 1, p, 1 - p)
}
dat <- dat |>
  mutate(sw0 = p_received(ps0_num, A0) / p_received(ps0_den, A0),
         sw1 = p_received(ps1_num, A1) / p_received(ps1_den, A1),
         sw  = sw0 * sw1)
sw_summary <- data.frame(t(as.matrix(summary(dat$sw))), check.names = FALSE)
sw_summary
save_table(sw_summary, "m13_msm_weight_summary")

# Balance check: after weighting, L1 should no longer predict A1 given A0
bal <- rbind(
  unweighted = coef(summary(glm(A1 ~ A0 + L1, family = binomial,
                                data = dat)))["L1", 1:2],
  weighted   = coef(summary(glm(A1 ~ A0 + L1, family = quasibinomial,
                                data = dat, weights = sw)))["L1", 1:2])
bal

# ---- the MSM -----------------------------------------------------------------
# Weighted regression of Y on the treatment history only. geeglm with one
# row per cluster gives a robust (sandwich) standard error, which accounts
# for the estimated weights conservatively.
fit_msm <- geeglm(Y ~ A0 + A1, data = dat, id = id, weights = sw,
                  corstr = "independence")
summary(fit_msm)

grab <- function(m, label, robust = FALSE) {
  s <- coef(summary(m))
  data.frame(method = label, term = c("A0", "A1"),
             estimate = s[c("A0", "A1"), 1], se = s[c("A0", "A1"), 2],
             truth = msm_truth, row.names = NULL)
}
msm_tab <- bind_rows(
  grab(fit_unadj, "Unadjusted regression"),
  grab(fit_L0,    "Adjusted for L0"),
  grab(fit_L01,   "Adjusted for L0 and L1"),
  grab(fit_msm,   "MSM, stabilised IPW")
) |>
  mutate(lower = estimate - 1.96 * se, upper = estimate + 1.96 * se,
         bias = estimate - truth)
msm_tab
save_table(msm_tab, "m13_msm_estimates")

p_msm <- ggplot(msm_tab, aes(estimate, method)) +
  geom_vline(aes(xintercept = truth), linetype = 2, colour = "#6B6B6B") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  facet_wrap(~ term, labeller = labeller(term = c(
    A0 = "Effect of treatment at time 0 (A0)",
    A1 = "Effect of treatment at time 1 (A1)"))) +
  labs(x = "Estimated effect on Y (dashed line: true MSM value)", y = NULL)
save_fig(p_msm, "m13_msm_estimates", width = 8, height = 3.6)

# ---- repeated simulation: is the difference systematic? ----------------------
# 200 datasets of 1,000 participants each. Bias = mean estimate - truth.
rep_est <- bind_rows(lapply(1:200, function(r) {
  d <- sim_tvc(1000, seed = r)
  w <- p_received(glm(A0 ~ 1, binomial, d), d$A0) /
       p_received(glm(A0 ~ L0, binomial, d), d$A0) *
       p_received(glm(A1 ~ A0, binomial, d), d$A1) /
       p_received(glm(A1 ~ A0 + L0 + L1, binomial, d), d$A1)
  data.frame(rep = r,
             method = rep(c("Adjusted for L0 and L1", "MSM, stabilised IPW"),
                          each = 2),
             term = c("A0", "A1"),
             estimate = c(coef(lm(Y ~ A0 + A1 + L0 + L1, d))[c("A0", "A1")],
                          coef(lm(Y ~ A0 + A1, d, weights = w))[c("A0", "A1")]))
}))
rep_summary <- rep_est |>
  group_by(method, term) |>
  summarise(mean_estimate = mean(estimate), empirical_se = sd(estimate),
            mcse_mean = sd(estimate) / sqrt(n())) |>
  ungroup() |>
  mutate(truth = msm_truth[term], bias = mean_estimate - truth)
rep_summary
save_table(rep_summary, "m13_msm_repeated_simulation")

message("13e done")
