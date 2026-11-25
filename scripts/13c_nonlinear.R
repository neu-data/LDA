# =============================================================================
# scripts/13c_nonlinear.R
# Module 13 (ADVANCED), section 3: nonlinear mixed-effects model.
#
# Purpose : fit an asymptotic exponential model to WALK SBP with nlme::nlme,
#             sbp_ij = base_i - A_i * (1 - exp(-month_ij / tau)) + e_ij,
#           with an arm-specific mean plateau change A and a common time
#           constant tau.
# Inputs  : data/walk_long.csv (via load_walk()), data/walk_truth.rds
# Outputs : outputs/tables/m13_nlme_*.csv, outputs/figures/m13_nlme_*.png,
#           outputs/models/m13_nlme*.rds
# Module  : 13 Advanced extensions
# =============================================================================

source(here::here("R", "setup.R"))
walk    <- load_walk()
walk_df <- as.data.frame(walk)          # nlme works with plain data frames
truth   <- readRDS(here("data", "walk_truth.rds"))

# ---- starting values --------------------------------------------------------
# base: the mean baseline SBP. A: the observed fall from month 0 to month 12
# in each arm. tau: most of the change is present by month 3, and tau is the
# time by which about 63% of the plateau change has happened, so 2 months is
# a sensible guess.
obs_means <- walk |>
  group_by(arm, month_sched) |>
  summarise(sbp = mean(sbp)) |>
  ungroup()
start_A <- obs_means |>
  group_by(arm) |>
  summarise(A = sbp[month_sched == 0] - sbp[month_sched == 12]) |>
  pull(A)
start_vals <- c(base = mean(walk$sbp[walk$visit == 1]),
                A = round(start_A, 1), tau = 2)
start_vals

# ---- model 1: random baseline only -----------------------------------------
# A ~ arm - 1 gives one plateau parameter per arm (cell-means coding), which
# reads more easily than a reference arm plus differences.
m_nl1 <- nlme(sbp ~ base - A * (1 - exp(-month / tau)),
              data   = walk_df,
              fixed  = list(base ~ 1, A ~ arm - 1, tau ~ 1),
              random = list(id = pdSymm(base ~ 1)),
              start  = start_vals)
summary(m_nl1)

# ---- model 2: random baseline and random plateau ---------------------------
m_nl2 <- nlme(sbp ~ base - A * (1 - exp(-month / tau)),
              data   = walk_df,
              fixed  = list(base ~ 1, A ~ arm - 1, tau ~ 1),
              random = list(id = pdDiag(base + A ~ 1)),
              start  = start_vals)
summary(m_nl2)
nl_lrt <- anova(m_nl1, m_nl2)  # ML fits; the null is on the boundary
nl_lrt
save_model(m_nl2, "m13_nlme_asymptotic")

ci_nl2 <- intervals(m_nl2, which = "fixed")$fixed
nl_tab <- data.frame(parameter = rownames(ci_nl2), ci_nl2, row.names = NULL)
# Truth: A excludes the part of the programme effect that acts through steps.
# The total plateau change adds the steps-mediated part (-0.6 x increase in
# steps), which a model without steps estimates.
nl_tab$truth_direct <- c(NA, truth$A, truth$tau)
nl_tab$truth_total  <- c(NA, truth$A - truth$beta_steps_within *
                           truth$steps_inc, truth$tau)
nl_tab
save_table(nl_tab, "m13_nlme_fixed_vs_truth")

# Differences in plateau change vs usual care, as change in SBP
# (negative = larger fall than usual care); delta method for a linear
# combination is exact here
fe <- nlme::fixef(m_nl2)
V  <- vcov(m_nl2)
L  <- rbind(SP_vs_UC = c(0, -1, 1, 0, 0), EP_vs_UC = c(0, -1, 0, 1, 0))
d_est <- -as.vector(L %*% fe)
d_se  <- sqrt(diag(L %*% V %*% t(L)))
plateau_diff <- data.frame(contrast = rownames(L), estimate = d_est,
                           se = d_se, lower = d_est - 1.96 * d_se,
                           upper = d_est + 1.96 * d_se)
plateau_diff
save_table(plateau_diff, "m13_nlme_plateau_differences")

# ---- checks -----------------------------------------------------------------
res_df <- walk_df |>
  mutate(fitted = fitted(m_nl2),
         resid  = resid(m_nl2, type = "normalized"))
p_res <- ggplot(res_df, aes(fitted, resid)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = course_colours[["orange"]]) +
  labs(x = "Fitted SBP (mmHg)", y = "Normalised residual")
save_fig(p_res, "m13_nlme_residuals")

# Population curves against observed means
curve <- expand.grid(month = seq(0, 12, by = 0.1),
                     arm = factor(c("UC", "SP", "EP"),
                                  levels = c("UC", "SP", "EP"))) |>
  mutate(A   = fe[paste0("A.arm", arm)],
         sbp = fe["base"] - A * (1 - exp(-month / fe["tau"])))
p_curve <- ggplot() +
  geom_line(data = curve, aes(month, sbp, colour = arm), linewidth = 0.9) +
  geom_point(data = obs_means, aes(month_sched, sbp, colour = arm),
             size = 2) +
  scale_colour_arm(name = NULL, labels = walk_arm_labels) +
  scale_x_continuous(breaks = c(0, 1, 3, 6, 9, 12)) +
  labs(x = "Months since randomisation", y = "SBP (mmHg)")
save_fig(p_curve, "m13_nlme_population_curves")

message("13c done")
