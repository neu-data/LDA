# =============================================================================
# 06_interpretation_prediction.R
# Module 6. Interpretation, contrasts and prediction
#
# Purpose : From the primary WALK model, compute estimated marginal means,
#           planned contrasts (difference in mean change from baseline to
#           month 12), multiplicity adjustments, population-average and
#           subject-specific predictions, confidence and prediction intervals,
#           a results table and a publication-style figure. Also fits the
#           constrained-baseline (cLDA) parameterisation.
# Inputs  : data/walk_long.csv (via load_walk())
# Outputs : outputs/figures/m06_*.png, outputs/tables/m06_*.csv,
#           outputs/models/m06_*.rds
# Module  : 06-interpretation-prediction
#
# SIMULATED DATA: the WALK study does not exist.
# =============================================================================

source(here::here("R", "setup.R"))
walk  <- load_walk()
truth <- readRDS(here("data", "walk_truth.rds"))
sched <- c(0, 1, 3, 6, 9, 12)

# ---- primary model -----------------------------------------------------------
m_prim <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk,
               control = lmer_ctrl)
save_model(m_prim, "m06_primary_lmer")

# ---- main effects, interactions and centring ----------------------------------
# In the linear-time model the arm "main effects" are differences at month 0.
# Centring time at month 12 turns them into differences at month 12; the fit
# is identical.
walk <- walk |> mutate(month12 = month - 12,
                       age10 = (age - 60) / 10, bmi_c = bmi - 29)
m_lin   <- lmer(sbp ~ month * arm + (1 + month | id), data = walk,
                control = lmer_ctrl)
m_lin12 <- lmer(sbp ~ month12 * arm + (1 + month12 | id), data = walk,
                control = lmer_ctrl)
centring <- data.frame(term = names(fixef(m_lin)),
                       time_from_0 = fixef(m_lin),
                       time_from_12 = fixef(m_lin12), row.names = NULL)
print(centring, digits = 4)
save_table(centring, "m06_centring_linear_time")

# ---- estimated marginal means -------------------------------------------------
emm <- emmeans(m_prim, ~ visit_f * arm)
emm_tab <- as.data.frame(summary(emm)) |>
  mutate(month = as.numeric(as.character(visit_f)))
save_table(emm_tab, "m06_emmeans")

# ---- planned contrasts --------------------------------------------------------
# Coefficient vectors over the 18 cells (visit varies fastest within arm)
cells  <- expand.grid(visit_f = levels(walk$visit_f), arm = levels(walk$arm))
cell_w <- function(a, v) as.numeric(cells$arm == a & cells$visit_f == v)
chg    <- function(a, v = "12") cell_w(a, v) - cell_w(a, "0")
L <- list("SP - UC" = chg("SP") - chg("UC"),
          "EP - UC" = chg("EP") - chg("UC"),
          "EP - SP" = chg("EP") - chg("SP"))

change_within <- summary(contrast(emm, list(UC = chg("UC"), SP = chg("SP"),
                                            EP = chg("EP"))), infer = TRUE)
dic_unadj <- summary(contrast(emm, L), infer = TRUE)
dic_holm  <- summary(contrast(emm, L[1:2], adjust = "holm"), infer = TRUE)
dic_dunn  <- summary(contrast(emm, L[1:2], adjust = "mvt"), infer = TRUE)
print(change_within); print(dic_unadj); print(dic_holm); print(dic_dunn)

# interaction-contrast approach: every visit vs baseline, every arm vs UC
ic <- summary(contrast(emm, interaction = c(visit_f = "trt.vs.ctrl1",
                                           arm = "trt.vs.ctrl1")),
              infer = TRUE)
print(ic)

# ---- results table ---------------------------------------------------------------
tc <- truth$true_change
results <- data.frame(
  quantity = c(paste("Mean change,", walk_arm_labels),
               "Difference SP - UC", "Difference EP - UC", "Difference EP - SP"),
  estimate = c(change_within$estimate, dic_unadj$estimate),
  lower    = c(change_within$lower.CL, dic_unadj$lower.CL),
  upper    = c(change_within$upper.CL, dic_unadj$upper.CL),
  p        = c(rep(NA, 3), dic_unadj$p.value)
) |>
  mutate(estimate_ci = fmt_ci(estimate, lower, upper),
         p_text = ifelse(is.na(p), "", fmt_p(p)))
print(results)
save_table(results, "m06_results_table")

# ---- population-average vs subject-specific predictions ---------------------------
walk$pred_pa <- predict(m_prim, re.form = NA)   # fixed effects only
walk$pred_ss <- predict(m_prim)                  # fixed effects + BLUPs
# On a complete grid (every participant at every scheduled month) the average of
# the subject-specific predictions equals the population-average prediction
# up to the (near-zero) average of the BLUPs: identity link.
grid <- walk |> distinct(id, arm) |>
  tidyr::crossing(visit_f = factor(sched, levels = sched)) |>
  mutate(month = as.numeric(as.character(visit_f)))
grid$pred_ss <- predict(m_prim, newdata = grid)
grid$pred_pa <- predict(m_prim, newdata = grid, re.form = NA)
pa_vs_ss <- grid |> group_by(arm, visit_f) |>
  summarise(mean_ss = mean(pred_ss), pa = mean(pred_pa)) |>
  mutate(diff = mean_ss - pa)
print(pa_vs_ss)
save_table(pa_vs_ss, "m06_pa_vs_ss_grid")

# ---- confidence vs prediction intervals ----------------------------------------
G  <- as.matrix(VarCorr(m_prim)$id)
s2 <- sigma(m_prim)^2
pi_tab <- emm_tab |>
  mutate(var_re   = G[1, 1] + 2 * G[1, 2] * month + G[2, 2] * month^2,
         se_newpt = sqrt(SE^2 + var_re + s2),       # new participant, one obs
         pi_lower = emmean - qnorm(0.975) * se_newpt,
         pi_upper = emmean + qnorm(0.975) * se_newpt)
save_table(pi_tab, "m06_ci_pi")

# simulation check of the prediction interval for a new participant
set.seed(606)
n_sim <- 20000
beta  <- MASS::mvrnorm(n_sim, fixef(m_prim), as.matrix(vcov(m_prim)))
b     <- MASS::mvrnorm(n_sim, c(0, 0), G)
# EP arm at month 12: intercept + visit_f12 + armEP + visit_f12:armEP
x_ep12 <- as.numeric(names(fixef(m_prim)) %in%
                       c("(Intercept)", "visit_f12", "armEP", "visit_f12:armEP"))
y_new <- as.vector(beta %*% x_ep12) + b[, 1] + 12 * b[, 2] +
  rnorm(n_sim, 0, sigma(m_prim))
print(quantile(y_new, c(0.025, 0.975)))
print(filter(pi_tab, arm == "EP", month == 12) |> select(pi_lower, pi_upper))

# ---- presentation figure ---------------------------------------------------------
obs_means <- walk |> group_by(arm, month_sched) |>
  summarise(mean = mean(sbp), n = n())
p_traj <- ggplot(emm_tab, aes(month, emmean, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.15,
              colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(data = obs_means, aes(month_sched, mean), shape = 1, size = 2) +
  scale_colour_arm(labels = walk_arm_labels, name = NULL) +
  scale_fill_arm(labels = walk_arm_labels, name = NULL) +
  scale_x_continuous(breaks = sched) +
  labs(x = "Months since randomisation",
       y = "Systolic blood pressure (mmHg)",
       caption = "Lines and bands: model-estimated means with 95% CIs. Open circles: observed means of attended visits.")
save_fig(p_traj, "m06_fitted_trajectories")

dic_time <- ic |> as.data.frame() |>
  mutate(month = as.numeric(sub(" - 0", "", visit_f_trt.vs.ctrl1)),
         arm = sub(" - UC", "", arm_trt.vs.ctrl1)) |>
  bind_rows(data.frame(month = 0, arm = c("SP", "EP"), estimate = 0,
                       lower.CL = 0, upper.CL = 0)) |>
  mutate(arm = factor(arm, levels = c("SP", "EP")))
p_dic <- ggplot(dic_time, aes(month, estimate, colour = arm)) +
  geom_hline(yintercept = 0, colour = "#9A9A9A") +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),
                  position = position_dodge(width = 0.5), size = 0.3) +
  geom_line(position = position_dodge(width = 0.5)) +
  scale_colour_arm(labels = walk_arm_labels, name = NULL) +
  scale_x_continuous(breaks = sched) +
  labs(x = "Months since randomisation",
       y = "Difference in mean change\nvs usual care (mmHg)")
save_fig(p_dic, "m06_difference_in_change")

p_pi <- ggplot(filter(pi_tab, arm == "EP"), aes(month, emmean)) +
  geom_ribbon(aes(ymin = pi_lower, ymax = pi_upper), fill = "#D55E00",
              alpha = 0.12) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), fill = "#D55E00",
              alpha = 0.35) +
  geom_line(colour = "#D55E00", linewidth = 0.8) +
  scale_x_continuous(breaks = sched) +
  labs(x = "Months since randomisation",
       y = "Systolic blood pressure (mmHg)",
       caption = "Enhanced programme. Dark band: 95% CI for the mean. Light band: 95% prediction interval for one measurement of a new participant.")
save_fig(p_pi, "m06_ci_vs_pi")

# ---- constrained baseline (cLDA) ---------------------------------------------------
# One common mean at baseline for all arms (randomisation makes the true
# baseline means equal); arm-specific means only after randomisation.
walk <- walk |>
  mutate(grp = ifelse(visit == 1, "base", paste(arm, visit_f, sep = "_")),
         grp = relevel(factor(grp), ref = "base"))
m_clda <- lmer(sbp ~ grp + (1 + month | id), data = walk, control = lmer_ctrl)
save_model(m_clda, "m06_clda_lmer")
e_clda <- emmeans(m_clda, ~ grp)
lv <- levels(walk$grp)
gw <- function(g) as.numeric(lv == g)
L_clda <- list("SP - UC" = gw("SP_12") - gw("UC_12"),
               "EP - UC" = gw("EP_12") - gw("UC_12"),
               "EP - SP" = gw("EP_12") - gw("SP_12"))
clda <- summary(contrast(e_clda, L_clda), infer = TRUE)
print(clda)
clda_vs_lda <- data.frame(contrast = names(L),
                          lda_est = dic_unadj$estimate, lda_se = dic_unadj$SE,
                          clda_est = clda$estimate, clda_se = clda$SE)
print(clda_vs_lda)
save_table(clda_vs_lda, "m06_clda_vs_lda")

message("Module 6 script finished.")
