# =============================================================================
# 09_time_varying_nonlinear.R
# Module 9. Time-varying predictors and nonlinear change.
#
# Purpose
#   (1) Separate within-person and between-person associations of daily
#       steps with SBP in WALK (person-mean centring, contextual effect),
#       and show how adjusting for steps changes the meaning of the arm
#       effect (steps is partly a mediator).
#   (2) Unequal measurement intervals: actual vs scheduled time, visit
#       number vs time, and the irregular KIDNEY visits.
#   (3) Polynomial, piecewise linear and natural spline trajectories, and
#       generalized additive mixed models (mgcv) for WALK and KIDNEY.
#
# Inputs   data/walk_long.csv, data/walk_truth.rds, data/kidney_long.csv
# Outputs  outputs/figures/m09_*.png, outputs/tables/m09_*.csv,
#          outputs/models/m09_*.rds
# =============================================================================

source(here::here("R", "setup.R"))

walk  <- load_walk()
truth <- readRDS(here("data", "walk_truth.rds"))

# ---- person-mean centring --------------------------------------------------
# steps_pm: the participant's mean steps over their observed visits
#           (between-person part). steps_gm: that mean centred at the grand
#           mean of the person means. steps_wc: the visit's deviation from the
#           participant's own mean (within-person part).
walk <- walk |>
  group_by(id) |>
  mutate(steps_pm = mean(steps)) |>
  ungroup() |>
  mutate(steps_wc = steps - steps_pm,
         steps_gm = steps_pm - mean(steps_pm[!duplicated(id)]),
         age10 = (age - 60) / 10,
         bmi_c = bmi - 29,
         month_p3 = pmax(month - 3, 0))

p_steps <- walk |>
  group_by(arm, visit_f) |>
  summarise(mean = mean(steps), se = sd(steps) / sqrt(n())) |>
  ggplot(aes(as.numeric(as.character(visit_f)), mean, colour = arm)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                width = 0.2) +
  scale_colour_arm(name = NULL, labels = walk_arm_labels) +
  labs(x = "Scheduled month", y = "Mean daily steps (thousands)")
save_fig(p_steps, "m09_steps_by_arm")

# ---- within and between associations ---------------------------------------
# All models use the primary mean structure (visit x arm) so that the steps
# coefficients are estimated net of the arm trajectories.
m_naive   <- lmer(sbp ~ steps + visit_f * arm + (1 + month | id),
                  data = walk, control = lmer_ctrl)
m_wb      <- lmer(sbp ~ steps_wc + steps_gm + visit_f * arm +
                    (1 + month | id), data = walk, control = lmer_ctrl)
m_wb_adj  <- lmer(sbp ~ steps_wc + steps_gm + visit_f * arm +
                    age10 + sex + bmi_c + (1 + month | id),
                  data = walk, control = lmer_ctrl)
# contextual-effect parameterisation: raw steps plus the person mean
m_ctx_adj <- lmer(sbp ~ steps + steps_gm + visit_f * arm +
                    age10 + sex + bmi_c + (1 + month | id),
                  data = walk, control = lmer_ctrl)

pick <- function(m, terms, model) {
  fixef_table(m) |> filter(term %in% terms) |> mutate(model = model)
}
tab_wb <- bind_rows(
  pick(m_naive, "steps", "Naive (single coefficient)"),
  pick(m_wb, c("steps_wc", "steps_gm"), "Within/between"),
  pick(m_wb_adj, c("steps_wc", "steps_gm"),
       "Within/between + age, sex, BMI"),
  pick(m_ctx_adj, c("steps", "steps_gm"),
       "Contextual parameterisation + age, sex, BMI")
)
save_table(tab_wb, "m09_within_between")
save_model(m_wb_adj, "m09_within_between_adj")

# ---- arm effects with and without steps -----------------------------------
m_primary <- lmer(sbp ~ visit_f * arm + (1 + month | id),
                  data = walk, control = lmer_ctrl)
m_steps   <- lmer(sbp ~ steps_wc + visit_f * arm + (1 + month | id),
                  data = walk, control = lmer_ctrl)

arm_contrast <- function(m, label) {
  e  <- emmeans(m, ~ visit_f | arm, at = list(visit_f = c("0", "12")))
  ch <- contrast(e, list(change = c(-1, 1)), by = "arm")
  summary(contrast(ch, "trt.vs.ctrl", by = NULL), infer = TRUE,
          adjust = "none") |>
    as.data.frame() |>
    mutate(model = label)
}
tab_med <- bind_rows(
  arm_contrast(m_primary, "Primary model (no steps)"),
  arm_contrast(m_steps, "Adjusted for within-person steps")
)
save_table(tab_med, "m09_arm_with_without_steps")

# ---- unequal intervals -------------------------------------------------------
m_actual <- lmer(sbp ~ month * arm + (1 + month | id), data = walk,
                 control = lmer_ctrl)
m_sched  <- lmer(sbp ~ month_sched * arm + (1 + month_sched | id),
                 data = walk, control = lmer_ctrl)
m_visit  <- lmer(sbp ~ I(visit - 1) * arm + (1 + I(visit - 1) | id),
                 data = walk, control = lmer_ctrl)
tab_time <- bind_rows(
  fixef_table(m_actual) |> mutate(model = "Actual month"),
  fixef_table(m_sched) |> mutate(model = "Scheduled month"),
  fixef_table(m_visit) |> mutate(model = "Visit number")
)
save_table(tab_time, "m09_time_scales")

kidney <- load_kidney() |>
  group_by(id) |> arrange(years, .by_group = TRUE) |>
  mutate(visit_no = row_number() - 1) |> ungroup() |>
  mutate(age10 = (age - 60) / 10)
mk_years <- lmer(egfr ~ years + (1 + years | id), data = kidney,
                 control = lmer_ctrl)
mk_visit <- lmer(egfr ~ visit_no + (1 + visit_no | id), data = kidney,
                 control = lmer_ctrl)
save_table(bind_rows(fixef_table(mk_years) |> mutate(model = "Years"),
                     fixef_table(mk_visit) |> mutate(model = "Visit number")),
           "m09_kidney_time_scales")

# ---- polynomial, piecewise and spline trajectories ------------------------
# ML (REML = FALSE) because the models differ in their fixed effects and are
# compared by AIC/BIC.
fit_ml <- function(f) lmer(f, data = walk, control = lmer_ctrl, REML = FALSE)
traj <- list(
  Linear      = fit_ml(sbp ~ month * arm + (1 + month | id)),
  Quadratic   = fit_ml(sbp ~ (month + I(month^2)) * arm + (1 + month | id)),
  Piecewise   = fit_ml(sbp ~ (month + month_p3) * arm + (1 + month | id)),
  `Natural spline` = fit_ml(sbp ~ splines::ns(month, df = 3) * arm +
                              (1 + month | id)),
  Categorical = fit_ml(sbp ~ visit_f * arm + (1 + month | id))
)
ns_knots <- attr(splines::ns(walk$month, df = 3), "knots")

nd <- expand.grid(month = seq(0, 12, by = 0.1),
                  arm = factor(c("UC", "SP", "EP"),
                               levels = c("UC", "SP", "EP"))) |>
  mutate(month_p3 = pmax(month - 3, 0))
curves <- bind_rows(lapply(names(traj)[1:4], function(nm) {
  p <- predict(traj[[nm]], newdata = nd, re.form = NA)
  nd |> mutate(fit = p, model = nm) |>
    group_by(arm) |> mutate(change = fit - fit[month == 0]) |> ungroup()
})) |>
  mutate(model = factor(model, levels = names(traj)[1:4]))
true_curve <- nd |>
  mutate(change = -truth$A[as.character(arm)] *
           (1 - exp(-month / truth$tau)) +
           truth$beta_steps_within * truth$steps_inc[as.character(arm)] *
           (1 - exp(-month / 1)))

p_traj <- ggplot(curves, aes(month, change, colour = arm)) +
  geom_line(data = true_curve, linetype = "dashed", colour = "black",
            aes(group = arm)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ model, nrow = 1) +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12)) +
  scale_colour_arm(name = NULL, labels = walk_arm_labels) +
  labs(x = "Month since randomisation",
       y = "Mean change in SBP (mmHg)",
       caption = "Dashed: true mean change (known only because the data are simulated)")
save_fig(p_traj, "m09_trajectories", width = 10, height = 4)

ic <- data.frame(model = names(traj),
                 n_par = sapply(traj, function(m) attr(logLik(m), "df")),
                 AIC = sapply(traj, AIC), BIC = sapply(traj, BIC))
m12 <- lapply(names(traj)[1:4], function(nm) {
  p <- predict(traj[[nm]], newdata = filter(nd, month %in% c(0, 12)),
               re.form = NA)
  ch <- p[c(2, 4, 6)] - p[c(1, 3, 5)]
  data.frame(model = nm, SP = ch[2] - ch[1], EP = ch[3] - ch[1])
}) |> bind_rows()
save_table(left_join(ic, m12, by = "model"), "m09_trajectory_comparison")

# ---- generalized additive mixed models --------------------------------------
# gamm(): random intercept and slope through nlme. Fast.
g_walk <- gamm(sbp ~ arm + s(month, by = arm, k = 5),
               random = list(id = ~ 1 + month), data = walk, method = "REML")
save_model(g_walk, "m09_gamm_walk")

# The same model written with random-effect smooths. gam() accepts it but
# takes several minutes here (each of the 720 random coefficients is a
# penalised term); bam() fits the identical formula in under a minute.
walk$idf <- factor(walk$id)
g_walk_bam <- bam(sbp ~ arm + s(month, by = arm, k = 5) +
                    s(idf, bs = "re") + s(idf, month, bs = "re"),
                  data = walk, method = "fREML")
save_model(g_walk_bam, "m09_bam_walk")

kidney$dm <- factor(kidney$diabetes, levels = 0:1, labels = c("No", "Yes"))
g_kid <- gamm(egfr ~ dm + age10 + s(years, by = dm, k = 6),
              random = list(id = ~ 1 + years), data = kidney,
              method = "REML")
save_model(g_kid, "m09_gamm_kidney")

ndk <- expand.grid(years = seq(0, 5, by = 0.05), dm = c("No", "Yes"),
                   age10 = 0)
pk <- predict(g_kid$gam, newdata = ndk, se.fit = TRUE)
ndk <- ndk |> mutate(fit = pk$fit, lo = fit - 1.96 * pk$se.fit,
                     hi = fit + 1.96 * pk$se.fit)
p_kid <- ggplot(ndk, aes(years, fit, colour = dm, fill = dm)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(No = "#0072B2", Yes = "#D55E00"),
                      name = "Diabetes") +
  scale_fill_manual(values = c(No = "#0072B2", Yes = "#D55E00"),
                    name = "Diabetes") +
  labs(x = "Years since entry", y = "Mean eGFR (mL/min/1.73m²), age 60")
save_fig(p_kid, "m09_gamm_kidney")

print(tab_wb); print(tab_med); print(left_join(ic, m12, by = "model"))
print(summary(g_walk$gam)); print(summary(g_kid$gam))
