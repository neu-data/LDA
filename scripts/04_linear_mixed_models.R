# =============================================================================
# scripts/04_linear_mixed_models.R
# Module 4. Linear mixed-effects models
#
# Purpose: fit and interpret random-intercept and random-slope models for the
#          WALK trial; variance components, ICC and the implied marginal
#          covariance; coding of time (continuous, centred, categorical,
#          quadratic); the reference models and the primary estimand
#          (difference in mean change at month 12) with emmeans; ML versus
#          REML (including a small simulation); BLUPs and shrinkage; the same
#          model in nlme::lme; basic residual checks.
# Inputs:  data/walk_long.csv (via load_walk()), data/walk_truth.rds
# Outputs: outputs/figures/m04_*.png, outputs/tables/m04_*.csv,
#          outputs/models/m04_*.rds
# Module:  04-linear-mixed-models
# =============================================================================

source(here::here("R", "setup.R"))
walk  <- load_walk()
truth <- readRDS(here("data", "walk_truth.rds"))
sched <- c(0, 1, 3, 6, 9, 12)

# ---- 1. Individual trajectories ---------------------------------------------
# A mixed model starts from the idea that every participant has a line of
# their own. Plot a random subset so the lines can be seen.
p_spag <- plot_spaghetti(walk, month, sbp, id, group = arm, n_show = 60) +
  labs(x = "Months since randomisation", y = "SBP (mmHg)")
save_fig(p_spag, "m04_spaghetti", width = 8, height = 4)

# ---- 2. Random-intercept model ----------------------------------------------
# Each participant has their own level (b0i) but all share the arm's slope.
m04_ri <- lmer(sbp ~ month * arm + (1 | id), data = walk, control = lmer_ctrl)
summary(m04_ri)
VarCorr(m04_ri)
icc_ri <- icc_lmer(m04_ri)   # share of total variance that is between-person
icc_ri
save_model(m04_ri, "m04_random_intercept")

# ---- 3. Random-slope model: the linear-time reference model -------------------
m04_lin <- lmer(sbp ~ month * arm + (1 + month | id), data = walk,
                control = lmer_ctrl)
summary(m04_lin)
save_model(m04_lin, "m04_linear_time")
save_table(fixef_table(m04_lin), "m04_linear_time_fixef")
save_table(varcomp_table(m04_lin), "m04_linear_time_varcomp")

# Is the random slope needed? REML LRT (same fixed effects), with the
# boundary-corrected p-value: the slope variance is tested at zero and the
# covariance is free, so the reference is a 50:50 mixture of chi2(1), chi2(2).
lrt_slope <- lr_test(m04_ri, m04_lin)
lrt_slope$p_mixture <- p_boundary_mixture(lrt_slope$chisq, df_mix = c(1, 2))
lrt_slope
save_table(lrt_slope, "m04_lrt_random_slope")

# Implied marginal covariance and correlation at the scheduled times.
# Random slopes make the variance and the correlation depend on time.
V_lin <- implied_cov_lmer(m04_lin, sched)
V_ri  <- implied_cov_lmer(m04_ri, sched)
round(V_lin, 1)
cov_to_cor(V_lin)
cov_to_cor(V_ri)
obs_cor <- visit_correlation(walk, id, month_sched, sbp)
round(obs_cor, 2)
save_table(round(V_lin, 2), "m04_implied_cov_linear_time")

# Variance as a smooth function of time, for a figure
G   <- as.matrix(VarCorr(m04_lin)$id)
s2  <- sigma(m04_lin)^2
tt  <- seq(0, 12, by = 0.25)
var_t <- data.frame(month = tt,
                    total = G[1, 1] + 2 * tt * G[1, 2] + tt^2 * G[2, 2] + s2,
                    model = "Random slope")
var_t <- rbind(var_t, data.frame(month = tt,
  total = as.numeric(VarCorr(m04_ri)$id) + sigma(m04_ri)^2,
  model = "Random intercept"))
obs_var <- walk |> group_by(month_sched) |>
  summarise(v = var(sbp))
p_var <- ggplot(var_t, aes(month, total, linetype = model)) +
  geom_line() +
  geom_point(data = obs_var, aes(month_sched, v), inherit.aes = FALSE,
             colour = course_colours["orange"], size = 2) +
  labs(x = "Months since randomisation", y = "Variance of SBP (mmHg²)",
       linetype = NULL)
save_fig(p_var, "m04_implied_variance")

# ---- 4. Coding of time --------------------------------------------------------
# 4a. Centred at month 12: arm "main effects" become differences at month 12
walk <- walk |> mutate(month12 = month - 12)
m04_c12 <- lmer(sbp ~ month12 * arm + (1 + month12 | id), data = walk,
                control = lmer_ctrl)
fixef(m04_c12)
save_model(m04_c12, "m04_linear_time_centred12")

# 4b. Quadratic time: a first look at curvature (Module 9 goes further)
m04_quad <- lmer(sbp ~ (month + I(month^2)) * arm + (1 + month | id),
                 data = walk, control = lmer_ctrl)
anova(m04_lin, m04_quad)   # refits with ML: fixed effects differ

# 4c. Categorical time: the primary model
m04_cat <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk,
                control = lmer_ctrl)
summary(m04_cat)
save_model(m04_cat, "m04_primary_categorical_time")
save_table(fixef_table(m04_cat), "m04_primary_fixef")
save_table(varcomp_table(m04_cat), "m04_primary_varcomp")

# Fitted mean trajectories from the three mean models against the observed
# means, to show what each coding of time assumes.
nd <- expand.grid(month = seq(0, 12, 0.25), arm = levels(walk$arm))
nd$lin  <- predict(m04_lin,  newdata = nd, re.form = NA)
nd$quad <- predict(m04_quad, newdata = nd, re.form = NA)
cat_means <- as.data.frame(emmeans(m04_cat, ~ visit_f | arm)) |>
  mutate(month = as.numeric(as.character(visit_f)))
p_time <- ggplot(nd, aes(month, colour = arm)) +
  geom_line(aes(y = lin, linetype = "Linear")) +
  geom_line(aes(y = quad, linetype = "Quadratic")) +
  geom_point(data = cat_means, aes(month, emmean, shape = "Categorical"),
             size = 2) +
  scale_colour_arm(labels = walk_arm_labels, name = NULL) +
  labs(x = "Months since randomisation", y = "Estimated mean SBP (mmHg)",
       linetype = NULL, shape = NULL) +
  theme(legend.box = "vertical")
save_fig(p_time, "m04_time_codings")

# ---- 5. Primary estimand: difference in mean change at month 12 ---------------
# EMMs at each visit, change from baseline within arm, then each programme
# minus usual care. adjust = "none": each comparison is its own estimand
# (multiplicity is discussed in Module 6).
emm_cat <- emmeans(m04_cat, ~ visit_f | arm)
chg12   <- contrast(emm_cat, list(chg12 = c(-1, 0, 0, 0, 0, 1)))
dmc12   <- contrast(chg12, "trt.vs.ctrl", by = NULL, adjust = "none")
chg12_s <- summary(chg12, infer = TRUE)
dmc12_s <- summary(dmc12, infer = TRUE)
chg12_s
dmc12_s

# Comparison with the truth: possible only because the data are simulated
true12 <- truth$true_change |> filter(month == 12)
true_chg <- setNames(true12$change, true12$arm)
compare_truth <- data.frame(
  quantity = c("Change UC", "Change SP", "Change EP", "SP - UC", "EP - UC"),
  estimate = c(chg12_s$estimate, dmc12_s$estimate),
  lower    = c(chg12_s$lower.CL, dmc12_s$lower.CL),
  upper    = c(chg12_s$upper.CL, dmc12_s$upper.CL),
  truth    = c(true_chg[c("UC", "SP", "EP")],
               true_chg["SP"] - true_chg["UC"],
               true_chg["EP"] - true_chg["UC"]))
compare_truth
save_table(compare_truth, "m04_primary_vs_truth")

# The same estimand from the linear-time model (assumes straight lines)
emm_lin <- emmeans(m04_lin, ~ month | arm, at = list(month = c(0, 12)))
dmc12_lin <- contrast(contrast(emm_lin, list(chg12 = c(-1, 1))),
                      "trt.vs.ctrl", by = NULL, adjust = "none")
summary(dmc12_lin, infer = TRUE)

# ---- 6. ML versus REML -----------------------------------------------------
m04_cat_ml <- update(m04_cat, REML = FALSE)
reml_ml <- rbind(
  data.frame(method = "REML", varcomp_table(m04_cat)),
  data.frame(method = "ML", varcomp_table(m04_cat_ml)))
reml_ml
save_table(reml_ml, "m04_reml_vs_ml_varcomp")

# Simulation: small trial (16 participants), random intercept only.
# ML ignores the degrees of freedom used by the fixed effects, so it
# underestimates variance components; REML corrects for this.
# A few of the 1,000 fits give borderline convergence warnings (gradient just
# above tolerance); these are expected with such small datasets and are left
# visible on purpose (Module 11).
set.seed(404)
sim_reml <- t(replicate(500, {
  d <- sim_two_arm(n_per_arm = 8, times = c(0, 3, 6, 9, 12),
                   sd_int = 11, sd_slope = 0, sd_resid = 6)
  f_reml <- suppressMessages(lmer(y ~ time * trt + (1 | id), data = d))
  f_ml   <- suppressMessages(lmer(y ~ time * trt + (1 | id), data = d,
                                  REML = FALSE))
  c(tau2_reml = as.numeric(VarCorr(f_reml)$id),
    tau2_ml   = as.numeric(VarCorr(f_ml)$id),
    s2_reml   = sigma(f_reml)^2, s2_ml = sigma(f_ml)^2)
}))
sim_summary <- data.frame(
  parameter = c("tau0^2", "tau0^2", "sigma^2", "sigma^2"),
  method = c("REML", "ML", "REML", "ML"),
  truth = c(121, 121, 36, 36),
  mean_estimate = colMeans(sim_reml),
  mcse = apply(sim_reml, 2, sd) / sqrt(nrow(sim_reml)))
sim_summary
save_table(sim_summary, "m04_reml_ml_simulation")

# ---- 7. BLUPs and shrinkage ---------------------------------------------------
# Per-participant OLS lines (no pooling) versus the model's predicted lines
# (partial pooling). Participants with one observation have no OLS slope.
ols_fits <- lme4::lmList(sbp ~ month | id, data = walk)
ols_cf <- coef(ols_fits) |>
  tibble::rownames_to_column("id") |>
  rename(ols_int = `(Intercept)`, ols_slope = month)
blup_cf <- coef(m04_lin)$id |>
  tibble::rownames_to_column("id") |>
  transmute(id,
            # participant-specific slope = month + month:arm + b1i
            blup_int = `(Intercept)`, blup_month = month,
            blup_SP = `month:armSP`, blup_EP = `month:armEP`)
n_obs <- walk |> count(id, arm, name = "n_obs")
shrink <- n_obs |>
  left_join(ols_cf, by = "id") |>
  left_join(blup_cf, by = "id") |>
  mutate(blup_slope = blup_month + (arm == "SP") * blup_SP +
                      (arm == "EP") * blup_EP,
         n_group = cut(n_obs, c(0, 3, 5, 6),
                       labels = c("2-3 visits", "4-5 visits", "6 visits")))
save_table(shrink, "m04_ols_vs_blup_slopes")

p_shrink <- shrink |>
  filter(!is.na(ols_slope)) |>
  ggplot(aes(ols_slope, blup_slope, colour = n_group)) +
  geom_abline(slope = 1, intercept = 0, colour = "#AAAAAA") +
  geom_point(alpha = 0.7, size = 1.4) +
  scale_colour_manual(values = unname(course_colours[c("orange", "blue",
                                                       "grey")]),
                      name = NULL) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(x = "Separate OLS slope (mmHg/month)",
       y = "Model-based slope, BLUP (mmHg/month)")
save_fig(p_shrink, "m04_shrinkage_slopes")

# ---- 8. The same model with nlme::lme -------------------------------------------
m04_cat_lme <- lme(sbp ~ visit_f * arm, random = ~ 1 + month | id,
                   data = walk, method = "REML")
agree <- data.frame(
  term = names(fixef(m04_cat)),
  lmer = unname(fixef(m04_cat)),
  lme = unname(nlme::fixef(m04_cat_lme)),
  se_lmer = sqrt(diag(as.matrix(vcov(m04_cat)))),
  se_lme = sqrt(diag(vcov(m04_cat_lme))))
agree
save_table(agree, "m04_lmer_vs_lme_fixef")
save_model(m04_cat_lme, "m04_primary_lme")
dmc12_lme <- contrast(contrast(emmeans(m04_cat_lme, ~ visit_f | arm),
                               list(chg12 = c(-1, 0, 0, 0, 0, 1))),
                      "trt.vs.ctrl", by = NULL, adjust = "none")
summary(dmc12_lme, infer = TRUE)

# ---- 9. Basic residual checks (full diagnostics in Module 11) -------------------
diag_df <- data.frame(fitted = fitted(m04_cat), resid = resid(m04_cat))
re_df <- ranef(m04_cat)$id
p_rf <- ggplot(diag_df, aes(fitted, resid)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_hline(yintercept = 0, colour = "#888888") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x,
              colour = course_colours["orange"]) +
  labs(x = "Fitted value (mmHg)", y = "Conditional residual (mmHg)")
p_qq <- ggplot(diag_df, aes(sample = resid)) +
  stat_qq(alpha = 0.3, size = 0.8) + stat_qq_line(colour = "#888888") +
  labs(x = "Normal quantile", y = "Residual (mmHg)")
p_qq_b0 <- ggplot(re_df, aes(sample = `(Intercept)`)) +
  stat_qq(alpha = 0.4, size = 0.8) + stat_qq_line(colour = "#888888") +
  labs(x = "Normal quantile", y = "Predicted b0i (mmHg)")
p_qq_b1 <- ggplot(re_df, aes(sample = month)) +
  stat_qq(alpha = 0.4, size = 0.8) + stat_qq_line(colour = "#888888") +
  labs(x = "Normal quantile", y = "Predicted b1i (mmHg/month)")
p_diag <- (p_rf | p_qq) / (p_qq_b0 | p_qq_b1)
save_fig(p_diag, "m04_residual_checks", width = 8, height = 6)

sessionInfo()
