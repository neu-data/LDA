# =============================================================================
# 05_covariance_models.R
# Module 5. Modelling covariance and choosing a model
#
# Purpose : Fit the main residual correlation structures (compound symmetry,
#           AR(1), continuous-time exponential, unstructured) and random-effects
#           structures to the WALK data, compute the marginal covariance each
#           one implies, compare the models (LRT, AIC, BIC), illustrate boundary
#           tests and singular fits, and show how much (or how little) the
#           treatment-effect estimate depends on the covariance model.
# Inputs  : data/walk_long.csv (via load_walk())
# Outputs : outputs/figures/m05_*.png, outputs/tables/m05_*.csv,
#           outputs/models/m05_*.rds
# Module  : 05-covariance-models
#
# SIMULATED DATA: the WALK study does not exist.
# =============================================================================

source(here::here("R", "setup.R"))
walk <- load_walk()

sched <- c(0, 1, 3, 6, 9, 12)

# ---- helper: implied marginal covariance at scheduled times ------------------
# Computes V = Z G Z' + R for gls and lme fits at a given set of times. The
# residual part R = sigma^2 * S C S, where S holds the varIdent standard-deviation
# ratios and C is the correlation matrix implied by the corStruct. Writing this
# out by hand (rather than calling getVarCov()) makes the formula visible and
# evaluates every model at the same scheduled times.
implied_cov <- function(model, times = sched) {
  k <- length(times)
  # residual correlation C
  cs <- model$modelStruct$corStruct
  C <- diag(k)
  if (!is.null(cs)) {
    if (inherits(cs, "corCompSymm")) {
      rho <- coef(cs, unconstrained = FALSE)
      C <- matrix(rho, k, k); diag(C) <- 1
    } else if (inherits(cs, c("corAR1", "corARMA"))) {
      # corAR1() objects are stored with class corARMA (p = 1)
      phi <- coef(cs, unconstrained = FALSE)
      C <- phi^abs(outer(seq_len(k), seq_len(k), "-"))   # by visit index
    } else if (inherits(cs, "corExp")) {
      p <- coef(cs, unconstrained = FALSE)
      nug <- if ("nugget" %in% names(p)) p[["nugget"]] else 0
      C <- (1 - nug) * exp(-abs(outer(times, times, "-")) / p[["range"]])
      diag(C) <- 1
    } else if (inherits(cs, "corSymm")) {
      r <- coef(cs, unconstrained = FALSE)
      C[lower.tri(C)] <- r
      C <- C + t(C) - diag(k)
    }
  }
  # heterogeneous SD ratios (varIdent), in the order of the scheduled visits
  vs <- model$modelStruct$varStruct
  s <- rep(1, k)
  if (!is.null(vs)) {
    w <- coef(vs, unconstrained = FALSE, allCoef = TRUE)
    s <- w[as.character(times)]
  }
  R <- model$sigma^2 * diag(s) %*% C %*% diag(s)
  # random-effects part Z G Z'
  if (inherits(model, "lme")) {
    G <- as.matrix(getVarCov(model))
    Z <- cbind(1, times)[, seq_len(ncol(G)), drop = FALSE]
    V <- Z %*% G %*% t(Z) + R
  } else {
    V <- R
  }
  dimnames(V) <- list(times, times)
  V
}

# long data frame of a matrix, for heatmaps
mat_long <- function(M, label) {
  data.frame(row = factor(rep(rownames(M), ncol(M)), levels = rownames(M)),
             col = factor(rep(colnames(M), each = nrow(M)),
                          levels = colnames(M)),
             value = as.vector(M), model = label)
}

# ---- observed covariance of residuals -----------------------------------------
# Remove the arm-by-visit means first (an OLS fit of the saturated mean model),
# otherwise differences between arm means would inflate the covariances.
walk$res_ols <- residuals(lm(sbp ~ visit_f * arm, data = walk))
obs_cov <- walk |>
  select(id, visit_f, res_ols) |>
  pivot_wider(names_from = visit_f, values_from = res_ols) |>
  select(-id) |>
  cov(use = "pairwise.complete.obs")
obs_cor <- cov2cor(obs_cov)
save_table(round(as.data.frame(obs_cov), 1) |> mutate(month = sched, .before = 1),
           "m05_observed_cov")
save_table(round(as.data.frame(obs_cor), 3) |> mutate(month = sched, .before = 1),
           "m05_observed_cor")

# ---- residual correlation structures with gls (REML) -------------------------
# All models share the same saturated mean (visit_f * arm), so their REML
# likelihoods are comparable.
f_mean <- sbp ~ visit_f * arm

m_ind <- gls(f_mean, data = walk)
m_cs  <- gls(f_mean, data = walk,
             correlation = corCompSymm(form = ~ 1 | id))
m_ar1 <- gls(f_mean, data = walk,
             correlation = corAR1(form = ~ visit | id))
m_exp <- gls(f_mean, data = walk,
             correlation = corExp(form = ~ month | id, nugget = TRUE))
m_cs_het <- gls(f_mean, data = walk,
                correlation = corCompSymm(form = ~ 1 | id),
                weights = varIdent(form = ~ 1 | visit_f))
# unstructured: a free correlation for every pair of visits and a free
# variance at every visit (6 variances + 15 correlations = 21 parameters)
m_un  <- gls(f_mean, data = walk,
             correlation = corSymm(form = ~ visit | id),
             weights = varIdent(form = ~ 1 | visit_f))

# ---- random effects, alone and combined with serial correlation (lme, REML) --
m_ri     <- lme(f_mean, random = ~ 1 | id, data = walk)
m_rs     <- lme(f_mean, random = ~ 1 + month | id, data = walk)
m_ri_exp <- lme(f_mean, random = ~ 1 | id, data = walk,
                correlation = corExp(form = ~ month | id, nugget = TRUE))
m_rs_exp <- lme(f_mean, random = ~ 1 + month | id, data = walk,
                correlation = corExp(form = ~ month | id, nugget = TRUE))

mods <- list(
  "Independence"                    = m_ind,
  "Compound symmetry"               = m_cs,
  "AR(1) by visit"                  = m_ar1,
  "Exponential + nugget"            = m_exp,
  "CS + heterogeneous variances"    = m_cs_het,
  "Unstructured"                    = m_un,
  "Random intercept"                = m_ri,
  "Random intercept + slope"        = m_rs,
  "Random intercept + exponential"  = m_ri_exp,
  "Random int. + slope + exponential" = m_rs_exp
)
for (nm in names(mods)) save_model(mods[[nm]],
  paste0("m05_", gsub("[^a-z0-9]+", "_", tolower(nm))))

# ---- model comparison table --------------------------------------------------
cmp <- data.frame(
  model  = names(mods),
  n_cov  = sapply(mods, function(m) attr(logLik(m), "df")) - length(coef(m_ind)),
  logLik = sapply(mods, function(m) as.numeric(logLik(m))),
  AIC    = sapply(mods, AIC),
  BIC    = sapply(mods, BIC)
) |>
  mutate(dAIC = AIC - min(AIC), dBIC = BIC - min(BIC))
rownames(cmp) <- NULL
print(cmp, digits = 6)
save_table(cmp, "m05_model_comparison")

# nested LRTs (REML, same mean model)
print(anova(m_cs, m_cs_het))       # homogeneous vs heterogeneous variances
print(anova(m_cs, m_un))           # CS nested in UN
print(anova(m_ri, m_ri_exp))       # adding serial correlation to RI
print(anova(m_rs, m_rs_exp))       # adding serial correlation to RI+RS

# ---- implied vs observed correlation: heatmaps --------------------------------
heat_models <- c("Compound symmetry", "AR(1) by visit", "Exponential + nugget",
                 "Unstructured", "Random intercept + slope",
                 "Random int. + slope + exponential")
heat <- bind_rows(
  mat_long(obs_cor, "Observed (residuals)"),
  lapply(heat_models, function(nm) mat_long(cov2cor(implied_cov(mods[[nm]])), nm))
) |>
  mutate(model = factor(model, levels = c("Observed (residuals)", heat_models)))

p_heat <- ggplot(heat, aes(col, row, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", limits = c(0.3, 1), breaks = c(0.4, 0.6, 0.8, 1),
                      name = "Correlation") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~ model, ncol = 4) +
  labs(x = "Scheduled month", y = "Scheduled month") +
  theme(panel.grid = element_blank())
save_fig(p_heat, "m05_correlation_heatmaps", width = 11, height = 6)

# correlation as a function of lag: observed vs implied
lag_df <- bind_rows(lapply(c("Observed (residuals)", heat_models), function(nm) {
  M <- if (nm == "Observed (residuals)") obs_cor else cov2cor(implied_cov(mods[[nm]]))
  idx <- which(upper.tri(M), arr.ind = TRUE)
  data.frame(model = nm, lag = sched[idx[, 2]] - sched[idx[, 1]],
             from = sched[idx[, 1]], cor = M[idx])
}))
p_lag <- ggplot(filter(lag_df, model != "Observed (residuals)") |>
                  mutate(model = factor(model, levels = heat_models)),
                aes(lag, cor)) +
  geom_point(data = filter(lag_df, model == "Observed (residuals)") |>
               select(-model), colour = "#9A9A9A", size = 1.6) +
  geom_point(colour = "#0072B2", size = 1.6) +
  facet_wrap(~ model, ncol = 3) +
  coord_cartesian(ylim = c(0.3, 1)) +
  labs(x = "Lag between visits (months)", y = "Correlation",
       caption = "Grey: observed residual correlations. Blue: implied by the model.")
save_fig(p_lag, "m05_correlation_by_lag", width = 9, height = 5.5)

# ---- variance components of the combined model --------------------------------
print(summary(m_rs_exp))
print(intervals(m_rs_exp, which = "var-cov"))

# ---- semivariogram of normalized residuals from the random-slope model --------
# Values below 1 at short lags mean that correlation remains after the random
# intercept and slope: evidence for a serial correlation component.
vg <- Variogram(m_rs, form = ~ month | id, resType = "normalized",
                maxDist = 12)
p_vg <- ggplot(vg, aes(dist, variog)) +
  geom_hline(yintercept = 1, colour = "#9A9A9A", linetype = 2) +
  geom_point(colour = "#222222") +
  geom_smooth(method = "loess", se = FALSE, colour = "#0072B2",
              linewidth = 0.7, formula = y ~ x) +
  labs(x = "Distance between measurements (months)",
       y = "Semivariogram of normalized residuals")
save_fig(p_vg, "m05_variogram_rs")

# ---- treatment effect across covariance models ---------------------------------
# Difference in mean change from baseline to month 12, programme minus usual care
cells <- expand.grid(visit_f = levels(walk$visit_f), arm = levels(walk$arm))
cell_w <- function(a, v) as.numeric(cells$arm == a & cells$visit_f == v)
chg12  <- function(a) cell_w(a, "12") - cell_w(a, "0")
L_dic  <- list("SP - UC" = chg12("SP") - chg12("UC"),
               "EP - UC" = chg12("EP") - chg12("UC"))

effect_by_model <- bind_rows(lapply(names(mods), function(nm) {
  m <- mods[[nm]]
  e <- if (inherits(m, "gls")) emmeans(m, ~ visit_f * arm, mode = "df.error") else
    emmeans(m, ~ visit_f * arm)
  s <- summary(contrast(e, L_dic), infer = TRUE)
  data.frame(model = nm, contrast = s$contrast, estimate = s$estimate,
             se = s$SE, lower = s$lower.CL, upper = s$upper.CL)
}))
m_prim <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk,
               control = lmer_ctrl)
save_model(m_prim, "m05_primary_lmer")
s_prim <- summary(contrast(emmeans(m_prim, ~ visit_f * arm), L_dic),
                  infer = TRUE)
effect_by_model <- bind_rows(
  effect_by_model,
  data.frame(model = "Primary model (lmer)", contrast = s_prim$contrast,
             estimate = s_prim$estimate, se = s_prim$SE,
             lower = s_prim$lower.CL, upper = s_prim$upper.CL)
)
print(effect_by_model, digits = 3)
save_table(effect_by_model, "m05_effect_by_model")

p_forest <- effect_by_model |>
  mutate(model = factor(model, levels = rev(unique(model)))) |>
  ggplot(aes(estimate, model)) +
  geom_vline(xintercept = 0, colour = "#9A9A9A") +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.3) +
  facet_wrap(~ contrast) +
  labs(x = "Difference in mean change, baseline to month 12 (mmHg)", y = NULL)
save_fig(p_forest, "m05_effect_by_model", width = 9, height = 5)

# ---- REML and fixed effects: a warning in two lines ---------------------------
# The same model with time in months or weeks: identical fit, identical ML
# log-likelihood, but a different REML log-likelihood. REML likelihoods depend
# on the parameterisation of X, so they cannot compare different mean models.
walk$week <- walk$month * 52 / 12
m_lin_mo <- lmer(sbp ~ month * arm + (1 + month | id), data = walk,
                 control = lmer_ctrl)
m_lin_wk <- lmer(sbp ~ week * arm + (1 + month | id), data = walk,
                 control = lmer_ctrl)
reml_demo <- data.frame(
  model = c("time in months", "time in weeks"),
  REML_logLik = c(logLik(m_lin_mo), logLik(m_lin_wk)),
  ML_logLik   = c(logLik(refitML(m_lin_mo)), logLik(refitML(m_lin_wk)))
)
print(reml_demo)
save_table(reml_demo, "m05_reml_parameterisation")
# anova() on lmer fits refits with ML automatically for fixed-effect comparisons
print(anova(m_lin_mo, m_prim))

# ---- boundary tests ------------------------------------------------------------
m_ri_lmer <- lmer(sbp ~ visit_f * arm + (1 | id), data = walk,
                  control = lmer_ctrl)
lrt_slope <- lr_test(m_ri_lmer, m_prim)        # slope variance + covariance
lrt_slope$p_mixture <- p_boundary_mixture(lrt_slope$chisq, df_mix = c(1, 2))
m_clinic <- lmer(sbp ~ visit_f * arm + (1 + month | id) + (1 | clinic),
                 data = walk, control = lmer_ctrl)
lrt_clinic <- lr_test(m_prim, m_clinic)        # one variance on the boundary
lrt_clinic$p_mixture <- p_boundary_mixture(lrt_clinic$chisq)
boundary <- bind_rows(cbind(test = "random slope (variance and covariance)",
                            lrt_slope),
                      cbind(test = "clinic random intercept", lrt_clinic))
print(boundary)
save_table(boundary, "m05_boundary_tests")

# ---- singular fit ---------------------------------------------------------------
# A random slope for clinic: with only 8 clinics, the clinic slope variance and
# its correlation with the clinic intercept cannot be estimated.
m_clinic_slope <- lmer(sbp ~ visit_f * arm + (1 + month | id) +
                         (1 + month | clinic), data = walk,
                       control = lmer_ctrl)
print(VarCorr(m_clinic_slope))
print(isSingular(m_clinic_slope))
save_model(m_clinic_slope, "m05_singular_clinic_slope")

effect_clinic <- bind_rows(lapply(list(
  "Primary model"                   = m_prim,
  "+ clinic intercept"              = m_clinic,
  "+ clinic intercept and slope (singular)" = m_clinic_slope), function(m) {
    s <- summary(contrast(emmeans(m, ~ visit_f * arm), L_dic), infer = TRUE)
    data.frame(contrast = s$contrast, estimate = s$estimate, se = s$SE)
  }), .id = "model")
print(effect_clinic, digits = 3)
save_table(effect_clinic, "m05_effect_clinic_models")

message("Module 5 script finished.")
