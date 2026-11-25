# =============================================================================
# 08b_glmm_count.R
# Purpose : Poisson and negative binomial mixed models for the DIARY migraine
#           counts with an exposure offset; overdispersion checks; testing
#           zero inflation only where justified (DHARMa, ZINB comparison);
#           rate ratios and expected counts per 28 days with emmeans;
#           conditional versus marginal expected counts; comparison with GEE
#           and with the true values.
# Inputs  : data/diary_long.csv, data/diary_truth.rds
# Outputs : outputs/tables/m08_diary_*.csv, outputs/figures/m08_dharma_diary_*.png,
#           outputs/models/m08_diary_*.rds
# Module  : 8 (Generalized linear mixed models)
# All data are simulated.
# =============================================================================

source(here::here("R", "setup.R"))

diary <- load_diary() |> arrange(id, period)
truth <- readRDS(here("data", "diary_truth.rds"))

# ---- 1. Poisson and negative binomial GLMMs ----------------------------------
f_diary <- migraine_days ~ post * arm + offset(log(days)) + (1 | id)
pois_glmm <- glmmTMB(f_diary, data = diary, family = poisson)
nb_glmm   <- glmmTMB(f_diary, data = diary, family = nbinom2)
print(summary(nb_glmm))

# ---- 2. Overdispersion --------------------------------------------------------
# Pearson chi-square / residual df for the Poisson GLMM: values well above 1
# indicate overdispersion beyond what the random intercept explains.
pearson_ratio <- function(m) sum(residuals(m, type = "pearson")^2) / df.residual(m)
disp_tab <- data.frame(model = c("Poisson GLMM", "NB GLMM"),
                       pearson_ratio = c(pearson_ratio(pois_glmm),
                                         pearson_ratio(nb_glmm)))
print(disp_tab)

res_pois <- simulateResiduals(pois_glmm, n = 1000, seed = 8)
res_nb   <- simulateResiduals(nb_glmm, n = 1000, seed = 8)
dharma_tab <- data.frame(
  model = c("Poisson GLMM", "NB GLMM"),
  dispersion = c(testDispersion(res_pois, plot = FALSE)$statistic,
                 testDispersion(res_nb, plot = FALSE)$statistic),
  p_dispersion = c(testDispersion(res_pois, plot = FALSE)$p.value,
                   testDispersion(res_nb, plot = FALSE)$p.value),
  zero_ratio = c(testZeroInflation(res_pois, plot = FALSE)$statistic,
                 testZeroInflation(res_nb, plot = FALSE)$statistic),
  p_zero = c(testZeroInflation(res_pois, plot = FALSE)$p.value,
             testZeroInflation(res_nb, plot = FALSE)$p.value),
  p_uniform = c(testUniformity(res_pois, plot = FALSE)$p.value,
                testUniformity(res_nb, plot = FALSE)$p.value))
print(dharma_tab)
save_table(dharma_tab, "m08_diary_dharma_tests")
for (nm in c("pois", "nb")) {
  png(here("outputs", "figures", paste0("m08_dharma_diary_", nm, ".png")),
      width = 1600, height = 800, res = 200)
  plot(if (nm == "pois") res_pois else res_nb, quantreg = FALSE)
  dev.off()
}

# ---- 3. Is zero inflation needed? --------------------------------------------
# Only the NB model is a sensible baseline: the Poisson model's apparent
# excess of zeros is a symptom of overdispersion, not of a zero process.
zinb_glmm <- glmmTMB(f_diary, ziformula = ~ 1, data = diary, family = nbinom2)
print(summary(zinb_glmm)$coefficients$zi)
model_tab <- data.frame(
  model = c("Poisson GLMM", "NB GLMM", "ZINB GLMM"),
  df = c(attr(logLik(pois_glmm), "df"), attr(logLik(nb_glmm), "df"),
         attr(logLik(zinb_glmm), "df")),
  logLik = c(logLik(pois_glmm), logLik(nb_glmm), logLik(zinb_glmm)),
  AIC = c(AIC(pois_glmm), AIC(nb_glmm), AIC(zinb_glmm)))
print(model_tab)
save_table(model_tab, "m08_diary_model_comparison")

# ---- 4. Rate ratios and expected counts ---------------------------------------
est_tab <- data.frame(
  term = names(fixef(nb_glmm)$cond),
  poisson = fixef(pois_glmm)$cond,
  se_poisson = sqrt(diag(vcov(pois_glmm)$cond)),
  negbin = fixef(nb_glmm)$cond,
  se_negbin = sqrt(diag(vcov(nb_glmm)$cond)),
  truth = c(truth$beta0, truth$beta_post, 0, truth$beta_trt))
print(est_tab)
save_table(est_tab, "m08_diary_estimates")
save_table(fixef_table(nb_glmm, exponentiate = TRUE), "m08_diary_nb_rate_ratios")

# Expected migraine days per 28 days for a participant with b = 0
# (conditional) and averaged over participants (marginal: exp(sigma^2 / 2)
# larger for a log link; emmeans' bias adjustment is exact here).
sig_nb <- attr(VarCorr(nb_glmm)$cond$id, "stddev")
emm_c <- emmeans(nb_glmm, ~ arm | post, offset = log(28), type = "response")
emm_m <- emmeans(nb_glmm, ~ arm | post, offset = log(28), type = "response",
                 bias.adjust = TRUE, sigma = sig_nb)
print(emm_c); print(emm_m)
save_table(as.data.frame(emm_c), "m08_diary_per28_conditional")
save_table(as.data.frame(emm_m), "m08_diary_per28_marginal")
# Treatment effect: the ratio of rate ratios (interaction contrast)
rr_int <- confint(contrast(emmeans(nb_glmm, ~ post * arm),
                           interaction = c("revpairwise", "revpairwise")),
                  type = "response")
print(rr_int)

# Comparison with the Module 7 GEE: the slopes agree for a log link
gee_diary <- geeglm(migraine_days ~ post * arm + offset(log(days)),
                    id = factor(id), waves = period + 1, data = diary,
                    family = poisson, corstr = "exchangeable")
cmp <- data.frame(term = names(coef(gee_diary)),
                  gee = coef(gee_diary),
                  nb_glmm = fixef(nb_glmm)$cond,
                  nb_glmm_intercept_marginal =
                    fixef(nb_glmm)$cond + c(sig_nb^2 / 2, 0, 0, 0))
print(cmp)
save_table(cmp, "m08_diary_gee_vs_glmm")

# ---- save models -----------------------------------------------------------
save_model(pois_glmm, "m08_diary_poisson_glmm")
save_model(nb_glmm, "m08_diary_nb_glmm")
save_model(zinb_glmm, "m08_diary_zinb_glmm")
