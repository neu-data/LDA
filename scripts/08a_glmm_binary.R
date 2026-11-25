# =============================================================================
# 08a_glmm_binary.R
# Purpose : Logistic mixed models for the CLEAR binary outcome (itch).
#           Laplace versus adaptive Gauss-Hermite quadrature (glmer) and
#           glmmTMB; conditional interpretation; attenuation of marginal
#           (GEE) relative to conditional (GLMM) coefficients; conditional
#           versus marginal predicted probabilities with uncertainty;
#           DHARMa residual checks; comparison with the true values.
# Inputs  : data/clear_long.csv, data/clear_truth.rds
# Outputs : outputs/tables/m08_*.csv, outputs/figures/m08_*.png,
#           outputs/models/m08_*.rds
# Module  : 8 (Generalized linear mixed models)
# All data are simulated.
# =============================================================================

source(here::here("R", "setup.R"))

clear <- load_clear() |>
  arrange(id, week) |>
  mutate(age45 = (age - 45) / 10)          # truth is defined at age 45
truth <- readRDS(here("data", "clear_truth.rds"))

# ---- helpers -----------------------------------------------------------------
# Gauss-Hermite nodes and weights (Golub-Welsch). For b ~ N(0, s^2):
# E f(b) = sum_k w_k / sqrt(pi) * f(sqrt(2) * s * x_k).
gauss_hermite <- function(n) {
  i <- seq_len(n - 1)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}
gh <- gauss_hermite(30)
# Marginal (population-averaged) probability for linear predictor(s) eta
# when the random intercept has SD sigma
marginal_p <- function(eta, sigma) {
  sapply(eta, function(e)
    sum(gh$w / sqrt(pi) * plogis(e + sqrt(2) * sigma * gh$x)))
}

# ---- 1. Laplace versus adaptive Gauss-Hermite quadrature -------------------
f_clear <- itch ~ week * arm + (1 | id)
agq <- c(1, 5, 10, 25)
fits_q <- lapply(agq, function(q)
  glmer(f_clear, data = clear, family = binomial, nAGQ = q))
fit_tmb <- glmmTMB(f_clear, data = clear, family = binomial)
quad_tab <- rbind(
  do.call(rbind, lapply(seq_along(agq), function(k) {
    m <- fits_q[[k]]
    data.frame(method = paste0("glmer nAGQ = ", agq[k]), t(fixef(m)),
               sd_b = attr(VarCorr(m)$id, "stddev"),
               logLik = as.numeric(logLik(m)), check.names = FALSE)
  })),
  data.frame(method = "glmmTMB (Laplace)", t(fixef(fit_tmb)$cond),
             sd_b = attr(VarCorr(fit_tmb)$cond$id, "stddev"),
             logLik = as.numeric(logLik(fit_tmb)), check.names = FALSE))
rownames(quad_tab) <- NULL
print(quad_tab)
save_table(quad_tab, "m08_clear_laplace_vs_agq")

glmm_clear <- fits_q[[4]]                 # nAGQ = 25: the reference fit
print(summary(glmm_clear))
save_table(fixef_table(glmm_clear, exponentiate = TRUE), "m08_clear_glmm_or")

# ---- 2. Comparison with the true conditional values --------------------------
# Possible only because the data are simulated. The adjusted model has the
# same covariates as the generating model, so its intercept is comparable.
glmm_adj <- glmer(itch ~ week * arm + age45 + sex + (1 | id), data = clear,
                  family = binomial, nAGQ = 25)
truth_tab <- data.frame(
  parameter = c("(Intercept)", "week", "week:armActive", "age45",
                "sexFemale", "SD random intercept"),
  truth = c(-truth$cuts[2], truth$beta_week, truth$beta_trtweek,
            truth$beta_age, truth$beta_female, truth$sd_b),
  estimate = c(fixef(glmm_adj)[c("(Intercept)", "week", "week:armActive",
                                 "age45", "sexFemale")],
               attr(VarCorr(glmm_adj)$id, "stddev")),
  se = c(sqrt(diag(vcov(glmm_adj)))[c("(Intercept)", "week",
                                       "week:armActive", "age45",
                                       "sexFemale")], NA))
print(truth_tab)
save_table(truth_tab, "m08_clear_glmm_vs_truth")

# ---- 3. Conditional versus marginal coefficients ---------------------------
gee_clear <- geeglm(itch ~ week * arm, id = factor(id), waves = visit,
                    data = clear, family = binomial, corstr = "exchangeable")
sig <- attr(VarCorr(glmm_clear)$id, "stddev")
atten <- 1 / sqrt(1 + 0.346 * sig^2)
# Marginal coefficients implied exactly by the fitted GLMM: integrate over
# the random intercept at the observed weeks and fit the marginal model.
grid <- expand.grid(week = truth$weeks,
                    arm = factor(c("Placebo", "Active"),
                                 levels = c("Placebo", "Active")))
X <- model.matrix(~ week * arm, grid)
grid$p_marg <- marginal_p(X %*% fixef(glmm_clear), sig)
implied <- coef(suppressWarnings(glm(p_marg ~ week * arm, family = binomial,
                                     data = grid)))
true_att <- 1 / sqrt(1 + 0.346 * truth$sd_b^2)
atten_tab <- data.frame(
  term = names(fixef(glmm_clear)),
  glmm_conditional = fixef(glmm_clear),
  approx_marginal = fixef(glmm_clear) * atten,
  exact_marginal_from_glmm = implied,
  gee_exchangeable = coef(gee_clear),
  truth_conditional = c(0, truth$beta_week, 0, truth$beta_trtweek),
  truth_approx_marginal = c(0, truth$beta_week, 0, truth$beta_trtweek) * true_att)
rownames(atten_tab) <- NULL
print(atten_tab)
save_table(atten_tab, "m08_clear_attenuation")

# Figure: individual conditional curves versus the population-average curve
wk <- seq(0, 12, by = 0.25)
qs <- qnorm(c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95))
curves <- expand.grid(week = wk, b = qs * sig,
                      arm = factor(c("Placebo", "Active"),
                                   levels = c("Placebo", "Active")))
curves$p <- plogis(model.matrix(~ week * arm, curves) %*% fixef(glmm_clear) +
                     curves$b)
avg <- expand.grid(week = wk, arm = factor(c("Placebo", "Active"),
                                           levels = c("Placebo", "Active")))
Xa <- model.matrix(~ week * arm, avg)
avg$conditional_b0 <- plogis(Xa %*% fixef(glmm_clear))
avg$marginal <- marginal_p(Xa %*% fixef(glmm_clear), sig)
avg_long <- avg |>
  pivot_longer(c(conditional_b0, marginal), names_to = "curve",
               values_to = "p") |>
  mutate(curve = factor(curve, levels = c("conditional_b0", "marginal"),
                        labels = c("Participant with b = 0",
                                   "Population average")))
p_curves <- ggplot() +
  geom_line(data = curves, aes(week, p, group = b), colour = "grey70",
            linewidth = 0.4) +
  geom_line(data = avg_long, aes(week, p, linetype = curve),
            colour = "#222222", linewidth = 0.9) +
  facet_wrap(~ arm) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 12)) +
  labs(x = "Week", y = "Probability of moderate or severe itch",
       linetype = NULL)
save_fig(p_curves, "m08_conditional_vs_marginal_curves", width = 8, height = 4)

# ---- 4. Predicted probabilities with uncertainty ---------------------------
# Conditional at b = 0: emmeans, delta method on the logit scale
emm_cond <- emmeans(glmm_clear, ~ arm | week, at = list(week = c(0, 12)),
                    type = "response")
print(emm_cond)
# Marginal: integrate over b, with uncertainty from draws of (sigma, beta)
# from their joint approximate sampling distribution. glmer stores the
# Hessian of the deviance with respect to (theta, beta) for nAGQ > 0.
H <- glmm_clear@optinfo$derivs$Hessian
V_joint <- 2 * solve(H)
stopifnot(max(abs(V_joint[-1, -1] - as.matrix(vcov(glmm_clear)))) < 1e-6)
set.seed(808)
draws <- MASS::mvrnorm(4000, mu = c(getME(glmm_clear, "theta"),
                                    fixef(glmm_clear)), Sigma = V_joint)
pg <- expand.grid(arm = factor(c("Placebo", "Active"),
                               levels = c("Placebo", "Active")),
                  week = c(0, 12))
Xp <- model.matrix(~ week * arm, pg)
marg_draws <- apply(draws, 1, function(d) marginal_p(Xp %*% d[-1], abs(d[1])))
pg$marginal <- marginal_p(Xp %*% fixef(glmm_clear), sig)
pg$lower <- apply(marg_draws, 1, quantile, 0.025)
pg$upper <- apply(marg_draws, 1, quantile, 0.975)
print(pg)
rd_draws <- marg_draws[4, ] - marg_draws[3, ]      # week 12: Active - Placebo
rd_marg <- c(estimate = pg$marginal[4] - pg$marginal[3],
             quantile(rd_draws, c(0.025, 0.975)))
print(rd_marg)
save_table(pg, "m08_clear_marginal_probabilities")

# ---- 5. DHARMa residual checks ----------------------------------------------
res_clear <- simulateResiduals(glmm_clear, n = 1000, seed = 8)
print(testUniformity(res_clear, plot = FALSE))
print(testDispersion(res_clear, plot = FALSE))
png(here("outputs", "figures", "m08_dharma_clear.png"), width = 1600,
    height = 800, res = 200)
plot(res_clear, quantreg = FALSE)
dev.off()

# ---- save models -------------------------------------------------------------
save_model(glmm_clear, "m08_glmm_clear_agq25")
save_model(glmm_adj, "m08_glmm_clear_adjusted")
save_model(fit_tmb, "m08_glmmtmb_clear")
