# =============================================================================
# 11_diagnostics.R
# Module 11. Diagnostics and troubleshooting
#
# Purpose
#   Residual diagnostics, influence analysis, random-effects checks,
#   simulation-based residuals for GLMMs, GEE mean-model checks, convergence
#   troubleshooting and separation, using the simulated course datasets.
#
# Inputs
#   data/walk_long.csv, data/clear_long.csv, data/diary_long.csv (via loaders)
#
# Outputs
#   outputs/figures/m11_*.png      diagnostic figures
#   outputs/tables/m11_*.csv       summary tables
#   outputs/models/m11_*.rds       fitted models and cached simulation results
#
# The slow parts (360 leave-one-participant-out refits and two small
# simulation studies) are cached in outputs/models/. The notes, workbook and
# answer key read these caches when they exist and recompute them (with the
# same seeds, so with identical results) when they do not.
#
# All data are simulated. Results illustrate methods, not empirical findings.
# =============================================================================

source(here::here("R", "setup.R"))

# Run a computation once and cache the result in outputs/models/<name>.rds.
cache_or_run <- function(name, expr) {
  f <- here("outputs", "models", paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  x <- expr                     # evaluated only when there is no cache
  save_model(x, name)
  x
}

walk  <- load_walk()
clear <- load_clear()
diary <- load_diary()

# ---- 1. primary WALK model -------------------------------------------------
m_pri <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk,
              control = lmer_ctrl)
save_model(m_pri, "m11_walk_primary")

# ---- 2. conditional and marginal residuals ---------------------------------
# Conditional residual: y - X beta - Z b  (what is left after the participant's
# own trajectory). Marginal residual: y - X beta (deviation from the arm mean).
walk_res <- walk |>
  mutate(
    fit_cond  = fitted(m_pri),
    fit_marg  = predict(m_pri, re.form = NA),
    res_cond  = resid(m_pri),
    res_cond_std = resid(m_pri, scaled = TRUE),
    res_marg  = sbp - fit_marg
  )

# Cholesky-scaled marginal residuals: for each participant, L_i^{-1} r_i where
# V_i = L_i L_i'. If the implied covariance V_i is correct these are
# approximately independent N(0, 1), so they can be checked like OLS residuals.
chol_resid <- function(model, data, id = "id", time = "month") {
  G  <- as.matrix(VarCorr(model)$id)
  s2 <- sigma(model)^2
  r  <- data$sbp - predict(model, re.form = NA)
  out <- numeric(nrow(data))
  for (g in unique(data[[id]])) {
    rows <- which(data[[id]] == g)
    Z <- cbind(1, data[[time]][rows])
    V <- Z %*% G %*% t(Z) + diag(s2, length(rows))
    out[rows] <- forwardsolve(t(chol(V)), r[rows])
  }
  out
}
walk_res$res_chol <- chol_resid(m_pri, walk_res)

p_rf <- ggplot(walk_res, aes(fit_cond, res_cond_std)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "#B2182B") +
  labs(x = "Conditional fitted value (mmHg)",
       y = "Standardised conditional residual")
p_rt <- ggplot(walk_res, aes(month, res_cond_std)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "#B2182B") +
  labs(x = "Actual time (months)", y = "Standardised conditional residual")
p_arm <- ggplot(walk_res, aes(visit_f, res_cond_std, fill = arm)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
  scale_fill_arm(name = NULL) +
  labs(x = "Scheduled month", y = "Standardised conditional residual")
p_qq <- ggplot(walk_res, aes(sample = res_cond_std)) +
  stat_qq(size = 0.7, alpha = 0.5) + stat_qq_line(colour = "#B2182B") +
  labs(x = "Standard normal quantile", y = "Standardised conditional residual")
p_sl <- ggplot(walk_res, aes(fit_cond, sqrt(abs(res_cond_std)))) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "#B2182B") +
  labs(x = "Conditional fitted value (mmHg)",
       y = expression(sqrt("|standardised residual|")))
save_fig((p_rf | p_rt) / (p_qq | p_sl), "m11_conditional_residuals",
         width = 9, height = 7)
save_fig(p_arm, "m11_residuals_by_arm")

p_marg <- ggplot(walk_res, aes(visit_f, res_marg, fill = arm)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
  scale_fill_arm(name = NULL) +
  labs(x = "Scheduled month", y = "Marginal residual (mmHg)")
p_chol <- ggplot(walk_res, aes(sample = res_chol)) +
  stat_qq(size = 0.7, alpha = 0.5) + stat_qq_line(colour = "#B2182B") +
  labs(x = "Standard normal quantile", y = "Cholesky-scaled marginal residual")
save_fig(p_marg | p_chol, "m11_marginal_residuals", width = 9, height = 4)

# Serial correlation left in the residuals: correlation between consecutive
# within-participant Cholesky residuals (should be near 0 if V_i is right).
lag_tab <- walk_res |>
  arrange(id, month) |>
  group_by(id) |>
  mutate(lag_chol = lag(res_chol), lag_cond = lag(res_cond)) |>
  ungroup() |>
  summarise(r_chol = cor(res_chol, lag_chol, use = "complete.obs"),
            r_cond = cor(res_cond, lag_cond, use = "complete.obs"))
save_table(lag_tab, "m11_residual_lag_correlation")

# ---- 3. influence: leave one participant out --------------------------------
# For the primary model, the month-12 EP - UC difference in mean change equals
# the coefficient visit_f12:armEP, so fixef() is enough (no emmeans needed).
ids <- unique(walk$id)
loo <- cache_or_run("m11_loo", {
  b_full <- fixef(m_pri)
  V_inv  <- solve(as.matrix(vcov(m_pri)))
  theta0 <- getME(m_pri, "theta")
  res <- lapply(ids, function(g) {
    m_g <- lmer(sbp ~ visit_f * arm + (1 + month | id),
                data = walk[walk$id != g, ], control = lmer_ctrl,
                start = theta0)          # warm start from the full fit
    d <- b_full - fixef(m_g)
    data.frame(id = g,
               ep12 = unname(fixef(m_g)["visit_f12:armEP"]),
               sp12 = unname(fixef(m_g)["visit_f12:armSP"]),
               cook = as.numeric(t(d) %*% V_inv %*% d) / length(d))
  })
  bind_rows(res)
})
loo <- loo |>
  left_join(walk |> distinct(id, arm), by = "id") |>
  left_join(walk |> count(id, name = "n_obs"), by = "id") |>
  mutate(dfbeta_ep12 = fixef(m_pri)["visit_f12:armEP"] - ep12)
save_table(loo, "m11_loo")

se_ep12 <- sqrt(vcov(m_pri)["visit_f12:armEP", "visit_f12:armEP"])
top_loo <- loo |> arrange(desc(abs(dfbeta_ep12))) |> slice_head(n = 5)

p_loo <- ggplot(loo, aes(reorder(id, ep12), ep12, colour = arm)) +
  geom_point(size = 0.9) +
  geom_hline(yintercept = fixef(m_pri)["visit_f12:armEP"], linetype = 2) +
  scale_colour_arm(name = NULL) +
  labs(x = "Participant omitted (ordered)",
       y = "EP - UC difference in change at month 12 (mmHg)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_cook <- ggplot(loo, aes(seq_along(id), cook, colour = arm)) +
  geom_point(size = 0.9) +
  scale_colour_arm(name = NULL) +
  labs(x = "Participant (index)", y = "Cook's distance (fixed effects)")
save_fig(p_loo / p_cook, "m11_influence", width = 8, height = 6.5)

# ---- 4. random effects: BLUPs and shrinkage ---------------------------------
re <- ranef(m_pri)$id
names(re) <- c("b0", "b1")
p_blup <- (ggplot(re, aes(sample = b0)) + stat_qq(size = 0.7) +
             stat_qq_line(colour = "#B2182B") +
             labs(x = "Normal quantile", y = "Predicted random intercept (mmHg)")) |
  (ggplot(re, aes(sample = b1)) + stat_qq(size = 0.7) +
     stat_qq_line(colour = "#B2182B") +
     labs(x = "Normal quantile", y = "Predicted random slope (mmHg/month)"))
save_fig(p_blup, "m11_blup_qq", width = 8, height = 3.8)

# Shrinkage demonstration: true intercepts from a two-component mixture
# (clearly bimodal), only 2 noisy measurements each. BLUPs pull each value
# towards the mean by tau^2 / (tau^2 + sigma^2 / n_i), blurring the modes.
set.seed(1101)
n_sh <- 400
true_b <- ifelse(runif(n_sh) < 0.5, -10, 10) + rnorm(n_sh, 0, 3)
sh_dat <- data.frame(id = factor(rep(seq_len(n_sh), each = 2)),
                     y  = 140 + rep(true_b, each = 2) + rnorm(2 * n_sh, 0, 15))
m_sh <- lmer(y ~ 1 + (1 | id), data = sh_dat, control = lmer_ctrl)
sh_plot <- bind_rows(
  data.frame(value = true_b, what = "True random intercepts"),
  data.frame(value = ranef(m_sh)$id[, 1], what = "Predicted (BLUPs)"))
sh_plot$what <- factor(sh_plot$what, levels = c("True random intercepts",
                                                 "Predicted (BLUPs)"))
p_sh <- ggplot(sh_plot, aes(value)) +
  geom_histogram(bins = 35, fill = "grey60", colour = "white") +
  facet_wrap(~what) +
  labs(x = "Random intercept (mmHg)", y = "Participants")
save_fig(p_sh, "m11_shrinkage", width = 8, height = 3.5)

# ---- 5. simulation: fixed effects and non-normal random effects -------------
# Same design as sim_two_arm(), but random effects drawn from a normal, a
# right-skewed (centred, scaled chi-square with 2 df) or a bimodal
# distribution, always with the same variances.
draw_re <- function(n, shape) {
  z <- switch(shape,
              normal  = rnorm(n),
              skewed  = (rchisq(n, 2) - 2) / 2,
              bimodal = (sample(c(-1, 1), n, TRUE) * 0.9 +
                           rnorm(n, 0, sqrt(1 - 0.81))))
  z
}
sim_re_shape <- function(shape, n_per_arm = 60, times = c(0, 3, 6, 9, 12),
                         delta = -0.5, sd_int = 11, sd_slope = 0.4,
                         sd_resid = 6) {
  n <- 2 * n_per_arm
  trt <- rep(0:1, each = n_per_arm)
  u0 <- sd_int * draw_re(n, shape)
  u1 <- sd_slope * draw_re(n, shape)
  d <- expand.grid(time = times, id = seq_len(n))
  d$trt <- trt[d$id]
  d$y <- 150 - 0.2 * d$time + delta * d$trt * d$time +
    u0[d$id] + u1[d$id] * d$time + rnorm(nrow(d), 0, sd_resid)
  d$id <- factor(d$id)
  d
}
re_sim <- cache_or_run("m11_re_shape_sim", {
  set.seed(1102)
  n_rep <- 300
  bind_rows(lapply(c("normal", "skewed", "bimodal"), function(sh) {
    bind_rows(lapply(seq_len(n_rep), function(r) {
      d <- sim_re_shape(sh)
      m <- lmer(y ~ time * trt + (1 + time | id), data = d,
                control = lmerControl(optimizer = "bobyqa",
                                      calc.derivs = FALSE))
      s <- coef(summary(m))["time:trt", ]
      data.frame(shape = sh, est = s[["Estimate"]], se = s[["Std. Error"]],
                 sd_int_hat = attr(VarCorr(m)$id, "stddev")[[1]])
    }))
  }))
})
re_sum <- re_sim |>
  group_by(shape) |>
  summarise(n_sims = n(),
            mean_est = mean(est), bias = mean(est) - (-0.5),
            mcse_bias = sd(est) / sqrt(n()),
            emp_se = sd(est), mean_model_se = mean(se),
            coverage = mean(abs(est - (-0.5)) <= qnorm(0.975) * se),
            mcse_cov = sqrt(coverage * (1 - coverage) / n()))
save_table(re_sum, "m11_re_shape_sim")

# ---- 6. GLMM diagnostics with DHARMa ----------------------------------------
# CLEAR: logistic random-intercept model for moderate-or-severe itch
m_clear <- glmer(itch ~ week * arm + (1 | id), family = binomial,
                 data = clear)
save_model(m_clear, "m11_clear_glmm")
sim_clear <- simulateResiduals(m_clear, n = 1000, seed = 1103)
sim_clear_u <- simulateResiduals(m_clear, n = 1000, seed = 1103,
                                 simulateREs = "unconditional")
grp_clear_u <- recalculateResiduals(sim_clear_u, group = clear$id)
png(here("outputs", "figures", "m11_dharma_clear.png"), width = 1600,
    height = 800, res = 160)
plot(sim_clear)
dev.off()

# DIARY: Poisson versus negative binomial GLMM with exposure offset
m_pois <- glmmTMB(migraine_days ~ arm * post + offset(log(days)) + (1 | id),
                  family = poisson, data = diary)
m_nb   <- glmmTMB(migraine_days ~ arm * post + offset(log(days)) + (1 | id),
                  family = nbinom2, data = diary)
save_model(m_pois, "m11_diary_poisson")
save_model(m_nb, "m11_diary_nb")
sim_pois <- simulateResiduals(m_pois, n = 1000, seed = 1104)
sim_nb   <- simulateResiduals(m_nb, n = 1000, seed = 1104)
png(here("outputs", "figures", "m11_dharma_diary_poisson.png"), width = 1600,
    height = 800, res = 160)
plot(sim_pois)
dev.off()
png(here("outputs", "figures", "m11_dharma_diary_nb.png"), width = 1600,
    height = 800, res = 160)
plot(sim_nb)
dev.off()

dharma_tab <- data.frame(
  model = c("CLEAR logistic GLMM", "DIARY Poisson GLMM", "DIARY NB GLMM"),
  uniformity_p = c(testUniformity(sim_clear, plot = FALSE)$p.value,
                   testUniformity(sim_pois, plot = FALSE)$p.value,
                   testUniformity(sim_nb, plot = FALSE)$p.value),
  dispersion = c(testDispersion(sim_clear, plot = FALSE)$statistic,
                 testDispersion(sim_pois, plot = FALSE)$statistic,
                 testDispersion(sim_nb, plot = FALSE)$statistic),
  dispersion_p = c(testDispersion(sim_clear, plot = FALSE)$p.value,
                   testDispersion(sim_pois, plot = FALSE)$p.value,
                   testDispersion(sim_nb, plot = FALSE)$p.value),
  zero_ratio = c(NA,
                 testZeroInflation(sim_pois, plot = FALSE)$statistic,
                 testZeroInflation(sim_nb, plot = FALSE)$statistic),
  zero_p = c(NA,
             testZeroInflation(sim_pois, plot = FALSE)$p.value,
             testZeroInflation(sim_nb, plot = FALSE)$p.value))
save_table(dharma_tab, "m11_dharma_tests")

# Calibration check of the participant-level (grouped) DHARMa test: data
# simulated from a correctly specified logistic random-intercept model.
# Conditional simulation (DHARMa 0.5.0 default) versus unconditional.
grp_calib <- cache_or_run("m11_dharma_grouped_calibration", {
  set.seed(1105)
  bind_rows(lapply(1:20, function(r) {
    n <- 200; weeks <- c(0, 2, 4, 6, 8, 12)
    b <- rnorm(n, 0, 2); trt <- rep(0:1, each = n / 2)
    d <- expand.grid(week = weeks, id = seq_len(n))
    d$trt <- trt[d$id]
    d$itch <- rbinom(nrow(d), 1,
                     plogis(b[d$id] - 0.06 * d$week - 0.16 * d$trt * d$week))
    m <- glmer(itch ~ week * trt + (1 | id), family = binomial, data = d)
    out <- sapply(c("conditional", "unconditional"), function(how) {
      s <- simulateResiduals(m, n = 300, simulateREs = how, seed = r)
      g <- recalculateResiduals(s, group = d$id)
      testDispersion(g, plot = FALSE)$p.value
    })
    data.frame(rep = r, p_conditional = out[1], p_unconditional = out[2])
  }))
})
save_table(grp_calib, "m11_dharma_grouped_calibration")

# ---- 7. GEE: residuals by visit and a check of the mean model --------------
clear_gee <- clear |> mutate(id_f = factor(id)) |> arrange(id_f, week)
g_lin <- geeglm(itch ~ week * arm, id = id_f, family = binomial,
                corstr = "exchangeable", data = clear_gee)
g_cat <- geeglm(itch ~ week_f * arm, id = id_f, family = binomial,
                corstr = "exchangeable", data = clear_gee)
save_model(g_lin, "m11_clear_gee_linear")
clear_gee <- clear_gee |>
  mutate(fit_lin = fitted(g_lin)[, 1], fit_cat = fitted(g_cat)[, 1],
         pres = (itch - fit_lin) / sqrt(fit_lin * (1 - fit_lin)))
p_gee_res <- ggplot(clear_gee, aes(week_f, pres, fill = arm)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
  scale_fill_arm(name = NULL) +
  labs(x = "Week", y = "Pearson residual (linear-in-week GEE)")
gee_check <- clear_gee |>
  group_by(arm, week) |>
  summarise(observed = mean(itch), linear = mean(fit_lin),
            categorical = mean(fit_cat), n = n())
p_gee_mean <- ggplot(gee_check, aes(week, observed, colour = arm)) +
  geom_point(size = 2) +
  geom_line(aes(y = linear), linewidth = 0.7) +
  geom_line(aes(y = categorical), linetype = 2) +
  scale_colour_arm(name = NULL) +
  labs(x = "Week", y = "Proportion with moderate or severe itch")
save_fig(p_gee_res | p_gee_mean, "m11_gee_checks", width = 9, height = 4)
save_table(gee_check, "m11_gee_mean_check")

# ---- 8. convergence ---------------------------------------------------------
# (a) Poorly scaled model: time in days and age in years
walk_sc <- walk |> mutate(day = month * 30.44)
m_bad <- lmer(sbp ~ day * arm + age + (1 + day | id), data = walk_sc,
              control = lmer_ctrl)
m_good <- lmer(sbp ~ month * arm + I((age - 60) / 10) + (1 + month | id),
               data = walk_sc, control = lmer_ctrl)
# lme4 compares the Hessian-scaled gradient (or the raw gradient, if
# smaller) with its 0.002 tolerance. ML log-likelihoods are compared because
# the REML criterion changes when fixed-effect columns are rescaled.
conv_grad <- function(m) {
  d <- m@optinfo$derivs
  sc <- with(d, solve(chol(Hessian), gradient))
  max(pmin(abs(sc), abs(d$gradient)))
}
conv_tab <- data.frame(
  model = c("Time in days, age in years", "Time in months, age per 10 years"),
  scaled_grad = c(conv_grad(m_bad), conv_grad(m_good)),
  ep_slope_per_month = c(fixef(m_bad)["day:armEP"] * 30.44,
                         fixef(m_good)["month:armEP"]),
  loglik_ml = c(as.numeric(logLik(refitML(m_bad))),
                as.numeric(logLik(refitML(m_good)))))
save_table(conv_tab, "m11_convergence_scaling")

# (b) allFit: refit the primary model with every available optimiser.
# allFit() needs the control object written out in the call (it cannot
# evaluate the shared lmer_ctrl symbol), so the model is refitted here with an
# identical, explicit control argument.
m_pri_af <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk,
                 control = lmerControl(optimizer = "bobyqa"))
af <- allFit(m_pri_af, verbose = FALSE)
af_sum <- summary(af)
af_tab <- data.frame(
  optimiser = names(af_sum$llik),
  loglik = af_sum$llik,
  ep12 = af_sum$fixef[, "visit_f12:armEP"],
  sd_int = af_sum$sdcor[, 1], sd_slope = af_sum$sdcor[, 2],
  warning = vapply(af_sum$msgs, function(x) if (is.null(x)) "" else
    sub("\n.*", "", paste(x, collapse = "; ")), ""))
save_table(af_tab, "m11_allfit")

# The same model with lme4's default control settings: allFit() then shows
# a gradient warning from one optimiser alongside identical estimates.
m_pri_def <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk)
af_def <- summary(allFit(m_pri_def, verbose = FALSE))
save_table(data.frame(optimiser = names(af_def$llik), loglik = af_def$llik,
                      ep12 = af_def$fixef[, "visit_f12:armEP"],
                      warning = vapply(af_def$msgs, function(x)
                        !is.null(x), logical(1))),
           "m11_allfit_default")

# (c) Singular fit: random slope on month for 8 clinics
m_sing <- lmer(sbp ~ visit_f * arm + (1 + month | id) + (1 + month | clinic),
               data = walk, control = lmer_ctrl)

# ---- 9. separation in a sparse binary example ------------------------------
# A small trial with a rare outcome: no events after baseline in the active
# arm, so the arm-by-period interaction is (quasi-)completely separated.
set.seed(1106)
n_sep <- 40
sep <- expand.grid(visit = 0:3, id = seq_len(n_sep)) |>
  mutate(active = as.integer(id > n_sep / 2), post = as.integer(visit > 0))
b_sep <- rnorm(n_sep, 0, 0.5)
sep$event <- rbinom(nrow(sep), 1,
                    plogis(-2 + b_sep[sep$id] - 3 * sep$active * sep$post))
sep$event[sep$active == 1 & sep$post == 1] <- 0   # guarantee separation
sep$id <- factor(sep$id)
save_table(sep, "m11_separation_data")
sep_tab <- with(sep, table(group = ifelse(active == 1, "Active", "Control"),
                           period = ifelse(post == 1, "Post", "Baseline"),
                           event))
m_sep_glm <- glm(event ~ active * post, family = binomial, data = sep)

# Firth's penalised likelihood (heinze2002), written out by hand because the
# course does not load logistf/brglm2. Adds 0.5 * log|I(beta)| to the
# log-likelihood, which keeps estimates finite under separation.
firth_logit <- function(X, y, maxit = 100, tol = 1e-8) {
  b <- rep(0, ncol(X))
  for (it in seq_len(maxit)) {
    p <- as.vector(plogis(X %*% b))
    w <- p * (1 - p)
    XtWX <- crossprod(X, w * X)
    h <- w * rowSums((X %*% solve(XtWX)) * X)   # hat values
    U <- crossprod(X, y - p + h * (0.5 - p))     # modified score
    step <- solve(XtWX, U)
    b <- b + as.vector(step)
    if (max(abs(step)) < tol) break
  }
  p <- as.vector(plogis(X %*% b))
  XtWX <- crossprod(X, p * (1 - p) * X)
  data.frame(term = colnames(X), estimate = b,
             se = sqrt(diag(solve(XtWX))), iterations = it)
}
X_sep <- model.matrix(~ active * post, data = sep)
fit_firth <- firth_logit(X_sep, sep$event)
m_sep_glmer <- glmer(event ~ active * post + (1 | id), family = binomial,
                     data = sep)
m_sep_prior <- glmmTMB(event ~ active * post + (1 | id), family = binomial,
                       data = sep,
                       priors = data.frame(prior = "normal(0, 2.5)",
                                           class = "fixef"))
sep_compare <- data.frame(
  method = c("glm (ML)", "glmer (ML, Laplace)", "Firth (penalised)",
             "glmmTMB, normal(0, 2.5) prior"),
  estimate = c(coef(m_sep_glm)["active:post"],
               fixef(m_sep_glmer)["active:post"],
               fit_firth$estimate[4],
               fixef(m_sep_prior)$cond["active:post"]),
  se = c(sqrt(vcov(m_sep_glm)["active:post", "active:post"]),
         sqrt(vcov(m_sep_glmer)["active:post", "active:post"]),
         fit_firth$se[4],
         sqrt(vcov(m_sep_prior)$cond["active:post", "active:post"])))
save_table(sep_compare, "m11_separation")

# Sensitivity of a binary GLMM to the integration method (Laplace vs AGQ)
agq_tab <- bind_rows(lapply(c(1, 5, 10, 25), function(q) {
  m <- glmer(itch ~ week * arm + (1 | id), family = binomial, data = clear,
             nAGQ = q)
  data.frame(nAGQ = q, week_arm = fixef(m)["week:armActive"],
             se = sqrt(vcov(m)["week:armActive", "week:armActive"]),
             sd_int = attr(VarCorr(m)$id, "stddev")[[1]])
}))
save_table(agq_tab, "m11_agq_sensitivity")

message("Module 11 script finished.")
