# =============================================================================
# 07a_gee.R
# Purpose : Generalized estimating equations (GEE) for continuous (WALK),
#           binary (CLEAR) and count (DIARY) outcomes; working correlation
#           structures; model-based versus sandwich standard errors; the
#           Mancl-DeRouen bias-corrected sandwich; inverse-probability-
#           weighted GEE for dropout in CLEAR.
# Inputs  : data/walk_long.csv, data/clear_long.csv, data/diary_long.csv,
#           data/clear_truth.rds (via the R/ loaders)
# Outputs : outputs/tables/m07_*.csv, outputs/figures/m07_*.png,
#           outputs/models/m07_*.rds
# Module  : 7 (Generalized estimating equations)
# All data are simulated.
# =============================================================================

source(here::here("R", "setup.R"))

# geeglm() assumes that the rows of each cluster are contiguous and in time
# order, so every dataset is sorted by participant and time first.
walk  <- load_walk()  |> arrange(id, month)
clear <- load_clear() |> arrange(id, week)
diary <- load_diary() |> arrange(id, period)

# ---- helper: sandwich variance by hand, with optional Mancl-DeRouen ---------
# Reproduces geeglm's robust variance exactly when correction = "none"
# (checked below). With correction = "MD" each cluster's residual vector is
# inflated by (I - H_i)^{-1}, where H_i is the cluster "hat" matrix
# (Mancl and DeRouen 2001). Supports independence, exchangeable and ar1.
sandwich_gee <- function(fit, correction = c("none", "MD")) {
  correction <- match.arg(correction)
  X <- model.matrix(fit); y <- fit$y; mu <- fitted(fit)
  fam <- family(fit); eta <- fit$linear.predictors
  w <- fit$prior.weights; if (is.null(w)) w <- rep(1, length(y))
  alpha <- fit$geese$alpha
  cl <- split(seq_along(y), fit$id)
  p <- ncol(X); B <- matrix(0, p, p); parts <- vector("list", length(cl))
  for (k in seq_along(cl)) {
    r <- cl[[k]]; n <- length(r); R <- diag(n)
    if (fit$corstr == "exchangeable" && n > 1) { R[] <- alpha; diag(R) <- 1 }
    if (fit$corstr == "ar1" && n > 1) R <- alpha^abs(outer(1:n, 1:n, "-"))
    A <- diag(sqrt(fam$variance(mu[r]) / w[r]), n)
    Vinv <- solve(A %*% R %*% A)            # the scale phi cancels
    D <- fam$mu.eta(eta[r]) * X[r, , drop = FALSE]
    B <- B + t(D) %*% Vinv %*% D
    parts[[k]] <- list(D = D, Vinv = Vinv, e = y[r] - mu[r])
  }
  Binv <- solve(B); M <- matrix(0, p, p)
  for (pk in parts) {
    e <- pk$e
    if (correction == "MD") {
      H <- pk$D %*% Binv %*% t(pk$D) %*% pk$Vinv
      e <- solve(diag(nrow(H)) - H, e)
    }
    u <- t(pk$D) %*% pk$Vinv %*% e
    M <- M + u %*% t(u)
  }
  V <- Binv %*% M %*% Binv
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

# ---- helper: Gaussian GEE with an unstructured working correlation ----------
# geepack 1.3.13 terminates R (segmentation fault) when corstr =
# "unstructured" is combined with `waves` on incomplete clusters, so the
# identity-link case is coded directly: GLS given the working correlation,
# then a moment estimate of the correlation from Pearson residuals.
gee_unstructured <- function(formula, data, id, wave, maxit = 50, tol = 1e-8) {
  mf <- model.frame(formula, data); X <- model.matrix(formula, mf)
  y <- model.response(mf); idv <- data[[id]]; wv <- data[[wave]]
  K <- max(wv); cl <- split(seq_along(y), idv); p <- ncol(X)
  beta <- qr.solve(X, y)
  for (it in seq_len(maxit)) {
    e <- as.vector(y - X %*% beta)
    phi <- sum(e^2) / (length(y) - p); rs <- e / sqrt(mean(e^2))
    S <- matrix(0, K, K); Nn <- matrix(0, K, K)
    for (r in cl) {
      S[wv[r], wv[r]]  <- S[wv[r], wv[r]] + rs[r] %o% rs[r]
      Nn[wv[r], wv[r]] <- Nn[wv[r], wv[r]] + 1
    }
    Rw <- S / Nn; diag(Rw) <- 1
    B <- matrix(0, p, p); u <- rep(0, p)
    for (r in cl) {
      Ri <- solve(Rw[wv[r], wv[r], drop = FALSE])
      B <- B + t(X[r, , drop = FALSE]) %*% Ri %*% X[r, , drop = FALSE]
      u <- u + t(X[r, , drop = FALSE]) %*% Ri %*% y[r]
    }
    new <- as.vector(solve(B, u))
    done <- max(abs(new - beta)) < tol
    beta <- new
    if (done) break
  }
  e <- as.vector(y - X %*% beta); Binv <- solve(B); M <- matrix(0, p, p)
  for (r in cl) {
    s <- t(X[r, , drop = FALSE]) %*% solve(Rw[wv[r], wv[r], drop = FALSE], e[r])
    M <- M + s %*% t(s)
  }
  names(beta) <- colnames(X)
  list(coefficients = beta, R = Rw, phi = phi,
       vcov_robust = Binv %*% M %*% Binv, vcov_naive = phi * Binv,
       iterations = it)
}

# ---- 1. WALK: continuous outcome, four working correlations -----------------
walk_cs <- c("independence", "exchangeable", "ar1")
gee_walk <- lapply(walk_cs, function(s)
  geeglm(sbp ~ visit_f * arm, id = factor(id), waves = visit,
         data = walk, corstr = s))
names(gee_walk) <- walk_cs
gee_walk_un <- gee_unstructured(sbp ~ visit_f * arm, data = walk,
                                id = "id", wave = "visit")

# Check: on the complete teaching file (no missing visits), where geepack's
# unstructured option works, the hand-coded function reproduces geeglm.
walk_full <- load_walk("complete") |> arrange(id, visit)
chk_gp <- geeglm(sbp ~ month_sched * arm + age, id = factor(id),
                 data = walk_full, corstr = "unstructured")
chk_un <- gee_unstructured(sbp ~ month_sched * arm + age, data = walk_full,
                           id = "id", wave = "visit")
cat("Max |coef difference| =",
    signif(max(abs(coef(chk_gp) - chk_un$coefficients)), 3), "\n")

terms12 <- c("visit_f12:armSP", "visit_f12:armEP")
walk_wc <- bind_rows(
  lapply(walk_cs, function(s) {
    f <- gee_walk[[s]]
    data.frame(corstr = s, term = terms12, estimate = coef(f)[terms12],
               se_robust = sqrt(diag(vcov(f)))[terms12],
               se_naive = sqrt(diag(f$geese$vbeta.naiv))[match(terms12, names(coef(f)))])
  }),
  data.frame(corstr = "unstructured", term = terms12,
             estimate = gee_walk_un$coefficients[terms12],
             se_robust = sqrt(diag(gee_walk_un$vcov_robust))[match(terms12, names(gee_walk_un$coefficients))],
             se_naive = sqrt(diag(gee_walk_un$vcov_naive))[match(terms12, names(gee_walk_un$coefficients))])
)
rownames(walk_wc) <- NULL
print(walk_wc)
save_table(walk_wc, "m07_walk_working_correlation")

# QIC from geepack (smaller is better; see the notes for its limits)
walk_qic <- t(sapply(gee_walk, QIC))
print(walk_qic)
save_table(data.frame(corstr = rownames(walk_qic), walk_qic), "m07_walk_qic")

# Estimand: difference in mean change from baseline to month 12 (emmeans)
emm_walk <- emmeans(gee_walk$exchangeable, ~ visit_f | arm)
chg_walk <- contrast(emm_walk, list(change_m12 = c(-1, 0, 0, 0, 0, 1)))
diff_walk <- confint(contrast(chg_walk, "trt.vs.ctrl", by = NULL))
print(diff_walk)
save_table(diff_walk, "m07_walk_gee_diff_change")

# Comparison with the course primary linear mixed model (valid under MAR)
lmm_walk <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk,
                 control = lmer_ctrl)
diff_lmm <- confint(contrast(contrast(emmeans(lmm_walk, ~ visit_f | arm),
                                      list(change_m12 = c(-1, 0, 0, 0, 0, 1))),
                             "trt.vs.ctrl", by = NULL))
print(diff_lmm)

p_wc <- ggplot(walk_wc, aes(estimate, corstr)) +
  geom_pointrange(aes(xmin = estimate - 1.96 * se_robust,
                      xmax = estimate + 1.96 * se_robust)) +
  facet_wrap(~ term) +
  labs(x = "Difference in mean change at month 12 vs usual care (mmHg)",
       y = "Working correlation")
save_fig(p_wc, "m07_walk_working_correlation", width = 8, height = 3.5)

# ---- 2. CLEAR: binary outcome ----------------------------------------------
gee_clear <- geeglm(itch ~ week * arm, id = factor(id), waves = visit,
                    data = clear, family = binomial, corstr = "exchangeable")
print(summary(gee_clear))
se_tab <- data.frame(term = names(coef(gee_clear)),
                     estimate = coef(gee_clear),
                     se_naive = sqrt(diag(gee_clear$geese$vbeta.naiv)),
                     se_robust = sqrt(diag(vcov(gee_clear))),
                     se_md = sqrt(diag(sandwich_gee(gee_clear, "MD"))))
print(se_tab)
stopifnot(max(abs(sandwich_gee(gee_clear) - vcov(gee_clear))) < 1e-10)
save_table(se_tab, "m07_clear_se_comparison")
save_table(fixef_table(gee_clear, exponentiate = TRUE), "m07_clear_gee_or")

emm_clear <- emmeans(gee_clear, ~ arm | week, at = list(week = c(0, 12)),
                     type = "response")
print(emm_clear)
rd12 <- confint(contrast(regrid(emmeans(gee_clear, ~ arm | week,
                                        at = list(week = 12))),
                         "revpairwise"))
print(rd12)

pred_clear <- expand.grid(week = seq(0, 12, 0.5),
                          arm = factor(c("Placebo", "Active"),
                                       levels = c("Placebo", "Active")))
pred_clear$p <- predict(gee_clear, newdata = pred_clear, type = "response")
obs_clear <- clear |> group_by(arm, week) |> summarise(p = mean(itch), n = n())
p_clear <- ggplot(pred_clear, aes(week, p, colour = arm)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = obs_clear, aes(size = n), alpha = 0.7) +
  scale_colour_arm(name = NULL) + scale_size_area(max_size = 3, guide = "none") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 12)) +
  labs(x = "Week", y = "Proportion with moderate or severe itch")
save_fig(p_clear, "m07_clear_gee_fitted")

# ---- 3. DIARY: count outcome with exposure ---------------------------------
gee_diary <- geeglm(migraine_days ~ post * arm + offset(log(days)),
                    id = factor(id), waves = period + 1, data = diary,
                    family = poisson, corstr = "exchangeable")
print(summary(gee_diary))
save_table(fixef_table(gee_diary, exponentiate = TRUE), "m07_diary_gee_rr")
emm_diary <- emmeans(gee_diary, ~ arm | post, offset = log(28),
                     type = "response")
print(emm_diary)
save_table(as.data.frame(emm_diary), "m07_diary_gee_per28")

# ---- 4. CLEAR: inverse-probability-weighted GEE ----------------------------
# Dropout in CLEAR depends on the last observed severity (MAR). We model the
# probability of remaining in the study at each visit, given the previous
# observed severity and arm, and weight each observation by the inverse of
# the cumulative probability of having remained up to that visit.
build_ipw <- function(obs) {
  grid <- obs |> distinct(id, arm) |> crossing(visit = 1:6) |>
    left_join(select(obs, id, visit, severity), by = c("id", "visit")) |>
    group_by(id) |>
    mutate(last_visit = max(visit[!is.na(severity)]), last_sev = severity) |>
    tidyr::fill(last_sev, .direction = "down") |>
    mutate(prev_sev = lag(last_sev)) |> ungroup()
  risk <- grid |> filter(visit >= 2, visit <= pmin(last_visit + 1, 6)) |>
    mutate(stay = as.integer(visit <= last_visit))
  drop_model <- glm(stay ~ prev_sev + arm, family = binomial, data = risk)
  risk$p_stay <- fitted(drop_model)
  risk <- risk |> group_by(id) |> mutate(cum_p = cumprod(p_stay)) |> ungroup()
  out <- obs |>
    left_join(select(risk, id, visit, cum_p), by = c("id", "visit")) |>
    mutate(cum_p = ifelse(visit == 1, 1, cum_p), ipw = 1 / cum_p)
  list(data = out, model = drop_model, risk = risk)
}
ipw_clear <- build_ipw(clear)
print(summary(ipw_clear$model)$coefficients)
print(summary(ipw_clear$data$ipw))

gee_unw <- geeglm(itch ~ week * arm, id = factor(id), data = ipw_clear$data,
                  family = binomial, corstr = "independence")
gee_ipw <- geeglm(itch ~ week * arm, id = factor(id), data = ipw_clear$data,
                  family = binomial, corstr = "independence", weights = ipw)

# True population-average (marginal) coefficients, possible only because the
# data are simulated: integrate the true conditional model over the random
# intercept and the CLEAR participants' age and sex, then fit the marginal
# logistic model to the true probabilities.
truth_clear <- readRDS(here("data", "clear_truth.rds"))
true_marginal_clear <- function(truth, base, n_grid = 200) {
  z <- qnorm((seq_len(n_grid) - 0.5) / n_grid)
  g <- expand.grid(week = truth$weeks,
                   arm = factor(c("Placebo", "Active"),
                                levels = c("Placebo", "Active")))
  g$p <- sapply(seq_len(nrow(g)), function(k) {
    lp <- -truth$cuts[2] + truth$beta_week * g$week[k] +
      truth$beta_trtweek * (g$arm[k] == "Active") * g$week[k] +
      truth$beta_age * (base$age - 45) / 10 +
      truth$beta_female * (base$sex == "Female")
    mean(plogis(outer(lp, truth$sd_b * z, "+")))
  })
  g
}
base_clear <- distinct(clear, id, arm, age, sex)
true_p <- true_marginal_clear(truth_clear, base_clear)
true_marg_coef <- coef(suppressWarnings(
  glm(p ~ week * arm, family = binomial, data = true_p)))

ipw_tab <- data.frame(term = names(coef(gee_unw)),
                      truth_marginal = true_marg_coef,
                      unweighted = coef(gee_unw),
                      se_unweighted = sqrt(diag(vcov(gee_unw))),
                      ipw = coef(gee_ipw),
                      se_ipw = sqrt(diag(vcov(gee_ipw))))
print(ipw_tab)
save_table(ipw_tab, "m07_clear_ipw_comparison")

# ---- save models ------------------------------------------------------------
save_model(gee_walk$exchangeable, "m07_gee_walk_exch")
save_model(gee_clear, "m07_gee_clear_exch")
save_model(gee_diary, "m07_gee_diary_exch")
save_model(gee_ipw, "m07_gee_clear_ipw")

sessionInfo()
