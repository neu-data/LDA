# =============================================================================
# scripts/13a_ordinal.R
# Module 13 (ADVANCED), section 1: repeated ordinal outcomes.
#
# Purpose : fit a cumulative-logit mixed model (proportional odds with a
#           random intercept) to CLEAR itch severity, compare the estimates
#           with the true generating values, examine the proportional odds
#           assumption and compute predicted category probabilities.
# Inputs  : data/clear_long.csv (via load_clear()), data/clear_truth.rds
# Outputs : outputs/tables/m13_clmm_*.csv, outputs/figures/m13_clmm_*.png,
#           outputs/models/m13_clmm*.rds
# Module  : 13 Advanced extensions
# =============================================================================

source(here::here("R", "setup.R"))
suppressPackageStartupMessages(library(ordinal))   # advanced module only

clear <- load_clear() |>
  mutate(id = factor(id),
         age10 = (age - 45) / 10)      # centred where the simulation centred it
truth <- readRDS(here("data", "clear_truth.rds"))

# ---- description -----------------------------------------------------------
sev_by_week <- clear |>
  count(arm, week, severity_f) |>
  group_by(arm, week) |>
  mutate(prop = n / sum(n)) |>
  ungroup()
save_table(sev_by_week, "m13_clear_severity_by_week")

p_obs <- ggplot(sev_by_week, aes(factor(week), prop, fill = severity_f)) +
  geom_col(width = 0.75) +
  facet_wrap(~ arm) +
  scale_fill_manual(values = c("#F0F0F0", "#BDBDBD", "#6B6B6B", "#222222"),
                    name = "Itch severity") +
  labs(x = "Week", y = "Proportion of observed participants")
save_fig(p_obs, "m13_clear_severity_bars")

# ---- cumulative-logit mixed model --------------------------------------------
# nAGQ = 7: adaptive Gauss-Hermite quadrature, more accurate than the Laplace
# approximation (the default, nAGQ = 1) for a random-intercept SD near 2.
m_clmm <- clmm(severity_f ~ week * arm + (1 | id), data = clear, nAGQ = 7)
m_clmm_adj <- clmm(severity_f ~ week * arm + age10 + sex + (1 | id),
                   data = clear, nAGQ = 7)
summary(m_clmm)
summary(m_clmm_adj)
save_model(m_clmm, "m13_clmm")
save_model(m_clmm_adj, "m13_clmm_adj")

# ordinal parameterises logit P(Y <= k) = theta_k - eta, so a positive
# coefficient means a shift towards MORE severe categories, as in the
# simulation's latent-variable model.
est_tab <- function(m) {
  cf <- coef(summary(m))
  data.frame(term = rownames(cf), estimate = cf[, 1], se = cf[, 2],
             lower = cf[, 1] - 1.96 * cf[, 2], upper = cf[, 1] + 1.96 * cf[, 2],
             row.names = NULL)
}
cmp <- est_tab(m_clmm_adj) |>
  bind_rows(data.frame(term = "SD(random intercept)",
                       estimate = sqrt(as.numeric(VarCorr(m_clmm_adj)$id)))) |>
  mutate(truth = c(truth$cuts, truth$beta_week, 0, truth$beta_age,
                   truth$beta_female, truth$beta_trtweek, truth$sd_b))
cmp
save_table(cmp, "m13_clmm_vs_truth")

# ---- proportional odds check -------------------------------------------------
# Fit a binary random-intercept GLMM to each cumulative split (>= Mild,
# >= Moderate, = Severe). Under proportional odds the slopes should agree up
# to sampling error; only the intercepts should differ.
splits <- lapply(1:3, function(k) {
  g <- glmer(I(severity >= k) ~ week * arm + (1 | id), family = binomial,
             data = clear, nAGQ = 7)
  s <- coef(summary(g))
  data.frame(split = paste0("severity >= ", k), term = rownames(s),
             estimate = s[, 1], se = s[, 2], row.names = NULL)
})
po_check <- bind_rows(splits) |> filter(term != "(Intercept)")
po_check
save_table(po_check, "m13_po_check_splits")

# Formal comparison: let both week slopes differ by cut-point (a
# "nominal" effect, available in clmm2) and compare with a likelihood-ratio
# test. week_act carries the extra slope in the active arm, so the nominal
# week terms are not aliased with a location week term (with week * arm in
# the location part, one nominal term would be redundant).
clear_df <- as.data.frame(clear) |>
  mutate(week_act = week * (arm == "Active"))
m_po  <- clmm2(severity_f ~ week + arm + week_act, random = id,
               data = clear_df, Hess = TRUE, nAGQ = 7)
m_npo <- clmm2(severity_f ~ arm, nominal = ~ week + week_act, random = id,
               data = clear_df, Hess = TRUE, nAGQ = 7)
summary(m_npo)
po_lrt <- anova(m_po, m_npo)
po_lrt
save_table(as.data.frame(po_lrt), "m13_po_lrt")

# ---- predicted category probabilities ----------------------------------------
# Conditional (subject-specific): for a participant with random intercept 0.
# Marginal (population-averaged): average over the random-intercept
# distribution with Gauss-Hermite quadrature.
pred_probs <- function(m, newdata, marginal = FALSE, n_gh = 30) {
  th  <- m$Theta
  bet <- m$beta
  X   <- model.matrix(~ week * arm, newdata)[, -1, drop = FALSE]
  eta <- as.vector(X %*% bet[colnames(X)])
  sd_b <- sqrt(as.numeric(VarCorr(m)$id))
  if (marginal) {
    gh <- statmod_free_gh(n_gh)
    nodes <- sqrt(2) * sd_b * gh$x; w <- gh$w / sqrt(pi)
  } else { nodes <- 0; w <- 1 }
  cum <- sapply(seq_along(th), function(k)
    sapply(eta, function(e) sum(w * plogis(th[k] - e - nodes))))
  cum <- matrix(cum, nrow = length(eta))
  p <- cbind(cum, 1) - cbind(0, cum)
  colnames(p) <- levels(clear$severity_f)
  cbind(newdata, p)
}
# Gauss-Hermite nodes and weights (Golub-Welsch), avoiding extra packages
statmod_free_gh <- function(n) {
  i <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- sqrt(i / 2); J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}
grid <- expand.grid(week = c(0, 2, 4, 6, 8, 12),
                    arm = factor(c("Placebo", "Active"),
                                 levels = c("Placebo", "Active")))
pp <- bind_rows(
  pred_probs(m_clmm, grid) |> mutate(type = "Conditional (b = 0)"),
  pred_probs(m_clmm, grid, marginal = TRUE) |>
    mutate(type = "Marginal (averaged over b)")
) |>
  pivot_longer(None:Severe, names_to = "severity", values_to = "prob") |>
  mutate(severity = factor(severity, levels = levels(clear$severity_f)))
save_table(pp, "m13_clmm_pred_probs")

p_pp <- ggplot(pp, aes(week, prob, colour = arm)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  facet_grid(type ~ severity) +
  scale_colour_arm(name = NULL) +
  scale_x_continuous(breaks = c(0, 4, 8, 12)) +
  labs(x = "Week", y = "Predicted probability")
save_fig(p_pp, "m13_clmm_pred_probs", width = 8, height = 5)

message("13a done")
