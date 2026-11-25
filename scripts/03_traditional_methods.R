# =============================================================================
# scripts/03_traditional_methods.R
# Module 3. Statistical foundations and traditional approaches
#
# Purpose: regression review with an arm-by-time interaction, paired t-test,
#          repeated-measures ANOVA (aov with Error(), multivariate lm with
#          sphericity corrections), follow-up-only vs change-score vs ANCOVA
#          at month 12, a simulation of their precision, Lord's paradox, and
#          separate tests at each visit.
# Inputs:  data/walk_long.csv, data/walk_wide.csv (via load_walk())
# Outputs: outputs/figures/m03_*.png, outputs/tables/m03_*.csv,
#          outputs/models/m03_*.rds
# Module:  03-traditional-methods
# =============================================================================

source(here::here("R", "setup.R"))
walk <- load_walk()
wide <- load_walk("wide")

yv <- paste0("sbp_m", c(0, 1, 3, 6, 9, 12))   # wide SBP columns in visit order

# ---- 1. Regression review: OLS with an arm-by-time interaction --------------
# This model treats all 1,897 rows as independent. The point estimates are
# reasonable, but the standard errors ignore the correlation within each
# participant, which is why Module 4 replaces it with a mixed model.
m03_ols <- lm(sbp ~ month * arm, data = walk)
summary(m03_ols)
confint(m03_ols)
save_model(m03_ols, "m03_ols_month_arm")
save_table(fixef_table(m03_ols), "m03_ols_coefficients")

# Residuals from the same participant are strongly correlated across visits
walk$res_ols <- resid(m03_ols)
res_cor <- visit_correlation(walk, id, month_sched, res_ols)
round(res_cor, 2)
save_table(round(res_cor, 3), "m03_ols_residual_correlation")

p_res <- ggplot(walk, aes(fitted(m03_ols), res_ols)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_hline(yintercept = 0, colour = "#888888") +
  labs(x = "Fitted SBP (mmHg)", y = "Residual (mmHg)")
save_fig(p_res, "m03_ols_residuals")

# ---- 2. Paired t-test: baseline vs month 12 within the EP arm ---------------
ep <- wide |> filter(arm == "EP", !is.na(sbp_m0), !is.na(sbp_m12))
uc <- wide |> filter(arm == "UC", !is.na(sbp_m0), !is.na(sbp_m12))
tt_ep <- t.test(ep$sbp_m12, ep$sbp_m0, paired = TRUE)
tt_uc <- t.test(uc$sbp_m12, uc$sbp_m0, paired = TRUE)
tt_ep; tt_uc
save_table(data.frame(arm = c("EP", "UC"), n = c(nrow(ep), nrow(uc)),
                      mean_change = c(tt_ep$estimate, tt_uc$estimate),
                      lower = c(tt_ep$conf.int[1], tt_uc$conf.int[1]),
                      upper = c(tt_ep$conf.int[2], tt_uc$conf.int[2]),
                      p = c(tt_ep$p.value, tt_uc$p.value)),
           "m03_paired_t")

# ---- 3. Repeated-measures ANOVA on completers --------------------------------
cc <- wide[complete.cases(wide[, yv]), ]       # listwise deletion
nrow(cc)
lc <- cc |>
  select(id, arm, all_of(yv)) |>
  pivot_longer(all_of(yv), names_to = "visit", values_to = "sbp") |>
  mutate(visit = factor(visit, levels = yv), id = factor(id))

# (a) univariate approach: aov() with an Error() stratum for participants
m03_aov <- aov(sbp ~ arm * visit + Error(id / visit), data = lc)
summary(m03_aov)

# (b) multivariate lm on wide data: sphericity test and corrected p-values
Y <- as.matrix(cc[, yv])
m03_mlm <- lm(Y ~ arm, data = cc)
idata <- data.frame(visit = factor(yv, levels = yv))
rm_tab <- anova(m03_mlm, X = ~1, idata = idata, test = "Spherical")
rm_tab
mauchly.test(m03_mlm, X = ~1, idata = idata)
save_model(m03_mlm, "m03_rm_mlm")
save_table(cbind(term = rownames(rm_tab), as.data.frame(rm_tab)),
           "m03_rm_anova_spherical")

# Observed covariance matrix among completers: unequal variances and
# correlations that decline with the gap between visits break sphericity
S_cc <- cov(Y)
save_table(cbind(visit = yv, as.data.frame(round(S_cc, 1))), "m03_cov_completers")

cov_long <- as.data.frame(as.table(round(cov2cor(S_cc), 2)))
p_cor <- ggplot(cov_long, aes(Var1, Var2, fill = Freq)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 3.2) +
  scale_fill_gradient(low = "#F2F2F2", high = "#0072B2", name = "Correlation") +
  labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank())
save_fig(p_cor, "m03_observed_correlation", width = 6, height = 5)

# ---- 4. Month 12: follow-up only, change score, ANCOVA ----------------------
d12 <- wide |>
  filter(!is.na(sbp_m0), !is.na(sbp_m12)) |>
  mutate(chg = sbp_m12 - sbp_m0)
m03_fu    <- lm(sbp_m12 ~ arm, data = d12)
m03_chg   <- lm(chg ~ arm, data = d12)
m03_ancova <- lm(sbp_m12 ~ arm + sbp_m0, data = d12)
tab12 <- bind_rows(
  fixef_table(m03_fu)     |> mutate(analysis = "Follow-up only"),
  fixef_table(m03_chg)    |> mutate(analysis = "Change score"),
  fixef_table(m03_ancova) |> mutate(analysis = "ANCOVA")
) |> filter(term %in% c("armSP", "armEP"))
tab12
save_table(tab12, "m03_month12_three_analyses")
save_model(m03_ancova, "m03_ancova_m12")

# ---- 5. Simulation: precision of the three analyses -------------------------
# Randomised two-arm trial, true effect -5 mmHg, baseline-follow-up
# correlation rho. Variance of the estimator is proportional to 1 (follow-up),
# 2(1 - rho) (change) and 1 - rho^2 (ANCOVA).
sim_three <- function(rho, n = 100, effect = -5, sd = 13, reps = 1000) {
  S <- sd^2 * matrix(c(1, rho, rho, 1), 2)
  est <- replicate(reps, {
    y <- MASS::mvrnorm(2 * n, c(145, 145), S)
    trt <- rep(0:1, each = n)
    y[, 2] <- y[, 2] + effect * trt
    c(fu   = coef(lm(y[, 2] ~ trt))[2],
      chg  = coef(lm(I(y[, 2] - y[, 1]) ~ trt))[2],
      anc  = coef(lm(y[, 2] ~ trt + y[, 1]))[2])
  })
  data.frame(rho = rho, analysis = c("Follow-up only", "Change score", "ANCOVA"),
             mean_est = unname(rowMeans(est)), emp_se = unname(apply(est, 1, sd)))
}
set.seed(303)
sim_res <- bind_rows(lapply(c(0.2, 0.5, 0.7, 0.9), sim_three))
sim_res
save_table(sim_res, "m03_sim_three_analyses")

p_sim <- ggplot(sim_res, aes(rho, emp_se, colour = analysis)) +
  geom_line() + geom_point() +
  scale_colour_manual(values = c("Follow-up only" = "#6B6B6B",
                                 "Change score" = "#0072B2",
                                 "ANCOVA" = "#D55E00"), name = NULL) +
  labs(x = "Correlation between baseline and follow-up",
       y = "Empirical SE of the treatment effect (mmHg)")
save_fig(p_sim, "m03_sim_precision")

# ---- 6. Lord's paradox: two non-randomised groups ---------------------------
# Groups differ at baseline; within each group nothing changes on average.
set.seed(17)
n_l <- 500
grp <- rep(c("A", "B"), each = n_l)
true_level <- rnorm(2 * n_l, ifelse(grp == "A", 135, 150), 10)
lord <- data.frame(grp = grp,
                   y0 = true_level + rnorm(2 * n_l, 0, 7),
                   y1 = true_level + rnorm(2 * n_l, 0, 7))
coef(summary(lm(I(y1 - y0) ~ grp, data = lord)))["grpB", ]
coef(summary(lm(y1 ~ grp + y0, data = lord)))["grpB", ]

# ---- 7. Separate tests at each visit ----------------------------------------
per_visit <- bind_rows(lapply(yv, function(v) {
  d <- wide[!is.na(wide[[v]]), ]
  f <- fixef_table(lm(reformulate("arm", v), data = d))
  f |> filter(term != "(Intercept)") |> mutate(visit = v, n = nrow(d))
}))
per_visit
save_table(per_visit, "m03_per_visit_tests")

# Family-wise error under the null: 6 visits x 2 comparisons, correlated
# visits with the completers' covariance matrix, 120 per arm.
set.seed(3)
fwer <- mean(replicate(2000, {
  y <- MASS::mvrnorm(360, rep(145, 6), S_cc)
  arm <- factor(rep(c("UC", "SP", "EP"), each = 120), levels = c("UC", "SP", "EP"))
  p <- sapply(1:6, function(j) coef(summary(lm(y[, j] ~ arm)))[2:3, 4])
  any(p < 0.05)
}))
fwer
save_table(data.frame(tests = 12, fwer_sim = fwer,
                      fwer_independent = 1 - 0.95^12), "m03_fwer")

sessionInfo()
