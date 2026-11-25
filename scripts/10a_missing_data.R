# =============================================================================
# 10a_missing_data.R
# Module 10. Missing data and dropout: analyses of the WALK and KIDNEY data.
#
# Purpose
#   (1) Describe the missingness pattern in WALK and what predicts dropout.
#   (2) Compare the mixed model on all observed data (valid under MAR) with
#       complete-case and LOCF analyses, and with the complete-data answer
#       (available only because the data are simulated).
#   (3) Multiple imputation: (a) wide format, per visit, separately by arm;
#       (b) long format, multilevel (2l.norm). Pool with Rubin's rules.
#   (4) Delta-adjustment sensitivity analysis and tipping-point display.
#   (5) Informative visit timing in KIDNEY: naive window means vs a mixed
#       model.
#
# Inputs   data/walk_long.csv, walk_wide.csv, walk_baseline.csv,
#          walk_complete.csv, kidney_long.csv, kidney_events.csv,
#          kidney_truth.rds
# Outputs  outputs/figures/m10_*.png, outputs/tables/m10_*.csv,
#          outputs/models/m10_*.rds
# Simulations are in scripts/10b_missing_data_simulations.R.
# =============================================================================

source(here::here("R", "setup.R"))

walk     <- load_walk()
walk_b   <- load_walk("baseline")
walk_w   <- load_walk("wide")
walk_c   <- load_walk("complete")

# ---- describing the missingness ---------------------------------------------
pattern <- walk_c |>
  arrange(id, visit) |>
  group_by(id, arm) |>
  summarise(pattern = paste(ifelse(observed, "O", "."), collapse = ""),
            n_obs = sum(observed),
            last = max(visit[observed])) |>
  ungroup() |>
  mutate(type = case_when(
    n_obs == 6 ~ "Complete",
    n_obs == last ~ "Dropout only",
    last == 6 ~ "Intermittent only",
    TRUE ~ "Intermittent and dropout"))
save_table(count(pattern, arm, type), "m10_pattern_types")
save_table(count(pattern, pattern, sort = TRUE), "m10_patterns")

p_map <- plot_missing_map(walk_c, id, month_sched, observed,
                          order_by = as.numeric(observed) * visit) +
  facet_wrap(~ arm, scales = "free_y",
             labeller = as_labeller(walk_arm_labels))
save_fig(p_map, "m10_missing_map", width = 8, height = 5)

# Who drops out? One row per participant-interval at risk: was the
# participant gone by the next scheduled visit, given the last observed SBP?
risk <- walk_c |>
  arrange(id, visit) |>
  group_by(id) |>
  mutate(last_seen = cummax(ifelse(observed, visit, 0)),
         last_sbp = sapply(seq_along(visit), function(k)
           sbp[max(which(observed[1:k]))])) |>
  ungroup() |>
  left_join(walk_b |> select(id, last_month_sched), by = "id") |>
  group_by(id) |>
  mutate(prev_sbp = lag(last_sbp), prev_sched = lag(month_sched)) |>
  ungroup() |>
  filter(visit > 1, prev_sched <= last_month_sched) |>
  mutate(drop_now = as.integer(prev_sched == last_month_sched),
         prev_sbp10 = (prev_sbp - 150) / 10)
m_drop <- glm(drop_now ~ prev_sbp10 + arm, family = binomial, data = risk)
save_table(fixef_table(m_drop, exponentiate = TRUE), "m10_dropout_model")

# ---- MAR-valid mixed model vs complete cases vs LOCF -----------------------
diff_change <- function(m) {
  e  <- emmeans(m, ~ visit_f | arm, at = list(visit_f = c("0", "12")))
  ch <- contrast(e, list(change = c(-1, 1)), by = "arm")
  summary(contrast(ch, "trt.vs.ctrl", by = NULL), infer = TRUE,
          adjust = "none") |> as.data.frame()
}

m_lmm  <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk,
               control = lmer_ctrl)
completers <- walk_b$id[!walk_b$dropout]
m_cc   <- lmer(sbp ~ visit_f * arm + (1 + month | id),
               data = filter(walk, id %in% completers), control = lmer_ctrl)
m_full <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk_c,
               control = lmer_ctrl)

locf <- walk |>
  group_by(id, arm) |>
  summarise(base = sbp[visit == 1], last = sbp[which.max(visit)]) |>
  ungroup()
m_locf <- lm(I(last - base) ~ arm, data = locf)
locf_tab <- fixef_table(m_locf) |> filter(term != "(Intercept)")

comp <- bind_rows(
  diff_change(m_full) |> mutate(analysis = "Complete data (simulation only)"),
  diff_change(m_lmm)  |> mutate(analysis = "Mixed model, all observed data"),
  diff_change(m_cc)   |> mutate(analysis = "Mixed model, completers only"),
  data.frame(contrast = c("change SP - change UC", "change EP - change UC"),
             estimate = locf_tab$estimate, SE = locf_tab$se,
             lower.CL = locf_tab$lower, upper.CL = locf_tab$upper,
             analysis = "LOCF, change to last value")
) |> select(analysis, contrast, estimate, SE, lower.CL, upper.CL)
save_table(comp, "m10_analysis_comparison")
save_model(m_lmm, "m10_primary_lmm")

# ---- Rubin's rules --------------------------------------------------------
# est, se: vectors of estimates and standard errors from m imputations.
rubin <- function(est, se, n = Inf, k = 1) {
  p <- mice::pool.scalar(Q = est, U = se^2, n = n, k = k)
  q <- qt(0.975, p$df)
  data.frame(estimate = p$qbar, within = p$ubar, between = p$b,
             total = p$t, se = sqrt(p$t), df = p$df, fmi = p$fmi,
             lower = p$qbar - q * sqrt(p$t), upper = p$qbar + q * sqrt(p$t))
}

# ---- (a) wide-format MI, per visit, separately by arm ----------------------
# Imputing within each arm lets every visit-to-visit relation differ by arm,
# so the imputation model is at least as rich as the analysis model.
walk_w <- walk_w |> mutate(age10 = (age - 60) / 10, bmi_c = bmi - 29)
imp_vars <- c("age10", "sex", "bmi_c",
              paste0("sbp_m", c(0, 1, 3, 6, 9, 12)))
m_imp <- 50
imp_wide <- lapply(levels(walk_w$arm), function(a) {
  mice(walk_w |> filter(arm == a) |> select(id, all_of(imp_vars)),
       m = m_imp, method = "norm", maxit = 10, seed = 2026,
       predictorMatrix = {
         pm <- make.predictorMatrix(walk_w |> filter(arm == a) |>
                                      select(id, all_of(imp_vars)))
         pm[, "id"] <- 0
         pm
       },
       printFlag = FALSE)
})
names(imp_wide) <- levels(walk_w$arm)
save_model(imp_wide, "m10_mi_wide")

# k-th completed dataset (all arms stacked)
completed_wide <- function(k) {
  bind_rows(lapply(names(imp_wide), function(a)
    mice::complete(imp_wide[[a]], k) |> mutate(arm = a))) |>
    mutate(arm = factor(arm, levels = c("UC", "SP", "EP")))
}

# Analysis of each completed dataset: change to month 12 regressed on arm.
# With no missing values left, this targets the same estimand as the
# primary mixed model (difference in mean change at month 12).
post_drop <- walk_b |> select(id, last_month_sched)
delta_fit <- function(delta, arms_shift = c("SP", "EP")) {
  res <- lapply(seq_len(m_imp), function(k) {
    d <- completed_wide(k) |>
      left_join(post_drop, by = "id") |>
      # delta adjustment: post-dropout month-12 values in the programme
      # arms are made delta mmHg higher than MAR predicts
      mutate(sbp_m12 = sbp_m12 + ifelse(arm %in% arms_shift &
                                          last_month_sched < 12, delta, 0))
    f <- lm(I(sbp_m12 - sbp_m0) ~ arm, data = d)
    data.frame(term = c("SP", "EP"), est = coef(f)[2:3],
               se = sqrt(diag(vcov(f)))[2:3])
  }) |> bind_rows()
  bind_rows(lapply(c("SP", "EP"), function(tm) {
    r <- filter(res, term == tm)
    rubin(r$est, r$se, n = nrow(walk_w), k = 3) |> mutate(term = tm)
  })) |> mutate(delta = delta)
}
mi_wide_mar <- delta_fit(0)
save_table(mi_wide_mar, "m10_mi_wide_mar")

# ---- (b) long-format multilevel MI (2l.norm) -------------------------------
# Rows for every scheduled visit; the missing SBP values are imputed from a
# linear mixed model with visit x arm fixed effects and a random intercept
# and slope on time. For missed visits the actual time is unknown, so the
# scheduled month stands in for it.
long_grid <- walk_c |>
  mutate(sbp = ifelse(observed, sbp, NA),
         month = ifelse(observed, month, month_sched),
         idn = as.integer(factor(id)))
X <- model.matrix(~ visit_f * arm, long_grid)[, -1]
colnames(X) <- make.names(colnames(X))
ld <- data.frame(idn = long_grid$idn, sbp = long_grid$sbp, X,
                 month = long_grid$month, cons = 1)
meth <- make.method(ld); meth[] <- ""; meth["sbp"] <- "2l.norm"
pred <- make.predictorMatrix(ld); pred[, ] <- 0
pred["sbp", ] <- 1                        # fixed effects
pred["sbp", "sbp"] <- 0
pred["sbp", "idn"] <- -2                  # cluster identifier
pred["sbp", c("cons", "month")] <- 2      # random intercept and slope
m_imp_long <- 20
imp_long <- mice(ld, m = m_imp_long, method = meth, predictorMatrix = pred,
                 maxit = 1, seed = 2027, printFlag = FALSE)
save_model(imp_long, "m10_mi_long")

res_long <- lapply(seq_len(m_imp_long), function(k) {
  d <- long_grid |> mutate(sbp = mice::complete(imp_long, k)$sbp)
  m <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = d,
            control = lmer_ctrl)
  s <- coef(summary(m))
  data.frame(term = c("SP", "EP"),
             est = s[c("visit_f12:armSP", "visit_f12:armEP"), 1],
             se  = s[c("visit_f12:armSP", "visit_f12:armEP"), 2])
}) |> bind_rows()
mi_long <- bind_rows(lapply(c("SP", "EP"), function(tm) {
  r <- filter(res_long, term == tm)
  rubin(r$est, r$se, n = nrow(walk_w), k = 3) |> mutate(term = tm)
}))
save_table(mi_long, "m10_mi_long")

# ---- delta adjustment and tipping point ----------------------------------
deltas <- seq(0, 30, by = 2)
tip <- bind_rows(lapply(deltas, delta_fit))
save_table(tip, "m10_tipping_point")

p_tip <- ggplot(tip, aes(delta, estimate, colour = term, fill = term)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point() +
  scale_colour_arm(name = NULL, labels = c(SP = "SP vs UC",
                                           EP = "EP vs UC")) +
  scale_fill_arm(name = NULL, labels = c(SP = "SP vs UC",
                                         EP = "EP vs UC")) +
  labs(x = "δ: post-dropout SBP in programme arms above MAR prediction (mmHg)",
       y = "Difference in mean change at month 12 (mmHg)")
save_fig(p_tip, "m10_tipping_point")

# ---- informative visit timing in KIDNEY --------------------------------
kidney <- load_kidney() |> mutate(age10 = (age - 60) / 10)
k_ev   <- load_kidney("events") |> mutate(age10 = (age - 60) / 10)
k_tr   <- readRDS(here("data", "kidney_truth.rds"))

kidney <- kidney |>
  group_by(id) |> arrange(years, .by_group = TRUE) |>
  mutate(next_gap = lead(years) - years) |> ungroup()
gap_tab <- kidney |>
  filter(!is.na(next_gap)) |>
  mutate(last_low = ifelse(egfr < 30, "eGFR < 30", "eGFR >= 30")) |>
  group_by(last_low) |>
  summarise(n_gaps = n(), mean_gap_years = mean(next_gap))
save_table(gap_tab, "m10_kidney_gaps")

m_kid <- lmer(egfr ~ years + I(years^2) + diabetes + diabetes:years +
                age10 + (1 + years | id), data = kidney, control = lmer_ctrl)
save_model(m_kid, "m10_kidney_lmm")
k_tab <- data.frame(year = 1:4) |>
  rowwise() |>
  mutate(n_obs = sum(abs(kidney$years - year) <= 0.25),
         naive = mean(kidney$egfr[abs(kidney$years - year) <= 0.25]),
         lmm = mean(predict(m_kid, newdata = mutate(k_ev, years = year),
                            re.form = NA)),
         truth = mean(k_tr$b0 + (k_tr$b1 + k_tr$b1_dm * k_ev$diabetes) *
                        year + k_tr$b2 * year^2 +
                        k_tr$b_age * k_ev$age10)) |>
  ungroup()
save_table(k_tab, "m10_kidney_window_means")

print(count(pattern, type)); print(fixef_table(m_drop, exponentiate = TRUE))
print(comp); print(mi_wide_mar); print(mi_long); print(tip); print(gap_tab)
print(k_tab)
