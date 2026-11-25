# =============================================================================
# 00_refresher.R
# Module 0. Foundations refresher: R and regression
#
# Purpose : Cross-sectional regression refresher on the WALK baseline data:
#           dplyr summaries, linear regression, categorical predictors,
#           a centred age-by-sex interaction, logistic regression with odds
#           ratios and predicted probabilities, residual checks, confidence
#           versus prediction intervals, and a first look at why repeated
#           measurements break the independence assumption.
# Inputs  : data/walk_baseline.csv, data/walk_long.csv (via load_walk())
# Outputs : outputs/figures/m00_*.png, outputs/tables/m00_*.csv,
#           outputs/models/m00_*.rds
# Notes   : Uses the same model specifications as notes/00-refresher.qmd,
#           so the numbers agree. ALL DATA ARE SIMULATED.
# =============================================================================

source(here::here("R", "setup.R"))

# ---- data -------------------------------------------------------------------
# One row per participant. Centring age at 60 (per 10 years) and BMI at 29
# makes the intercept refer to a real kind of participant.
walk_base <- load_walk("baseline") |>
  mutate(
    age10  = (age - 60) / 10,
    bmi_c  = bmi - 29,
    hyp150 = as.integer(sbp0 >= 150)
  )

baseline_by_arm <- walk_base |>
  group_by(arm) |>
  summarise(n = n(), mean_age = mean(age),
            pct_women = 100 * mean(sex == "Female"),
            mean_sbp0 = mean(sbp0), sd_sbp0 = sd(sbp0),
            pct_150 = 100 * mean(hyp150))
print(baseline_by_arm)
save_table(baseline_by_arm, "m00_baseline_by_arm")

# ---- exploratory plots ------------------------------------------------------
p_hist <- ggplot(walk_base, aes(sbp0)) +
  geom_histogram(binwidth = 4, fill = "#9A9A9A", colour = "white") +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = "Baseline systolic blood pressure (mmHg)", y = "Participants")
save_fig(p_hist, "m00_sbp0_histogram")

p_scatter <- ggplot(walk_base, aes(bmi, sbp0, colour = sex)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  scale_colour_manual(values = c(Male = "#0072B2", Female = "#D55E00"),
                      name = NULL) +
  labs(x = "Body mass index (kg/m²)", y = "Baseline SBP (mmHg)")
save_fig(p_scatter, "m00_sbp0_bmi_sex")

# ---- linear regression ------------------------------------------------------
m_lin <- lm(sbp0 ~ age10 + sex + bmi_c, data = walk_base)
print(summary(m_lin))
lin_tab <- fixef_table(m_lin)
lin_tab$lower <- confint(m_lin)[, 1]   # t-based intervals, as in confint()
lin_tab$upper <- confint(m_lin)[, 2]
print(lin_tab)
save_table(lin_tab, "m00_linear_coefficients")
save_model(m_lin, "m00_linear")

# Uncentred version: same slopes, intercept refers to age 0 and BMI 0
m_raw <- lm(sbp0 ~ age + sex + bmi, data = walk_base)
print(round(coef(m_raw), 3))

# Categorical predictor with three levels; arm was randomised, so baseline
# differences between arms can only be chance
m_arm <- lm(sbp0 ~ arm, data = walk_base)
print(coef(summary(m_arm)))
print(anova(m_arm))
m_arm_ep <- lm(sbp0 ~ relevel(arm, ref = "EP"), data = walk_base)
print(round(coef(m_arm_ep), 2))

# ---- residual checks --------------------------------------------------------
diag_df <- data.frame(fitted = fitted(m_lin), resid = resid(m_lin),
                      std_resid = rstandard(m_lin))
p_rvf <- ggplot(diag_df, aes(fitted, resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "#0072B2") +
  labs(x = "Fitted SBP (mmHg)", y = "Residual (mmHg)")
p_qq <- ggplot(diag_df, aes(sample = std_resid)) +
  stat_qq(alpha = 0.5) +
  stat_qq_line(colour = "#0072B2") +
  labs(x = "Theoretical normal quantile", y = "Standardised residual")
save_fig(p_rvf + p_qq, "m00_residual_checks", width = 9, height = 4)

# ---- interaction: age by sex, age centred at 60 -------------------------------
m_int <- lm(sbp0 ~ age10 * sex + bmi_c, data = walk_base)
print(summary(m_int))
print(anova(m_lin, m_int))
print(emtrends(m_int, ~ sex, var = "age10"))
int_tab <- fixef_table(m_int)
int_tab$lower <- confint(m_int)[, 1]
int_tab$upper <- confint(m_int)[, 2]
save_table(int_tab, "m00_interaction_coefficients")
save_model(m_int, "m00_interaction")

# Without centring the sex "main effect" refers to age 0
m_int_raw <- lm(sbp0 ~ age * sex + bmi_c, data = walk_base)
print(coef(summary(m_int_raw))["sexFemale", ])

nd <- expand.grid(age = 40:80,
                  sex = factor(c("Male", "Female"), levels = c("Male", "Female")),
                  bmi_c = 0) |>
  mutate(age10 = (age - 60) / 10)
nd <- cbind(nd, predict(m_int, newdata = nd, interval = "confidence"))
p_int <- ggplot(nd, aes(age, fit, colour = sex, fill = sex)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c(Male = "#0072B2", Female = "#D55E00"),
                      name = NULL) +
  scale_fill_manual(values = c(Male = "#0072B2", Female = "#D55E00"),
                    name = NULL) +
  labs(x = "Age (years)", y = "Fitted mean SBP (mmHg)")
save_fig(p_int, "m00_interaction_age_sex")

p_int_check <- walk_base |>
  mutate(resid = resid(m_int)) |>
  ggplot(aes(age, resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "#0072B2") +
  facet_wrap(~ sex) +
  labs(x = "Age (years)", y = "Residual (mmHg)")
save_fig(p_int_check, "m00_interaction_residuals")

# ---- logistic regression: odds of baseline SBP >= 150 --------------------------
m_log <- glm(hyp150 ~ age10 + sex + bmi_c, family = binomial, data = walk_base)
print(summary(m_log))
# Profile-likelihood intervals on the log-odds scale, then exponentiated
or_tab <- data.frame(term = names(coef(m_log)),
                     exp(cbind(OR = coef(m_log), confint(m_log))),
                     p = coef(summary(m_log))[, 4], row.names = NULL)
print(or_tab)
save_table(or_tab, "m00_logistic_odds_ratios")
save_model(m_log, "m00_logistic")

# Predicted probabilities: interval on the link scale, then back-transform
profiles <- data.frame(
  label = c("Man, 60, BMI 29", "Woman, 60, BMI 29", "Woman, 70, BMI 35"),
  age10 = c(0, 0, 1),
  sex   = factor(c("Male", "Female", "Female"), levels = c("Male", "Female")),
  bmi_c = c(0, 0, 6)
)
lp <- predict(m_log, newdata = profiles, type = "link", se.fit = TRUE)
profiles <- profiles |>
  mutate(prob  = plogis(lp$fit),
         lower = plogis(lp$fit - 1.96 * lp$se.fit),
         upper = plogis(lp$fit + 1.96 * lp$se.fit))
print(profiles)
save_table(profiles, "m00_logistic_predicted_probabilities")

calib <- walk_base |>
  mutate(p_hat = fitted(m_log), group = ntile(p_hat, 5)) |>
  group_by(group) |>
  summarise(pred = mean(p_hat), obs = mean(hyp150), n = n()) |>
  mutate(se = sqrt(obs * (1 - obs) / n))
p_calib <- ggplot(calib, aes(pred, obs)) +
  geom_abline(linetype = "dashed") +
  geom_errorbar(aes(ymin = obs - 1.96 * se, ymax = obs + 1.96 * se),
                width = 0.01) +
  geom_point(size = 2) +
  coord_equal(xlim = c(0.1, 0.8), ylim = c(0.1, 0.8)) +
  labs(x = "Mean predicted probability", y = "Observed proportion")
save_fig(p_calib, "m00_logistic_calibration", width = 5, height = 5)

# ---- confidence vs prediction intervals ---------------------------------------
new_men <- data.frame(age10 = 0,
                      sex = factor("Male", levels = c("Male", "Female")),
                      bmi_c = c(0, 6))
cipi <- data.frame(bmi = c(29, 35),
                   predict(m_lin, newdata = new_men, interval = "confidence"),
                   predict(m_lin, newdata = new_men,
                           interval = "prediction")[, 2:3])
names(cipi) <- c("bmi", "fit", "ci_lower", "ci_upper", "pi_lower", "pi_upper")
print(cipi)
save_table(cipi, "m00_ci_vs_pi")

grid <- data.frame(age10 = 0,
                   sex = factor("Male", levels = c("Male", "Female")),
                   bmi_c = seq(19, 40, by = 0.5) - 29)
band <- cbind(grid,
              as.data.frame(predict(m_lin, grid, interval = "confidence")),
              setNames(as.data.frame(predict(m_lin, grid,
                                             interval = "prediction"))[, 2:3],
                       c("plo", "phi"))) |>
  mutate(bmi = bmi_c + 29)
p_band <- ggplot(band, aes(bmi)) +
  geom_ribbon(aes(ymin = plo, ymax = phi), fill = "#0072B2", alpha = 0.12) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#0072B2", alpha = 0.35) +
  geom_line(aes(y = fit), colour = "#0072B2", linewidth = 0.9) +
  geom_point(data = filter(walk_base, sex == "Male", between(age, 55, 65)),
             aes(bmi, sbp0), alpha = 0.6) +
  labs(x = "Body mass index (kg/m²)", y = "SBP (mmHg)")
save_fig(p_band, "m00_ci_vs_pi")

# ---- why repeated measurements break independence -----------------------------
# Fitting the same regression to every visit treats 1,897 correlated rows as
# independent. The BMI standard error shrinks, but BMI was measured once per
# person, so the extra precision is spurious.
walk <- load_walk()
m_once <- lm(sbp0 ~ bmi, data = walk_base)
m_all  <- lm(sbp ~ bmi, data = walk)
naive_tab <- rbind(baseline_only = coef(summary(m_once))["bmi", 1:2],
                   all_visits    = coef(summary(m_all))["bmi", 1:2])
print(naive_tab)
save_table(data.frame(model = rownames(naive_tab), naive_tab,
                      row.names = NULL), "m00_naive_all_visits")

m1_sbp  <- walk |> filter(month_sched == 1) |> select(id, sbp_m1 = sbp)
pairs01 <- inner_join(walk_base, m1_sbp, by = "id")
cat("Correlation of SBP at baseline and month 1:",
    round(cor(pairs01$sbp0, pairs01$sbp_m1), 2), "\n")
