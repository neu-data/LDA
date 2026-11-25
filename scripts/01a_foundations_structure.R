# =============================================================================
# 01a_foundations_structure.R
# Module 1. Foundations of longitudinal data
#
# Purpose : Describe the structure of the WALK and KIDNEY data (subjects,
#           observations, visits, scheduled and actual times, balance and
#           follow-up), show within-subject dependence (observed correlation
#           across visits) and estimate between- and within-subject variation,
#           and compare naive OLS with a mixed model on WALK.
# Inputs  : data/walk_long.csv, data/kidney_long.csv (via load_walk(),
#           load_kidney())
# Outputs : outputs/figures/m01_*.png, outputs/tables/m01_*.csv,
#           outputs/models/m01_*.rds
# Notes   : Same model specifications as notes/01-foundations.qmd.
#           ALL DATA ARE SIMULATED.
# =============================================================================

source(here::here("R", "setup.R"))
walk   <- load_walk()
kidney <- load_kidney()

# ---- design schematic: four kinds of data ---------------------------------------
# Rows are units, columns are time; filled tiles are measurements.
schematic <- bind_rows(
  expand.grid(unit = 1:6, time = 1) |> mutate(design = "Cross-sectional"),
  expand.grid(unit = 1:6, time = 1:5) |> mutate(design = "Longitudinal (balanced)"),
  data.frame(unit = rep(1:6, c(5, 3, 5, 2, 4, 5)),
             time = c(1:5, 1:3, 1:5, 1:2, c(1, 2, 4, 5), 1:5)) |>
    mutate(design = "Longitudinal (unbalanced)"),
  expand.grid(unit = 1:6, time = 1:3) |>
    mutate(cluster = rep(c("Clinic A", "Clinic B"), length.out = n()),
           design = "Clustered (no time order)")
) |>
  mutate(design = factor(design, levels = c("Cross-sectional",
                                            "Longitudinal (balanced)",
                                            "Longitudinal (unbalanced)",
                                            "Clustered (no time order)")))
p_schematic <- ggplot(schematic, aes(time, factor(unit))) +
  geom_tile(fill = "#0072B2", colour = "white", linewidth = 1, alpha = 0.8) +
  facet_wrap(~ design, nrow = 1) +
  scale_x_continuous(breaks = 1:5) +
  labs(x = "Measurement occasion (or member of cluster)", y = "Unit") +
  theme(panel.grid = element_blank())
save_fig(p_schematic, "m01_design_schematic", width = 10, height = 3.2)

# ---- structure checks -----------------------------------------------------------------
print(check_long(walk, id, month_sched, "sbp"))
obs_per_person <- walk |> count(id, name = "n_obs") |> count(n_obs, name = "participants")
print(obs_per_person)
save_table(obs_per_person, "m01_walk_obs_per_participant")

flow <- visit_flow(walk, id, month_sched, arm)
print(flow)
save_table(flow, "m01_walk_visit_flow")

timing <- walk |>
  group_by(month_sched) |>
  summarise(n = n(), min_actual = min(month), max_actual = max(month),
            sd_deviation = sd(month - month_sched))
print(timing)
save_table(timing, "m01_walk_timing")

p_timing <- ggplot(walk, aes(month, visit_f)) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.25, size = 0.8) +
  geom_point(data = distinct(walk, month_sched, visit_f),
             aes(month_sched, visit_f), shape = 124, size = 8,
             colour = "#D55E00") +
  scale_x_continuous(breaks = c(0, 1, 3, 6, 9, 12)) +
  labs(x = "Actual time of visit (months since randomisation)",
       y = "Scheduled visit (month)")
save_fig(p_timing, "m01_walk_actual_vs_scheduled")

print(check_long(kidney, id, years, "egfr"))
kidney_fu <- kidney |>
  group_by(id) |>
  summarise(n_visits = n(), follow_up = max(years))
print(summary(kidney_fu))

set.seed(4)
ids_w <- sample(unique(walk$id), 15)
ids_k <- sample(unique(kidney$id), 15)
timing_df <- bind_rows(
  walk |> filter(id %in% ids_w) |>
    transmute(study = "WALK (scheduled visits)", id, t = month / 12),
  kidney |> filter(id %in% ids_k) |>
    transmute(study = "KIDNEY (routine clinic visits)", id, t = years)
) |>
  mutate(study = factor(study, levels = c("WALK (scheduled visits)",
                                          "KIDNEY (routine clinic visits)")))
p_visits <- ggplot(timing_df, aes(t, reorder(id, t, max))) +
  geom_line(aes(group = id), colour = "#D9D9D9") +
  geom_point(size = 1.4) +
  facet_wrap(~ study, scales = "free") +
  labs(x = "Years since entry", y = "Participant") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_fig(p_visits, "m01_visit_patterns", width = 9, height = 4)

# Visits are more frequent after a low eGFR: informative visit timing
kidney_gaps <- kidney |>
  arrange(id, years) |>
  group_by(id) |>
  mutate(gap_months = 12 * (years - lag(years)), prev_egfr = lag(egfr)) |>
  ungroup() |>
  filter(!is.na(gap_months)) |>
  group_by(previous_egfr_below_30 = prev_egfr < 30) |>
  summarise(n_gaps = n(), mean_gap_months = mean(gap_months))
print(kidney_gaps)
save_table(kidney_gaps, "m01_kidney_visit_gaps")

# ---- within-subject dependence ------------------------------------------------------------
cor_walk <- visit_correlation(walk, id, month_sched, sbp)
print(round(cor_walk, 2))
save_table(data.frame(month = rownames(cor_walk), round(cor_walk, 3)),
           "m01_walk_correlation")

cor_long <- as.data.frame(as.table(cor_walk)) |>
  setNames(c("row", "col", "r")) |>
  mutate(row = factor(row, levels = rev(colnames(cor_walk))))
p_cor <- ggplot(cor_long, aes(col, row, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.5) +
  scale_fill_gradient(low = "#F2F2F2", high = "#0072B2", limits = c(0.5, 1),
                      name = "Correlation") +
  labs(x = "Scheduled month", y = "Scheduled month") +
  theme(panel.grid = element_blank())
save_fig(p_cor, "m01_walk_correlation_heatmap", width = 6, height = 5)

# ---- between- and within-person spread ----------------------------------------------------------
set.seed(21)
ids_b <- sample(unique(walk$id[walk$visit == 6]), 25)
spread_df <- walk |>
  filter(id %in% ids_b) |>
  group_by(id) |>
  mutate(person_mean = mean(sbp)) |>
  ungroup()
p_spread <- ggplot(spread_df, aes(reorder(id, person_mean), sbp)) +
  geom_hline(yintercept = mean(walk$sbp), linetype = "dashed") +
  geom_point(alpha = 0.5, size = 1.4) +
  geom_point(aes(y = person_mean), shape = 95, size = 8, colour = "#D55E00") +
  labs(x = "Participant (ordered by own mean SBP)", y = "SBP (mmHg)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p_spread, "m01_between_within_spread")

decomp <- walk |>
  group_by(id) |>
  mutate(person_mean = mean(sbp), deviation = sbp - person_mean) |>
  ungroup()
descriptive <- c(total_var = var(walk$sbp),
                 var_person_means = var(distinct(decomp, id, person_mean)$person_mean),
                 var_within_deviations = var(decomp$deviation))
print(descriptive)

# ---- variance components from random-intercept models ------------------------------------------------
m_ri0 <- lmer(sbp ~ 1 + (1 | id), data = walk, control = lmer_ctrl)
m_ri  <- lmer(sbp ~ visit_f * arm + (1 | id), data = walk, control = lmer_ctrl)
vc <- rbind(
  data.frame(model = "Intercept only", varcomp_table(m_ri0), icc = icc_lmer(m_ri0)),
  data.frame(model = "Visit x arm means", varcomp_table(m_ri), icc = icc_lmer(m_ri))
)
print(vc)
save_table(vc, "m01_variance_components")
save_model(m_ri, "m01_random_intercept")

# ---- naive OLS versus the linear-time mixed model on WALK ------------------------------------------------
m_ols <- lm(sbp ~ month * arm, data = walk)
m_lin <- lmer(sbp ~ month * arm + (1 + month | id), data = walk,
              control = lmer_ctrl)
se_compare <- data.frame(term = names(coef(m_ols)),
                         est_ols = coef(m_ols),
                         se_ols  = coef(summary(m_ols))[, 2],
                         est_lmm = fixef(m_lin),
                         se_lmm  = coef(summary(m_lin))[, 2],
                         row.names = NULL) |>
  mutate(se_ratio = se_ols / se_lmm)
print(se_compare)
save_table(se_compare, "m01_walk_ols_vs_lmm")
save_model(m_lin, "m01_linear_time")

# ---- preview: the primary estimand -------------------------------------------------------------------------
# Difference in mean change from baseline to month 12 versus usual care, from
# the course's primary model. Modules 4 and 6 develop this model fully.
m_primary <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk,
                  control = lmer_ctrl)
primary <- fixef_table(m_primary) |>
  filter(term %in% c("visit_f12:armSP", "visit_f12:armEP")) |>
  mutate(true_value = c(-4.6, -8.1))
print(primary)
save_table(primary, "m01_primary_estimand_preview")
save_model(m_primary, "m01_primary")

p_means <- plot_observed_means(walk, month_sched, sbp, arm) +
  scale_x_continuous(breaks = c(0, 1, 3, 6, 9, 12)) +
  labs(x = "Scheduled month", y = "Observed mean SBP (mmHg)")
save_fig(p_means, "m01_observed_means")

p_spag <- plot_spaghetti(walk, month, sbp, id, group = arm, n_show = 45) +
  labs(x = "Months since randomisation", y = "SBP (mmHg)")
save_fig(p_spag, "m01_spaghetti", width = 9, height = 4)
