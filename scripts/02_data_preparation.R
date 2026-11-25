# =============================================================================
# 02_data_preparation.R
# Module 2. Preparing and exploring longitudinal data in R
#
# Purpose  Reshape, join, clean and explore the WALK data. The script takes the
#          messy visit export (walk_visits_raw.csv) through a documented
#          cleaning workflow, checks the result against the analysis file
#          (walk_long.csv), describes missing visits and participant flow, and
#          draws the standard exploratory plots. It ends by comparing observed
#          means with complete-data means and mixed-model estimates.
# Inputs   data/walk_visits_raw.csv, data/walk_baseline.csv, data/walk_long.csv,
#          data/walk_wide.csv, data/walk_complete.csv (teaching only)
# Outputs  outputs/tables/m02_*.csv, outputs/figures/m02_*.png,
#          outputs/models/m02_primary_lmm.rds
# Module   02-data-preparation
#
# THE DATA ARE SIMULATED. Results illustrate methods, not empirical findings.
# =============================================================================

source(here::here("R", "setup.R"))

baseline <- load_walk("baseline")
walk     <- load_walk("long")
wide     <- load_walk("wide")

sched_months <- c(0, 1, 3, 6, 9, 12)
days_per_month <- 30.44   # WALK convention: 365.25/12 rounded to 2 decimals

# ---- 1. reshaping: wide to long and back -----------------------------------
# Wide data hold one column per scheduled visit. Pivoting longer gives one row
# per participant and scheduled visit; missed visits become explicit NA rows.
long_from_wide <- wide |>
  pivot_longer(cols = starts_with("sbp_m"),
               names_to = "month_sched", names_prefix = "sbp_m",
               names_transform = list(month_sched = as.numeric),
               values_to = "sbp")

# Dropping the NA rows should reproduce the attended visits in walk_long
check_reshape <- long_from_wide |>
  filter(!is.na(sbp)) |>
  full_join(walk |> select(id, month_sched, sbp_long = sbp),
            by = c("id", "month_sched"))
stopifnot(nrow(check_reshape) == nrow(walk),
          all(check_reshape$sbp == check_reshape$sbp_long))

# And back to wide: the result must equal the original wide file
wide_again <- walk |>
  select(id, clinic, arm, age, sex, bmi, month_sched, sbp) |>
  pivot_wider(names_from = month_sched, values_from = sbp,
              names_prefix = "sbp_m")
stopifnot(isTRUE(all.equal(as.data.frame(wide_again), as.data.frame(wide),
                           check.attributes = FALSE)))

# ---- 2. joining baseline data to visit data --------------------------------
# Join only the variables not already in the visit file, state the expected
# relationship, and check rows and unmatched keys in both directions.
walk_joined <- walk |>
  left_join(baseline |> select(id, rand_date, sbp0, steps0, dropout,
                               last_month_sched),
            by = "id", relationship = "many-to-one")
stopifnot(nrow(walk_joined) == nrow(walk),
          nrow(anti_join(walk, baseline, by = "id")) == 0,
          nrow(anti_join(baseline, walk, by = "id")) == 0)

# ---- 3. cleaning the raw visit export --------------------------------------
# Read every column as text first so that nothing is silently converted.
raw <- read_csv(here("data", "walk_visits_raw.csv"),
                col_types = cols(.default = col_character()))

# 3a. standardise identifiers and visit labels
step1 <- raw |>
  mutate(
    id          = toupper(trimws(participant)),
    month_sched = as.numeric(sub("^M", "", visit_label)),
    visit       = match(month_sched, sched_months)
  )
stopifnot(all(grepl("^W[0-9]{3}$", step1$id)), !anyNA(step1$visit))

# 3b. parse dates and convert text to numbers
step2 <- step1 |>
  mutate(
    visit_date = as.Date(visit_date, format = "%Y-%m-%d"),
    sbp        = parse_double(sbp_mmhg),
    steps      = parse_double(steps_k)
  ) |>
  select(id, visit, month_sched, visit_date, sbp, steps)
stopifnot(!anyNA(step2$visit_date), !anyNA(step2$steps))

# 3c. exact duplicates: identical in every column, so one copy carries no
#     information. Done after standardising IDs, or w004/W004 copies are missed.
n_exact_dup <- sum(duplicated(step2))
step3 <- distinct(step2)

# 3d. conflicting duplicates: same participant and visit, different values.
#     Rule (fixed before looking at outcomes or arms): keep the earliest-dated
#     record and flag the visit. A rule chosen after seeing which version
#     helps which arm would be a source of bias.
conflicts <- step3 |>
  add_count(id, visit, name = "n_records") |>
  filter(n_records > 1) |>
  arrange(id, visit, visit_date)

step4 <- step3 |>
  group_by(id, visit) |>
  mutate(conflict_flag = n() > 1) |>
  arrange(visit_date, .by_group = TRUE) |>
  slice(1) |>
  ungroup()
stopifnot(nrow(count(step4, id, visit) |> filter(n > 1)) == 0)

# 3e. attach randomisation date and derive actual time in months
step5 <- step4 |>
  left_join(baseline |> select(id, rand_date),
            by = "id", relationship = "many-to-one", unmatched = "error") |>
  mutate(
    days        = as.numeric(visit_date - rand_date),
    month       = round(days / days_per_month, 2),
    target_date = rand_date + round(month_sched * days_per_month),
    dev_days    = as.numeric(visit_date - target_date)
  )
stopifnot(all(step5$days >= 0), all(step5$month[step5$visit == 1] == 0))

# 3f. plausibility: SBP outside 70-250 mmHg is implausible for this population.
#     The implausible value is set to missing and logged as a data query.
sbp_limits <- c(70, 250)
step6 <- step5 |>
  mutate(
    sbp_flag = case_when(
      is.na(sbp) ~ "blank",
      sbp < sbp_limits[1] | sbp > sbp_limits[2] ~ "implausible",
      TRUE ~ "ok"),
    sbp_raw = sbp,
    sbp     = if_else(sbp_flag == "implausible", NA_real_, sbp)
  )

# 3g. visit windows from the (hypothetical) WALK protocol, in days
window_days <- c(`0` = 0, `1` = 7, `3` = 14, `6` = 28, `9` = 28, `12` = 28)
walk_clean <- step6 |>
  mutate(window = window_days[as.character(month_sched)],
         out_of_window = abs(dev_days) > window) |>
  select(id, visit, month_sched, visit_date, rand_date, month, dev_days,
         out_of_window, sbp, sbp_raw, sbp_flag, steps, conflict_flag) |>
  arrange(id, visit)

query_log <- bind_rows(
  conflicts |> transmute(id, visit, issue = "conflicting duplicate",
                         detail = paste0(visit_date, ": SBP ", sbp),
                         action = "kept earliest-dated record"),
  walk_clean |> filter(sbp_flag != "ok") |>
    transmute(id, visit, issue = paste("SBP", sbp_flag),
              detail = paste0("raw value: ", sbp_raw),
              action = "set to missing; query sent to clinic")
)

cleaning_log <- tibble(
  step = c("Rows read", "Exact duplicate rows removed",
           "Rows after removing exact duplicates",
           "Visits with conflicting records", "Rows after resolving conflicts",
           "Blank SBP values", "Implausible SBP values set to missing",
           "Out-of-window visits (kept, flagged)"),
  n = c(nrow(raw), n_exact_dup, nrow(step3),
        n_distinct(paste(conflicts$id, conflicts$visit)), nrow(walk_clean),
        sum(walk_clean$sbp_flag == "blank"),
        sum(walk_clean$sbp_flag == "implausible"),
        sum(walk_clean$out_of_window))
)
print(cleaning_log)
save_table(cleaning_log, "m02_cleaning_log")
save_table(query_log, "m02_query_log")

window_table <- walk_clean |>
  filter(visit > 1) |>
  group_by(month_sched) |>
  summarise(window_days = first(window_days[as.character(month_sched)]),
            n = n(), min_dev = min(dev_days), max_dev = max(dev_days),
            sd_dev = round(sd(dev_days), 1),
            n_out = sum(out_of_window))
save_table(window_table, "m02_visit_windows")

# 3h. agreement with the analysis file walk_long.csv
agreement <- walk_clean |>
  full_join(walk |> select(id, visit, month_ref = month, sbp_ref = sbp,
                           steps_ref = steps),
            by = c("id", "visit"))
agreement_summary <- tibble(
  check = c("Rows in cleaned data", "Rows in walk_long",
            "Keys only in cleaned data", "Keys only in walk_long",
            "SBP values that differ (incl. missing)",
            "Steps values that differ",
            "Rows where month differs", "Largest month difference"),
  value = c(nrow(walk_clean), nrow(walk),
            sum(is.na(agreement$month_ref)), sum(is.na(agreement$month)),
            sum(is.na(agreement$sbp) | is.na(agreement$sbp_ref) |
                  agreement$sbp != agreement$sbp_ref),
            sum(agreement$steps != agreement$steps_ref),
            sum(abs(agreement$month - agreement$month_ref) > 0.001),
            max(abs(agreement$month - agreement$month_ref)))
)
print(agreement_summary)
save_table(agreement_summary, "m02_agreement_walk_long")

# ---- 4. missing visits: the complete id x scheduled-visit grid --------------
# Built from the design (randomised participants x protocol visits), not from
# walk_complete.csv, which no real study would have.
grid <- baseline |>
  select(id, arm) |>
  crossing(visit = seq_along(sched_months)) |>
  mutate(month_sched = sched_months[visit]) |>
  left_join(walk |> select(id, visit, sbp), by = c("id", "visit")) |>
  mutate(observed = !is.na(sbp)) |>
  group_by(id) |>
  mutate(last_seen = max(visit[observed]),
         miss_type = case_when(observed ~ "observed",
                               visit < last_seen ~ "intermittent",
                               TRUE ~ "after dropout")) |>
  ungroup()

missing_summary <- grid |>
  count(arm, miss_type) |>
  pivot_wider(names_from = miss_type, values_from = n, values_fill = 0)
save_table(missing_summary, "m02_missing_types")

# ---- 5. summary tables --------------------------------------------------------
mean_sd <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"),
          mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}
baseline_table <- baseline |>
  group_by(arm) |>
  summarise(
    n = as.character(n()),
    `Age, years` = mean_sd(age),
    `Female, n (%)` = sprintf("%d (%.0f%%)", sum(sex == "Female"),
                              100 * mean(sex == "Female")),
    `BMI, kg/m2` = mean_sd(bmi),
    `SBP, mmHg` = mean_sd(sbp0),
    `Steps, 1000/day` = mean_sd(steps0)
  ) |>
  pivot_longer(-arm, names_to = "Characteristic") |>
  pivot_wider(names_from = arm, values_from = value)
print(baseline_table)
save_table(baseline_table, "m02_baseline_table")

visit_summary <- walk |>
  group_by(arm, month_sched) |>
  summarise(n = n(), mean = mean(sbp), sd = sd(sbp),
            median = median(sbp), mean_month = mean(month))
save_table(visit_summary, "m02_visit_summary")

flow <- visit_flow(walk, id, month_sched, arm) |>
  left_join(count(baseline, arm, name = "randomised"), by = "arm") |>
  mutate(cell = sprintf("%d (%.0f%%)", n_seen, 100 * n_seen / randomised)) |>
  select(arm, month_sched, cell) |>
  pivot_wider(names_from = arm, values_from = cell)
print(flow)
save_table(flow, "m02_visit_flow")

# ---- 6. exploratory plots -----------------------------------------------------
walk_lab <- walk |>
  mutate(arm_lab = factor(walk_arm_labels[as.character(arm)],
                          levels = walk_arm_labels))

p_spag_all <- plot_spaghetti(walk_lab, month, sbp, id, group = arm_lab,
                             alpha = 0.15) +
  labs(x = "Months since randomisation", y = "SBP (mmHg)")
save_fig(p_spag_all, "m02_spaghetti_all", width = 9, height = 4)

p_spag_sub <- plot_spaghetti(walk_lab, month, sbp, id, group = arm_lab,
                             n_show = 30, seed = 2, alpha = 0.7) +
  labs(x = "Months since randomisation", y = "SBP (mmHg)")
save_fig(p_spag_sub, "m02_spaghetti_subset", width = 9, height = 4)

set.seed(4)
panel_ids <- baseline |> group_by(arm) |> slice_sample(n = 4) |> pull(id)
p_panels <- walk_lab |>
  filter(id %in% panel_ids) |>
  ggplot(aes(month, sbp, colour = arm)) +
  geom_line() + geom_point(size = 1.2) +
  facet_wrap(~ id, ncol = 4) +
  scale_colour_arm(name = NULL, labels = walk_arm_labels) +
  labs(x = "Months since randomisation", y = "SBP (mmHg)")
save_fig(p_panels, "m02_individual_panels", width = 9, height = 6)

p_means <- plot_observed_means(walk_lab, month_sched, sbp, arm_lab) +
  scale_x_continuous(breaks = sched_months) +
  labs(x = "Scheduled month", y = "Observed mean SBP (mmHg), 95% CI")
save_fig(p_means, "m02_observed_means")

p_map <- plot_missing_map(grid, id, visit, observed,
                          order_by = last_seen) +
  facet_wrap(~ arm, scales = "free_y",
             labeller = as_labeller(walk_arm_labels)) +
  scale_x_discrete(labels = sched_months) +
  labs(x = "Scheduled month")
save_fig(p_map, "m02_missing_map", width = 9, height = 5)

dropout_by_arm <- baseline |>
  count(arm, last_month_sched) |>
  group_by(arm) |>
  mutate(pct = 100 * n / sum(n))
save_table(dropout_by_arm, "m02_dropout_by_arm")

retention <- grid |>
  group_by(arm, month_sched) |>
  summarise(pct_still_in = 100 * mean(last_seen >= visit))
p_retention <- ggplot(retention, aes(month_sched, pct_still_in,
                                     colour = arm)) +
  geom_step(linewidth = 0.8) + geom_point() +
  scale_colour_arm(name = NULL, labels = walk_arm_labels) +
  scale_x_continuous(breaks = sched_months) +
  labs(x = "Scheduled month",
       y = "Participants not yet dropped out (%)")
save_fig(p_retention, "m02_retention")

# Were participants who dropped out different just before they left? Compare
# each participant's SBP with the mean of their arm at the same visit.
dropout_last_sbp <- walk |>
  left_join(baseline |> select(id, last_month_sched), by = "id",
            relationship = "many-to-one") |>
  filter(month_sched < 12) |>
  group_by(arm, month_sched) |>
  mutate(dev_from_mean = sbp - mean(sbp)) |>
  ungroup() |>
  mutate(status = if_else(month_sched == last_month_sched,
                          "Last visit before dropout", "Seen again later")) |>
  group_by(status) |>
  summarise(n = n(), mean_dev = mean(dev_from_mean))
save_table(dropout_last_sbp, "m02_dropout_last_sbp")

pattern_means <- walk |>
  left_join(baseline |> select(id, last_month_sched), by = "id",
            relationship = "many-to-one") |>
  mutate(pattern = factor(last_month_sched, levels = sched_months,
                          labels = paste("Last seen month", sched_months))) |>
  group_by(arm, pattern, month_sched) |>
  summarise(mean_sbp = mean(sbp), n = n())
p_pattern <- pattern_means |>
  ggplot(aes(month_sched, mean_sbp, colour = pattern)) +
  geom_line() + geom_point(aes(size = n)) +
  facet_wrap(~ arm, labeller = as_labeller(walk_arm_labels)) +
  scale_size_area(max_size = 3, guide = "none") +
  scale_colour_manual(values = c("#D55E00", "#E69F00", "#56B4E9", "#009E73",
                                 "#CC79A7", "#222222"), name = NULL) +
  scale_x_continuous(breaks = sched_months) +
  labs(x = "Scheduled month", y = "Mean observed SBP (mmHg)")
save_fig(p_pattern, "m02_pattern_means", width = 9, height = 4.5)

# ---- 7. observed means versus model-estimated means ---------------------------
# Teaching only: walk_complete.csv shows what the means would have been with
# no missing data. A real study never has this file.
complete <- load_walk("complete")
m_primary <- lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk,
                  control = lmer_ctrl)
save_model(m_primary, "m02_primary_lmm")

em_visits <- as.data.frame(emmeans(m_primary, ~ visit_f | arm))
means_compare <- complete |>
  group_by(arm, visit_f) |>
  summarise(complete_mean = mean(sbp),
            observed_mean = mean(sbp[observed]),
            n_observed = sum(observed)) |>
  left_join(em_visits |> select(arm, visit_f, model_mean = emmean),
            by = c("arm", "visit_f")) |>
  mutate(observed_minus_complete = observed_mean - complete_mean,
         model_minus_complete = model_mean - complete_mean) |>
  ungroup()
print(means_compare)
save_table(means_compare, "m02_observed_vs_model_means")

p_compare <- means_compare |>
  select(arm, visit_f, `All scheduled visits (teaching only)` = complete_mean,
         `Observed visits only` = observed_mean,
         `Mixed model, observed data` = model_mean) |>
  pivot_longer(-c(arm, visit_f), names_to = "source", values_to = "mean") |>
  mutate(month = as.numeric(as.character(visit_f))) |>
  ggplot(aes(month, mean, linetype = source, shape = source)) +
  geom_line() + geom_point() +
  facet_wrap(~ arm, labeller = as_labeller(walk_arm_labels)) +
  scale_x_continuous(breaks = sched_months) +
  labs(x = "Scheduled month", y = "Mean SBP (mmHg)", linetype = NULL,
       shape = NULL)
save_fig(p_compare, "m02_observed_vs_model", width = 9, height = 4.5)

message("Module 2 script finished: ", nrow(walk_clean), " cleaned visits.")
