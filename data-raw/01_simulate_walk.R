# =============================================================================
# 01_simulate_walk.R
# Simulates the WALK study: the core continuous-outcome dataset of the course.
#
# THIS IS SIMULATED DATA. The WALK study does not exist. All "results" obtained
# from these data are illustrations of methods, not empirical findings.
#
# Design (hypothetical)
#   A three-arm randomised trial of a community walking programme for adults
#   with raised blood pressure. 8 clinics, 360 participants (120 per arm).
#   Arms: Usual care (UC), Standard programme (SP), Enhanced programme (EP).
#   Outcome: systolic blood pressure (SBP, mmHg).
#   Scheduled visits: months 0 (baseline), 1, 3, 6, 9, 12.
#   Actual visit times vary around the schedule (irregular timing).
#   Time-varying covariate: average daily steps (thousands) at each visit.
#   Missingness: intermittent missed visits (completely at random) and
#   dropout that depends on the most recent observed SBP (missing at random).
#
# Outputs (written to data/)
#   walk_baseline.csv     one row per participant (baseline covariates)
#   walk_long.csv         one row per attended visit (analysis dataset)
#   walk_wide.csv         one row per participant, SBP by scheduled visit
#   walk_complete.csv     all scheduled visits before missingness (TEACHING
#                         ONLY - never available in a real study)
#   walk_visits_raw.csv   messy visit-level export for Module 2 cleaning
#   walk_truth.rds        true generating parameters
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(here)

# Seed 3 was chosen from 15 candidate seeds (20261026, 1-14) because its
# estimated treatment effects are typical of the sampling distribution rather
# than unusually close to or far from the true values listed in walk_truth.rds.
set.seed(3)

# ---- design ----------------------------------------------------------------
n_clinic     <- 8
n_per_arm    <- 120
arms         <- c("UC", "SP", "EP")
arm_labels   <- c(UC = "Usual care", SP = "Standard programme",
                  EP = "Enhanced programme")
sched_months <- c(0, 1, 3, 6, 9, 12)
n            <- n_per_arm * length(arms)

# ---- true parameters -------------------------------------------------------
truth <- list(
  intercept      = 150,     # mean SBP at baseline for reference covariates
  age_per10      = 3.0,     # mmHg per 10 years (age centred at 60)
  female         = -4.0,    # women vs men
  bmi            = 0.6,     # mmHg per kg/m2 (BMI centred at 29)
  # treatment trajectory: mean change = -A * (1 - exp(-t / tau))
  A              = c(UC = 2, SP = 6, EP = 9),
  tau            = 2,       # months; how quickly the change develops
  sd_clinic      = 3.0,     # between-clinic SD of intercept
  sd_b0          = 11.0,    # between-participant SD of intercept
  sd_b1          = 0.35,    # SD of participant-specific linear slopes (/month)
  cor_b0b1       = -0.25,
  sd_serial      = 5.0,     # SD of the serially correlated residual process
  phi_serial     = 3.0,     # months; exp(-|t - t'| / phi) correlation
  sd_meas        = 4.0,     # independent measurement error SD
  # steps (thousands/day): person mean, arm-induced increase, fluctuation
  steps_mean     = 6.0,
  steps_sd_pers  = 1.8,
  steps_inc      = c(UC = 0.2, SP = 1.2, EP = 2.0),
  steps_sd_occ   = 1.2,
  beta_steps_between = -1.5, # mmHg per 1000 steps difference between people
  beta_steps_within  = -0.6, # mmHg per 1000 steps change within a person
  # missingness
  p_miss_visit   = 0.05,    # intermittent missed visit, MCAR
  drop_int       = -3.6,    # dropout logit intercept (per interval)
  drop_sbp       = 0.045,   # dropout logit per mmHg of last observed SBP - 150
  drop_ep        = 0.4      # extra dropout in the enhanced programme
)

# ---- participants ----------------------------------------------------------
clinic_eff <- rnorm(n_clinic, 0, truth$sd_clinic)

base <- tibble(
  id     = sprintf("W%03d", seq_len(n)),
  clinic = sample(rep(seq_len(n_clinic), length.out = n)),
  arm    = factor(rep(arms, each = n_per_arm), levels = arms),
  age    = round(pmin(pmax(rnorm(n, 60, 9), 40), 80)),
  sex    = factor(ifelse(runif(n) < 0.52, "Female", "Male"),
                  levels = c("Male", "Female")),
  bmi    = round(pmin(pmax(rnorm(n, 29, 4.5), 19), 45), 1)
) |>
  # randomisation within clinic is approximated by shuffling arm order
  group_by(clinic) |>
  mutate(arm = sample(arm)) |>
  ungroup() |>
  arrange(id)

# participant random effects (intercept, slope)
S <- matrix(c(truth$sd_b0^2,
              truth$cor_b0b1 * truth$sd_b0 * truth$sd_b1,
              truth$cor_b0b1 * truth$sd_b0 * truth$sd_b1,
              truth$sd_b1^2), 2)
b <- MASS::mvrnorm(n, c(0, 0), S)

# habitual steps: lower in people with higher BMI (a between-person confounder)
base <- base |>
  mutate(
    b0         = b[, 1],
    b1         = b[, 2],
    steps_pers = truth$steps_mean - 0.12 * (bmi - 29) +
                 rnorm(n, 0, truth$steps_sd_pers),
    rand_date  = as.Date("2024-01-08") + sample(0:330, n, replace = TRUE)
  )

# ---- visits (complete data) ------------------------------------------------
# actual time = scheduled time + jitter; jitter grows with scheduled time
jitter_sd <- c(0, 0.12, 0.25, 0.35, 0.40, 0.45)

visits <- base |>
  select(id, arm, clinic, age, sex, bmi, b0, b1, steps_pers) |>
  tidyr::crossing(visit = seq_along(sched_months)) |>
  mutate(
    month_sched = sched_months[visit],
    month       = pmax(month_sched + rnorm(n(), 0, jitter_sd[visit]), 0),
    month       = ifelse(visit == 1, 0, month),
    month       = round(month, 2)
  )

# serially correlated residuals (exponential correlation in actual time)
serial <- visits |>
  group_by(id) |>
  group_modify(function(d, key) {
    R <- exp(-abs(outer(d$month, d$month, "-")) / truth$phi_serial)
    e <- as.vector(t(chol(truth$sd_serial^2 * R)) %*% rnorm(nrow(d)))
    tibble(visit = d$visit, e_serial = e)
  }) |>
  ungroup()

visits <- visits |>
  left_join(serial, by = c("id", "visit")) |>
  mutate(
    steps_inc = truth$steps_inc[as.character(arm)] * (1 - exp(-month / 1)),
    steps     = round(pmax(steps_pers + steps_inc +
                             rnorm(n(), 0, truth$steps_sd_occ), 0.5), 1),
    mu_trt    = -truth$A[as.character(arm)] * (1 - exp(-month / truth$tau)),
    sbp_true  = truth$intercept +
      truth$age_per10 * (age - 60) / 10 +
      truth$female * (sex == "Female") +
      truth$bmi * (bmi - 29) +
      clinic_eff[clinic] +
      b0 + b1 * month +
      mu_trt +
      truth$beta_steps_between * (steps_pers - truth$steps_mean) +
      truth$beta_steps_within * (steps - steps_pers) +
      e_serial,
    sbp       = round(sbp_true + rnorm(n(), 0, truth$sd_meas))
  ) |>
  # the steps effect above makes the arm-induced rise in steps part of the
  # total treatment effect; truth$A is the effect *not* acting through steps
  select(id, clinic, arm, age, sex, bmi, visit, month_sched, month,
         steps, sbp)

# ---- missingness -----------------------------------------------------------
# (1) dropout: before each post-baseline visit, the hazard of leaving depends
#     on the most recent observed SBP (MAR) and on arm (observed)
# (2) intermittent missed visits among those still enrolled (MCAR)
miss <- visits |>
  arrange(id, visit) |>
  group_by(id) |>
  group_modify(function(d, key) {
    n_v     <- nrow(d)
    present <- rep(TRUE, n_v)
    dropped <- FALSE
    last_obs <- d$sbp[1]
    for (k in 2:n_v) {
      if (!dropped) {
        lp <- truth$drop_int + truth$drop_sbp * (last_obs - 150) +
              truth$drop_ep * (d$arm[1] == "EP")
        if (runif(1) < plogis(lp)) dropped <- TRUE
      }
      if (dropped) {
        present[k] <- FALSE
      } else if (runif(1) < truth$p_miss_visit) {
        present[k] <- FALSE            # missed this visit, still enrolled
      } else {
        last_obs <- d$sbp[k]
      }
    }
    # last visit attended defines the dropout point
    last_att <- max(which(present))
    tibble(visit = d$visit, present = present,
           dropout = last_att < n_v,
           last_visit = d$visit[last_att])
  }) |>
  ungroup()

complete <- visits |>
  left_join(miss, by = c("id", "visit")) |>
  mutate(arm = factor(arm, levels = arms))

walk_long <- complete |>
  filter(present) |>
  select(id, clinic, arm, age, sex, bmi, visit, month_sched, month, steps, sbp)

dropout_info <- complete |>
  distinct(id, dropout, last_visit)

walk_baseline <- base |>
  select(id, clinic, arm, age, sex, bmi, rand_date) |>
  left_join(walk_long |> filter(visit == 1) |>
              select(id, sbp0 = sbp, steps0 = steps), by = "id") |>
  left_join(dropout_info, by = "id") |>
  mutate(last_month_sched = sched_months[last_visit]) |>
  select(-last_visit)

walk_wide <- walk_long |>
  select(id, clinic, arm, age, sex, bmi, month_sched, sbp) |>
  pivot_wider(names_from = month_sched, values_from = sbp,
              names_prefix = "sbp_m")

walk_complete <- complete |>
  select(id, clinic, arm, age, sex, bmi, visit, month_sched, month, steps,
         sbp, observed = present)

# ---- messy raw export for Module 2 ----------------------------------------
# Visit dates instead of months, text visit labels, a few duplicated rows,
# conflicting duplicates, one implausible value and inconsistent id case.
raw <- walk_long |>
  left_join(base |> select(id, rand_date), by = "id") |>
  mutate(
    visit_label = paste0("M", month_sched),
    visit_date  = rand_date + round(month * 30.44),
    sbp_raw     = as.character(sbp),
    steps_raw   = as.character(steps)
  )

dup_exact <- raw |> slice_sample(n = 6)
dup_conf  <- raw |> slice_sample(n = 3) |>
  mutate(sbp_raw = as.character(as.numeric(sbp_raw) + c(4, -6, 10)),
         visit_date = visit_date + 2)
bad_value <- sample(seq_len(nrow(raw)), 2)
raw$sbp_raw[bad_value[1]] <- paste0(raw$sbp_raw[bad_value[1]], "0")  # x10 typo
raw$sbp_raw[bad_value[2]] <- ""                                       # blank
lower_ids <- sample(unique(raw$id), 4)
raw$id[raw$id %in% lower_ids] <- tolower(raw$id[raw$id %in% lower_ids])

walk_visits_raw <- bind_rows(raw, dup_exact, dup_conf) |>
  slice_sample(prop = 1) |>
  transmute(participant = id, visit_label,
            visit_date = format(visit_date, "%Y-%m-%d"),
            sbp_mmhg = sbp_raw, steps_k = steps_raw)

# ---- write ------------------------------------------------------------------
dir.create(here("data"), showWarnings = FALSE)
write_csv(walk_baseline,   here("data", "walk_baseline.csv"))
write_csv(walk_long,       here("data", "walk_long.csv"))
write_csv(walk_wide,       here("data", "walk_wide.csv"))
write_csv(walk_complete,   here("data", "walk_complete.csv"))
write_csv(walk_visits_raw, here("data", "walk_visits_raw.csv"))
# True mean change from baseline at each scheduled month, by arm (the
# estimands that the course analyses try to recover). Includes the part of
# the programme effect that acts through increased walking.
true_change <- expand.grid(arm = arms, month = sched_months,
                           stringsAsFactors = FALSE) |>
  mutate(change = -truth$A[arm] * (1 - exp(-month / truth$tau)) +
           truth$beta_steps_within * truth$steps_inc[arm] *
           (1 - exp(-month / 1))) |>
  as_tibble()

saveRDS(c(truth, list(clinic_eff = clinic_eff, arm_labels = arm_labels,
                      sched_months = sched_months,
                      true_change = true_change)),
        here("data", "walk_truth.rds"))

message("WALK: ", n, " participants; ", nrow(walk_long), " observed visits of ",
        nrow(walk_complete), " scheduled; dropout ",
        round(100 * mean(dropout_info$dropout)), "%")
