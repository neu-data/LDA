# =============================================================================
# 05_simulate_mood.R
# Simulates the MOOD study: the dataset for the FINAL PROJECT.
#
# THIS IS SIMULATED DATA. The MOOD study does not exist.
#
# Design (hypothetical)
#   Two-arm randomised trial of a guided digital self-help programme versus
#   a waiting-list control for adults with depressive symptoms. 4 sites,
#   240 participants (120 per arm).
#   Scheduled assessments: weeks 0 (baseline), 2, 4, 8, 12 and 24.
#   Actual assessment times vary around the schedule.
#   Outcome: depressive symptom score, 0-27 (higher = worse), recorded as an
#   integer. Floor effects appear as participants improve.
#   Derived binary outcome: response = score < 10.
#   Time-varying covariate: average nightly sleep (hours) in the past week.
#   Missingness: intermittent MCAR and dropout that depends on the last
#   observed score AND, in part, on the current unobserved score (so MAR is
#   only an approximation; sensitivity analysis is worthwhile).
#
# The true values are saved in data/mood_truth.rds, which is for INSTRUCTORS.
# Participants should not open it before submitting the project.
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(here)

set.seed(20261030)

n_per_arm <- 120
n         <- 2 * n_per_arm
weeks     <- c(0, 2, 4, 8, 12, 24)

truth <- list(
  intercept   = 16,     # mean baseline score (latent, before rounding)
  age_per10   = -0.4,
  female      = 0.8,
  site_sd     = 0.8,
  # mean change: control  -C * (1 - exp(-week / 6))
  #              programme -(C + D) * (1 - exp(-week / 4)), partial waning
  C           = 3.0,
  D           = 4.0,
  wane        = 0.04,   # loss of programme benefit per week after week 12
  sd_b0       = 3.5,
  sd_b1       = 0.12,   # per week
  cor_b       = 0.1,
  sd_e        = 2.6,
  rho_e       = 0.5,    # AR(1)-type correlation between adjacent visits
  sleep_mean  = 6.4, sleep_sd_pers = 0.8, sleep_sd_occ = 0.6,
  beta_sleep_within  = -0.6,   # per hour, within person
  beta_sleep_between = -1.2,   # per hour, between persons
  p_miss      = 0.05,
  drop_int    = -3.0, drop_last = 0.10, drop_curr = 0.08
)

site_eff <- rnorm(4, 0, truth$site_sd)

part <- tibble(
  id   = sprintf("M%03d", seq_len(n)),
  site = sample(rep(1:4, length.out = n)),
  arm  = factor(sample(rep(c("Waitlist", "Programme"), each = n_per_arm)),
                levels = c("Waitlist", "Programme")),
  age  = round(pmin(pmax(rnorm(n, 38, 12), 18), 70)),
  sex  = factor(ifelse(runif(n) < 0.62, "Female", "Male"),
                levels = c("Male", "Female")),
  sleep_pers = rnorm(n, truth$sleep_mean, truth$sleep_sd_pers)
)
B <- MASS::mvrnorm(n, c(0, 0),
                   matrix(c(truth$sd_b0^2,
                            truth$cor_b * truth$sd_b0 * truth$sd_b1,
                            truth$cor_b * truth$sd_b0 * truth$sd_b1,
                            truth$sd_b1^2), 2))
part$b0 <- B[, 1]; part$b1 <- B[, 2]

jit <- c(0, 0.3, 0.4, 0.6, 0.8, 1.2)

mean_change <- function(arm, w) {
  ctrl <- -truth$C * (1 - exp(-w / 6))
  prog <- -(truth$C + truth$D) * (1 - exp(-w / 4)) +
    truth$wane * pmax(w - 12, 0)
  ifelse(arm == "Programme", prog, ctrl)
}

complete <- part |>
  crossing(visit = seq_along(weeks)) |>
  arrange(id, visit) |>
  group_by(id) |>
  mutate(
    week_sched = weeks[visit],
    week   = round(pmax(week_sched + rnorm(n(), 0, jit[visit]), 0), 1),
    week   = ifelse(visit == 1, 0, week),
    e      = as.numeric(arima.sim(list(ar = truth$rho_e), n = n(),
                                  sd = truth$sd_e * sqrt(1 - truth$rho_e^2))),
    sleep  = round(sleep_pers + 0.3 * (arm == "Programme") *
                     (1 - exp(-week / 4)) +
                     rnorm(n(), 0, truth$sleep_sd_occ), 1)
  ) |>
  ungroup() |>
  mutate(
    latent = truth$intercept + truth$age_per10 * (age - 38) / 10 +
      truth$female * (sex == "Female") + site_eff[site] +
      b0 + b1 * week + mean_change(arm, week) +
      truth$beta_sleep_between * (sleep_pers - truth$sleep_mean) +
      truth$beta_sleep_within * (sleep - sleep_pers) + e,
    score  = pmin(pmax(round(latent), 0), 27)
  )

# dropout: depends on last observed score (MAR part) and current score (MNAR)
complete <- complete |>
  group_by(id) |>
  mutate(observed = {
    obs <- rep(TRUE, n()); gone <- FALSE; last <- score[1]
    for (k in 2:n()) {
      if (!gone && runif(1) < plogis(truth$drop_int +
                                     truth$drop_last * (last - 15) +
                                     truth$drop_curr * (score[k] - 15)))
        gone <- TRUE
      if (gone) obs[k] <- FALSE
      else if (runif(1) < truth$p_miss) obs[k] <- FALSE
      else last <- score[k]
    }
    obs
  }) |>
  ungroup()

mood_long <- complete |>
  filter(observed) |>
  mutate(response = as.integer(score < 10)) |>
  select(id, site, arm, age, sex, visit, week_sched, week, sleep, score,
         response)

# true mean change from baseline at each scheduled week (latent scale)
true_change <- expand.grid(arm = c("Waitlist", "Programme"), week = weeks,
                           stringsAsFactors = FALSE) |>
  mutate(change = mean_change(arm, week) +
           truth$beta_sleep_within * 0.3 * (arm == "Programme") *
           (1 - exp(-week / 4)))

dir.create(here("data"), showWarnings = FALSE)
write_csv(mood_long, here("data", "mood_long.csv"))
saveRDS(c(truth, list(site_eff = site_eff, true_change = true_change)),
        here("data", "mood_truth.rds"))

message("MOOD: ", n, " participants; ", nrow(mood_long), " observed of ",
        nrow(complete), " scheduled assessments")
