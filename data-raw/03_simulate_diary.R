# =============================================================================
# 03_simulate_diary.R
# Simulates the DIARY study: repeated COUNT outcomes with varying exposure.
#
# THIS IS SIMULATED DATA. The DIARY study does not exist.
#
# Design (hypothetical)
#   Two-arm randomised trial of a preventive treatment for migraine.
#   160 participants (80 per arm). A 28-day baseline diary period (period 0)
#   is followed by six 28-day treatment periods (1-6).
#   Outcome: number of migraine days recorded in the period.
#   Exposure: number of days on which the diary was completed (14-28), so the
#   natural quantity is a RATE (migraine days per diary day).
#
# Generating model (subject-specific negative binomial with random intercept)
#   log E[Y_ij | b_i] = log(days_ij) + beta0 + b_i
#                       + beta_post * post_ij
#                       + beta_trt  * active_i * post_ij
#                       + beta_age  * (age - 40)/10
#   Y_ij | b_i ~ NegBin(mean, size = theta), truncated at days_ij.
#   No structural (excess) zeros are generated: zeros arise from low rates.
#
# Outputs: data/diary_long.csv, data/diary_truth.rds
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(here)

set.seed(12)

n_per_arm <- 80
n         <- 2 * n_per_arm
periods   <- 0:6

truth <- list(
  beta0     = log(0.18),  # about 5 migraine days per 28 days at baseline
  beta_post = -0.10,      # placebo / period change after baseline
  beta_trt  = -0.40,      # active vs placebo after baseline (rate ratio 0.67)
  beta_age  = -0.05,      # per 10 years
  sd_b      = 0.55,       # random-intercept SD (log scale)
  theta     = 5,          # negative binomial size (overdispersion)
  p_short   = 0.30,       # chance a diary period is incomplete
  drop_p    = 0.035       # per-period dropout probability (MCAR)
)

part <- tibble(
  id  = sprintf("D%03d", seq_len(n)),
  arm = factor(sample(rep(c("Placebo", "Active"), each = n_per_arm)),
               levels = c("Placebo", "Active")),
  age = round(pmin(pmax(rnorm(n, 40, 11), 18), 70)),
  sex = factor(ifelse(runif(n) < 0.75, "Female", "Male"),
               levels = c("Male", "Female")),
  b   = rnorm(n, 0, truth$sd_b),
  drop_after = {
    d <- rgeom(n, truth$drop_p) + 1     # last completed period
    pmin(d, max(periods))
  }
)

diary_long <- part |>
  crossing(period = periods) |>
  filter(period <= drop_after) |>
  mutate(
    days = ifelse(runif(n()) < truth$p_short,
                  sample(14:27, n(), replace = TRUE), 28L),
    post = as.integer(period > 0),
    mu   = days * exp(truth$beta0 + b + truth$beta_post * post +
                        truth$beta_trt * (arm == "Active") * post +
                        truth$beta_age * (age - 40) / 10),
    migraine_days = pmin(rnbinom(n(), mu = mu, size = truth$theta), days)
  ) |>
  select(id, arm, age, sex, period, post, days, migraine_days)

dir.create(here("data"), showWarnings = FALSE)
write_csv(diary_long, here("data", "diary_long.csv"))
saveRDS(c(truth, list(periods = periods)), here("data", "diary_truth.rds"))

message("DIARY: ", n, " participants; ", nrow(diary_long), " diary periods")
