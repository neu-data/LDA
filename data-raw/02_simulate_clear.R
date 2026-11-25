# =============================================================================
# 02_simulate_clear.R
# Simulates the CLEAR study: repeated BINARY and ORDINAL outcomes.
#
# THIS IS SIMULATED DATA. The CLEAR study does not exist.
#
# Design (hypothetical)
#   Two-arm randomised trial of a topical treatment for a chronic itchy skin
#   condition. 200 participants (100 per arm).
#   Scheduled visits: weeks 0, 2, 4, 6, 8, 12.
#   Ordinal outcome: itch severity (0 = none, 1 = mild, 2 = moderate,
#   3 = severe). Binary outcome: moderate-or-severe itch (severity >= 2).
#
# Generating model (subject-specific, cumulative logit with random intercept)
#   latent Y*_ij = b_i + beta_week * week + beta_trtweek * active * week
#                  + beta_age * (age - 45)/10 + logistic error
#   severity = k when Y* falls between cut-points c_k and c_(k+1).
#   It follows that the binary indicator (severity >= 2) obeys a logistic
#   random-intercept model with intercept -c_2 and the same slopes.
#
# Outputs: data/clear_long.csv, data/clear_truth.rds
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(here)

set.seed(20261027)

n_per_arm <- 100
weeks     <- c(0, 2, 4, 6, 8, 12)
n         <- 2 * n_per_arm

truth <- list(
  cuts         = c(-2.5, 0, 2.5),  # cut-points between 0|1, 1|2, 2|3
  beta_week    = -0.06,            # log-odds per week, placebo (conditional)
  beta_trtweek = -0.16,            # additional log-odds per week, active
  beta_age     = 0.20,             # per 10 years
  beta_female  = 0.30,
  sd_b         = 2.0,              # random-intercept SD (logit scale)
  p_miss_visit = 0.04,
  drop_int     = -3.9,             # dropout logit per interval
  drop_sev     = 0.5               # per unit of last observed severity
)

part <- tibble(
  id  = sprintf("C%03d", seq_len(n)),
  arm = factor(sample(rep(c("Placebo", "Active"), each = n_per_arm)),
               levels = c("Placebo", "Active")),
  age = round(pmin(pmax(rnorm(n, 45, 13), 18), 80)),
  sex = factor(ifelse(runif(n) < 0.55, "Female", "Male"),
               levels = c("Male", "Female")),
  b   = rnorm(n, 0, truth$sd_b)
)

complete <- part |>
  crossing(week = weeks) |>
  mutate(
    lp = b + truth$beta_week * week +
      truth$beta_trtweek * (arm == "Active") * week +
      truth$beta_age * (age - 45) / 10 +
      truth$beta_female * (sex == "Female"),
    ystar    = lp + rlogis(n()),
    severity = findInterval(ystar, truth$cuts),
    itch     = as.integer(severity >= 2)
  ) |>
  arrange(id, week)

# dropout depends on the last observed severity (MAR); intermittent MCAR
complete <- complete |>
  group_by(id) |>
  mutate(present = {
    pres <- rep(TRUE, n()); dropped <- FALSE; last <- severity[1]
    for (k in 2:n()) {
      if (!dropped && runif(1) < plogis(truth$drop_int +
                                         truth$drop_sev * last)) dropped <- TRUE
      if (dropped) pres[k] <- FALSE
      else if (runif(1) < truth$p_miss_visit) pres[k] <- FALSE
      else last <- severity[k]
    }
    pres
  }) |>
  ungroup()

clear_long <- complete |>
  filter(present) |>
  mutate(visit = match(week, weeks)) |>
  select(id, arm, age, sex, visit, week, severity, itch)

dir.create(here("data"), showWarnings = FALSE)
write_csv(clear_long, here("data", "clear_long.csv"))
saveRDS(c(truth, list(weeks = weeks)), here("data", "clear_truth.rds"))

message("CLEAR: ", n, " participants; ", nrow(clear_long),
        " observed visits of ", nrow(complete))
