# =============================================================================
# 04_simulate_kidney.R
# Simulates the KIDNEY cohort: IRREGULAR observation times, INFORMATIVE visit
# timing and a TIME-TO-EVENT outcome (used in Modules 9, 10 and 13).
#
# THIS IS SIMULATED DATA. The KIDNEY cohort does not exist.
#
# Design (hypothetical observational cohort)
#   300 adults with chronic kidney disease followed for up to 5 years.
#   eGFR (mL/min/1.73m2) is measured at routine clinic visits whose timing is
#   irregular: visits are more frequent when the previous eGFR was low, so the
#   visit process depends on the outcome history.
#   Event: kidney failure (start of dialysis or transplant); measurements stop
#   at the event. Censoring: administrative at 5 years or loss to follow-up.
#
# Generating model
#   True eGFR: m_i(t) = (b0 + u0_i) + (b1 + b1_dm*diabetes + u1_i) t + b2 t^2
#              + b_age * (age - 60)/10
#   Observed eGFR = m_i(t) + error (SD 4).
#   Hazard of kidney failure: h_i(t) = h0 * exp(alpha * (m_i(t) - 30)
#                                               + gamma_dm * diabetes)
#
# Outputs: data/kidney_long.csv, data/kidney_events.csv, data/kidney_truth.rds
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(here)

set.seed(20261029)

n <- 300
truth <- list(
  b0 = 50, b1 = -1.5, b1_dm = -1.5, b2 = -0.15, b_age = -4,
  sd_u0 = 12, sd_u1 = 1.5, cor_u = 0.2, sd_e = 4,
  h0 = 0.05, alpha = -0.13, gamma_dm = 0.3,
  gap_mean_hi = 0.50,  # mean years between visits when last eGFR >= 30
  gap_mean_lo = 0.25,  # ... when last eGFR < 30 (informative timing)
  cens_rate = 0.05,    # loss to follow-up hazard (per year)
  max_fu = 5
)

part <- tibble(
  id       = sprintf("K%03d", seq_len(n)),
  age      = round(pmin(pmax(rnorm(n, 62, 10), 30), 85)),
  sex      = factor(ifelse(runif(n) < 0.45, "Female", "Male"),
                    levels = c("Male", "Female")),
  diabetes = rbinom(n, 1, 0.4)
)
U <- MASS::mvrnorm(n, c(0, 0),
                   matrix(c(truth$sd_u0^2,
                            truth$cor_u * truth$sd_u0 * truth$sd_u1,
                            truth$cor_u * truth$sd_u0 * truth$sd_u1,
                            truth$sd_u1^2), 2))
part$u0 <- U[, 1]
part$u1 <- U[, 2]

m_true <- function(p, t) {
  (truth$b0 + p$u0) + (truth$b1 + truth$b1_dm * p$diabetes + p$u1) * t +
    truth$b2 * t^2 + truth$b_age * (p$age - 60) / 10
}

# event time by inversion on a fine grid (cumulative hazard)
grid <- seq(0.01, truth$max_fu, by = 0.01)
ev <- lapply(seq_len(n), function(i) {
  p   <- part[i, ]
  h   <- truth$h0 * exp(truth$alpha * (m_true(p, grid) - 30) +
                          truth$gamma_dm * p$diabetes)
  H   <- cumsum(h) * 0.01
  e   <- -log(runif(1))
  t_e <- if (any(H >= e)) grid[which(H >= e)[1]] else Inf
  t_c <- min(rexp(1, truth$cens_rate), truth$max_fu)
  tibble(id = p$id, time = round(min(t_e, t_c), 2),
         failure = as.integer(t_e <= t_c))
})
events <- bind_rows(ev)

# irregular visits up to the event/censoring time
obs <- lapply(seq_len(n), function(i) {
  p    <- part[i, ]
  tend <- events$time[i]
  t    <- 0
  out  <- list()
  repeat {
    y <- m_true(p, t) + rnorm(1, 0, truth$sd_e)
    out[[length(out) + 1]] <- tibble(id = p$id, years = round(t, 2),
                                     egfr = round(max(y, 4), 1))
    gap_mean <- if (y < 30) truth$gap_mean_lo else truth$gap_mean_hi
    t <- t + max(rexp(1, 1 / gap_mean), 0.05)
    if (t >= tend) break
  }
  bind_rows(out)
})

kidney_long <- bind_rows(obs) |>
  left_join(part |> select(id, age, sex, diabetes), by = "id") |>
  select(id, age, sex, diabetes, years, egfr)

kidney_events <- part |>
  select(id, age, sex, diabetes) |>
  left_join(events, by = "id")

dir.create(here("data"), showWarnings = FALSE)
write_csv(kidney_long,   here("data", "kidney_long.csv"))
write_csv(kidney_events, here("data", "kidney_events.csv"))
saveRDS(truth, here("data", "kidney_truth.rds"))

message("KIDNEY: ", n, " participants; ", nrow(kidney_long), " visits; ",
        sum(events$failure), " kidney failures")
