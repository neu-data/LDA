# =============================================================================
# 10b_missing_data_simulations.R
# Module 10. Missing data and dropout: simulation studies.
#
# Purpose
#   (1) Impose MCAR, MAR and MNAR dropout on the complete WALK data
#       (walk_complete.csv, available only because the data are simulated)
#       and compare four analyses: observed means, complete-case mixed
#       model, last observation carried forward (LOCF) and the mixed model
#       using all observed data.
#   (2) Show the consequences of informative visit timing, using the
#       generating model of the KIDNEY cohort (kidney_truth.rds).
#
# Inputs   data/walk_complete.csv, data/walk_baseline.csv,
#          data/kidney_truth.rds
# Outputs  outputs/tables/m10_sim_mechanisms_raw.csv
#          outputs/tables/m10_sim_mechanisms.csv
#          outputs/tables/m10_sim_visit_timing_raw.csv
#          outputs/tables/m10_sim_visit_timing.csv
#          outputs/figures/m10_sim_mechanisms.png
#          outputs/figures/m10_sim_visit_timing.png
# Run time about 10-15 minutes. The notes, slides and answer key read the
# saved tables; they run this script only if the tables are missing.
# =============================================================================

source(here::here("R", "setup.R"))

# Faster lmer settings for repeated fits: the derivative check is skipped
# because convergence of these well-behaved models was checked on the real
# data; it changes nothing about the estimates.
sim_ctrl <- lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)

# ---- part 1: missing-data mechanisms on the complete WALK data ------------

wc <- load_walk("complete")

# Difference in mean change (baseline to month 12) vs usual care from the
# primary model. Computed from the fixed effects directly (faster than
# emmeans and identical for this parameterisation).
primary_contrast <- function(d) {
  m <- suppressMessages(suppressWarnings(
    lmer(sbp ~ visit_f * arm + (1 + month | id), data = d,
         control = sim_ctrl)))
  b <- fixef(m)
  c(SP = unname(b["visit_f12:armSP"]), EP = unname(b["visit_f12:armEP"]))
}

# The benchmark: the same model on the complete data (no missingness).
bench <- primary_contrast(wc)

# Monotone dropout. Before each post-baseline visit k a participant still in
# the study leaves with probability plogis(a + b * x). The story: people whose
# blood pressure has risen since baseline become discouraged and leave.
#   MCAR: x = 0 (leaving is unrelated to blood pressure);
#   MAR : x = change from baseline at the last *observed* visit (k - 1),
#         which the analyst can see;
#   MNAR: x = change from baseline at visit k itself, the value that goes
#         unrecorded because the participant has left.
impose_dropout <- function(d, mechanism, a, b) {
  d <- arrange(d, id, visit)
  sbp <- matrix(d$sbp, ncol = 6, byrow = TRUE)
  keep <- matrix(TRUE, nrow(sbp), 6)
  for (k in 2:6) {
    at_risk <- keep[, k - 1]
    x <- switch(mechanism,
                MCAR = 0,
                MAR  = sbp[, k - 1] - sbp[, 1],
                MNAR = sbp[, k] - sbp[, 1])
    leave <- at_risk & (runif(nrow(sbp)) < plogis(a + b * x))
    keep[, k] <- at_risk & !leave
  }
  d$observed_sim <- as.vector(t(keep))
  d
}

# The same intercept and slope give roughly a third of participants
# dropping out by month 12 under each mechanism (see the dropout column).
mech_par <- list(MCAR = c(a = -2.3, b = 0),
                 MAR  = c(a = -2.3, b = 0.15),
                 MNAR = c(a = -2.3, b = 0.15))

analyse_once <- function(d) {
  obs <- filter(d, observed_sim)
  # (a) naive: observed mean at month 12 minus observed mean at baseline
  om <- obs |> group_by(arm, month_sched) |> summarise(m = mean(sbp)) |>
    tidyr::pivot_wider(names_from = month_sched, values_from = m)
  ch <- om$`12` - om$`0`
  naive <- c(SP = ch[2] - ch[1], EP = ch[3] - ch[1])
  # (b) complete cases: participants observed at month 12
  cc_ids <- obs$id[obs$month_sched == 12]
  cc <- primary_contrast(filter(obs, id %in% cc_ids))
  # (c) LOCF: last observed SBP used as the month-12 value
  locf_d <- obs |> group_by(id, arm) |>
    summarise(base = sbp[visit == 1], last = sbp[which.max(visit)]) |>
    ungroup()
  lf <- coef(lm(I(last - base) ~ arm, data = locf_d))
  locf <- c(SP = unname(lf["armSP"]), EP = unname(lf["armEP"]))
  # (d) mixed model on all observed data (valid under MAR)
  lmm <- primary_contrast(obs)
  dropout <- mean(!tapply(d$observed_sim, d$id, function(z) z[6]))
  rbind(
    data.frame(method = "Observed means", arm = names(naive), est = naive),
    data.frame(method = "Complete cases", arm = names(cc), est = cc),
    data.frame(method = "LOCF", arm = names(locf), est = locf),
    data.frame(method = "Mixed model, all data", arm = names(lmm), est = lmm)
  ) |> mutate(dropout = dropout)
}

set.seed(20261101)
n_rep <- 200
sim1 <- lapply(names(mech_par), function(mech) {
  p <- mech_par[[mech]]
  bind_rows(lapply(seq_len(n_rep), function(r) {
    d <- impose_dropout(wc, mech, p["a"], p["b"])
    analyse_once(d) |> mutate(rep = r)
  })) |> mutate(mechanism = mech)
}) |> bind_rows()
rownames(sim1) <- NULL
save_table(sim1, "m10_sim_mechanisms_raw")

# Bias is measured against the complete-data estimate (the answer we would
# have obtained with no missing data), so it isolates the effect of the
# missingness from ordinary sampling error. Monte Carlo SEs follow
# Morris et al. (2019).
sim1_sum <- sim1 |>
  mutate(benchmark = bench[arm]) |>
  group_by(mechanism, method, arm) |>
  summarise(n_rep = n(), mean_est = mean(est),
            bias = mean(est - benchmark),
            mcse_bias = sd(est) / sqrt(n()),
            emp_se = sd(est), dropout = mean(dropout)) |>
  ungroup() |>
  mutate(benchmark = bench[arm],
         mechanism = factor(mechanism, levels = c("MCAR", "MAR", "MNAR")),
         method = factor(method, levels = c("Observed means", "LOCF",
                                            "Complete cases",
                                            "Mixed model, all data"))) |>
  arrange(mechanism, method, arm)
save_table(sim1_sum, "m10_sim_mechanisms")

p1 <- ggplot(sim1_sum,
             aes(bias, method, colour = arm)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = bias - 1.96 * mcse_bias,
                      xmax = bias + 1.96 * mcse_bias),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ mechanism) +
  scale_colour_arm(name = "Contrast vs usual care",
                   labels = c(SP = "SP vs UC", EP = "EP vs UC")) +
  labs(x = "Bias in difference in mean change at month 12 (mmHg)",
       y = NULL)
save_fig(p1, "m10_sim_mechanisms", width = 9, height = 4)

# ---- part 2: informative visit timing (KIDNEY generating model) ----------

kt <- readRDS(here("data", "kidney_truth.rds"))

# Simulate one cohort. informative = TRUE uses the KIDNEY visit rule (mean
# gap 0.25 years after an eGFR < 30, otherwise 0.5 years); FALSE uses a
# fixed mean gap of 0.4 years regardless of eGFR. failure = TRUE stops
# measurement at kidney failure, as in the real KIDNEY data.
sim_kidney <- function(n = 300, informative = TRUE, failure = FALSE) {
  age10 <- (round(pmin(pmax(rnorm(n, 62, 10), 30), 85)) - 60) / 10
  dm <- rbinom(n, 1, 0.4)
  S <- matrix(c(kt$sd_u0^2, kt$cor_u * kt$sd_u0 * kt$sd_u1,
                kt$cor_u * kt$sd_u0 * kt$sd_u1, kt$sd_u1^2), 2)
  U <- MASS::mvrnorm(n, c(0, 0), S)
  m_true <- function(i, t) {
    kt$b0 + U[i, 1] + (kt$b1 + kt$b1_dm * dm[i] + U[i, 2]) * t +
      kt$b2 * t^2 + kt$b_age * age10[i]
  }
  grid <- seq(0.01, kt$max_fu, by = 0.01)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    t_end <- kt$max_fu
    if (failure) {
      H <- cumsum(kt$h0 * exp(kt$alpha * (m_true(i, grid) - 30) +
                                kt$gamma_dm * dm[i])) * 0.01
      e <- -log(runif(1))
      if (any(H >= e)) t_end <- grid[which(H >= e)[1]]
    }
    t <- 0; ts <- numeric(0); ys <- numeric(0)
    repeat {
      y <- m_true(i, t) + rnorm(1, 0, kt$sd_e)
      ts <- c(ts, t); ys <- c(ys, y)
      gap <- if (!informative) 0.4 else if (y < 30) kt$gap_mean_lo else
        kt$gap_mean_hi
      t <- t + max(rexp(1, 1 / gap), 0.05)
      if (t >= t_end) break
    }
    out[[i]] <- data.frame(id = i, years = ts, egfr = pmax(ys, 4),
                           diabetes = dm[i], age10 = age10[i])
  }
  d <- do.call(rbind, out)
  # Target: mean eGFR trajectory in the whole cohort at years 1-4, i.e. the
  # population mean of the true trajectories (had everyone been followed).
  truth <- sapply(1:4, function(t)
    mean(kt$b0 + (kt$b1 + kt$b1_dm * dm) * t + kt$b2 * t^2 +
           kt$b_age * age10))
  list(d = d, truth = truth, pop = data.frame(diabetes = dm, age10 = age10))
}

analyse_kidney <- function(s) {
  d <- s$d
  # naive: mean of all eGFR values recorded within 3 months of each year
  naive <- sapply(1:4, function(t) mean(d$egfr[abs(d$years - t) <= 0.25]))
  m <- suppressMessages(suppressWarnings(
    lmer(egfr ~ years + I(years^2) + diabetes + diabetes:years + age10 +
           (1 + years | id), data = d, control = sim_ctrl)))
  lmm <- sapply(1:4, function(t)
    mean(predict(m, newdata = transform(s$pop, years = t), re.form = NA)))
  data.frame(year = rep(1:4, 2),
             method = rep(c("Naive window mean", "Mixed model"), each = 4),
             error = c(naive, lmm) - rep(s$truth, 2))
}

set.seed(20261102)
n_rep_k <- 100
scen <- expand.grid(informative = c(FALSE, TRUE), failure = c(FALSE, TRUE))
sim2 <- lapply(seq_len(nrow(scen)), function(j) {
  bind_rows(lapply(seq_len(n_rep_k), function(r) {
    analyse_kidney(sim_kidney(informative = scen$informative[j],
                              failure = scen$failure[j])) |>
      mutate(rep = r)
  })) |>
    mutate(visits = ifelse(scen$informative[j],
                           "Visits depend on last eGFR",
                           "Visits unrelated to eGFR"),
           follow_up = ifelse(scen$failure[j],
                              "Follow-up stops at kidney failure",
                              "No kidney failure"))
}) |> bind_rows()
save_table(sim2, "m10_sim_visit_timing_raw")

sim2_sum <- sim2 |>
  group_by(visits, follow_up, method, year) |>
  summarise(n_rep = n(), bias = mean(error),
            mcse = sd(error) / sqrt(n())) |>
  ungroup()
save_table(sim2_sum, "m10_sim_visit_timing")

p2 <- ggplot(sim2_sum, aes(year, bias, colour = method)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = bias - 1.96 * mcse, ymax = bias + 1.96 * mcse),
                  position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  facet_grid(follow_up ~ visits) +
  scale_colour_manual(values = c("Naive window mean" = "#D55E00",
                                 "Mixed model" = "#0072B2"), name = NULL) +
  labs(x = "Years since entry", y = "Bias in mean eGFR (mL/min/1.73m²)")
save_fig(p2, "m10_sim_visit_timing", width = 8, height = 5.5)

message("Benchmark (complete data): SP ", round(bench["SP"], 2),
        ", EP ", round(bench["EP"], 2))
print(sim1_sum, n = Inf)
print(sim2_sum, n = Inf)
