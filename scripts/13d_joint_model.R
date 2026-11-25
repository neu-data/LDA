# =============================================================================
# scripts/13d_joint_model.R
# Module 13 (ADVANCED), section 4: joint model for longitudinal eGFR and time
# to kidney failure (KIDNEY cohort), fitted with JMbayes2.
#
# Purpose : (1) linear mixed model for eGFR, natural spline in time;
#           (2) Cox model for kidney failure;
#           (3) joint model linking them through the current true eGFR;
#           (4) comparison with a naive time-dependent Cox model that uses the
#               last observed eGFR. True association: -0.13 per eGFR unit.
# Inputs  : data/kidney_long.csv, data/kidney_events.csv, kidney_truth.rds
# Outputs : outputs/tables/m13_jm_*.csv, outputs/figures/m13_jm_*.png,
#           outputs/models/m13_jm.rds
# Module  : 13 Advanced extensions
# =============================================================================

source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(survival)                     # advanced module only
  library(JMbayes2)
})
# JMbayes2 masks a few functions; be explicit where it matters
lme   <- nlme::lme
coxph <- survival::coxph

truth <- readRDS(here("data", "kidney_truth.rds"))
kid_long <- load_kidney("long") |>
  mutate(age10 = (age - 60) / 10) |>
  as.data.frame()
kid_ev <- load_kidney("events") |>
  mutate(age10 = (age - 60) / 10) |>
  as.data.frame()

# ---- description -----------------------------------------------------------
kid_desc <- data.frame(
  participants = nrow(kid_ev),
  failures     = sum(kid_ev$failure),
  visits       = nrow(kid_long),
  median_visits = median(table(kid_long$id)),
  median_fu    = median(kid_ev$time)
)
kid_desc
save_table(kid_desc, "m13_kidney_description")

set.seed(13)
show_ids <- c(sample(kid_ev$id[kid_ev$failure == 1], 6),
              sample(kid_ev$id[kid_ev$failure == 0], 6))
p_traj <- kid_long |>
  filter(id %in% show_ids) |>
  left_join(kid_ev |> select(id, failure), by = "id") |>
  mutate(status = ifelse(failure == 1, "Kidney failure", "No failure")) |>
  ggplot(aes(years, egfr, group = id)) +
  geom_line(colour = "#555555", alpha = 0.7) +
  geom_point(size = 0.9, colour = "#555555") +
  facet_wrap(~ status) +
  labs(x = "Years since entry", y = "eGFR (mL/min/1.73m2)")
save_fig(p_traj, "m13_kidney_trajectories")

# ---- submodels --------------------------------------------------------------
# Longitudinal: a natural cubic spline in time (one interior knot at 2.5
# years, fixed boundary knots so predictions use the same basis), faster
# linear decline with diabetes, random intercept and random linear slope.
# A quadratic (years + I(years^2)) gives almost the same fit, but its two
# time terms are strongly correlated and the MCMC mixed poorly for them
# (Rhat about 1.1-1.4); the spline basis is closer to orthogonal.
lme_egfr <- lme(egfr ~ ns(years, knots = 2.5, Boundary.knots = c(0, 5)) +
                  diabetes + diabetes:years + age10,
                random = ~ years | id, data = kid_long,
                control = lmeControl(opt = "optim"))
summary(lme_egfr)

# The same submodel with a quadratic, for comparison with the truth
lme_quad <- lme(egfr ~ years + I(years^2) + diabetes + diabetes:years + age10,
                random = ~ years | id, data = kid_long,
                control = lmeControl(opt = "optim"))
save_table(data.frame(term = names(nlme::fixef(lme_quad)),
                      estimate = nlme::fixef(lme_quad),
                      truth = c(truth$b0, truth$b1, truth$b2, 0,
                                truth$b_age, truth$b1_dm), row.names = NULL),
           "m13_kidney_lme_quadratic_vs_truth")

# Survival: baseline covariate only; the eGFR link is added by jm().
cox_fail <- coxph(Surv(time, failure) ~ diabetes, data = kid_ev)

# ---- joint model -----------------------------------------------------------
# MCMC settings: 3 chains, 25000 iterations each, first 5000 discarded.
# The association parameter mixes slowly: with 10000 iterations (3000
# burn-in) its Rhat was about 1.12, with 25000 (5000 burn-in) about 1.005.
# This takes several minutes. For a publication, also refit with a different
# seed and check that conclusions do not change. The notes load the saved fit.
jm_settings <- list(n_chains = 3L, n_iter = 25000L, n_burnin = 5000L,
                    seed = 2026L)
jm_fit <- jm(cox_fail, lme_egfr, time_var = "years",
             functional_forms = ~ value(egfr),
             n_chains = jm_settings$n_chains, n_iter = jm_settings$n_iter,
             n_burnin = jm_settings$n_burnin, seed = jm_settings$seed)
summary(jm_fit)
save_model(jm_fit, "m13_jm")

s_jm <- summary(jm_fit)
jm_surv <- data.frame(term = rownames(s_jm$Survival), s_jm$Survival,
                      row.names = NULL, check.names = FALSE)
jm_long <- data.frame(term = rownames(s_jm$Outcome1), s_jm$Outcome1,
                      row.names = NULL, check.names = FALSE)
jm_surv; jm_long
save_table(jm_surv, "m13_jm_survival_part")
save_table(jm_long, "m13_jm_longitudinal_part")

# ---- convergence diagnostics -------------------------------------------------
# Gelman-Rubin potential scale reduction factors (Rhat); values close to 1
# (below about 1.05-1.1) suggest the chains agree.
rhat <- data.frame(parameter = c(rownames(s_jm$Survival),
                                 rownames(s_jm$Outcome1)),
                   Rhat = c(s_jm$Survival$Rhat, s_jm$Outcome1$Rhat))
rhat
save_table(rhat, "m13_jm_rhat")

# Trace plot for the association parameter, built from the stored chains
alpha_chains <- jm_fit$mcmc$alphas
trace_df <- bind_rows(lapply(seq_along(alpha_chains), function(k)
  data.frame(iteration = seq_len(nrow(alpha_chains[[k]])),
             value = as.numeric(alpha_chains[[k]][, 1]),
             chain = factor(k))))
p_trace <- ggplot(trace_df, aes(iteration, value, colour = chain)) +
  geom_line(linewidth = 0.25, alpha = 0.8) +
  geom_hline(yintercept = truth$alpha, linetype = 2) +
  scale_colour_manual(values = unname(course_colours[c("blue", "orange",
                                                       "green")])) +
  labs(x = "Iteration after burn-in",
       y = "Association: log HR per eGFR unit", colour = "Chain")
save_fig(p_trace, "m13_jm_trace_alpha", height = 3.5)

# ---- naive comparisons -------------------------------------------------------
# Time-dependent Cox with last-observed eGFR carried forward (counting
# process data built with survival::tmerge)
td <- tmerge(kid_ev |> select(id, diabetes, age10, time, failure),
             kid_ev |> select(id, time, failure),
             id = id, fail = event(time, failure))
td <- tmerge(td, kid_long |> select(id, years, egfr), id = id,
             egfr = tdc(years, egfr))
cox_naive <- coxph(Surv(tstart, tstop, fail) ~ egfr + diabetes, data = td)
summary(cox_naive)

# Baseline eGFR only
base_egfr <- kid_long |> filter(years == 0) |> select(id, egfr0 = egfr)
cox_base <- coxph(Surv(time, failure) ~ egfr0 + diabetes,
                  data = left_join(kid_ev, base_egfr, by = "id"))

ci_cox <- function(m, term, label) {
  s <- summary(m)$coefficients
  data.frame(model = label, estimate = s[term, "coef"],
             lower = s[term, "coef"] - 1.96 * s[term, "se(coef)"],
             upper = s[term, "coef"] + 1.96 * s[term, "se(coef)"])
}
assoc_tab <- bind_rows(
  data.frame(model = "Joint model (current true eGFR)",
             estimate = jm_surv$Mean[jm_surv$term == "value(egfr)"],
             lower = jm_surv$`2.5%`[jm_surv$term == "value(egfr)"],
             upper = jm_surv$`97.5%`[jm_surv$term == "value(egfr)"]),
  ci_cox(cox_naive, "egfr", "Time-dependent Cox (last observed eGFR)"),
  ci_cox(cox_base, "egfr0", "Cox with baseline eGFR only")
) |>
  mutate(truth = truth$alpha,
         hr_per_10_lower = exp(-10 * estimate))
assoc_tab
save_table(assoc_tab, "m13_jm_vs_naive")

# ---- dynamic prediction for one participant ---------------------------------
# Use the first 2 years of data of a participant who later failed, and
# predict the eGFR path and the probability of remaining failure-free.
pid <- kid_ev |> filter(failure == 1, time > 3) |> slice(1) |> pull(id)
nd  <- kid_long |> filter(id == pid, years <= 2)
nd_ev <- kid_ev |> filter(id == pid) |> mutate(time = 2, failure = 0)
nd_all <- left_join(nd, nd_ev |> select(id, time, failure), by = "id")
pred_long <- predict(jm_fit, newdata = nd_all, process = "longitudinal",
                     times = seq(2, 5, by = 0.1), return_newdata = TRUE)
pred_surv <- predict(jm_fit, newdata = nd_all, process = "event",
                     times = seq(2, 5, by = 0.1), return_newdata = TRUE)
dyn <- list(id = pid, true_time = kid_ev$time[kid_ev$id == pid],
            long = pred_long, surv = pred_surv)
save_model(dyn, "m13_jm_dynamic_prediction")

message("13d done")
