# =============================================================================
# scripts/13b_three_level.R
# Module 13 (ADVANCED), section 2: three-level longitudinal data.
#
# Purpose : visits within participants within clinics (WALK, 8 clinics).
#           Fit a three-level LMM, partition the variance, and compare
#           clinic as a random effect with clinic as a fixed effect and with
#           ignoring clinic, for the primary treatment contrast.
# Inputs  : data/walk_long.csv (via load_walk()), data/walk_truth.rds
# Outputs : outputs/tables/m13_threelevel_*.csv, outputs/figures/m13_*.png,
#           outputs/models/m13_threelevel*.rds
# Module  : 13 Advanced extensions
# =============================================================================

source(here::here("R", "setup.R"))
walk  <- load_walk()
truth <- readRDS(here("data", "walk_truth.rds"))

# ---- the clinic level in the data ------------------------------------------
clinic_summary <- walk |>
  filter(visit == 1) |>
  group_by(clinic) |>
  summarise(n_participants = n(), mean_sbp0 = mean(sbp),
            n_UC = sum(arm == "UC"), n_SP = sum(arm == "SP"),
            n_EP = sum(arm == "EP"))
clinic_summary
save_table(clinic_summary, "m13_clinic_summary")

p_clin <- walk |>
  filter(visit == 1) |>
  ggplot(aes(clinic, sbp)) +
  geom_boxplot(outlier.size = 0.6, fill = "#F0F0F0") +
  labs(x = "Clinic", y = "Baseline SBP (mmHg)")
save_fig(p_clin, "m13_clinic_baseline_box")

# ---- variance partition: random intercepts only ------------------------------
# Random intercepts for clinic and for participant. Participant ids are
# unique across clinics, so (1 | id) is automatically nested in clinic.
m3_ri <- lmer(sbp ~ visit_f * arm + (1 | clinic) + (1 | id),
              data = walk, control = lmer_ctrl)
vc <- as.data.frame(VarCorr(m3_ri))
vpart <- data.frame(level = c("Clinic", "Participant (within clinic)",
                              "Visit (residual)"),
                    variance = vc$vcov[match(c("clinic", "id", "Residual"),
                                             vc$grp)]) |>
  mutate(sd = sqrt(variance), share = variance / sum(variance))
vpart
save_table(vpart, "m13_threelevel_variance_partition")
# Correlation between two visits of the same participant, and between two
# different participants at the same clinic
icc_participant <- (vpart$variance[1] + vpart$variance[2]) / sum(vpart$variance)
icc_clinic      <- vpart$variance[1] / sum(vpart$variance)
c(icc_participant = icc_participant, icc_clinic = icc_clinic)

# ---- three-level version of the primary model -------------------------------
m2_primary <- lmer(sbp ~ visit_f * arm + (1 + month | id),
                   data = walk, control = lmer_ctrl)
m3_primary <- lmer(sbp ~ visit_f * arm + (1 + month | id) + (1 | clinic),
                   data = walk, control = lmer_ctrl)
m3_fixed   <- lmer(sbp ~ visit_f * arm + clinic + (1 + month | id),
                   data = walk, control = lmer_ctrl)
VarCorr(m3_primary)
save_model(m3_primary, "m13_threelevel_primary")

# LRT for the clinic variance (both REML with the same fixed effects; the
# null value is on the boundary, so use the 50:50 chi-square mixture)
lrt_clinic <- lr_test(m2_primary, m3_primary)
lrt_clinic$p_boundary <- p_boundary_mixture(lrt_clinic$chisq)
lrt_clinic

# Profile-likelihood interval for the clinic SD: wide with only 8 clinics
ci_prof <- confint(m3_primary, parm = "sd_(Intercept)|clinic", method = "profile",
                   signames = FALSE)
ci_prof
save_table(data.frame(param = rownames(ci_prof), ci_prof, row.names = NULL),
           "m13_threelevel_profile_ci")

# ---- the treatment contrast under three ways of handling clinic -------------
contrast_m12 <- function(m, label) {
  emm <- emmeans(m, ~ visit_f | arm, at = list(visit_f = c("0", "12")))
  chg <- contrast(emm, list(change = c(-1, 1)), by = "arm")
  dif <- contrast(chg, "trt.vs.ctrl", by = NULL, ref = 1)
  s <- summary(dif, infer = TRUE)
  data.frame(model = label, contrast = s$contrast, estimate = s$estimate,
             se = s$SE, df = s$df, lower = s$lower.CL, upper = s$upper.CL)
}
# For m3_fixed, emmeans averages the fixed clinic effects over the clinics
trt_tab <- bind_rows(
  contrast_m12(m2_primary, "Two-level (clinic ignored)"),
  contrast_m12(m3_primary, "Three-level (clinic random)"),
  contrast_m12(m3_fixed,   "Clinic as fixed effect")
)
trt_tab
save_table(trt_tab, "m13_threelevel_treatment_contrasts")

# ---- predicted clinic effects vs the (simulated) truth -----------------------
clin_re <- data.frame(clinic = factor(1:8),
                      predicted = ranef(m3_primary)$clinic[, 1],
                      true = truth$clinic_eff)
clin_re
save_table(clin_re, "m13_clinic_effects_vs_truth")
p_re <- ggplot(clin_re, aes(true, predicted, label = clinic)) +
  geom_abline(linetype = 2, colour = "#9A9A9A") +
  geom_text(size = 3.5) +
  coord_equal() +
  labs(x = "True clinic effect (mmHg)",
       y = "Predicted clinic effect (mmHg)")
save_fig(p_re, "m13_clinic_blup_vs_truth", width = 5, height = 5)

message("13b done")
