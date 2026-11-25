# WALK trial: statistical analysis plan (worked example)

> **Simulated study.** WALK is a fictional trial created for this course. This SAP is a teaching example that shows how the template in `final-project/sap_template.md` is filled in for a longitudinal trial. It contains no results. The design-stage assumptions in Section 3.3 are invented for the example.

## Section 1. Administrative information

1.1 **Title.** Statistical analysis plan for the WALK trial: a three-arm randomised trial of community walking programmes for adults with raised blood pressure.

1.2 **Registration.** Not applicable (fictional study).

1.3 **SAP version and history.**

| Version | Date | Changes | Made before unblinded data were seen? |
|---|---|---|---|
| 0.9 | (draft) | First draft for the trial steering committee | Yes |
| 1.0 | (final) | Visit windows added (3.7); fitting-problem rules added (6.3, item 6) | Yes |

1.4 **Protocol version.** 2.0.

1.5 **Roles.** Written by the trial statistician; approved by the senior statistician. The primary analysis is run by the trial statistician, who remains blind to allocation (arms coded A, B and C) until the database is locked and this SAP is signed.

1.6 **Signatures.** Trial statistician, senior statistician and chief investigator, each with a date.

## Section 2. Introduction

2.1 **Background.** Regular walking lowers blood pressure, but it is uncertain whether a structured community programme achieves sustained reductions, and whether an enhanced version with coaching and step-count feedback adds to a standard programme.

2.2 **Objectives.**

- *Primary:* to estimate the effect of the enhanced programme (EP) and of the standard programme (SP), each compared with usual care (UC), on the change in systolic blood pressure (SBP) from baseline to 12 months.
- *Secondary:* (a) to estimate the effects on change in SBP at months 1, 3, 6 and 9, describing how quickly any effect develops; (b) to compare EP with SP at 12 months; (c) to estimate the effect on daily step count.

## Section 3. Study methods

3.1 **Design.** Parallel-group, three-arm, individually randomised superiority trial with 1:1:1 allocation, in 8 clinics.

3.2 **Randomisation.** Computer-generated, stratified by clinic, permuted blocks of 6 and 9.

3.3 **Sample size.** At the design stage a difference in mean change at 12 months of 5 mmHg between a programme arm and usual care was considered clinically important. Pilot data suggested an SD of SBP of about 13 mmHg at each visit and a correlation of about 0.7 between baseline and 12 months, so the SD of the 12-month change is about $\sqrt{2 \times 13^2 \times 0.3} \approx 10$ mmHg. With two primary comparisons, each tested at two-sided $\alpha = 0.025$ (Bonferroni) with 90% power, a two-sample comparison of change needs $n = 2(z_{0.9875} + z_{0.90})^2 \sigma^2/\Delta^2 = 2(2.241 + 1.282)^2 \times 100/25 \approx 100$ participants per arm. Allowing for 15% loss by month 12 gives $100/0.85 \approx 118$, rounded up to **120 per arm (360 in total)**. This calculation is conservative because the mixed model also uses the intermediate visits and the participants who drop out. A simulation-based check (Module 12, `scripts/12_planning.R`) can be used to confirm it under the planned model.

3.4 **Framework.** Superiority.

3.5 **Interim analyses.** None.

3.6 **Timing of final analysis.** After the last participant's month-12 visit (or withdrawal), database lock, and sign-off of this SAP.

3.7 **Schedule of measurement occasions.**

| Visit | Scheduled month | Window (actual months) | Measurements |
|---|---|---|---|
| 1 (baseline) | 0 | before randomisation | SBP, steps, age, sex, BMI |
| 2 | 1 | 0.5 to 1.9 | SBP, steps |
| 3 | 3 | 2.0 to 4.4 | SBP, steps |
| 4 | 6 | 4.5 to 7.4 | SBP, steps |
| 5 | 9 | 7.5 to 10.4 | SBP, steps |
| 6 | 12 | 10.5 to 14.0 | SBP, steps |

The schedule is front-loaded (visits at 1 and 3 months) because the programmes were expected to act within the first three months; the later visits show whether the effect is sustained. Measurements taken outside a window are kept and assigned to the nearest scheduled visit; their actual time is used in the random-effects part of the model (see 3.8).

3.8 **Time scale.** Scheduled visit (categorical) for the mean; actual months since randomisation for the random slope.

## Section 4. Statistical principles

4.1 **Confidence intervals and p-values.** Two-sided 95% CIs for all estimates. The two primary comparisons (EP vs UC, SP vs UC) are tested with the Holm procedure at familywise $\alpha = 0.05$. Secondary outcomes are reported with 95% CIs and without p-values.

4.2 **Analysis populations.** Intention-to-treat: all randomised participants, analysed in the arm to which they were allocated, using all available SBP measurements. Per-protocol: participants who attended at least 75% of programme sessions (programme arms) or who did not join an external walking group (UC).

4.3 **Adherence.** Session attendance summarised by arm as median (IQR) and the percentage attending at least 75%.

## Section 5. Study population

5.1 **Screening.** Numbers screened, eligible and randomised.

5.2 **Participant flow.** A CONSORT diagram, plus a table of the number of participants seen at each scheduled visit by arm.

5.3 **Withdrawal.** Date and reason for withdrawal or loss to follow-up; the last visit attended.

5.4 **Baseline characteristics.** Age, sex, BMI, clinic, SBP and steps by arm: mean (SD) or n (%). No significance tests.

## Section 6. Analysis

### 6.1 Outcome definitions

- **SBP** (mmHg): mean of the second and third of three seated readings taken with a validated automatic device. Continuous.
- **Steps** (thousands/day): mean daily step count over the 7 days before the visit, from a sealed pedometer. Continuous.

### 6.2 Primary estimand

| Attribute | Entry |
|---|---|
| Treatment conditions | Assignment to EP (or SP) versus assignment to UC |
| Population | Adults with raised SBP meeting the trial eligibility criteria |
| Variable | Change in SBP from baseline to month 12 |
| Intercurrent events and strategies | (1) Stopping programme attendance: treatment policy (the outcome is used whatever the attendance). (2) Starting or changing antihypertensive medication: treatment policy for the primary estimand; a hypothetical strategy (as if no medication change) is a supplementary estimand. (3) Death: expected to be rare; handled by the treatment policy strategy, with deaths listed. |
| Population-level summary | Difference in mean change (programme minus UC), in mmHg. For a continuous outcome and an identity link the conditional and marginal differences coincide. |

### 6.3 Analysis methods

1. **Primary model.** A linear mixed-effects model fitted by REML:

   $$ y_{ij} = \beta_0 + \sum_{v} \gamma_v \,\mathbb{1}(\text{visit}_{ij} = v) + \beta_{\text{SP}}\text{SP}_i + \beta_{\text{EP}}\text{EP}_i + \sum_{v}\left(\delta_{v}^{\text{SP}}\text{SP}_i + \delta_{v}^{\text{EP}}\text{EP}_i\right)\mathbb{1}(\text{visit}_{ij} = v) + b_{0i} + b_{1i} t_{ij} + \varepsilon_{ij} $$

   with $v \in \{1, 3, 6, 9, 12\}$, $(b_{0i}, b_{1i}) \sim N(\mathbf{0}, G)$ with unstructured $G$, and $\varepsilon_{ij} \sim N(0, \sigma^2)$ independent. In R:
   `lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk, control = lmerControl(optimizer = "bobyqa"))` (packages lme4 and lmerTest; Satterthwaite degrees of freedom).
2. **Estimand from the model.** For each arm, the estimated marginal mean at month 12 minus that at month 0, computed with emmeans; the primary contrasts are the differences of these changes, EP minus UC and SP minus UC. In this model they equal the coefficients $\delta_{12}^{\text{EP}}$ and $\delta_{12}^{\text{SP}}$.
3. **Baseline.** Baseline SBP is modelled as an outcome (visit 1), so participants with any SBP measurement contribute. Because of randomisation the arms have the same expected baseline mean.
4. **Covariates.** The primary model is unadjusted. A supplementary adjusted analysis adds age (per 10 years, centred at 60), sex and BMI (centred at 29 kg/m²), chosen before the data were seen. Clinic (a stratification factor) is added as a random intercept in a sensitivity analysis.
5. **Assumption checks** (recorded in the analysis log, reported in the supplement): standardised conditional residuals against fitted values, actual time and arm; normal Q–Q and scale-location plots; Q–Q plots of predicted random effects (interpreted with the shrinkage caveat); leave-one-participant-out influence on the primary contrasts. **Pre-specified responses:** if the residual variance clearly differs between visits, a model with visit-specific residual variances is reported as a sensitivity analysis. If a single participant changes a primary contrast by more than half of its standard error, the data for that participant are checked against the source documents and the analysis is repeated without them as a sensitivity analysis. The primary model is not changed on the basis of diagnostics.
6. **Fitting problems.** If a convergence warning occurs: (a) check the gradient and refit with all available optimisers (`allFit()`), and if the estimates agree to the reported precision, keep the fit; (b) otherwise remove the intercept-slope correlation; (c) then remove the random slope, using a random intercept with an exponential serial correlation fitted with `nlme::lme` instead.
7. **Alternative method.** If the mixed model cannot be fitted, an unstructured-covariance model for the scheduled visits (`nlme::gls` with `corSymm` and `varIdent`) will be used.

### 6.4 Missing data

About 5% of visits are expected to be missed intermittently, and 15% of participants are expected to drop out by month 12. The primary analysis uses all observed SBP measurements under the assumption that the data are missing at random (MAR) given arm and the observed SBP history, which the likelihood-based mixed model accommodates without imputation. Dropout related to high SBP is plausible and is consistent with MAR when it depends on previously observed values. No baseline covariates are expected to be missing in the primary model.

### 6.5 Sensitivity analyses

| Analysis | Assumption relaxed |
|---|---|
| Delta-adjusted multiple imputation: after dropout, imputed SBP values in the programme arms shifted by +2, +4 and +6 mmHg (50 imputations) | MAR (allows worse outcomes after dropout in the programme arms) |
| Exponential serial correlation added to the residuals (`nlme::lme` with `corExp`) | Independent residual errors |
| Clinic random intercept | No clustering by clinic |
| Per-protocol population | Treatment-policy strategy for non-attendance (gives a different estimand; interpreted with caution) |
| Leave-one-participant-out | No single participant drives the result |

### 6.6 Subgroups

Sex, and age below or above 60 years, through three-way interactions (subgroup × arm × visit) with the month-12 contrast estimated within each subgroup. These are exploratory.

### 6.7 Secondary outcomes

Change in SBP at months 1, 3, 6 and 9, and EP minus SP at month 12: contrasts from the primary model. Steps: the same model structure as the primary analysis.

### 6.8 Harms

Falls and musculoskeletal injuries by arm: counts and percentages of participants.

### 6.9 Software and reproducibility

R 4.6.1 with lme4, lmerTest, nlme, emmeans and mice. The package versions are recorded in an `renv.lock` file committed to the study repository with the analysis code; `sessionInfo()` output is saved with the results.

## Section 7. Reporting

7.1 **Tables and figures.** Baseline table; participant flow by visit; observed mean SBP by arm and visit; estimated mean change by arm and visit with 95% CIs (figure); table of primary and secondary contrasts; sensitivity analyses next to the primary result (forest-style figure); diagnostic figures in the supplement.

7.2 **Content of the statistical methods section.** The model, estimation method and software; how the contrasts are computed; the missing-data assumption and the sensitivity analyses; the diagnostics and any deviations from this SAP.

7.3 **Analysis log.** All decisions made during the analysis are recorded in `analysis_log.md` (based on `admin/templates/analysis_log_template.md`) and archived with the code.
