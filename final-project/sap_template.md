# Statistical analysis plan (SAP) template for a longitudinal study

This template follows the structure recommended by Gamble et al. (2017, *JAMA*) for statistical analysis plans in clinical trials, and it uses the estimand framework of the ICH E9(R1) addendum (2019). It adds the items that longitudinal analyses need: the schedule of measurement occasions, the handling of time, the covariance model, intercurrent events and dropout, and diagnostics.

How to use it:

- Replace every instruction in *italics* with study-specific text. Delete a section only if it genuinely does not apply, and then say why in one line rather than leaving it blank.
- Finalise and date the SAP **before** the treatment comparison is unblinded or, in an observational study, before outcome data are analysed. Later changes go into the version history with reasons.
- Write so that a statistician who has never seen the study could run the analysis from this document alone.

Course references: Module 12 notes (planning, SAP, reporting); Module 11 notes (diagnostics, and the analysis log in `admin/templates/analysis_log_template.md`). A filled example for the WALK study is in `admin/templates/walk_sap_example.md`.

---

## Section 1. Administrative information

1.1 **Title** *(include "statistical analysis plan" and the study name).*

1.2 **Trial or study registration** *(registry and number).*

1.3 **SAP version and date**, with a **version history**: *for each version, the date, the changes made, the reason, and whether the change was made before or after any unblinded or outcome data were seen.*

1.4 **Protocol version** that this SAP refers to.

1.5 **Roles and responsibilities**: *who wrote the SAP, the senior statistician responsible, and who will carry out the analysis (and whether they are blind to allocation).*

1.6 **Signatures**: *author, senior statistician, chief investigator, with dates.*

## Section 2. Introduction

2.1 **Background and rationale** *(a short paragraph; refer to the protocol for detail).*

2.2 **Objectives**: *the primary objective and each secondary objective, each written as a question about change over time or about a difference in change between groups.*

## Section 3. Study methods

3.1 **Design**: *for example, a parallel-group randomised trial with allocation ratio 1:1:1, or a prospective cohort.*

3.2 **Randomisation**: *method, stratification factors and block sizes (stratification factors usually need to appear as covariates).*

3.3 **Sample size**: *the target difference; the outcome SD or the variance components assumed, and their source; the correlation structure assumed; alpha and power; the expected dropout and how it was allowed for. If the calculation was done by simulation, give the simulation model, the number of simulated trials, the Monte Carlo standard error of the power estimate, and where the code is stored.*

3.4 **Framework**: *superiority, non-inferiority (with margin) or equivalence.*

3.5 **Statistical interim analyses and stopping guidance**: *if there are none, say so.*

3.6 **Timing of final analysis**: *for example, when all participants have completed the month-12 visit or withdrawn.*

3.7 **Schedule of measurement occasions**: *a table of scheduled visits, the time of each since randomisation or baseline, the outcomes measured at each, and the visit windows (the range of actual times allowed for each scheduled visit). State how measurements outside the windows will be handled: kept with their actual time, assigned to the nearest scheduled visit, or excluded.*

3.8 **Time scale in the analysis**: *whether the analysis uses scheduled time, actual time, or both (for example scheduled visit as a factor for the mean and actual time for a random slope).*

## Section 4. Statistical principles

4.1 **Confidence intervals and p-values**: *two-sided level; whether p-values are reported for secondary outcomes; adjustment for multiplicity (which comparisons and which method, or why none).*

4.2 **Analysis populations**: *intention-to-treat (all randomised, analysed as randomised), per-protocol, safety. Define each precisely.*

4.3 **Adherence and protocol deviations**: *how adherence is defined and summarised, and which deviations are major.*

## Section 5. Study population

5.1 **Screening and eligibility** data to be summarised.

5.2 **Recruitment and participant flow**: *a CONSORT-style diagram with the numbers randomised and the numbers seen at each scheduled visit, by arm.*

5.3 **Withdrawal and follow-up**: *how withdrawal, loss to follow-up and death are recorded, with the timing of each.*

5.4 **Baseline characteristics**: *list the variables and the summaries (mean and SD, median and IQR, counts and percentages). Do not test baseline differences in randomised trials.*

## Section 6. Analysis

### 6.1 Outcome definitions

*For each outcome: the variable, its units, how and when it is measured, any derivation (for example the mean of two readings), and the outcome type (continuous, binary, ordinal, count with exposure, time to event).*

### 6.2 Estimands (ICH E9(R1))

For the primary outcome, and for each key secondary outcome, give the five attributes.

| Attribute | Entry |
|---|---|
| Treatment condition(s) of interest | *for example "enhanced programme versus usual care, as assigned"* |
| Population | *the participants targeted by the question* |
| Variable (endpoint) | *for example "change in systolic blood pressure from baseline to month 12"* |
| Intercurrent events and the strategy for each | *events after randomisation that affect the existence or interpretation of the outcome (for example stopping the programme, starting a new medication, death). For each, the strategy: treatment policy, hypothetical, composite, while on treatment, or principal stratum.* |
| Population-level summary | *for example "difference in mean change between arms", "odds ratio", "rate ratio"* |

*State explicitly whether the summary is subject-specific (conditional) or population-average (marginal). This matters for binary, ordinal and count outcomes.*

### 6.3 Analysis methods

For the primary analysis give the following.

1. **The model**, written as an equation or as unambiguous software syntax, with:
   - the mean model: how time enters (categorical scheduled visit, linear, piecewise or spline), the group-by-time interaction, and the covariates with their coding and centring;
   - the random effects (which ones, and their covariance) and/or the residual correlation structure; for GEE, the working correlation and the use of robust standard errors;
   - the distribution and link, for outcomes that are not continuous;
   - the estimation method (REML or ML; Laplace or adaptive quadrature for GLMMs) and the software, package and version;
   - the denominator degrees-of-freedom method (for example Satterthwaite or Kenward–Roger).
2. **How the estimand is obtained from the model**: *for example "the contrast of estimated marginal means for change from baseline to month 12, EP minus UC, computed with emmeans".*
3. **Whether baseline is modelled as an outcome or used as a covariate**, and why.
4. **Covariate adjustment**: *list the covariates, including stratification factors, and state that they were chosen before outcome data were seen.*
5. **Assumption checks** to be performed (see Module 11): *residual plots, influence analysis, random-effects checks, simulation-based residuals for GLMMs, the GEE mean-model check*. Give **the pre-specified response** if an assumption is clearly violated, for example: "if the residual variance increases with time, a model with visit-specific residual variances will be fitted and reported as a sensitivity analysis; the primary model will not be changed".
6. **Fitting problems**: *what will be done, in order, if the model fails to converge or is singular. For example: rescale time and covariates; compare optimisers with `allFit()`; then simplify the random-effects structure in a stated order (remove the intercept-slope correlation, then the random slope).*
7. **Alternative methods** if the planned model cannot be fitted.

### 6.4 Missing data

- *The expected pattern (intermittent, dropout) and the anticipated amount.*
- *The assumption of the primary analysis (usually MAR given the variables in the model) and why it is plausible.*
- *How the primary model uses incomplete data (for example a likelihood-based mixed model using all observed visits).*
- *Auxiliary variables or multiple imputation, if any: the imputation model, the number of imputations and the software.*
- *Missing baseline covariates.*

### 6.5 Sensitivity analyses

*Each sensitivity analysis must have a stated purpose, that is, the assumption it relaxes. Examples: delta-adjusted multiple imputation over a range of deltas (MNAR); an alternative covariance structure; analysis on actual rather than scheduled time; the per-protocol population; the influence of individual participants (leave-one-participant-out).*

### 6.6 Pre-specified subgroup analyses

*The subgroups, the interaction test used, and a statement that subgroup analyses are exploratory.*

### 6.7 Secondary outcomes

*For each one, the model, estimand and summary, in the same detail as for the primary outcome or by reference to it.*

### 6.8 Harms

*How adverse events are summarised.*

### 6.9 Statistical software and reproducibility

*R version, packages and versions; how the computational environment is recorded (for example an `renv.lock` file under version control); where the analysis code is archived.*

### 6.10 References

*Methodological references for any non-standard method.*

## Section 7. Reporting

7.1 **Tables and figures to be produced** *(shells or a list): baseline table; participant flow by visit; observed means by arm and visit; model-based estimated means with 95% CIs; the primary contrast table; diagnostic figures (supplement).*

7.2 **What will be reported for every model**: *the model specification, the estimation method, the covariance assumptions, the numbers of participants and observations used, the estimate of the estimand with its 95% CI, the variance components, the assumption checks and their outcome, and all deviations from this plan.*

7.3 **Analysis log**: *decisions made during the analysis are recorded in the analysis log (`admin/templates/analysis_log_template.md`), which is archived with the code.*

## Appendix A. Reporting checklist for a longitudinal analysis

Use this list when writing the statistical methods and results.

**Design and data**

- [ ] Measurement schedule and visit windows; scheduled versus actual time
- [ ] Numbers of participants and observations used, by arm and visit
- [ ] Amount and pattern of missing data by arm and visit; reasons for dropout

**Estimand**

- [ ] Treatment, population, variable, intercurrent events with strategies, population-level summary
- [ ] Conditional (subject-specific) or marginal (population-average) interpretation stated

**Model**

- [ ] Mean model: how time is modelled, interactions, covariates and their coding
- [ ] Random effects and/or residual correlation structure, or GEE working correlation with robust SEs
- [ ] Distribution and link (for outcomes that are not continuous)
- [ ] Estimation method, degrees-of-freedom method, software and versions
- [ ] How the reported contrast is computed from the model

**Assumptions and diagnostics**

- [ ] Which checks were made and what they showed
- [ ] Convergence or singularity problems and how they were resolved
- [ ] Influential participants and their effect on the main estimate

**Missing data**

- [ ] Assumption (MCAR, MAR, MNAR) and its justification
- [ ] Method (likelihood, multiple imputation, weighting) with details
- [ ] Sensitivity analyses to departures from the assumption

**Results**

- [ ] Estimate, 95% CI and (if reported) p-value for each estimand, in the units of the outcome
- [ ] Estimated means or probabilities by arm and visit, with a figure
- [ ] Variance components (or working correlation) with interpretation
- [ ] Results of sensitivity analyses next to the primary result
- [ ] Deviations from the SAP, with reasons

**Reproducibility**

- [ ] Code archived and referenced; package versions recorded (`renv.lock` or `sessionInfo()`)
