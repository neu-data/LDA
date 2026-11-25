# Analysis log

An analysis log is a dated, factual record of what was done to the data and why. Write each entry when the decision is made, not afterwards. A reader (a reviewer, a co-author, or you in a year's time) should be able to reconstruct every result from the log, the scripts and the data.

Rules for using this template:

- Keep one log per study analysis, under version control next to the scripts.
- Never delete or rewrite an entry. If a decision is reversed, add a new entry that refers to the old one.
- Say whether each decision was **pre-specified** (in the protocol or statistical analysis plan, SAP), **data-driven but blind to the treatment comparison**, or **made after seeing results**. The last category must be reported as such.
- Record the effect of each decision on the primary estimate when it is known.

## 1. Study and analysis details

| Item | Entry |
|---|---|
| Study | |
| Analysis | (for example: primary analysis, 12-month outcomes) |
| SAP version and date | |
| Data version (file name, date, checksum or commit) | |
| Analyst(s) | |
| Reviewer | |
| R version and package snapshot (`renv.lock` commit or `sessionInfo()` file) | |
| Code repository and commit at sign-off | |

## 2. Decision log

Copy the block below for each entry. Number entries consecutively.

### Entry NNN

| Field | Entry |
|---|---|
| Date | YYYY-MM-DD |
| Analyst | |
| Script and line / commit | (for example `scripts/04_primary.R`, commit `a1b2c3d`) |
| Category | Data management / Model specification / Convergence or fitting / Diagnostics / Missing data / Sensitivity analysis / Reporting / Other |
| Status | Pre-specified / Data-driven, blind to comparison / Post hoc (after seeing results) |
| What was observed | (the fact that prompted the decision: a warning, a plot, a count) |
| Decision | (what was done) |
| Reason | (why; cite the SAP section if relevant) |
| Alternatives considered | (and why they were not chosen) |
| Effect on primary estimate | (for example "EP − UC at month 12 changed from −9.0 to −8.9 mmHg"; "not yet known"; "none") |
| Reported where | (main paper, supplement, or not reported and why) |
| Reviewed by / date | |

### Example entry (from the course, WALK primary analysis)

| Field | Entry |
|---|---|
| Date | 2026-10-26 |
| Analyst | Course team |
| Script and line / commit | `scripts/11_diagnostics.R`, section 8 |
| Category | Convergence or fitting |
| Status | Data-driven, blind to comparison |
| What was observed | Linear-time model with time in days and age in years gave the warning "Model is nearly unidentifiable: very large eigenvalue. Rescale variables?" |
| Decision | Time expressed in months and age per 10 years (centred at 60). The primary model was also refitted with all available optimisers using `allFit()` |
| Reason | Rescaling is a reparameterisation and does not change the model. All five optimisers gave the same log-likelihood and estimates for the primary model |
| Alternatives considered | Simplifying the random effects (not needed: the fit is not singular) |
| Effect on primary estimate | None. EP slope difference −0.723 mmHg/month before and after rescaling; month-12 EP − UC estimate −8.97 mmHg with every optimiser |
| Reported where | Supplement, statistical methods |
| Reviewed by / date | |

## 3. Running checklists

### Data
- [ ] Number of participants and observations agrees with the participant-flow diagram
- [ ] Duplicate subject-time rows checked and resolved (entry number: )
- [ ] Out-of-range values checked (entry number: )
- [ ] Scheduled versus actual times checked; visit windows applied as in the SAP

### Model fitting
- [ ] Primary model fitted exactly as specified in the SAP
- [ ] Convergence: warnings recorded, gradient checked, optimisers compared if needed
- [ ] Singular fits recorded and the response documented
- [ ] Software versions recorded

### Diagnostics
- [ ] Conditional residuals against fitted values, time and arm
- [ ] Normal Q–Q plot of residuals; scale-location plot
- [ ] Influence of individual participants on the primary estimate
- [ ] Random-effects distribution (with the shrinkage caveat)
- [ ] For GLMMs: simulation-based residuals, dispersion, zero inflation
- [ ] For GEE: residuals by visit, fitted against observed means

### Missing data and sensitivity
- [ ] Missing-data pattern summarised by arm and visit
- [ ] Assumption of the primary analysis (for example MAR) stated
- [ ] Pre-specified sensitivity analyses run (for example delta adjustment)

## 4. Deviations from the SAP

| No. | SAP section | Planned | Done instead | Reason | Log entry |
|---|---|---|---|---|---|
| 1 | | | | | |

## 5. Sign-off

| Role | Name | Date | Comment |
|---|---|---|---|
| Analyst | | | |
| Independent reviewer | | | |
