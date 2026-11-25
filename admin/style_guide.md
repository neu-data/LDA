# Authoring and style guide

This guide keeps every module of the course consistent. Read it before you write or edit any material. It covers file names, structure, terminology, notation, writing style, code conventions and verification.

## 1. Project layout and file names

| Folder | Content | File name pattern (NN = two-digit module number) |
|---|---|---|
| `R/` | Shared helper functions, loaded by `source(here::here("R", "setup.R"))` | Do **not** edit existing files. If a module needs its own helpers, define them inside the module script/notes. |
| `data/`, `data-raw/` | Simulated datasets and their simulation scripts | Do not modify. |
| `scripts/` | Runnable analysis scripts that save to `outputs/` | `NN_topic.R`; several per module allowed: `08a_glmm_binary.R`, `08b_glmm_count.R` |
| `notes/` | Course notes (Quarto, HTML) | `NN-topic.qmd` |
| `slides/` | Reveal.js slides with instructor notes | `NN-topic-slides.qmd` |
| `workbook/` | Participant exercises | `wbNN-topic.qmd` |
| `answer-key/` | Instructor solutions with executed output | `akNN-topic.qmd` |
| `assessments/` | Knowledge checks; practical assessments (Modules 2, 4, 6, 8, 10) | `kcNN-topic.qmd`, `paNN-topic.qmd` |
| `outputs/figures`, `outputs/tables`, `outputs/models` | Saved by scripts via `save_fig()`, `save_table()`, `save_model()` | Prefix with module: `m04_spaghetti.png` |
| `admin/cheatsheet_input/`, `admin/glossary_input/` | Inputs for the cheat sheets and glossary | `mNN.md` |

Module topics (use these short names): `00-refresher`, `01-foundations`, `02-data-preparation`, `03-traditional-methods`, `04-linear-mixed-models`, `05-covariance-models`, `06-interpretation-prediction`, `07-gee`, `08-glmm`, `09-time-varying-nonlinear`, `10-missing-data`, `11-diagnostics`, `12-planning-reporting`, `13-advanced`.

Render from the project root with Quarto:

```
"/c/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe" render notes/04-linear-mixed-models.qmd
```

Run scripts with `"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" scripts/04_linear_mixed_models.R` from the project root. Folder-level settings live in each folder's `_metadata.yml`: formats, CSS, MathML math, embedded resources, `comment: "#>"`. Do not override them without reason.

## 2. Front matter and shared code

Every `.qmd` starts like this:

```yaml
---
title: "Module 4. Linear mixed-effects models"
subtitle: "Course notes"          # or "Slides", "Workbook", "Answer key", "Knowledge check"
bibliography: ../references/references.bib
---
```

Its first code chunk (label `setup`, `#| include: true` in notes, `#| echo: false` in slides) is:

```r
source(here::here("R", "setup.R"))
walk <- load_walk()      # or load_clear(), load_diary(), load_kidney()
```

`setup.R` loads the packages: here, readr, dplyr, tidyr, ggplot2, patchwork, lme4, lmerTest, nlme, glmmTMB, geepack, emmeans, mgcv, mice, DHARMa. It sets the ggplot theme and defines `lmer_ctrl` (use `control = lmer_ctrl` in every `lmer()` call) and these helpers:

- **Loading data:** `load_walk("long"|"baseline"|"wide"|"complete")`, `load_clear()`, `load_diary()`, `load_kidney("long"|"events")`
- **Checking long data:** `check_long(data, id, time, "outcome")`, `visit_correlation(data, id, visit, outcome)`, `visit_flow(data, id, visit, group)`
- **Formatting and extracting results:** `fmt_ci(est, lo, hi, digits)`, `fmt_p(p)`, `fixef_table(model, exponentiate = FALSE)`, `varcomp_table(m)`, `icc_lmer(m)`, `implied_cov_lmer(m, times)`, `cov_to_cor(V)`, `lr_test(m0, m1)`, `p_boundary_mixture(chisq)`
- **Plotting:** `plot_spaghetti(data, x, y, id, group, n_show)`, `plot_observed_means(data, x, y, group)`, `plot_missing_map(grid, id, visit, observed, order_by)`, `scale_colour_arm()`, `scale_fill_arm()`, `arm_colours`, `walk_arm_labels`
- **Saving outputs:** `save_fig(p, name, width, height)`, `save_table(x, name)`, `save_model(m, name)`
- **Simulation:** `sim_two_arm(...)`, `power_two_arm(...)`

Only the advanced module may load `ordinal`, `survival` and `JMbayes2`. Do not add other packages. If you think one is essential, explain why in your report rather than adding it. `MASS::` and `rlang::` calls are acceptable because they ship with R or are already dependencies.

## 3. Datasets and the true values

Read `data/dictionaries/data_dictionary.md`. **Every dataset is simulated.** Say so at the start of each document with this callout:

```markdown
::: {.callout-warning title="Simulated data"}
All examples use simulated data from fictional studies. Results illustrate methods; they are not empirical findings.
:::
```

Because the data are simulated, you may compare estimates with the true values in `walk_truth.rds`, `clear_truth.rds`, `diary_truth.rds` and `kidney_truth.rds`. Label such comparisons clearly as possible only because the data are simulated. Truth for WALK, change from baseline at month 12: usual care −2.1, standard programme −6.7, enhanced programme −10.2 mmHg. The differences in mean change versus usual care are −4.6 (SP) and −8.1 (EP).

Reference models for WALK, to be used consistently across modules:

- **Linear-time model:** `lmer(sbp ~ month * arm + (1 + month | id), data = walk, control = lmer_ctrl)`. Time is actual months since randomisation.
- **Primary model (categorical time):** `lmer(sbp ~ visit_f * arm + (1 + month | id), data = walk, control = lmer_ctrl)`. Scheduled visit is a factor for the mean, and there is a random slope on actual time. The primary estimand is the **difference in mean change from baseline to month 12** between each programme and usual care, obtained with `emmeans`.
- A covariate-adjusted version adds `+ age10 + sex + bmi_c`, defined as `age10 = (age - 60)/10` and `bmi_c = bmi - 29`.

You may of course fit other models for teaching purposes (gls with corExp, splines and so on). Call them by descriptive names.

## 4. Structure of course notes

A module's notes are a substantial teaching text, not an outline: typically 6,000–10,000 words including code. Use connected paragraphs. Use bullet lists only where the content really is a list. Open with:

1. A short paragraph saying what the module is about and why it matters. Do not start with "In this module we will delve…".
2. **Learning outcomes**, copied from `syllabus/syllabus.qmd` with their numbers.
3. **Prerequisites** (which modules).
4. The simulated-data callout.

For **every method**, cover all ten elements below in this order. Use these as `###` subheadings (adapt the wording slightly if needed) so that readers find them in the same place every time:

1. *The research question and estimand*
2. *Data structure and outcome types*
3. *Intuition*
4. *The model*: the equation in display maths, followed by an explanation of every symbol
5. *Assumptions and limitations*
6. *Fitting the model in R*
7. *Interpreting the output*: coefficients, uncertainty, predictions, with annotated output
8. *Checking the model*
9. *Reporting*: a short model paragraph in the style of a paper, with real numbers from inline R
10. *Common mistakes*: a callout listing errors and how to avoid them

End every module with:

- **Summary**: a short paragraph, not a repetition of the outcomes.
- **Further reading**: 2–5 verified references with `[@key]`.
- **Exercises**: point to the workbook.

Label material:

- `## Heading <span class="core-tag">CORE</span>` for core sections.
- `## Heading <span class="advanced-tag">ADVANCED</span>` for advanced sections. Add one sentence explaining when the specialist method is needed.

Callout conventions:

- `.callout-important` titled "Common mistakes"
- `.callout-tip` titled "Reporting example"
- `.callout-note` titled "Advanced extension" or "Why this matters"
- `.callout-warning` for the simulated-data notice and for genuine pitfalls

**Annotated R output.** Show the real output of a model (for example `summary()`). Then explain it piece by piece in prose that refers to specific numbers, using inline R such as `` `r round(fixef(m)["month:armEP"], 2)` `` so that the text always matches the output. Never type a number into the text by hand if it comes from a model. A useful pattern is a short table listing "Output line → What it means → Value here".

**Figures.** Use ggplot2 with the course theme (already set). Use `scale_colour_arm()` for arms. Label axes with units. Give each figure a `#| fig-cap:` and keep the default sizes from `_metadata.yml` unless a figure needs more space. Prefer simulated demonstrations and plots to abstract description. **Diagrams** may be drawn with ggplot2 (for example data-structure schematics or covariance heatmaps) or written as Mermaid code blocks (```` ```{mermaid} ````). Mermaid renders in HTML.

## 5. Notation (use exactly this throughout)

| Symbol | Meaning |
|---|---|
| $i = 1,\dots,N$ | Subject (participant); $N$ subjects |
| $j = 1,\dots,n_i$ | Measurement occasion within subject $i$ |
| $y_{ij}$ | Outcome for subject $i$ at occasion $j$; $\mathbf{y}_i$ the vector for subject $i$ |
| $t_{ij}$ | Time of occasion $j$ for subject $i$ (actual time unless stated) |
| $\mathbf{x}_{ij}$, $X_i$ | Covariates / fixed-effects design matrix |
| $\boldsymbol\beta$, $\beta_k$ | Fixed effects |
| $b_{0i}$, $b_{1i}$, $\mathbf{b}_i$ | Random intercept, random slope, vector of random effects |
| $Z_i$ | Random-effects design matrix |
| $G$ | Covariance matrix of $\mathbf{b}_i$; elements $\tau_0^2$, $\tau_1^2$, $\tau_{01}$ |
| $\varepsilon_{ij}$, $\sigma^2$ | Residual (within-subject) error and its variance |
| $R_i$ | Residual covariance matrix for subject $i$ |
| $V_i = Z_i G Z_i^\top + R_i$ | Marginal covariance of $\mathbf{y}_i$ |
| $\rho$ | Correlation parameter (for example AR(1)); ICC written $\rho_I$ or "ICC" |
| $\mu_{ij}$ | Marginal mean $E(Y_{ij}\mid \mathbf{x}_{ij})$ (GEE) |
| $g(\cdot)$ | Link function |
| $\phi$ | Dispersion/scale parameter |
| $R(\boldsymbol\alpha)$ | GEE working correlation matrix |

Treatment indicators are written $\text{SP}_i$, $\text{EP}_i$ (WALK) or $\text{Active}_i$. Use $\log$ for natural logarithms.

## 6. Terminology (use these terms consistently)

- **participant** for a person in a study; **subject** in general statistical statements; **cluster** in GEE contexts (and explain that a participant is a cluster).
- **observation** = one row/one measurement; **visit** = a scheduled contact; **measurement occasion** = the time at which a measurement was taken. Distinguish **scheduled time** (`month_sched`) from **actual time** (`month`).
- **long format** / **wide format**.
- **within-subject (within-person)** vs **between-subject (between-person)**.
- **fixed effects**, **random effects**, **residual error**, **variance components**, **intraclass correlation (ICC)**.
- **random intercept**, **random slope**.
- **residual correlation structure** (compound symmetry, AR(1), exponential/continuous-time, unstructured); **marginal covariance**.
- **subject-specific (conditional)** vs **population-average (marginal)**. Use both words the first time in each module.
- **working correlation**, **robust (sandwich) standard error**, **model-based standard error**.
- **estimated marginal means (EMMs)**; **difference in mean change**; **contrast**.
- **confidence interval** vs **prediction interval**.
- **REML** (restricted maximum likelihood) and **ML** (maximum likelihood).
- **MCAR, MAR, MNAR**; **dropout** (monotone) vs **intermittent missingness**; **delta adjustment**.
- WALK arms: **usual care (UC)**, **standard programme (SP)**, **enhanced programme (EP)**. CLEAR/DIARY arms: **placebo**, **active**.
- British spelling: programme, randomised, modelling, analyse, centring, behaviour. Write "generalized" in model names (generalized linear mixed model, generalized estimating equations), as is standard in the literature.

## 7. Writing style

Write naturally and precisely for an intelligent reader who is not yet an expert.

- Lead with the point. Explain why before how.
- Use plain verbs and short-to-medium sentences.
- Develop ideas in paragraphs, not fragmented bullets.
- Explain every technical term at first use.
- **Do not use em dashes (—)**. Use full stops, commas, colons or parentheses. En dashes are fine in ranges (Modules 4–6), and minus signs in numbers.
- Avoid stock phrases ("delve", "pivotal", "leverage", "seamless", "it is important to note", "robust" except for robust standard errors, "unlock", "landscape").
- Avoid "not just X but Y" constructions, rhetorical questions that you then answer, and closing sentences that merely restate the paragraph.
- No emoji, no decorative arrows. Use bold sparingly.
- Never invent results, outputs or citations. Cite only keys in `references/references.bib`, which have been verified. If you need a reference that is not there, mention it in your report instead of citing it.

Available citation keys:

- **Books and handbooks:** verbeke2000, molenberghs2005, molenberghs2017 (Wiley StatsRef article "Longitudinal Data Analysis"), fitzmaurice2008 (handbook), fitzmaurice2011, diggle2002, singer2003, hedeker2006, pinheiro2000, gelman2007, wood2017, rizopoulos2012, hernan2020
- **Missing data:** little2019, vanbuuren2018, carpenter2013, molenberghs2007, daniels2008, nrc2010
- **Classic papers:** laird1982, liang1986, zeger1988, neuhaus1991, rubin1976, kenward1997, mancl2001, robins1995, self1987, stram1994, lindstrom1990, greenhouse1959, huynh1976, heinze2002
- **R software:** bates2015 (lme4), kuznetsova2017 (lmerTest), brooks2017 (glmmTMB), halekoh2006 (geepack), vanbuuren2011 (mice), lenth2026 (emmeans), hartig2026 (DHARMa), christensen2026 (ordinal), rizopoulos2026 (JMbayes2), ushey2026 (renv), wickham2014 (tidy data)
- **Applied guidance:** curran2011, enders2007, barr2013, matuschek2017, senn2006, vickers2001, twisk2018, mallinckrodt2008, pullenayegum2016, morris2019, arnold2011, gamble2017, ich2019

## 8. Slides

- **Format and length:** reveal.js, 16:9 (set in `_metadata.yml`). About 25–45 slides per module.
- **One message per slide:** state it in the slide title as a sentence where possible, for example "Random slopes let each participant change at their own rate".
- **Content:** concise text, figures, small tables and short code. Put extended explanation in instructor notes: `::: {.notes}` … `:::` on every content slide. Notes should tell the instructor what to say, what to point at, likely questions and timing-free activity prompts.
- **Readability:** keep code chunks short (no more than about 12 lines) and make sure figures are legible.
- **Sections:** separate them with `#` section slides. Use `.small` class text sparingly. Label advanced content with the advanced tag.
- **Exercises:** end each section with a "Try it" slide pointing to the matching workbook exercise.
- **Consistency:** slides must execute their own code or reuse the same model specifications as the notes, so the numbers match the notes.

## 9. Workbook and answer key

- **Workbook (`wbNN`):** a short introduction, then **Guided exercises** with step-by-step instructions, partially completed code (`#| eval: false`, using `___` blanks) and questions to answer in words. Then **Independent exercises**, which give a task with no scaffolding. Typically there are 4–6 guided and 3–5 independent exercises per module. Every exercise lists the learning outcomes it addresses (for example "LO 4.3").
- **Answer key (`akNN`):** exercise for exercise, with complete executed code, the real output, and model answers to the written questions that use inline R values. Add marking notes on acceptable alternatives and common errors. The answer key must render without errors, which is the verification that the outputs are real.

## 10. Assessments

- **Knowledge check (`kcNN`):** 10–15 questions, a mix of multiple choice, true/false-with-justification, and short interpretation of a given output (real output produced in the document). Map each question to a learning outcome. Put answers with explanations in a final section headed "Answers and explanations (for instructors)".
- **Practical assessment (`paNN`, Modules 2, 4, 6, 8, 10 only):** a task brief, the data, the deliverables, and a marking guide with points per criterion. A model solution can go at the end in a section marked for instructors, executed in the document.

## 11. Scripts

Each script is runnable top to bottom from the project root. It starts with a header comment (purpose, inputs, outputs, module), `source(here::here("R", "setup.R"))`, and clearly titled sections (`# ---- section ----`). Comment the *why*, not just the *what*. It saves the key figures, tables and models with `save_fig()`, `save_table()` and `save_model()`, using names prefixed `mNN_`. Use relative paths via `here()` only. The scripts and notes must use the same model specifications so that the numbers agree.

## 12. Inputs for cheat sheets and glossary

For each module, write:

- `admin/cheatsheet_input/mNN.md`. For each core method in the module, give:
  - its purpose
  - key assumptions
  - a short R example (no more than 6 lines)
  - the quantities to report
  - one common error

  Also include any decision points relevant to the method-selection diagrams, in plain text.
- `admin/glossary_input/mNN.md`: terms introduced in the module, each with a one- or two-sentence definition consistent with this guide.

## 13. Verification (required)

1. Run your scripts and render every `.qmd` you wrote. Fix all errors.
2. Read the rendered output. Check that every interpretive statement (direction, size, significance) agrees with the actual numbers, and that inline R values appear correctly.
3. Search your files for "—" (em dash) and remove any you introduced.
4. Keep R warnings visible when they teach something (for example a singular fit). Otherwise resolve them, or explain them in the text.
5. In your final report, list:
   - the files created;
   - the render and run status of each;
   - anything that could not be run or verified, and why;
   - any deviation from this guide.
