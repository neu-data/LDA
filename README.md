# Longitudinal Data Analysis: Concepts, Methods and Practical Applications Using R

These are complete course materials for researchers, analysts and learners with basic statistical knowledge. The course runs from the structure of longitudinal data through linear mixed models, covariance modelling, GEE and GLMMs to time-varying predictors, nonlinear change, missing data, diagnostics, study planning and reporting. An advanced module covers ordinal outcomes, three-level data, nonlinear mixed models, joint models and time-varying confounding.

The course is organised into progressive modules (0–13) without dates or fixed teaching hours. Start with `syllabus/syllabus.html` for the learning outcomes and suggested route, or open `index.html` for links to every material.

**All datasets are simulated.** The WALK, CLEAR, DIARY, KIDNEY and MOOD studies are fictional. Results illustrate methods and are not empirical findings.

## Online materials and the book

The rendered materials are published with GitHub Pages at <https://neu-data.github.io/LDA/>.

The course notes are also published as a book, *Longitudinal Data Analysis in R: Concepts, Methods and Practical Applications* (English edition), online at <https://neu-data.github.io/LDA/book/en/_book/> with a PDF download. The book's chapters are generated from `notes/`, so the notes remain the single source. To rebuild the book after editing the notes:

```
Rscript book/make_book_sources.R
cd book/en
quarto render
```

Chapter numbers match module numbers. The `book/en/` folder is the English edition; other language editions can be added as `book/<language>/`.

## Quick start

1. Install R (version 4.4 or later; the materials were run with **R 4.6.1** on Windows 11) and RStudio. RStudio includes Quarto. Quarto 1.4 or later is needed to re-render documents; the materials were rendered with Quarto 1.10.18.
2. Open `LDA_course.Rproj` in RStudio. This sets the working directory to the project root, and every path in the course is relative to it (via `here::here()`).
3. Install the packages (once):

   ```r
   install.packages(c(
     "here", "readr", "dplyr", "tidyr", "ggplot2", "patchwork",
     "lme4", "lmerTest", "pbkrtest", "nlme", "glmmTMB", "geepack", "emmeans",
     "mgcv", "mice", "DHARMa", "knitr", "rmarkdown",
     "ordinal", "survival", "JMbayes2"   # advanced module only
   ))
   ```

   To reproduce the exact package versions used to build the materials, use the lockfile instead: `install.packages("renv")`, then `renv::restore()`.
4. Check the set-up by running `source(here::here("R", "setup.R"))`. It should load without errors.
5. Open any rendered `.html` file in a web browser. Rendered files are self-contained, so they need no internet connection and can be shared as single files.

## Why each package is needed

| Package | Purpose in the course |
|---|---|
| here, readr, dplyr, tidyr | Relative paths, reading data, data manipulation, reshaping long and wide |
| ggplot2, patchwork | All figures; multi-panel layouts |
| lme4, lmerTest, pbkrtest | Linear and generalized linear mixed models; Satterthwaite and Kenward–Roger degrees of freedom |
| nlme | Residual correlation structures (`gls`, `lme`), heterogeneous variances, nonlinear mixed models |
| glmmTMB | Negative binomial and zero-inflated mixed models |
| geepack | Generalized estimating equations |
| emmeans | Estimated marginal means, contrasts and predictions |
| mgcv | Generalized additive (mixed) models |
| mice | Multiple imputation, including multilevel imputation |
| DHARMa | Simulation-based residual diagnostics for GLMMs |
| ordinal, survival, JMbayes2 | Advanced module only: ordinal mixed models, survival models, joint models |
| knitr, rmarkdown | Rendering the Quarto documents |

`R/setup.R` loads the core packages, sets the figure theme and defines the shared helper functions. Every script and document starts with `source(here::here("R", "setup.R"))`.

## File organisation

```
LDA_course.Rproj        RStudio project (open this first)
index.html              Course home page with links to all materials
README.md               This file
renv.lock, DESCRIPTION  Exact package versions used; package list for renv
_quarto.yml             Quarto project settings (each folder has a _metadata.yml)

syllabus/               Syllabus and measurable learning outcomes
notes/                  Course notes, one per module (NN-topic.qmd and .html)
slides/                 Reveal.js slides with instructor notes (press S for speaker view)
workbook/               Participant workbook: guided and independent exercises (wbNN)
answer-key/             Instructor answer key with executed, verified output (akNN)
assessments/            Knowledge checks (kcNN) and practical assessments (paNN)
cheatsheets/            Four printable reference sheets (cs01-cs04)
diagrams/               Six method-selection diagrams: editable Mermaid (.mmd) + rendered page
glossary/               Glossary of terms (268 entries)
pdf/                    Printable PDFs of notes, slides, cheat sheets, diagrams and reference documents
references/             Verified references (.bib), verification report, rendered list
final-project/          Final project brief and rubric, SAP template, instructor model solution

R/                      Reusable functions: setup, theme, data, model, plot, simulation helpers
data-raw/               Simulation scripts for every dataset (fixed seeds)
data/                   Simulated datasets (.csv), true values (.rds), data dictionary
scripts/                Numbered analysis scripts, one or more per module
outputs/                Figures, tables and fitted models saved by the scripts
admin/                  Style guide, completion checklist, templates (analysis log, example SAP)
```

## Execution order

The datasets are already in `data/`. To rebuild everything from scratch, run the steps below from the project root (in RStudio, or with `Rscript` in a terminal). The time given for each step is approximate, on the laptop used to build the course.

1. **Data** (optional; it reproduces `data/` exactly): `scripts/00_make_data.R` (under 1 min).
2. **Module scripts**, in numerical order. Each writes to `outputs/` and can also be run on its own.

| Script | Module | Time |
|---|---|---|
| `00_refresher.R` | 0 | < 1 min |
| `01a_foundations_structure.R`, `01b_independence_simulation.R` | 1 | < 1 min; about 5 min |
| `02_data_preparation.R` | 2 | < 1 min |
| `03_traditional_methods.R` | 3 | < 1 min |
| `04_linear_mixed_models.R` | 4 | about 2 min |
| `05_covariance_models.R` | 5 | < 1 min |
| `06_interpretation_prediction.R` | 6 | < 1 min |
| `07a_gee.R`, `07b_gee_simulations.R` | 7 | 1 min; about 4 min |
| `08a_glmm_binary.R`, `08b_glmm_count.R` | 8 | about 2 min each |
| `09_time_varying_nonlinear.R` | 9 | about 1–2 min |
| `10a_missing_data.R`, `10b_missing_data_simulations.R` | 10 | about 6 min; 15–20 min |
| `11_diagnostics.R` | 11 | about 10–15 min (360 leave-one-out refits) |
| `12_planning.R` | 12 | about 10 min (power simulations) |
| `13a_ordinal.R` … `13e_msm.R` | 13 | a few minutes; the joint model (13d) is the slowest |

3. **Documents.** Render with Quarto from the project root, for example `quarto render notes/04-linear-mixed-models.qmd`, or render everything with `quarto render`. Rendering every file takes a long time (about an hour), because the notes and answer keys run their analyses. Several documents reuse simulation results saved in `outputs/models/` by the module scripts, or re-run the simulation if those files are missing. Run the scripts first to keep rendering fast.

## Materials per module

Each module has course notes, slides, a workbook, an answer key and a knowledge check. Modules 2, 4, 6, 8 and 10 also have practical assessments. Material is labelled <span>CORE</span> or <span>ADVANCED</span> throughout. Module 0 is an optional refresher of R and regression. Module 13 is entirely advanced.

## Verification and known limitations

Every script was run and every document rendered without errors on the build machine. Numbers in the notes and answer keys are computed when the document is rendered (inline R), so the prose always matches the output. `admin/completion_checklist.md` lists each deliverable with its status and records anything that remains untested, together with known issues and decisions.

## Licence and attribution

Prepared for Neudata Consultancy training. Package and method references are in `references/`.
