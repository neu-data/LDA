# Completion checklist

**Status as of 23 November 2026.** "Rendered" means the Quarto document executed all its R code and produced HTML without errors on the build machine (Windows 11, R 4.6.1, Quarto 1.10.18). "Runs" means the script ran from top to bottom with `Rscript` from the project root.

## Deliverables

| Deliverable | Location | Status |
|---|---|---|
| Syllabus with measurable module-level learning outcomes | `syllabus/syllabus.qmd` | Rendered |
| Course home page linking all materials | `index.qmd` | Rendered; all 87 links checked |
| README (installation, organisation, execution order) | `README.md` | Written |
| Style and terminology guide | `admin/style_guide.md` | Written; used by all modules |
| Simulated datasets with fixed seeds | `data-raw/01`–`05_*.R`, `data/*.csv` | Run; regeneration verified byte-identical (`scripts/00_make_data.R`) |
| Data dictionaries | `data/dictionaries/data_dictionary.md` | Written (WALK, CLEAR, DIARY, KIDNEY, MOOD) |
| Reusable functions | `R/*.R` | Used throughout |
| Package versions | `renv.lock`, `DESCRIPTION`, `outputs/*session_info*.txt` | Recorded |
| Verified references | `references/references.bib`, `verification_report.md`, `references.qmd` | 57 entries checked against Crossref/publisher/CRAN; rendered |
| Glossary | `glossary/glossary.qmd` | Rendered; 268 terms merged from module inputs |
| Cheat sheets (4) | `cheatsheets/cs01`–`cs04` | Rendered; every R example executed during rendering |
| Method-selection diagrams (6) | `diagrams/D1`–`D6.mmd`, `decision-diagrams.qmd` | Editable Mermaid sources; page rendered and inspected |
| Final project brief and rubric | `final-project/final-project-brief.qmd` | Rendered |
| SAP template | `final-project/sap_template.md` | Written |
| Final project instructor model solution | `final-project/instructor-solution.qmd` | Rendered |
| Analysis log template; example WALK SAP | `admin/templates/` | Written |
| Printable PDFs (notes, slides, cheat sheets, diagrams, glossary, syllabus, brief, references) | `pdf/` | Printed from the rendered HTML with headless Edge |

## Modules

Every module has notes, slides with instructor notes, a workbook, an answer key and a knowledge check. Practical assessments are provided for Modules 2, 4, 6, 8 and 10.

| Module | Scripts | Notes | Slides | Workbook | Answer key | Knowledge check | Practical |
|---|---|---|---|---|---|---|---|
| 0 Refresher | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | n/a |
| 1 Foundations | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | n/a |
| 2 Data preparation | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | Rendered |
| 3 Traditional methods | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | n/a |
| 4 Linear mixed models | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | Rendered |
| 5 Covariance models | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | n/a |
| 6 Interpretation and prediction | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | Rendered |
| 7 GEE | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | n/a |
| 8 GLMM | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | Rendered |
| 9 Time-varying and nonlinear | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | n/a |
| 10 Missing data | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | Rendered |
| 11 Diagnostics | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | n/a |
| 12 Planning and reporting | Runs | Rendered | Rendered | Rendered | Rendered | Rendered | n/a |
| 13 Advanced | Runs (13d joint model about 11 min) | Rendered | Rendered | Rendered | Rendered | Rendered | n/a |

## Verification carried out

* **Numbers in prose.** The notes, answer keys and assessments compute their numbers with inline R, so the text matches the output whenever a document is rendered. Module authors read the rendered text and corrected interpretive statements that did not match the numbers.
* **Consistency of the primary analysis across modules.** Modules 1, 4, 5, 6, 9, 11 and 12 report the same WALK primary-model results. The difference in mean change at month 12 is −5.38 (95% CI −8.17 to −2.58) for SP vs UC and −8.97 (−11.79 to −6.15) for EP vs UC. The true values are −4.6 and −8.1. The GEE estimates in Module 7 (−5.40, −8.98) and the multiple-imputation estimates in Module 10 (−5.33, −8.98) agree.
* **Slide layout.** Every deck was printed to PDF at 1280 × 720 (`pdf/slides/`), and each printed with exactly one page per slide, so no slide content overflows onto an extra page. Before that, an automated in-browser measurement flagged three figures or tables that were too tall (Module 3 sphericity heat map, Module 5 semivariogram and AIC table), and these were resized.
* **Style.**
  * Searched all course files for em dashes: none remain.
  * Searched for placeholders: one in the Module 5 answer key was replaced by a computed statement.
* **References.**
  * Only keys in the verified `.bib` file are cited.
  * Package citations reflect CRAN as of 26 October 2026.
* **Reproducibility.**
  * Data regeneration reproduces identical files.
  * All paths are relative (`here::here()`).

## Untested or partly verified

1. **Other platforms.** Everything was run on one Windows 11 machine with R 4.6.1. Not tested on macOS, Linux or earlier R versions.
2. **Slide appearance.** Overflow was checked automatically and through the PDF page counts. Individual slides were not all inspected by eye for aesthetics (for example the legibility of small axis labels), so review dense slides before teaching.
3. **PDF copies.** The PDFs in `pdf/` were produced by printing the HTML in headless Edge. Their page count and file size were checked, but they were not proof-read page by page.
4. **Chapter pointers.** Chapter numbers in the "Further reading" sections (for example "Diggle et al. 2002, chapters 1–2") come from the books' known structure but were not re-checked against each book's table of contents. The references themselves are verified.
5. **Two bibliographic details.**
   * The page range of Wickham (2014) is omitted.
   * The last page of Neuhaus et al. (1991) comes from secondary sources.
6. **Whole-project render.** Documents were rendered one by one. A single `quarto render` of the whole project from empty caches was not timed.
7. **Saved model outputs.** Module 13 slides, notes and the knowledge check read the joint model saved by `scripts/13d_joint_model.R` (about 11 minutes of MCMC). Run that script before rendering the Module 13 slides.
8. **Simulation-based results.** Several documents read simulation results cached in `outputs/models/` by the scripts, or re-run them if they are missing (Modules 1, 10, 11, 12). The knowledge check for Module 1 (`kc01`) needs `outputs/tables/m01_sim_summary.csv`, so run `scripts/01b_independence_simulation.R` first.

## Known software issues documented in the materials

* **geepack 1.3.13** crashes R (segmentation fault) with `corstr = "unstructured"` plus `waves` on incomplete clusters. Module 7 uses a hand-coded unstructured GEE for the Gaussian case instead.
* **A character cluster `id` in `geeglm()`** is silently coerced to `NA`, so all rows become one cluster. All materials use `factor(id)`, and Module 7 teaches this pitfall.
* **mice 3.19.0.**
  * `2l.lmer` fails with a random slope when some participants have one observed visit, so Module 10 uses `2l.norm`.
  * `pool()` for `lmer` fits needs broom.mixed, which is not a course package. Rubin's rules are applied with `pool.scalar()`.
* **DHARMa 0.5.0** simulates conditional on fitted random effects by default. Participant-level grouped tests then reject correct models, so the materials use `simulateREs = "unconditional"` for grouped checks (Modules 8 and 11).
* **`lme4::allFit()`** cannot refit models whose call uses the `lmer_ctrl` object. Write `control = lmerControl(optimizer = "bobyqa")` in full when you plan to use `allFit()` (noted in `R/setup.R` and Module 11).
* **`mgcv::gam()` with participant random effects as `bs = "re"` smooths** is very slow on WALK (about 9 minutes). Module 9 uses `gamm()` or `bam()`.

## Decisions and notes for instructors

* "Molenberghs and Verbeke 2017" was identified as their article "Longitudinal Data Analysis" in *Wiley StatsRef: Statistics Reference Online* (doi:10.1002/9781118445112.stat00291.pub2). No 2017 book edition by these authors was found. Their books (2000, 2005) and the 2008 handbook they co-edited are also cited. If a different 2017 source was intended, add it to `references/references.bib`.
* **Folder size.** The project is about 1 GB, almost all of it `outputs/models/m13_jm.rds` (about 800 MB of MCMC samples for the Module 13 joint model). To share a smaller copy, leave this file out; recipients regenerate it by running `scripts/13d_joint_model.R` (about 11 minutes) before rendering the Module 13 slides.
* The WALK simulation seed was chosen from 15 candidates so that the estimates are typical of their sampling distribution. The choice is documented in `data-raw/01_simulate_walk.R`.
* **Practical assessments.** The practical assessments `pa02` and `pa10` include the instructor model solution in the same file. Remove the sections marked for instructors before giving the brief to participants. Both write their own datasets to `assessments/pa02-data/` and `assessments/pa10-data/`.
* **Observed correlations.** Module 5 computes observed correlations as `cov2cor()` of the pairwise-complete covariance. Other modules use `visit_correlation()` (pairwise `cor()`). The two differ slightly (for example 0.85 vs 0.88 for months 0 and 1). Both are legitimate summaries.
* **Hand-typed numbers.** A few numbers in slide instructor notes (Modules 3 and 4) are typed by hand. They match the output at the time of writing, but they will not update automatically if the data are regenerated with a different seed.
* **Final project data.** The MOOD data include dropout that depends partly on unobserved outcomes, so a MAR analysis is somewhat optimistic. This is deliberate (see the instructor guidance in the brief).
