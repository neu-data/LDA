# =============================================================================
# R/setup.R
# Loads the course packages and helper functions. Every script, note and
# exercise starts with:  source(here::here("R", "setup.R"))
#
# The course uses a deliberately small set of packages. Why each is needed:
#
#   Data handling and graphics
#     here        builds paths relative to the project root (no setwd()).
#     readr       reads and writes CSV files with predictable column types.
#     dplyr       filters, groups and summarises data frames.
#     tidyr       reshapes between long and wide formats.
#     ggplot2     all figures.
#     patchwork   places several ggplot figures side by side.
#
#   Models
#     lme4        linear and generalised linear mixed models (lmer, glmer).
#     lmerTest    Satterthwaite degrees of freedom and p-values for lmer.
#     nlme        residual correlation structures (gls, lme) and nonlinear
#                 mixed models (nlme).
#     glmmTMB     negative binomial and zero-inflated mixed models.
#     geepack     generalised estimating equations (geeglm).
#     emmeans     estimated marginal means, contrasts and predictions.
#     mgcv        generalised additive (mixed) models.
#     mice        multiple imputation, including multilevel methods.
#     DHARMa      simulation-based residual checks for GLMMs.
#
#   Advanced module only (loaded inside the scripts that need them)
#     ordinal     cumulative-link mixed models for ordinal outcomes.
#     survival    survival objects and Cox models.
#     JMbayes2    joint models for longitudinal and time-to-event data.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(lme4)
  library(lmerTest)
  library(nlme, exclude = "lmList")   # avoid masking lme4::lmList
  library(glmmTMB)
  library(geepack)
  library(emmeans)
  library(mgcv)
  library(mice, exclude = c("filter", "complete"))
  library(DHARMa)
})

# nlme and lme4 both export some generics; make the lme4 versions explicit
fixef   <- lme4::fixef
ranef   <- lme4::ranef
VarCorr <- lme4::VarCorr

options(
  dplyr.summarise.inform = FALSE,
  digits = 4,
  pillar.sigfig = 4,
  emmeans = list(lmer.df = "satterthwaite")
)

# Optimiser settings used for every lmer() fit in the course. The minqa
# BOBYQA optimiser converges cleanly for the WALK models; see Module 11 for
# how to investigate convergence warnings rather than simply silencing them.
# Note: lme4::allFit() cannot re-evaluate a model whose call refers to
# `lmer_ctrl`; when you plan to use allFit(), write the control argument out
# in full: control = lmerControl(optimizer = "bobyqa").
lmer_ctrl <- lmerControl(optimizer = "bobyqa")

for (f in c("theme_course.R", "data_helpers.R", "model_helpers.R",
            "plot_helpers.R", "sim_helpers.R")) {
  source(here("R", f))
}
rm(f)
