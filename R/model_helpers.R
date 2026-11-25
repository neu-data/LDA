# =============================================================================
# R/model_helpers.R
# Small functions that extract and format quantities the course reports
# repeatedly. They wrap standard package functions; they do not replace them.
# =============================================================================

# Format an estimate with its confidence interval: "-6.1 (-8.4 to -3.8)"
fmt_ci <- function(est, lo, hi, digits = 1) {
  f <- function(x) formatC(x, format = "f", digits = digits)
  paste0(f(est), " (", f(lo), " to ", f(hi), ")")
}

# Format a p-value for text: "p < 0.001" or "p = 0.034"
fmt_p <- function(p) {
  ifelse(p < 0.001, "p < 0.001",
         paste0("p = ", formatC(p, format = "f", digits = 3)))
}

# Fixed-effects table for lmer (lmerTest), lme, gls, glmer, glmmTMB, geeglm.
# Wald confidence intervals; for lmerTest models the Satterthwaite df are
# used for the t-based interval.
fixef_table <- function(model, level = 0.95, exponentiate = FALSE) {
  z <- qnorm(1 - (1 - level) / 2)
  if (inherits(model, "lmerModLmerTest")) {
    s  <- coef(summary(model))
    df <- s[, "df"]
    q  <- qt(1 - (1 - level) / 2, df)
    out <- data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2],
                      df = df, lower = s[, 1] - q * s[, 2],
                      upper = s[, 1] + q * s[, 2], p = s[, 5])
  } else if (inherits(model, "glmmTMB")) {
    s <- coef(summary(model))$cond
    out <- data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2],
                      lower = s[, 1] - z * s[, 2], upper = s[, 1] + z * s[, 2],
                      p = s[, 4])
  } else if (inherits(model, "geeglm")) {
    s <- coef(summary(model))
    out <- data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2],
                      lower = s[, 1] - z * s[, 2], upper = s[, 1] + z * s[, 2],
                      p = s[, 4])
  } else if (inherits(model, c("lme", "gls"))) {
    s <- summary(model)$tTable
    dfs <- if ("DF" %in% colnames(s)) s[, "DF"] else
      rep(model$dims$N - model$dims$p, nrow(s))
    q <- qt(1 - (1 - level) / 2, dfs)
    out <- data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2],
                      df = dfs, lower = s[, 1] - q * s[, 2],
                      upper = s[, 1] + q * s[, 2], p = s[, ncol(s)])
  } else {
    s <- coef(summary(model))
    out <- data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2],
                      lower = s[, 1] - z * s[, 2], upper = s[, 1] + z * s[, 2],
                      p = s[, ncol(s)])
  }
  if (exponentiate) {
    out$estimate <- exp(out$estimate)
    out$lower    <- exp(out$lower)
    out$upper    <- exp(out$upper)
  }
  rownames(out) <- NULL
  out
}

# Variance components of an lmer fit as a tidy data frame
varcomp_table <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  vc$type <- ifelse(is.na(vc$var2), "SD", "correlation")
  vc$value <- ifelse(is.na(vc$var2), vc$sdcor, vc$sdcor)
  vc[, c("grp", "var1", "var2", "type", "value")]
}

# Intraclass correlation from a random-intercept lmer fit:
# between-subject variance / (between-subject + residual variance)
icc_lmer <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  v_b <- vc$vcov[vc$var1 == "(Intercept)" & is.na(vc$var2)][1]
  v_e <- vc$vcov[vc$grp == "Residual"]
  v_b / (v_b + v_e)
}

# Implied marginal covariance of a random-intercept/random-slope model with
# independent residuals, at times t:  V = Z G Z' + sigma^2 I
implied_cov_lmer <- function(model, times, group = 1) {
  G  <- as.matrix(lme4::VarCorr(model)[[group]])
  s2 <- sigma(model)^2
  Z  <- cbind(1, times)[, seq_len(ncol(G)), drop = FALSE]
  V  <- Z %*% G %*% t(Z) + diag(s2, length(times))
  dimnames(V) <- list(times, times)
  V
}

# Convert a covariance matrix to a correlation matrix with rounded output
cov_to_cor <- function(V, digits = 2) round(cov2cor(V), digits)

# Likelihood-ratio test between two nested fits, returning a one-row summary
lr_test <- function(m0, m1) {
  l0 <- logLik(m0); l1 <- logLik(m1)
  stat <- as.numeric(2 * (l1 - l0))
  df   <- attr(l1, "df") - attr(l0, "df")
  data.frame(logLik_0 = as.numeric(l0), logLik_1 = as.numeric(l1),
             chisq = stat, df = df,
             p = pchisq(stat, df, lower.tail = FALSE))
}

# p-value for a single variance component tested on the boundary
# (equal mixture of chi-square(0) and chi-square(1); Self & Liang 1987)
p_boundary_mixture <- function(chisq, df_mix = c(0, 1)) {
  if (identical(df_mix, c(0, 1))) {
    0.5 * pchisq(chisq, 1, lower.tail = FALSE)
  } else {
    0.5 * pchisq(chisq, df_mix[1], lower.tail = FALSE) +
      0.5 * pchisq(chisq, df_mix[2], lower.tail = FALSE)
  }
}
