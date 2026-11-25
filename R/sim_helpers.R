# =============================================================================
# R/sim_helpers.R
# Simulation tools used for demonstrations (Modules 1, 3, 5) and for
# simulation-based power and sample-size planning (Module 12).
# =============================================================================

# Simulate a simple two-arm longitudinal trial with a random intercept and
# random slope, linear mean trajectories and optional monotone dropout.
#
#   y_ij = b0 + (b_time + delta * trt_i) t_j + u0_i + u1_i t_j + e_ij
#
# Arguments
#   n_per_arm   participants per arm
#   times       measurement occasions (e.g. months)
#   delta       difference in mean slope between arms (the estimand)
#   sd_int, sd_slope, cor_int_slope, sd_resid  variance parameters
#   dropout     probability of dropping out before each post-baseline visit
sim_two_arm <- function(n_per_arm = 50, times = c(0, 3, 6, 9, 12),
                        b0 = 150, b_time = -0.2, delta = -0.5,
                        sd_int = 11, sd_slope = 0.4, cor_int_slope = -0.2,
                        sd_resid = 6, dropout = 0) {
  n   <- 2 * n_per_arm
  trt <- rep(0:1, each = n_per_arm)
  G   <- matrix(c(sd_int^2, cor_int_slope * sd_int * sd_slope,
                  cor_int_slope * sd_int * sd_slope, sd_slope^2), 2)
  u   <- MASS::mvrnorm(n, c(0, 0), G)
  d <- expand.grid(time = times, id = seq_len(n))
  d$trt <- trt[d$id]
  d$y <- b0 + (b_time + delta * d$trt) * d$time +
    u[d$id, 1] + u[d$id, 2] * d$time + rnorm(nrow(d), 0, sd_resid)
  if (dropout > 0) {
    # visit index at which each subject is last seen
    last <- 1 + rgeom(n, dropout)
    d <- d[match(d$time, times) <= last[d$id], ]
  }
  d$id <- factor(d$id)
  d
}

# Estimate power for the arm-by-time interaction by repeated simulation.
# Returns the proportion of simulated trials with p < alpha and the Monte
# Carlo standard error of that proportion (Morris et al. 2019).
power_two_arm <- function(n_sims = 200, alpha = 0.05, seed = 1, ...) {
  set.seed(seed)
  p <- vapply(seq_len(n_sims), function(i) {
    d <- sim_two_arm(...)
    m <- suppressMessages(suppressWarnings(
      lmerTest::lmer(y ~ time * trt + (time | id), data = d,
                     control = lme4::lmerControl(calc.derivs = FALSE))))
    coef(summary(m))["time:trt", "Pr(>|t|)"]
  }, numeric(1))
  pw <- mean(p < alpha)
  c(power = pw, mcse = sqrt(pw * (1 - pw) / n_sims), n_sims = n_sims)
}
