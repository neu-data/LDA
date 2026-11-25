# =============================================================================
# 07b_gee_simulations.R
# Purpose : Two simulation studies for Module 7.
#           (1) Few clusters: how well do model-based, sandwich and
#               Mancl-DeRouen-corrected standard errors work with 10, 40 and
#               200 clusters? Reports Monte Carlo standard errors.
#           (2) Dropout in a CLEAR-like trial with 20,000 participants:
#               bias of unweighted GEE under MAR dropout, and its removal by
#               inverse-probability weighting.
#           Also checks the hand-coded Mancl-DeRouen variance against the
#           leave-one-cluster-out jackknife (exact for a linear model).
# Inputs  : data/walk_long.csv, data/clear_truth.rds
# Outputs : outputs/tables/m07_small_cluster_sim.csv,
#           outputs/tables/m07_small_cluster_summary.csv,
#           outputs/tables/m07_ipw_large_sample.csv,
#           outputs/figures/m07_small_cluster_coverage.png
# Module  : 7 (Generalized estimating equations)
# Run time: about 2-3 minutes.
# =============================================================================

source(here::here("R", "setup.R"))

# ---- helper: sandwich variance with optional Mancl-DeRouen correction -------
# (identical to the function in scripts/07a_gee.R)
sandwich_gee <- function(fit, correction = c("none", "MD")) {
  correction <- match.arg(correction)
  X <- model.matrix(fit); y <- fit$y; mu <- fitted(fit)
  fam <- family(fit); eta <- fit$linear.predictors
  w <- fit$prior.weights; if (is.null(w)) w <- rep(1, length(y))
  alpha <- fit$geese$alpha
  cl <- split(seq_along(y), fit$id)
  p <- ncol(X); B <- matrix(0, p, p); parts <- vector("list", length(cl))
  for (k in seq_along(cl)) {
    r <- cl[[k]]; n <- length(r); R <- diag(n)
    if (fit$corstr == "exchangeable" && n > 1) { R[] <- alpha; diag(R) <- 1 }
    if (fit$corstr == "ar1" && n > 1) R <- alpha^abs(outer(1:n, 1:n, "-"))
    A <- diag(sqrt(fam$variance(mu[r]) / w[r]), n)
    Vinv <- solve(A %*% R %*% A)
    D <- fam$mu.eta(eta[r]) * X[r, , drop = FALSE]
    B <- B + t(D) %*% Vinv %*% D
    parts[[k]] <- list(D = D, Vinv = Vinv, e = y[r] - mu[r])
  }
  Binv <- solve(B); M <- matrix(0, p, p)
  for (pk in parts) {
    e <- pk$e
    if (correction == "MD") {
      H <- pk$D %*% Binv %*% t(pk$D) %*% pk$Vinv
      e <- solve(diag(nrow(H)) - H, e)
    }
    u <- t(pk$D) %*% pk$Vinv %*% e
    M <- M + u %*% t(u)
  }
  V <- Binv %*% M %*% Binv
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

# ---- check: Mancl-DeRouen equals the jackknife for a linear model -----------
# For an identity link with independence working correlation, the MD
# variance equals sum_i (beta_(-i) - beta)(beta_(-i) - beta)', where
# beta_(-i) leaves out cluster i. This verifies the hand-coded formula.
walk_sub <- load_walk() |> filter(clinic %in% 1:2) |> arrange(id, month)
g_ind <- geeglm(sbp ~ month * arm, id = factor(id), data = walk_sub,
                corstr = "independence")
jack <- sapply(unique(walk_sub$id), function(i)
  coef(lm(sbp ~ month * arm, data = filter(walk_sub, id != i))) - coef(g_ind))
check_md <- max(abs(sandwich_gee(g_ind, "MD") - jack %*% t(jack)))
check_sw <- max(abs(sandwich_gee(g_ind) - vcov(g_ind)))
cat("Max |MD - jackknife| =", signif(check_md, 3),
    "; max |hand sandwich - geeglm| =", signif(check_sw, 3), "\n")

# ---- simulation 1: few clusters ---------------------------------------------
# Continuous outcome, 5 equally spaced occasions, half the clusters active.
# The truth has a random intercept AND a random slope, so the exchangeable
# working correlation is (mildly) wrong: model-based SEs are not trustworthy,
# and only the sandwich is justified. Target: time:active = -0.5.
sim_small_clusters <- function(n_clusters, times = 0:4, effect = -0.5,
                               sd_b0 = 2, sd_b1 = 0.5, sd_e = 1.5) {
  d <- expand.grid(time = times, id = seq_len(n_clusters))
  d$active <- as.integer(d$id > n_clusters / 2)
  b0 <- rnorm(n_clusters, 0, sd_b0); b1 <- rnorm(n_clusters, 0, sd_b1)
  d$y <- -0.5 * d$time + effect * d$active * d$time +
    b0[d$id] + b1[d$id] * d$time + rnorm(nrow(d), 0, sd_e)
  fit <- geeglm(y ~ time * active, id = id, data = d, corstr = "exchangeable")
  k <- "time:active"
  c(estimate = unname(coef(fit)[k]),
    se_model = sqrt(fit$geese$vbeta.naiv[4, 4]),
    se_sandwich = sqrt(vcov(fit)[k, k]),
    se_md = sqrt(sandwich_gee(fit, "MD")[k, k]))
}

set.seed(2027)
n_sim <- 1000
sims <- bind_rows(lapply(c(10, 40, 200), function(N)
  data.frame(N = N, t(replicate(n_sim, sim_small_clusters(N))))))
save_table(sims, "m07_small_cluster_sim")

covers <- function(est, se, q = qnorm(0.975)) abs(est - (-0.5)) <= q * se
sim_summary <- sims |>
  group_by(N) |>
  summarise(empirical_sd = sd(estimate),
            mean_se_model = mean(se_model),
            mean_se_sandwich = mean(se_sandwich),
            mean_se_md = mean(se_md),
            cover_model = mean(covers(estimate, se_model)),
            cover_sandwich = mean(covers(estimate, se_sandwich)),
            cover_md = mean(covers(estimate, se_md)),
            cover_md_t = mean(covers(estimate, se_md, qt(0.975, N - 4)))) |>
  mutate(mcse_cover = sqrt(0.95 * 0.05 / n_sim),
         mcse_empirical_sd = empirical_sd / sqrt(2 * (n_sim - 1)))
print(as.data.frame(sim_summary))
save_table(sim_summary, "m07_small_cluster_summary")

cov_long <- sim_summary |>
  select(N, cover_model, cover_sandwich, cover_md, cover_md_t) |>
  pivot_longer(-N, names_to = "method", values_to = "coverage") |>
  mutate(method = factor(method,
                         levels = c("cover_model", "cover_sandwich",
                                    "cover_md", "cover_md_t"),
                         labels = c("Model-based", "Sandwich",
                                    "Mancl-DeRouen (z)", "Mancl-DeRouen (t, N-4 df)")),
         lo = coverage - 1.96 * sqrt(coverage * (1 - coverage) / n_sim),
         hi = coverage + 1.96 * sqrt(coverage * (1 - coverage) / n_sim))
p_cov <- ggplot(cov_long, aes(factor(N), coverage, colour = method)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = unname(course_colours[c("grey", "blue",
                                                       "orange", "green")]),
                      name = NULL) +
  labs(x = "Number of clusters", y = "Coverage of nominal 95% CI")
save_fig(p_cov, "m07_small_cluster_coverage")

# ---- simulation 2: IPW-GEE in a large CLEAR-like trial ----------------------
truth_clear <- readRDS(here("data", "clear_truth.rds"))
simulate_clear_like <- function(n, truth) {
  part <- tibble(id = seq_len(n), active = rep(0:1, length.out = n),
                 age = round(pmin(pmax(rnorm(n, 45, 13), 18), 80)),
                 female = as.integer(runif(n) < 0.55),
                 b = rnorm(n, 0, truth$sd_b))
  d <- part |> crossing(visit = 1:6) |> arrange(id, visit) |>
    mutate(week = truth$weeks[visit],
           lp = -truth$cuts[2] + b + truth$beta_week * week +
             truth$beta_trtweek * active * week +
             truth$beta_age * (age - 45) / 10 + truth$beta_female * female,
           severity = findInterval(lp + rlogis(n()), truth$cuts),
           itch = as.integer(severity >= 2),
           arm = factor(ifelse(active == 1, "Active", "Placebo"),
                        levels = c("Placebo", "Active")))
  S <- matrix(d$severity, ncol = 6, byrow = TRUE)
  present <- matrix(TRUE, n, 6); dropped <- rep(FALSE, n); last <- S[, 1]
  for (k in 2:6) {
    dropped <- dropped |
      (runif(n) < plogis(truth$drop_int + truth$drop_sev * last))
    present[, k] <- !dropped & runif(n) >= truth$p_miss_visit
    last <- ifelse(present[, k], S[, k], last)
  }
  d$present <- as.vector(t(present))
  d
}
build_ipw <- function(obs) {
  grid <- obs |> distinct(id, arm) |> crossing(visit = 1:6) |>
    left_join(select(obs, id, visit, severity), by = c("id", "visit")) |>
    group_by(id) |>
    mutate(last_visit = max(visit[!is.na(severity)]), last_sev = severity) |>
    tidyr::fill(last_sev, .direction = "down") |>
    mutate(prev_sev = lag(last_sev)) |> ungroup()
  risk <- grid |> filter(visit >= 2, visit <= pmin(last_visit + 1, 6)) |>
    mutate(stay = as.integer(visit <= last_visit))
  drop_model <- glm(stay ~ prev_sev + arm, family = binomial, data = risk)
  risk$p_stay <- fitted(drop_model)
  risk <- risk |> group_by(id) |> mutate(cum_p = cumprod(p_stay)) |> ungroup()
  out <- obs |>
    left_join(select(risk, id, visit, cum_p), by = c("id", "visit")) |>
    mutate(cum_p = ifelse(visit == 1, 1, cum_p), ipw = 1 / cum_p)
  list(data = out, model = drop_model, risk = risk)
}

set.seed(99)
big <- simulate_clear_like(20000, truth_clear)
g_full <- geeglm(itch ~ week * arm, id = id, data = big, family = binomial,
                 corstr = "independence")
big_ipw <- build_ipw(filter(big, present))
g_unw <- geeglm(itch ~ week * arm, id = id, data = big_ipw$data,
                family = binomial, corstr = "independence")
g_ipw <- geeglm(itch ~ week * arm, id = id, data = big_ipw$data,
                family = binomial, corstr = "independence", weights = ipw)
ipw_large <- data.frame(term = names(coef(g_full)),
                        complete_data = coef(g_full),
                        se_complete = sqrt(diag(vcov(g_full))),
                        observed_unweighted = coef(g_unw),
                        observed_ipw = coef(g_ipw))
print(ipw_large)
save_table(ipw_large, "m07_ipw_large_sample")
