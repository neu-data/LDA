# =============================================================================
# R/plot_helpers.R
# Reusable longitudinal plots and a consistent way of saving figures.
# =============================================================================

# Spaghetti plot: one line per subject, optionally faceted by group.
# `n_show` draws a random subset so that individual lines remain visible.
plot_spaghetti <- function(data, x, y, id, group = NULL, n_show = NULL,
                           seed = 1, alpha = 0.35) {
  d <- data
  if (!is.null(n_show)) {
    set.seed(seed)
    keep <- sample(unique(pull(d, {{ id }})), n_show)
    d <- filter(d, {{ id }} %in% keep)
  }
  p <- ggplot(d, aes({{ x }}, {{ y }}, group = {{ id }})) +
    geom_line(alpha = alpha, linewidth = 0.4, colour = "#555555") +
    geom_point(alpha = alpha, size = 0.8, colour = "#555555")
  if (!rlang::quo_is_null(rlang::enquo(group))) {
    p <- p + facet_wrap(vars({{ group }}))
  }
  p
}

# Observed means with 95% CI by group and visit (observed data only).
plot_observed_means <- function(data, x, y, group) {
  s <- data |>
    group_by({{ group }}, {{ x }}) |>
    summarise(mean = mean({{ y }}, na.rm = TRUE),
              se   = sd({{ y }}, na.rm = TRUE) / sqrt(sum(!is.na({{ y }}))),
              n    = sum(!is.na({{ y }})), .groups = "drop") |>
    mutate(lower = mean - qt(0.975, n - 1) * se,
           upper = mean + qt(0.975, n - 1) * se)
  ggplot(s, aes({{ x }}, mean, colour = {{ group }}, group = {{ group }})) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25,
                  linewidth = 0.5) +
    scale_colour_arm(name = NULL)
}

# Missingness map: subjects (rows) by scheduled visits (columns).
plot_missing_map <- function(complete_grid, id, visit, observed, order_by) {
  ggplot(complete_grid,
         aes(factor({{ visit }}), reorder({{ id }}, {{ order_by }}),
             fill = {{ observed }})) +
    geom_tile() +
    scale_fill_manual(values = c(`TRUE` = "#4D4D4D", `FALSE` = "#F0F0F0"),
                      labels = c(`TRUE` = "Observed", `FALSE` = "Missing"),
                      name = NULL) +
    labs(x = "Scheduled visit", y = "Participant") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank())
}

# Save a figure to outputs/figures with standard dimensions.
# Widths are chosen to stay legible on A4 notes and 16:9 slides.
save_fig <- function(plot, name, width = 7, height = 4.5, dpi = 200) {
  dir.create(here("outputs", "figures"), showWarnings = FALSE,
             recursive = TRUE)
  ggsave(here("outputs", "figures", paste0(name, ".png")), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  invisible(here("outputs", "figures", paste0(name, ".png")))
}

save_table <- function(x, name) {
  dir.create(here("outputs", "tables"), showWarnings = FALSE,
             recursive = TRUE)
  write_csv(as.data.frame(x), here("outputs", "tables", paste0(name, ".csv")))
  invisible(x)
}

save_model <- function(model, name) {
  dir.create(here("outputs", "models"), showWarnings = FALSE,
             recursive = TRUE)
  saveRDS(model, here("outputs", "models", paste0(name, ".rds")))
  invisible(model)
}
