# =============================================================================
# R/data_helpers.R
# Reading the course datasets with consistent types, and checking the
# structure of long-format longitudinal data.
# =============================================================================

walk_arm_labels <- c(UC = "Usual care", SP = "Standard programme",
                     EP = "Enhanced programme")

# ---- loaders ----------------------------------------------------------------
# Each loader returns a tibble with factors set to the course reference
# levels (Usual care / Placebo / Male first) and useful derived variables.

load_walk <- function(which = c("long", "baseline", "wide", "complete")) {
  which <- match.arg(which)
  d <- read_csv(here("data", paste0("walk_", which, ".csv")),
                show_col_types = FALSE)
  d <- d |>
    mutate(
      arm    = factor(arm, levels = c("UC", "SP", "EP")),
      sex    = factor(sex, levels = c("Male", "Female")),
      clinic = factor(clinic)
    )
  if ("month_sched" %in% names(d)) {
    d <- d |> mutate(visit_f = factor(month_sched,
                                      levels = c(0, 1, 3, 6, 9, 12)))
  }
  d
}

load_clear <- function() {
  read_csv(here("data", "clear_long.csv"), show_col_types = FALSE) |>
    mutate(arm      = factor(arm, levels = c("Placebo", "Active")),
           sex      = factor(sex, levels = c("Male", "Female")),
           week_f   = factor(week),
           severity_f = factor(severity, levels = 0:3, ordered = TRUE,
                               labels = c("None", "Mild", "Moderate",
                                          "Severe")))
}

load_diary <- function() {
  read_csv(here("data", "diary_long.csv"), show_col_types = FALSE) |>
    mutate(arm = factor(arm, levels = c("Placebo", "Active")),
           sex = factor(sex, levels = c("Male", "Female")),
           period_f = factor(period))
}

load_kidney <- function(which = c("long", "events")) {
  which <- match.arg(which)
  read_csv(here("data", paste0("kidney_", which, ".csv")),
           show_col_types = FALSE) |>
    mutate(sex = factor(sex, levels = c("Male", "Female")))
}

# ---- structure checks ------------------------------------------------------
# check_long(): a quick audit of a long-format dataset. It does not change
# the data; it reports problems so the analyst can decide what to do.

check_long <- function(data, id, time, outcome = NULL) {
  id <- rlang::ensym(id); time <- rlang::ensym(time)
  d  <- data
  n_rows <- nrow(d)
  n_ids  <- dplyr::n_distinct(dplyr::pull(d, !!id))
  dups   <- d |> count(!!id, !!time) |> filter(n > 1)
  per_id <- d |> count(!!id)
  out <- list(
    rows = n_rows,
    subjects = n_ids,
    obs_per_subject = summary(per_id$n),
    duplicated_id_time = nrow(dups),
    missing_time = sum(is.na(pull(d, !!time)))
  )
  if (!is.null(outcome)) {
    out$missing_outcome <- sum(is.na(d[[outcome]]))
  }
  class(out) <- c("long_check", "list")
  out
}

print.long_check <- function(x, ...) {
  cat("Rows:                       ", x$rows, "\n")
  cat("Subjects:                   ", x$subjects, "\n")
  cat("Observations per subject:   ",
      paste(names(x$obs_per_subject), round(x$obs_per_subject, 2),
            collapse = "  "), "\n")
  cat("Duplicated subject-time rows:", x$duplicated_id_time, "\n")
  cat("Missing time values:        ", x$missing_time, "\n")
  if (!is.null(x$missing_outcome))
    cat("Missing outcome values:     ", x$missing_outcome, "\n")
  invisible(x)
}

# Observed correlation of the outcome between scheduled visits (wide data)
visit_correlation <- function(data, id, visit, outcome) {
  data |>
    select({{ id }}, {{ visit }}, {{ outcome }}) |>
    arrange({{ visit }}) |>
    pivot_wider(names_from = {{ visit }}, values_from = {{ outcome }}) |>
    select(-{{ id }}) |>
    cor(use = "pairwise.complete.obs")
}

# Participant flow: how many are seen at each scheduled visit
visit_flow <- function(data, id, visit, group = NULL) {
  data |>
    group_by(across(c({{ group }}, {{ visit }}))) |>
    summarise(n_seen = n_distinct({{ id }}), .groups = "drop")
}
