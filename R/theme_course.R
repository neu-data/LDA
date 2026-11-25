# =============================================================================
# R/theme_course.R
# One restrained visual style for every figure in notes, slides and outputs.
# Colours are from the Okabe-Ito colour-blind-safe palette.
# =============================================================================

course_colours <- c(
  ink     = "#222222",
  grey    = "#6B6B6B",
  light   = "#D9D9D9",
  blue    = "#0072B2",
  orange  = "#D55E00",
  green   = "#009E73",
  purple  = "#8E5A9F"
)

# Arm colours are fixed so a given arm looks the same everywhere
arm_colours <- c(
  "UC" = "#6B6B6B", "Usual care" = "#6B6B6B",
  "SP" = "#0072B2", "Standard programme" = "#0072B2",
  "EP" = "#D55E00", "Enhanced programme" = "#D55E00",
  "Placebo" = "#6B6B6B", "Active" = "#0072B2",
  "Waitlist" = "#6B6B6B", "Programme" = "#0072B2"
)

theme_course <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = "") +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "#ECECEC", linewidth = 0.4),
      axis.line         = element_line(colour = "#9A9A9A", linewidth = 0.4),
      axis.ticks        = element_line(colour = "#9A9A9A", linewidth = 0.4),
      axis.title        = element_text(colour = "#222222"),
      axis.text         = element_text(colour = "#444444"),
      strip.text        = element_text(face = "bold", hjust = 0,
                                       colour = "#222222"),
      legend.position   = "bottom",
      legend.title      = element_text(size = rel(0.9)),
      plot.title        = element_text(face = "bold", size = rel(1.05),
                                       colour = "#222222"),
      plot.title.position = "plot",
      plot.subtitle     = element_text(colour = "#555555", size = rel(0.9)),
      plot.caption      = element_text(colour = "#777777", size = rel(0.75),
                                       hjust = 0)
    )
}

scale_colour_arm <- function(...) {
  scale_colour_manual(values = arm_colours, ...)
}
scale_fill_arm <- function(...) {
  scale_fill_manual(values = arm_colours, ...)
}

theme_set(theme_course())
update_geom_defaults("point", list(colour = "#222222"))
update_geom_defaults("line",  list(colour = "#222222"))
