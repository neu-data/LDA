# =============================================================================
# scripts/00_make_data.R
# Regenerates every simulated course dataset from data-raw/ (fixed seeds) and
# records the R session used. Run from the project root:
#   Rscript scripts/00_make_data.R
# The datasets are already supplied in data/, so running this is optional;
# it demonstrates that they are fully reproducible.
# =============================================================================

library(here)

sims <- c("01_simulate_walk.R", "02_simulate_clear.R",
          "03_simulate_diary.R", "04_simulate_kidney.R",
          "05_simulate_mood.R")

for (s in sims) {
  message("Running data-raw/", s)
  # each script runs in its own environment so objects do not leak between
  source(here("data-raw", s), local = new.env())
}

# Record the computational environment used to create the data
dir.create(here("outputs"), showWarnings = FALSE)
writeLines(capture.output(sessionInfo()),
           here("outputs", "session_info_data.txt"))
message("Done. Data written to data/; session info in outputs/.")
