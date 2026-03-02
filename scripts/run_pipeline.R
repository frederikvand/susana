# =============================================================================
# run_pipeline.R
# Master pipeline: sources all analysis scripts in sequence.
# Run from the project root:
#   cd d:/projects/susana
#   Rscript scripts/run_pipeline.R
# =============================================================================

cat("\n##############################################################\n")
cat("  SUSANA Analysis Pipeline\n")
cat("  Running all scripts in sequence ...\n")
cat("##############################################################\n\n")

t0 <- Sys.time()

scripts <- c(
  "scripts/00_data_prep.R",
  "scripts/01_substrate_analysis.R",
  "scripts/02_aboveground_analysis.R",
  "scripts/03_belowground_analysis.R",
  "scripts/04_germination_analysis.R",
  "scripts/05_derive_som_glm.R",
  "scripts/06_copy_figures_to_manuscript.R"
)

for (s in scripts) {
  cat("\n>>> Sourcing:", s, "\n")
  tryCatch(
    source(s, local = new.env(parent = globalenv())),
    error = function(e) {
      cat("\n!!! ERROR in", s, ":\n", conditionMessage(e), "\n")
      stop(paste("Pipeline halted at", s), call. = FALSE)
    }
  )
}

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)

cat("\n##############################################################\n")
cat("  Pipeline complete in", elapsed, "minutes\n")
cat("  All outputs written to output/ and evaluation/manuscript/\n")
cat("##############################################################\n")
