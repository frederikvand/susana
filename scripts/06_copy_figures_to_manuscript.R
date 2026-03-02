# =============================================================================
# 06_copy_figures_to_manuscript.R
# Copies all publication figures and tables from output/ to manuscript/figures/
# and updates \includegraphics paths in the .tex files
# =============================================================================

suppressPackageStartupMessages(library(here))

cat("\n======================================================\n")
cat("  06_copy_figures_to_manuscript.R\n")
cat("======================================================\n\n")

fig_dir <- here("evaluation", "manuscript", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Define the figure mapping: source -> manuscript filename ---
figure_map <- list(
  # Substrate figures
  list(src  = here("output", "substrate", "substrate_toc_d50.pdf"),
       dest = "substrate_toc_d50.pdf"),
  list(src  = here("output", "substrate", "substrate_in_out_comparison.pdf"),
       dest = "substrate_in_out_comparison.pdf"),
  list(src  = here("output", "substrate", "substrate_tn.pdf"),
       dest = "substrate_tn.pdf"),
  list(src  = here("output", "substrate", "soil_moisture_boxplots.pdf"),
       dest = "soil_moisture_boxplots.pdf"),

  # Aboveground figures
  list(src  = here("output", "common_garden", "aboveground_biomass", "aboveground_growth_main.pdf"),
       dest = "aboveground_growth_main.pdf"),
  list(src  = here("output", "common_garden", "aboveground_biomass", "marram_height_all.pdf"),
       dest = "marram_height_all.pdf"),
  list(src  = here("output", "common_garden", "aboveground_biomass", "number_of_leaves.pdf"),
       dest = "number_of_leaves.pdf"),

  # Belowground figures
  list(src  = here("output", "common_garden", "belowground_biomass", "belowground_main.pdf"),
       dest = "belowground_main.pdf"),
  list(src  = here("output", "common_garden", "belowground_biomass", "root_depth_profile.pdf"),
       dest = "root_depth_profile.pdf"),
  list(src  = here("output", "common_garden", "belowground_biomass", "root_vs_shoot_scatter.pdf"),
       dest = "root_vs_shoot_scatter.pdf"),
  list(src  = here("output", "common_garden", "belowground_biomass", "height_vs_shoot_biomass.pdf"),
       dest = "height_vs_shoot_biomass.pdf"),

  # Germination figures
  list(src  = here("output", "germination", "germination_main.pdf"),
       dest = "germination_main.pdf"),
  list(src  = here("output", "germination", "soil_temperature_timeseries.pdf"),
       dest = "soil_temperature_timeseries.pdf"),
  list(src  = here("output", "germination", "soil_moisture_timeseries.pdf"),
       dest = "soil_moisture_timeseries.pdf"),
  list(src  = here("output", "germination", "microclimate_boxplots.pdf"),
       dest = "microclimate_boxplots.pdf"),
  list(src  = here("output", "germination", "diurnal_temperature.pdf"),
       dest = "diurnal_temperature.pdf"),

  # SOM response curves
  list(src  = here("output", "som_response", "som_response_R_derivation.pdf"),
       dest = "som_response_curves.pdf")
)

# --- Copy figures ---
cat("--- Copying figures ---\n")
copied <- 0
missing <- 0
for (item in figure_map) {
  if (file.exists(item$src)) {
    file.copy(item$src, file.path(fig_dir, item$dest), overwrite = TRUE)
    cat("  OK:", item$dest, "\n")
    copied <- copied + 1
  } else {
    # Try PNG fallback
    png_src <- sub("\\.pdf$", ".png", item$src)
    png_dest <- sub("\\.pdf$", ".png", item$dest)
    if (file.exists(png_src)) {
      file.copy(png_src, file.path(fig_dir, png_dest), overwrite = TRUE)
      cat("  OK (PNG):", png_dest, "\n")
      copied <- copied + 1
    } else {
      cat("  MISSING:", item$src, "\n")
      missing <- missing + 1
    }
  }
}
cat(sprintf("\n  Copied: %d  Missing: %d\n", copied, missing))

# --- Update LaTeX paths ---
cat("\n--- Updating LaTeX \\includegraphics paths ---\n")

update_tex_paths <- function(tex_file) {
  lines <- readLines(tex_file, warn = FALSE)
  original <- lines

  # Map of old path patterns -> new figure-relative paths
  replacements <- list(
    # Substrate
    c("../../output/substrate/TOC.PNG",
      "figures/substrate_toc_d50.pdf"),
    c("../../output/substrate/substrate_toc_d50.pdf",
      "figures/substrate_toc_d50.pdf"),
    c("../../output/substrate/substrate_toc_d50.png",
      "figures/substrate_toc_d50.pdf"),
    c("../../output/substrate/boxplots substrates.PNG",
      "figures/substrate_in_out_comparison.pdf"),
    c("../../output/substrate/boxplots in-out.PNG",
      "figures/substrate_in_out_comparison.pdf"),
    c("../../output/substrate/substrate_in_out_comparison.pdf",
      "figures/substrate_in_out_comparison.pdf"),
    c("../../output/substrate/substrate_tn.pdf",
      "figures/substrate_tn.pdf"),

    # Soil moisture
    c("../../output/germination/soil moisture boxplots.png",
      "figures/soil_moisture_boxplots.pdf"),
    c("../../output/germination/soil_moisture_boxplots.pdf",
      "figures/soil_moisture_boxplots.pdf"),
    c("../../output/substrate/soil_moisture_boxplots.pdf",
      "figures/soil_moisture_boxplots.pdf"),

    # Aboveground
    c("../../output/common_garden/aboveground_biomass/marram_height_loess.png",
      "figures/aboveground_growth_main.pdf"),
    c("../../output/common_garden/aboveground_biomass/marram_height_loess.pdf",
      "figures/aboveground_growth_main.pdf"),
    c("../../output/common_garden/aboveground_biomass/leaves_number_loess.png",
      "figures/aboveground_growth_main.pdf"),
    c("../../output/common_garden/aboveground_biomass/aboveground_growth_main.pdf",
      "figures/aboveground_growth_main.pdf"),
    c("../../output/common_garden/aboveground_biomass/marram_height_all.png",
      "figures/marram_height_all.pdf"),
    c("../../output/common_garden/aboveground_biomass/marram_height_all.pdf",
      "figures/marram_height_all.pdf"),
    c("../../output/common_garden/aboveground_biomass/number_of_leaves.png",
      "figures/number_of_leaves.pdf"),
    c("../../output/common_garden/aboveground_biomass/number_of_leaves.pdf",
      "figures/number_of_leaves.pdf"),

    # Belowground
    c("../../output/common_garden/belowground_biomass/marram height vs shoot biomass.PNG",
      "figures/height_vs_shoot_biomass.pdf"),
    c("../../output/common_garden/belowground_biomass/height_vs_shoot_biomass.pdf",
      "figures/height_vs_shoot_biomass.pdf"),
    c("../../output/common_garden/root_shoot_ratio/root-shoot ratio.PNG",
      "figures/belowground_main.pdf"),
    c("../../output/common_garden/belowground_biomass/belowground_main.pdf",
      "figures/belowground_main.pdf"),

    # Germination
    c("../../output/germination/percentage of germinated seeds.png",
      "figures/germination_main.pdf"),
    c("../../output/germination/length of seedlings.png",
      "figures/germination_main.pdf"),
    c("../../output/germination/germination_main.pdf",
      "figures/germination_main.pdf"),

    # Microclimate
    c("../../output/germination/climate variables.png",
      "figures/microclimate_boxplots.pdf"),
    c("../../output/germination/air temperature.png",
      "figures/soil_temperature_timeseries.pdf"),
    c("../../output/germination/soil temperature.png",
      "figures/soil_temperature_timeseries.pdf"),
    c("../../output/germination/surface temperature.png",
      "figures/soil_temperature_timeseries.pdf"),
    c("../../output/germination/temperature 4am vs 4 pm substrate.png",
      "figures/diurnal_temperature.pdf"),
    c("../../output/germination/moisture ifo temperature.png",
      "figures/soil_moisture_timeseries.pdf"),

    # SOM response
    c("../../output/som_response_R_derivation.png",
      "figures/som_response_curves.pdf"),
    c("../../output/som_response/som_response_R_derivation.png",
      "figures/som_response_curves.pdf"),
    c("../../output/som_response/som_response_R_derivation.pdf",
      "figures/som_response_curves.pdf")
  )

  for (rep in replacements) {
    # Escape special regex characters in the old path
    old_pat <- gsub("([.()\\[\\]])", "\\\\\\1", rep[1])
    lines <- gsub(old_pat, rep[2], lines, fixed = FALSE)
  }

  n_changed <- sum(lines != original)
  if (n_changed > 0) {
    writeLines(lines, tex_file)
    cat("  Updated", tex_file, ":", n_changed, "lines changed\n")
  } else {
    cat("  No changes needed in", tex_file, "\n")
  }
}

update_tex_paths(here("evaluation", "manuscript", "SOM_dune_development_draft.tex"))
update_tex_paths(here("evaluation", "manuscript", "SOM_dune_development_appendix.tex"))


cat("\n======================================================\n")
cat("  06_copy_figures_to_manuscript.R COMPLETE\n")
cat("======================================================\n")
