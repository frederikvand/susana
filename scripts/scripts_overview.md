# SUSANA Analysis Scripts Overview

Pipeline for the SUSANA manuscript: *Optimizing Alternative Sediments for Dune Construction*.
All scripts are in `d:/projects/susana/scripts/`.

## Execution Order

Run scripts sequentially; each depends on output from the previous steps.

| # | Script | Purpose | Key Statistics | Outputs |
|---|--------|---------|---------------|---------|
| — | `shared_theme.R` | Centralised ggplot2 theme, viridis substrate palette, `save_pub_plot()` helper | — | Sourced by all scripts |
| 0 | `00_data_prep.R` | Reads raw Excel → standardised `.rds` files in `output/clean_data/` | — | 11 `.rds` files |
| 1 | `01_substrate_analysis.R` | Substrate characterisation: TOC, TN, D50 | Wilcoxon rank-sum (Alt vs Ctrl, In vs Out) | `substrate_toc_d50`, `substrate_in_out_comparison`, `substrate_tn`, `soil_moisture_boxplots` (PDF+PNG); `substrate_statistics.tex`, `substrate_summary.tex` |
| 2 | `02_aboveground_analysis.R` | Common-garden height & leaf count over time | Gaussian LMM `~ substrate * poly(time, 2) + (1\|plant_id)` for height; Poisson GLMM (glmmTMB, log link) for leaf count; EMMs with FDR-adjusted CLD; `ggeffects::ggpredict()` for model trajectories | `aboveground_growth_main`, `marram_height_loess`, `leaves_number_loess`, `marram_height_all`, `number_of_leaves` (PDF+PNG); `aboveground_lmm_results.tex` |
| 3 | `03_belowground_analysis.R` | Root biomass, root:shoot ratio, depth profiles | Gamma GLM (log link) primary + KW fallback; Pearson correlation height vs shoot | `belowground_main`, `root_depth_profile`, `root_vs_shoot_scatter`, `height_vs_shoot_biomass` (PDF+PNG); `belowground_biomass_summary.tex`, `belowground_statistics.tex` |
| 4 | `04_germination_analysis.R` | Germination rate, seedling growth, microclimate | Quasibinomial GLM for germination (overdispersion corrected); Poisson GLMM (glmmTMB) `~ substrate * poly(time, 2) + (1\|pot_id)` for seedling count; KW for moisture/temp; `ggeffects::ggpredict()` for trajectories | `germination_main`, `soil_temperature_timeseries`, `soil_moisture_timeseries`, `microclimate_boxplots`, `diurnal_temperature` (PDF+PNG); `germination_statistics.tex`, `microclimate_summary.tex` |
| 5 | `05_derive_som_glm.R` | SOM response curves for Living Dunes parameterisation | Weighted NLS (Levenberg–Marquardt) on treatment means; Binomial GLMM (glmmTMB) for mortality | `som_response_R_derivation` (PDF+PNG); `growth_model_coefficients.tex`, `mortality_model_coefficients.tex`, `som_derived_parameters.tex`, `model_fit_statistics.tex` |
| 6 | `06_copy_figures_to_manuscript.R` | Copies PDFs to `evaluation/manuscript/figures/`; updates `\includegraphics` paths in `.tex` | — | 17 figures copied; 2 `.tex` files patched |

## Statistical Methods Summary

| Method | Package | Used in |
|--------|---------|---------|
| Gaussian LMM (quadratic time) | `lme4::lmer` | 02 (height) |
| Poisson GLMM (log link) | `glmmTMB` | 02 (leaf count), 04 (seedling count) |
| Quasibinomial GLM | base R `glm` | 04 (germination) |
| Gamma GLM (log link) | base R `glm` | 03 (biomass) |
| Estimated marginal means + CLD | `emmeans`, `multcompView` | 02, 03, 04 |
| Model-predicted trajectories | `ggeffects::ggpredict` | 02, 04 |
| Kruskal–Wallis (parallel reference) | base R `kruskal.test` | 03, 04 |
| Wilcoxon rank-sum | base R `wilcox.test` | 01 |
| Weighted NLS (Levenberg–Marquardt) | `minpack.lm::nlsLM` | 05 (SOM response) |
| Binomial GLMM | `glmmTMB` | 05 (SOM mortality) |

## Output Directories

| Directory | Content |
|-----------|---------|
| `output/clean_data/` | Standardised `.rds` data files |
| `output/substrate/` | Substrate figures + tables |
| `output/common_garden/aboveground_biomass/` | Aboveground figures + table |
| `output/common_garden/belowground_biomass/` | Belowground figures + tables |
| `output/germination/` | Germination & microclimate figures + tables |
| `output/som_response/` | SOM dose–response figures + tables |
| `evaluation/manuscript/figures/` | Publication-ready copies |
| `evaluation/manuscript/tables/` | LaTeX table fragments (\input) |
