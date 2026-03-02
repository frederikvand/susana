# SUSANA — Optimizing Alternative Sediments for Dune Construction

This repository contains the analysis pipeline and manuscript for the SUSANA project, maintained by the **Terrestrial Ecology Unit (TEREC)** at Ghent University.

## Manuscript

The primary output is a scientific manuscript: `evaluation/manuscript/SOM_dune_development_draft.tex`.

It synthesizes findings from two core work package deliverables (WP 1.1 and WP 1.4) concerning boundary conditions for dune vegetation development and the fundamental/realised niche of *Ammophila arenaria* (marram grass) on alternative sediments. Empirical findings on biological dune development are combined with geophysical storm erosion modeling (KU Leuven) and **dune development modelling using the Living Dunes model**, simulating dune restoration across 10 different substrate types.

## Analysis Pipeline

All analysis is orchestrated by numbered R scripts in `scripts/`. Run the full pipeline from the project root:

```bash
cd d:/projects/susana
Rscript scripts/run_pipeline.R
```

Or run individual scripts:

```bash
Rscript scripts/02_aboveground_analysis.R
```

### Script Sequence

| # | Script | Purpose |
|---|--------|---------|
| — | `shared_theme.R` | Centralised ggplot2 theme (`theme_susana()`), viridis substrate palette, `save_pub_plot()` helper |
| 00 | `00_data_prep.R` | Reads raw Excel data → standardised `.rds` files in `output/clean_data/` |
| 01 | `01_substrate_analysis.R` | Substrate characterisation: TOC, TN, D50 (Wilcoxon rank-sum tests) |
| 02 | `02_aboveground_analysis.R` | Common-garden height & leaf count: Gaussian LMM (height), Poisson GLMM (leaf count), EMMs with CLD |
| 03 | `03_belowground_analysis.R` | Root biomass, root:shoot ratio: Gamma GLM + Kruskal-Wallis, Pearson correlation |
| 04 | `04_germination_analysis.R` | Germination & seedling growth: quasibinomial GLM (germination), Poisson GLMM (seedlings), microclimate KW |
| 05 | `05_derive_som_glm.R` | SOM dose–response curves: weighted NLS (Levenberg–Marquardt) for growth, binomial GLMM for mortality |
| 06 | `06_copy_figures_to_manuscript.R` | Copies 17 publication-ready PDFs to `evaluation/manuscript/figures/` |

For a full statistical methods overview, see `scripts/scripts_overview.md`.

### Reproducibility

- All scripts use `here::here()` for portable paths — no hardcoded absolute paths
- All analysis scripts set `set.seed(2024)` for reproducibility
- `run_pipeline.R` runs all scripts in sequence with error handling and timing
- The pipeline takes approximately 1 minute on a standard machine
- Outputs are deterministic: re-running produces identical figures and tables

### Dependencies

R packages (install with `install.packages()`):

```
ggplot2, dplyr, tidyr, readxl, lme4, glmmTMB, emmeans, multcomp,
multcompView, patchwork, viridis, scales, FSA, minpack.lm, ggeffects, here
```

LaTeX: TeX Live 2024 (or equivalent) with `latexmk`.

## Directory Structure

```
susana/
├── data/                          # Raw experimental data (Excel, CSV)
│   ├── geomorphology/             # Point cloud scan data
│   ├── photos/                    # Photographic records
│   ├── raw/                       # Unprocessed data dumps
│   ├── substrate/                 # Soil analysis data (TOC, TN, grain size)
│   └── vegetation_development/    # Plant growth measurements
│
├── scripts/                       # Analysis pipeline (R)
│   ├── shared_theme.R             # Centralised plotting theme & helpers
│   ├── 00_data_prep.R             # Data ingestion & cleaning
│   ├── 01_substrate_analysis.R    # Substrate characterisation
│   ├── 02_aboveground_analysis.R  # Height & leaf count LMMs
│   ├── 03_belowground_analysis.R  # Root biomass Gamma GLMs
│   ├── 04_germination_analysis.R  # Germination & microclimate
│   ├── 05_derive_som_glm.R        # SOM response curves for Living Dunes
│   ├── 06_copy_figures_to_manuscript.R  # Copy figures to manuscript dir
│   ├── run_pipeline.R             # Master pipeline runner
│   ├── scripts_overview.md        # Detailed statistical methods
│   ├── original/                  # Original student scripts (archived)
│   └── z_previous_version/        # Intermediate script versions (archived)
│
├── output/                        # Generated outputs (figures, tables, data)
│   ├── clean_data/                # Standardised .rds files from script 00
│   ├── substrate/                 # Substrate figures (PDF + PNG)
│   ├── common_garden/
│   │   ├── aboveground_biomass/   # Height & leaf count figures
│   │   └── belowground_biomass/   # Root biomass figures
│   ├── germination/               # Germination & microclimate figures
│   ├── som_response/              # SOM dose–response figures
│   └── z_old/                     # Archived debug logs & old figures
│
├── evaluation/
│   ├── manuscript/
│   │   ├── SOM_dune_development_draft.tex     # Main manuscript
│   │   ├── SOM_dune_development_appendix.tex  # Supplementary appendix
│   │   ├── references.bib                     # Bibliography
│   │   ├── figures/               # Publication-ready PDFs (copied by script 06)
│   │   └── tables/                # LaTeX table fragments (\input)
│   ├── analysis/                  # Additional analysis projects
│   ├── deliverables/              # WP reports
│   └── presentation/              # Slides, posters
│
├── docs/                          # Documentation & literature
│   ├── dune_modelling/            # Dune model documentation
│   ├── literature/                # Reference papers
│   ├── proposal/                  # Grant proposals
│   └── protocol/                  # Experimental protocols
│
└── z_old/                         # Top-level archived files
```

## Compiling the Manuscript

```bash
cd evaluation/manuscript
latexmk -pdf SOM_dune_development_draft.tex
```

The manuscript uses `\graphicspath{{figures/}}` so all figure paths are relative. Tables are included via `\input{tables/filename.tex}` and are regenerated by the pipeline.
