# Figure Integration Analysis

Based on the available figures in `c:\Users\frevdael\Documents\01_DATA\susana\output\`, here is a recommended strategy for integrating them into the `SOM_dune_development_draft.tex` manuscript versus the Appendix/Supplemental Materials.

The strategy focuses on supporting the core narrative of the paper: SOM acts as a **"double-edged sword"**, enhancing vegetation parameters (germination, height, leaf count) and eco-hydrology ("sponge effect"), but suppressing aeolian sand supply via cohesion. We want to prioritize the clearest, most synthesized visual evidence for the main text.

---

## 1. Main Manuscript Figures
*These figures provide direct visual evidence for the key findings reported in the Results and support the Discussion's primary arguments.*

### A. Substrate Characterization
- **Current Table 1:** Currently, you summarize substrate TOC, TN, and $D_{50}$ in Table 1.
- **Figures:** If a visual representation of the gradient is preferred over the table, you could use a multi-panel figure combining:
  - `output/substrate/TOC.PNG`
  - `output/substrate/TN.PNG`
  - `output/substrate/d(0.5).PNG`
- *Recommendation:* Keep Table 1 for exact numbers, but potentially include `TOC.PNG` as a small panel in the Experimental Setup figure to visually establish the SOM gradient across the 10 mixtures.

### B. Germination and Initial Establishment (Germination Experiment)
- `output/germination/percentage of germinated seeds.png`
- `output/germination/length of seedlings.png`
- *Relevance:* These directly support Section 3.1.2, clearly showing that germination success and initial seedling length are significantly lower on reference beach sand (R) compared to the alternative substrates. A classic 2-panel figure.

### C. Aboveground Growth Trajectories (Common-Garden)
- `output/common_garden/aboveground_biomass/leaves_number_loess.png`
- `output/common_garden/aboveground_biomass/marram_height_loess.png`
- *Relevance:* These smoothed (LOESS) trajectories are perfect for the main text (Section 3.1.1). They show the clean divergence over time: pure D producing the most leaves, and moderate mixtures (D50, D37.5) producing the optimal height compared to the poor performance of R. The smoothed curves are easier to read in a main paper than dense boxplots.

### D. Belowground Allocation (Common-Garden)
- `output/common_garden/belowground_biomass/marram grass biomass allocation on alternative substrates.PNG`
- *Relevance:* While aboveground metrics show one side of the story, root biomass allocation shows how the plant invests resources based on substrate nutrients. This is very relevant to eco-geomorphology.

### E. Microclimate & The "Sponge Effect" Proof
- `output/germination/soil moisture boxplots.png` (or `Soil moisture Vol.png`)
- *Relevance:* This provides essential empirical evidence for the SOM "sponge effect" (dynamic field capacity) mentioned in the Methods and Discussion. Showing that soil moisture is consistently higher in D, D75, D50 than R validates the Saxton & Rawls pedotransfer hypothesis applied in the Living Dunes model.

---

## 2. Appendix / Supplemental Material Figures
*These figures provide necessary detailed statistical backing, raw data distributions, or secondary analyses that would otherwise clutter the main narrative flow.*

### A. Detailed Aboveground Distributions
- `output/common_garden/aboveground_biomass/marram_height_all.png` (or `Marram_height.png`)
- `output/common_garden/aboveground_biomass/number_of_leaves.png`
- *Why Appendix:* Time-series boxplots contain too much visual 'noise' for the main text, but they serve as rigorous proof of the ANOVA/Tukey HSD statistics (showing medians, quartiles, and outliers) that back the LOESS curves in the main text.

### B. Plant Correlations
- `output/common_garden/belowground_biomass/marram height vs shoot biomass.PNG`
- `output/common_garden/root_shoot_ratio/root-shoot ratio.PNG`
- *Why Appendix:* Correlation checks (e.g., verifying that max height correlates with root density, or that root-shoot ratios shift) are valuable context for reviewers but secondary to the main treatment effect (Substrate) on absolute performance.

### C. Detailed Microclimate Data
- `output/germination/climate variables.png`
- `output/germination/air temperature.png`
- `output/germination/soil temperature.png`
- `output/germination/surface temperature.png`
- `output/germination/moisture ifo temperature.png`
- `output/germination/temperature 4am vs 4 pm substrate.png`
- *Why Appendix:* While soil moisture averages are critical for the main text, the full suite of high-resolution temperature and climate variable logs validates the uniformity of the experimental conditions at Zeebrugge.

### D. Detailed Substrate Variance
- `output/substrate/boxplots substrates.PNG`
- `output/substrate/boxplots in-out.PNG`
- *Why Appendix:* If reviewers ask to see the variance among the 5 substrate samples per type tested, these boxplots show spread around the means reported in Table 1.

---

## Summary of Next Steps
1. **Combine for Main Text:** I suggest assembling 3 robust multi-panel figures for the experimental section of the paper:
   - **Fig 1.** Substrate properties (TOC gradient) & Microclimate validation (soil moisture boxplots highlighting the sponge effect).
   - **Fig 2.** Germination performance (% germinated and seedling length).
   - **Fig 3.** Aboveground growth over time (LOESS curves for height and leaf count).
2. **Deposit to Appendix:** Group the high-resolution boxplots (height, leaves, climate logs) into a supplementary PDF file referencing the main text ANOVAs.
