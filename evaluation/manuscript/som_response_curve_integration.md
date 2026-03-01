# Deriving and Integrating SOM Response Curves for *Ammophila arenaria*

This report outlines the methodology for deriving Soil Organic Matter (SOM) response curves from the SUSANA common garden experimental data and the plan for integrating these environmental responses into the Living Dunes model.

## Part 1: Deriving SOM Response Curves from Experimental Data

The common garden experiment provides data on *Ammophila arenaria* growth across a gradient of alternative sediments (combinations of dune sand, beach sand, and dredged material). To translate this data into model parameters, we must isolate the specific effects of SOM on plant vitality and demography.

According to the "double-edged sword" concept and the analysis of the experimental data, SOM influences *Ammophila arenaria* through three primary mechanisms:
1.  **Nutrient Boost (Abiotic Positive):** Initial surge in growth due to nitrogen availability.
2.  **Biotic Penalty (Pathogens):** Increased mortality and reduced vigor at higher SOM levels due to soil-borne pathogens (fungi, nematodes).
3.  **Abiotic Penalty (Sponge Effect):** Increased soil moisture retention at high SOM fractions leads to hypoxia.

### 1. Data Preparation and Independent Variables

The primary independent variable for these curves is the **Soil Organic Matter (SOM) mass fraction** ($f_{som}$). This can be derived directly from the Total Organic Carbon (TOC) analysis in the dataset, often using a standard conversion factor (e.g., SOM = TOC $\times$ 1.724) or using the TOC values directly if the model normalizes them.

Secondary environmental variables available in the dataset that act as intermediaries or modulators include:
*   **$d_{50}$ (Median Grain Size):** Influences drainage and cohesion.
*   **Nitrogen (N) Content:** Drives the initial nutrient boost.
*   **Volumetric Soil Moisture ($\theta$):** The primary driver of the abiotic penalty (hypoxia).

### 2. Deriving the Response Functions

We need to establish empirical relationships between SOM (or its consequences) and the core demographic rates: **Growth Rate ($G$)** and **Mortality Rate ($M$)**.

#### A. The "Nutrient Boost" Function (Growth Multiplier)
The data shows that moderate SOM (e.g., the 50/50 D/W mixture) yields the highest aboveground biomass and leaf length. This acts as a multiplier on the base growth rate.
*   **Data Source:** Leaf length ($cm$) and Number of leaves over time.
*   **Function Shape:** A parabolic or log-normal curve peaking at an optimal SOM concentration, or a saturation curve based on Nitrogen if we separate the variables. Given the negative effects at very high SOM, a single combined function is often simplest.
*   **Derivation:**
    1.  Normalize growth metrics (e.g., total leaf length increase over the first season) relative to the reference Dune Sand (D).
    2.  Plot normalized growth against initial soil Nitrogen content or SOM fraction.
    3.  Fit a non-linear regression function. A suitable form might be a modified Gaussian or a generic beta function bounded between 0 and maximum observed SOM.
    4.  $Growth\_Multiplier(SOM) = \dots$

#### B. The "Biotic Penalty" Function (Vitality Reduction/Mortality Increase)
High SOM environments promote pathogenic activity, leading to reduced vigor and increased mortality, a widely documented phenomenon for *A. arenaria*.
*   **Data Source:** Survival rates, below-ground biomass allocation (root/shoot ratio), and late-season vigor decline.
*   **Function Shape:** A logistic or threshold function turning on sharply above a critical SOM concentration (often cited as ~1.0% by mass).
*   **Derivation:**
    1.  Identify a vitality index (e.g., final survival probability or a scaled combined health metric).
    2.  Plot survival/vitality against SOM.
    3.  Fit a reverse logistic function: $Penalty(SOM) = \frac{1}{1 + e^{k(SOM - SOM_{crit})}}$
    4.  This penalty can directly increase the basal mortality rate or decrease the maximum potential growth rate in the model.

#### C. The "Sponge Effect" Function (Abiotic Mortality)
Fine-grained, organic-rich sediments retain excessive water. *A. arenaria* is highly sensitive to waterlogging (hypoxia).
*   **Data Source:** Volumetric soil moisture data correlated with substrate type and plant survival.
*   **Function Shape:** A threshold mortality function that activates when soil moisture exceeds field capacity or a critical anoxic threshold.
*   **Derivation:**
    1.  Utilize established pedotransfer functions (e.g., Saxton and Rawls) to link SOM and $d_{50}$ to soil moisture retention capabilities at field capacity ($\theta_{fc}$).
    2.  Correlate continuous periods of soil moisture exceeding a critical threshold ($\theta_{crit}$, potentially defined around 4-6% for optimal, but much higher for hypoxia) with elevated mortality in the experimental setup or literature.
    3.  Define an abiotic survival modifier: $S_{abiotic}(\theta) = \dots$ (declining rapidly as $\theta$ exceeds $\theta_{crit}$).

### 3. Statistical Fitting Procedure

1.  **Extract Matrices:** Form feature matrices ($X$) comprising SOM, N, $d_{50}$, etc., and target vectors ($y$) for growth (e.g., $\Delta$ height) and mortality from the `data/` directory.
2.  **Model Selection:** Use `scipy.optimize.curve_fit` or general linear mixed models (GLMMs) in R/Python (`statsmodels`) to fit the proposed functional shapes (logistic, parabolic).
3.  **Evaluate Fit:** Use AIC/BIC and $R^2$ to select the most parsimonious curves that capture the observed dynamics without overfitting the noise.

### 4. Implementation and Results

We performed standard General Linear Modeling (GLM) and Non-Linear Least Squares (NLS) regression in **R** (using `minpack.lm` and base stats) to derive strict ecological parameters for *Ammophila arenaria*, decoupling the SOM response across four independent demographic vectors:

*   **Mortality Rate Modifier (quasibinomial GLM):**
    We fitted a strict quasibinomial (pseudobinomial) model linking SOM to survival probability to account for overdispersion in the experimental survival ratios.
    *   **Fit:** Significant intercept ($p = 0.003$), yielding an `inverse_logistic` equivalent defined by a steepness of $-0.044$ and a midpoint threshold at very high theoretical limits (due to the relatively slow linear decline across the tested gradient).
*   **Height Growth Multiplier (Gaussian NLS):**
    Height exhibited the most classic parabolic "Nutrient Boost" before hitting the "Sponge Effect".
    *   **Fit:** $\mu = 1.36\%$ SOM, $\sigma = 0.89\%$, and a multiplier peak of $+30\%$.
*   **Density Growth Multiplier (Gaussian NLS):**
    Shoot density exhibited a much narrower optimal window, preferring lower SOM mixtures.
    *   **Fit:** Sharp peak at $\mu = 0.62\%$ SOM, $\sigma = 0.067\%$.
*   **Root Growth Multiplier (Gaussian NLS):**
    Belowground biomass (root growth) showed extreme hyper-sensitivity to specific intermediate SOM states.
    *   **Fit:** Peak at $\mu = 0.83\%$ SOM, $\sigma = 0.067\%$.

These statistically robust curves are visualized below:

![SOM Response R Derivations](file:///d:/projects/susana/output/som_response_R_derivation.png)

These distinctive parameters must be individually injected into `ammophila_arenaria.json` (e.g., nesting the height parameters under `"height_growth"`, density under `"density_growth"`, and the quasibinomial GLM output under `"mortality"`). This decoupling ensures the complex "Double-Edged Sword" manifests correctly down to the morphological level within Living Dunes.

---

## Part 2: Integrating SOM Response into Living Dunes

Currently, `demography_tools.py` in Living Dunes calculates demographic responses based primarily on hydrodynamics (inundation) and morphodynamics (burial/erosion). We need to introduce the **edaphic environment** (specifically SOM and its derivatives) as a first-class driver.

### 1. Extending `soil_organic_matter.py`

This module currently tracks spatial SOM distribution over time (production and decay). It must be updated to export localized soil properties that plant demography depends on.

*   **Changes Required:** Ensure the `SoilOrganicMatter` class provides easy public access to spatial arrays of `som_fraction` and derived `soil_moisture` (perhaps coupling with a simplified moisture state if not fully resolved by a groundwater module).
*   **Output:** The state array of SOM concentration needs to be passed to the vegetation update step.

### 2. Modifying `demography_tools.py`

We must integrate the derived mathematical functions into the core demographic calculations. The primary entry point for environmental stressors is `calculate_demographic_response`.

*   **Current State:**
    ```python
    def calculate_demographic_response(
        demography: np.ndarray,
        elevation_change: np.ndarray,
        inundation_depth: np.ndarray,
        salinity: np.ndarray,
        species_params: SpeciesParameters
    ) -> DemographicResponse:
    ```
*   **Proposed State:** Add `som_fraction` (and potentially `soil_moisture`) to the signature.
    ```python
    def calculate_demographic_response(
        demography: np.ndarray,
        elevation_change: np.ndarray,
        inundation_depth: np.ndarray,
        salinity: np.ndarray,
        som_fraction: np.ndarray, # NEW INJECTION
        soil_moisture: np.ndarray, # NEW INJECTION
        species_params: SpeciesParameters
    ) -> DemographicResponse:
    ```

*   **Internal Calculation Updates:**
    Inside `demography_tools.py`, we add new modifier functions equivalent to existing burial/erosion modifiers:
    1.  `calculate_som_growth_modifier(som_fraction, species_params)`: Implements the "Nutrient Boost" parabolic curve. Returns a multiplier array $[0, \infty)$ (typically $0.5$ to $1.5$) applied to `base_growth_rate`.
    2.  `calculate_som_mortality(som_fraction, soil_moisture, species_params)`: Combines the "Biotic Penalty" (logistic depending on SOM) and "Sponge Effect" (threshold depending on moisture). Returns an additive mortality rate array applied to the total mortality.

### 3. Updating the Configuration Registry (`ammophila_arenaria.json`)

The mathematical derivations from Part 1 yield specific parameter values (optimal SOM, critical thresholds, logistic steepness). These must be defined in the species configuration schema.

*   **Changes to `living_dunes/config/default/atlantic_europe/species/ammophila_arenaria.json`:**
    Add a new sub-dictionary under the `environmental_response` block:
    ```json
    "edaphic_response": {
        "som_growth_optimal": 0.005,      // Output of Part 1 (e.g., 0.5% SOM)
        "som_growth_multiplier_max": 1.4, // Output of Part 1 (40% boost)
        "som_biotic_penalty_k": 15.0,     // Steepness of logistic penalty
        "som_critical_threshold": 0.01,   // 1.0% SOM threshold for pathogens
        "hypoxia_moisture_threshold": 0.15 // Volumetric moisture threshold
    }
    ```
*   **Schema Update:** Ensure `config/schemas/species_schema.json` is updated to accept these new keys.

### 4. Updating Documentation (`demography_tools.md`)

The comprehensive markdown documentation for the module must reflect the new physics.

*   **Changes to `docs/modules/vegetation/technical/demography_tools.md`:**
    1.  **Environmental Interaction Models Section:** Add a new subsection titled **"Edaphic Factors (Soil Organic Matter and Moisture)"**.
    2.  **Mathematical Formulations:** Document the exact equations derived in Part 1 (the parabolic growth multiplier and logistic mortality penalty).
    3.  **Variable Definitions:** Clearly define all new species parameters (`som_growth_optimal`, etc.) introduced in the JSON configuration.
    4.  **Flowchart/Diagram:** If there is a mermaid diagram showing data flow into `calculate_demographic_response`, update it to include the `SoilOrganicMatter` state arrays.
