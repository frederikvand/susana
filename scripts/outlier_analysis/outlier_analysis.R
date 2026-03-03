# =============================================================================
# outlier_analysis.R
# Analyse density growth outlier for SOM response curve
# =============================================================================

suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(here)
})

# Load data
substrate_path <- here("data", "substrate", "calculated_substrate_properties.xlsx")
ag_path <- here("data", "vegetation_development", "common_garden", "aboveground_biomass", "aboveground_biomass_zeebrugge.xlsx")

df_sub <- read_excel(substrate_path)
som_vals <- pull(df_sub, `TOC (%)`)
lookup <- data.frame(
    substrate = c("R", "W", "D", "12.5", "25", "37.5", "50", "62.5", "75", "87.5"),
    som = som_vals
)

cat("\n--- Substrate TOC values ---\n")
print(lookup)

df_ag <- read_excel(ag_path)
df_ag$length <- suppressWarnings(as.numeric(df_ag$length))
df_ag$nr     <- suppressWarnings(as.numeric(df_ag$nr))
df_ag$substrate <- as.character(df_ag$substrate)
df_ag <- df_ag %>%
    left_join(lookup, by = "substrate") %>%
    filter(!is.na(som))

# Get peak leaf count per individual plant
ag_max <- df_ag %>%
    group_by(ID, column, som) %>%
    summarize(max_nr = max(nr, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(max_nr))

ref_n <- mean(ag_max$max_nr[ag_max$som == min(ag_max$som)], na.rm = TRUE)
ag_max$n_norm <- ag_max$max_nr / ref_n

# Treatment means for density
means_density <- ag_max %>%
    group_by(som) %>%
    summarize(mean_n = mean(n_norm), se_n = sd(n_norm) / sqrt(n()), n = n(), .groups = "drop")
cat("\n--- Density treatment means ---\n")
print(as.data.frame(means_density))

# Show individual measurements for the peak SOM
peak_som <- means_density$som[which.max(means_density$mean_n)]
cat("\nPeak at TOC =", peak_som, "\n")
cat("--- Individual data at peak SOM ---\n")
peak_data <- ag_max %>% filter(som == peak_som)
print(as.data.frame(peak_data))

# Show boxplot statistics for each treatment
cat("\n--- Summary stats per treatment ---\n")
stats <- ag_max %>% group_by(som) %>%
    summarize(
        n = n(), Mean = mean(n_norm), Median = median(n_norm), SD = sd(n_norm),
        Min = min(n_norm), Max = max(n_norm),
        Q25 = quantile(n_norm, 0.25), Q75 = quantile(n_norm, 0.75),
        IQR = IQR(n_norm),
        Upper_fence = quantile(n_norm, 0.75) + 1.5 * IQR(n_norm),
        .groups = "drop")
print(as.data.frame(stats), digits = 3)

# Identify statistical outliers (Tukey fences)
cat("\n--- Tukey outliers (>1.5*IQR above Q3) per treatment ---\n")
for (s in sort(unique(ag_max$som))) {
    sub <- ag_max$n_norm[ag_max$som == s]
    q3 <- quantile(sub, 0.75)
    iqr <- IQR(sub)
    fence <- q3 + 1.5 * iqr
    outliers <- sub[sub > fence]
    if (length(outliers) > 0) {
        cat("TOC =", s, ": Outliers:", outliers, "(fence =", fence, ")\n")
        ids <- ag_max$ID[ag_max$som == s & ag_max$n_norm > fence]
        cat("  Plant IDs:", ids, "\n")
    }
}

# Compare NLS fit with and without the outlier
suppressPackageStartupMessages(library(minpack.lm))
lm_control <- nls.lm.control(maxiter = 200)

cat("\n============================================================\n")
cat("  NLS FIT COMPARISON: WITH vs WITHOUT OUTLIER\n")
cat("============================================================\n")

# 1) Fit WITH all data
fit_all <- nlsLM(
    mean_n ~ b0 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
    data    = means_density,
    start   = list(b0 = 1.0, amp = 0.7, mu = 0.7, sigma = 0.3),
    lower   = c(b0 = 0.5, amp = 0, mu = 0, sigma = 0.1),
    upper   = c(b0 = 1.5, amp = 2, mu = 3, sigma = 3.0),
    weights = 1 / means_density$se_n^2,
    control = lm_control
)
cat("\n--- With all data ---\n")
print(coef(fit_all))

# Find the outlier plant and remove it, recalculate means
# First, identify the outlier at the peak SOM treatment
peak_plants <- ag_max %>% filter(som == peak_som) %>% arrange(desc(n_norm))
cat("\nPlants at peak SOM (", peak_som, ") sorted by normalised leaf count:\n")
print(as.data.frame(peak_plants))

outlier_id <- peak_plants$ID[1]  # highest value
cat("\nCandidate outlier plant:", outlier_id, " with n_norm =", peak_plants$n_norm[1], "\n")

# Check if this plant is a Tukey outlier
q3_peak <- quantile(peak_plants$n_norm, 0.75)
iqr_peak <- IQR(peak_plants$n_norm)
fence_peak <- q3_peak + 1.5 * iqr_peak
cat("Tukey fence at peak SOM:", fence_peak, "\n")
cat("Is outlier by Tukey criterion:", peak_plants$n_norm[1] > fence_peak, "\n")

# 2) Remove outlier, recalculate treatment means, refit
ag_max_clean <- ag_max %>% filter(ID != outlier_id)
means_density_clean <- ag_max_clean %>%
    group_by(som) %>%
    summarize(mean_n = mean(n_norm), se_n = sd(n_norm) / sqrt(n()), n = n(), .groups = "drop")

cat("\n--- Density treatment means WITHOUT outlier ---\n")
print(as.data.frame(means_density_clean))

fit_clean <- nlsLM(
    mean_n ~ b0 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
    data    = means_density_clean,
    start   = list(b0 = 1.0, amp = 0.7, mu = 0.7, sigma = 0.3),
    lower   = c(b0 = 0.5, amp = 0, mu = 0, sigma = 0.1),
    upper   = c(b0 = 1.5, amp = 2, mu = 3, sigma = 3.0),
    weights = 1 / means_density_clean$se_n^2,
    control = lm_control
)
cat("\n--- Without outlier ---\n")
print(coef(fit_clean))

cat("\n--- Comparison ---\n")
cf_all <- coef(fit_all)
cf_clean <- coef(fit_clean)
comp <- data.frame(
    Parameter = names(cf_all),
    With_all = as.numeric(cf_all),
    Without_outlier = as.numeric(cf_clean),
    Pct_change = round(100 * (as.numeric(cf_clean) - as.numeric(cf_all)) / as.numeric(cf_all), 1)
)
print(comp)

# Convert to SOM (Van Bemmelen x1.724)
cat("\n--- Living Dunes parameters (SOM units, Van Bemmelen x1.724) ---\n")
cat("WITH outlier:    mean =", cf_all["mu"] * 1.724, " std =", abs(cf_all["sigma"]) * 1.724, " amplitude =", cf_all["amp"], "\n")
cat("WITHOUT outlier: mean =", cf_clean["mu"] * 1.724, " std =", abs(cf_clean["sigma"]) * 1.724, " amplitude =", cf_clean["amp"], "\n")

# Also check: what if we use median instead of mean for the peak treatment
cat("\n--- Sensitivity: median-based mean at peak treatment ---\n")
peak_median <- median(ag_max$n_norm[ag_max$som == peak_som])
peak_mean_all <- mean(ag_max$n_norm[ag_max$som == peak_som])
peak_mean_clean <- mean(ag_max_clean$n_norm[ag_max_clean$som == peak_som])
cat("Peak SOM treatment (", peak_som, "):\n")
cat("  Mean (all):          ", peak_mean_all, "\n")
cat("  Mean (no outlier):   ", peak_mean_clean, "\n")
cat("  Median (all):        ", peak_median, "\n")
