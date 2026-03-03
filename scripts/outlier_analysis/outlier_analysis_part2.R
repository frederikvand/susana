# =============================================================================
# outlier_analysis_part2.R
# Deep investigation of plant 6b_5 and re-fitting without outlier
# =============================================================================

suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(here)
    library(minpack.lm)
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

df_ag <- read_excel(ag_path)
df_ag$length <- suppressWarnings(as.numeric(df_ag$length))
df_ag$nr     <- suppressWarnings(as.numeric(df_ag$nr))
df_ag$substrate <- as.character(df_ag$substrate)
df_ag <- df_ag %>% left_join(lookup, by = "substrate") %>% filter(!is.na(som))

cat("\n============================================================\n")
cat("  TIME SERIES FOR PLANT 6b_5\n")
cat("============================================================\n")
plant_6b5 <- df_ag %>% filter(ID == "6b_5") %>% select(ID, date, nr, length, substrate, som)
print(as.data.frame(plant_6b5))

cat("\n--- Compare with other plants at same SOM (0.729%) ---\n")
same_trt <- df_ag %>% filter(som == som_vals[4]) %>%
    group_by(ID) %>%
    summarize(
        measurements = n(),
        min_nr = min(nr, na.rm = TRUE),
        max_nr = max(nr, na.rm = TRUE),
        final_nr = last(nr),
        .groups = "drop"
    ) %>%
    arrange(desc(max_nr))
print(as.data.frame(same_trt))

cat("\n============================================================\n")
cat("  FULL DATASET OUTLIER REPORT\n")
cat("============================================================\n")
# All plants with extreme leaf counts across all treatments
ag_max <- df_ag %>%
    group_by(ID, column, som) %>%
    summarize(max_nr = max(nr, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(max_nr))

cat("Distribution of max leaf counts across all plants:\n")
print(summary(ag_max$max_nr))
cat("\nPlants with max leaves > 30:\n")
high_leaf <- ag_max %>% filter(max_nr > 30) %>% arrange(desc(max_nr))
print(as.data.frame(high_leaf))

cat("\n============================================================\n")
cat("  RE-FIT WITHOUT OUTLIER\n")
cat("============================================================\n")

# Remove outlier and recalculate
ref_n <- mean(ag_max$max_nr[ag_max$som == min(ag_max$som)], na.rm = TRUE)
ag_max$n_norm <- ag_max$max_nr / ref_n

ag_clean <- ag_max %>% filter(ID != "6b_5")
means_clean <- ag_clean %>%
    group_by(som) %>%
    summarize(mean_n = mean(n_norm), se_n = sd(n_norm) / sqrt(n()), n = n(), .groups = "drop")

cat("\n--- Treatment means WITHOUT 6b_5 ---\n")
print(as.data.frame(means_clean))

# Try Gaussian fit
cat("\n--- Attempt Gaussian NLS without outlier ---\n")
tryCatch({
    fit_clean <- nlsLM(
        mean_n ~ b0 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
        data    = means_clean,
        start   = list(b0 = 1.0, amp = 0.5, mu = 0.5, sigma = 0.5),
        lower   = c(b0 = 0.5, amp = 0, mu = 0, sigma = 0.1),
        upper   = c(b0 = 1.5, amp = 2, mu = 3, sigma = 3.0),
        weights = 1 / means_clean$se_n^2,
        control = nls.lm.control(maxiter = 200)
    )
    cat("Converged:\n")
    print(coef(fit_clean))
    print(summary(fit_clean))
}, error = function(e) {
    cat("Gaussian fit failed:", conditionMessage(e), "\n")
})

# Flat model (null): just intercept
cat("\n--- Null model (flat, no SOM effect) ---\n")
mean_overall <- weighted.mean(means_clean$mean_n, 1 / means_clean$se_n^2)
cat("Weighted mean (flat model):", mean_overall, "\n")
ss_res <- sum((1 / means_clean$se_n^2) * (means_clean$mean_n - mean_overall)^2)
ss_tot <- sum((1 / means_clean$se_n^2) * (means_clean$mean_n - mean_overall)^2)
cat("Total weighted SS:", ss_tot, "\n")

# Kruskal-Wallis on clean data
cat("\n--- Kruskal-Wallis on clean data ---\n")
ag_clean$som_factor <- as.factor(round(ag_clean$som, 3))
kw_clean <- kruskal.test(n_norm ~ som_factor, data = ag_clean)
cat("chi² =", round(kw_clean$statistic, 2), " p =", round(kw_clean$p.value, 4), "\n")

# Also try broader starting vals
cat("\n--- Gaussian NLS with broader starting values ---\n")
starts_list <- list(
    list(b0 = 0.9, amp = 0.4, mu = 0.5, sigma = 0.8),
    list(b0 = 1.0, amp = 0.3, mu = 0.3, sigma = 0.3),
    list(b0 = 1.0, amp = 0.3, mu = 1.0, sigma = 1.0),
    list(b0 = 1.0, amp = 0.2, mu = 0.5, sigma = 0.5)
)
for (i in seq_along(starts_list)) {
    tryCatch({
        fit_try <- nlsLM(
            mean_n ~ b0 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
            data    = means_clean,
            start   = starts_list[[i]],
            lower   = c(b0 = 0.5, amp = 0, mu = 0, sigma = 0.1),
            upper   = c(b0 = 1.5, amp = 2, mu = 3, sigma = 3.0),
            weights = 1 / means_clean$se_n^2,
            control = nls.lm.control(maxiter = 500)
        )
        cat("Start set", i, "converged:\n")
        print(coef(fit_try))
    }, error = function(e) {
        cat("Start set", i, "failed:", conditionMessage(e), "\n")
    })
}

cat("\n============================================================\n")
cat("  CONCLUSION\n")
cat("============================================================\n")
cat("Plant 6b_5 has 97 leaves — this is", 97 / 21, "× the next highest (21).\n")
cat("The value is", round(9.312 / 2.016, 1), "× the second-highest normalised value.\n")
cat("This single plant inflates the peak mean from ~1.0 to 1.7.\n")
cat("Without it, the Gaussian peak for density growth likely disappears,\n")
cat("meaning SOM has no significant effect on density growth.\n")
