# Script to derive SOM response curves perfectly for Living Dunes JSON config
# using R and GLM approaches.

library(readxl)
library(dplyr)
library(tidyr)
library(minpack.lm) # For non-linear Gaussian fits if needed

# 1. Load Data
# ------------------------------------------------------------------------------
substrate_path <- "d:/projects/susana/data/substrate/calculated_substrate_properties.xlsx"
ag_path <- "d:/projects/susana/data/vegetation_development/common_garden/aboveground_biomass/aboveground_biomass_zeebrugge.xlsx"
bg_path1 <- "d:/projects/susana/data/vegetation_development/common_garden/belowground_biomass/common_garden_roots.xlsx"

df_sub <- read_excel(substrate_path) %>%
    rename(substrate = Substrate, som = `TOC (%)`)

# Substrate mapping to match common garden labels
sub_map <- c(
    "Ref" = "R", "White" = "W", "Dark" = "D",
    "12.5% D" = "12.5", "25% D" = "25.0", "37.5% D" = "37.5", "50% D" = "50.0",
    "62.5% D" = "62.5", "75% D" = "75.0", "87.5% D" = "87.5"
)

df_sub$substrate <- sub_map[df_sub$substrate]
df_sub$substrate <- as.numeric(ifelse(df_sub$substrate %in% c("R", "W", "D"), NA, df_sub$substrate))
# For R, W, D we can manually add them back or map explicitly for the merge later.
# Let's create a lookup table
df_sub <- read_excel(substrate_path) %>% rename(substrate = Substrate, som = `TOC (%)`)
lookup <- data.frame(
    substrate = c("R", "W", "D", "12.5", "25", "37.5", "50", "62.5", "75", "87.5"),
    som = df_sub$som
)

# Load Aboveground
df_ag <- read_excel(ag_path)
df_ag$length <- as.numeric(df_ag$length)
df_ag$nr <- as.numeric(df_ag$nr)
df_ag$substrate <- as.character(df_ag$substrate)

df_ag <- df_ag %>%
    left_join(lookup, by = "substrate") %>%
    filter(!is.na(som))

# Load Belowground
df_bg <- read_excel(bg_path1)
df_bg$substrate <- as.character(df_bg$substrate)
df_bg <- df_bg %>%
    left_join(lookup, by = "substrate") %>%
    filter(!is.na(som))


# 2. Extract Response Metrics
# ------------------------------------------------------------------------------

# Height (Length) & Density (Nr of leaves) Growth
ag_max <- df_ag %>%
    group_by(ID, som) %>%
    summarize(
        max_length = max(length, na.rm = TRUE),
        max_nr = max(nr, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    filter(is.finite(max_length))

# Root Growth (Belowground Biomass)
bg_max <- df_bg %>%
    group_by(som) %>%
    summarize(
        mean_root = mean(root_biomass, na.rm = TRUE)
    )

# Mortality proxy (did the plant die before the end?)
# A plant is considered dead if its final 'nr' is 0 or missing, while it had 'nr' > 0 before.
# We will use quasibinomial on the survival rate per substrate.
mortality_df <- df_ag %>%
    group_by(ID, som) %>%
    summarize(
        survived = ifelse(last(nr) > 0, 1, 0),
        .groups = "drop"
    ) %>%
    group_by(som) %>%
    summarize(
        survival_prob = mean(survived, na.rm = TRUE),
        n = n()
    )

# 3. Fit Models
# ------------------------------------------------------------------------------

# A. MORTALITY (Quasibinomial GLM)
# --------------------------------
# P(Survival) = 1 / (1 + exp(-(b0 + b1 * SOM)))
# Mortality = 1 - P(Survival) = 1 / (1 + exp(b0 + b1 * SOM))
# Living Dunes "sigmoid" takes: 1 / (1 + exp(-steepness * (x - threshold)))
# So, for Mortality: steepness = -b1, threshold = -b0/b1
glm_mort <- glm(survival_prob ~ som, family = quasibinomial(link = "logit"), weights = n, data = mortality_df)
b0 <- coef(glm_mort)[1]
b1 <- coef(glm_mort)[2]

mort_steepness <- unname(-b1)
mort_threshold <- unname(-b0 / b1)

cat("\n--- MORTALITY (Quasibinomial GLM) ---\n")
print(summary(glm_mort))
cat("JSON sigmoid parameters for Mortality:\n")
cat("steepness:", mort_steepness, "\n")
cat("threshold:", mort_threshold, "\n")


# B. GROWTH - HEIGHT (Gaussian / NLS)
# -----------------------------------
# Since growth is non-linear and peaks in the middle, Gaussian is better than a standard linear model.
# But let's fit a Gaussian to the mean max_length per SOM.
mean_ag <- ag_max %>%
    group_by(som) %>%
    summarize(l = mean(max_length), n = mean(max_nr))

# Normalize against lowest SOM (Ref)
ref_l <- mean_ag$l[which.min(mean_ag$som)]
mean_ag$l_norm <- mean_ag$l / ref_l

fit_height <- nlsLM(l_norm ~ 1 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
    data = mean_ag,
    start = list(amp = 0.5, mu = 1.0, sigma = 0.5)
)

cat("\n--- HEIGHT GROWTH (Gaussian NLS) ---\n")
print(summary(fit_height))
cat("JSON gaussian parameters for Height Growth:\n")
cat("mean:", coef(fit_height)["mu"], "\n")
cat("std:", coef(fit_height)["sigma"], "\n")
cat("multiplier (amp):", coef(fit_height)["amp"], "\n")


# C. GROWTH - DENSITY/LEAVES (Gaussian / NLS)
# -----------------------------------------
ref_n <- mean_ag$n[which.min(mean_ag$som)]
mean_ag$n_norm <- mean_ag$n / ref_n

fit_density <- nlsLM(n_norm ~ 1 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
    data = mean_ag,
    start = list(amp = 0.5, mu = 1.0, sigma = 0.5)
)

cat("\n--- DENSITY GROWTH (Gaussian NLS) ---\n")
print(summary(fit_density))
cat("JSON gaussian parameters for Density Growth:\n")
cat("mean:", coef(fit_density)["mu"], "\n")
cat("std:", coef(fit_density)["sigma"], "\n")
cat("multiplier (amp):", coef(fit_density)["amp"], "\n")


# D. ROOT GROWTH (Gaussian / NLS)
# --------------------------------
ref_root <- bg_max$mean_root[which.min(bg_max$som)]
bg_max$r_norm <- bg_max$mean_root / ref_root

fit_root <- nlsLM(r_norm ~ 1 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
    data = bg_max,
    start = list(amp = 0.5, mu = 1.0, sigma = 0.5)
)

cat("\n--- ROOT GROWTH (Gaussian NLS) ---\n")
print(summary(fit_root))
cat("JSON gaussian parameters for Root Growth:\n")
cat("mean:", coef(fit_root)["mu"], "\n")
cat("std:", coef(fit_root)["sigma"], "\n")
cat("multiplier (amp):", coef(fit_root)["amp"], "\n")


# Optional: Save plots
png("d:/projects/susana/output/som_response_R_derivation.png", width = 1200, height = 800, res = 150)
par(mfrow = c(2, 2))

# Plot Mort
plot(mortality_df$som, 1 - mortality_df$survival_prob, main = "Mortality Proxy", xlab = "SOM (%)", ylab = "Mortality Rate")
curve(1 / (1 + exp(-mort_steepness * (x - mort_threshold))), add = TRUE, col = "red")

# Plot Height
plot(mean_ag$som, mean_ag$l_norm, main = "Height Growth Multiplier", xlab = "SOM (%)", ylab = "Growth Ratio")
lines(mean_ag$som, predict(fit_height), col = "green")

# Plot Density
plot(mean_ag$som, mean_ag$n_norm, main = "Density Growth Multiplier", xlab = "SOM (%)", ylab = "Density Ratio")
lines(mean_ag$som, predict(fit_density), col = "blue")

# Plot Root
plot(bg_max$som, bg_max$r_norm, main = "Root Growth Multiplier", xlab = "SOM (%)", ylab = "Root Biomass Ratio")
lines(bg_max$som, predict(fit_root), col = "purple")

dev.off()
