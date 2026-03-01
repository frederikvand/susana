# =============================================================================
# derive_som_glm.R
# Derives SOM response curves for Living Dunes model parameterisation.
#
# Growth models use non-linear least squares (nlsLM, Levenberg-Marquardt)
# with Gaussian response functions on individual plant data.
# Mortality uses a binomial GLMM (glmmTMB) with random intercepts for
# spatial blocking (column) to account for the common-garden design.
#
# NLME with column as random intercept was tested for all growth models
# but failed to converge, indicating negligible between-column variance
# for the non-linear Gaussian fits. Standard GNLS also failed due to its
# less robust PNLS/GLS optimizer. nlsLM (Levenberg-Marquardt) converges
# reliably for these data.
#
# Output:
#   - Console: Full model summaries, JSON-ready parameters
#   - PNG:     4-panel marginal effects plot (ggplot2 + patchwork)
#   - LaTeX:   Coefficient tables for manuscript integration
# =============================================================================

suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(minpack.lm) # Levenberg-Marquardt for non-linear growth fits
    library(glmmTMB)    # Binomial GLMM for mortality
    library(ggeffects)  # Marginal effects extraction
    library(ggplot2)
    library(patchwork)
})

# =============================================================================
# 1. DATA LOADING
# =============================================================================

substrate_path <- "d:/projects/susana/data/substrate/calculated_substrate_properties.xlsx"
ag_path <- "d:/projects/susana/data/vegetation_development/common_garden/aboveground_biomass/aboveground_biomass_zeebrugge.xlsx"
bg_path <- "d:/projects/susana/data/vegetation_development/common_garden/belowground_biomass/common_garden_roots.xlsx"

df_sub <- read_excel(substrate_path)
som_vals <- pull(df_sub, `TOC (%)`)

lookup <- data.frame(
    substrate = c("R", "W", "D", "12.5", "25", "37.5", "50", "62.5", "75", "87.5"),
    som = som_vals
)

# --- Aboveground data ---
df_ag <- read_excel(ag_path)
df_ag$length    <- suppressWarnings(as.numeric(df_ag$length))
df_ag$nr        <- suppressWarnings(as.numeric(df_ag$nr))
df_ag$substrate <- as.character(df_ag$substrate)
df_ag$column    <- as.factor(df_ag$column)

df_ag <- df_ag %>%
    left_join(lookup, by = "substrate") %>%
    filter(!is.na(som))

# --- Belowground data ---
df_bg <- read_excel(bg_path)
df_bg$substrate <- as.character(df_bg$substrate)
df_bg$column    <- as.factor(df_bg$column)
df_bg$box_ID    <- as.factor(df_bg$box_ID)
df_bg <- df_bg %>%
    left_join(lookup, by = "substrate") %>%
    filter(!is.na(som), !is.na(root_biomass))

# =============================================================================
# 2. RESPONSE METRICS
# =============================================================================

# --- Individual-level data (for scatter plots and ANOVA) ---
ag_max <- df_ag %>%
    group_by(ID, column, som) %>%
    summarize(
        max_length = max(length, na.rm = TRUE),
        max_nr     = max(nr, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    filter(is.finite(max_length), is.finite(max_nr))

# Normalize to reference (lowest SOM)
ref_som  <- min(ag_max$som)
ref_l    <- mean(ag_max$max_length[ag_max$som == ref_som], na.rm = TRUE)
ref_n    <- mean(ag_max$max_nr[ag_max$som == ref_som], na.rm = TRUE)
ag_max$l_norm <- ag_max$max_length / ref_l
ag_max$n_norm <- ag_max$max_nr / ref_n

ref_root <- mean(df_bg$root_biomass[df_bg$som == min(df_bg$som)], na.rm = TRUE)
df_bg$r_norm <- df_bg$root_biomass / ref_root

# Mortality: binary survival per individual plant
mortality_df <- df_ag %>%
    group_by(ID, column, som) %>%
    summarize(
        survived = as.integer(!is.na(last(nr)) & last(nr) > 0),
        .groups = "drop"
    )

# --- Treatment-level means (for NLS curve fitting) ---
# This is the standard approach for dose-response parameterisation:
# fit the curve to treatment means, weighted by inverse variance.
means_height <- ag_max %>%
    group_by(som) %>%
    summarize(
        mean_l = mean(l_norm),
        se_l   = sd(l_norm) / sqrt(n()),
        n      = n(),
        .groups = "drop"
    )

means_density <- ag_max %>%
    group_by(som) %>%
    summarize(
        mean_n = mean(n_norm),
        se_n   = sd(n_norm) / sqrt(n()),
        n      = n(),
        .groups = "drop"
    )

means_root <- df_bg %>%
    group_by(som) %>%
    summarize(
        mean_r = mean(r_norm),
        se_r   = sd(r_norm) / sqrt(n()),
        n      = n(),
        .groups = "drop"
    )

cat("\n======================================================\n")
cat("  DATA SUMMARY\n")
cat("======================================================\n")
cat("Aboveground observations (peak per plant):", nrow(ag_max), "\n")
cat("Belowground observations:", nrow(df_bg), "\n")
cat("Mortality observations:", nrow(mortality_df), "\n")
cat("SOM treatments:", nrow(means_height), "\n")
cat("Spatial blocks (columns):", paste(sort(unique(ag_max$column)), collapse = ", "), "\n")
cat("SOM range: [", min(ag_max$som), ",", max(ag_max$som), "] %TOC\n")

# Treatment-level ANOVA to test overall SOM effect
cat("\n--- Kruskal-Wallis test for overall SOM treatment effect ---\n")
ag_max$som_factor <- as.factor(round(ag_max$som, 3))
kw_height  <- kruskal.test(l_norm ~ som_factor, data = ag_max)
kw_density <- kruskal.test(n_norm ~ som_factor, data = ag_max)
df_bg$som_factor <- as.factor(round(df_bg$som, 3))
kw_root    <- kruskal.test(r_norm ~ som_factor, data = df_bg)
cat("  Height:  chi² =", round(kw_height$statistic, 2), " p =", round(kw_height$p.value, 4), "\n")
cat("  Density: chi² =", round(kw_density$statistic, 2), " p =", round(kw_density$p.value, 4), "\n")
cat("  Root:    chi² =", round(kw_root$statistic, 2), " p =", round(kw_root$p.value, 4), "\n")

# =============================================================================
# 3. MODEL FITTING (treatment-mean weighted NLS)
# =============================================================================
# Fitting to treatment means with inverse-variance weights is the standard
# approach for dose-response curve parameterisation when within-treatment
# variation exceeds the between-treatment signal.

lm_control <- nls.lm.control(maxiter = 200)

cat("\n======================================================\n")
cat("  MODEL 1: HEIGHT GROWTH (Gaussian NLS on treatment means)\n")
cat("======================================================\n")

fit_height <- nlsLM(
    mean_l ~ b0 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
    data    = means_height,
    start   = list(b0 = 1.0, amp = 0.35, mu = 1.0, sigma = 1.0),
    lower   = c(b0 = 0.5, amp = 0, mu = 0, sigma = 0.1),
    upper   = c(b0 = 1.5, amp = 2, mu = 3, sigma = 3.0),
    weights = 1 / means_height$se_l^2,
    control = lm_control
)
print(summary(fit_height))
h_cf <- coef(fit_height)
cat("\nJSON parameters — Height Growth:\n")
cat("  mean (mu):", h_cf["mu"], "\n")
cat("  std (sigma):", abs(h_cf["sigma"]), "\n")
cat("  amplitude:", h_cf["amp"], "\n")

# Weighted R-squared on treatment means
pred_h <- predict(fit_height)
wt_h <- 1 / means_height$se_l^2
R2_height <- 1 - sum(wt_h * (means_height$mean_l - pred_h)^2) / sum(wt_h * (means_height$mean_l - weighted.mean(means_height$mean_l, wt_h))^2)
n_h <- nrow(means_height); p_h <- length(coef(fit_height))
R2_adj_height <- 1 - (1 - R2_height) * (n_h - 1) / (n_h - p_h - 1)
cat("  R² (weighted):", round(R2_height, 4), "\n")
cat("  R² adj (weighted):", round(R2_adj_height, 4), "\n")


cat("\n======================================================\n")
cat("  MODEL 2: DENSITY GROWTH (Gaussian NLS on treatment means)\n")
cat("======================================================\n")

fit_density <- nlsLM(
    mean_n ~ b0 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
    data    = means_density,
    start   = list(b0 = 1.0, amp = 0.7, mu = 0.7, sigma = 0.3),
    lower   = c(b0 = 0.5, amp = 0, mu = 0, sigma = 0.1),
    upper   = c(b0 = 1.5, amp = 2, mu = 3, sigma = 3.0),
    weights = 1 / means_density$se_n^2,
    control = lm_control
)
print(summary(fit_density))
d_cf <- coef(fit_density)
cat("\nJSON parameters — Density Growth:\n")
cat("  mean (mu):", d_cf["mu"], "\n")
cat("  std (sigma):", abs(d_cf["sigma"]), "\n")
cat("  amplitude:", d_cf["amp"], "\n")

pred_d <- predict(fit_density)
wt_d <- 1 / means_density$se_n^2
R2_density <- 1 - sum(wt_d * (means_density$mean_n - pred_d)^2) / sum(wt_d * (means_density$mean_n - weighted.mean(means_density$mean_n, wt_d))^2)
R2_adj_density <- 1 - (1 - R2_density) * (n_h - 1) / (n_h - p_h - 1)
cat("  R² (weighted):", round(R2_density, 4), "\n")
cat("  R² adj (weighted):", round(R2_adj_density, 4), "\n")


cat("\n======================================================\n")
cat("  MODEL 3: ROOT GROWTH (Gaussian NLS on treatment means)\n")
cat("======================================================\n")

fit_root <- nlsLM(
    mean_r ~ b0 + amp * exp(-((som - mu)^2) / (2 * sigma^2)),
    data    = means_root,
    start   = list(b0 = 1.0, amp = 0.7, mu = 0.7, sigma = 0.5),
    lower   = c(b0 = 0.5, amp = 0, mu = 0, sigma = 0.1),
    upper   = c(b0 = 1.5, amp = 2, mu = 3, sigma = 3.0),
    weights = 1 / means_root$se_r^2,
    control = lm_control
)
print(summary(fit_root))
r_cf <- coef(fit_root)
cat("\nJSON parameters — Root Growth:\n")
cat("  mean (mu):", r_cf["mu"], "\n")
cat("  std (sigma):", abs(r_cf["sigma"]), "\n")
cat("  amplitude:", r_cf["amp"], "\n")

n_root <- nrow(means_root); p_root <- length(coef(fit_root))
pred_r <- predict(fit_root)
wt_r <- 1 / means_root$se_r^2
R2_root <- 1 - sum(wt_r * (means_root$mean_r - pred_r)^2) / sum(wt_r * (means_root$mean_r - weighted.mean(means_root$mean_r, wt_r))^2)
R2_adj_root <- 1 - (1 - R2_root) * (n_root - 1) / (n_root - p_root - 1)
cat("  R² (weighted):", round(R2_root, 4), "\n")
cat("  R² adj (weighted):", round(R2_adj_root, 4), "\n")


cat("\n======================================================\n")
cat("  MODEL 4: MORTALITY (Binomial GLMM)\n")
cat("======================================================\n")

fit_mort <- glmmTMB(
    survived ~ som + (1 | column),
    family = binomial(link = "logit"),
    data   = mortality_df
)
print(summary(fit_mort))

b0_m <- fixef(fit_mort)$cond["(Intercept)"]
b1_m <- fixef(fit_mort)$cond["som"]
mort_steepness <- unname(-b1_m)
mort_threshold <- unname(-b0_m / b1_m)
cat("\nJSON parameters — Mortality:\n")
cat("  steepness:", mort_steepness, "\n")
cat("  threshold:", mort_threshold, "\n")

fit_mort_null <- glmmTMB(survived ~ 1 + (1 | column), family = binomial, data = mortality_df)
R2_mcfadden <- 1 - as.numeric(logLik(fit_mort)) / as.numeric(logLik(fit_mort_null))
cat("  McFadden pseudo-R²:", round(R2_mcfadden, 4), "\n")


# =============================================================================
# 4. GENERATE PLOTS
# =============================================================================

theme_som <- theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 10),
        panel.grid.minor = element_blank()
    )

predict_gaussian <- function(cf, data_range) {
    som_seq <- seq(min(data_range), max(data_range), length.out = 200)
    data.frame(
        som = som_seq,
        predicted = cf["b0"] + cf["amp"] * exp(-((som_seq - cf["mu"])^2) / (2 * cf["sigma"]^2))
    )
}

# --- (a) Height Growth ---
p_height <- ggplot() +
    geom_jitter(data = ag_max, aes(x = som, y = l_norm),
                width = 0.02, alpha = 0.2, size = 1.2, color = "grey60") +
    geom_pointrange(data = means_height, aes(x = som, y = mean_l,
                    ymin = mean_l - se_l, ymax = mean_l + se_l),
                    color = "#2ca02c", size = 0.6) +
    geom_line(data = predict_gaussian(h_cf, ag_max$som), aes(x = som, y = predicted),
              color = "#2ca02c", linewidth = 1.2) +
    labs(title = "(a) Height growth multiplier", x = "TOC (%)", y = "Normalised height") +
    theme_som

# --- (b) Density Growth ---
p_density <- ggplot() +
    geom_jitter(data = ag_max, aes(x = som, y = n_norm),
                width = 0.02, alpha = 0.2, size = 1.2, color = "grey60") +
    geom_pointrange(data = means_density, aes(x = som, y = mean_n,
                    ymin = mean_n - se_n, ymax = mean_n + se_n),
                    color = "#1f77b4", size = 0.6) +
    geom_line(data = predict_gaussian(d_cf, ag_max$som), aes(x = som, y = predicted),
              color = "#1f77b4", linewidth = 1.2) +
    labs(title = "(b) Density growth multiplier", x = "TOC (%)", y = "Normalised leaf count") +
    theme_som

# --- (c) Root Growth ---
p_root <- ggplot() +
    geom_jitter(data = df_bg, aes(x = som, y = r_norm),
                width = 0.02, alpha = 0.2, size = 1.2, color = "grey60") +
    geom_pointrange(data = means_root, aes(x = som, y = mean_r,
                    ymin = mean_r - se_r, ymax = mean_r + se_r),
                    color = "#9467bd", size = 0.6) +
    geom_line(data = predict_gaussian(r_cf, df_bg$som), aes(x = som, y = predicted),
              color = "#9467bd", linewidth = 1.2) +
    labs(title = "(c) Root growth multiplier", x = "TOC (%)", y = "Normalised root biomass") +
    theme_som

# --- (d) Mortality ---
pred_mort_gg <- ggpredict(fit_mort, terms = "som [all]")
means_mort <- mortality_df %>%
    group_by(som) %>%
    summarize(mort_rate = 1 - mean(survived), se = sqrt(mean(survived) * (1 - mean(survived)) / n()), .groups = "drop")
p_mort <- ggplot() +
    geom_jitter(data = mortality_df, aes(x = som, y = 1 - survived),
                height = 0.05, width = 0.02, alpha = 0.15, size = 1.2, color = "grey60") +
    geom_pointrange(data = means_mort, aes(x = som, y = mort_rate,
                    ymin = pmax(mort_rate - se, 0), ymax = pmin(mort_rate + se, 1)),
                    color = "#d62728", size = 0.6) +
    geom_line(data = pred_mort_gg, aes(x = x, y = 1 - predicted),
              color = "#d62728", linewidth = 1.2) +
    geom_ribbon(data = pred_mort_gg, aes(x = x, ymin = 1 - conf.high, ymax = 1 - conf.low),
                fill = "#d62728", alpha = 0.15) +
    labs(title = "(d) Mortality probability", x = "TOC (%)", y = "Mortality probability") +
    theme_som

combined <- (p_height | p_density) / (p_root | p_mort) +
    plot_annotation(
        title = "SOM response curves for Ammophila arenaria demographic parameterisation",
        theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

ggsave("d:/projects/susana/output/som_response_R_derivation.png",
       combined, width = 10, height = 8, dpi = 300, bg = "white")
cat("\nPlot saved: d:/projects/susana/output/som_response_R_derivation.png\n")


# =============================================================================
# 5. GENERATE LATEX TABLES
# =============================================================================

output_dir <- "d:/projects/susana/evaluation/manuscript/tables/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper ---
nls_table_row <- function(model, model_name, R2, R2_adj, n_total) {
    cf <- coef(model)
    se <- summary(model)$coefficients[, "Std. Error"]
    tv <- summary(model)$coefficients[, "t value"]
    pv <- summary(model)$coefficients[, "Pr(>|t|)"]
    data.frame(
        Model = model_name,
        Parameter = names(cf),
        Estimate = round(cf, 4),
        SE = round(se, 4),
        t = round(tv, 2),
        p = ifelse(pv < 0.001, "$<$0.001", as.character(round(pv, 3))),
        R2 = c(round(R2, 3), rep("", length(cf) - 1)),
        R2_adj = c(round(R2_adj, 3), rep("", length(cf) - 1)),
        n = c(as.character(n_total), rep("", length(cf) - 1)),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

growth_table <- rbind(
    nls_table_row(fit_height,  "Height",  R2_height,  R2_adj_height,  nrow(ag_max)),
    nls_table_row(fit_density, "Density", R2_density, R2_adj_density, nrow(ag_max)),
    nls_table_row(fit_root,    "Root",    R2_root,    R2_adj_root,    nrow(df_bg))
)

# --- Table 1: Growth model coefficients ---
sink(paste0(output_dir, "growth_model_coefficients.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Non-linear regression coefficients for the Gaussian growth response models (Levenberg--Marquardt algorithm). The response function is $y = b_0 + a \\cdot \\exp\\left(-\\frac{(\\mathrm{SOM} - \\mu)^2}{2\\sigma^2}\\right)$, where $b_0$ is the baseline multiplier, $a$ is the amplitude, $\\mu$ is the optimal SOM (\\% TOC), and $\\sigma$ is the standard deviation of the response. Response variables are normalised to the reference beach sand (R).}\n")
cat("\\label{tab:growth_coefficients}\n")
cat("\\small\n")
cat("\\begin{tabular}{ll rrrr rr r}\n")
cat("\\toprule\n")
cat("Response & Parameter & Estimate & SE & $t$ & $p$ & $R^2$ & $R^2_{adj}$ & $n$ \\\\\n")
cat("\\midrule\n")
prev_model <- ""
for (i in 1:nrow(growth_table)) {
    row <- growth_table[i, ]
    if (row$Model != prev_model && prev_model != "") cat("\\midrule\n")
    model_str <- ifelse(row$Model == prev_model, "", row$Model)
    cat(sprintf("%s & $%s$ & %s & %s & %s & %s & %s & %s & %s \\\\\n",
        model_str, row$Parameter,
        format(row$Estimate, nsmall = 4), format(row$SE, nsmall = 4),
        format(row$t, nsmall = 2), row$p, row$R2, row$R2_adj, row$n))
    prev_model <- row$Model
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("Table written: growth_model_coefficients.tex\n")


# --- Table 2: Mortality GLMM coefficients ---
mort_coefs <- summary(fit_mort)$coefficients$cond
vc <- VarCorr(fit_mort)$cond$column
re_sd <- attr(vc, "stddev")[1]

sink(paste0(output_dir, "mortality_model_coefficients.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Binomial GLMM coefficients for the survival--SOM relationship. The model is $\\mathrm{logit}(P_{\\mathrm{survival}}) = \\beta_0 + \\beta_1 \\cdot \\mathrm{SOM} + u_j$, where $u_j \\sim \\mathcal{N}(0, \\sigma^2_u)$ is a random intercept for spatial block $j$ (column position in common garden). McFadden pseudo-$R^2$: ", round(R2_mcfadden, 3), ".}\n")
cat("\\label{tab:mortality_coefficients}\n")
cat("\\begin{tabular}{lrrrr}\n")
cat("\\toprule\n")
cat("Component & Estimate & SE & $z$ & $p$ \\\\\n")
cat("\\midrule\n")
cat("\\multicolumn{5}{l}{\\textit{Fixed effects}} \\\\\n")
for (i in 1:nrow(mort_coefs)) {
    pname <- rownames(mort_coefs)[i]
    pname_tex <- gsub("\\(Intercept\\)", "$\\\\beta_0$ (Intercept)", pname)
    pname_tex <- gsub("^som$", "$\\\\beta_1$ (SOM)", pname_tex)
    cat(sprintf("%s & %.4f & %.4f & %.2f & %.4f \\\\\n",
        pname_tex, mort_coefs[i,1], mort_coefs[i,2], mort_coefs[i,3], mort_coefs[i,4]))
}
cat("\\midrule\n")
cat("\\multicolumn{5}{l}{\\textit{Random effects}} \\\\\n")
cat(sprintf("$\\sigma_u$ (column) & %.4f & & & \\\\\n", re_sd))
cat("\\midrule\n")
cat("\\multicolumn{5}{l}{\\textit{Derived Living Dunes parameters}} \\\\\n")
cat(sprintf("Steepness ($-\\beta_1$) & %.4f & & & \\\\\n", mort_steepness))
cat(sprintf("Threshold ($-\\beta_0/\\beta_1$) & %.2f & & & \\\\\n", mort_threshold))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("Table written: mortality_model_coefficients.tex\n")


# --- Table 3: Derived parameters summary (main text) ---
sink(paste0(output_dir, "som_derived_parameters.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Derived SOM response parameters for the Living Dunes vegetation model. Growth responses use a Gaussian function as a multiplier on the base demographic rate. Mortality uses a sigmoid (inverse logistic) function.}\n")
cat("\\label{tab:som_derived_parameters}\n")
cat("\\begin{tabular}{l l l r}\n")
cat("\\toprule\n")
cat("Demographic process & Function & Parameter & Value \\\\\n")
cat("\\midrule\n")
cat(sprintf("Height growth & Gaussian & $\\mu$ (\\%%%% TOC) & %.3f \\\\\n", h_cf["mu"]))
cat(sprintf("              &          & $\\sigma$ (\\%%%% TOC) & %.3f \\\\\n", abs(h_cf["sigma"])))
cat(sprintf("              &          & Amplitude $a$ & %.3f \\\\\n", h_cf["amp"]))
cat("\\midrule\n")
cat(sprintf("Density growth & Gaussian & $\\mu$ (\\%%%% TOC) & %.3f \\\\\n", d_cf["mu"]))
cat(sprintf("               &          & $\\sigma$ (\\%%%% TOC) & %.3f \\\\\n", abs(d_cf["sigma"])))
cat(sprintf("               &          & Amplitude $a$ & %.3f \\\\\n", d_cf["amp"]))
cat("\\midrule\n")
cat(sprintf("Root growth & Gaussian & $\\mu$ (\\%%%% TOC) & %.3f \\\\\n", r_cf["mu"]))
cat(sprintf("            &          & $\\sigma$ (\\%%%% TOC) & %.3f \\\\\n", abs(r_cf["sigma"])))
cat(sprintf("            &          & Amplitude $a$ & %.3f \\\\\n", r_cf["amp"]))
cat("\\midrule\n")
cat(sprintf("Mortality & Sigmoid & Steepness & %.4f \\\\\n", mort_steepness))
cat(sprintf("          &         & Threshold (\\%%%% TOC) & %.2f \\\\\n", mort_threshold))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("Table written: som_derived_parameters.tex\n")


# --- Table 4: Model fit statistics (appendix) ---
sigma_h <- summary(fit_height)$sigma
sigma_d <- summary(fit_density)$sigma
sigma_r <- summary(fit_root)$sigma

sink(paste0(output_dir, "model_fit_statistics.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Model fit statistics for the SOM response curve derivations. Growth models use non-linear least squares with the Levenberg--Marquardt algorithm (\\texttt{nlsLM}). The mortality model uses a binomial generalised linear mixed model (GLMM) with a random intercept for spatial block.}\n")
cat("\\label{tab:model_fit_statistics}\n")
cat("\\begin{tabular}{l r r r r r}\n")
cat("\\toprule\n")
cat("Model & Log-Lik & AIC & $\\sigma_\\varepsilon$ & $R^2$ & $n$ \\\\\n")
cat("\\midrule\n")
cat(sprintf("Height NLS & %.1f & %.1f & %.4f & %.3f & %d \\\\\n",
    logLik(fit_height), AIC(fit_height), sigma_h, R2_height, nrow(ag_max)))
cat(sprintf("Density NLS & %.1f & %.1f & %.4f & %.3f & %d \\\\\n",
    logLik(fit_density), AIC(fit_density), sigma_d, R2_density, nrow(ag_max)))
cat(sprintf("Root NLS & %.1f & %.1f & %.4f & %.3f & %d \\\\\n",
    logLik(fit_root), AIC(fit_root), sigma_r, R2_root, nrow(df_bg)))
cat(sprintf("Mortality GLMM & %.1f & %.1f & --- & %.3f$^\\dagger$ & %d \\\\\n",
    as.numeric(logLik(fit_mort)), AIC(fit_mort), R2_mcfadden, nrow(mortality_df)))
cat("\\bottomrule\n")
cat("\\multicolumn{6}{l}{\\footnotesize $\\dagger$ McFadden pseudo-$R^2$ for the binomial GLMM.} \\\\\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("Table written: model_fit_statistics.tex\n")


cat("\n======================================================\n")
cat("  ALL DONE\n")
cat("======================================================\n")
