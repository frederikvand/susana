# =============================================================================
# 02_aboveground_analysis.R
# Common garden aboveground biomass: height & leaf count over time
#
# Statistical approach:
#   Height  — Gaussian LMM with quadratic time (poly(time, 2))
#   Leaves  — Negative Binomial GLMM (glmmTMB, nbinom2) with quadratic time
#   Both    — plant_id random intercept for repeated measures
# Visualisation:
#   Model-predicted trajectories via ggeffects (NOT raw means or LOESS)
# =============================================================================

cat("\n======================================================\n")
cat("  02_aboveground_analysis.R\n")
cat("======================================================\n\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(glmmTMB)
  library(emmeans)
  library(ggeffects)
  library(patchwork)
  library(scales)
  library(multcomp)
  library(here)
})

source(here("scripts", "shared_theme.R"))
set.seed(2024)

# --- Load clean data ---
ag <- readRDS(here("output", "clean_data", "aboveground_planted.rds"))

out_dir <- here("output", "common_garden", "aboveground_biomass")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- here("evaluation", "manuscript", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# =====================================================
# DATA PREPARATION
# =====================================================
cat("--- Data summary ---\n")

ag_growth <- ag %>%
  filter(!is.na(length), length > 0, !is.na(leaf_count))

cat("  Observations (alive):", nrow(ag_growth), "\n")
cat("  Unique plants:", length(unique(ag_growth$plant_id)), "\n")
cat("  Timepoints:", length(unique(ag_growth$timepoint)), "\n")
cat("  Substrates:", length(unique(ag_growth$substrate)), "\n")
cat("  Leaf count  — mean:", round(mean(ag_growth$leaf_count), 2),
    " var:", round(var(ag_growth$leaf_count), 2),
    " (overdispersed:", var(ag_growth$leaf_count) > mean(ag_growth$leaf_count),
    ")\n")

# glmmTMB / emmeans require unordered factors
ag_growth$substrate_f <- factor(ag_growth$substrate, ordered = FALSE)


# =====================================================
# MODEL 1: Height — Gaussian LMM with quadratic time
# =====================================================
cat("\n--- Gaussian LMM: height ~ substrate * poly(time, 2) + (1|plant_id) ---\n")

lmm_height <- lmer(
  length ~ substrate_f * poly(time_numeric, 2) + (1 | plant_id),
  data = ag_growth,
  REML = TRUE
)

# Compare linear vs quadratic (ML for model comparison)
lmm_h_lin  <- lmer(length ~ substrate_f * time_numeric +
                      (1 | plant_id), data = ag_growth, REML = FALSE)
lmm_h_quad <- lmer(length ~ substrate_f * poly(time_numeric, 2) +
                      (1 | plant_id), data = ag_growth, REML = FALSE)
lr_test <- anova(lmm_h_lin, lmm_h_quad)
cat("  Linear AIC:", round(AIC(lmm_h_lin), 1),
    "  Quadratic AIC:", round(AIC(lmm_h_quad), 1), "\n")
cat("  LR test (linear vs quad): Chisq =",
    round(lr_test$Chisq[2], 2), " p =",
    format.pval(lr_test$`Pr(>Chisq)`[2], digits = 3), "\n")

# EMMs averaged over time
emm_height <- emmeans(lmm_height, "substrate_f")
pairs_height <- pairs(emm_height, adjust = "fdr")
cld_height <- cld(emm_height, adjust = "fdr", Letters = letters)
cld_height <- as.data.frame(cld_height) %>%
  mutate(.group = trimws(.group))

cat("\n  CLD for height:\n")
print(cld_height[, c("substrate_f", "emmean", "SE", ".group")])


# =====================================================
# MODEL 2: Leaf count — Negative Binomial GLMM
# =====================================================
cat("\n--- NB GLMM: leaf_count ~ substrate * poly(time, 2) + (1|plant_id) ---\n")

pois_leaves <- glmmTMB(
  leaf_count ~ substrate_f * poly(time_numeric, 2) + (1 | plant_id),
  family = poisson,
  data = ag_growth
)

nb_ok <- tryCatch({
  nb_leaves <- glmmTMB(
    leaf_count ~ substrate_f * poly(time_numeric, 2) + (1 | plant_id),
    family = nbinom2,
    data = ag_growth
  )
  # Check for degenerate NB (theta -> Inf means Poisson is sufficient)
  theta <- sigma(nb_leaves)
  hess_ok <- !is.na(AIC(nb_leaves))
  if (!hess_ok || theta > 1e6) {
    cat("  NB dispersion degenerate (theta =", round(theta, 0),
        ") -> random intercept absorbs overdispersion\n")
    cat("  Falling back to Poisson GLMM\n")
    FALSE
  } else {
    cat("  Poisson AIC:", round(AIC(pois_leaves), 1),
        "  NB AIC:", round(AIC(nb_leaves), 1), "\n")
    cat("  NB dispersion:", round(theta, 2), "\n")
    TRUE
  }
}, error = function(e) {
  cat("  NB GLMM failed:", conditionMessage(e), "\n")
  cat("  Using Poisson GLMM instead\n")
  FALSE
})

# Select best leaf model
if (nb_ok) {
  leaf_model <- nb_leaves
  leaf_family_label <- "Negative Binomial"
} else {
  leaf_model <- pois_leaves
  leaf_family_label <- "Poisson"
}
cat("  Selected leaf model:", leaf_family_label, "GLMM\n")
cat("  AIC:", round(AIC(leaf_model), 1), "\n")

# EMMs on the response (count) scale
emm_leaves <- emmeans(leaf_model, "substrate_f", type = "response")
cld_leaves <- cld(emm_leaves, adjust = "fdr", Letters = letters)
cld_leaves <- as.data.frame(cld_leaves) %>%
  mutate(.group = trimws(.group))

# Normalise column name (Poisson uses 'rate', NB uses 'response')
resp_col <- intersect(names(cld_leaves), c("response", "rate", "emmean"))
if (length(resp_col) > 0) {
  names(cld_leaves)[names(cld_leaves) == resp_col[1]] <- "emm_response"
} else {
  stop("Cannot find response column in CLD output")
}

cat("\n  CLD for leaf count (response scale):\n")
print(cld_leaves[, c("substrate_f", "emm_response", "SE", ".group")])


# =====================================================
# MODEL PREDICTIONS via ggeffects
# =====================================================
cat("\n--- Model predictions (ggeffects) ---\n")

# Unique time values in data
time_vals <- sort(unique(ag_growth$time_numeric))

# Dense prediction grid for smooth curves
time_grid <- seq(min(time_vals), max(time_vals), length.out = 50)

pred_height <- ggpredict(lmm_height,
                         terms = c("time_numeric [time_grid]", "substrate_f"))
pred_height_df <- as.data.frame(pred_height) %>%
  rename(time_numeric = x, substrate = group)

pred_leaves <- ggpredict(leaf_model,
                         terms = c("time_numeric [time_grid]", "substrate_f"),
                         bias_correction = TRUE)
pred_leaves_df <- as.data.frame(pred_leaves) %>%
  rename(time_numeric = x, substrate = group)


# =====================================================
# LaTeX TABLE
# =====================================================
cat("\n--- Writing LaTeX table ---\n")

sink(file.path(table_dir, "aboveground_lmm_results.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Estimated marginal means (EMMs) for aboveground performance metrics.")
cat(" Height was modelled with a Gaussian LMM (quadratic time);")
cat(sprintf(" leaf count with a %s GLMM (\\texttt{glmmTMB}, quadratic time);", leaf_family_label))
cat(" both with plant identity as random intercept.")
cat(" Groups sharing the same letter are not significantly different")
cat(" (FDR-adjusted pairwise contrasts, $p < 0.05$).}\n")
cat("\\label{tab:aboveground_lmm}\n")
cat("\\small\n")
cat("\\begin{tabular}{l rr cc}\n")
cat("\\toprule\n")
cat("Substrate & Height EMM (cm) & Leaf count EMM & Height group & Leaf group \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(cld_height))) {
  s <- as.character(cld_height$substrate_f[i])
  h_emm <- round(cld_height$emmean[i], 1)
  h_grp <- cld_height$.group[i]
  l_row <- cld_leaves[as.character(cld_leaves$substrate_f) == s, ]
  l_emm <- round(l_row$emm_response[1], 1)
  l_grp <- l_row$.group[1]
  cat(sprintf("%-8s & %.1f & %.1f & %s & %s \\\\\n",
              s, h_emm, l_emm, h_grp, l_grp))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("  Table written: aboveground_lmm_results.tex\n")


# =====================================================
# PLOTS: Model-predicted trajectories
# =====================================================
cat("\n--- Plots ---\n")

# --- 1. Height: model fit + raw data (main text) ---
p_height <- ggplot() +
  geom_point(data = ag_growth,
             aes(x = time_numeric, y = length, colour = substrate),
             alpha = 0.08, size = 0.5, show.legend = FALSE) +
  geom_ribbon(data = pred_height_df,
              aes(x = time_numeric, ymin = conf.low, ymax = conf.high,
                  fill = substrate),
              alpha = 0.15) +
  geom_line(data = pred_height_df,
            aes(x = time_numeric, y = predicted, colour = substrate),
            linewidth = 0.8) +
  scale_colour_substrate() +
  scale_fill_substrate() +
  scale_x_continuous(breaks = time_vals,
                     labels = paste0("t", time_vals)) +
  labs(x = "Measurement period", y = "Predicted leaf height (cm)",
       colour = "Substrate", fill = "Substrate") +
  theme_susana() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 2. Leaf count: NB model fit + raw data (main text) ---
p_leaves <- ggplot() +
  geom_point(data = ag_growth,
             aes(x = time_numeric, y = leaf_count, colour = substrate),
             alpha = 0.08, size = 0.5, show.legend = FALSE) +
  geom_ribbon(data = pred_leaves_df,
              aes(x = time_numeric, ymin = conf.low, ymax = conf.high,
                  fill = substrate),
              alpha = 0.15) +
  geom_line(data = pred_leaves_df,
            aes(x = time_numeric, y = predicted, colour = substrate),
            linewidth = 0.8) +
  scale_colour_substrate() +
  scale_fill_substrate() +
  scale_x_continuous(breaks = time_vals,
                     labels = paste0("t", time_vals)) +
  labs(x = "Measurement period", y = "Predicted leaf count",
       colour = "Substrate", fill = "Substrate") +
  theme_susana() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Combined: stacked vertically ---
p_ag_main <- p_height / p_leaves +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "bottom")
save_pub_plot(p_ag_main,
              file.path(out_dir, "aboveground_growth_main.png"),
              width = 180, height = 240)

# Save individual panels for backward compatibility
save_pub_plot(p_height,
              file.path(out_dir, "marram_height_loess.png"),
              width = 180, height = 120)
save_pub_plot(p_leaves,
              file.path(out_dir, "leaves_number_loess.png"),
              width = 180, height = 120)

# --- 3. Detailed boxplots: Height (appendix) ---
p_height_all <- ggplot(ag_growth,
                       aes(x = substrate, y = length, fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_substrate() +
  facet_wrap(~ timepoint, nrow = 2) +
  labs(x = "Substrate", y = "Leaf height (cm)") +
  theme_susana() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
save_pub_plot(p_height_all,
              file.path(out_dir, "marram_height_all.png"),
              width = 260, height = 160)

# --- 4. Detailed boxplots: Leaf count (appendix) ---
p_leaves_all <- ggplot(ag_growth,
                       aes(x = substrate, y = leaf_count, fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_substrate() +
  facet_wrap(~ timepoint, nrow = 2) +
  labs(x = "Substrate", y = "Number of leaves") +
  theme_susana() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
save_pub_plot(p_leaves_all,
              file.path(out_dir, "number_of_leaves.png"),
              width = 260, height = 160)


cat("\n======================================================\n")
cat("  02_aboveground_analysis.R COMPLETE\n")
cat("======================================================\n")
