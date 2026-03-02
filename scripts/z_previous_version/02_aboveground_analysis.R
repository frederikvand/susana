# =============================================================================
# 02_aboveground_analysis.R
# Common garden aboveground biomass: height & leaf count over time, LMMs
# =============================================================================

cat("\n======================================================\n")
cat("  02_aboveground_analysis.R\n")
cat("======================================================\n\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(emmeans)
  library(patchwork)
  library(scales)
})

source("d:/projects/susana/scripts/shared_theme.R")

# --- Load clean data ---
ag <- readRDS("d:/projects/susana/output/clean_data/aboveground_planted.rds")

out_dir <- "d:/projects/susana/output/common_garden/aboveground_biomass"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- "d:/projects/susana/evaluation/manuscript/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# =====================================================
# DATA PREPARATION
# =====================================================
cat("--- Data summary ---\n")

# Filter out dead plants (length == 0 or NA) for growth analysis
ag_growth <- ag %>%
  filter(!is.na(length), length > 0, !is.na(leaf_count))

cat("  Observations (alive):", nrow(ag_growth), "\n")
cat("  Unique plants:", length(unique(ag_growth$plant_id)), "\n")
cat("  Timepoints:", length(unique(ag_growth$timepoint)), "\n")
cat("  Substrates:", length(unique(ag_growth$substrate)), "\n")


# =====================================================
# STATISTICS: Linear Mixed-Effects Models
# =====================================================
cat("\n--- LMM: Height ---\n")

# Height LMM: substrate + time as fixed, plant_id as random intercept
lmm_height <- lmer(length ~ substrate * time_numeric + (1 | plant_id),
                    data = ag_growth,
                    REML = TRUE)
cat("  Height LMM summary:\n")
print(summary(lmm_height)$coefficients[1:5, ])

# Estimated marginal means for substrate effect (averaged over time)
emm_height <- emmeans(lmm_height, "substrate")
cat("\n  EMMs for height by substrate:\n")
print(summary(emm_height))

# Pairwise comparisons with FDR correction
pairs_height <- pairs(emm_height, adjust = "fdr")
cat("\n  Significant pairwise height comparisons (FDR < 0.05):\n")
sig_h <- as.data.frame(pairs_height) %>% filter(p.value < 0.05)
if (nrow(sig_h) > 0) {
  print(sig_h[, c("contrast", "estimate", "p.value")])
} else {
  cat("  None\n")
}

# Compact Letter Display for plots
cld_height <- multcomp::cld(emm_height, adjust = "fdr", Letters = letters)
cld_height <- as.data.frame(cld_height) %>%
  mutate(.group = trimws(.group))
cat("\n  CLD for height:\n")
print(cld_height[, c("substrate", "emmean", ".group")])


cat("\n--- LMM: Leaf count ---\n")

lmm_leaves <- lmer(leaf_count ~ substrate * time_numeric + (1 | plant_id),
                    data = ag_growth,
                    REML = TRUE)
cat("  Leaf LMM summary:\n")
print(summary(lmm_leaves)$coefficients[1:5, ])

emm_leaves <- emmeans(lmm_leaves, "substrate")
pairs_leaves <- pairs(emm_leaves, adjust = "fdr")
cld_leaves <- multcomp::cld(emm_leaves, adjust = "fdr", Letters = letters)
cld_leaves <- as.data.frame(cld_leaves) %>%
  mutate(.group = trimws(.group))
cat("\n  CLD for leaf count:\n")
print(cld_leaves[, c("substrate", "emmean", ".group")])


# =====================================================
# LaTeX TABLE: LMM results
# =====================================================
cat("\n--- Writing LaTeX table ---\n")

sink(file.path(table_dir, "aboveground_lmm_results.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Estimated marginal means (EMMs) for aboveground performance metrics")
cat(" from linear mixed-effects models (LMM) with substrate and time as fixed effects")
cat(" and plant identity as random intercept. Groups sharing the same letter are not")
cat(" significantly different (pairwise contrasts, FDR-adjusted $p < 0.05$).}\n")
cat("\\label{tab:aboveground_lmm}\n")
cat("\\small\n")
cat("\\begin{tabular}{l rr cc}\n")
cat("\\toprule\n")
cat("Substrate & Height EMM (cm) & Leaf count EMM & Height group & Leaf group \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(cld_height))) {
  s <- as.character(cld_height$substrate[i])
  h_emm <- round(cld_height$emmean[i], 1)
  h_grp <- cld_height$.group[i]
  # Match substrate in leaf CLD
  l_row <- cld_leaves[as.character(cld_leaves$substrate) == s, ]
  l_emm <- round(l_row$emmean[1], 1)
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
# PLOTS
# =====================================================
cat("\n--- Plots ---\n")

# Compute means and SE per substrate × timepoint for ribbon plots
ag_summary <- ag_growth %>%
  group_by(substrate, time_numeric) %>%
  summarise(
    mean_height = mean(length, na.rm = TRUE),
    se_height   = sd(length, na.rm = TRUE) / sqrt(n()),
    mean_leaves = mean(leaf_count, na.rm = TRUE),
    se_leaves   = sd(leaf_count, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 1. Height LOESS over time (main text)
p_height_loess <- ggplot(ag_summary,
                         aes(x = time_numeric, y = mean_height,
                             colour = substrate, fill = substrate)) +
  geom_ribbon(aes(ymin = mean_height - se_height,
                  ymax = mean_height + se_height),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  scale_colour_substrate() +
  scale_fill_substrate() +
  scale_x_continuous(breaks = 1:13, labels = paste0("t", 1:13)) +
  labs(x = "Measurement period", y = "Mean leaf height (cm)",
       colour = "Substrate", fill = "Substrate") +
  theme_susana() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Leaf count LOESS over time (main text)
p_leaves_loess <- ggplot(ag_summary,
                          aes(x = time_numeric, y = mean_leaves,
                              colour = substrate, fill = substrate)) +
  geom_ribbon(aes(ymin = mean_leaves - se_leaves,
                  ymax = mean_leaves + se_leaves),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  scale_colour_substrate() +
  scale_fill_substrate() +
  scale_x_continuous(breaks = 1:13, labels = paste0("t", 1:13)) +
  labs(x = "Measurement period", y = "Mean leaf count",
       colour = "Substrate", fill = "Substrate") +
  theme_susana() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combined main text figure
p_ag_main <- p_height_loess + p_leaves_loess +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "bottom")
save_pub_plot(p_ag_main,
              file.path(out_dir, "aboveground_growth_main.png"),
              width = 260, height = 120)

# Also save individual plots for backward compatibility
save_pub_plot(p_height_loess,
              file.path(out_dir, "marram_height_loess.png"),
              width = 180, height = 120)
save_pub_plot(p_leaves_loess,
              file.path(out_dir, "leaves_number_loess.png"),
              width = 180, height = 120)

# 3. Detailed boxplots: Height per substrate per timepoint (appendix)
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

# 4. Detailed boxplots: Leaf count per substrate per timepoint (appendix)
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
