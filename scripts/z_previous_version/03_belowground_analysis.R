# =============================================================================
# 03_belowground_analysis.R
# Common garden belowground biomass: root distribution, root:shoot ratio
# Uses Kruskal-Wallis + Dunn's test (small n per group, non-normal)
# =============================================================================

cat("\n======================================================\n")
cat("  03_belowground_analysis.R\n")
cat("======================================================\n\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(FSA)
  library(patchwork)
  library(scales)
  library(emmeans)
})

source("d:/projects/susana/scripts/shared_theme.R")

# --- Load clean data ---
bg <- readRDS("d:/projects/susana/output/clean_data/belowground.rds")

out_dir <- "d:/projects/susana/output/common_garden/belowground_biomass"
rs_dir  <- "d:/projects/susana/output/common_garden/root_shoot_ratio"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rs_dir,  recursive = TRUE, showWarnings = FALSE)
table_dir <- "d:/projects/susana/evaluation/manuscript/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# =====================================================
# DATA PREPARATION
# =====================================================
cat("--- Data summary ---\n")
cat("  Total rows:", nrow(bg), "\n")
cat("  Unique plants:", length(unique(bg$plant_id)), "\n")
cat("  Layers:", paste(unique(as.character(bg$layer)), collapse = ", "), "\n")

# --- Total root biomass per plant (sum across layers, excluding shoots) ---
root_total <- bg %>%
  filter(layer != "shoots", !is.na(root_biomass)) %>%
  group_by(plant_id, substrate) %>%
  summarise(
    total_root_biomass = sum(root_biomass, na.rm = TRUE),
    n_layers           = n(),
    .groups = "drop"
  )

# --- Shoot biomass per plant ---
shoot_total <- bg %>%
  filter(layer == "shoots" | (layer == "top" & !is.na(shoot_biomass))) %>%
  group_by(plant_id, substrate) %>%
  summarise(
    total_shoot_biomass = sum(shoot_biomass, na.rm = TRUE),
    .groups = "drop"
  )

# --- Merge for root:shoot ratio ---
rs_df <- root_total %>%
  left_join(shoot_total, by = c("plant_id", "substrate")) %>%
  filter(!is.na(total_shoot_biomass), total_shoot_biomass > 0) %>%
  mutate(root_shoot_ratio = total_root_biomass / total_shoot_biomass)

cat("  Plants with root data:", nrow(root_total), "\n")
cat("  Plants with root:shoot ratio:", nrow(rs_df), "\n")
cat("  Root:shoot range:",
    round(min(rs_df$root_shoot_ratio, na.rm = TRUE), 2), "-",
    round(max(rs_df$root_shoot_ratio, na.rm = TRUE), 2), "\n")

# --- Root biomass per layer (depth profile) ---
root_layers <- bg %>%
  filter(layer != "shoots", !is.na(root_biomass)) %>%
  mutate(
    depth_label = factor(layer,
                         levels = c("top", "1", "2", "3", "4"),
                         labels = c("0-6 cm", "6-16 cm", "16-26 cm",
                                    "26-36 cm", "36-46 cm"),
                         ordered = TRUE)
  ) %>%
  filter(!is.na(depth_label))


# =====================================================
# STATISTICS
# =====================================================
cat("\n--- Statistics ---\n")

# --- Kruskal-Wallis: total root biomass ---
kw_root <- kruskal.test(total_root_biomass ~ substrate, data = root_total)
cat("  KW total root biomass: chi2 =", round(kw_root$statistic, 2),
    "  df =", kw_root$parameter, "  p =",
    format.pval(kw_root$p.value, digits = 3), "\n")

# Dunn's post-hoc if significant
if (kw_root$p.value < 0.05) {
  dunn_root <- dunnTest(total_root_biomass ~ substrate,
                        data = root_total, method = "bh")
  cat("  Dunn's test (BH-adjusted) significant pairs:\n")
  sig_pairs <- dunn_root$res %>% filter(P.adj < 0.05)
  if (nrow(sig_pairs) > 0) {
    print(sig_pairs[, c("Comparison", "Z", "P.adj")])
  } else {
    cat("    None (all FDR-adjusted p > 0.05)\n")
  }
} else {
  dunn_root <- NULL
  cat("  No post-hoc needed (KW p >= 0.05)\n")
}

# --- Kruskal-Wallis: root:shoot ratio ---
kw_rs <- kruskal.test(root_shoot_ratio ~ substrate, data = rs_df)
cat("  KW root:shoot ratio: chi2 =", round(kw_rs$statistic, 2),
    "  df =", kw_rs$parameter, "  p =",
    format.pval(kw_rs$p.value, digits = 3), "\n")

if (kw_rs$p.value < 0.05) {
  dunn_rs <- dunnTest(root_shoot_ratio ~ substrate,
                      data = rs_df, method = "bh")
  cat("  Dunn's test (BH-adjusted) significant pairs:\n")
  sig_pairs_rs <- dunn_rs$res %>% filter(P.adj < 0.05)
  if (nrow(sig_pairs_rs) > 0) {
    print(sig_pairs_rs[, c("Comparison", "Z", "P.adj")])
  } else {
    cat("    None (all FDR-adjusted p > 0.05)\n")
  }
} else {
  dunn_rs <- NULL
  cat("  No post-hoc needed (KW p >= 0.05)\n")
}

# --- Kruskal-Wallis: shoot biomass ---
kw_shoot <- kruskal.test(total_shoot_biomass ~ substrate,
                         data = shoot_total %>% filter(!is.na(total_shoot_biomass)))
cat("  KW shoot biomass: chi2 =", round(kw_shoot$statistic, 2),
    "  df =", kw_shoot$parameter, "  p =",
    format.pval(kw_shoot$p.value, digits = 3), "\n")

# --- Summary statistics per substrate ---
bg_summary <- rs_df %>%
  group_by(substrate) %>%
  summarise(
    n = n(),
    root_mean    = mean(total_root_biomass, na.rm = TRUE),
    root_sd      = sd(total_root_biomass, na.rm = TRUE),
    shoot_mean   = mean(total_shoot_biomass, na.rm = TRUE),
    shoot_sd     = sd(total_shoot_biomass, na.rm = TRUE),
    rs_mean      = mean(root_shoot_ratio, na.rm = TRUE),
    rs_sd        = sd(root_shoot_ratio, na.rm = TRUE),
    .groups = "drop"
  )
cat("\n  Summary per substrate:\n")
print(as.data.frame(bg_summary))


# =====================================================
# LaTeX TABLES
# =====================================================
cat("\n--- Writing LaTeX tables ---\n")

# Biomass summary table
sink(file.path(table_dir, "belowground_biomass_summary.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Belowground biomass and root:shoot ratio summary per substrate.")
cat(" Values are mean $\\pm$ SD. Kruskal-Wallis tests:")
cat(sprintf(" total root biomass $\\chi^2(%d)$ = %.2f, $p$ = %s;",
            kw_root$parameter, kw_root$statistic,
            format.pval(kw_root$p.value, digits = 3)))
cat(sprintf(" root:shoot ratio $\\chi^2(%d)$ = %.2f, $p$ = %s.}\n",
            kw_rs$parameter, kw_rs$statistic,
            format.pval(kw_rs$p.value, digits = 3)))
cat("\\label{tab:belowground_summary}\n")
cat("\\small\n")
cat("\\begin{tabular}{l r rr rr rr}\n")
cat("\\toprule\n")
cat("Substrate & $n$ & Root (g) & & Shoot (g) & & R:S ratio & \\\\\n")
cat("          &     & Mean & SD & Mean & SD & Mean & SD \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(bg_summary))) {
  r <- bg_summary[i, ]
  cat(sprintf("%-8s & %d & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n",
              as.character(r$substrate), r$n,
              r$root_mean, r$root_sd,
              r$shoot_mean, r$shoot_sd,
              r$rs_mean, r$rs_sd))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("  Table written: belowground_biomass_summary.tex\n")

# Statistics table
sink(file.path(table_dir, "belowground_statistics.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Non-parametric statistical tests for belowground biomass metrics")
cat(" across the SOM gradient. Kruskal-Wallis test with Dunn's post-hoc")
cat(" (Benjamini-Hochberg correction).}\n")
cat("\\label{tab:belowground_stats}\n")
cat("\\small\n")
cat("\\begin{tabular}{l r r r}\n")
cat("\\toprule\n")
cat("Response variable & $\\chi^2$ & df & $p$-value \\\\\n")
cat("\\midrule\n")
cat(sprintf("Total root biomass & %.2f & %d & %s \\\\\n",
            kw_root$statistic, kw_root$parameter,
            format.pval(kw_root$p.value, digits = 3)))
cat(sprintf("Root:shoot ratio & %.2f & %d & %s \\\\\n",
            kw_rs$statistic, kw_rs$parameter,
            format.pval(kw_rs$p.value, digits = 3)))
cat(sprintf("Aboveground biomass & %.2f & %d & %s \\\\\n",
            kw_shoot$statistic, kw_shoot$parameter,
            format.pval(kw_shoot$p.value, digits = 3)))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("  Table written: belowground_statistics.tex\n")


# =====================================================
# PLOTS
# =====================================================
cat("\n--- Plots ---\n")

# 1. Total root biomass per substrate (main text)
p_root_total <- ggplot(root_total,
                       aes(x = substrate, y = total_root_biomass,
                           fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, shape = 16) +
  scale_fill_substrate(guide = "none") +
  labs(x = "Substrate", y = "Total root biomass (g)") +
  theme_susana()

# 2. Root:shoot ratio per substrate (main text)
p_rs <- ggplot(rs_df,
               aes(x = substrate, y = root_shoot_ratio,
                   fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, shape = 16) +
  scale_fill_substrate(guide = "none") +
  labs(x = "Substrate", y = "Root:shoot ratio") +
  theme_susana()

# Combined main text figure
p_bg_main <- p_root_total + p_rs +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))
save_pub_plot(p_bg_main,
              file.path(out_dir, "belowground_main.png"),
              width = 260, height = 120)

# 3. Root depth profile - mean biomass per layer per substrate (appendix)
depth_summary <- root_layers %>%
  group_by(substrate, depth_label) %>%
  summarise(
    mean_biomass = mean(root_biomass, na.rm = TRUE),
    se_biomass   = sd(root_biomass, na.rm = TRUE) / sqrt(n()),
    n            = n(),
    .groups      = "drop"
  )

p_depth <- ggplot(depth_summary,
                  aes(x = depth_label, y = mean_biomass,
                      fill = substrate)) +
  geom_col(position = position_dodge(width = 0.8),
           alpha = 0.7, width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - se_biomass,
                    ymax = mean_biomass + se_biomass),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.3) +
  scale_fill_substrate() +
  labs(x = "Soil depth", y = "Mean root biomass (g)",
       fill = "Substrate") +
  theme_susana() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
save_pub_plot(p_depth,
              file.path(out_dir, "root_depth_profile.png"),
              width = 200, height = 130)

# 4. Shoot biomass vs root biomass scatter (appendix)
p_scatter <- ggplot(rs_df,
                    aes(x = total_shoot_biomass, y = total_root_biomass,
                        colour = substrate)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = TRUE, alpha = 0.1, linewidth = 0.5,
              colour = "grey30") +
  scale_colour_substrate() +
  labs(x = "Shoot biomass (g)", y = "Root biomass (g)",
       colour = "Substrate") +
  theme_susana() +
  theme(legend.position = "bottom")
save_pub_plot(p_scatter,
              file.path(out_dir, "root_vs_shoot_scatter.png"),
              width = 160, height = 130)

# 5. Height vs shoot biomass (was in old student output)
# Load aboveground data for cross-reference
ag <- readRDS("d:/projects/susana/output/clean_data/aboveground_planted.rds")

# Get final height per plant
ag_final <- ag %>%
  filter(!is.na(length), length > 0) %>%
  group_by(plant_id, substrate) %>%
  summarise(final_height = last(length), .groups = "drop")

# Merge with shoot biomass
height_shoot <- ag_final %>%
  inner_join(shoot_total, by = c("plant_id", "substrate")) %>%
  filter(!is.na(total_shoot_biomass))

if (nrow(height_shoot) > 0) {
  cor_hs <- cor.test(height_shoot$final_height,
                     height_shoot$total_shoot_biomass,
                     method = "pearson")
  cat("  Height vs shoot biomass: r =", round(cor_hs$estimate, 3),
      "  p =", format.pval(cor_hs$p.value, digits = 3), "\n")

  p_hs <- ggplot(height_shoot,
                 aes(x = final_height, y = total_shoot_biomass,
                     colour = substrate)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x,
                se = TRUE, alpha = 0.1, colour = "grey30",
                linewidth = 0.5) +
    scale_colour_substrate() +
    labs(x = "Final leaf height (cm)", y = "Shoot biomass (g)",
         colour = "Substrate") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f, p = %s",
                             cor_hs$estimate,
                             format.pval(cor_hs$p.value, digits = 2)),
             hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic") +
    theme_susana() +
    theme(legend.position = "bottom")
  save_pub_plot(p_hs,
                file.path(out_dir, "height_vs_shoot_biomass.png"),
                width = 160, height = 130)
}


cat("\n======================================================\n")
cat("  03_belowground_analysis.R COMPLETE\n")
cat("======================================================\n")
