# =============================================================================
# 04_germination_analysis.R
# Germination experiment: seed germination and seedling growth over time
# Uses LMMs for repeated measures, microclimate logger analysis
# =============================================================================

cat("\n======================================================\n")
cat("  04_germination_analysis.R\n")
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
germ <- readRDS("d:/projects/susana/output/clean_data/germination.rds")
loggers <- readRDS("d:/projects/susana/output/clean_data/loggers.rds")

out_dir_germ   <- "d:/projects/susana/output/germination"
out_dir_logger <- "d:/projects/susana/output/germination"
dir.create(out_dir_germ,   recursive = TRUE, showWarnings = FALSE)
table_dir <- "d:/projects/susana/evaluation/manuscript/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Germination substrate palette (subset of 6)
GERM_LEVELS <- c("R", "W", "D25", "D50", "D75", "D")
germ_palette <- substrate_palette[GERM_LEVELS]


# =====================================================
# DATA PREPARATION
# =====================================================
cat("--- Data summary ---\n")
cat("  Total observations:", nrow(germ), "\n")
cat("  Unique pots:", length(unique(germ$pot_id)), "\n")
cat("  Substrates:", paste(GERM_LEVELS, collapse = ", "), "\n")
cat("  Timepoints:", length(unique(germ$timepoint)), "\n")

# Filter to timepoints with actual data (t1, t2, t6, t11)
valid_tp <- c("t1", "t2", "t6", "t11")
germ_valid <- germ %>%
  filter(timepoint %in% valid_tp, !is.na(leaf_count)) %>%
  mutate(timepoint = factor(timepoint, levels = valid_tp, ordered = TRUE))

cat("  Observations with data:", nrow(germ_valid), "\n")

# Initial germination counts (t1): 10 seeds per pot
germ_t1 <- germ_valid %>%
  filter(timepoint == "t1") %>%
  mutate(pct_germinated = (leaf_count / 10) * 100)

cat("  Mean germination at t1:", round(mean(germ_t1$pct_germinated), 1), "%\n")


# =====================================================
# STATISTICS: GERMINATION
# =====================================================
cat("\n--- Statistics: Germination (t1) ---\n")

# Kruskal-Wallis for germination percentage (non-normal count data)
kw_germ <- kruskal.test(pct_germinated ~ substrate, data = germ_t1)
cat("  KW germination rate: chi2 =", round(kw_germ$statistic, 2),
    "  df =", kw_germ$parameter,
    "  p =", format.pval(kw_germ$p.value, digits = 3), "\n")

if (kw_germ$p.value < 0.05) {
  # dunnTest requires unordered factor
  germ_t1_dunn <- germ_t1 %>% mutate(substrate = factor(substrate, ordered = FALSE))
  dunn_germ <- FSA::dunnTest(pct_germinated ~ substrate,
                             data = germ_t1_dunn, method = "bh")
  cat("  Dunn's test significant pairs (FDR < 0.05):\n")
  sig_g <- dunn_germ$res %>% filter(P.adj < 0.05)
  if (nrow(sig_g) > 0) {
    print(sig_g[, c("Comparison", "Z", "P.adj")])
  } else {
    cat("    None\n")
  }
} else {
  dunn_germ <- NULL
  cat("  No post-hoc needed (KW p >= 0.05)\n")
}

# Germination summary per substrate
germ_summary <- germ_t1 %>%
  group_by(substrate) %>%
  summarise(
    n = n(),
    mean_pct = mean(pct_germinated, na.rm = TRUE),
    sd_pct   = sd(pct_germinated, na.rm = TRUE),
    .groups  = "drop"
  )
cat("\n  Germination % at t1 per substrate:\n")
print(as.data.frame(germ_summary))


# =====================================================
# STATISTICS: SEEDLING GROWTH LMM
# =====================================================
cat("\n--- Statistics: Seedling count over time (LMM) ---\n")

# LMM: leaf count ~ substrate * time + (1 | pot_id)
lmm_germ <- lmer(leaf_count ~ substrate * time_numeric + (1 | pot_id),
                  data = germ_valid,
                  REML = TRUE)

# EMMs for substrate effect
emm_germ <- emmeans(lmm_germ, "substrate")
cat("  EMMs for seedling count by substrate:\n")
print(summary(emm_germ))

# Pairwise comparisons
pairs_germ <- pairs(emm_germ, adjust = "fdr")
sig_germ_lmm <- as.data.frame(pairs_germ) %>% filter(p.value < 0.05)
cat("\n  Significant pairwise seedling count differences (FDR < 0.05):\n")
if (nrow(sig_germ_lmm) > 0) {
  print(sig_germ_lmm[, c("contrast", "estimate", "p.value")])
} else {
  cat("  None\n")
}

# CLD for plots
cld_germ <- multcomp::cld(emm_germ, adjust = "fdr", Letters = letters)
cld_germ <- as.data.frame(cld_germ) %>%
  mutate(.group = trimws(.group))
cat("\n  CLD for seedling count:\n")
print(cld_germ[, c("substrate", "emmean", ".group")])


# =====================================================
# STATISTICS: MICROCLIMATE
# =====================================================
cat("\n--- Statistics: Microclimate ---\n")

# Daily summaries
logger_daily <- loggers %>%
  filter(!is.na(T1_soil), !is.na(Vol_moisture)) %>%
  group_by(substrate, date) %>%
  summarise(
    T_soil_mean = mean(T1_soil, na.rm = TRUE),
    T_soil_max  = max(T1_soil, na.rm = TRUE),
    T_soil_min  = min(T1_soil, na.rm = TRUE),
    T_air_mean  = mean(T3_air, na.rm = TRUE),
    moisture_mean = mean(Vol_moisture, na.rm = TRUE),
    .groups = "drop"
  )

# Kruskal-Wallis for soil temperature
kw_temp <- kruskal.test(T_soil_mean ~ substrate, data = logger_daily)
cat("  KW mean soil temp: chi2 =", round(kw_temp$statistic, 2),
    "  p =", format.pval(kw_temp$p.value, digits = 3), "\n")

# Kruskal-Wallis for moisture
kw_moist <- kruskal.test(moisture_mean ~ substrate, data = logger_daily)
cat("  KW mean moisture: chi2 =", round(kw_moist$statistic, 2),
    "  p =", format.pval(kw_moist$p.value, digits = 3), "\n")

# Microclimate summary
micro_summary <- logger_daily %>%
  group_by(substrate) %>%
  summarise(
    n_days = n(),
    T_soil = sprintf("%.1f +/- %.1f",
                     mean(T_soil_mean), sd(T_soil_mean)),
    T_air  = sprintf("%.1f +/- %.1f",
                     mean(T_air_mean), sd(T_air_mean)),
    Moisture = sprintf("%.1f +/- %.1f",
                       mean(moisture_mean), sd(moisture_mean)),
    .groups = "drop"
  )
cat("\n  Microclimate summary:\n")
print(as.data.frame(micro_summary))


# =====================================================
# LaTeX TABLES
# =====================================================
cat("\n--- Writing LaTeX tables ---\n")

# Germination statistics table
sink(file.path(table_dir, "germination_statistics.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Germination and seedling establishment across substrates.")
cat(" Initial germination rate at $t_1$ (10 seeds per pot, $n$ = 16 pots per substrate)")
cat(" tested with Kruskal-Wallis;")
cat(" seedling count over time tested with LMM (substrate $\\times$ time, random intercept per pot).")
cat(sprintf(" Germination KW: $\\chi^2(%d)$ = %.2f, $p$ = %s.}\n",
            kw_germ$parameter, kw_germ$statistic,
            format.pval(kw_germ$p.value, digits = 3)))
cat("\\label{tab:germination_stats}\n")
cat("\\small\n")
cat("\\begin{tabular}{l r rr rr c}\n")
cat("\\toprule\n")
cat("Substrate & $n$ & Germ. \\% & SD & EMM count & SE & Group \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(germ_summary))) {
  s <- as.character(germ_summary$substrate[i])
  g <- germ_summary[i, ]
  c_row <- cld_germ[as.character(cld_germ$substrate) == s, ]
  cat(sprintf("%-6s & %d & %.1f & %.1f & %.1f & %.1f & %s \\\\\n",
              s, g$n, g$mean_pct, g$sd_pct,
              c_row$emmean[1], c_row$SE[1], c_row$.group[1]))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("  Table written: germination_statistics.tex\n")

# Microclimate table
sink(file.path(table_dir, "microclimate_summary.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Microclimate conditions in germination experiment substrates,")
cat(" measured by TMS-4 loggers. Values are daily means $\\pm$ SD.}\n")
cat("\\label{tab:microclimate}\n")
cat("\\small\n")
cat("\\begin{tabular}{l r l l l}\n")
cat("\\toprule\n")
cat("Substrate & Days & Soil temp. ($^{\\circ}$C) & Air temp. ($^{\\circ}$C) &")
cat(" Vol. moisture (\\%) \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(micro_summary))) {
  r <- micro_summary[i, ]
  cat(sprintf("%-6s & %d & %s & %s & %s \\\\\n",
              as.character(r$substrate), r$n_days,
              r$T_soil, r$T_air, r$Moisture))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("  Table written: microclimate_summary.tex\n")


# =====================================================
# PLOTS
# =====================================================
cat("\n--- Plots ---\n")

# 1. Germination rate at t1 (main text)
p_germ_bar <- ggplot(germ_summary,
                     aes(x = substrate, y = mean_pct, fill = substrate)) +
  geom_col(alpha = 0.7, width = 0.7) +
  geom_errorbar(aes(ymin = mean_pct - sd_pct,
                    ymax = mean_pct + sd_pct),
                width = 0.2, linewidth = 0.3) +
  geom_jitter(data = germ_t1,
              aes(x = substrate, y = pct_germinated),
              width = 0.1, alpha = 0.4, size = 1.5,
              shape = 16, inherit.aes = FALSE) +
  scale_fill_manual(values = germ_palette, guide = "none") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
  labs(x = "Substrate", y = "Germination (%)") +
  theme_susana()

# 2. Seedling count over time (main text)
germ_time_summary <- germ_valid %>%
  group_by(substrate, timepoint, time_numeric) %>%
  summarise(
    mean_count = mean(leaf_count, na.rm = TRUE),
    se_count   = sd(leaf_count, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

p_germ_time <- ggplot(germ_time_summary,
                      aes(x = time_numeric, y = mean_count,
                          colour = substrate, fill = substrate)) +
  geom_ribbon(aes(ymin = mean_count - se_count,
                  ymax = mean_count + se_count),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = germ_palette) +
  scale_fill_manual(values = germ_palette) +
  scale_x_continuous(breaks = c(1, 2, 6, 11),
                     labels = c("t1", "t2", "t6", "t11")) +
  labs(x = "Measurement period", y = "Mean seedling count",
       colour = "Substrate", fill = "Substrate") +
  theme_susana()

# Combined germination figure
p_germ_main <- p_germ_bar + p_germ_time +
  plot_layout(ncol = 2, widths = c(1, 1.3),
              guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "bottom")
save_pub_plot(p_germ_main,
              file.path(out_dir_germ, "germination_main.png"),
              width = 260, height = 120)

# 3. Soil temperature time series (appendix)
p_temp <- ggplot(logger_daily,
                 aes(x = date, y = T_soil_mean,
                     colour = substrate)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  scale_colour_manual(values = germ_palette) +
  labs(x = "Date", y = expression("Mean soil temperature ("*degree*"C)"),
       colour = "Substrate") +
  theme_susana() +
  theme(legend.position = "bottom")
save_pub_plot(p_temp,
              file.path(out_dir_logger, "soil_temperature_timeseries.pdf"),
              width = 200, height = 100)

# 4. Soil moisture time series (appendix)
p_moisture <- ggplot(logger_daily,
                     aes(x = date, y = moisture_mean,
                         colour = substrate)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  scale_colour_manual(values = germ_palette) +
  labs(x = "Date", y = "Volumetric moisture (%)",
       colour = "Substrate") +
  theme_susana() +
  theme(legend.position = "bottom")
save_pub_plot(p_moisture,
              file.path(out_dir_logger, "soil_moisture_timeseries.pdf"),
              width = 200, height = 100)

# 5. Temperature + moisture boxplots per substrate (appendix)
p_temp_box <- ggplot(logger_daily,
                     aes(x = substrate, y = T_soil_mean,
                         fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_manual(values = germ_palette, guide = "none") +
  labs(x = "Substrate",
       y = expression("Daily mean soil temp. ("*degree*"C)")) +
  theme_susana()

p_moist_box <- ggplot(logger_daily,
                      aes(x = substrate, y = moisture_mean,
                          fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_manual(values = germ_palette, guide = "none") +
  labs(x = "Substrate", y = "Daily mean moisture (%)") +
  theme_susana()

p_micro_box <- p_temp_box + p_moist_box +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))
save_pub_plot(p_micro_box,
              file.path(out_dir_logger, "microclimate_boxplots.pdf"),
              width = 260, height = 110)

# 6. Diurnal temperature pattern (appendix)
logger_diurnal <- loggers %>%
  filter(!is.na(T1_soil)) %>%
  group_by(substrate, hour) %>%
  summarise(
    T_mean = mean(T1_soil, na.rm = TRUE),
    T_se   = sd(T1_soil, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_diurnal <- ggplot(logger_diurnal,
                    aes(x = hour, y = T_mean,
                        colour = substrate, fill = substrate)) +
  geom_ribbon(aes(ymin = T_mean - T_se, ymax = T_mean + T_se),
              alpha = 0.1, colour = NA) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = germ_palette) +
  scale_fill_manual(values = germ_palette) +
  scale_x_continuous(breaks = seq(0, 23, 4)) +
  labs(x = "Hour of day", y = expression("Soil temperature ("*degree*"C)"),
       colour = "Substrate", fill = "Substrate") +
  theme_susana() +
  theme(legend.position = "bottom")
save_pub_plot(p_diurnal,
              file.path(out_dir_logger, "diurnal_temperature.pdf"),
              width = 180, height = 110)


cat("\n======================================================\n")
cat("  04_germination_analysis.R COMPLETE\n")
cat("======================================================\n")
