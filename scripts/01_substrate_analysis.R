# =============================================================================
# 01_substrate_analysis.R
# Substrate characterisation: TOC, TN, D50 across alternative vs control
# =============================================================================

cat("\n======================================================\n")
cat("  01_substrate_analysis.R\n")
cat("======================================================\n\n")

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(here)
})

source(here("scripts", "shared_theme.R"))
set.seed(2024)

# --- Load data directly (substrate_sampling is the main CN/GS dataset) ---
sub_raw <- read_excel(here("data", "substrate", "substrate_sampling.xlsx"))

# Clean substrate names: map Dutch variants to English
sub_df <- sub_raw %>%
  mutate(
    substrate_clean = case_when(
      grepl("^control",     substrate, ignore.case = TRUE) ~ "Control",
      grepl("^alternative", substrate, ignore.case = TRUE) ~ "Alternative",
      grepl("^alternatief", substrate, ignore.case = TRUE) ~ "Alternative",
      TRUE ~ NA_character_
    ),
    substrate = factor(substrate_clean,
                       levels = c("Control", "Alternative")),
    in_out    = case_when(
      in_out == "in"  ~ "Inside",
      in_out == "out" ~ "Outside",
      TRUE ~ NA_character_
    ),
    in_out    = factor(in_out, levels = c("Inside", "Outside")),
    pct_TOC   = as.numeric(TOC),
    pct_TN    = as.numeric(TN),
    d50       = as.numeric(d0.5)
  ) %>%
  filter(!is.na(substrate))

cat("  Samples:", nrow(sub_df), "\n")
cat("  Substrate:", paste(table(sub_df$substrate), collapse = " / "), "\n")
cat("  In/out:", paste(table(sub_df$in_out, useNA = "always"), collapse = " / "), "\n")

# Load 10-substrate gradient properties
calc_props <- readRDS(here("output", "clean_data", "substrate_properties.rds"))

out_dir <- here("output", "substrate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- here("evaluation", "manuscript", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Colour palette for substrate type
substrate_type_cols <- c("Alternative" = "#E69F00", "Control" = "#56B4E9")


# =====================================================
# STATISTICS
# =====================================================
cat("\n--- Statistics ---\n")

# Wilcoxon rank-sum tests: alternative vs control
w_toc <- wilcox.test(pct_TOC ~ substrate, data = sub_df)
w_tn  <- wilcox.test(pct_TN  ~ substrate, data = sub_df)
w_d50 <- wilcox.test(d50     ~ substrate, data = sub_df %>% filter(!is.na(d50)))
cat("  Wilcoxon TOC (Alt vs Ctrl): W =", round(w_toc$statistic, 1),
    "  p =", format.pval(w_toc$p.value, digits = 3), "\n")
cat("  Wilcoxon TN  (Alt vs Ctrl): W =", round(w_tn$statistic, 1),
    "  p =", format.pval(w_tn$p.value, digits = 3), "\n")
cat("  Wilcoxon D50 (Alt vs Ctrl): W =", round(w_d50$statistic, 1),
    "  p =", format.pval(w_d50$p.value, digits = 3), "\n")

# In vs out comparison (within alternative substrate)
sub_io <- sub_df %>% filter(!is.na(in_out), substrate == "Alternative")
w_toc_io <- wilcox.test(pct_TOC ~ in_out, data = sub_io)
w_tn_io  <- wilcox.test(pct_TN  ~ in_out, data = sub_io)
w_d50_io <- wilcox.test(d50     ~ in_out, data = sub_io %>% filter(!is.na(d50)))
cat("  Wilcoxon TOC in/out: p =", format.pval(w_toc_io$p.value, digits = 3), "\n")
cat("  Wilcoxon TN  in/out: p =", format.pval(w_tn_io$p.value, digits = 3), "\n")
cat("  Wilcoxon D50 in/out: p =", format.pval(w_d50_io$p.value, digits = 3), "\n")

# Summary statistics
summary_stats <- sub_df %>%
  group_by(substrate) %>%
  summarise(
    n = n(),
    TOC_mean = mean(pct_TOC, na.rm = TRUE),
    TOC_sd   = sd(pct_TOC, na.rm = TRUE),
    TN_mean  = mean(pct_TN, na.rm = TRUE),
    TN_sd    = sd(pct_TN, na.rm = TRUE),
    D50_mean = mean(d50, na.rm = TRUE),
    D50_sd   = sd(d50, na.rm = TRUE),
    .groups = "drop"
  )
cat("\n  Summary:\n")
print(as.data.frame(summary_stats))


# =====================================================
# STATISTICS TABLE (LaTeX)
# =====================================================
cat("\n--- Writing LaTeX tables ---\n")

sink(file.path(table_dir, "substrate_statistics.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Wilcoxon rank-sum test results comparing sediment properties")
cat(" between alternative (dredged) and control (beach) substrates,")
cat(" and between inside and outside sampling positions within the alternative substrate.}\n")
cat("\\label{tab:substrate_stats}\n")
cat("\\small\n")
cat("\\begin{tabular}{l l r r}\n")
cat("\\toprule\n")
cat("Comparison & Property & $W$ & $p$-value \\\\\n")
cat("\\midrule\n")
cat(sprintf("Alternative vs Control & TOC (\\%%) & %.0f & %s \\\\\n",
            w_toc$statistic, format.pval(w_toc$p.value, digits = 3)))
cat(sprintf("Alternative vs Control & TN (\\%%) & %.0f & %s \\\\\n",
            w_tn$statistic, format.pval(w_tn$p.value, digits = 3)))
cat(sprintf("Alternative vs Control & $D_{50}$ ($\\mu$m) & %.0f & %s \\\\\n",
            w_d50$statistic, format.pval(w_d50$p.value, digits = 3)))
cat("\\midrule\n")
cat(sprintf("Inside vs Outside & TOC (\\%%) & %.0f & %s \\\\\n",
            w_toc_io$statistic, format.pval(w_toc_io$p.value, digits = 3)))
cat(sprintf("Inside vs Outside & TN (\\%%) & %.0f & %s \\\\\n",
            w_tn_io$statistic, format.pval(w_tn_io$p.value, digits = 3)))
cat(sprintf("Inside vs Outside & $D_{50}$ ($\\mu$m) & %.0f & %s \\\\\n",
            w_d50_io$statistic, format.pval(w_d50_io$p.value, digits = 3)))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("  Table written: substrate_statistics.tex\n")

# 10-substrate gradient table
sink(file.path(table_dir, "substrate_summary.tex"))
cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Physical and chemical properties of the 10 experimental")
cat(" substrates used in the common garden experiment.}\n")
cat("\\label{tab:substrate_summary}\n")
cat("\\small\n")
cat("\\begin{tabular}{l r r r l}\n")
cat("\\toprule\n")
cat("Substrate & TOC (\\%) & TN (\\%) & $D_{50}$ ($\\mu$m) & Description \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(calc_props))) {
  row <- calc_props[i, ]
  desc <- ifelse(is.null(row$Description) || is.na(row$Description),
                 "", as.character(row$Description))
  cat(sprintf("%-12s & %.2f & %.3f & %.0f & %s \\\\\n",
              as.character(row$substrate_label),
              row$pct_TOC, row$pct_TN, row$d50_um, desc))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
sink()
cat("  Table written: substrate_summary.tex\n")


# =====================================================
# PLOTS
# =====================================================
cat("\n--- Plots ---\n")

# 1. TOC: alternative vs control (main text)
p_toc <- ggplot(sub_df, aes(x = substrate, y = pct_TOC, fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1.5, shape = 16) +
  scale_fill_manual(values = substrate_type_cols, guide = "none") +
  labs(x = "Sediment type", y = "TOC (%)") +
  annotate("text", x = 1.5, y = max(sub_df$pct_TOC, na.rm = TRUE) * 1.05,
           label = paste0("p = ", format.pval(w_toc$p.value, digits = 2)),
           size = 3, fontface = "italic") +
  theme_susana()

# 2. D50 (main text)
p_d50 <- ggplot(sub_df %>% filter(!is.na(d50)),
                aes(x = substrate, y = d50, fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1.5, shape = 16) +
  scale_fill_manual(values = substrate_type_cols, guide = "none") +
  labs(x = "Sediment type",
       y = expression(D[50]~"("*mu*"m)")) +
  annotate("text", x = 1.5,
           y = max(sub_df$d50, na.rm = TRUE) * 1.05,
           label = paste0("p = ", format.pval(w_d50$p.value, digits = 2)),
           size = 3, fontface = "italic") +
  theme_susana()

p_main <- p_toc + p_d50 +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))
save_pub_plot(p_main,
              file.path(out_dir, "substrate_toc_d50.png"),
              width = 180, height = 100)

# 3. In/Out comparison (appendix)
sub_io_long <- sub_df %>%
  filter(!is.na(in_out)) %>%
  pivot_longer(cols = c(pct_TOC, pct_TN),
               names_to = "property", values_to = "value") %>%
  mutate(property = recode(property,
                           pct_TOC = "TOC (%)", pct_TN = "TN (%)"))

p_in_out <- ggplot(sub_io_long,
                   aes(x = in_out, y = value, fill = substrate)) +
  geom_boxplot(alpha = 0.7, linewidth = 0.3,
               position = position_dodge(width = 0.75)) +
  geom_point(alpha = 0.4, size = 1, shape = 16,
             position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.75)) +
  scale_fill_manual(values = substrate_type_cols) +
  facet_wrap(~ property, scales = "free_y") +
  labs(x = "Sampling position", y = NULL, fill = "Sediment") +
  theme_susana()
save_pub_plot(p_in_out,
              file.path(out_dir, "substrate_in_out_comparison.png"),
              width = 200, height = 110)

# 4. TN plot (appendix)
p_tn <- ggplot(sub_df, aes(x = substrate, y = pct_TN, fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1.5, shape = 16) +
  scale_fill_manual(values = substrate_type_cols, guide = "none") +
  labs(x = "Sediment type", y = "TN (%)") +
  annotate("text", x = 1.5, y = max(sub_df$pct_TN, na.rm = TRUE) * 1.05,
           label = paste0("p = ", format.pval(w_tn$p.value, digits = 2)),
           size = 3, fontface = "italic") +
  theme_susana()
save_pub_plot(p_tn,
              file.path(out_dir, "substrate_tn.png"),
              width = 120, height = 100)

# 5. Soil moisture boxplots (from logger data)
loggers <- tryCatch(
  readRDS(here("output", "clean_data", "loggers.rds")),
  error = function(e) NULL
)
if (!is.null(loggers)) {
  p_moisture <- ggplot(loggers %>% filter(!is.na(Vol_moisture)),
                       aes(x = substrate, y = Vol_moisture,
                           fill = substrate)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, linewidth = 0.3) +
    scale_fill_viridis_d(option = "D", guide = "none") +
    labs(x = "Substrate", y = "Volumetric soil moisture (%)") +
    theme_susana()
  save_pub_plot(p_moisture,
                file.path(out_dir, "soil_moisture_boxplots.png"),
                width = 160, height = 100)
}


cat("\n======================================================\n")
cat("  01_substrate_analysis.R COMPLETE\n")
cat("======================================================\n")
