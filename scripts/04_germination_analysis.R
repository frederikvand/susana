# =============================================================================
# 04_germination_analysis.R
# Germination experiment + microclimate logger analysis
#
# Statistical approach:
#   Germination at t1 — Binomial GLM (cbind successes/failures)
#   Seedling count     — Negative Binomial GLMM (glmmTMB) with poly(time, 2)
#   Microclimate       — Kruskal-Wallis (non-parametric, single logger/substrate)
# Visualisation:
#   Model-predicted trajectories via ggeffects
# =============================================================================

cat("\n======================================================\n")
cat("  04_germination_analysis.R\n")
cat("======================================================\n\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
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
germ    <- readRDS(here("output", "clean_data", "germination.rds"))
loggers <- readRDS(here("output", "clean_data", "loggers.rds"))

out_dir   <- here("output", "germination")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- here("evaluation", "manuscript", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Germination substrates (6 of 10)
GERM_LEVELS <- c("R", "W", "D25", "D50", "D75", "D")
germ_palette <- substrate_palette[GERM_LEVELS]


# =====================================================
# DATA PREPARATION
# =====================================================
cat("--- Data summary ---\n")
cat("  Total obs:", nrow(germ), "\n")
cat("  Pots:", length(unique(germ$pot_id)), "\n")

valid_tp <- c("t1", "t2", "t6", "t11")
germ_valid <- germ %>%
  filter(timepoint %in% valid_tp, !is.na(leaf_count)) %>%
  mutate(timepoint = factor(timepoint, levels = valid_tp))

cat("  Valid obs:", nrow(germ_valid), "\n")
cat("  Seedling count — mean:", round(mean(germ_valid$leaf_count), 2),
    " var:", round(var(germ_valid$leaf_count), 2), "\n")

# Germination at t1: 10 seeds per pot
germ_t1 <- germ_valid %>%
  filter(timepoint == "t1") %>%
  mutate(
    n_germinated = as.integer(leaf_count),
    n_not        = as.integer(10L - n_germinated),
    pct_germinated = (n_germinated / 10) * 100,
    substrate_f  = factor(substrate, ordered = FALSE)
  )

# Unordered factor for seedling models
germ_valid$substrate_f <- factor(germ_valid$substrate, ordered = FALSE)

cat("  Mean germination at t1:", round(mean(germ_t1$pct_germinated), 1), "%\n")


# =====================================================
# MODEL 1: Germination — Binomial GLM
# =====================================================
cat("\n--- Binomial GLM: germination ~ substrate ---\n")

glm_germ <- glm(
  cbind(n_germinated, n_not) ~ substrate_f,
  family = binomial(link = "logit"),
  data   = germ_t1
)

cat("  AIC:", round(AIC(glm_germ), 1), "\n")
cat("  Deviance:", round(deviance(glm_germ), 2),
    " (df =", df.residual(glm_germ), ")\n")

# Check overdispersion (residual deviance / df)
disp_ratio <- deviance(glm_germ) / df.residual(glm_germ)
cat("  Dispersion ratio:", round(disp_ratio, 2))
if (disp_ratio > 1.5) {
  cat(" — overdispersed, switching to quasibinomial\n")
  glm_germ <- glm(
    cbind(n_germinated, n_not) ~ substrate_f,
    family = quasibinomial(link = "logit"),
    data   = germ_t1
  )
} else {
  cat(" — acceptable\n")
}

# EMMs on probability scale → percentages
emm_germ <- emmeans(glm_germ, "substrate_f", type = "response")
cld_germ <- cld(emm_germ, adjust = "fdr", Letters = letters)
cld_germ <- as.data.frame(cld_germ) %>%
  mutate(.group = trimws(.group),
         pct    = prob * 100,
         pct_lo = asymp.LCL * 100,
         pct_hi = asymp.UCL * 100)
cat("\n  Binomial GLM EMMs (germination probability):\n")
print(cld_germ[, c("substrate_f", "prob", "SE", ".group")])

# Model-predicted germination for plot
pred_germ <- ggpredict(glm_germ, terms = "substrate_f")
pred_germ_df <- as.data.frame(pred_germ) %>%
  rename(substrate = x) %>%
  mutate(pct = predicted * 100,
         pct_lo = conf.low * 100,
         pct_hi = conf.high * 100)


# =====================================================
# MODEL 2: Seedling count — NB GLMM
# =====================================================
cat("\n--- NB GLMM: seedling_count ~ substrate * poly(time, 2) + (1|pot_id) ---\n")

nb_seedling <- glmmTMB(
  leaf_count ~ substrate_f * poly(time_numeric, 2) + (1 | pot_id),
  family = nbinom2,
  data   = germ_valid
)

# Compare with Poisson
pois_seedling <- glmmTMB(
  leaf_count ~ substrate_f * poly(time_numeric, 2) + (1 | pot_id),
  family = poisson,
  data   = germ_valid
)
cat("  Poisson AIC:", round(AIC(pois_seedling), 1),
    "  NB AIC:", round(AIC(nb_seedling), 1), "\n")
cat("  NB dispersion:", round(sigma(nb_seedling), 2), "\n")

# EMMs on response scale
emm_seedling <- emmeans(nb_seedling, "substrate_f", type = "response")
cld_seedling <- cld(emm_seedling, adjust = "fdr", Letters = letters)
cld_seedling <- as.data.frame(cld_seedling) %>%
  mutate(.group = trimws(.group))
cat("\n  CLD for seedling count (response scale):\n")
print(cld_seedling[, c("substrate_f", "response", "SE", ".group")])

# Model-predicted trajectories
time_vals_germ <- sort(unique(germ_valid$time_numeric))
pred_seedling <- ggpredict(nb_seedling,
                           terms = c("time_numeric [time_vals_germ]",
                                     "substrate_f"))
pred_seedling_df <- as.data.frame(pred_seedling) %>%
  rename(time_numeric = x, substrate = group)


# =====================================================
# MICROCLIMATE STATISTICS
# =====================================================
cat("\n--- Microclimate statistics ---\n")

logger_daily <- loggers %>%
  filter(!is.na(T1_soil), !is.na(Vol_moisture)) %>%
  group_by(substrate, date) %>%
  summarise(
    T_soil_mean   = mean(T1_soil, na.rm = TRUE),
    T_soil_max    = max(T1_soil, na.rm = TRUE),
    T_soil_min    = min(T1_soil, na.rm = TRUE),
    T_air_mean    = mean(T3_air, na.rm = TRUE),
    moisture_mean = mean(Vol_moisture, na.rm = TRUE),
    .groups = "drop"
  )

kw_temp  <- kruskal.test(T_soil_mean ~ substrate, data = logger_daily)
kw_moist <- kruskal.test(moisture_mean ~ substrate, data = logger_daily)
cat("  KW soil temp: chi2 =", round(kw_temp$statistic, 2),
    " p =", format.pval(kw_temp$p.value, digits = 3), "\n")
cat("  KW moisture: chi2 =", round(kw_moist$statistic, 2),
    " p =", format.pval(kw_moist$p.value, digits = 3), "\n")

micro_summary <- logger_daily %>%
  group_by(substrate) %>%
  summarise(
    n_days   = n(),
    T_soil   = sprintf("%.1f +/- %.1f", mean(T_soil_mean), sd(T_soil_mean)),
    T_air    = sprintf("%.1f +/- %.1f", mean(T_air_mean), sd(T_air_mean)),
    Moisture = sprintf("%.1f +/- %.1f", mean(moisture_mean), sd(moisture_mean)),
    .groups  = "drop"
  )
cat("\n  Microclimate summary:\n")
print(as.data.frame(micro_summary))


# =====================================================
# LaTeX TABLES
# =====================================================
cat("\n--- Writing LaTeX tables ---\n")

# --- Germination statistics ---
germ_summary <- germ_t1 %>%
  group_by(substrate) %>%
  summarise(n = n(),
            mean_pct = mean(pct_germinated),
            sd_pct   = sd(pct_germinated),
            .groups  = "drop")

sink(file.path(table_dir, "germination_statistics.tex"))
cat("\\begin{table}[htbp]\n\\centering\n")
cat("\\caption{Germination and seedling establishment across substrates.")
cat(" Germination at $t_1$ (10 seeds per pot, $n=16$)")
cat(" modelled with a binomial GLM")
if (disp_ratio > 1.5) cat(" (quasibinomial due to overdispersion)")
cat("; seedling count over time modelled with a Negative Binomial GLMM")
cat(" (\\texttt{glmmTMB}, quadratic time, random intercept per pot).")
cat(" Groups sharing the same letter are not significantly different")
cat(" (FDR-adjusted pairwise contrasts, $p < 0.05$).}\n")
cat("\\label{tab:germination_stats}\n\\small\n")
cat("\\begin{tabular}{l r rr rr cc}\n\\toprule\n")
cat("Substrate & $n$ & Germ. \\% & SD & EMM count & SE & Germ. group & Count group \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(germ_summary))) {
  s <- as.character(germ_summary$substrate[i])
  g <- germ_summary[i, ]
  c_germ <- cld_germ[as.character(cld_germ$substrate_f) == s, ]
  c_seed <- cld_seedling[as.character(cld_seedling$substrate_f) == s, ]
  cat(sprintf("%-6s & %d & %.1f & %.1f & %.1f & %.1f & %s & %s \\\\\n",
              s, g$n, g$mean_pct, g$sd_pct,
              c_seed$response[1], c_seed$SE[1],
              c_germ$.group[1], c_seed$.group[1]))
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
sink()
cat("  Table: germination_statistics.tex\n")

# --- Microclimate ---
sink(file.path(table_dir, "microclimate_summary.tex"))
cat("\\begin{table}[htbp]\n\\centering\n")
cat("\\caption{Microclimate conditions in germination experiment substrates,")
cat(" measured by TMS-4 loggers. Values are daily means $\\pm$ SD.}\n")
cat("\\label{tab:microclimate}\n\\small\n")
cat("\\begin{tabular}{l r l l l}\n\\toprule\n")
cat("Substrate & Days & Soil temp. ($^{\\circ}$C) & Air temp. ($^{\\circ}$C) &")
cat(" Vol. moisture (\\%) \\\\\n\\midrule\n")
for (i in seq_len(nrow(micro_summary))) {
  r <- micro_summary[i, ]
  cat(sprintf("%-6s & %d & %s & %s & %s \\\\\n",
              as.character(r$substrate), r$n_days,
              r$T_soil, r$T_air, r$Moisture))
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
sink()
cat("  Table: microclimate_summary.tex\n")


# =====================================================
# PLOTS
# =====================================================
cat("\n--- Plots ---\n")

# --- 1. Germination: model-predicted probabilities (main text) ---
p_germ_bar <- ggplot() +
  # Raw data points
  geom_jitter(data = germ_t1,
              aes(x = substrate, y = pct_germinated),
              width = 0.15, alpha = 0.35, size = 1.5, shape = 16,
              colour = "grey40") +
  # Model-predicted means + CI from binomial GLM
  geom_pointrange(data = cld_germ,
                  aes(x = substrate_f, y = pct, ymin = pct_lo, ymax = pct_hi,
                      colour = substrate_f),
                  size = 0.6, linewidth = 0.8, show.legend = FALSE) +
  # CLD letters
  geom_text(data = cld_germ,
            aes(x = substrate_f, y = pct_hi + 5, label = .group),
            fontface = "bold", size = 3.5) +
  scale_colour_manual(values = germ_palette) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
  labs(x = "Substrate", y = "Germination (%)  [binomial GLM]") +
  theme_susana()

# --- 2. Seedling count: NB model trajectories (main text) ---
p_germ_time <- ggplot() +
  # Raw data as faint jittered points
  geom_point(data = germ_valid,
             aes(x = time_numeric, y = leaf_count, colour = substrate),
             alpha = 0.1, size = 0.8, show.legend = FALSE,
             position = position_jitter(width = 0.15, height = 0.15)) +
  # Model CI ribbons
  geom_ribbon(data = pred_seedling_df,
              aes(x = time_numeric, ymin = conf.low, ymax = conf.high,
                  fill = substrate),
              alpha = 0.15) +
  # Model predicted lines
  geom_line(data = pred_seedling_df,
            aes(x = time_numeric, y = predicted, colour = substrate),
            linewidth = 0.8) +
  scale_colour_manual(values = germ_palette) +
  scale_fill_manual(values = germ_palette) +
  scale_x_continuous(breaks = time_vals_germ,
                     labels = paste0("t", time_vals_germ)) +
  labs(x = "Measurement period",
       y = "Predicted seedling count  [NB GLMM]",
       colour = "Substrate", fill = "Substrate") +
  theme_susana()

# --- Combined: stacked vertically ---
p_germ_main <- p_germ_bar / p_germ_time +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "bottom")
save_pub_plot(p_germ_main,
              file.path(out_dir, "germination_main.png"),
              width = 180, height = 240)


# --- 3. Soil temperature time series (appendix) ---
p_temp <- ggplot(logger_daily,
                 aes(x = date, y = T_soil_mean, colour = substrate)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  scale_colour_manual(values = germ_palette) +
  labs(x = "Date", y = expression("Mean soil temperature ("*degree*"C)"),
       colour = "Substrate") +
  theme_susana() +
  theme(legend.position = "bottom")
save_pub_plot(p_temp,
              file.path(out_dir, "soil_temperature_timeseries.pdf"),
              width = 200, height = 100)

# --- 4. Soil moisture time series (appendix) ---
p_moisture <- ggplot(logger_daily,
                     aes(x = date, y = moisture_mean, colour = substrate)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  scale_colour_manual(values = germ_palette) +
  labs(x = "Date", y = "Volumetric moisture (%)",
       colour = "Substrate") +
  theme_susana() +
  theme(legend.position = "bottom")
save_pub_plot(p_moisture,
              file.path(out_dir, "soil_moisture_timeseries.pdf"),
              width = 200, height = 100)

# --- 5. Microclimate boxplots (appendix, stacked) ---
p_temp_box <- ggplot(logger_daily,
                     aes(x = substrate, y = T_soil_mean, fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_manual(values = germ_palette, guide = "none") +
  labs(x = "Substrate",
       y = expression("Daily mean soil temp. ("*degree*"C)")) +
  theme_susana()

p_moist_box <- ggplot(logger_daily,
                      aes(x = substrate, y = moisture_mean, fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_manual(values = germ_palette, guide = "none") +
  labs(x = "Substrate", y = "Daily mean moisture (%)") +
  theme_susana()

p_micro_box <- p_temp_box / p_moist_box +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))
save_pub_plot(p_micro_box,
              file.path(out_dir, "microclimate_boxplots.pdf"),
              width = 180, height = 220)

# --- 6. Diurnal temperature (appendix) ---
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
              file.path(out_dir, "diurnal_temperature.pdf"),
              width = 180, height = 110)


cat("\n======================================================\n")
cat("  04_germination_analysis.R COMPLETE\n")
cat("======================================================\n")
