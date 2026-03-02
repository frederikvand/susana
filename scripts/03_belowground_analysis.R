# =============================================================================
# 03_belowground_analysis.R
# Common garden belowground biomass: root distribution, root:shoot ratio
#
# Statistical approach:
#   Primary  — Gamma GLM (log link) for positive continuous biomass data
#   Fallback — Kruskal-Wallis + Dunn (BH) when Gamma fails
#   Both reported for completeness
# =============================================================================

cat("\n======================================================\n")
cat("  03_belowground_analysis.R\n")
cat("======================================================\n\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(FSA)
  library(emmeans)
  library(patchwork)
  library(scales)
  library(multcomp)
  library(here)
})

source(here("scripts", "shared_theme.R"))
set.seed(2024)

# --- Load clean data ---
bg <- readRDS(here("output", "clean_data", "belowground.rds"))

out_dir <- here("output", "common_garden", "belowground_biomass")
rs_dir  <- here("output", "common_garden", "root_shoot_ratio")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rs_dir,  recursive = TRUE, showWarnings = FALSE)
table_dir <- here("evaluation", "manuscript", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# =====================================================
# DATA PREPARATION
# =====================================================
cat("--- Data summary ---\n")
cat("  Total rows:", nrow(bg), "\n")
cat("  Unique plants:", length(unique(bg$plant_id)), "\n")

# Total root biomass per plant
root_total <- bg %>%
  filter(layer != "shoots", !is.na(root_biomass)) %>%
  group_by(plant_id, substrate) %>%
  summarise(
    total_root_biomass = sum(root_biomass, na.rm = TRUE),
    n_layers           = n(),
    .groups = "drop"
  )

# Shoot biomass per plant
shoot_total <- bg %>%
  filter(layer == "shoots" | (layer == "top" & !is.na(shoot_biomass))) %>%
  group_by(plant_id, substrate) %>%
  summarise(
    total_shoot_biomass = sum(shoot_biomass, na.rm = TRUE),
    .groups = "drop"
  )

# Root:shoot ratio
rs_df <- root_total %>%
  left_join(shoot_total, by = c("plant_id", "substrate")) %>%
  filter(!is.na(total_shoot_biomass), total_shoot_biomass > 0) %>%
  mutate(root_shoot_ratio = total_root_biomass / total_shoot_biomass)

cat("  Plants with root data:", nrow(root_total), "\n")
cat("  Plants with R:S ratio:", nrow(rs_df), "\n")
cat("  R:S range:",
    round(min(rs_df$root_shoot_ratio, na.rm = TRUE), 2), "--",
    round(max(rs_df$root_shoot_ratio, na.rm = TRUE), 2), "\n")

# Root biomass per depth layer
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

# Unordered substrate factor for GLM
root_total$substrate_f <- factor(root_total$substrate, ordered = FALSE)
rs_df$substrate_f      <- factor(rs_df$substrate, ordered = FALSE)
shoot_total$substrate_f <- factor(shoot_total$substrate, ordered = FALSE)


# =====================================================
# STATISTICS: Gamma GLM + Kruskal-Wallis
# =====================================================
cat("\n--- Statistics ---\n")

# ------------------------------------------------------------------
# Helper: try Gamma GLM, report results, fall back to KW
# ------------------------------------------------------------------
fit_gamma_or_kw <- function(formula_gamma, formula_kw, data,
                            response_label) {
  results <- list()

  # --- KW (always computed for reference) ---
  kw <- kruskal.test(formula_kw, data = data)
  cat(sprintf("  KW %s: chi2 = %.2f, df = %d, p = %s\n",
              response_label, kw$statistic, kw$parameter,
              format.pval(kw$p.value, digits = 3)))
  results$kw <- kw

  # Dunn post-hoc if KW significant
  if (kw$p.value < 0.05) {
    data_uo <- data
    data_uo$substrate <- factor(data_uo$substrate, ordered = FALSE)
    dunn <- dunnTest(formula_kw, data = data_uo, method = "bh")
    sig <- dunn$res %>% filter(P.adj < 0.05)
    if (nrow(sig) > 0) {
      cat("    Dunn significant pairs:\n")
      print(sig[, c("Comparison", "Z", "P.adj")])
    } else {
      cat("    No significant Dunn pairs\n")
    }
    results$dunn <- dunn
  }

  # --- Gamma GLM (only on strictly positive data) ---
  resp_vals <- model.response(model.frame(formula_kw, data))
  data_pos <- data[resp_vals > 0, ]
  n_zero <- sum(resp_vals == 0)
  if (n_zero > 0) {
    cat(sprintf("    %d zero values excluded for Gamma GLM\n", n_zero))
  }

  gamma_ok <- tryCatch({
    gm <- glm(formula_gamma, family = Gamma(link = "log"), data = data_pos)
    cat(sprintf("  Gamma GLM %s: AIC = %.1f, deviance = %.2f\n",
                response_label, AIC(gm), deviance(gm)))

    # Likelihood-ratio test vs null
    gm_null <- glm(update(formula_gamma, . ~ 1),
                    family = Gamma(link = "log"), data = data_pos)
    lr_test <- anova(gm_null, gm, test = "F")
    cat(sprintf("    Gamma F-test vs null: F = %.2f, p = %s\n",
                lr_test$F[2],
                format.pval(lr_test$`Pr(>F)`[2], digits = 3)))

    # EMMs on response scale
    emm <- emmeans(gm, "substrate_f", type = "response")
    cld_df <- cld(emm, adjust = "fdr", Letters = letters)
    cld_df <- as.data.frame(cld_df) %>% mutate(.group = trimws(.group))
    cat(sprintf("    Gamma EMMs for %s:\n", response_label))
    print(cld_df[, c("substrate_f", "response", "SE", ".group")])

    results$gamma  <- gm
    results$emm    <- emm
    results$cld    <- cld_df
    results$method <- "Gamma GLM"
    TRUE
  }, error = function(e) {
    cat(sprintf("    Gamma GLM failed for %s: %s\n",
                response_label, conditionMessage(e)))
    FALSE
  })

  if (!gamma_ok) {
    results$method <- "Kruskal-Wallis"
  }

  results
}

# --- Total root biomass ---
cat("\n  [Total root biomass]\n")
res_root <- fit_gamma_or_kw(
  total_root_biomass ~ substrate_f,
  total_root_biomass ~ substrate,
  root_total, "total root biomass"
)

# --- Root:shoot ratio ---
cat("\n  [Root:shoot ratio]\n")
res_rs <- fit_gamma_or_kw(
  root_shoot_ratio ~ substrate_f,
  root_shoot_ratio ~ substrate,
  rs_df, "root:shoot ratio"
)

# --- Shoot biomass ---
cat("\n  [Shoot biomass]\n")
shoot_pos <- shoot_total %>% filter(!is.na(total_shoot_biomass))
res_shoot <- fit_gamma_or_kw(
  total_shoot_biomass ~ substrate_f,
  total_shoot_biomass ~ substrate,
  shoot_pos, "shoot biomass"
)

# --- Height vs shoot biomass correlation ---
ag <- readRDS(here("output", "clean_data", "aboveground_planted.rds"))
ag_final <- ag %>%
  filter(!is.na(length), length > 0) %>%
  group_by(plant_id, substrate) %>%
  summarise(final_height = last(length), .groups = "drop")
height_shoot <- ag_final %>%
  inner_join(shoot_total, by = c("plant_id", "substrate")) %>%
  filter(!is.na(total_shoot_biomass))
cor_hs <- NULL
if (nrow(height_shoot) > 0) {
  cor_hs <- cor.test(height_shoot$final_height,
                     height_shoot$total_shoot_biomass, method = "pearson")
  cat("\n  Height vs shoot biomass: r =", round(cor_hs$estimate, 3),
      "  p =", format.pval(cor_hs$p.value, digits = 3), "\n")
}

# --- Per-substrate summary ---
bg_summary <- rs_df %>%
  group_by(substrate) %>%
  summarise(
    n = n(),
    root_mean  = mean(total_root_biomass), root_sd = sd(total_root_biomass),
    shoot_mean = mean(total_shoot_biomass), shoot_sd = sd(total_shoot_biomass),
    rs_mean    = mean(root_shoot_ratio), rs_sd = sd(root_shoot_ratio),
    .groups = "drop"
  )
cat("\n  Summary per substrate:\n")
print(as.data.frame(bg_summary))


# =====================================================
# LaTeX TABLES
# =====================================================
cat("\n--- Writing LaTeX tables ---\n")

# --- Biomass summary ---
sink(file.path(table_dir, "belowground_biomass_summary.tex"))
cat("\\begin{table}[htbp]\n\\centering\n")
cat("\\caption{Belowground biomass and root:shoot ratio per substrate (mean $\\pm$ SD).")
cat(" Statistical tests: Gamma GLM (log link) where feasible;")
cat(" Kruskal-Wallis rank-sum otherwise. See \\cref{tab:belowground_stats} for test details.}\n")
cat("\\label{tab:belowground_summary}\n\\small\n")
cat("\\begin{tabular}{l r rr rr rr}\n\\toprule\n")
cat("Substrate & $n$ & Root (g) & & Shoot (g) & & R:S ratio & \\\\\n")
cat("          &     & Mean & SD & Mean & SD & Mean & SD \\\\\n\\midrule\n")
for (i in seq_len(nrow(bg_summary))) {
  r <- bg_summary[i, ]
  cat(sprintf("%-8s & %d & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n",
              as.character(r$substrate), r$n,
              r$root_mean, r$root_sd,
              r$shoot_mean, r$shoot_sd,
              r$rs_mean, r$rs_sd))
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
sink()
cat("  Table: belowground_biomass_summary.tex\n")

# --- Statistics table ---
sink(file.path(table_dir, "belowground_statistics.tex"))
cat("\\begin{table}[htbp]\n\\centering\n")
cat("\\caption{Statistical tests for belowground biomass metrics across the SOM gradient.")
cat(" Gamma GLM (log link) was fitted where data were strictly positive;")
cat(" Kruskal-Wallis is reported in parallel as a non-parametric reference.")
cat(" Benjamini-Hochberg correction applied to pairwise comparisons.}\n")
cat("\\label{tab:belowground_stats}\n\\small\n")
cat("\\begin{tabular}{l l r r r}\n\\toprule\n")
cat("Response & Method & Statistic & df & $p$-value \\\\\n\\midrule\n")
# Root
cat(sprintf("Total root biomass & KW & $\\chi^2$ = %.2f & %d & %s \\\\\n",
            res_root$kw$statistic, res_root$kw$parameter,
            format.pval(res_root$kw$p.value, digits = 3)))
if (!is.null(res_root$gamma)) {
  lr <- anova(glm(total_root_biomass ~ 1, Gamma(link="log"),
                  root_total[root_total$total_root_biomass>0,]),
              res_root$gamma, test = "F")
  cat(sprintf("                   & Gamma & F = %.2f & %d & %s \\\\\n",
              lr$F[2], lr$Df[2],
              format.pval(lr$`Pr(>F)`[2], digits = 3)))
}
# R:S ratio
cat(sprintf("Root:shoot ratio & KW & $\\chi^2$ = %.2f & %d & %s \\\\\n",
            res_rs$kw$statistic, res_rs$kw$parameter,
            format.pval(res_rs$kw$p.value, digits = 3)))
if (!is.null(res_rs$gamma)) {
  lr2 <- anova(glm(root_shoot_ratio ~ 1, Gamma(link="log"),
                   rs_df[rs_df$root_shoot_ratio>0,]),
               res_rs$gamma, test = "F")
  cat(sprintf("                 & Gamma & F = %.2f & %d & %s \\\\\n",
              lr2$F[2], lr2$Df[2],
              format.pval(lr2$`Pr(>F)`[2], digits = 3)))
}
# Shoot
cat(sprintf("Shoot biomass & KW & $\\chi^2$ = %.2f & %d & %s \\\\\n",
            res_shoot$kw$statistic, res_shoot$kw$parameter,
            format.pval(res_shoot$kw$p.value, digits = 3)))
if (!is.null(res_shoot$gamma)) {
  lr3 <- anova(glm(total_shoot_biomass ~ 1, Gamma(link="log"),
                   shoot_pos[shoot_pos$total_shoot_biomass>0,]),
               res_shoot$gamma, test = "F")
  cat(sprintf("              & Gamma & F = %.2f & %d & %s \\\\\n",
              lr3$F[2], lr3$Df[2],
              format.pval(lr3$`Pr(>F)`[2], digits = 3)))
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
sink()
cat("  Table: belowground_statistics.tex\n")


# =====================================================
# PLOTS (stacked vertically)
# =====================================================
cat("\n--- Plots ---\n")

# 1. Total root biomass boxplot
p_root_total <- ggplot(root_total,
                       aes(x = substrate, y = total_root_biomass,
                           fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, shape = 16) +
  scale_fill_substrate(guide = "none") +
  labs(x = "Substrate", y = "Total root biomass (g)") +
  theme_susana()

# 2. Root:shoot ratio boxplot
p_rs <- ggplot(rs_df,
               aes(x = substrate, y = root_shoot_ratio, fill = substrate)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, shape = 16) +
  scale_fill_substrate(guide = "none") +
  labs(x = "Substrate", y = "Root:shoot ratio") +
  theme_susana()

# Combined: stacked vertically
p_bg_main <- p_root_total / p_rs +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))
save_pub_plot(p_bg_main,
              file.path(out_dir, "belowground_main.png"),
              width = 180, height = 240)

# 3. Root depth profile
depth_summary <- root_layers %>%
  group_by(substrate, depth_label) %>%
  summarise(
    mean_biomass = mean(root_biomass, na.rm = TRUE),
    se_biomass   = sd(root_biomass, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_depth <- ggplot(depth_summary,
                  aes(x = depth_label, y = mean_biomass, fill = substrate)) +
  geom_col(position = position_dodge(0.8), alpha = 0.7, width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - se_biomass,
                    ymax = mean_biomass + se_biomass),
                position = position_dodge(0.8),
                width = 0.2, linewidth = 0.3) +
  scale_fill_substrate() +
  labs(x = "Soil depth", y = "Mean root biomass (g)", fill = "Substrate") +
  theme_susana() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
save_pub_plot(p_depth,
              file.path(out_dir, "root_depth_profile.png"),
              width = 200, height = 130)

# 4. Root vs shoot scatter
p_scatter <- ggplot(rs_df,
                    aes(x = total_shoot_biomass, y = total_root_biomass,
                        colour = substrate)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = TRUE, alpha = 0.1, linewidth = 0.5, colour = "grey30") +
  scale_colour_substrate() +
  labs(x = "Shoot biomass (g)", y = "Root biomass (g)",
       colour = "Substrate") +
  theme_susana() +
  theme(legend.position = "bottom")
save_pub_plot(p_scatter,
              file.path(out_dir, "root_vs_shoot_scatter.png"),
              width = 160, height = 130)

# 5. Height vs shoot biomass
if (nrow(height_shoot) > 0 && !is.null(cor_hs)) {
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
