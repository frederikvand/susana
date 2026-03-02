# =============================================================================
# 00_data_prep.R
# Reads all raw Excel data, cleans, standardises, and saves as .rds
# =============================================================================

cat("\n======================================================\n")
cat("  00_data_prep.R — Data ingestion and cleaning\n")
cat("======================================================\n\n")

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(here)
})

source(here("scripts", "shared_theme.R"))
set.seed(2024)

# --- Paths ---
data_root   <- here("data")
out_root    <- here("output", "clean_data")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)


# =====================================================
# 1. SUBSTRATE CHARACTERISATION
# =====================================================
cat("--- Substrate data ---\n")

# CN soil: alternative vs control, sampled in/out of tidal zone
cn_raw <- read_excel(file.path(data_root, "substrate/CN_soil.xlsx"), sheet = 1)
cn_df <- cn_raw %>%
  rename(pct_TN = `%TN`, pct_TOC = `%TOC`, in_out = in_uit) %>%
  mutate(
    substrate = factor(substrate, levels = c("control", "alternative")),
    in_out    = factor(in_out, levels = c("in", "uit"),
                       labels = c("Inside", "Outside")),
    pct_TOC   = as.numeric(pct_TOC),
    pct_TN    = as.numeric(pct_TN)
  ) %>%
  filter(!is.na(substrate), !is.na(pct_TOC))

cat("  CN samples:", nrow(cn_df), "\n")

# Grain size data
gs_raw <- read_excel(file.path(data_root, "substrate/GS_soil.xlsx"), sheet = 1)
# The grain size data has d(0.5) = D50. Link via sample ID
cat("  GS samples:", nrow(gs_raw), "\n")

# Substrate sampling: includes TOC, TN, d0.5, substrate, in_out
sub_samp <- tryCatch({
  df <- read_excel(file.path(data_root, "substrate/substrate_sampling.xlsx"))
  df %>% mutate(
    substrate = factor(substrate, levels = c("control", "alternative")),
    in_out    = factor(in_out, levels = c("in", "out"),
                       labels = c("Inside", "Outside")),
    d50       = as.numeric(d0.5)
  ) %>%
  filter(!is.na(substrate))
}, error = function(e) {
  cat("  WARNING: Could not read substrate_sampling.xlsx:", e$message, "\n")
  NULL
})
if (!is.null(sub_samp)) {
  cat("  Substrate sampling:", nrow(sub_samp), "\n")
}

# Calculated substrate properties (10-substrate gradient summary)
calc_props <- read_excel(
  file.path(data_root, "substrate/calculated_substrate_properties.xlsx")
) %>%
  rename(substrate_label = Substrate,
         pct_TOC = `TOC (%)`,
         pct_TN  = `N (%)`,
         d50_um  = `D50 (µm)`) %>%
  mutate(substrate = standardise_substrate(substrate_label))

cat("  Calculated properties:", nrow(calc_props), "substrates\n")

saveRDS(cn_df,       file.path(out_root, "substrate_cn.rds"))
saveRDS(gs_raw,      file.path(out_root, "substrate_gs.rds"))
saveRDS(calc_props,  file.path(out_root, "substrate_properties.rds"))
if (!is.null(sub_samp)) {
  saveRDS(sub_samp,  file.path(out_root, "substrate_sampling.rds"))
}
cat("  -> saved substrate .rds files\n\n")


# =====================================================
# 2. ABOVEGROUND BIOMASS (common garden)
# =====================================================
cat("--- Aboveground biomass ---\n")

ag_raw <- read_excel(
  file.path(data_root,
  "vegetation_development/common_garden/aboveground_biomass/aboveground_biomass_zeebrugge.xlsx"),
  sheet = 1
)

ag_df <- ag_raw %>%
  mutate(
    substrate    = standardise_substrate(substrate),
    plant_status = factor(plant,
                          levels = c("controle", "plant", "replant")),
    length       = as.numeric(length),
    leaf_count   = as.numeric(nr),
    timepoint    = factor(t, levels = paste0("t", 1:13), ordered = TRUE),
    time_numeric = as.numeric(gsub("t", "", t)),
    plant_id     = paste(column, row, sep = "_")
  ) %>%
  filter(!is.na(substrate))

# Separate controls and planted
ag_controls <- ag_df %>% filter(plant_status == "controle")
ag_planted  <- ag_df %>% filter(plant_status %in% c("plant", "replant"))

cat("  Total observations:", nrow(ag_df), "\n")
cat("  Planted plants:", nrow(ag_planted),
    "  Controls:", nrow(ag_controls), "\n")
cat("  Timepoints:", length(unique(ag_df$timepoint)), "\n")
cat("  Unique plants:", length(unique(ag_planted$plant_id)), "\n")

saveRDS(ag_df,       file.path(out_root, "aboveground_all.rds"))
saveRDS(ag_planted,  file.path(out_root, "aboveground_planted.rds"))
saveRDS(ag_controls, file.path(out_root, "aboveground_controls.rds"))
cat("  -> saved aboveground .rds files\n\n")


# =====================================================
# 3. BELOWGROUND BIOMASS (roots)
# =====================================================
cat("--- Belowground biomass ---\n")

bg_raw <- read_excel(
  file.path(data_root,
  "vegetation_development/common_garden/belowground_biomass/common_garden_roots.xlsx"),
  sheet = 1
)

bg_df <- bg_raw %>%
  mutate(
    substrate     = standardise_substrate(substrate),
    layer         = factor(layer,
                           levels = c("shoots", "top", "1", "2", "3", "4"),
                           ordered = TRUE),
    total_biomass = as.numeric(total_biomass),
    biomass_per_L = as.numeric(biomass_per_L),
    shoot_biomass = as.numeric(shoot_biomass),
    root_biomass  = as.numeric(root_biomass),
    plant_id      = paste(column, row, sep = "_")
  ) %>%
  filter(!is.na(substrate))

cat("  Observations:", nrow(bg_df), "\n")
cat("  Unique plants:", length(unique(bg_df$plant_id)), "\n")

saveRDS(bg_df, file.path(out_root, "belowground.rds"))
cat("  -> saved belowground .rds file\n\n")


# =====================================================
# 4. GERMINATION EXPERIMENT
# =====================================================
cat("--- Germination data ---\n")

germ_raw <- read_excel(
  file.path(data_root,
  "vegetation_development/germination/kiemingsexperimenten data.xlsx"),
  sheet = 1
)

# Map germination substrate names (D, D75, D50, D25, W, R)
germ_substrate_map <- c(
  "D"   = "D",
  "D75" = "D75",
  "D50" = "D50",
  "D25" = "D25",
  "W"   = "W",
  "R"   = "R"
)
# The germination experiment uses only 6 of the 10 substrates
GERM_SUBSTRATE_LEVELS <- c("R", "W", "D25", "D50", "D75", "D")

germ_df <- germ_raw %>%
  mutate(
    substrate_raw = as.character(substrate),
    substrate     = factor(substrate_raw, levels = GERM_SUBSTRATE_LEVELS,
                           ordered = TRUE),
    leaf_count    = as.numeric(nr),
    length        = as.numeric(length),
    timepoint     = factor(t, levels = paste0("t", 1:13), ordered = TRUE),
    time_numeric  = as.numeric(gsub("t", "", t)),
    pot_id        = paste(column, row, sep = "_")
  ) %>%
  filter(!is.na(substrate))

# Compute germination percentage at t1 (10 seeds planted per pot)
germ_t1 <- germ_df %>%
  filter(timepoint == "t1") %>%
  mutate(pct_germinated = leaf_count * 10)

cat("  Observations:", nrow(germ_df), "\n")
cat("  Unique pots:", length(unique(germ_df$pot_id)), "\n")
cat("  Substrates:", paste(levels(germ_df$substrate), collapse = ", "), "\n")

saveRDS(germ_df,  file.path(out_root, "germination.rds"))
saveRDS(germ_t1,  file.path(out_root, "germination_t1.rds"))
cat("  -> saved germination .rds files\n\n")


# =====================================================
# 5. MICROCLIMATE LOGGER DATA
# =====================================================
cat("--- Logger data ---\n")

logger_dir <- file.path(data_root,
                         "vegetation_development/germination/logger data")

# Logger ID -> substrate mapping
logger_map <- data.frame(
  file      = c("l_Dc.csv", "l_75c.csv", "l_50c.csv",
                "l_25c.csv", "l_Wc.csv", "l_Rc.csv"),
  substrate = c("D", "D75", "D50", "D25", "W", "R"),
  stringsAsFactors = FALSE
)

logger_list <- list()
for (i in seq_len(nrow(logger_map))) {
  fpath <- file.path(logger_dir, logger_map$file[i])
  if (!file.exists(fpath)) {
    cat("  WARNING: missing", fpath, "\n")
    next
  }
  tmp <- read.csv(fpath, header = FALSE, sep = ";",
                   stringsAsFactors = FALSE)
  # TMS-4 logger format: index; datetime; T1; T2; T3; moisture_raw; Vol; flag
  if (ncol(tmp) >= 7) {
    colnames(tmp)[1:7] <- c("index", "datetime", "T1_soil", "T2_surface",
                             "T3_air", "moisture_raw", "Vol_moisture")
  }
  tmp$substrate_raw <- logger_map$substrate[i]
  tmp <- tmp %>%
    mutate(
      datetime     = as.POSIXct(datetime, format = "%Y.%m.%d %H:%M"),
      T1_soil      = as.numeric(T1_soil),
      T2_surface   = as.numeric(T2_surface),
      T3_air       = as.numeric(T3_air),
      Vol_moisture = as.numeric(Vol_moisture),
      substrate    = factor(substrate_raw,
                            levels = GERM_SUBSTRATE_LEVELS, ordered = TRUE)
    ) %>%
    filter(!is.na(datetime))
  logger_list[[i]] <- tmp
}

loggers_all <- bind_rows(logger_list) %>%
  arrange(substrate, datetime) %>%
  mutate(
    date = as.Date(datetime),
    hour = as.integer(format(datetime, "%H"))
  )

cat("  Logger records:", nrow(loggers_all), "\n")
cat("  Date range:", as.character(min(loggers_all$date, na.rm = TRUE)),
    "to", as.character(max(loggers_all$date, na.rm = TRUE)), "\n")

saveRDS(loggers_all, file.path(out_root, "loggers.rds"))
cat("  -> saved loggers .rds file\n\n")


cat("======================================================\n")
cat("  DATA PREP COMPLETE\n")
cat("======================================================\n")
