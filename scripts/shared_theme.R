# =============================================================================
# shared_theme.R
# Publication-quality plotting theme and helpers for the SUSANA manuscript
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(patchwork)
})

# --- Substrate ordering and TOC lookup ---
SUBSTRATE_LEVELS <- c("R", "W", "D12.5", "D25", "D37.5", "D50",
                       "D62.5", "D75", "D87.5", "D")
SUBSTRATE_LABELS <- c("R", "W", "12.5", "25", "37.5", "50",
                       "62.5", "75", "87.5", "D")

# Lookup table: substrate -> mean TOC (%)
SUBSTRATE_TOC <- c(R = 0.06, W = 0.25, D12.5 = 0.49, D25 = 0.73,
                   D37.5 = 0.97, D50 = 1.21, D62.5 = 1.45,
                   D75 = 1.68, D87.5 = 1.92, D = 2.16)

# --- Publication theme ---
theme_susana <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Panel
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(fill = NA, colour = "grey30", linewidth = 0.4),
      # Axes
      axis.title     = element_text(size = rel(1.0), face = "plain"),
      axis.text       = element_text(size = rel(0.85), colour = "grey30"),
      axis.ticks      = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks.length = unit(1.5, "mm"),
      # Legend
      legend.title    = element_text(size = rel(0.9), face = "bold"),
      legend.text     = element_text(size = rel(0.85)),
      legend.key.size = unit(5, "mm"),
      legend.position = "right",
      # Strip (facet labels)
      strip.text      = element_text(size = rel(0.95), face = "bold",
                                     margin = margin(3, 3, 3, 3)),
      strip.background = element_rect(fill = "grey95", colour = "grey30",
                                      linewidth = 0.3),
      # Plot title
      plot.title       = element_text(size = rel(1.15), face = "bold",
                                      hjust = 0, margin = margin(b = 6)),
      plot.subtitle    = element_text(size = rel(0.95), hjust = 0,
                                      margin = margin(b = 6)),
      plot.margin      = margin(8, 8, 8, 8)
    )
}


# --- Viridis-based substrate colour/fill scales ---
# Ordered D-heavy (dark purple) to R (yellow) using viridis option D
substrate_palette <- viridis(n = length(SUBSTRATE_LEVELS), option = "D",
                              begin = 0.05, end = 0.95)
names(substrate_palette) <- SUBSTRATE_LEVELS

scale_fill_substrate <- function(...) {
  scale_fill_manual(values = substrate_palette, ...)
}

scale_colour_substrate <- function(...) {
  scale_colour_manual(values = substrate_palette, ...)
}


# --- Helper: save publication figure as PDF and PNG ---
save_pub_plot <- function(plot, filename, width = 180, height = 120,
                          units = "mm", dpi = 600) {
  # Save as PDF (vector) for typesetting
  pdf_file <- sub("\\.[^.]+$", ".pdf", filename)
  ggsave(pdf_file, plot = plot, width = width, height = height,
         units = units, device = cairo_pdf)

  # Save as PNG (raster) for quick preview / supplementary
  png_file <- sub("\\.[^.]+$", ".png", filename)
  ggsave(png_file, plot = plot, width = width, height = height,
         units = units, dpi = dpi)

  cat("  Saved:", pdf_file, "\n")
  cat("  Saved:", png_file, "\n")
}


# --- Helper: standardise substrate factor from raw data ---
standardise_substrate <- function(x) {
  # Map common raw notations to the canonical form
  # Handles both dot and comma decimals (European format)
  mapping <- c(
    "R"    = "R",
    "W"    = "W",
    "12.5" = "D12.5", "12,5" = "D12.5", "12"   = "D12.5", "D12.5" = "D12.5",
    "25"   = "D25",   "D25"  = "D25",
    "37.5" = "D37.5", "37,5" = "D37.5", "37"   = "D37.5", "D37.5" = "D37.5",
    "50"   = "D50",   "D50"  = "D50",
    "62.5" = "D62.5", "62,5" = "D62.5", "62"   = "D62.5", "D62.5" = "D62.5",
    "75"   = "D75",   "D75"  = "D75",
    "87.5" = "D87.5", "87,5" = "D87.5", "87"   = "D87.5", "D87.5" = "D87.5",
    "D"    = "D"
  )
  mapped <- mapping[as.character(x)]
  factor(mapped, levels = SUBSTRATE_LEVELS, ordered = TRUE)
}


# --- Helper: add compact letter display (significance grouping) ---
add_significance_letters <- function(df, y_col, group_col,
                                     y_nudge = 0.05, size = 3) {
  # Expects a data.frame with columns for the value and group,
  # plus a column 'cld' containing the compact letter display
  geom_text(data = df, aes(x = .data[[group_col]],
                            y = .data[[y_col]], label = cld),
            vjust = -0.5, size = size, fontface = "bold")
}


cat("  shared_theme.R loaded\n")
