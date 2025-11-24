# ============================================================================
# COMPLETE PLS-DA ANALYSIS PIPELINE
# ============================================================================
# From data import through PLS-DA to multivariate discrimination visualization
# Generates biplots and loadings plots for garden environment classification
# Generated for scientific publication
# ============================================================================

# --- LOAD REQUIRED LIBRARIES ---
library(mixOmics)
library(dplyr)
library(xlsx)
library(patchwork)
library(ggrepel)
library(ggplot2)

# ============================================================================
# SECTION 1: DATA IMPORT (MATCHING PERMANOVA PIPELINE)
# ============================================================================

# --- Read metabolite data ---
metab <- read.xlsx("Table_S2.xlsx", 1)

# --- Read phenotypic data ---
pheno <- read.xlsx("Table_S3.xlsx", 1)

# --- Rename phenotypic trait columns ---
colnames(pheno)[-(1:6)] <- c("Leaves_mean", "Leaves_increase", "Rosettes_mean", "Rosettes_increase",
                             "Volume_mean", "Volume_increase", "Total_ripe_fruits",
                             "Total_runners", "Fruit_volume_mean", "Fertilized_seeds_mean", "Leaf_damage_mean",
                             "Leaf_damage_CV", "Leaf_damage_max")

# ============================================================================
# SECTION 2: PREPARE SUBSETS (MATCHING PERMANOVA PIPELINE)
# ============================================================================

# --- Extract metabolite data ---
metab_subset <- metab[, 7:ncol(metab)]
metab_subset <- as.data.frame(metab_subset)
metab_subset$GARDENID <- metab$GARDENID
metab_subset$GENOID <- metab$GENOID

# --- Extract phenotypic data ---
pheno_subset <- pheno[, 7:ncol(pheno)]
pheno_subset <- as.data.frame(pheno_subset)
pheno_subset$GARDENID <- pheno$GARDENID
pheno_subset$GENOID <- pheno$GENOID

# ============================================================================
# SECTION 3: IDENTIFY COMMON SAMPLES (MATCHING PERMANOVA PIPELINE)
# ============================================================================

# --- Create unique sample identifiers ---
metab$ID <- paste(metab$GARDENID, metab$GENOID, metab$BLOCK, sep = "_")
pheno$ID <- paste(pheno$GARDENID, pheno$GENOID, pheno$BLOCK, sep = "_")

# --- Find common samples ---
common_samples <- intersect(metab$ID, pheno$ID)

# --- Filter to common samples ---
metab_subset <- metab_subset[metab$ID %in% common_samples, ]
pheno_subset <- pheno_subset[pheno$ID %in% common_samples, ]

# --- Sort by same order ---
metab_subset <- metab_subset[order(metab$ID[metab$ID %in% common_samples]), ]
pheno_subset <- pheno_subset[order(pheno$ID[pheno$ID %in% common_samples]), ]

# ============================================================================
# SECTION 4: COMBINE DATA (MATCHING PERMANOVA PIPELINE)
# ============================================================================

# --- Merge metabolites and phenotypes ---
combined_df <- cbind(metab_subset[, (37:38)], metab_subset[, -(37:38)], pheno_subset[, -(14:15)])
rownames(combined_df) <- common_samples

# ============================================================================
# SECTION 5: LOG2 TRANSFORMATION (MATCHING PERMANOVA PIPELINE)
# ============================================================================

# --- Handle negative values ---
min_value <- min(combined_df[, -(1:2)], na.rm = TRUE)

if (min_value < 0) {
  combined_log2 <- combined_df %>%
    mutate_at(vars(-(1:2)), ~ log2(. + abs(min_value) + 1))
} else {
  combined_log2 <- combined_df %>%
    mutate_at(vars(-(1:2)), log2)
}

# ============================================================================
# SECTION 6: Z-SCORE SCALING (MATCHING PERMANOVA PIPELINE)
# ============================================================================

# --- Apply z-score scaling ---
combined_df <- combined_log2 %>%
  mutate_at(vars(-(1:2)), scale)

# ============================================================================
# SECTION 7: PREPARE DATA FOR PLS-DA
# ============================================================================

# --- Define predictor matrix (X) and response vector (Y) ---
X <- combined_df[, -(1:2)]
Y <- as.factor(combined_df$GARDENID)

# ============================================================================
# SECTION 8: PERFORM PLS-DA WITH 3 COMPONENTS
# ============================================================================

# --- Run PLS-DA (data already scaled) ---
plsda_result <- plsda(X, Y, ncomp = 3, scale = FALSE)

# ============================================================================
# SECTION 9: VISUALIZATION PARAMETERS
# ============================================================================

# --- Parameters for biplot ---
n_top_vars <- 15      # Number of top variables to display
arrow_scale <- 7      # Scaling factor for arrow length

# --- Dictionary for pretty names ---
ugly_names <- colnames(X)

pretty_metab_names <- c(
  "Pyruvic acid", "Valine", "Isoleucine", "Glycine", 
  "Phosphoric acid", "Proline", "Urea", "Glyceric acid", 
  "Alanine", "Serine", "Succinic acid", "Threonine", 
  "Fumaric acid", "Nicotinic acid", "Erythritol", "Malic acid", 
  "GABA", "Aspartic acid", "Threonic acid", "Xylose", 
  "Xylitol", "Pyroglutamic acid", "Glutamic acid", "Arginine", 
  "Tryptophan", "Trehalose", "Maltitol", "Galactinol", 
  "Quinic acid", "Fructose", "Glucose", "Citric acid", 
  "Methyl α-D-glucopyranoside", "Dehydroascorbic acid", "myo-Inositol", "Sucrose"
)

pretty_pheno_names <- c(
  "Mean leaf number", "Leaf number change", "Mean rosette number",
  "Rosette number change", "Mean plant volume", "Plant volume change",
  "Total ripe fruits", "Total runners", "Mean fruit volume",
  "Fertilized achenes (%)", "Mean leaf damage (%)", "Leaf damage (CV)",
  "Max leaf damage (%)"
)

all_pretty_names <- c(pretty_metab_names, pretty_pheno_names)
name_mapper <- setNames(all_pretty_names, ugly_names)

# --- Define garden colors ---
# --- Define garden colors (matching RColorBrewer Set2 palette) ---
garden_colors <- c(
  "ALN_21" = "#66C2A5",
  "GON_21" = "#FC8D62",
  "GON_22" = "#8DA0CB",
  "KEV_22" = "#E78AC3",
  "RUI_21" = "#A6D854"
)


# ============================================================================
# SECTION 10: PREPARE DATA FOR BIPLOT VISUALIZATION
# ============================================================================

# --- Extract sample coordinates (scores) ---
ind_coords <- as.data.frame(plsda_result$variates$X)
ind_coords$GARDENID <- Y
levels(ind_coords$GARDENID) <- c("ALN_21", "GON_21", "GON_22", "KEV_22", "RUI_21")

# --- Calculate explained variance for axis labels ---
explained_variance <- plsda_result$prop_expl_var$X * 100
xlab <- paste0("LV 1 (", round(explained_variance[1], 2), "%)")
ylab_2 <- paste0("LV 2 (", round(explained_variance[2], 2), "%)")
ylab_3 <- paste0("LV 3 (", round(explained_variance[3], 2), "%)")

# ============================================================================
# SECTION 11: PREPARE LOADINGS FOR LV1 vs LV2
# ============================================================================

# --- Extract loadings and select top variables ---
loadings_12 <- as.data.frame(plsda_result$loadings$X[, 1:2])
colnames(loadings_12) <- c("comp1", "comp2")
loadings_12$varname <- rownames(loadings_12)

loadings_12 <- loadings_12 %>%
  mutate(magnitude = sqrt(comp1^2 + comp2^2)) %>%
  arrange(desc(magnitude)) %>%
  slice_head(n = n_top_vars) %>%
  mutate(pretty_label = name_mapper[varname],
         comp1_scaled = comp1 * arrow_scale,
         comp2_scaled = comp2 * arrow_scale)

# ============================================================================
# SECTION 12: PREPARE LOADINGS FOR LV1 vs LV3
# ============================================================================

# --- Extract loadings and select top variables ---
loadings_13 <- as.data.frame(plsda_result$loadings$X[, c(1, 3)])
colnames(loadings_13) <- c("comp1", "comp3")
loadings_13$varname <- rownames(loadings_13)

loadings_13 <- loadings_13 %>%
  mutate(magnitude = sqrt(comp1^2 + comp3^2)) %>%
  arrange(desc(magnitude)) %>%
  slice_head(n = n_top_vars) %>%
  mutate(pretty_label = name_mapper[varname],
         comp1_scaled = comp1 * arrow_scale,
         comp3_scaled = comp3 * arrow_scale)

# ============================================================================
# SECTION 13: CREATE BIPLOT LV1 vs LV2
# ============================================================================

LV1_2 <- ggplot(data = ind_coords, aes(x = comp1, y = comp2, color = GARDENID)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = GARDENID), type = "norm", level = 0.95, alpha = 0.2, 
               geom = "polygon", show.legend = FALSE) +
  geom_segment(data = loadings_12, aes(x = 0, y = 0, xend = comp1_scaled, yend = comp2_scaled),
               arrow = arrow(length = unit(0.2, 'cm')), color = 'black', alpha = 0.7) +
  geom_text_repel(data = loadings_12, aes(x = comp1_scaled, y = comp2_scaled, label = pretty_label),
                  size = 4, color = 'black', max.overlaps = Inf, fontface = "bold") +
  scale_color_manual(values = garden_colors) +
  scale_fill_manual(values = garden_colors) +
  labs(title = "", subtitle = "LV1 vs LV2", x = xlab, y = ylab_2) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

# ============================================================================
# SECTION 14: CREATE BIPLOT LV1 vs LV3
# ============================================================================

LV1_3 <- ggplot(data = ind_coords, aes(x = comp1, y = comp3, color = GARDENID)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = GARDENID), type = "norm", level = 0.95, alpha = 0.2, 
               geom = "polygon", show.legend = FALSE) +
  geom_segment(data = loadings_13, aes(x = 0, y = 0, xend = comp1_scaled, yend = comp3_scaled),
               arrow = arrow(length = unit(0.2, 'cm')), color = 'black', alpha = 0.7) +
  geom_text_repel(data = loadings_13, aes(x = comp1_scaled, y = comp3_scaled, label = pretty_label),
                  size = 4, color = 'black', max.overlaps = Inf, fontface = "bold") +
  scale_color_manual(values = garden_colors) +
  scale_fill_manual(values = garden_colors) +
  labs(title = "", subtitle = "LV1 vs LV3", x = xlab, y = ylab_3) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

# ============================================================================
# SECTION 15: COMBINE PLOTS WITH PATCHWORK AND SHARED LEGEND
# ============================================================================

# --- Combine both biplots with shared legend at the bottom ---
PLSDA_FINAL <- (LV1_2 + LV1_3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key = element_rect(fill = "white", color = NA)
  ) &
  guides(color = guide_legend(title = "Common Garden", nrow = 1),
         fill = guide_legend(title = "Common Garden", nrow = 1))

# ============================================================================
# SECTION 16: DISPLAY FINAL PLOT
# ============================================================================

# --- Display final plot ---
print(PLSDA_FINAL)

# ============================================================================
# END OF SCRIPT
# ============================================================================