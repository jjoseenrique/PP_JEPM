# ============================================================================
# PERFORMANCE-RELATED DATA ANALYSIS WITH HEATMAP VISUALIZATION
# ============================================================================
# From data import through normalization to plant performance-climatic correlations
# Generates heatmaps and correlations with bioclimatic variables
# Generated for scientific publication
# ============================================================================

# --- LOAD REQUIRED LIBRARIES ---
library(circlize)
library(dplyr)
library(xlsx)
library(ComplexHeatmap)
library(RColorBrewer)
library(nasapower)
library(tidyverse)
library(grid)

# ============================================================================
# SECTION 1: DATA PREPARATION AND NORMALIZATION
# ============================================================================

# --- Helper function to clean NaN values ---
clean_nan <- function(df) {
  as.data.frame(lapply(df, function(x) {
    x[is.nan(x)] <- NA
    return(x)
  }))
}

# --- Read phenotypic data from Excel ---
pheno <- read.xlsx("Table_S3.xlsx", 1)

# --- Rename phenotypic trait columns ---
colnames(pheno)[-(1:6)] <- c("Leaves_mean", "Leaves_increase", "Rosettes_mean", "Rosettes_increase",
                             "Volume_mean", "Volume_increase", "Total_ripe_fruits",
                             "Total_runners", "Fruit_volume_mean", "Fertilized_seeds_mean", "Leaf_damage_mean",
                             "Leaf_damage_CV", "Leaf_damage_max")

# --- Select phenotypic trait columns (from column 7 onwards) ---
pheno_subset <- pheno[, 7:ncol(pheno)]
pheno_subset <- as.data.frame(pheno_subset)

# --- Add grouping columns: GARDENID and GENOID ---
pheno_subset$GARDENID <- pheno$GARDENID
pheno_subset$GENOID <- pheno$GENOID

# ============================================================================
# SECTION 2: NORMALITY TEST AND LOG2 TRANSFORMATION WITH OFFSET
# ============================================================================

# --- Perform Shapiro-Wilk normality test (excluding ID columns) ---
shapiro <- apply(pheno_subset[, -((ncol(pheno_subset) - 1):ncol(pheno_subset))], 2, shapiro.test)

# --- Extract p-values from Shapiro-Wilk test ---
p_values <- sapply(shapiro, function(x) x$p.value)

# --- Filter metabolites with p-value > 0.05 ---
significant <- p_values[p_values > 0.05]

# --- Find minimum value to handle negative numbers ---
min_value <- min(pheno_subset[, -c((ncol(pheno_subset) - 1):ncol(pheno_subset))], na.rm = TRUE)

# --- Apply log2 transformation with offset for negative values ---
if (min_value < 0) {
  pheno_log2 <- pheno_subset %>%
    mutate_at(vars(-c((ncol(pheno_subset) - 1):ncol(pheno_subset))), ~ log2(. + abs(min_value) + 1))
} else {
  pheno_log2 <- pheno_subset %>%
    mutate_at(vars(-c((ncol(pheno_subset) - 1):ncol(pheno_subset))), log2)
}

# --- Apply z-score scaling to normalized data ---
pheno_subset <- pheno_log2 %>%
  mutate_at(vars(-c((ncol(pheno_subset) - 1):ncol(pheno_subset))), scale)

# ============================================================================
# SECTION 3: CALCULATE MEAN VALUES BY GARDEN AND GENOTYPE
# ============================================================================

# --- Calculate mean phenotypic values by GARDENID and GENOID ---
pheno_mean <- pheno_subset %>%
  group_by(GARDENID, GENOID) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = 'drop')

# --- Standardize garden IDs with year information ---
pheno_mean$GARDENID <- gsub("ALN", "ALN_21", pheno_mean$GARDENID)
pheno_mean$GARDENID <- gsub("RUI", "RUI_21", pheno_mean$GARDENID)
pheno_mean$GARDENID <- gsub("KEV", "KEV_22", pheno_mean$GARDENID)

# --- Create row names combining GENOID and GARDENID ---
pheno_mean$RowName <- paste(pheno_mean$GENOID, pheno_mean$GARDENID, sep = "_")
pheno_mean <- as.data.frame(pheno_mean)
rownames(pheno_mean) <- pheno_mean$RowName

# --- Store filtered phenotypic data for later merging ---
pheno_filtered <- pheno_mean

# --- Remove ID columns and convert to data frame ---
pheno_mean <- pheno_mean[, -c(1:2, ncol(pheno_mean))]
pheno_mean <- as.data.frame(pheno_mean)
row_names_pheno <- row.names(pheno_mean)

# --- Clean NaN values ---
pheno_mean <- clean_nan(pheno_mean)
row.names(pheno_mean) <- row_names_pheno

# ============================================================================
# SECTION 4: PREPARE DATA FOR HEATMAP VISUALIZATION
# ============================================================================

# --- Transpose data for heatmap (traits as rows, samples as columns) ---
heatmap_data <- t(pheno_mean)

# ============================================================================
# SECTION 5: RETRIEVE CLIMATIC DATA FROM NASA POWER API
# ============================================================================

# --- Define garden locations ---
gardens_data <- data.frame(
  name = c("GON", "ALN", "RUI", "KEV"),
  lat = c(50.98405, 55.65765, 60.433333, 69.750),
  lon = c(3.79869, 13.082, 22.16666, 27.017),
  lab = c("GON", "ALN", "RUI", "KEV")
)

# --- Download climatic data from NASA POWER for each garden ---
# GON garden data
GON_21 <- get_power(
  community = "AG",
  lonlat = c(3.79869, 50.98405),
  pars = c("T2M", "T2M_MIN", "T2M_MAX", "PRECTOTCORR"),
  dates = c("2021", "2022"),
  temporal_api = "monthly"
)

GON_22 <- get_power(
  community = "AG",
  lonlat = c(3.79869, 50.98405),
  pars = c("T2M", "T2M_MIN", "T2M_MAX", "PRECTOTCORR"),
  dates = c("2022", "2023"),
  temporal_api = "monthly"
)

# ALN garden data
ALN_21 <- get_power(
  community = "AG",
  lonlat = c(13.082, 55.65765),
  pars = c("T2M", "T2M_MIN", "T2M_MAX", "PRECTOTCORR"),
  dates = c("2021", "2022"),
  temporal_api = "monthly"
)

# RUI garden data
RUI_21 <- get_power(
  community = "AG",
  lonlat = c(22.16666, 60.433333),
  pars = c("T2M", "T2M_MIN", "T2M_MAX", "PRECTOTCORR"),
  dates = c("2021", "2022"),
  temporal_api = "monthly"
)

# KEV garden data
KEV_22 <- get_power(
  community = "AG",
  lonlat = c(27.017, 69.750),
  pars = c("T2M", "T2M_MIN", "T2M_MAX", "PRECTOTCORR"),
  dates = c("2021", "2022"),
  temporal_api = "monthly"
)

# ============================================================================
# SECTION 6: PROCESS CLIMATIC DATA
# ============================================================================

# --- Extract annual data for each garden and year ---
GON_21 <- GON_21[GON_21$YEAR == 2021, c("LON", "LAT", "PARAMETER", "YEAR", "ANN")]
GON_22 <- GON_22[GON_22$YEAR == 2022, c("LON", "LAT", "PARAMETER", "YEAR", "ANN")]
ALN_21 <- ALN_21[ALN_21$YEAR == 2021, c("LON", "LAT", "PARAMETER", "YEAR", "ANN")]
RUI_21 <- RUI_21[RUI_21$YEAR == 2021, c("LON", "LAT", "PARAMETER", "YEAR", "ANN")]
KEV_22 <- KEV_22[KEV_22$YEAR == 2022, c("LON", "LAT", "PARAMETER", "YEAR", "ANN")]

# --- Add garden identifiers ---
GON_21 <- cbind(GARDEN = c("GON_21", "GON_21"), GON_21)
GON_22 <- cbind(GARDEN = c("GON_22", "GON_22"), GON_22)
ALN_21 <- cbind(GARDEN = c("ALN_21", "ALN_21"), ALN_21)
RUI_21 <- cbind(GARDEN = c("RUI_21", "RUI_21"), RUI_21)
KEV_22 <- cbind(GARDEN = c("KEV_22", "KEV_22"), KEV_22)

# --- Combine all climatic data ---
Clim_info <- bind_rows(GON_21, GON_22, ALN_21, RUI_21, KEV_22)
Clim_info <- Clim_info[, c("GARDEN", "LON", "LAT", "PARAMETER", "YEAR", "ANN")]

# --- Save climatic information to Excel ---
write.xlsx2(Clim_info, "Climatic_info.xlsx", row.names = FALSE, na.rm = TRUE)

# ============================================================================
# SECTION 7: SIMPLIFY CLIMATIC DATA
# ============================================================================

# --- Remove NAs and shift values to the left ---
shift_left <- function(x) {
  x[!is.na(x)]
}

Clim_info_simp <- t(apply(Clim_info, 1, shift_left))
Clim_info_simp <- as.data.frame(Clim_info_simp)
names(Clim_info_simp) <- c("GARDEN", "LON", "LAT", "PARAMETER", "YEAR", "ANN")

# --- Convert ANN to numeric ---
Clim_info_simp$ANN <- as.numeric(Clim_info_simp$ANN)

# --- Pivot climatic data to have parameters as columns ---
Clim_info_pivot <- Clim_info_simp %>%
  pivot_wider(names_from = PARAMETER, values_from = ANN)

Clim_info_pivot <- Clim_info_pivot[, c("GARDEN", "PRECTOTCORR", "T2M", "T2M_MAX", "T2M_MIN")]

# ============================================================================
# SECTION 8: CALCULATE CORRELATIONS BETWEEN PHENOTYPES AND CLIMATIC VARIABLES
# ============================================================================

# --- Clean phenotypic data ---
pheno_filtered <- pheno_filtered[, !names(pheno_filtered) %in% "RowName"]
pheno_filtered <- clean_nan(pheno_filtered)

# --- Merge phenotypic and climatic data ---
combined_data <- merge(pheno_filtered, Clim_info_pivot, by.x = "GARDENID", by.y = "GARDEN")

# --- Convert climatic variables to numeric ---
combined_data_numeric <- combined_data %>%
  mutate(across(c(PRECTOTCORR, T2M, T2M_MAX, T2M_MIN), as.numeric))

# --- Extract phenotypic variables and climatic variables ---
variables <- combined_data_numeric[, 3:(ncol(combined_data_numeric) - 4)]
last_4_vars <- combined_data_numeric[, (ncol(combined_data_numeric) - 3):ncol(combined_data_numeric)]

# --- Verify data types ---
if (!all(sapply(variables, is.numeric))) {
  stop("Some columns in 'variables' are not numeric")
}

if (!all(sapply(last_4_vars, is.numeric))) {
  stop("Some columns in 'last_4_vars' are not numeric")
}

# --- Initialize correlation and p-value matrices ---
correlation_matrix <- matrix(NA, nrow = ncol(variables), ncol = ncol(last_4_vars))
p_value_matrix <- matrix(NA, nrow = ncol(variables), ncol = ncol(last_4_vars))

# --- Calculate correlation coefficients and p-values ---
for (i in 1:ncol(variables)) {
  for (j in 1:ncol(last_4_vars)) {
    test <- cor.test(variables[, i], last_4_vars[, j], use = "pairwise.complete.obs")
    correlation_matrix[i, j] <- test$estimate
    p_value_matrix[i, j] <- test$p.value
  }
}

# --- Convert matrices to data frames and assign row/column names ---
correlation_matrix <- as.data.frame(correlation_matrix)
p_value_matrix <- as.data.frame(p_value_matrix)
rownames(correlation_matrix) <- colnames(variables)
colnames(correlation_matrix) <- colnames(last_4_vars)
rownames(p_value_matrix) <- colnames(variables)
colnames(p_value_matrix) <- colnames(last_4_vars)

# --- Identify non-significant correlations (p > 0.05) ---
significance_threshold <- 0.05
non_significant_traits <- rownames(p_value_matrix)[apply(p_value_matrix, 1, function(row) all(row >= significance_threshold))]

# --- Rename climate variables for better readability ---
colnames(correlation_matrix) <- c("Precipitation", "T_mean", "T_max", "T_min")

# ============================================================================
# SECTION 9: CREATE MAIN HEATMAP (Z-SCORE PHENOTYPES)
# ============================================================================

# --- Define color palette for phenotype heatmap ---
colores <- colorRamp2(c(-2, 0, 2), c("navy", "ivory", "firebrick3"))

# --- Extract garden identifiers from column names ---
garden_ids <- gsub("G\\d+_", "", colnames(heatmap_data))
unique_gardens <- unique(garden_ids)
garden_colors <- brewer.pal(length(unique_gardens), "Set2")
names(garden_colors) <- unique_gardens

# --- Create garden annotation for top of heatmap ---
col_annotation <- HeatmapAnnotation(
  Garden = garden_ids,
  col = list(Garden = garden_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# --- Custom function to visualize NA values ---
na_cell_fun <- function(j, i, x, y, width, height, fill) {
  if (is.na(heatmap_data[i, j])) {
    grid.rect(x, y, width, height, gp = gpar(fill = "white", col = NA))
    grid.lines(
      x = unit(c(x - width / 2, x + width / 2), "npc"),
      y = unit(c(y - height / 2, y + height / 2), "npc"),
      gp = gpar(col = "black", lwd = 0.5)
    )
  } else {
    grid.rect(x, y, width, height, gp = gpar(fill = fill, col = NA))
  }
}

# --- Perform hierarchical clustering on rows ---
row_clusters <- cutree(hclust(dist(heatmap_data)), k = 3)

# --- Create main heatmap ---
main_heatmap <- Heatmap(
  heatmap_data,
  name = "Z-score",
  col = colores,
  show_heatmap_legend = FALSE,
  top_annotation = col_annotation,
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_title = NULL,
  na_col = "white",
  cell_fun = na_cell_fun,
  row_split = row_clusters,
  border = TRUE,
  border_gp = gpar(col = "black"),
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  cluster_columns = TRUE
)

# ============================================================================
# SECTION 10: PREPARE CORRELATION DATA WITH FORCED ROW ORDERING
# ============================================================================

# --- Reorder correlation matrix to match phenotype order ---
names_mets <- colnames(combined_data_numeric)[3:(ncol(combined_data_numeric) - 4)]
df <- correlation_matrix
df_reordered <- df[names_mets, ]
df_reordered

# --- Define pretty phenotypic trait names in the correct order ---
pretty_names <- c(
  "Mean leaf number",              # Leaves_mean
  "Leaf number change",            # Leaves_increase
  "Mean rosette number",           # Rosettes_mean
  "Rosette number change",         # Rosettes_increase
  "Mean plant volume",             # Volume_mean
  "Plant volume change",           # Volume_increase
  "Total ripe fruits",             # Total_ripe_fruits
  "Total runners",                 # Total_runners
  "Mean fruit volume",             # Fruit_volume_mean
  "Fertilized achenes (%)",        # Fertilized_seeds_mean
  "Mean leaf damage (%)",          # Leaf_damage_mean
  "Leaf damage (CV)",              # Leaf_damage_CV
  "Max leaf damage (%)"            # Leaf_damage_max
)

if (length(names_mets) != length(pretty_names)) {
  stop(sprintf("Mismatch: %d traits vs %d names", length(names_mets), length(pretty_names)))
}

# --- Apply pretty names to reordered correlation data ---
rownames(df_reordered) <- pretty_names

# --- Calculate circle sizes based on correlation magnitude ---
max_corr <- max(abs(correlation_matrix))
circle_scale <- abs(correlation_matrix) / max_corr
fixed_max_radius <- 0.125

# --- Define color palette for correlation heatmap ---
col_inverted_2 <- colorRamp2(c(-1, 0, 1), c("#3b1b79", "ivory", "#204f16"))

# ============================================================================
# SECTION 11: CREATE CORRELATION HEATMAP WITH CIRCLE SCALING AND FORCED ORDERING
# ============================================================================

# --- Custom function to draw circles proportional to correlation strength ---
correlation_cell_fun <- function(j, i, x, y, width, height, fill) {
  # Draw cell border
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA))

  # Calculate circle radius based on correlation magnitude
  radius <- circle_scale[i, j] * fixed_max_radius

  # Draw circle with color corresponding to correlation value
  grid.circle(
    x = x, y = y, r = radius,
    gp = gpar(fill = col_inverted_2(correlation_matrix[i, j]), col = NA)
  )
}

# --- Create correlation heatmap with same row_clusters as main heatmap ---
correlation_heatmap <- Heatmap(
  df_reordered,
  name = "Correlation",
  col = col_inverted_2,
  column_names_side = "top",
  row_names_rot = 45,
  row_names_gp = gpar(fontsize = 15),
  show_heatmap_legend = FALSE,
  rect_gp = gpar(type = "none"),
  cell_fun = correlation_cell_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_split = row_clusters
)

# ============================================================================
# SECTION 12: CREATE LEGENDS AND COMBINE HEATMAPS
# ============================================================================

# --- Create legend for Z-score heatmap ---
lgd_zscore <- Legend(
  title = "Z-score",
  col_fun = colores,
  direction = "horizontal"
)

# --- Create legend for correlation heatmap ---
lgd_cor <- Legend(
  title = "Correlation",
  col_fun = col_inverted_2,
  direction = "horizontal"
)

# --- Create legend for garden annotations ---
lgd_garden <- Legend(
  title = "Common Garden",
  labels = names(garden_colors),
  legend_gp = gpar(fill = garden_colors),
  direction = "horizontal",
  nrow = 1
)

# --- Combine all legends ---
packed_legends <- packLegend(
  lgd_zscore, lgd_cor, lgd_garden,
  direction = "horizontal",
  gap = unit(8, "mm")
)

# ============================================================================
# SECTION 13: COMBINE AND VISUALIZE FINAL FIGURE
# ============================================================================

# --- Combine main and correlation heatmaps side by side ---
combined_heatmap <- main_heatmap + correlation_heatmap

# --- Draw the complete figure with legends ---
draw(
  combined_heatmap,
  annotation_legend_list = list(packed_legends),
  annotation_legend_side = "bottom"
)

# ============================================================================
# END OF SCRIPT
# ============================================================================
