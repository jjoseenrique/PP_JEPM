# ============================================================================
# COMPLETE ICC-LATITUDE REGRESSION ANALYSIS PIPELINE
# ============================================================================
# From data import through ICC calculation to regression plots with BH correction
# Generates plasticity vs latitude analysis with statistical annotation
# Generated for scientific publication
# ============================================================================

# --- LOAD REQUIRED LIBRARIES ---
library(dplyr)
library(xlsx)
library(lme4)
library(ggplot2)
library(cowplot)
library(ggtext)

# ============================================================================
# SECTION 1: DATA IMPORT
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
# SECTION 2: PREPARE SUBSETS
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
# SECTION 3: IDENTIFY COMMON SAMPLES
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
# SECTION 4: COMBINE DATA
# ============================================================================

# --- Merge metabolites and phenotypes ---
combined_df <- cbind(metab_subset[, (37:38)], metab_subset[, -(37:38)], pheno_subset[, -(14:15)])
rownames(combined_df) <- common_samples

# ============================================================================
# SECTION 5: LOG2 TRANSFORMATION
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
# SECTION 6: Z-SCORE SCALING
# ============================================================================

# --- Apply z-score scaling ---
combined_df <- combined_log2 %>%
  mutate_at(vars(-(1:2)), scale)

# ============================================================================
# SECTION 7: PREPARE DATA FOR ICC CALCULATION
# ============================================================================

# --- Store combined data for ICC calculation ---
filtered_traits <- combined_df

# --- Extract genotypes ---
genos <- unique(filtered_traits$GENOID)

# ============================================================================
# SECTION 8: CALCULATE ICC (INTRACLASS CORRELATION COEFFICIENT)
# ============================================================================

# --- Create subsets by genotype and calculate ICC for each trait ---
lista_traits <- list()

for (trait_col in colnames(filtered_traits)[-(1:2)]) {
  ICCs_vector <- NULL
  
  for (geno in genos) {
    # --- Subset data for specific genotype ---
    geno_data <- as.data.frame(subset(filtered_traits, GENOID == geno))
    
    # --- Fit linear mixed model ---
    formula <- as.formula(paste(trait_col, "~ 1 + (1|GARDENID)"))
    modelo <- lmer(formula, data = geno_data, REML = TRUE)
    
    # --- Extract variance components ---
    varcomp <- as.data.frame(VarCorr(modelo))
    
    # --- Calculate ICC ---
    ICC <- varcomp[1, "vcov"] / (varcomp[1, "vcov"] + varcomp[2, "vcov"])
    
    ICCs_vector <- c(ICCs_vector, ICC)
    names(ICCs_vector)[length(ICCs_vector)] <- paste0("subset_", geno)
  }
  
  lista_traits[[trait_col]] <- ICCs_vector
}

# ============================================================================
# SECTION 9: DEFINE PRETTY NAME MAPPINGS
# ============================================================================

# --- Pretty names for metabolites ---
pretty_metab_names <- c(
  "Pyruvic acid", "Valine", "Isoleucine", "Glycine", "Phosphoric acid", 
  "Proline", "Urea", "Glyceric acid", "Alanine", "Serine", "Succinic acid", 
  "Threonine", "Fumaric acid", "Nicotinic acid", "Erythritol", "Malic acid", 
  "GABA", "Aspartic acid", "Threonic acid", "Xylose", "Xylitol", 
  "Pyroglutamic acid", "Glutamic acid", "Arginine", "Tryptophan", 
  "Trehalose", "Maltitol", "Galactinol", "Quinic acid", "Fructose", 
  "Glucose", "Citric acid", "Methyl α-D-glucopyranoside", 
  "Dehydroascorbic acid", "myo-Inositol", "Sucrose"
)

# --- Pretty names for phenotypes ---
pretty_pheno_names <- c(
  "Mean leaf number", "Leaf number change", "Mean rosette number", 
  "Rosette number change", "Mean plant volume", "Plant volume change", 
  "Total ripe fruits", "Total runners", "Mean fruit volume", 
  "Fertilized achenes (%)", "Mean leaf damage (%)", "Leaf damage (CV)", 
  "Max leaf damage (%)"
)

# --- Create name mapper dictionary ---
ugly_names <- names(lista_traits)
all_pretty_names <- c(pretty_metab_names, pretty_pheno_names)

if (length(all_pretty_names) == length(ugly_names)) {
  name_mapper <- setNames(all_pretty_names, ugly_names)
} else {
  stop("Number of pretty names does not match number of traits")
}

# ============================================================================
# SECTION 10: PREPARE DATA FOR REGRESSION ANALYSIS
# ============================================================================

# --- Genotype coordinates ---
latitudes <- data.frame(
  Genotype = c("G01","G02","G03","G04","G05","G06","G07","G09","G10","G11","G12","G13","G14","G15","G16"),
  Latitude = c(37.7796, 40.2938, 43.1344, 45.94, 47.96666667, 48.801407, 50.015654, 
               54.5729, 55.5703, 60.1018, 60.2292, 60.4051, 70.1671, 70.0301, 70.0324)
)

# --- Convert lista_traits to data frame ---
icc_df <- as.data.frame(do.call(cbind, lista_traits))
rownames(icc_df) <- gsub("subset_", "", names(lista_traits[[1]]))

# --- Combine with latitude data ---
temp_data <- cbind(Latitude = latitudes$Latitude, icc_df)

# ============================================================================
# SECTION 11: REGRESSION ANALYSIS WITH BH CORRECTION
# ============================================================================

plot_list <- list()
p_values_raw <- numeric()
model_summaries <- list()
trait_vars <- colnames(temp_data)[-1]

# --- Loop 1: Collect p-values and create base plots ---
for (i in trait_vars) {
  formula <- as.formula(paste0("`", i, "` ~ Latitude"))
  mreg <- lm(formula, data = temp_data)
  s_mreg <- summary(mreg)
  
  p_values_raw[i] <- s_mreg$coefficients[2, 4]
  model_summaries[[i]] <- list(r_squared = s_mreg$r.squared)
  
  p <- ggplot(data = temp_data, aes(x = Latitude, y = .data[[i]])) +
    geom_point(color = "turquoise3", size = 2) +
    stat_smooth(method = "lm", formula = y ~ x, color = "black", linetype = "dashed", linewidth = 0.8) +
    labs(x = "", y = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          plot.title = element_markdown(hjust = 0.5, size = 12),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
  
  plot_list[[i]] <- p
}

# --- Apply Benjamini-Hochberg correction ---
p_adjusted <- p.adjust(p_values_raw, method = "BH")

# ============================================================================
# SECTION 12: ADD TITLES WITH STATISTICS AND PRETTY NAMES
# ============================================================================

# --- Loop 2: Update plot titles with statistics and pretty names ---
for (i in trait_vars) {
  r_squared <- round(model_summaries[[i]]$r_squared, 3)
  p_val_adj <- p_adjusted[i]
  
  p_formatted <- if (p_val_adj <= 0.01) {
    sprintf("%.2e", p_val_adj)
  } else {
    round(p_val_adj, 3)
  }
  
  # --- Get pretty name from mapper ---
  trait_name_pretty <- name_mapper[i]
  
  # --- Format title with subscript and conditional bolding ---
  title_text <- ifelse(p_val_adj < 0.05,
                       paste0("<b>", trait_name_pretty, "  R² = ", r_squared, ",  p-value<sub>BH</sub> = ", p_formatted, "</b>"),
                       paste0(trait_name_pretty, "  R² = ", r_squared, ",  p-value<sub>BH</sub> = ", p_formatted))
  
  plot_list[[i]] <- plot_list[[i]] + ggtitle(title_text)
}

# ============================================================================
# SECTION 13: ORDER PLOTS BY ADJUSTED P-VALUE
# ============================================================================

sorted_vars <- names(sort(p_adjusted, decreasing = TRUE))
plot_list_ordered <- plot_list[sorted_vars]

# ============================================================================
# SECTION 14: ADD LATITUDE LABELS TO BOTTOM ROW
# ============================================================================

num_plots <- length(plot_list_ordered)
num_cols <- 5
num_rows <- ceiling(num_plots / num_cols)
start_of_last_row <- (num_rows - 1) * num_cols + 1

for (k in seq_along(plot_list_ordered)) {
  if (k >= start_of_last_row) {
    plot_list_ordered[[k]] <- plot_list_ordered[[k]] + labs(x = "Latitude (°N)")
  }
}

# ============================================================================
# SECTION 15: CREATE AND DISPLAY COMPLETE GRID
# ============================================================================

final_grid <- plot_grid(plotlist = plot_list_ordered, ncol = 5)

print(final_grid)

# ============================================================================
# SECTION 16: SELECT AND DISPLAY SPECIFIC TRAITS
# ============================================================================

# --- Define traits to select (by original column names) ---
trait_row1 <- c("Glucose", "Fructose", "Leaves_mean")
trait_row2 <- c("Phosphoric_acid", "Sucrose", "Leaf_damage_CV")
trait_row3 <- c("Succinic_acid", "Citric_acid", "Volume_mean")

all_selected_traits <- c(trait_row1, trait_row2, trait_row3)

# --- Verify all selected traits exist ---
missing_traits <- all_selected_traits[!all_selected_traits %in% names(plot_list)]

if (length(missing_traits) > 0) {
  warning("The following traits were not found: ", paste(missing_traits, collapse = ", "))
  # Remove missing traits from selection
  trait_row1 <- trait_row1[trait_row1 %in% names(plot_list)]
  trait_row2 <- trait_row2[trait_row2 %in% names(plot_list)]
  trait_row3 <- trait_row3[trait_row3 %in% names(plot_list)]
}

# ============================================================================
# SECTION 18: EXTRACT SELECTED PLOT OBJECTS
# ============================================================================

# --- Extract plots for each row ---
plots_row1 <- lapply(trait_row1, function(trait) plot_list[[trait]])
plots_row2 <- lapply(trait_row2, function(trait) plot_list[[trait]])
plots_row3 <- lapply(trait_row3, function(trait) plot_list[[trait]])

# --- Add latitude label to last row ---
for (i in 1:length(plots_row3)) {
  plots_row3[[i]] <- plots_row3[[i]] + labs(x = "Latitude (°N)")
}

# ============================================================================
# SECTION 19: CREATE SELECTED TRAITS GRID
# ============================================================================
final_selection_grid <- plot_grid(
  plotlist = c(plots_row1, plots_row2, plots_row3),
  ncol = 3,
  nrow = 3
)

print(final_selection_grid)

# ============================================================================
# SECTION 21: EXPORT STATISTICS TABLE
# ============================================================================

# --- Create summary table ---
stats_table <- data.frame(
  Trait = names(p_values_raw),
  TraitPretty = sapply(names(p_values_raw), function(x) name_mapper[x]),
  P_Value_Raw = round(p_values_raw, 4),
  P_Value_BH = round(p_adjusted, 4),
  Significant = ifelse(p_adjusted < 0.05, "Yes", "No")
)

stats_table <- stats_table[order(stats_table$P_Value_BH), ]

print(stats_table)

# ============================================================================
# END OF SCRIPT
# ============================================================================