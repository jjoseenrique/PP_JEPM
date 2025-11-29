# ============================================================================
# VARIANCE PARTITIONING PIPELINE - CLIMATIC ORIGIN VS GARDEN ENVIRONMENT
# ============================================================================
# From data import through distance matrix calculation to sequential db-RDA
# Quantifies adaptive signals independent of immediate environment effects
# Generated for scientific publication
# ============================================================================

# --- LOAD REQUIRED LIBRARIES ---
library(mixOmics)
library(dplyr)
library(xlsx)
library(cowplot)
library(gridExtra)
library(lme4)
library(geodata)
library(ggtext)
library(vegan)
library(VIM)
library(terra)
library(raster)
library(ggplot2)
library(reshape2)
library(FactoMineR)
library(factoextra)

# --- SET SEED FOR REPRODUCIBILITY ---
set.seed(20)

# ============================================================================
# SECTION 1: DATA IMPORT
# ============================================================================

# --- Read Excel files ---
metab <- read.xlsx("Table_S2.xlsx", 1)
pheno <- read.xlsx("Table_S3.xlsx", 1)

# --- Rename phenotypic columns ---
colnames(pheno)[-(1:6)] <- c("Leaves_mean", "Leaves_increase", "Rosettes_mean", 
                             "Rosettes_increase", "Volume_mean", "Volume_increase", 
                             "Total_ripe_fruits", "Total_runners", "Fruit_volume_mean", 
                             "Fertilized_seeds_mean", "Leaf_damage_mean",
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
combined_df <- cbind(metab_subset[, (37:38)], metab_subset[, -(37:38)], 
                     pheno_subset[, -(14:15)])
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
# SECTION 7: KNN IMPUTATION
# ============================================================================

# --- Apply KNN imputation ---
traits_for_imputation <- combined_df[, -c(1, 2)]
imputed_traits <- kNN(traits_for_imputation, k = 5, imp_var = FALSE)
filtered_traits_imputed <- cbind(combined_df[, c(1, 2)], imputed_traits)

# ============================================================================
# SECTION 8: ADD BIOCLIMATIC VARIABLES
# ============================================================================

# --- Genotype coordinates ---
points_data <- data.frame(
  name = c("G01", "G02", "G03", "G04", "G05", "G06", "G07", "G09", "G10", 
           "G11", "G12", "G13", "G14", "G15", "G16"),
  lat = c(37.7796, 40.2938, 43.1344, 45.94, 47.96666667, 48.801407, 50.015654, 
          54.5729, 55.5703, 60.1018, 60.2292, 60.4051, 70.1671, 70.0301, 70.0324),
  lon = c(-3.7849, -5.0091, -4.888, 10.81, 7.83333333, 2.130122, 2.6973567, 
          24.6722, 9.7466, 18.3433333, 25.0233, 25.1562, 24.7561, 22.0653, 23.4012)
)

# --- Extract bioclimatic data ---
var <- c("bio")
valores <- list()

for (i in 1:15) {
  coords <- data.frame(x = points_data$lon[i], y = points_data$lat[i])
  
  for (n in var) {
    Clim_tiff <- worldclim_tile(var = n, lon = points_data$lon[i], 
                                lat = points_data$lat[i], path = tempdir())
    values <- terra::extract(Clim_tiff, coords)
    
    if (is.na(values[1, 2])) {
      Clim_tiff <- worldclim_tile(var = n, lon = points_data$lon[i], 
                                  lat = points_data$lat[i], path = tempdir())
      values <- terra::extract(Clim_tiff, coords)
    }
    
    valores[[paste0(points_data$name[i])]] <- values[, -1]
  }
}

# --- Format bioclimatic data ---
test1 <- data.frame()
for (i in valores) {
  test1 <- rbind(test1, unlist(i))
}

test1 <- cbind.data.frame(test1[, 1], test1[, (12:ncol(test1))], test1[, (2:11)])
colnames(test1) <- c(paste0("BIO_", 1:19))
rownames(test1) <- c(paste0("G0", 1:7), "G09", paste0("G1", 0:6))

# --- Add GENOID and merge ---
test1$GENOID <- rownames(test1)
full_data_with_bio <- merge(filtered_traits_imputed, 
                            test1[, c("GENOID", "BIO_1", "BIO_12")], 
                            by = "GENOID")

# ============================================================================
# SECTION 9: DEFINE TRAIT GROUPS (SEPARATED)
# ============================================================================

metabolic_traits <- c("Pyruvic_acid", "Valine", "Isoleucine", "Glycine", 
                      "Phosphoric_acid", "Proline", "Urea", "Glyceric_acid", 
                      "Alanine", "Serine", "Succinic_acid", "Threonine", 
                      "Fumaric_acid", "Nicotinic_acid", "Erythritol", "Malic_acid", 
                      "GABA", "Aspartic_acid", "Threonic_acid", "Xylose", "Xylitol", 
                      "Pyroglutamic_acid", "Glutamic_acid", "Arginine", "Tryptophan", 
                      "Trehalose", "Maltitol", "Galactinol", "Quinic_acid", 
                      "Fructose", "Glucose", "Citric_acid", 
                      "Methyl_α.D.glucopyranoside", "Dehydroascorbic_acid", 
                      "Inositol_myo", "Sucrose")

biomass_traits <- c("Leaves_mean", "Leaves_increase", "Rosettes_mean", 
                    "Rosettes_increase", "Volume_mean", "Volume_increase")

reproduction_traits <- c("Total_ripe_fruits", "Total_runners", "Fruit_volume_mean", 
                         "Fertilized_seeds_mean")

herbivore_damage_traits <- c("Leaf_damage_mean", "Leaf_damage_CV", "Leaf_damage_max")

# --- List with separated groups ---
list_of_trait_groups <- list(
  Metabolic = metabolic_traits,
  Biomass = biomass_traits,
  Reproduction = reproduction_traits,
  Herbivore_Damage = herbivore_damage_traits
)

# ============================================================================
# SECTION 10: PERMANOVA BY TRAIT GROUPS (MARGINAL EFFECTS FOR VERIFICATION)
# ============================================================================

predictors <- full_data_with_bio[, c("BIO_1", "BIO_12", "GARDENID", "GENOID")]
permanova_results_by_group <- list()

for (group_name in names(list_of_trait_groups)) {
  
  # --- Select traits for this group ---
  current_traits <- full_data_with_bio[, list_of_trait_groups[[group_name]]]
  
  # --- Calculate distance matrix ---
  dist_matrix_group <- vegdist(current_traits, method = "euclidean")
  
  # --- Temperature model (MARGINAL EFFECTS) ---
  permanova_temp_group <- adonis2(
    formula = dist_matrix_group ~ GARDENID + BIO_1,
    data = predictors,
    permutations = 9999,
    by = "margin"
  )
  
  # --- Precipitation model (MARGINAL EFFECTS) ---
  permanova_precip_group <- adonis2(
    formula = dist_matrix_group ~ GARDENID + BIO_12,
    data = predictors,
    permutations = 9999,
    by = "margin"
  )
  
  # --- Store results ---
  permanova_results_by_group[[group_name]] <- list(
    Temp_Model = permanova_temp_group,
    Precip_Model = permanova_precip_group
  )
}

# ============================================================================
# SECTION 11: db-RDA SEQUENTIAL ANALYSIS (PRIMARY RESULTS FOR TABLE 2)
# ============================================================================

# --- Function for db-RDA sequential analysis using capscale() ---
dbRDA_sequential <- function(group_name, group_traits, bio_var, bio_var_name) {
  
  # --- Prepare data ---
  traits_data <- full_data_with_bio[, group_traits]
  predictors_env <- full_data_with_bio[, c("GARDENID", bio_var)]
  
  # --- Model 1: GARDENID only ---
  rda_garden <- capscale(traits_data ~ GARDENID, 
                         data = predictors_env, 
                         distance = "euclidean")
  anova_garden <- anova.cca(rda_garden, permutations = 9999)
  r2_garden <- RsquareAdj(rda_garden)$adj.r.squared
  
  # --- Model 2: Bioclimatic variable conditioned on GARDENID ---
  # This is the KEY model: isolates variance independent of GARDENID
  rda_bioclim_cond <- capscale(traits_data ~ get(bio_var) + Condition(GARDENID), 
                               data = predictors_env, 
                               distance = "euclidean")
  anova_bioclim_cond <- anova.cca(rda_bioclim_cond, permutations = 9999)
  r2_bioclim <- RsquareAdj(rda_bioclim_cond)$adj.r.squared
  f_bioclim <- anova_bioclim_cond$F[1]
  p_bioclim <- anova_bioclim_cond$`Pr(>F)`[1]
  
  # --- Model 3: Full model (for verification) ---
  rda_full <- capscale(traits_data ~ GARDENID + get(bio_var), 
                       data = predictors_env, 
                       distance = "euclidean")
  anova_full <- anova.cca(rda_full, by = "terms", permutations = 9999)
  r2_full <- RsquareAdj(rda_full)$adj.r.squared
  
  # --- Return results ---
  return(list(
    group = group_name,
    bioclim_var = bio_var_name,
    r2_garden = r2_garden,
    r2_bioclim_cond = r2_bioclim,
    r2_full = r2_full,
    f_bioclim = f_bioclim,
    p_bioclim = p_bioclim,
    p_garden = anova_garden$`Pr(>F)`[1],
    f_garden = anova_garden$F[1]
  ))
}

# --- Execute db-RDA for all combinations ---
dbrda_results <- list()

for (group_name in names(list_of_trait_groups)) {
  
  # Temperature (BIO_1)
  result_temp <- dbRDA_sequential(
    group_name = group_name,
    group_traits = list_of_trait_groups[[group_name]],
    bio_var = "BIO_1",
    bio_var_name = "Annual Mean Temperature (BIO_1)"
  )
  dbrda_results[[paste(group_name, "BIO_1", sep = "_")]] <- result_temp
  
  # Precipitation (BIO_12)
  result_precip <- dbRDA_sequential(
    group_name = group_name,
    group_traits = list_of_trait_groups[[group_name]],
    bio_var = "BIO_12",
    bio_var_name = "Annual Precipitation (BIO_12)"
  )
  dbrda_results[[paste(group_name, "BIO_12", sep = "_")]] <- result_precip
}

# ============================================================================
# SECTION 12: CREATE PUBLICATION-READY SUMMARY TABLE (TABLE 2)
# ============================================================================
# This table reports db-RDA results: R² values + F-statistics + p-values
# from conditional models. These are the PRIMARY results for testing 
# climate-mediated adaptation

publication_table <- data.frame(
  `Trait Group` = character(),
  `Climatic Driver` = character(),
  `R² (%) (Climatic Covariate)` = numeric(),
  `R² (%) (Garden Environment)` = numeric(),
  `F-statistic` = numeric(),
  `P-value` = character(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# --- Loop through each trait group ---
for (group_name in names(list_of_trait_groups)) {
  
  # --- Temperature (BIO1) model ---
  result_temp_key <- paste(group_name, "BIO_1", sep = "_")
  result_temp <- dbrda_results[[result_temp_key]]
  
  # Extract conditional R² (from Model 2: BIO_1 | Condition(GARDENID))
  r2_bio1_cond <- result_temp$r2_bioclim_cond * 100
  # Extract GARDENID-only R² (from Model 1)
  r2_garden_temp <- result_temp$r2_garden * 100
  # Extract F-statistic and p-value
  f_bio1 <- result_temp$f_bioclim
  p_bio1 <- result_temp$p_bioclim
  
  # Format p-value
  p_bio1_formatted <- if (p_bio1 < 0.0001) {
    "< 10^-4 ***"
  } else if (p_bio1 < 0.001) {
    paste0(format(p_bio1, scientific = TRUE, digits = 2), " **")
  } else if (p_bio1 < 0.01) {
    paste0(round(p_bio1, 3), " *")
  } else {
    round(p_bio1, 3)
  }
  
  # --- Precipitation (BIO12) model ---
  result_precip_key <- paste(group_name, "BIO_12", sep = "_")
  result_precip <- dbrda_results[[result_precip_key]]
  
  # Extract conditional R² (from Model 2: BIO_12 | Condition(GARDENID))
  r2_bio12_cond <- result_precip$r2_bioclim_cond * 100
  # Extract GARDENID-only R² (from Model 1)
  r2_garden_precip <- result_precip$r2_garden * 100
  # Extract F-statistic and p-value
  f_bio12 <- result_precip$f_bioclim
  p_bio12 <- result_precip$p_bioclim
  
  # Format p-value
  p_bio12_formatted <- if (p_bio12 < 0.0001) {
    "< 10^-4 ***"
  } else if (p_bio12 < 0.001) {
    paste0(format(p_bio12, scientific = TRUE, digits = 2), " **")
  } else if (p_bio12 < 0.01) {
    paste0(round(p_bio12, 3), " *")
  } else {
    round(p_bio12, 3)
  }
  
  # --- Add Temperature row ---
  publication_table <- rbind(publication_table, data.frame(
    `Trait Group` = group_name,
    `Climatic Driver` = "Temperature (BIO1)",
    `R² (%) (Climatic Covariate)` = round(r2_bio1_cond, 2),
    `R² (%) (Garden Environment)` = round(r2_garden_temp, 2),
    `F-statistic` = round(f_bio1, 2),
    `P-value` = p_bio1_formatted,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  
  # --- Add Precipitation row ---
  publication_table <- rbind(publication_table, data.frame(
    `Trait Group` = group_name,
    `Climatic Driver` = "Precipitation (BIO12)",
    `R² (%) (Climatic Covariate)` = round(r2_bio12_cond, 2),
    `R² (%) (Garden Environment)` = round(r2_garden_precip, 2),
    `F-statistic` = round(f_bio12, 2),
    `P-value` = p_bio12_formatted,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
}

print(publication_table)

# ============================================================================
# END OF SCRIPT
# ============================================================================
