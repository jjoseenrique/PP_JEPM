# ============================================================================
# COMPLETE ICC-PCA ANALYSIS PIPELINE WITH db-RDA AND PRETTY NAMES
# ============================================================================
# From data import through ICC calculation to PCA visualization and db-RDA
# Generated for scientific publication
# ============================================================================

# --- LOAD REQUIRED LIBRARIES ---
library(dplyr)
library(xlsx)
library(lme4)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(ggplot2)
library(vegan)

# --- SET SEED FOR REPRODUCIBILITY ---
set.seed(20)

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
# SECTION 9: PREPARE DATA FOR PCA (MATCHING ORIGINAL EXACTLY)
# ============================================================================

# --- EXACTLY as in original script ---
lista_traits_df <- as.data.frame(do.call(cbind, lista_traits))
original_names <- names(lista_traits_df)

# --- Map each trait to a number (1 to n) ---
num_labels <- as.character(seq_along(original_names))
names(num_labels) <- original_names
colnames(lista_traits_df) <- num_labels[original_names]
rownames(lista_traits_df) <- gsub("subset_", "", rownames(lista_traits_df))

# ============================================================================
# SECTION 10: FULL PCA (ALL TRAITS) - scale.unit = FALSE
# ============================================================================

# --- Perform PCA without scaling ---
pca_full <- PCA(lista_traits_df, scale.unit = FALSE, graph = FALSE)
explained_full <- round(pca_full$eig[1:2, 2], 2)

# ============================================================================
# SECTION 11: SELECT TOP 20 TRAITS BY PC1 LOADING
# ============================================================================

# --- Extract loadings and select top 20 ---
loadings <- pca_full$var$coord
top20_nums <- rownames(loadings)[order(abs(loadings[, 1]), decreasing = TRUE)[1:20]]
top20_names <- original_names[as.integer(top20_nums)]

# --- Subset with numbered trait columns ---
df_top20 <- lista_traits_df[, top20_nums]

# ============================================================================
# SECTION 12: DEFINE PRETTY NAME MAPPINGS (MATCHING ORIGINAL)
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

# --- Function to convert original names to pretty names (MATCHING ORIGINAL) ---
get_pretty_name <- function(original_name) {
  # --- Metabolite mapping with UNDERSCORE LOWERCASE names (as in original) ---
  metabolite_mapping <- setNames(pretty_metab_names, 
                                 c("pyruvic_acid", "valine", "isoleucine", "glycine", "phosphoric_acid",
                                   "proline", "urea", "glyceric_acid", "alanine", "serine", "succinic_acid",
                                   "threonine", "fumaric_acid", "nicotinic_acid", "erythritol", "malic_acid",
                                   "gaba", "aspartic_acid", "threonic_acid", "xylose", "xylitol",
                                   "pyroglutamic_acid", "glutamic_acid", "arginine", "tryptophan",
                                   "trehalose", "maltitol", "galactinol", "quinic_acid", "fructose",
                                   "glucose", "citric_acid", "methyl_alpha_d_glucopyranoside",
                                   "dehydroascorbic_acid", "myo_inositol", "sucrose"))
  
  if (original_name %in% names(metabolite_mapping)) {
    return(metabolite_mapping[original_name])
  }
  
  # --- Phenotype mapping ---
  phenotype_mapping <- setNames(pretty_pheno_names,
                                c("mean_leaf_number", "leaf_number_change", "mean_rosette_number",
                                  "rosette_number_change", "mean_plant_volume", "plant_volume_change",
                                  "total_ripe_fruits", "total_runners", "mean_fruit_volume",
                                  "fertilized_achenes_percent", "mean_leaf_damage_percent", 
                                  "leaf_damage_cv", "max_leaf_damage_percent"))
  
  if (original_name %in% names(phenotype_mapping)) {
    return(phenotype_mapping[original_name])
  }
  
  # --- Fallback ---
  pretty_name <- gsub("_", " ", original_name)
  pretty_name <- tools::toTitleCase(pretty_name)
  return(pretty_name)
}

# --- Apply pretty names to top20 traits ---
top20_pretty_names <- sapply(top20_names, get_pretty_name)

# ============================================================================
# SECTION 13: DEFINE PCA THEME
# ============================================================================

# --- Consistent theme for both plots ---
tema_pca <- function() {
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    aspect.ratio = 1,
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 13),
    legend.position = "none"
  )
}

# ============================================================================
# SECTION 14: CREATE BIPLOT FOR ALL TRAITS
# ============================================================================

plot_all_traits <- fviz_pca_biplot(
  pca_full,
  repel = TRUE,
  labelsize = 5,
  geom.var = c("text", "point"),
  geom.ind = c("point", "text"),
  title = "",
  xlab = paste0("PC1 (", explained_full[1], "%)"),
  ylab = paste0("PC2 (", explained_full[2], "%)")
) +
  tema_pca()

# ============================================================================
# SECTION 15: PCA ON TOP 20 TRAITS WITH SCALE.UNIT = TRUE AND PRETTY NAMES
# ============================================================================

# --- IMPORTANT: Assign pretty names BEFORE PCA (as in original) ---
colnames(df_top20) <- top20_pretty_names

# --- Perform PCA WITH scaling (scale.unit = TRUE) as in original ---
pca_selected <- PCA(as.data.frame(df_top20), scale.unit = F, graph = FALSE)
explained_top20 <- round(pca_selected$eig[1:2, 2], 2)

# ============================================================================
# SECTION 16: CREATE BIPLOT FOR TOP 20 TRAITS
# ============================================================================

plot_selected_traits <- fviz_pca_biplot(
  pca_selected,
  labelsize = 5,
  repel = TRUE,
  geom.var = c("point", "text"),
  geom.ind = c("point", "text"),
  title = "",
  xlab = paste0("PC1 (", explained_top20[1], "%)"),
  ylab = paste0("PC2 (", explained_top20[2], "%)")
) +
  tema_pca()

# ============================================================================
# SECTION 17: COMBINE PLOTS
# ============================================================================

# --- Combine both PCA biplots with equal sizes ---
combined_pca <- plot_grid(
  plot_all_traits,
  plot_selected_traits,
  ncol = 2,
  align = "hv",
  axis = "lrtb"
)

# --- Display combined plot ---
print(combined_pca)

# ============================================================================
# SECTION 18: CREATE TRAIT LEGEND WITH PRETTY NAMES
# ============================================================================

# --- Create mapping between IDs, original names, and pretty names ---
legend_df <- data.frame(
  TraitID = seq_along(original_names),
  TraitName = original_names,
  PrettyName = sapply(original_names, get_pretty_name)
)

print(legend_df)

# ============================================================================
# SECTION 19: PREPARE DATA FOR db-RDA ANALYSIS
# ============================================================================

# --- Convert lista_traits to ICC data frame (genotypes as rows) ---
icc_df <- as.data.frame(do.call(cbind, lista_traits))
rownames(icc_df) <- gsub("subset_", "", names(lista_traits[[1]]))

# --- Genotype latitude coordinates ---
latitudes <- data.frame(
  Genotype = c("G01","G02","G03","G04","G05","G06","G07","G09","G10","G11","G12","G13","G14","G15","G16"),
  Latitude = c(37.7796, 40.2938, 43.1344, 45.94, 47.96666667, 48.801407, 50.015654, 
               54.5729, 55.5703, 60.1018, 60.2292, 60.4051, 70.1671, 70.0301, 70.0324)
)

# ============================================================================
# SECTION 20: db-RDA ANALYSIS - ICC PLASTICITY PROFILES vs LATITUDE
# ============================================================================

# --- Calculate distance matrix from ICC values ---
dist_icc <- vegdist(icc_df, method = "euclidean")

# --- Prepare data for db-RDA ---
latitude_vector <- latitudes$Latitude[match(rownames(icc_df), latitudes$Genotype)]
dbRDA_data <- data.frame(Latitude = latitude_vector)

# --- Run db-RDA: ICC matrix ~ Latitude ---
rda_result <- capscale(dist_icc ~ Latitude, data = dbRDA_data, distance = "euclidean")

# --- Test significance with permutations ---
anova_result <- anova.cca(rda_result, permutations = 9999)
print(anova_result)

# --- Extract statistics ---
r2_adj <- RsquareAdj(rda_result)$adj.r.squared
f_value <- anova_result$F[1]
p_value <- anova_result$`Pr(>F)`[1]

# ============================================================================
# END OF SCRIPT
# ============================================================================
