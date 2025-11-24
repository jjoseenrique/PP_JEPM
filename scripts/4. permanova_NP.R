# ============================================================================
# DATA PROCESSING TO PERMANOVA ANALYSIS
# ============================================================================
# From data import through normalization to multivariate PERMANOVA analysis
# Evaluates garden and genotype effects with interaction terms
# Generated for scientific publication
# ============================================================================

# --- LOAD REQUIRED LIBRARIES ---
library(dplyr)
library(xlsx)
library(vegan)

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
# SECTION 5: NORMALITY TEST
# ============================================================================

# --- Shapiro-Wilk test ---
shapiro <- apply(combined_df[, -(1:2)], 2, shapiro.test)
p_values <- sapply(shapiro, function(x) x$p.value)

# ============================================================================
# SECTION 6: LOG2 TRANSFORMATION
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
# SECTION 7: Z-SCORE SCALING
# ============================================================================

# --- Apply z-score scaling ---
combined_df <- combined_log2 %>%
  mutate_at(vars(-(1:2)), scale)

# ============================================================================
# SECTION 8: PREPARE DATA FOR PERMANOVA
# ============================================================================

# --- Use combined_df directly as filtered_traits (NO imputation) ---
filtered_traits <- combined_df

# ============================================================================
# SECTION 9: PERMANOVA - SEQUENTIAL (GARDENID first)
# ============================================================================

set.seed(20)

permanova_result_full <- adonis2(
  filtered_traits[, -c(1, 2)] ~ GARDENID * GENOID,
  data = filtered_traits,
  permutations = 9999,
  method = "euclidean",
  by = "term",
  na.rm = TRUE
)

print(permanova_result_full)

# ============================================================================
# SECTION 10: PERMANOVA - REVERSED (GENOID first, check if it changes results)
# ============================================================================

permanova_result_full_2 <- adonis2(
  filtered_traits[, -c(1, 2)] ~ GENOID * GARDENID,
  data = filtered_traits,
  permutations = 9999,
  method = "euclidean",
  by = "term",
  na.rm = TRUE
)

print(permanova_result_full_2)

# ============================================================================
# END OF SCRIPT
# ============================================================================