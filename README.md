# Reproducible Analysis Code

## Overview

This repository contains all R scripts and documentation necessary to reproduce the statistical analyses described in our manuscript on phenotypic plasticity and local adaptation in *Fragaria vesca* along a latitudinal gradient in Europe.

## Repository Structure

```
├── README.md                                    # This file
├── data/
│   ├── Table_S2.xlsx                          # Raw metabolomic data (36 metabolites)
│   └── Table_S3.xlsx                          # Raw phenotypic data (13 traits)
├── scripts/
    ├── 1.-metabolite_analysis_NP.R            # Metabolite data processing and heatmap visualization
    ├── 2.-phenotype_analysis_NP.R             # Phenotype data processing and heatmap visualization
    ├── 3.-plsda_biplot_NP.R                   # PLS-DA multivariate analysis and biplots
    ├── 4.-permanova_NP.R                      # PERMANOVA: garden and genotype effects
    ├── 5.-ICC-PCA-DBRDA_NP.R                  # ICC calculation, PCA, and db-RDA analysis
    ├── 6.-ICC_latitude_regression_NP.R        # Regression: ICC vs latitude with Benjamini-Hochberg correction
    ├── 7.-PERMANOVA-DBRDA_Plasticity_NP.R    # Variance partitioning: climatic vs garden effects
    └── README_Scripts.md                      # Detailed script descriptions

```

## Requirements

### R Version
- R ≥ 4.3.1

### Required Packages

```r
# Install required packages
required_packages <- c(
  "dplyr", "xlsx", "lme4", "vegan", "mixOmics", "FactoMineR", 
  "factoextra", "cowplot", "ggplot2", "ggtext", "geodata", 
  "terra", "raster", "reshape2", "VIM", "ComplexHeatmap", 
  "RColorBrewer", "nasapower", "tidyverse", "circlize", 
  "patchwork", "ggrepel"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}
```

## Data Files

### Table_S2.xlsx (Metabolomic Data)
- **Dimensions:** 581 samples × 36 metabolites + metadata
- **Columns:** GARDENID, GENOID, BLOCK, [36 metabolite columns (columns 7-42)]
- **Metabolites:** Primary metabolites identified via GC-MS profiling
- **Metabolite classes:** Organic acids, amino acids, sugars, sugar alcohols
- **Processing:** Log₂-transformed, z-scored, k-NN imputed (k=5)
- **Format:** Excel worksheet 1

### Table_S3.xlsx (Performance-related Data)
- **Dimensions:** 581 samples × 13 performance-related traits + metadata
- **Columns:** GARDENID, GENOID, BLOCK, [13 phenotype columns (columns 7-19)]
- **Traits (categories):**
  - **Morphological (6):** Leaves_mean, Leaves_increase, Rosettes_mean, Rosettes_increase, Volume_mean, Volume_increase
  - **Reproductive (4):** Total_ripe_fruits, Total_runners, Fruit_volume_mean, Fertilized_seeds_mean
  - **Defense/Herbivory (3):** Leaf_damage_mean, Leaf_damage_CV, Leaf_damage_max
- **Processing:** Log₂-transformed, z-scored, k-NN imputed (k=5)
- **Format:** Excel worksheet 1

Both datasets are filtered to **common samples** (N = 581) with complete data across all analyses.

## Scripts Overview

### 1. `metabolite_analysis_NP.R`
**Purpose:** Data import, normalization, and exploratory analysis of metabolite data (Table_S2)

**Key Steps:**
- Data import and Shapiro-Wilk normality testing
- Log₂ transformation and z-score scaling
- Mean calculation by GARDEN_ID and GENOID
- Heatmap visualization grouped by common garden
- Correlation with climatic data (monthly sampling period)
- NASA POWER API data extraction for 5 gardens (GON_21, GON_22, ALN_21, RUI_21, KEV_22)

**Outputs:**
- Normalized metabolite matrix
- Heatmaps with z-scores and correlation values
- Climatic_info.xlsx with garden-specific climate data

**Key Functions:**
- `ComplexHeatmap::Heatmap()` for visualization
- `nasapower::get_power()` for climate data retrieval
- `RColorBrewer::colorRamp2()` for color scaling

---

### 2. `phenotype_analysis_NP.R`
**Purpose:** Data import, normalization, and exploratory analysis of phenotypic data (Table_S3)

**Key Steps:**
- Data import and Shapiro-Wilk normality testing
- Log₂ transformation and z-score scaling
- Mean calculation by GARDEN_ID and GENOID
- Heatmap visualization grouped by common garden
- Correlation with climatic data (annual averages)
- NASA POWER API data extraction for climate variables

**Outputs:**
- Normalized phenotype matrix
- Heatmaps with z-scores and correlation values
- Integration with climate data from 4 gardens

**Note:** Very similar structure to script 1 but applied to phenotypic data instead of metabolites

---

### 3. `plsda_biplot_NP.R`
**Purpose:** Discriminate sample groups by garden-year combinations using multivariate analysis

**Key Steps:**
- Combined data import (Table_S2 + Table_S3) with identical preprocessing to PERMANOVA
- PLS-DA with 3 latent variables (LV1, LV2, LV3)
- Extract top 15 variables by loading magnitude for each biplot
- Generate biplots: LV1 vs LV2 and LV1 vs LV3
- Add 95% confidence ellipses by garden-year group
- Apply pretty name mappings to metabolites and phenotypes

**Outputs:**
- PLS-DA biplots with sample scores and variable loadings
- Explained variance for each latent variable
- Sample grouping visualization by GARDEN_ID

**Key Functions:**
- `mixOmics::plsda()` with 3 components
- `ggplot2::stat_ellipse()` for confidence regions
- `ggrepel::geom_text_repel()` for non-overlapping labels

---

### 4. `permanova_NP.R`
**Purpose:** Test for significant effects of GARDEN_ID, genotype, and their interaction on combined metabolite + phenotype data

**Models Tested:**
1. Sequential: GARDENID + GENOID + GARDENID:GENOID (Type I SS)
2. Reversed: GENOID + GARDENID + GENOID:GARDENID (to check order effects)

**Outputs:**
- **Table 1:** PERMANOVA results (R², F-statistics, p-values)
- Variance partitioning for all 49 traits combined
- Comparison of sequential vs reversed models

**Key Functions:**
- `vegan::adonis2()` with by = "term"
- Euclidean distances on combined metabolite + phenotype matrix
- 9,999 permutations

**Note:** Simplest script - no imputation, just log₂ + z-score on combined data

---

### 5. `ICC-PCA-DBRDA_NP.R`
**Purpose:** Quantify trait-specific plasticity via ICC, perform dimensionality reduction, and test latitudinal structuring

**Key Steps:**

#### Part A: ICC Calculation
- Fit linear mixed models for each trait within each genotype:
  ```
  Trait ~ 1 + (1|GARDENID)
  ```
- Calculate ICC as: var(GARDENID) / [var(GARDENID) + var(residual)]
- Produces ICC matrix (15 genotypes × 49 traits)

**ICC Interpretation:**
- ICC = 0: No plasticity (trait identical across gardens)
- ICC = 1: Maximum plasticity (all variation is environmental)
- Values represent proportion of trait variance attributable to environment

#### Part B: PCA on ICC Values
- Perform PCA on full ICC matrix (scale.unit = FALSE)
- Select top 20 traits by PC1 loading magnitude
- Generate biplots: All traits + Top 20 traits
- Map trait IDs to biological names

#### Part C: db-RDA Analysis (ICC vs Latitude)
- Calculate Euclidean distance matrix from ICC values
- Test: ICC distance matrix ~ Latitude
- Model: `capscale(dist_icc ~ Latitude, distance = "euclidean")`
- Adjusted R², F-statistic, p-value (9,999 permutations)

**Outputs:**
- ICC matrix (15 genotypes × 49 traits)
- PCA biplots and explained variance
- db-RDA results with latitude structuring test
- Trait legend (ID ↔ original name ↔ pretty name)

**Key Functions:**
- `lme4::lmer()` for mixed model fitting
- `lme4::VarCorr()` for variance component extraction
- `FactoMineR::PCA()` and `factoextra::fviz_pca_biplot()`
- `vegan::vegdist()`, `vegan::capscale()`, `vegan::anova.cca()`

---

### 6. `ICC_latitude_regression_NP.R`
**Purpose:** Test for latitudinal structuring of individual trait plasticity

**Key Steps:**
- Linear regression for each of 49 traits: ICC ~ Latitude
- Calculate Pearson correlation (r) and p-values
- Apply Benjamini-Hochberg correction for multiple testing
- Generate scatter plots with regression lines, statistics, confidence intervals
- Display results ordered by adjusted p-value
- Create full grid of all 49 traits and selection of key traits

**Outputs:**
- Scatter plots for all 49 traits (grid visualization)
- Selection grid (3×3) with key traits
- Summary statistics table (trait, r, p_raw, p_BH, significance)
- Figure highlighting significant traits (p_BH < 0.05)

**Key Statistics:**
- Raw p-values from individual regressions
- Benjamini-Hochberg corrected p-values
- R² and Pearson correlation coefficients
- F-statistics for each trait

**Interpretation:**
- Positive r: Plasticity increases toward northern (higher latitude) genotypes
- Negative r: Plasticity decreases toward northern genotypes
- p_BH < 0.05: Significant latitudinal structuring after multiple testing correction

---

### 7. `PERMANOVA-DBRDA_Plasticity_NP.R`
**Purpose:** Partition phenotypic variance by immediate environment vs. climatic origin using separate trait groups

**Design:**
- Separate analyses for 4 trait groups:
  1. Metabolic (36 metabolites)
  2. Biomass (6 morphological traits)
  3. Reproduction (4 reproductive traits)
  4. Herbivore_Damage (3 damage-related traits)

- Two bioclimatic variables tested separately:
  - BIO_1: Annual Mean Temperature (from WorldClim)
  - BIO_12: Annual Precipitation (from WorldClim)

**db-RDA Models (for each trait group × bioclimatic variable):**

1. **Model 1 (Environment effect):**
   ```
   capscale(traits ~ GARDENID, distance = "euclidean")
   ```
   R² = variance explained by immediate garden environment

2. **Model 2 (Climate effect, conditional - MAIN RESULT):**
   ```
   capscale(traits ~ BIO + Condition(GARDENID), distance = "euclidean")
   ```
   R² = variance explained by climate AFTER removing garden effects
   Isolates adaptive/heritable signal independent of immediate environment

3. **Model 3 (Full model, for verification):**
   ```
   capscale(traits ~ GARDENID + BIO, distance = "euclidean")
   ```
   R² = total variance explained by both factors

**Outputs:**
- **Table 2 (PRIMARY - db-RDA conditional):**
  - R² (%) for Climatic Covariate (Model 2)
  - R² (%) for Garden Environment (Model 1)
  - F-statistics and p-values from permutation tests (9,999 permutations)
  
- **Supplementary Table:** PERMANOVA marginal results (8 trait groups × 2 bioclimate vars) for robustness verification
  - Note: Values differ slightly from db-RDA due to different statistical approaches

**Key Statistics:**
- Adjusted R² (accounts for number of predictors using `RsquareAdj()`)
- F-statistics from 9,999 permutation tests
- p-values (< 0.0001 highly significant; formatted as "< 10^-4 ***")

**Interpretation:**
- R² Garden: Magnitude of phenotypic plasticity (immediate environment effect)
- R² Climate: Magnitude of local adaptation signal (climatic origin effect)
- Ratio: Relative importance of plasticity vs. adaptation
- Negative R² Climate: Indicates no explanatory power (noise/type I error)

**Data Processing:**
- 581 samples × 49 traits combined
- Log₂ transformation with offset handling
- Z-score scaling
- k-NN imputation (k=5)
- WorldClim v2.1 bioclimatic variables at 2.5 arc-minute resolution

---

## Execution Instructions

### Step 1: Prepare Your Environment

```r
# Set working directory to repo folder
setwd("path/to/repository")

# Verify data files exist in data/ subdirectory
file.exists("data/Table_S2.xlsx")
file.exists("data/Table_S3.xlsx")
```

### Step 2: Run Scripts Sequentially

```r
# Execute scripts in numerical order
# Scripts 1-3 are exploratory and independent
# Scripts 4-7 build on common data preprocessing but can run independently

source("scripts/1.-metabolite_analysis_NP.R")      
source("scripts/2.-phenotype_analysis_NP.R")       
source("scripts/3.-plsda_biplot_NP.R")             
source("scripts/4.-permanova_NP.R")                
source("scripts/5.-ICC-PCA-DBRDA_NP.R")            
source("scripts/6.-ICC_latitude_regression_NP.R")  
source("scripts/7.-PERMANOVA-DBRDA_Plasticity_NP.R") 
```

## Reproducibility Notes

- **Random Seed:** Set to 20 across all scripts using `set.seed(20)`
- **Permutations:** 9,999 permutations for all permutation tests (scripts 4, 5, 7)
- **Correction Method:** Benjamini-Hochberg for multiple testing (script 6)
- **Distance Metric:** Euclidean distances on log₂-transformed, z-scored data
- **Imputation:** k-NN with k=5 for missing values (scripts 5, 7)
- **Bioclimatic Data:** WorldClim v2.1 at 2.5 arc-minute resolution via `geodata::worldclim_tile()`

## Key Statistical Methods

### Variance Partitioning (db-RDA with Conditioning)
- **Conditional ordination** removes immediate environment effects
- Isolates climate-driven (adaptive/heritable) signals from plastic responses
- Uses Euclidean distances on log₂-transformed, z-scored traits
- More conservative than marginal approaches (PERMANOVA by = "margin")
- Implemented with `vegan::capscale()` and `Condition()` syntax

### ICC (Intraclass Correlation Coefficient)
- **Linear mixed model:** Trait ~ 1 + (1|GARDENID)
- **Calculation:** ICC = var(GARDENID) / total_variance
- **Interpretation:** Proportion of trait variance due to environment (plasticity)
- **Genotype-specific:** Calculated separately for each of 15 genotypes
- Provides genotype-by-environment interaction quantification

### Multiple Testing Correction
- **Benjamini-Hochberg:** Controls false discovery rate at α = 0.05
- Applied when performing multiple independent tests (49 trait-by-trait regressions in script 6)
- Uses `p.adjust(..., method = "BH")` from base R `stats` package

## Troubleshooting

### Error: "Table_S2.xlsx not found" or "Table_S3.xlsx not found"
- Verify file paths are correct relative to working directory
- Ensure filenames match exactly (case-sensitive on Linux/Mac)
- Files should be in `data/` subdirectory

### Error in NASA POWER data download (scripts 1-2)
- Requires internet connection
- Some dates/locations may have missing data
- If download fails, manually create `Climatic_info.xlsx` or skip climate analysis

### Error: "Object 'X' not found" or model fitting issues (scripts 5-7)
- Ensure previous scripts executed successfully
- Check for NA values in data after imputation
- Verify lme4 and vegan packages are updated

### Memory Issues or Slow Execution
- Current analyses use ~581 samples × ~50 traits; typical memory < 1 GB
- Permutation tests (9,999 iterations) are computationally intensive
- If execution is slow, reduce permutations from 9,999 to 999 (for testing only)

## Citation

If you use these scripts, please cite:

> Pérez-Martín JE, Batsleer F, Vandegehuchte ML, et al. (Year). Phenotypic plasticity and local adaptation in *Fragaria vesca* along a European latitudinal gradient. *New Phytologist*, XX(X), XXX-XXX.

And reference this repository and data archive:

> Reproducible code and data available at: https://github.com/jjoseenrique/NP_JEPM
> 
> Raw data deposited in Zenodo: https://doi.org/[ZENODO_DOI]

## Contact & Support

For questions regarding these analyses, contact:
- **Primary Contact:** José E. Pérez-Martín (jepm@uma.es)
- **Corresponding Authors:** Sonia Osorio & David Posé
- **Supervision:** Dries Bonte, Femke Batsleer & Martijn L. Vandegehuchte

## License

This code is provided under the [MIT License](LICENSE)

---

**Last Updated:** November 2025  
**R Version Used:** 4.3.1  
**Status:** Complete and tested  
**Total Scripts:** 7  
**Total Data Files:** 2 (Table_S2.xlsx, Table_S3.xlsx)  
**Reproducibility:** Seeds, permutations, and corrections all specified
