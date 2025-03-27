# CardiometabolicExWASMR: Exposome-Wide Association Study (ExWAS) and Mendelian Randomization (MR) Analysis Repository

We created this repository to contain the code for the project **"Biobank-scale exposome-wide risk factors in CAD & T2D: observational, predictive, and causal evidence."**


Access to UK Biobank data is granted by following the steps described at the UK Biobank website [https://www.ukbiobank.ac.uk/principles-of-access/]

This repository includes scripts for:

- **Exposome-Wide Association Studies (ExWAS)** to assess exposome-wide factors for Coronary Artery Disease (CAD) and Type 2 Diabetes (T2D).
- **Mendelian Randomization (MR)** to evaluate causal relationships between exposures and outcomes.


---

## 📂 Directory Overview
```
.
├── ExWAS/                      # Contains ExWAS scripts and helper functions
│   ├── Baseline_PEWAS_Cox_Functions_script.R         # Helper function for Cox models
│   ├── Baseline_PEWAS_Logistic_Functions_script.R    # Helper function for Logistic models
│   ├── Incident_analysis/
│   │   ├── binary/
│   │   │   ├── Cox_CAD_binary_exposure_ExWAS.R
│   │   │   ├── Cox_T2D_binary_exposure_ExWAS.R
│   │   │   ├── Incident_CAD_binary_exposure_Logistic_ExWAS.R
│   │   │   └── Incident_T2D_binary_exposure_Logistic_ExWAS.R
│   │   └── untransformed/
│   │       ├── Cox_CAD_untransformed_exposure_ExWAS.R
│   │       ├── Cox_T2D_untransformed_exposure_ExWAS.R
│   │       ├── Incident_CAD_untransformed_exposure_Logistic_ExWAS.R
│   │       └── Incident_T2D_untransformed_exposure_Logistic_ExWAS.R
│   ├── Prevalent_analysis/
│   │   ├── binary/
│   │   │   ├── Prevalent_CAD_binary_exposure_ExWAS.R
│   │   │   └── Prevalent_T2D_binary_exposure_ExWAS.R
│   │   └── untransformed/
│   │       ├── Prevalent_CAD_untransformed_exposure_ExWAS.R
│   │       └── Prevalent_T2D_untransformed_exposure_ExWAS.R
│   └── README.md
├── MR/                         # Contains MR analysis scripts and helper functions
│   ├── compute_mr.R
│   ├── CAD_all_exposures_bidirectional_MR_analysis.R
│   ├── t2d_all_exposures_bidirectional_MR_analysis.R
│   └── README.md
├── exposures_id_name_mapping_updated_10_06_2022.csv  # Input CSV with exposure IDs
└── results/                    # Directory for storing results
```

---

## 🚀 Getting Started

### 1. **Prerequisites**
Ensure that R is installed and that all required packages are available. The following R packages are required for both ExWAS and MR analyses:

```r
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,      # Data manipulation and visualization
  TwoSampleMR,    # Mendelian Randomization analysis
  data.table,     # Efficient data manipulation
  purrr,          # Functional programming
  MRInstruments,  # MR instruments for analysis
  ggrepel,        # Text and label repelling in ggplot2
  survival,       # Cox proportional hazards models
  glmnet,         # Regularized regression
  gdata,          # Data manipulation utilities
  RNOmni,         # Rank-based normalization
  biglm,          # Bounded memory linear regression
  broom,          # Convert models into tidy data frames
  pROC            # ROC curve analysis
)
```

---

### 2. **Cloning the Repository**
```bash
git clone https://github.com/your-username/CardiometabolicExWASMR.git
cd CardiometabolicExWASMR
```

---

## 📊 Exposome-Wide Association Study (ExWAS)

### 📂 ExWAS Directory Structure
```
ExWAS/
├── Baseline_PEWAS_Cox_Functions_script.R         # Helper function for Cox models
├── Baseline_PEWAS_Logistic_Functions_script.R    # Helper function for Logistic models
├── Incident_analysis/
│   ├── binary/
│   │   ├── Cox_CAD_binary_exposure_ExWAS.R
│   │   ├── Cox_T2D_binary_exposure_ExWAS.R
│   │   ├── Incident_CAD_binary_exposure_Logistic_ExWAS.R
│   │   └── Incident_T2D_binary_exposure_Logistic_ExWAS.R
│   └── untransformed/
│       ├── Cox_CAD_untransformed_exposure_ExWAS.R
│       ├── Cox_T2D_untransformed_exposure_ExWAS.R
│       ├── Incident_CAD_untransformed_exposure_Logistic_ExWAS.R
│       └── Incident_T2D_untransformed_exposure_Logistic_ExWAS.R
├── Prevalent_analysis/
│   ├── binary/
│   │   ├── Prevalent_CAD_binary_exposure_ExWAS.R
│   │   └── Prevalent_T2D_binary_exposure_ExWAS.R
│   └── untransformed/
│       ├── Prevalent_CAD_untransformed_exposure_ExWAS.R
│       └── Prevalent_T2D_untransformed_exposure_ExWAS.R
└── README.md
```

---

### ▶️ **Running ExWAS Analyses**

#### 1. **Incident Analysis**
- Scripts for binary and untransformed (continuous) variables are located under:
  - `ExWAS/Incident_analysis/binary/`
  - `ExWAS/Incident_analysis/untransformed/`

#### Run Binary Models
```r
# Cox model for CAD with binary exposure
source("ExWAS/Incident_analysis/binary/Cox_CAD_binary_exposure_ExWAS.R")
```

#### Run Untransformed Models
```r
# Cox model for CAD with untransformed exposure
source("ExWAS/Incident_analysis/untransformed/Cox_CAD_untransformed_exposure_ExWAS.R")
```

---

#### 2. **Prevalent Analysis**
- Scripts for binary and untransformed (continuous) variables are located under:
  - `ExWAS/Prevalent_analysis/binary/`
  - `ExWAS/Prevalent_analysis/untransformed/`

#### Run Binary Models
```r
# Prevalent model for CAD with binary exposure
source("ExWAS/Prevalent_analysis/binary/Prevalent_CAD_binary_exposure_ExWAS.R")
```

#### Run Untransformed Models
```r
# Prevalent model for CAD with untransformed exposure
source("ExWAS/Prevalent_analysis/untransformed/Prevalent_CAD_untransformed_exposure_ExWAS.R")
```

---

## 🔬 Mendelian Randomization (MR) Analysis

### 📂 MR Directory Structure
```
MR/
├── compute_mr.R                            # Core MR function used in both analyses
├── CAD_all_exposures_bidirectional_MR_analysis.R   # MR analysis for CAD exposures
├── t2d_all_exposures_bidirectional_MR_analysis.R   # MR analysis for T2D exposures
└── README.md
```

---

### ▶️ **Running MR Analyses**

#### 1. **Load the `compute_mr` Function**
The `compute_mr.R` file is sourced by both MR analysis scripts. Ensure it is loaded before running any MR script.

```r
source("MR/compute_mr.R")
```

#### 2. **Run CAD MR Analysis**
```r
source("MR/CAD_all_exposures_bidirectional_MR_analysis.R")
```

#### 3. **Run T2D MR Analysis**
```r
source("MR/t2d_all_exposures_bidirectional_MR_analysis.R")
```

---

## 📈 Results and Output
- Results will be saved in the `results/` directory.
- Each exposure will generate results files with MR estimates, heterogeneity tests, and sensitivity analyses.

---

## 🛠️ Troubleshooting

- **Missing Packages:** Make sure all required packages are installed.
- **Path Issues:** Set the working directory correctly using `setwd()` before running any script.
- **API Access:** If using external APIs for GWAS data in MR, ensure access is correctly configured.

---

## 📄 License
This project is licensed under the MIT License.

---



![Figure_1](https://github.com/user-attachments/assets/44e0cf60-ac61-4433-bd97-c0689ebd4dce)
