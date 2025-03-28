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
## ▶️ Running Scripts from Command Line

### 1. **Run ExWAS Scripts from Command Line**
- Example: Run Cox model for CAD with binary exposure:
```bash
Rscript ExWAS/Incident_analysis/binary/Cox_CAD_binary_exposure_ExWAS.R
```

- Example: Run prevalent model for T2D with untransformed exposure:
```bash
Rscript ExWAS/Prevalent_analysis/untransformed/Prevalent_T2D_untransformed_exposure_ExWAS.R
```

---

### 2. **Run MR Scripts from Command Line**
- Example: Run CAD MR analysis:
```bash
Rscript MR/CAD_all_exposures_bidirectional_MR_analysis.R
```

- Example: Run T2D MR analysis:
```bash
Rscript MR/t2d_all_exposures_bidirectional_MR_analysis.R
```


---

## 📈 Results and Output
- Results will be saved in the `results/` directory.
- Each exposure will generate results files with MR estimates, heterogeneity tests, and sensitivity analyses.

---

## 📄 License
This project is licensed under the MIT License.

---



![Figure_1](https://github.com/user-attachments/assets/44e0cf60-ac61-4433-bd97-c0689ebd4dce)
