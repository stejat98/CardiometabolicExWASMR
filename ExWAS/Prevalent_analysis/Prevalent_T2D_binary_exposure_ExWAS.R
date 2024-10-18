library(tidyverse)
library(gdata)
library(RNOmni)
library(biglm)
library(broom)
library(fst)

source("/home/st320/Baseline_PEWAS_Logistic_Functions_script.R") 



load("/home/st320/UKB_PEWAS/adjustments_survival_analysis.Rdata")



data <- read.fst("/n/groups/patel/sivateja/UKB/prevalent_t2d_ExWAS_input_data_08_31_2024.fst")

data <- data %>% filter(f.21000.0.0 %in% c("British","Irish", "Any other white background"))


load("/home/st320/UKB_PEWAS/binary_exposures.Rdata")

exposures <- binvals[-length(binvals)]

for (exposureCol in exposures){
  data[[exposureCol]] <- as.factor(data[[exposureCol]])
}


## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data, depvar = "prevalent_t2d"
             ,adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/Prevalent_T2D_Binary_Exposures_Baseline_Actual_03_29_22_Baseline_Pheno%sOverall_biglm_results", "T2D"))
