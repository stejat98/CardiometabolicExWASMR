library(tidyverse)
library(gdata)
library(RNOmni)
library(biglm)
library(broom)
library(fst)

source("Baseline_PEWAS_Logistic_Functions_script.R") 


## load list of adjustment variables
load("adjustments_survival_analysis.Rdata")


## load input data file containing exposures and prevalent disease outcome data
data <- read.fst("prevalent_disease_ExWAS_input_data_08_31_2024.fst")
data <- data %>% filter(f.21000.0.0 %in% c("British","Irish", "Any other white background"))


load("binary_exposures.Rdata")

exposures <- binvals[-length(binvals)]

for (exposureCol in exposures){
  data[[exposureCol]] <- as.factor(data[[exposureCol]])
}


## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data, depvar = "prevalent_cad"
             ,adjustments = adjustments,exposures=exposures,outFileName = sprintf("Prevalent_CAD_Binary_Exposures_Baseline_Actual_03_29_22_Baseline_Pheno%sOverall_biglm_results", "CAD"))
