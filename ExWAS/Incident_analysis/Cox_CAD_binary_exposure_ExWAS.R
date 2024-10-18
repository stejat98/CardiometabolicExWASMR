library(tidyverse)
library(gdata)
library(RNOmni)
library(biglm)
library(broom)
library(fst)
library(survival)



source("Baseline_PEWAS_Cox_Functions_script.R") 


load("adjustments_survival_analysis.Rdata")


data <- read.fst("incident_ExWAS_input_data.fst")
data <- data %>% filter(f.21000.0.0 %in% c("British","Irish", "Any other white background"))

adjustments <- c(adjustments, "x.738")

load("binary_exposures.Rdata")

exposures <- binvals[-length(binvals)]

for (exposureCol in exposures){
  data[[exposureCol]] <- as.factor(data[[exposureCol]])
}

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASCoxph(data=data, depvar = "cad"
          ,adjustments = adjustments,exposures=exposures,outFileName = sprintf("Coxph_Binary_Exposures_Baseline_Actual_06_28_23_Baseline_Pheno%sOverall_biglm_results", "CAD"))
