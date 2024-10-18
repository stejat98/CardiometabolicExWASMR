library(tidyverse)
library(gdata)
library(RNOmni)
library(biglm)
library(broom)
library(fst)


source("Baseline_PEWAS_Logistic_Functions_script.R") 



load("adjustments_survival_analysis.Rdata")


data <- read.fst("incident_ExWAS_input_data.fst")
data <- data %>% filter(f.21000.0.0 %in% c("British","Irish", "Any other white background"))

adjustments <- c(adjustments, "x.738")

load("continuous_ordinal_exposures.Rdata")

exposures <- contin_ord_vars



## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data, depvar = "cad"
             ,adjustments = adjustments,exposures=exposures,outFileName = sprintf("Untransformed_Continuous_Exposures_Baseline_Actual_06_28_23_Baseline_Pheno%sOverall_biglm_results", "CAD"))
