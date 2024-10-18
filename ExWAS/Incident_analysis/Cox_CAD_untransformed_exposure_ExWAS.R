library(tidyverse)
library(gdata)
library(RNOmni)
library(biglm)
library(broom)
library(fst)
library(survival)



source("/home/st320/Baseline_PEWAS_Cox_Functions_script.R") 



load("/home/st320/UKB_PEWAS/adjustments_survival_analysis.Rdata")


data <- read.fst("/n/groups/patel/sivateja/UKB/all_cause_mort_data_exposures_full.fst")
data <- data %>% filter(f.21000.0.0 %in% c("British","Irish", "Any other white background"))

adjustments <- c(adjustments, "x.738")

load("/home/st320/UKB_PEWAS/continuous_ordinal_exposures.Rdata")

exposures <- contin_ord_vars


surv_object <- Surv(time = data$cad_time, event = data$cad_censored)

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASCoxph(data=data, depvar = "cad"
             ,adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/Cox_Untransformed_Continuous_Exposures_Baseline_Actual_06_28_23_Baseline_Pheno%sOverall_biglm_results", "CAD"))
