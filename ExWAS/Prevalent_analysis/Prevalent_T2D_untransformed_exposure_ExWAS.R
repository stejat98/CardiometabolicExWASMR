library(tidyverse)
library(gdata)
library(RNOmni)
library(biglm)
library(broom)
library(fst)


source("Baseline_PEWAS_Logistic_Functions_script.R") 



load("adjustments_survival_analysis.Rdata")



data <- read.fst("prevalent_t2d_ExWAS_input_data.fst")
data <- data %>% filter(f.21000.0.0 %in% c("British","Irish", "Any other white background"))

adjustments <- c(adjustments, "x.738")

load("continuous_ordinal_exposures.Rdata")

exposures_all_mappings <- read_csv("exposures_id_name_mapping_updated_10_06_2022.csv")

nutrients_mappings_filtered <- exposures_all_mappings %>% filter(Category.Specific %in% c("Estimated food nutrients yesterday","Total weight by food group yesterday"))


new_nutrient_exposures <- nutrients_mappings_filtered$variableName

new_nutrient_exposures <- sprintf("x.%s",new_nutrient_exposures)

exposures <- c(contin_ord_vars, new_nutrient_exposures)


## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data, depvar = "prevalent_t2d"
             ,adjustments = adjustments,exposures=exposures,outFileName = sprintf("Prevalent_T2D_Untransformed_Continuous_Exposures_Baseline_Actual_06_28_23_Baseline_Pheno%sOverall_biglm_results", "T2D"))
