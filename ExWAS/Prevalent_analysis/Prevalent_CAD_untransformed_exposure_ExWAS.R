library(tidyverse)
library(gdata)
library(RNOmni)
library(biglm)
library(broom)
library(fst)


source("/home/st320/Baseline_PEWAS_Logistic_Functions_script.R") 



load("/home/st320/UKB_PEWAS/adjustments_survival_analysis.Rdata")



data <- read.fst("/n/groups/patel/sivateja/UKB/prevalent_disease_ExWAS_input_data_08_31_2024.fst")
data <- data %>% filter(f.21000.0.0 %in% c("British","Irish", "Any other white background"))

adjustments <- c(adjustments, "x.738")

load("/home/st320/UKB_PEWAS/continuous_ordinal_exposures.Rdata")

exposures_all_mappings <- read_csv("exposures_id_name_mapping_updated_10_06_2022.csv")

nutrients_mappings_filtered <- exposures_all_mappings %>% filter(Category.Specific %in% c("Estimated food nutrients yesterday","Total weight by food group yesterday"))


new_nutrient_exposures <- nutrients_mappings_filtered$variableName

new_nutrient_exposures <- sprintf("x.%s",new_nutrient_exposures)

exposures <- c(contin_ord_vars, new_nutrient_exposures)


## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data, depvar = "prevalent_cad"
             ,adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/Prevalent_CAD_Untransformed_Continuous_Exposures_Baseline_Actual_06_28_23_Baseline_Pheno%sOverall_biglm_results", "CAD"))
