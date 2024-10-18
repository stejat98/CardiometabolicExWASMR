# Install necessary packages if not already installed
required_packages <- c("tidyverse", "TwoSampleMR", "data.table", "purrr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)

library(TwoSampleMR)

library(data.table)

library(purrr)


options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

# setwd to desired path that contains input files and where the files will be outputted to
#setwd("")

exposures_all_mappings <- read_csv("exposures_id_name_mapping_updated_10_06_2022.csv")

mr_id_vector <- na.omit(unique(exposures_all_mappings$mr_id))


compute_mr <- function(exposure,pheno){
  tryCatch(
    {
      exposure_dat <- extract_instruments(outcomes = c(exposure))
      exposure_dat <- clump_data(exposure_dat)
      outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(pheno), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
      dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
      mr_results <- mr(dat)
      mr_results <- generate_odds_ratios(mr_results)
      mr_results_temp <- mr_results
      mr_results_temp$mrID <- exposure
      mr_results_temp$direction <- "forward"
      
      ## testing reversal
      exposure_dat2 <- extract_instruments(outcomes = c(pheno))
      exposure_dat2 <- clump_data(exposure_dat2)
      outcome_dat2 <- extract_outcome_data(exposure_dat2$SNP, c(exposure), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
      dat2 <- harmonise_data(exposure_dat2, outcome_dat2, action = 2)
      mr_results_2 <- mr(dat2)
      mr_results_2_temp <- mr_results_2
      mr_results_2_temp$mrID <- exposure
      mr_results_2_temp$direction <- "reverse"
      colnames(mr_results_2_temp) <- map_chr(colnames(mr_results_2_temp),function(x){sprintf("reverse_%s",x)})
      mr_results_temp_full <- cbind(mr_results_temp,mr_results_2_temp)
      
      return(mr_results_temp_full)
      
    },
    error = function(e) {
      cat(sprintf("MR Failed on: exposure: %s ,  pheno: %s \n", exposure, pheno))
    })
  
}




mr_results_filtered <- map_dfr(mr_id_vector,function(x){compute_mr(x,pheno='finn-a-I9_CHD')})

saveRDS(mr_results_filtered,"all_exposures_MR_CAD_results.RDS")

all_exposures_MR_CAD_results_IVW <- mr_results_filtered %>% filter(method == "Inverse variance weighted" & reverse_method == "Inverse variance weighted")

write_csv(all_exposures_MR_CAD_results_IVW, "all_exposures_MR_CAD_results_IVW.csv")

