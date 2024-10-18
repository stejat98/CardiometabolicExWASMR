## Sivateja Tangirala
## 10/18/2024
## Baseline Functions (coxph regression)

library(tidyverse)
library(gdata)
library(RNOmni)
library(biglm)
library(broom)
library(pROC)
library(survival)

constructFormulaCoxph <- function (exposure, data, depvar, adjustments) {
  cat(sprintf("constructing formula for %s regressed on %s and adjusted for %s \n",
              depvar, exposure, paste(adjustments, collapse = ' + ')))
  formula <- NULL
  if (!is.null(adjustments)) {
    formula <-   {`if` (is.numeric(data[,exposure]),
                        sprintf("surv_object ~ I(scale(%s)) + %s", exposure, adjustments = paste(adjustments, collapse = ' + ')),
                        sprintf("surv_object ~ %s + %s", exposure, adjustments = paste(adjustments, collapse = ' + ')))}
  } else {
    formula <-   {`if` (is.numeric(data[,exposure]),
                        sprintf("surv_object ~ I(scale(%s))", exposure),
                        sprintf("surv_object ~ %s", exposure))}
  }
  cat(sprintf("Formula: %s \n", formula))
  return (list(formula,exposure,adjustments))
}


executeModelCoxph <- function (formula,data,exposure,depvar,adjustments) {
  mod <- NULL
  results <- NULL
  depvar_time <- sprintf("%s_time",depvar)
  depvar_censored <- sprintf("%s_censored",depvar)
  df_full <- data[,c(exposure,adjustments,depvar_time, depvar_censored)]
  
  df_filter_NA_full <- df_full %>% drop_na()
  
  df_filter_NA_2_full <- df_filter_NA_full[, !duplicated(colnames(df_filter_NA_full), fromLast = TRUE)]
  
  # df_filter_NA_3_full <- drop.levels(df_filter_NA_2_full)
  
  data1 <- df_filter_NA_2_full
  
  surv_object <- Surv(time = data1[,depvar_time], event = data1[,depvar_censored])
  
  cat(sprintf("Executing formula: %s \n", formula))
  cat(exposure)
  print(sprintf("\n The dimension of data1: %i",dim(data1)))
  tryCatch(
    mod1 <- coxph(as.formula(formula), data = data1),
    error = function(e) {
      cat(sprintf("Failed on: %s \n", formula))
      rm(formula)
    }) 
  tryCatch(
    mod2 <- coxph(as.formula(sprintf("surv_object ~ %s", adjustments = paste(adjustments, collapse = ' + '))),data=data1) ,
    error = function(e) {
      cat(sprintf("mod2 Failed on: %s \n", formula))
      rm(formula)
    })
  if(exists("mod1") & exists("mod2")){
    if (!is.null(mod1) & !is.null(mod2)){
      tryCatch({
        results1 <- mod1 %>% tidy()
      },
      error = function(e) {
        cat(sprintf("Failed on: %s  \n", formula))
        
      })
      if (!is.null(results1)) {
        tryCatch({results2 <- mod2 %>% tidy()
        
        results_bind_1 <- cbind(results1,Cstat = rep(summary(mod1)$concordance[1],times=nrow(results1)), SE_Cstat = rep(summary(mod1)$concordance[2],times=nrow(results1)),Cstat_adjVariables = rep(summary(mod2)$concordance[1],times=nrow(results1)),SE_Cstat_adjVariables = rep(summary(mod2)$concordance[2],times=nrow(results1)),Phenotype= rep(depvar,times=nrow(results1)),SampleSize = rep(nrow(data1),times=nrow(results1)),Exposure=rep(exposure,times=nrow(results1)))
        rm(results1,results2)},error = function(e){
          cat(sprintf("Failed on: coxph --  %s or %s  \n", sprintf("%s (full model)", depvar),sprintf("%s (baseline covariates only model)", depvar)))
          results_bind_1 <- cbind(results1,Cstat = rep(NA,times=nrow(results1)), SE_Cstat = rep(NA,times=nrow(results1)),Cstat_adjVariables = rep(NA,times=nrow(results1)),SE_Cstat_adjVariables = rep(NA,times=nrow(results1)),Phenotype= rep(depvar,times=nrow(results1)),SampleSize = rep(nrow(data1),times=nrow(results1)),Exposure=rep(exposure,times=nrow(results1)))
          
        })
      }
      
    }
    if(exists("results_bind_1")){
      return(as.data.frame(results_bind_1))
    }
  }
  
}

executeEWASCoxph.map <- function  (data, depvar,adjustments,
                                      exposures) {
  
  
  return_df  <- exposures %>%
    map(constructFormulaCoxph, data=data, depvar, adjustments) %>%
    map_dfr(.f= function(x){executeModelCoxph(x[[1]],data,exposure =x[[2]],depvar,adjustments=x[[3]])})
  
  return(return_df)
}



EWASCoxph <- function  (data, depvar,adjustments,exposures,
                           outFileName) {
  
  cat(sprintf("executing EWAS (Coxph) for %s", depvar))
  
  cols <- ncol(data)
  
  executeEWASCoxph (data, depvar, adjustments,
                       exposures, outFileName)
  
  
  
}


executeEWASCoxph <- function (data, depvar, adjustments,
                                 exposures, outFileName) {
  ## call executeEWASCoxph.map and save results
  
  results <- executeEWASCoxph.map(data,depvar=depvar,
                                     adjustments=adjustments,
                                     exposures=exposures)
  
  saveRDS(results,sprintf("%s.RDS",outFileName))
  
  
}





