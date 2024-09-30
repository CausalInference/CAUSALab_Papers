################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_setting_restriction.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-24
# This R program is to set the parametric g-formula settings

library(gfoRmula)
library(Hmisc)
# this line is to make sure our customized functions have been loaded.
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_cum.R")

#################################################################
# new additions for the customized time variable
# this maps period values of 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 to
# 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7

coarsen <- function(x){
  return(ifelse(x < 2, 2, floor(x / 2) + 1))
}

fit_custom <- function(covparams, covname, obs_data, j){
  # Note that the following code is never used to estimate the counterfactual mean, etc.
  # It is just here because the package expects some fitted model to be returned
  return(lm(obs_data$period_coarse ~ 1))
}

predict_custom <- function(obs_data, newdf, fitcov, time_name, t, condition, covname, ...){
  coarsen <- function(x){
    return(ifelse(x < 2, 2, floor(x / 2) + 1))
  }
  return(rep(coarsen(t), times = nrow(newdf)))
}
covfits_custom <- c(fit_custom, NA, NA, NA)
covpredict_custom <- c(predict_custom, NA, NA, NA)



id <- "id"
time_name <- "period"
outcome_name <- "event"
outcome_type <- "survival"
covnames <- c(
  "period_coarse",
  "dia", "hbp", "chl",
  "cal", "bmi", "xcond2", 
  "act", "totfv", "totgrn", "junk", 
  "totred", "totproc", "sodajuice",
  "totalc")
covtypes <- c(
  'custom',
  "absorbing", "absorbing", "absorbing",  
  "normal","normal", "absorbing",
  "normal", "normal", "normal", "normal", 
  "normal", "zero-inflated normal", "zero-inflated normal", 
  "zero-inflated normal")



basecovs <- c("age_base", "fhxpca90", "smkhx", "hlthy_base", "fhxmi",
              # pre-baseline variables
              "pre_act", "pre_totfv", "pre_totgrn","pre_junk", "pre_totred",
              "pre_totproc", "pre_sodajuice", "pre_totalc", "pre_cal",
              "base_hbp", "base_chl", "base_dia", "pre_bmi", "pre_hlthy"
)
histories <- c(lagged, cum, tsswitch)
histvars <- list(
  c("dia", "hbp", "chl", "xcond2", 
    "act", "totfv", "totgrn", "junk", 
    "totred", "totproc", "sodajuice",
    "totalc", "cal", "bmi",
    # to use "abosrbing" type, we need to make sure Lag1 term of these variables is included
    "dia", "hbp", "chl", "xcond2"),
  c("dia", "hbp", "chl", "xcond2"),
  c("dia", "hbp", "chl", "xcond2"))

covparams <- list(
  covmodels = c(
    period_coarse ~ 1, 
    # first we model 4 histories diabetes, hbp, chl, and xcond2 based on values one lag before
    # once a history variable has been modeled, itself and related time since switch will be included in the following models
    dia ~
      
      #lag 1 term of time-varying variables
      lag1_dia +  
      lag1_hbp + 
      lag1_chl +  
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0))   +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0, 4.5)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0)) +
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) + as.numeric(lag1_sodajuice >= 2.0) +
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) + as.numeric(lag1_totalc >= 3.0)  +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) +
      as.numeric(lag1_bmi < 2.917) +
      as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +
      as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      # baseline variables
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period),
    
    
    hbp ~   
      # concurrent term of time-varying variables 
      dia +  ts_int_dia +
      
      # lag 1 term
      lag1_dia + 
      lag1_hbp + 
      lag1_chl + 
      lag1_xcond2 +  
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0))  +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0)) +
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +        
      as.numeric(lag1_sodajuice >= 2.0) +
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        
      as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) +
      as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +
      as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+
      fhxmi+
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period),
    
    
    chl ~   
      dia + ts_int_dia +
      hbp + ts_int_hbp +
      
      lag1_dia+ 
      lag1_hbp + 
      lag1_chl + 
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) +
      as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +
      as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period),
    
    
    cal~
      # concurrent terms
      dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) +
      as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +
      as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period_coarse),
    
    bmi~  
      # concurrent terms
      dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) + cal:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12))+
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 +  
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) +
      as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +
      as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period),
    
    
    
    xcond2 ~   
      dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) + cal:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12))+
      as.numeric(bmi < 2.917) + as.numeric(bmi >= 2.917 & bmi <3.219) +as.numeric(bmi >= 3.219 & bmi < 3.401) +
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 +  
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) + as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx+
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90 + 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period),
    
    
    # make models for the intervention variables 
    
    act ~ dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      xcond2 + ts_int_xcond2 +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) + cal:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12))+
      as.numeric(bmi < 2.917) + as.numeric(bmi >= 2.917 & bmi <3.219) +as.numeric(bmi >= 3.219 & bmi < 3.401) +
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) + as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period),
    
    totfv ~  dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      xcond2 + ts_int_xcond2 +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046)+
      as.numeric(bmi < 2.917) + as.numeric(bmi >= 2.917 & bmi <3.219) +as.numeric(bmi >= 3.219 & bmi < 3.401) +
      Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) + as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period_coarse) ,
    
    totgrn ~  dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      xcond2 + ts_int_xcond2 +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) + 
      as.numeric(bmi < 2.917) + as.numeric(bmi >= 2.917 & bmi <3.219) +as.numeric(bmi >= 3.219 & bmi < 3.401) +
      Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(totfv, knots = c(2.4, 5.0, 9.0)) + 
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) + as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period_coarse),
    
    junk ~ dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      xcond2 + ts_int_xcond2 +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) + 
      as.numeric(bmi < 2.917) + as.numeric(bmi >= 2.917 & bmi <3.219) +as.numeric(bmi >= 3.219 & bmi < 3.401) +
      Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(totfv, knots = c(2.4, 5.0, 9.0)) +  
      Hmisc::rcspline.eval(totgrn, knots = c(0.5, 1.5, 3.0)) + 
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) + as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period_coarse),
    
    totred ~   dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      xcond2 + ts_int_xcond2 +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) + 
      as.numeric(bmi < 2.917) + as.numeric(bmi >= 2.917 & bmi <3.219) +as.numeric(bmi >= 3.219 & bmi < 3.401) +
      Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(totfv, knots = c(2.4, 5.0, 9.0)) +  
      Hmisc::rcspline.eval(totgrn, knots = c(0.5, 1.5, 3.0)) +  
      Hmisc::rcspline.eval(junk, knots = c(0.1, 1.0, 2.5, 5.0)) + 
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      as.numeric(lag1_bmi < 2.917) + as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      fhxpca90+ 
      fhxmi+
      as.factor(hlthy_base) +
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period_coarse),
    
    totproc ~ dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      xcond2 + ts_int_xcond2 +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) + 
      as.numeric(bmi < 2.917) + as.numeric(bmi >= 2.917 & bmi <3.219) +as.numeric(bmi >= 3.219 & bmi < 3.401) +
      Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(totfv, knots = c(2.4, 5.0, 9.0)) +  
      Hmisc::rcspline.eval(totgrn, knots = c(0.5, 1.5, 3.0)) + 
      Hmisc::rcspline.eval(junk, knots = c(0.1, 1.0, 2.5, 5.0)) + 
      Hmisc::rcspline.eval(totred, knots = c(1.0, 3.0, 7.0))+ 
      
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period_coarse),
    
    sodajuice ~   dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      xcond2 + ts_int_xcond2 +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) + 
      as.numeric(bmi < 2.917) + as.numeric(bmi >= 2.917 & bmi <3.219) +as.numeric(bmi >= 3.219 & bmi < 3.401) +
      Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(totfv, knots = c(2.4, 5.0, 9.0)) +  
      Hmisc::rcspline.eval(totgrn, knots = c(0.5, 1.5, 3.0)) +  
      Hmisc::rcspline.eval(junk, knots = c(0.1, 1.0, 2.5, 5.0)) + 
      Hmisc::rcspline.eval(totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(totproc, knots = c(0.1, 1.0, 5.0)) +
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period_coarse),
    
    totalc~ dia + ts_int_dia + 
      hbp + ts_int_hbp + 
      chl + ts_int_chl +
      xcond2 + ts_int_xcond2 +
      as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) + 
      as.numeric(bmi < 2.917) + as.numeric(bmi >= 2.917 & bmi <3.219) +as.numeric(bmi >= 3.219 & bmi < 3.401) +
      Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(totfv, knots = c(2.4, 5.0, 9.0)) +  
      Hmisc::rcspline.eval(totgrn, knots = c(0.5, 1.5, 3.0)) + 
      Hmisc::rcspline.eval(junk, knots = c(0.1, 1.0, 2.5, 5.0)) + 
      Hmisc::rcspline.eval(totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(sodajuice >= 0.1 & sodajuice < 1.0) + as.numeric(sodajuice >= 1.0 & sodajuice < 2.0) +   as.numeric(sodajuice >= 2.0)+
      
      
      # lag 1 terms of time-varying variables
      lag1_dia + 
      lag1_hbp + 
      lag1_chl +
      lag1_xcond2 + 
      Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
      Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
      Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
      Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
      Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
      Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
      as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
      as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
      as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
      
      smkhx + 
      as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
      as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
      as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
      as.numeric(age_base >=70 ) + 
      as.factor(hlthy_base) +
      fhxpca90+ 
      fhxmi+ 
      # pre-baseline variables
      as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
      as.numeric(pre_act >= 50) + 
      as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
      as.numeric(pre_totfv >= 9.0)+
      as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
      as.numeric(pre_totgrn >= 3.0) + 
      as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
      as.numeric(pre_junk >= 5.0) + 
      as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
      as.numeric(pre_totred >= 7.0) + 
      as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
      as.numeric(pre_totproc >= 5.0) +
      as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
      as.numeric(pre_sodajuice >= 2.0) +
      as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
      as.numeric(pre_totalc >= 3.0) +
      as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
      as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
      as.numeric(pre_bmi>=30) +
      base_dia + base_chl + base_hbp + 
      as.factor(period_coarse)))


compevent_name <- "dead_other"
compevent_model <- dead_other ~ 
  Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0)) + 
  Hmisc::rcspline.eval(totfv, knots = c(2.4, 5.0, 9.0)) + 
  Hmisc::rcspline.eval(totgrn, knots = c(0.5, 1.5, 3.0))+
  Hmisc::rcspline.eval(junk, knots = c(0.1, 1.0, 2.5, 5.0))+ 
  Hmisc::rcspline.eval(totred, knots = c(1.0, 3.0, 7.0))+ 
  Hmisc::rcspline.eval(totproc, knots = c(0.1, 1.0, 5.0))+ 
  as.numeric(sodajuice >= 0.1 & sodajuice < 1.0) + as.numeric(sodajuice >= 1.0 & sodajuice < 2.0) +               as.numeric(sodajuice >= 2.0)+ 
  as.numeric(totalc >= 0.1 & totalc < 0.5) + as.numeric(totalc >= 0.5 & totalc < 3.0) +        as.numeric(totalc >= 3.0) + 
  as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) +
  
  totfv:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) +
  totgrn:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  junk:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) +
  totred:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  totproc:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  sodajuice:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  totalc:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  cal:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12))+ 
  
  as.numeric(bmi < 2.917) +
  as.numeric(bmi >= 2.917 & bmi <3.219) +
  as.numeric(bmi >= 3.219 & bmi < 3.401) +
  dia + ts_int_dia+
  hbp + ts_int_hbp+
  chl + ts_int_chl+
  xcond2 + ts_int_xcond2+
  
  
  lag1_dia + 
  lag1_hbp + 
  lag1_chl +
  lag1_xcond2 + 
  Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
  Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
  Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
  Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
  Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
  Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
  as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
  as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
  as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
  as.numeric(lag1_bmi < 2.917) + as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
  
  
  smkhx +
  as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
  as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
  as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
  as.numeric(age_base >=70 ) + 
  as.factor(hlthy_base) +
  fhxpca90 +
  fhxmi+
  # pre-baseline variables
  as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
  as.numeric(pre_act >= 50) + 
  as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
  as.numeric(pre_totfv >= 9.0)+
  as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
  as.numeric(pre_totgrn >= 3.0) + 
  as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
  as.numeric(pre_junk >= 5.0) + 
  as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
  as.numeric(pre_totred >= 7.0) + 
  as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
  as.numeric(pre_totproc >= 5.0) +
  as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
  as.numeric(pre_sodajuice >= 2.0) +
  as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
  as.numeric(pre_totalc >= 3.0) +
  as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
  as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
  as.numeric(pre_bmi>=30) +
  base_dia + base_chl + base_hbp + 
  as.factor(period)

# add censoring models
censor_name <- "censor"
censor_model <- censor ~ 
  Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0)) + 
  Hmisc::rcspline.eval(totfv, knots = c(2.4, 5.0, 9.0)) + 
  Hmisc::rcspline.eval(totgrn, knots = c(0.5, 1.5, 3.0, 4.5))+
  Hmisc::rcspline.eval(junk, knots = c(1.0, 2.5, 5.0))+ 
  Hmisc::rcspline.eval(totred, knots = c(1.0, 3.0, 7.0))+ 
  Hmisc::rcspline.eval(totproc, knots = c(0.1, 1.0, 5.0))+ 
  as.numeric(sodajuice >= 0.1 & sodajuice < 1.0) + as.numeric(sodajuice >= 1.0 & sodajuice < 2.0) +
  as.numeric(sodajuice >= 2.0)+ 
  as.numeric(totalc >= 0.1 & totalc < 0.5) + as.numeric(totalc >= 0.5 & totalc < 3.0) +       
  as.numeric(totalc >= 3.0)  +
  as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) +
  
  totfv:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) +
  totgrn:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  junk:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) +
  totred:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  totproc:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  sodajuice:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  totalc:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  cal:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12))+ 
  
  as.factor(bmi < 2.917) +
  as.factor(bmi >= 2.917 & bmi <3.219) +
  as.factor(bmi >= 3.219 & bmi < 3.401) +
  dia + ts_int_dia+
  hbp + ts_int_hbp+
  chl + ts_int_chl+
  xcond2 + ts_int_xcond2+
  
  
  lag1_dia + 
  lag1_hbp + 
  lag1_chl +
  lag1_xcond2 + 
  Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
  Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
  Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
  Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
  Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
  Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
  as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
  as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
  as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
  as.numeric(lag1_bmi < 2.917) + as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
  smkhx +
  as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
  as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
  as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
  as.numeric(age_base >=70 ) + 
  as.factor(hlthy_base) + 
  fhxpca90 +
  fhxmi+ 
  as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
  as.numeric(pre_act >= 50) + 
  as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
  as.numeric(pre_totfv >= 9.0)+
  as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
  as.numeric(pre_totgrn >= 3.0) + 
  as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
  as.numeric(pre_junk >= 5.0) + 
  as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
  as.numeric(pre_totred >= 7.0) + 
  as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
  as.numeric(pre_totproc >= 5.0) +
  as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
  as.numeric(pre_sodajuice >= 2.0) +
  as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
  as.numeric(pre_totalc >= 3.0) +
  as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
  as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
  as.numeric(pre_bmi>=30) +
  as.numeric(base_dia) + as.numeric(base_chl) + as.numeric(base_hbp) + 
  as.factor(period)


ymodel <- event ~ 
  Hmisc::rcspline.eval(act, knots = c(0.1, 7.5, 20.0, 50.0)) + 
  Hmisc::rcspline.eval(totfv, knots = c(2.4, 5.0, 9.0)) + 
  Hmisc::rcspline.eval(totgrn, knots = c(0.5, 1.5, 3.0))+
  Hmisc::rcspline.eval(junk, knots = c(0.1, 1.0, 2.5, 5.0))+ 
  Hmisc::rcspline.eval(totred, knots = c(1.0, 3.0, 7.0))+ 
  Hmisc::rcspline.eval(totproc, knots = c(0.1, 1.0, 5.0))+ 
  as.numeric(sodajuice >= 0.1 & sodajuice < 1.0) + as.numeric(sodajuice >= 1.0 & sodajuice < 2.0) +               as.numeric(sodajuice >= 2.0)+ 
  as.numeric(totalc >= 0.1 & totalc < 0.5) + as.numeric(totalc >= 0.5 & totalc < 3.0) +        as.numeric(totalc >= 3.0) + 
  as.numeric(cal>= 7.495542 & cal < 7.824046) + as.numeric(cal>= 7.824046) +
  
  totfv:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) +
  totgrn:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  junk:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) +
  totred:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  totproc:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  sodajuice:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  totalc:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12)) + 
  cal:factor(period %in% c(-1, 0, 2, 4, 6, 8, 10, 12))+ 
  
  as.numeric(bmi < 2.917) +
  as.numeric(bmi >= 2.917 & bmi <3.219) +
  as.numeric(bmi >= 3.219 & bmi < 3.401) +
  dia + ts_int_dia+
  hbp + ts_int_hbp+
  chl + ts_int_chl+
  xcond2 + ts_int_xcond2+
  
  lag1_dia + 
  lag1_hbp + 
  lag1_chl +
  lag1_xcond2 + 
  Hmisc::rcspline.eval(lag1_act, knots = c(0.1, 7.5, 20.0, 50.0))  +
  Hmisc::rcspline.eval(lag1_totfv, knots = c(2.4, 5.0, 9.0)) +
  Hmisc::rcspline.eval(lag1_totgrn, knots = c(0.5, 1.5, 3.0)) +
  Hmisc::rcspline.eval(lag1_junk, knots = c(0.1, 1.0, 2.5, 5.0)) +
  Hmisc::rcspline.eval(lag1_totred, knots = c(1.0, 3.0, 7.0))+
  Hmisc::rcspline.eval(lag1_totproc, knots = c(0.1, 1.0, 5.0)) +
  as.numeric(lag1_sodajuice >= 0.1 & lag1_sodajuice < 1.0) + as.numeric(lag1_sodajuice >= 1.0 & lag1_sodajuice < 2.0) +               as.numeric(lag1_sodajuice >= 2.0)+
  as.numeric(lag1_totalc >= 0.1 & lag1_totalc < 0.5) + as.numeric(lag1_totalc >= 0.5 & lag1_totalc < 3.0) +        as.numeric(lag1_totalc >= 3.0) +
  as.numeric(lag1_cal>= 7.495542 & lag1_cal < 7.824046) + as.numeric(lag1_cal>= 7.824046) + 
  as.numeric(lag1_bmi < 2.917) + as.numeric(lag1_bmi >= 2.917 & lag1_bmi <3.219) +as.numeric(lag1_bmi >= 3.219 & lag1_bmi < 3.401) +
  smkhx +
  as.numeric(age_base < 45) + as.numeric(age_base >= 45 & age_base < 50) +
  as.numeric(age_base >= 50 & age_base < 55) + as.numeric(age_base >= 55 & age_base < 60) +
  as.numeric(age_base >= 60 & age_base < 65) + as.numeric(age_base >= 65& age_base < 70) + 
  as.numeric(age_base >=70 ) + 
  as.factor(hlthy_base) +
  fhxpca90 +
  fhxmi+
  # pre-baseline variables
  as.numeric(pre_act >= 7.5 & pre_act < 20) + as.numeric(pre_act >= 20 & pre_act < 50) + 
  as.numeric(pre_act >= 50) + 
  as.numeric(pre_totfv >= 2.4 & pre_totfv < 5.0) + as.numeric(pre_totfv >= 5.0 & pre_totfv < 9.0) +
  as.numeric(pre_totfv >= 9.0)+
  as.numeric(pre_totgrn >= 0.5 & pre_totgrn <1.5) + as.numeric(pre_totgrn >=1.5 & pre_totgrn < 3.0) +
  as.numeric(pre_totgrn >= 3.0) + 
  as.numeric(pre_junk >= 1.0 & pre_junk < 2.5) + as.numeric(pre_junk >= 2.5 & pre_junk < 5.0) + 
  as.numeric(pre_junk >= 5.0) + 
  as.numeric(pre_totred >= 1.0 & pre_totred < 3.0) + as.numeric(pre_totred >= 3.0 & pre_totred < 7.0) +
  as.numeric(pre_totred >= 7.0) + 
  as.numeric(pre_totproc >= 0.1 & pre_totproc < 1.0) + as.numeric(pre_totproc >= 1.0 & pre_totproc < 5.0) +
  as.numeric(pre_totproc >= 5.0) +
  as.numeric(pre_sodajuice >= 0.1 & pre_sodajuice < 1.0) + as.numeric(pre_sodajuice >= 1.0 & pre_sodajuice < 2.0) +
  as.numeric(pre_sodajuice >= 2.0) +
  as.numeric(pre_totalc >= 0.1 & pre_totalc < 1.0) + as.numeric(pre_totalc >= 1.0 & pre_totalc < 2.0) +
  as.numeric(pre_totalc >= 3.0) +
  as.numeric(pre_cal >= 1800 & pre_cal < 2500) + as.numeric(pre_cal >= 2500) +
  as.numeric(pre_bmi>= 18.5 & pre_bmi < 25) + as.numeric(pre_bmi >= 25 & pre_bmi <30) + 
  as.numeric(pre_bmi>=30) +
  base_dia + base_chl + base_hbp + 
  as.factor(period)
###################
# made restrictions to diets because we know that these values were only measured in
# 1986, 1990, 1994, 1998, 2002, 2006, 2010, 2014
# period 1(1992), 3(1996), 5(2000), 7(2004), 9(2008), 11(2012), 13(2016) are carried from previous survey
# at period 0,2,4,6,8,10 the intake should be modeled
restrictions <- list(c('totfv', 'period %in% c(0, 2, 4, 6, 8, 10, 12)', carry_forward),
                     c('totgrn', 'period %in% c(0, 2, 4, 6, 8, 10, 12)', carry_forward), 
                     c('junk', 'period %in% c(0, 2, 4, 6, 8, 10, 12)', carry_forward),
                     c('totred', 'period %in% c(0, 2, 4, 6, 8, 10, 12)', carry_forward),
                     c('totproc', 'period %in% c(0, 2, 4, 6, 8, 10, 12)', carry_forward),
                     c('sodajuice', 'period %in% c(0, 2, 4, 6, 8, 10, 12)', carry_forward),
                     c('totalc', 'period %in% c(0, 2, 4, 6, 8, 10, 12)', carry_forward),
                     c('cal', 'period %in% c(0, 2, 4, 6, 8, 10, 12)', carry_forward))

###################################################
# This line is a new intervention function

xcond_int <- function(newdf, pool, intvar, intvals, time_name, t) {
  if (nrow(newdf[newdf[[intvar]] < intvals[[1]] & newdf[["xcond2"]] == 0]) != 0) {
    classtmp <- class(newdf[[intvar]])
    myclass <- paste("as.", classtmp, sep = "")
    newdf[newdf[[intvar]] < intvals[[1]] & newdf[["xcond2"]] == 0, `:=`((intvar), 
                                                                        get(myclass)(intvals[[1]]))]
  }
  if (nrow(newdf[newdf[[intvar]] > intvals[[2]] & newdf[["xcond2"]] == 0]) != 0) {
    classtmp <- class(newdf[[intvar]])
    myclass <- paste("as.", classtmp, sep = "")
    newdf[newdf[[intvar]] > intvals[[2]] & newdf[["xcond2"]] == 0, `:=`((intvar), 
                                                                        get(myclass)(intvals[[2]]))]
  }
} 

###########################
