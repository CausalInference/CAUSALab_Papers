################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/sens_lag1_term/gformula_lag1_term_tot_adv.R
# Programer: Fuyu Guo
# Last updated: 2024-Mar-24
# Please note: this R code used the same model settings as the one for the primary analysis
# we only changed the input dataset whereby treatments and covariates were lagged by 1 term 
# except for the composite variable for cancer diagnoses 
# (except prostate cancer and nonmelanoma skin cancer),
# myocardial infarction, stroke, congestive heart failure, and functional impairment. 

# 1) Purpose: this R program used the parametric g-formula to conduct the sensitivity analysis
# for minimum latency period of 1 term
# 2) Outcome: total prostate cancer; advanced prostate cancer
# 3) input files:
#     - SAS dataset: 
#         "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/lag_mlp/wcrf_fg_1115_tot.sas7bdat"
#         "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/lag_mlp/wcrf_fg_1115_adv.sas7bdat"
#     - R codes:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_analyses_clean.R"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_cum.R"
#        ""/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_setting_restriction.R""
# 4) Study design: Prospective cohort
#
# 5) Follow-up: HPFS: 1990-2016 (baseline : 1990 ; pre-baseline: 1986)


###############################################
# R settings
options(warn = -1)

# Make sure the library path is correct and load pacakges
library(tidyr)
library(readr)
library(dplyr)
library(haven)
library(data.table)
library(gfoRmulaMLP)
library(parallel)
setwd("/n/hpnh/Users/fyguo/R_proj/code_R3/sens_lag1_term")


# Use supporting R codes
# This R program is used to load functions "extract_prebase_base()" and "trans_var()" which
# will be used to clean data and makes it ready for g-formula
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_analyses_clean.R")

# This R program aims to create a customized function to be used by g-formula
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_cum.R")


fit <- list(2L)
######################################
# total prostate cancer
# read in newly generated total data
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/lag_mlp/wcrf_fg_1115_tot.sas7bdat")

#  extract prebaseline, baseline variables, and
# conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)

source("/n/hpnh/Users/fyguo/R_proj/code_R3/sens_lag1_term/gformula_setting_restriction.R")
ncores <- 20

baseline_poprisk =0.01142210
baseline_prodp0 = 1- baseline_poprisk
baseline_prodd0 = 0.994

fit_total <- gformula(obs_data =dta,
                      id = id,
                      time_name = time_name,
                      outcome_name = outcome_name,
                      outcome_type = outcome_type,
                      covnames = covnames,
                      covtypes = covtypes,
                      histories = histories,
                      histvars = histvars,
                      covparams = covparams,
                      compevent_name = "dead_other",
                      compevent_model = compevent_model,
                      compevent_cens = FALSE,
                      censor_name = censor_name,
                      censor_model = censor_model,
                      ymodel = ymodel,
                      baseline_prodp0 = baseline_prodp0,
                      baseline_poprisk = baseline_poprisk,
                      baseline_prodd0 = baseline_prodd0,
                      basecovs = basecovs,
                      seed = 1,
                      #sim_data_b = T,
                      parallel= TRUE,
                      ncores = ncores,
                      nsamples = 500,
                      
                      intvars = list(
                        c("act","totfv", "totgrn", "junk", "totred", "totproc",
                          "sodajuice", "totalc")
                      ),
                      
                      int_descript = c(
                        'Activity and diet intervention'
                      ),
                      interventions = list(
                        list(c(xcond_int, 7.5, Inf),
                             c(threshold, 5, Inf),
                             c(threshold, 3, Inf),
                             c(threshold, -Inf, 1),
                             c(threshold, -Inf, 3),
                             c(threshold, -Inf, 0.999),
                             c(threshold, -Inf, 0),
                             c(threshold, -Inf, 0))
                      ),
                      restrictions = restrictions,
                      int_visit_type = rep(FALSE, 10),
                      covfits_custom = covfits_custom,
                      covpredict_custom = covpredict_custom
)
fit[[1]] <- fit_total
rm(list = "fit_total")


##########################################
##########################################
######################################
# advanced prostate cancer
# read in advanced data
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/lag_mlp/wcrf_fg_1115_adv.sas7bdat")

#  extract prebaseline, baseline variables, and
# conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)

ncores <- 20

baseline_poprisk <- 0.001831560
baseline_prodp0 <- 1- baseline_poprisk
fit_adv <- gformula(obs_data =dta,
                    id = id,
                    time_name = time_name,
                    outcome_name = outcome_name,
                    outcome_type = outcome_type,
                    covnames = covnames,
                    covtypes = covtypes,
                    histories = histories,
                    histvars = histvars,
                    covparams = covparams,
                    compevent_name = "dead_other",
                    compevent_model = compevent_model,
                    compevent_cens = FALSE,
                    censor_name = censor_name,
                    censor_model = censor_model,
                    ymodel = ymodel,
                    baseline_prodp0 = baseline_prodp0,
                    baseline_poprisk = baseline_poprisk,
                    baseline_prodd0 = baseline_prodd0,

                    basecovs = basecovs,
                    seed = 1,
                    #sim_data_b = T,
                    parallel= TRUE,
                    ncores = ncores,
                    nsamples = 500,
                    
                    intvars = list(
                      c("act","totfv", "totgrn", "junk", "totred", "totproc",
                        "sodajuice", "totalc")
                    ),
                    
                    int_descript = c(
                      'Activity and diet intervention'
                    ),
                    interventions = list(
                      list(c(xcond_int, 7.5, Inf),
                           c(threshold, 5, Inf),
                           c(threshold, 3, Inf),
                           c(threshold, -Inf, 1),
                           c(threshold, -Inf, 3),
                           c(threshold, -Inf, 0.999),
                           c(threshold, -Inf, 0),
                           c(threshold, -Inf, 0))
                    ),
                    restrictions = restrictions,
                    int_visit_type = rep(FALSE, 10),
                    covfits_custom = covfits_custom,
                    covpredict_custom = covpredict_custom
)
fit[[2]] <- fit_adv
rm(list = "fit_adv")


saveRDS(fit, "gformula_lag1_term_tot_adv.rds")


