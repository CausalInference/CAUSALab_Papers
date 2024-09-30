################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/sens_dynamic_strategy/gformula_dynamic_tot_adv.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-27
# Please note: this R code used the same model settings as the one for the primary analysis
# we only changed the intervention setting for the dietary interventions from static to 
# dynamic which is denoted by xcond_int() function

# 1) Purpose: this R program used the parametric g-formula to estimate 26-yr
# risks under dynamic interventions based on the WCRF 2018 recommendations 
# 2) Outcome: total prostate cancer; advanced prostate cancer
# 3) input files:
#     - SAS dataset: 
#         "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_tot.sas7bdat"
#         "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_adv.sas7bdat"
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
.libPaths("/n/home00/fyguo/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyr)
library(readr)
library(dplyr)
library(haven)
library(data.table)
library(gfoRmula)
library(parallel)
setwd("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_dynamic_strategy")


# Use supporting R codes
# This R program is used to load functions "extract_prebase_base()" and "trans_var()" which
# will be used to clean data and makes it ready for g-formula
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_analyses_clean.R")

# This R program aims to create a customized function to be used by g-formula
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_cum.R")


# make a list to store g-formula results
fit <- list(2L)
######################################
# total prostate cancer
# read in total data
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_tot.sas7bdat")

#  extract prebaseline, baseline variables, and
# conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)

# This R program stated the detailed g-formula settings used in the analysis
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_setting_restriction.R")
ncores <- 50
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
                      # here we set interventions to be xcond_int() for dietary interventions,
                      # which are the only difference from the primary analysis
                      interventions = list(
                        list(c(xcond_int, 7.5, Inf),
                             c(xcond_int, 5, Inf),
                             c(xcond_int, 3, Inf),
                             c(xcond_int, -Inf, 1),
                             c(xcond_int, -Inf, 3),
                             c(xcond_int, -Inf, 0.999),
                             c(xcond_int, -Inf, 0),
                             c(xcond_int, -Inf, 0))
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
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_adv.sas7bdat")

#  extract prebaseline, baseline variables, and
# conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)

# This R program stated the detailed g-formula settings used in the analysis
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_setting_restriction.R")

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
                             c(xcond_int, 5, Inf),
                             c(xcond_int, 3, Inf),
                             c(xcond_int, -Inf, 1),
                             c(xcond_int, -Inf, 3),
                             c(xcond_int, -Inf, 0.999),
                             c(xcond_int, -Inf, 0),
                             c(xcond_int, -Inf, 0))
                      ),
                      restrictions = restrictions,
                      int_visit_type = rep(FALSE, 10),
                      covfits_custom = covfits_custom,
                      covpredict_custom = covpredict_custom
)
fit[[2]] <- fit_adv
rm(list = "fit_adv")

Sys.time()

saveRDS(fit, "gformula_dynamic_tot_adv.rds")


