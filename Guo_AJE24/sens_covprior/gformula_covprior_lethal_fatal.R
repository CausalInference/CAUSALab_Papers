################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/gformula_covprior_lethal_fatal.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-30

# 1) Purpose: this R program used the parametric g-formula to estimate 26-yr
# risks on the WCRF 2018 recommendations for the sensitivity analysis that 
# adjusted for covariates from the previous questionnaire rather than 
# the questionnaire concurrent with diet and physical activity 
# 2) Outcome: lethal prostate cancer; fatal prostate cancer
# 3) input files:
#     - SAS dataset: 
#         "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/covprior/wcrf_fg_1115_lethal.sas7bdat"
#         "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/covprior/wcrf_fg_1115_fatal.sas7bdat"
#     - R codes:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/gformula_analyses_clean.R"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_cum.R"
#        ""/n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/gformula_setting_restriction.R""
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
setwd("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior")


# Use supporting R codes
# This R program is used to load functions "extract_prebase_base()" and "trans_var()" which
# will be used to clean data and makes it ready for g-formula
source("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/gformula_analyses_clean.R")

# This R program aims to create a customized function to be used by g-formula
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_cum.R")

fit <- list(2L)
######################################
# lethal prostate cancer
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/covprior/wcrf_fg_1115_lethal.sas7bdat")

#  extract prebaseline, baseline variables, and
# conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)


source("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/gformula_setting_restriction.R")
ncores <- 50
fit_lethal <- gformula(obs_data =dta,
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
fit[[1]] <- fit_lethal
rm(list = "fit_lethal")


##########################################
##########################################
######################################
# fatal prostate cancer
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/covprior/wcrf_fg_1115_fatal.sas7bdat")

#  extract prebaseline, baseline variables, and
# conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)

ncores <- 50
fit_fatal <- gformula(obs_data =dta,
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
fit[[2]] <- fit_fatal
rm(list = "fit_fatal")

Sys.time()
saveRDS(fit, "gformula_sens_covprior_lethal_fatal.rds")


