################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/sens_low_adherence_total_effects/gformula_low_adherence_tot.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-27
# Please note: the majority of the R programs were approximately the same among analyses for
# total prostate cancer, advanced prostate cancer, lethal prostate cancer, and fatal prostate cancer
# We loaded just different SAS datasets for each outcome. Thus, we will have four R programs whose
# names marked their outcomes.


# 1) Purpose: this R program used the parametric g-formula to estimate 26-yr
# risks under low-adherence interventions based on the WCRF 2018 recommendations 
# 2) Outcome: lethal prostate cancer (primary analysis for total effects comparing the risks under the 
# intervention to the risks under the low-adherence intervention) 
# 3) input files:
#     - SAS dataset: "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_lethal.sas7bdat"
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
setwd("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y")


# Use supporting R codes
# This R program is used to load functions "extract_prebase_base()" and "trans_var()" which
# will be used to clean data and makes it ready for g-formula
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_analyses_clean.R")

# This R program aims to create a customized function to be used by g-formula
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_cum.R")

#############################################################################
# read in data 
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_lethal.sas7bdat")

#  extract prebaseline, baseline variables, and
# conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)


ncores <- 50
# This R program stated the detailed g-formula settings used in the analysis
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_setting_restriction.R")
###################
# All the above settings are the same as the codes for the primary analyses.
# Bellow we set intervention and reference group differently.
# We ran each intervention strategy sequentially.

# Physical activity
fit_act <- gformula(obs_data =dta,
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
                    parallel= TRUE,
                    ncores = ncores,
                    nsamples = 500,
                    
                    intvars = list(c("act"),
                                   c("act")
                                   
                    ),
                    
                    int_descript = c("total activity",
                                     "total activity extreme"
                    ),
                    
                    interventions = list(list(c(xcond_int, 7.5, Inf)),
                                         list(c(threshold, 0, 7.49))
                    ),
                    ref_int = 2,
                    restrictions = restrictions,
                    int_visit_type = rep(FALSE, 10),
                    covfits_custom = covfits_custom,
                    covpredict_custom = covpredict_custom
)

###################################
#Fruits and vegetables

fit_fv <- gformula(obs_data =dta,
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
                   parallel= TRUE,
                   ncores = ncores,
                   nsamples = 500,
                   
                   intvars = list(
                     c("totfv"),
                     c("totfv")
                   ),
                   
                   int_descript = c(
                     "total fruit and vegetables",
                     "total fruit and vegetables"
                   ),
                   
                   interventions = list(
                     list(c(threshold, 5, Inf)),
                     list(c(threshold, 0, 4.99))
                   ),
                   ref_int = 2,
                   restrictions = restrictions,
                   int_visit_type = rep(FALSE, 10),
                   covfits_custom = covfits_custom,
                   covpredict_custom = covpredict_custom
)

###########################################
# Whole grains and legumes

fit_grn <- gformula(obs_data =dta,
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
                    parallel= TRUE,
                    ncores = ncores,
                    nsamples = 500,
                    
                    intvars = list(
                      c("totgrn"),
                      c("totgrn")
                    ),
                    
                    int_descript = c(
                      "total grains",
                      "total grains extreme"
                    ),
                    
                    interventions = list(
                      list(c(threshold, 3, Inf)),
                      list(c(threshold, 0, 2.99))
                    ),
                    ref_int = 2,
                    
                    restrictions = restrictions,
                    int_visit_type = rep(FALSE, 10),
                    covfits_custom = covfits_custom,
                    covpredict_custom = covpredict_custom
)

##############################################
# Processed food

fit_junk <- gformula(obs_data =dta,
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
                     parallel= TRUE,
                     ncores = ncores,
                     nsamples = 500,
                     
                     intvars = list(
                       c("junk"),
                       c("junk")
                       
                     ),
                     
                     int_descript = c(
                       "junk food",
                       "junk food extreme"
                       
                     ),
                     
                     interventions = list(
                       
                       list(c(threshold, -Inf, 1)),
                       list(c(threshold, 1.01, Inf))
                       
                     ),
                     ref_int = 2,
                     
                     restrictions = restrictions,
                     int_visit_type = rep(FALSE, 10),
                     covfits_custom = covfits_custom,
                     covpredict_custom = covpredict_custom
)

###########################################
# Red meat

fit_red <-  gformula(obs_data =dta,
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
                     parallel= TRUE,
                     ncores = ncores,
                     nsamples = 500,
                     
                     intvars = list(
                       c("totred"),
                       c("totred")
                     ),
                     
                     int_descript = c(
                       "total red meat",
                       "tot red meat extreme"
                     ),
                     
                     interventions = list(
                       list(c(threshold, -Inf, 3)),
                       list(c(threshold, 3.01, Inf))
                     ),
                     ref_int = 2,
                     
                     restrictions = restrictions,
                     int_visit_type = rep(FALSE, 10),
                     covfits_custom = covfits_custom,
                     covpredict_custom = covpredict_custom
)

#########################################
# Processed meat

fit_proc <- gformula(obs_data =dta,
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
                     parallel= TRUE,
                     ncores = ncores,
                     nsamples = 500,
                     
                     intvars = list(
                       c("totproc"),
                       c("totproc")
                     ),
                     
                     int_descript = c(
                       "total processed meat",
                       "total processed meat extreme"
                     ),
                     
                     interventions = list(
                       
                       list(c(threshold, -Inf, 0.99)),
                       list(c(threshold, 1, Inf))
                       
                     ),
                     ref_int = 2,
                     
                     restrictions = restrictions,
                     int_visit_type = rep(FALSE, 10),
                     covfits_custom = covfits_custom,
                     covpredict_custom = covpredict_custom
)

###########################################
# Sugar-sweetened beverages

fit_soda <- gformula(obs_data =dta,
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
                     parallel= TRUE,
                     ncores = ncores,
                     nsamples = 500,
                     
                     intvars = list(
                       c("sodajuice"),
                       c("sodajuice")
                     ),
                     
                     int_descript = c(
                       "total sodajuice",
                       "total sodajuice"
                       
                     ),
                     
                     interventions = list(
                       
                       list(c(threshold, -Inf, 0)),
                       list(c(threshold, 0.01, Inf))
                       
                       
                     ),
                     ref_int = 2,
                     
                     restrictions = restrictions,
                     int_visit_type = rep(FALSE, 10),
                     covfits_custom = covfits_custom,
                     covpredict_custom = covpredict_custom
)

#####################################
# Total alcohol 

fit_alc <- gformula(obs_data =dta,
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
                    parallel= TRUE,
                    ncores = ncores,
                    nsamples = 500,
                    
                    intvars = list(
                      c("totalc"),
                      c("totalc")
                    ),
                    
                    int_descript = c(
                      "total alcohol",
                      "total alcohol extreme"
                    ),
                    
                    interventions = list(
                      
                      list(c(threshold, -Inf, 0)),
                      list(c(threshold, 0.01, Inf))
                      
                    ),
                    ref_int = 2,
                    
                    
                    restrictions = restrictions,
                    int_visit_type = rep(FALSE, 10),
                    covfits_custom = covfits_custom,
                    covpredict_custom = covpredict_custom
)

fit <- list(fit_act, fit_fv, fit_grn, fit_junk,
            fit_red, fit_proc, fit_soda, fit_alc)

saveRDS(fit, "gformula_lethal_low_adherence.rds")
Sys.time()

