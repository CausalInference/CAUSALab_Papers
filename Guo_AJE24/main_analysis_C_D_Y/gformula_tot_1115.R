################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/main_analysis_C_D_Y/gformula_tot_1115.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-24
# Please note: the majority of the R programs were approximately the same among analyses for
# total prostate cancer, advanced prostate cancer, lethal prostate cancer, and fatal prostate cancer
# We loaded just different SAS datasets for each outcome. Thus, we will have four R programs whose
# names marked their outcomes.


# 1) Purpose: this R program used the parametric g-formula to estimate 26-yr
# risks under interventions based on the WCRF 2018 recommendations 
# 2) Outcome: total prostate cancer (primary analysis for total effects) 
# 3) input files:
#     - SAS dataset: "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_tot.sas7bdat"
#     - R codes:
#        "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/main_analysis_C_D_Y/gformula_analyses_clean.R"
#        "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/main_analysis_C_D_Y/gformula_cum.R"
#        "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/main_analysis_C_D_Y/gformula_setting_restriction.R""
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
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_tot.sas7bdat")
#  extract prebaseline, baseline variables, and conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)
############################################################

# set number of cores used in the parallel computations
ncores <- 50

# This R program stated the detailed g-formula settings used in the analysis
source("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_setting_restriction.R")


# Apply g-formula
fit <- gformula(obs_data =dta,
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
                      threads = 1,
                      
                      intvars = list(c("act"), 
                                     c("totfv"),
                                     c("totgrn"),
                                     c("junk"),
                                     c("totred"), 
                                     c("totproc"),
                                     c("sodajuice"),
                                     c("totalc"),
                                     c("totfv", "totgrn", "junk", "totred", "totproc",
                                       "sodajuice", "totalc"),
                                     c("act","totfv", "totgrn", "junk", "totred", "totproc",
                                       "sodajuice", "totalc")
                      ),
                      
                      int_descript = c("total activity",
                                       "total fruit and vegetables",
                                       "total grains",
                                       "junk food",
                                       "total red meat", 
                                       "total processed meat",
                                       "total sodajuice",
                                       "total alchol",
                                       'Diet intervention', 
                                       'Activity and diet intervention'
                      ),
                      interventions = list(list(c(xcond_int, 7.5, Inf)),
                                           
                                           list(c(threshold, 5, Inf)),
                                           
                                           list(c(threshold, 3, Inf)),
                                           
                                           list(c(threshold, -Inf, 1)),
                                           
                                           list(c(threshold, -Inf, 3)),
                                           
                                           list(c(threshold, -Inf, 0.999)),
                                           
                                           list(c(threshold, -Inf, 0)),
                                           
                                           list(c(threshold, -Inf, 0)),
                                           
                                           list(c(threshold, 5, Inf),
                                                c(threshold, 3, Inf),
                                                c(threshold, -Inf, 1),
                                                c(threshold, -Inf, 3),
                                                c(threshold, -Inf, 0.999),
                                                c(threshold, -Inf, 0),
                                                c(threshold, -Inf, 0)),
                                           
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


# show systematic time
Sys.time()

# Save the g-formula results for the future use
saveRDS(fit, "gformula_tot_1115.rds")
