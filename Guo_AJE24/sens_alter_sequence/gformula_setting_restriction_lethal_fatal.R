#################################################
# CAUSALAB
# this file is to do sensitivity analyses for changing time-varying treatments' sequence
# for the joint intervention


library(tidyverse)
library(haven)
library(tictoc)
library(gfoRmula)
library(data.table)
setwd("/n/hpnh/Users/fyguo/R_proj/analysis_week_1115/sens_alter_sequence")
source("/n/hpnh/Users/fyguo/R_proj/analysis_week_1115/main_analysis_C_D_Y/gformula_analyses_clean.R")
source("/n/hpnh/Users/fyguo/R_proj/analysis_week_1115/main_analysis_C_D_Y/gformula_cum.R")

fit <- list(2L)
######################################
# total prostate cancer
# read in total data
dta <- read_sas("/n/hpnh/Users/fyguo/data/generate_data_1115/competing_data_structure_C_D_Y/wcrf_fg_1115_lethal.sas7bdat")


#  extract prebaseline, baseline variables, and
# conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)


source("/n/hpnh/Users/fyguo/R_proj/analysis_week_1115/sens_alter_sequence/gformula_setting_restriction.R")
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
                      sim_data_b = T,
                      #parallel= TRUE,
                      #ncores = ncores,
                      #nsamples = 500,
                      
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
toc()
fit[[1]] <- fit_lethal
rm(list = "fit_lethal")

##########################################
##########################################
######################################
# advanced prostate cancer
# read in advanced data
dta <- read_sas("/n/hpnh/Users/fyguo/data/generate_data_1115/competing_data_structure_C_D_Y/wcrf_fg_1115_fatal.sas7bdat")

#  extract prebaseline, baseline variables, and
# conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)



tic()
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
                    sim_data_b = T,
                    #parallel= TRUE,
                    #ncores = ncores,
                    #nsamples = 500,
                    
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
toc()
fit[[2]] <- fit_fatal
rm(list = "fit_fatal")



saveRDS(fit, "gformula_alter_sequence_lethal_fatal.rds")


