################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/check_%_to_be_intervened_on.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-26

# 1) Purpose: this R program estimated the proportion of people to be intervened on during
# the 26-yr follow-up
# 2) Outcome: 1. total prostate cancer (primary analysis for fatal prostate cancer risks) 

# 3) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_tot_1115.rds"# 4) Study design: Prospective cohort
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
# Here we calculate the % to be intervened on at baseline
# As the four outcome datasets were the same at baseline, we read in total prostate cancer file to 
# calculate the baseline % be intervened on
# read in data 
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_tot.sas7bdat")
#  extract prebaseline, baseline variables, and conduct transformations 
dta <- dta  %>%
  extract_prebase_base() %>%
  trans_var() 
dta$period_coarse <- coarsen(dta$period)
dta <- as.data.table(dta)

# extract baseline data
dta_base <- dta[period == 0,,]

# percentange of single components
dta_base[,.(per_activity = mean(act < 7.5),
            per_fruitAndvegetable = mean(totfv < 5),
            per_grainAndlegume = mean(totgrn < 3),
            per_processed_food = mean(junk > 1),
            per_red_meat = mean(totred > 3),
            per_processed_meat = mean(totproc > 1),
            per_SSB = mean(sodajuice > 0),
            per_alcohol = mean(totalc > 0)),
         ]

# % to be intervened on for the joint dietary intervention
cat("% for joint dietary intervention",
    mean(dta_base$totfv < 5|
           dta_base$totgrn < 3|
           dta_base$junk > 1|
           dta_base$totred > 3|
           dta_base$totproc > 1|
           dta_base$sodajuice > 0|
           dta_base$totalc > 0))

# % to be intervened on for the joint dietary intervention
cat("% for joint dietary intervention",
    mean(dta_base$act < 7.5|
           dta_base$totfv < 5|
           dta_base$totgrn < 3|
           dta_base$junk > 1|
           dta_base$totred > 3|
           dta_base$totproc > 1|
           dta_base$sodajuice > 0|
           dta_base$totalc > 0))

##########################################################
##########################################################
# Average % intervened on and cumulative % intervened on can
# be read directly from the gformula results
fit_tot <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_tot_1115.rds")
fit_tot