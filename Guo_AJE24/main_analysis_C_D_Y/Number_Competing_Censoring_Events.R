################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/Number_Competing_Censoring_Events.R
# Programer: Fuyu Guo
# Last updated: 2023-Aug-14

# 1) Purpose: this R program calcualted the number of competing events, cenosring, and outcomes during the
# 26-yr follow-up among 4 prostate cancer outcomes.
# 2) Outcome: lethal prostate cancer (primary analysis for lethal prostate cancer risks) 
# 3) input files:
#     - SAS dataset: 
#     "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_tot.sas7bdat"
#     "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_adv.sas7bdat"
#     "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_lethal.sas7bdat"
#     "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_fatal.sas7bdat"
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

#########################################################
# Dataset for total prostate cancer
# read in data 
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_tot.sas7bdat")
# number of follow-up years = number of rows * 2 as each row represents 2 years
cat("Number of follow-up years for total prostate cancer:", dim(dta)[1]*2)

# number of events
cat("Number of outcome events for total prostate cancer:", sum(dta$event, na.rm = T))

# number of competing events
cat("Number of competing events for total prostate cancer:", sum(dta$dead_other, na.rm = T))

# number of censoring
cat("Number of censoring for total prostate cancer:", sum(dta$censor, na.rm = T))

# number of administrative censoring
cat("Number of administrative censoring for total prostate cancer:", 
    27859-sum(dta$event, na.rm = T)-sum(dta$dead_other, na.rm = T) -  sum(dta$censor, na.rm = T))


#########################################################
#########################################################
# Dataset for advanced prostate cancer
# read in data 
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_adv.sas7bdat")
# number of follow-up years = number of rows * 2 as each row represents 2 years
cat("Number of follow-up years for advanced prostate cancer:", dim(dta)[1]*2)

# number of events
cat("Number of outcome events for advanced prostate cancer:", sum(dta$event, na.rm = T))

# number of competing events
cat("Number of competing events for advanced prostate cancer:", sum(dta$dead_other, na.rm = T))

# number of censoring
cat("Number of competing events for advanced prostate cancer:", sum(dta$censor, na.rm = T))

# number of administrative censoring
cat("Number of administrative censoring for advanced prostate cancer:", 
    27859-sum(dta$event, na.rm = T)-sum(dta$dead_other, na.rm = T) -  sum(dta$censor, na.rm = T))

#########################################################
#########################################################
# Dataset for lethal prostate cancer
# read in data 
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_lethal.sas7bdat")
# number of follow-up years = number of rows * 2 as each row represents 2 years
cat("Number of follow-up years for lethal prostate cancer:", dim(dta)[1]*2)

# number of events
cat("Number of outcome events for lethal prostate cancer:", sum(dta$event, na.rm = T))

# number of competing events
cat("Number of competing events for lethal prostate cancer:", sum(dta$dead_other, na.rm = T))

# number of censoring
cat("Number of competing events for lethal prostate cancer:", sum(dta$censor, na.rm = T))
cat("Number of administrative censoring for advanced prostate cancer:", 
    27859-sum(dta$event, na.rm = T)-sum(dta$dead_other, na.rm = T) -  sum(dta$censor, na.rm = T))
# number of administrative censoring

#########################################################
#########################################################
# Dataset for fatal prostate cancer
# read in data 
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_fatal.sas7bdat")
# number of follow-up years = number of rows * 2 as each row represents 2 years
cat("Number of follow-up years for lethal prostate cancer:", dim(dta)[1]*2)

# number of events
cat("Number of outcome events for lethal prostate cancer:", sum(dta$event, na.rm = T))

# number of competing events
cat("Number of competing events for lethal prostate cancer:", sum(dta$dead_other, na.rm = T))

# number of censoring
cat("Number of competing events for lethal prostate cancer:", sum(dta$censor, na.rm = T))

# number of administrative censoring
cat("Number of administrative censoring for lethal prostate cancer:", 
    27859-sum(dta$event, na.rm = T)-sum(dta$dead_other, na.rm = T) -  sum(dta$censor, na.rm = T))
