################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/sens_lag1_term/read_gformula_results.R
# Programer: Fuyu Guo
# Last updated: 2024-Mar-24

# 1) Purpose: this R program helped read g-formula results for the sensitivity analysis 
# that conducted a minimum latency period analysis of 2 years 
# 2) Outcome: 1. total prostate cancer 
#             2. advanced prostate cancer 
#             3. lethal prostate cancer
#             4. fatal prostate cancer 

# 3) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_lag1_term/gformula_dynamic_tot_adv.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_lag1_term/gformula_dynamic_lethal_fatal.rds"
# 4) Study design: Prospective cohort
#
# 5) Follow-up: HPFS: 1992-2016

options(warn = -1,
        scipen = 9999,
        digits = 4)
library(tidyr)
library(dplyr)
library(readr)
library(gfoRmulaMLP)

#read total and advance pca g-formula results
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_R3/sens_lag1_term/gformula_lag1_term_tot_adv.rds")
fit

# read lethal and fatal pca g-formula results
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_R3/sens_lag1_term/gformula_lag1_term_lethal_fatal.rds")
fit
