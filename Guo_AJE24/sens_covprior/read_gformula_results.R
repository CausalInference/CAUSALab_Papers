################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/read_gformula_results
# Programer: Fuyu Guo
# Last updated: 2023-Jul-1


# 1) Purpose: this R program helped read g-formula results for the sensitivity analysis that
# adjusted for covariates from the previous questionnaire 
# rather than the questionnaire concurrent with diet and physical activity
# 2) Outcome: 1. total prostate cancer 
#             2. advanced prostate cancer 
#             3. lethal prostate cancer
#             4. fatal prostate cancer 

# 3) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/gformula_sens_covprior_tot_adv.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/gformula_sens_covprior_lethal_fatal.rds"
# 4) Study design: Prospective cohort
#
# 5) Follow-up: HPFS: 1990-2016 (baseline : 1990 ; pre-baseline: 1986)


options(warn = -1,
        scipen = 9999,
        digits = 4)
.libPaths("/n/home00/fyguo/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyr)
library(dplyr)
library(readr)
library(gfoRmula)
###############################################
# total prostate cancer and advanced
#read total and advance pca g-formula results
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/gformula_sens_covprior_tot_adv.rds")
fit

# read lethal and fatal pca g-formula results
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_covprior/gformula_sens_covprior_lethal_fatal.rds")
fit

Sys.time()
