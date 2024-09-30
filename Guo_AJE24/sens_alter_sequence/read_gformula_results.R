################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/sens_alter_sequence/read_gformula_results
# Programer: Fuyu Guo
# Last updated: 2023-Jul-1


# 1) Purpose: this R program helped read g-formula results for the sensitivity analysis that
# altered the order of time-varying treatment covariates 
# 2) Outcome: 1. total prostate cancer 
#             2. advanced prostate cancer 
#             3. lethal prostate cancer
#             4. fatal prostate cancer 

# 3) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_alter_sequence/gformula_alter_sequence_tot_adv.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_alter_sequence/gformula_alter_sequence_lethal_fatal.rds"
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

#read total and advance pca g-formula results
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_alter_sequence/gformula_alter_sequence_tot_adv.rds")
fit

# read lethal and fatal pca g-formula results
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_alter_sequence/gformula_alter_sequence_lethal_fatal.rds")
fit

Sys.time()