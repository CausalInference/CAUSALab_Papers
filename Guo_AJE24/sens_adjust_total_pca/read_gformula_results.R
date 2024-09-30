################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/sens_adjust_total_pca/read_gformula_results
# Programer: Fuyu Guo
# Last updated: 2023-Jun-30


# 1) Purpose: this R program helped read g-formula results for the sensitivity analysis that
# additionally adjusted for PSA testing
# 2) Outcome: 1. lethal prostate cancer
#             2. fatal prostate cancer 

# 3) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_adjust_total_pca/gformula_lethal_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_adjust_total_pca/gformula_fatal_1115.rds"
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
# lethal prostate cancer
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_adjust_total_pca/gformula_lethal_1115.rds")
fit

##############################################
# fatal prostate cancer
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_adjust_total_pca/gformula_fatal_1115.rds")
fit

Sys.time()