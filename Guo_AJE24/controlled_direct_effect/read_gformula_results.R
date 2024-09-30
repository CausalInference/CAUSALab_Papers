################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/controlled_direct_effect/read_gformula_results
# Programer: Fuyu Guo
# Last updated: 2023-Jun-27


# 1) Purpose: this R program helped read g-formula results from the secondary analysis
# to estimate controlled direct effects
# 2) Outcome: 1. total prostate cancer 
#             2. advanced prostate cancer 
#             3. lethal prostate cancer 
#             4. fatal prostate cancer 

# 3) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/controlled_direct_effect/gformula_tot_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/controlled_direct_effect/gformula_adv_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/controlled_direct_effect/gformula_lethal_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/controlled_direct_effect/gformula_fatal_1115.rds"
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


# total pca
fit_tot <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/controlled_direct_effect/gformula_tot_1115.rds")
fit_tot


# advanced pca
fit_adv <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/controlled_direct_effect/gformula_adv_1115.rds")
fit_adv


# lethal pca
fit_lethal <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/controlled_direct_effect/gformula_lethal_1115.rds")
fit_lethal

# fatal pca
fit_fatal <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/controlled_direct_effect/gformula_fatal_1115.rds")
fit_fatal
 
Sys.time()