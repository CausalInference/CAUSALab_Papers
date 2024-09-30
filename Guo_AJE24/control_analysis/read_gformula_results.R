################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/control_analysis/read_gformula_results
# Programer: Fuyu Guo
# Last updated: 2023-Jul-1


# 1) Purpose: this R program helped read g-formula results for the secondary analysis
# using negative outcome control as well as positive outcome control
# 2) Outcome: 1. accidental or injury-related deaths 
#             2. caridovascular disease deaths 


# 3) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_fun_splines/gformula_fun_splines_tot_adv.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_fun_splines/gformula_fun_splines_lethal_fatal.rds"
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

# Negative outcome control: accidental or injury-related death

fit_acd <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/control_analysis/gformula_acd_1115.rds")
fit_acd

# Positive outcome control: cardiovascular disease deaths

fit_cvd <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/control_analysis/gformula_cvd_1115.rds")
fit_cvd

Sys.time()
