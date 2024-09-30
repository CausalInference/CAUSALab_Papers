################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/sens_xcond3/read_gformula_results.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-28


# 1) Purpose: this R program helped read g-formula results for the sensitivity analysis 
# that expands the set of health conditions that excused individuals from following 
# the physical activity intervention to 
# include cancer (except prostate and nonmelanoma skin cancer), myocardial infarction, stroke, 
# congestive heart failure, functional impairment, angina pectoris, pulmonary embolism, diabetes, 
# chronic renal failure, rheumatoid arthritis, gout, ulcerative colitis/Crohn disease, and Parkinson disease.
# 2) Outcome: 1. total prostate cancer 
#             2. advanced prostate cancer 
#             3. lethal prostate cancer
#             4. fatal prostate cancer 

# 3) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_xcond3/gformula_xcond3_tot_adv.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/sens_xcond3/gformula_xcond3_lethal_fatal.rds"
# 4) Study design: Prospective cohort
#
# 5) Follow-up: HPFS: 1990-2016 (baseline : 1990 ; pre-baseline: 1986)


library(tidyverse)
library(gfoRmula)
options(scipen = 9999,
        digits = 4)

# total pca and advanced pca
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_xcond3/gformula_xcond3_tot_adv.rds")
fit

##################################
# lethal pca and fatal pca
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/sens_xcond3/gformula_xcond3_lethal_fatal.rds")
fit 
