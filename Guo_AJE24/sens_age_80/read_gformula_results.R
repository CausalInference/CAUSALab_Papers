# CAUSALAB
library(tidyverse)
library(gfoRmula)
options(scipen = 9999,
        digits = 4)

###############################################
# total prostate cancer and advanced
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/analysis_week_1115/sens_age_80/gformula_sens_age_80_tot_adv.rds")
fit

##############################################
# lethal and fatal prostate cancer 
fit <- read_rds("/n/hpnh/Users/fyguo/R_proj/analysis_week_1115/sens_age_80/gformula_sens_age_80_lethal_fatal.rds")
fit
