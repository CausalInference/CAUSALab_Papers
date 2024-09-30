################################################################
# CAUSALAB 
# table2_baseline_characteristics
# Programer: Fuyu Guo
# Last updated: 2023-Jun-25

# 1) Purpose: this R program will describe baseline demographic features
# 2) input files:
#     - SAS dataset: "/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_fatal.sas7bdat"

################################################################
# R settings
.libPaths("/n/home00/fyguo/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyr)
library(dplyr)
library(haven)
library(data.table)
library(tableone)

# Read in SAS data
#  - We use total prostate cancer to describe baseline characteristics.
#  - All the four datasets have no difference at period 0.
dta <- read_sas("/n/hpnh/Users/fyguo/data/SAS_code_review_21Jun2023/competing_data_structure_C_D_Y/wcrf_fg_1115_tot.sas7bdat")


# dta is the baseline data
dta <- dta %>% filter(period == 0) %>%
  as.data.table()

# Make a demographic table

myVars <- c("age", "bmi", "act", 
            "totfv", "totgrn", "junk", "totred", "totproc", "sodajuice", "totalc", 
            "cal", "dia", "hbp", "chl", "fhxpca90", "fhxmi", "smkhx", "psa", "hlthy")
catVars <- c("dia", "hbp", "chl", "fhxpca90", "fhxmi", "smkhx", "psa", "hlthy")
tab <- CreateTableOne(vars = myVars, 
                      data = dta, 
                      factorVars = catVars)
print(tab, nonnormal = T, formatOptions = list(big.mark = ","))