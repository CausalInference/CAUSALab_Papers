
########################################################################################################################################
#   
# CAUSALab 
# 
# File...........: Baseline characteristics and standardised mean differences of people included in the analysis 
#                  in ART-naive and ART-experienced individuals
# 
# Project........: Integrase strand-transfer inhibitor use and cardiovascular events in adults with HIV: 
#                  an emulation of target trials in the HIV-CAUSAL Collaboration 
#                  and the Antiretroviral Therapy Cohort Collaboration
# 
# Author.........: Sophia Rein
# 
# Date Created...: 2nd of May 2023
# 
# *------------------------------------------------------------------------------------------------------------
#   
#   Purpose........: Analysis for paper / CAUSALab transparency initiative
# 
# *------------------------------------------------------------------------------------------------------------
# Modification: Added Standardized mean differences
# 
# Date: 27th of July 2023  		Author: Sophia Rein	
# 
########################################################################################################################################
if (!require("data.table")) install.packages("data.table")
library(data.table)
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if (!require("tableone")) install.packages("tableone")
library(tableone)

insti <- fread(file="path to file")

################################################################################################################
## ART-naive individuals
naive <- insti %>%
  filter(elig_init==1)
## list of covariates
vars1 <- c("sex", "age", "new_ethnic", "new_mode", "cd4_v", "rna_v", "aids", "hcv_cat", "hbv_cat",
           "bmi_cat2", "hypertension", "total_chol2", "current_smoke_v", "diabetes_new", "ckd_v",
           "has_ABC", "bldate", "cohort")
# list of categorical covariates
cats1 <- c("sex", "new_ethnic", "new_mode", "aids", "hcv_cat", "hbv_cat", 
           "bmi_cat2", "hypertension", "total_chol2", "current_smoke_v",
           "diabetes_new", "has_ABC", "ckd_v", "has_ABC", "bldate", "cohort")
## construct the table
tab1 <- CreateTableOne(vars = vars1, strata = "has_insti", data = naive, factorVars=cats1, test = FALSE)
print(tab1,showAllLevels = TRUE, nonnormal=TRUE, smd = TRUE)

################################################################################################################
## ART-experienced individuals
experienced <- insti %>%
  filter(elig_switch==1)
## list of covariates
vars2 <- c("sex", "age", "new_ethnic", "new_mode", "cd4_v", "tsinceart", "aids", "hcv_cat", "hbv_cat",
           "bmi_cat2", "hypertension", "total_chol2", "current_smoke_v", "diabetes_new", "ckd_v", 
           "abc_6mo", "bldate", "cohort")
# list of categorical covariates
cats2 <- c("sex", "new_ethnic", "new_mode", "aids", "hcv_cat", "hbv_cat", 
           "bmi_cat2", "hypertension", "total_chol2", "current_smoke_v",
           "diabetes_new", "ckd_v", "abc_6mo", "bldate", "cohort")
## construct the table
tab2 <- CreateTableOne(vars = vars2, strata = "has_insti", data = experienced, factorVars=cats2, test = FALSE)
print(tab2,showAllLevels = TRUE, nonnormal=TRUE, smd = TRUE)
