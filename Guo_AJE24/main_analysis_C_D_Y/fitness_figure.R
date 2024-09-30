################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_fitness_figure.R
# Programer: Fuyu Guo
# Last updated: 2023-Jul-1


# 1) Purpose: this R program read in g-formula results and made fittness plots to check any
# model misspecification and its influence for four outcomes.

# 2) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_tot_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_adv_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_lethal_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_fatal_1115.rds"
# 3) Study design: Prospective cohort
#
# 4) Follow-up: HPFS: 1990-2016 (baseline : 1990 ; pre-baseline: 1986)

# load packages 
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(ggsci)
library(data.table)
library(scales)
setwd("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y")



###########################
# Total prostate figure 
fit_tot <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_tot_1115.rds")
plot(fit_tot,
     covnames = fit_tot$covnames[-1],
     ylab_cov = c("Diabetes (%)",
                  "Hypertension (%)",
                  "Hypercholesterolemia \n (%)",
                  "Total energy intake \n (kcal/day)",
                  "Body-mass index \n (kg/m2)",
                  "Serious health condition \n (%)",
                  "Physical activity \n (MET-hours/week)",
                  "Fruits and non-starchy \n vegetables (servings/day)",
                  "Whole grain and legumes \n (servings/day)",
                  "Processed foods \n (servings/day)",
                  "Red meat \n (servings/week)",
                  "Processed meat \n (servings/week)",
                  "Sugar-sweetened beverages \n (servings/day)",
                  "Alcohol \n (servings/day)"),
     ylab_risk = "Risk (%)",
     ncol = 4)
ggsave("fitness_total.png",
    width = 300,
    height = 200,
    dpi = 300,
    bg = "white",
    units = "mm")


###########################
# Advanced prostate figure 

fit_adv <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_adv_1115.rds")
plot(fit_adv,
     covnames = fit_tot$covnames[-1],
     ylab_cov = c("Diabetes (%)",
                  "Hypertension (%)",
                  "Hypercholesterolemia \n (%)",
                  "Total energy intake \n (kcal/day)",
                  "Body-mass index \n (kg/m2)",
                  "Serious health condition \n (%)",
                  "Physical activity \n (MET-hours/week)",
                  "Fruits and non-starchy \n vegetables (servings/day)",
                  "Whole grain and legumes \n (servings/day)",
                  "Processed foods \n (servings/day)",
                  "Red meat \n (servings/week)",
                  "Processed meat \n (servings/week)",
                  "Sugar-sweetened beverages \n (servings/day)",
                  "Alcohol \n (servings/day)"),
     ylab_risk = "Risk (%)",
     ncol = 4)
ggsave("fitness_adv.png",
       width = 300,
       height = 200,
       dpi = 300,
       bg = "white",
       units = "mm")
###########################
# Lethal prostate figure 

fit_lethal <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_lethal_1115.rds")
plot(fit_lethal,
     covnames = fit_tot$covnames[-1],
     ylab_cov = c("Diabetes (%)",
                  "Hypertension (%)",
                  "Hypercholesterolemia \n (%)",
                  "Total energy intake \n (kcal/day)",
                  "Body-mass index \n (kg/m2)",
                  "Serious health condition \n (%)",
                  "Physical activity \n (MET-hours/week)",
                  "Fruits and non-starchy \n vegetables (servings/day)",
                  "Whole grain and legumes \n (servings/day)",
                  "Processed foods \n (servings/day)",
                  "Red meat \n (servings/week)",
                  "Processed meat \n (servings/week)",
                  "Sugar-sweetened beverages \n (servings/day)",
                  "Alcohol \n (servings/day)"),
     ylab_risk = "Risk (%)",
     ncol = 4)
ggsave("fitness_lethal.png",
       width = 300,
       height = 200,
       dpi = 300,
       bg = "white",
       units = "mm")
###########################
# Fatal prostate figure 

fit_fatal <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_fatal_1115.rds")
plot(fit_fatal,
     covnames = fit_tot$covnames[-1],
     ylab_cov = c("Diabetes (%)",
                  "Hypertension (%)",
                  "Hypercholesterolemia \n (%)",
                  "Total energy intake \n (kcal/day)",
                  "Body-mass index \n (kg/m2)",
                  "Serious health condition \n (%)",
                  "Physical activity \n (MET-hours/week)",
                  "Fruits and non-starchy \n vegetables (servings/day)",
                  "Whole grain and legumes \n (servings/day)",
                  "Processed foods \n (servings/day)",
                  "Red meat \n (servings/week)",
                  "Processed meat \n (servings/week)",
                  "Sugar-sweetened beverages \n (servings/day)",
                  "Alcohol \n (servings/day)"),
     ylab_risk = "Risk (%)",
     ncol = 4)
ggsave("fitness_fatal.png",
       width = 300,
       height = 200,
       dpi = 300,
       bg = "white",
       units = "mm")

Sys.time()
