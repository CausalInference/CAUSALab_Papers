################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_figure_RD.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-25


# 1) Purpose: this R program read in g-formula results and made forest plots to show 26-ur risk
# differences for total, advanced, lethal, and fatal prostate cancer.
# 2) Outcome: 1. total prostate cancer (primary analysis for fatal prostate cancer risks) 
#             2. advanced prostate cancer (primary analysis for fatal prostate cancer risks)
#             3. lethal prostate cancer (primary analysis for fatal prostate cancer risks)
#             4. fatal prostate cancer (primary analysis for fatal prostate cancer risks)

# 3) input files:
#     - R data storing the g-formula results:
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_tot_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_adv_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_lethal_1115.rds"
#        "/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_fatal_1115.rds"
# 4) Study design: Prospective cohort
#
# 5) Follow-up: HPFS: 1990-2016 (baseline : 1990 ; pre-baseline: 1986)

# load packages 
.libPaths("/n/home00/fyguo/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(ggsci)
library(gfoRmula)
library(data.table)
library(scales)
setwd("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y")

# read in total pca results
fit_tot <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_tot_1115.rds")

# extract essential results
df_tot <- fit_tot$result[k==12 & `Interv.` %in% 1:10,.(`Risk difference`,  `RD lower 95% CI`, `RD upper 95% CI`)]
df_tot[, `:=`(intervention = 
                c("Physical activity",
                  "Fruits and non-starchy vegetables",
                  "Whole grain and legumes",
                  "Processed foods high in fat, starches, or sugar",
                  "Red meat",
                  "Processed meat",
                  "Sugar-sweetened beverages",
                  "Alcohol",
                  "Joint intervention on diet",
                  "Joint intervention on diet and physical activity")) ]

colnames(df_tot)[c(1,2, 3)] <- c('RD', "RD_low_ci", "RD_up_ci")
df_tot <- df_tot[,c(4,1,2,3)]
df_tot$outcome <- "Total Prostate Cancer"

#########################################
# advanced pca
fit_adv <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_adv_1115.rds")
df_adv <- fit_adv$result[k==12 & Interv. %in% 1:10,.(`Risk difference`,  `RD lower 95% CI`, `RD upper 95% CI`)]
df_adv[, `:=`(intervention = 
                c("Physical activity",
                  "Fruits and non-starchy vegetables",
                  "Whole grain and legumes",
                  "Processed foods high in fat, starches, or sugar",
                  "Red meat",
                  "Processed meat",
                  "Sugar-sweetened beverages",
                  "Alcohol",
                  "Joint intervention on diet",
                  "Joint intervention on diet and physical activity")) ]

colnames(df_adv)[c(1,2, 3)] <- c('RD', "RD_low_ci", "RD_up_ci")
df_adv <- df_adv[,c(4,1,2,3)]
df_adv$outcome <- "Advanced Prostate Cancer"

##################################################
# lethal pca
fit_lethal <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_lethal_1115.rds")
fit_lethal
df_lethal <- fit_lethal$result[k==12 & Interv. %in% 1:10,.(`Risk difference`,  `RD lower 95% CI`, `RD upper 95% CI`)]
df_lethal[, `:=`(intervention = 
                   c("Physical activity",
                     "Fruits and non-starchy vegetables",
                     "Whole grain and legumes",
                     "Processed foods high in fat, starches, or sugar",
                     "Red meat",
                     "Processed meat",
                     "Sugar-sweetened beverages",
                     "Alcohol",
                     "Joint intervention on diet",
                     "Joint intervention on diet and physical activity")) ]

colnames(df_lethal)[c(1,2, 3)] <- c('RD', "RD_low_ci", "RD_up_ci")
df_lethal <- df_lethal[,c(4,1,2,3)]
df_lethal$outcome <- "Lethal Prostate Cancer"

##################################################
# fatal pca
fit_fatal <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_fatal_1115.rds")
df_fatal <- fit_fatal$result[k==12 & Interv. %in% 1:10, .(`Risk difference`,  `RD lower 95% CI`, `RD upper 95% CI`)]
df_fatal[, `:=`(intervention = 
                  c("Physical activity",
                    "Fruits and non-starchy vegetables",
                    "Whole grain and legumes",
                    "Processed foods high in fat, starches, or sugar",
                    "Red meat",
                    "Processed meat",
                    "Sugar-sweetened beverages",
                    "Alcohol",
                    "Joint intervention on diet",
                    "Joint intervention on diet and physical activity")
                )]

colnames(df_fatal)[c(1,2, 3)] <- c('RD',"RD_low_ci", "RD_up_ci")
df_fatal <- df_fatal[,c(4,1,2,3)]
df_fatal$outcome <- "Fatal Prostate Cancer"

#######################################
# combined all the estimated data
df <- rbind(df_tot, df_adv, df_lethal, df_fatal)
df$outcome <- factor(df$outcome,
                     levels = c(
                       "Fatal Prostate Cancer",
                       "Lethal Prostate Cancer",
                       "Advanced Prostate Cancer",
                       "Total Prostate Cancer") %>%
                       rev())

# rename the intervention variable
df$intervention <- factor(df$intervention,
                          levels = c("Physical activity",
                                     "Fruits and non-starchy vegetables",
                                     "Whole grain and legumes",
                                     "Processed foods high in fat, starches, or sugar",
                                     "Red meat",
                                     "Processed meat",
                                     "Sugar-sweetened beverages",
                                     "Alcohol",
                                     "Joint intervention on diet",
                                     "Joint intervention on diet and physical activity") %>%
                            rev()) 
# draw the plot
ggplot(df)  +
  geom_point(aes(x=RD*100, 
                 y =intervention,
                 color = intervention),
             position = position_dodge2(width = 0.5))+
  geom_linerange(aes(x = RD*100, 
                     y = intervention, 
                     xmin = RD_low_ci*100, 
                     xmax = RD_up_ci*100, 
                     color = intervention),
                 position = position_dodge2(width = 0.5)) +
  facet_wrap(~outcome, 
             nrow = 2)+
  scale_y_discrete(labels = c("Physical activity: \u2265 7.5 MET-hours ",
                              "Fruits and non-starchy vegetables: \u2265 5 servings/day",
                              "Whole grain and legumes: \u2265 3 servings/day",
                              "Processed foods high in fat, starches, or sugar: \n \u2264 3 servings/day",
                              "Red meat: \u2264 3 servings/day",
                              "Processed meat: <1 servings/day",
                              "Sugar-sweetened beverages: 0 serving/day",
                              "Alcohol: 0 serving/day",
                              "Joint intervention on diet",
                              "Joint intervention on diet and physical activity") %>%
                     rev()) +
  scale_color_d3(name = "Intervention",
                     labels = c("Physical activity",
                                "Fruits and non-starchy vegetables",
                                "Whole grain and legumes",
                                "Processed foods high in fat, starches, or sugar",
                                "Red meat",
                                "Processed meat",
                                "Sugar-sweetened beverages",
                                "Alcohol",
                                "Joint intervention on diet",
                                "Joint intervention on diet and physical activity") %>%
                       rev(),
                     guide = guide_legend(reverse=TRUE)) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = "none")+
  geom_vline(aes(xintercept = 0.0),
             linetype = "dotted") +
  scale_x_continuous(labels = number_format(accuracy = 0.1)) + 
  labs(x = "26-Year Risk Difference (%)",
       y = "")

ggsave("gformula_figure_RD.pdf",
       dpi = 300,
       units = "cm",
       width = 25,
       height = 15)


