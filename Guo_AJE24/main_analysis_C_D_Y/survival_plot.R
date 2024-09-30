################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_survival_plot.R
# Programer: Fuyu Guo
# Last updated: 2023-Jul-1


# 1) Purpose: this R program read in g-formula results and made survival curves

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
.libPaths("/n/home00/fyguo/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(ggsci)
library(gfoRmula, lib.loc = "/n/home00/fyguo/R/x86_64-pc-linux-gnu-library/4.1")
library(data.table)
library(scales)
setwd("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y")

# total pca
fit_tot <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_tot_1115.rds")
df_tot <- fit_tot$result[`Interv.` %in% c(0, 1, 9, 10), .(`Interv.`, `k`,`g-form risk`, `Risk lower 95% CI`, `Risk upper 95% CI`)]
colnames(df_tot) <- c("Interventions",
                      "Time period",
                      "Risk",
                      "Risk_low_CI",
                      "Risk_up_CI")
df_tot[,`:=`(
  Survival = 1-Risk,
  Survival_low_CI = 1-Risk_up_CI,
  Survival_up_CI = 1-Risk_low_CI
)]

df_tot$outcome <- "Total Prostate Cancer"


# advanced pca
fit_adv <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_adv_1115.rds")
df_adv <- fit_adv$result[`Interv.` %in% c(0, 1, 9, 10), .(`Interv.`, `k`,`g-form risk`, `Risk lower 95% CI`, `Risk upper 95% CI`)]
colnames(df_adv) <- c("Interventions",
                      "Time period",
                      "Risk",
                      "Risk_low_CI",
                      "Risk_up_CI")
df_adv[,`:=`(
  Survival = 1-Risk,
  Survival_low_CI = 1-Risk_up_CI,
  Survival_up_CI = 1-Risk_low_CI
)]


df_adv$outcome <- "Advanced Prostate Cancer"

# lethal pca
fit_lethal <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_lethal_1115.rds")
df_lethal <- fit_lethal$result[`Interv.` %in% c(0, 1, 9, 10), .(`Interv.`, `k`,`g-form risk`, `Risk lower 95% CI`, `Risk upper 95% CI`)]
colnames(df_lethal) <- c("Interventions",
                         "Time period",
                         "Risk",
                         "Risk_low_CI",
                         "Risk_up_CI")
df_lethal[,`:=`(
  Survival = 1-Risk,
  Survival_low_CI = 1-Risk_up_CI,
  Survival_up_CI = 1-Risk_low_CI
)]


df_lethal$outcome <- "Lethal Prostate Cancer"

# fatal pca
fit_fatal <- read_rds("/n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_fatal_1115.rds")
df_fatal <- fit_fatal$result[`Interv.` %in% c(0, 1, 9, 10), .(`Interv.`, `k`,`g-form risk`, `Risk lower 95% CI`, `Risk upper 95% CI`)]
colnames(df_fatal) <- c("Interventions",
                        "Time period",
                        "Risk",
                        "Risk_low_CI",
                        "Risk_up_CI")
df_fatal[,`:=`(
  Survival = 1-Risk,
  Survival_low_CI = 1-Risk_up_CI,
  Survival_up_CI = 1-Risk_low_CI
)]

df_fatal$outcome <- "Fatal Prostate Cancer"

############################################
# combine all the dataframes
df <- rbind(df_tot, df_adv, df_lethal, df_fatal)

df <- df %>%
  mutate(Interventions = case_when(
    Interventions == 0 ~ "No intervention (Natural course)",
    Interventions == 1 ~ "Intervention on physical activity",
    Interventions == 9 ~ "Intervention on diet",
    Interventions == 10 ~ "Joint intervention"
  ))


df$Interventions <- factor(df$Interventions,
                           levels = c(
                             "No intervention (Natural course)",
                             "Intervention on physical activity",
                             "Intervention on diet",
                             "Joint intervention"
                           ))


df$outcome <- factor(df$outcome,
                     levels = c("Total Prostate Cancer",
                                "Advanced Prostate Cancer",
                                "Lethal Prostate Cancer",
                                "Fatal Prostate Cancer"))

# transfer the time period to actual year
df$`Time period` <- (df$`Time period`+1)*2 




# at the beginning of 1990 all the Survival probability should be 1
tmp <- df[1:16,]
tmp <- tmp %>% mutate(`Time period` = 0,
                      outcome = rep(c("Total Prostate Cancer",
                                      "Advanced Prostate Cancer",
                                      "Lethal Prostate Cancer",
                                      "Fatal Prostate Cancer"),
                                    each = 4),
                      Risk = 0,
                      Risk_low_CI = 0,
                      Risk_up_CI = 0,
                      Survival = 1,
                      Survival_low_CI = 1,
                      Survival_up_CI = 1)
df <- rbind(tmp, df)

library(ggsci)
plot_tot <- ggplot(df[outcome == "Total Prostate Cancer"]) + 
  geom_line(aes(x = `Time period`,
                y = Survival*100,
                color = Interventions,
                group = Interventions)) + 
  geom_point(aes(x = `Time period`,
                 y = Survival*100,
                 shape = Interventions,
                 color = Interventions,
                 group = Interventions)) +
  scale_color_jama()+
  scale_fill_jama()+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_x_continuous(breaks = seq(0, 26.5, 2))+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())+
  scale_y_continuous(labels = number_format(accuracy = 0.1)) +
  labs(title = "Total Prostate Cancer",
       x = "Years",
       y = "Survival (%)")


plot_adv <- ggplot(df[outcome == "Advanced Prostate Cancer"]) + 
  geom_line(aes(x = `Time period`,
                y = Survival*100,
                color = Interventions,
                group = Interventions)) + 
  geom_point(aes(x = `Time period`,
                 y = Survival*100,
                 shape = Interventions,
                 color = Interventions,
                 group = Interventions)) +
  scale_color_jama()+
  scale_fill_jama()+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_x_continuous(breaks = seq(0, 26.5, 2))+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())+
  labs(title = "Advanced Prostate Cancer",
       x = "Years",
       y = NULL)+
  coord_cartesian(ylim = c(97.5, 100)) 


plot_lethal <- ggplot(df[outcome == "Lethal Prostate Cancer"]) + 
  geom_line(aes(x = `Time period`,
                y = Survival*100,
                color = Interventions,
                group = Interventions)) + 
  geom_point(aes(x = `Time period`,
                 y = Survival*100,
                 shape = Interventions,
                 color = Interventions,
                 group = Interventions)) +
  scale_color_jama()+
  scale_fill_jama()+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_x_continuous(breaks = seq(0, 26.5, 2))+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Lethal Prostate Cancer",
       x = "Years",
       y = "Survival (%)")+
  coord_cartesian(ylim = c(97.5, 100))


plot_fatal <- ggplot(df[outcome == "Fatal Prostate Cancer"]) + 
  geom_line(aes(x = `Time period`,
                y = Survival*100,
                color = Interventions,
                group = Interventions)) + 
  geom_point(aes(x = `Time period`,
                 y = Survival*100,
                 shape = Interventions,
                 color = Interventions,
                 group = Interventions)) +
  scale_color_jama()+
  scale_fill_jama()+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_x_continuous(breaks = seq(0, 26.5, 2))+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Fatal Prostate Cancer",
       x = "Years",
       y = NULL)+
  coord_cartesian(ylim = c(97.5, 100))


########################################
# combine all the plots together
library(ggpubr)
png("gformula_figure_survival.png",
    width = 20,
    height = 14,
    units = "cm",
    res = 300)
ggarrange(plot_tot,plot_adv,plot_lethal,plot_fatal,
          nrow = 2,
          ncol = 2,
          common.legend = T,
          legend = "bottom")

dev.off()
Sys.time()

