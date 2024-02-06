########################################################################################################################################
#   
# CAUSALab 
# 
# File...........: Estimated cumulative incidence of cardiovascular events in ART-naive and ART-experienced individuals [creation of figure]
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
#
########################################################################################################################################

if (!require("data.table")) install.packages("data.table")
library(data.table)
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!require("gridExtra")) install.packages("gridExtra")
library(gridExtra)


#### PLOT IN ART-NAIVE INDIVIDUALS
risk_naive_final <- fread(file="path to file")

plot_naive <- 
  ggplot(risk_naive_final,
         aes(x=month, y=risk)) +
  geom_line(aes(x=month, y=risk, color=group2, linetype=group2)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group2), alpha=0.1) +
  xlab("Months") +
  scale_x_continuous(breaks=seq(0, 48, by=4)) +
  ylab("Cumulative Incidence (%)") + 
  scale_y_continuous(limits=c(0, 0.012),
                     breaks=c(0, 0.0024, 0.0048, 0.0072, 0.0096, 0.012),
                     labels=c("0.0%", "0.24%", "0.48%",
                              "0.72%", "0.96%", "1.2%")) +
  theme_bw() +
  ggtitle("a.) ART-naive individuals") +
  theme(axis.text.y = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.15, 0.9), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key=element_blank(),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.2),
        plot.title = element_text(size = 11, face="bold")) +
  scale_linetype_manual(values=c("solid", "dotdash")) +
  scale_color_manual(values=c("indianred4", "steelblue4"))


#### PLOT IN ART-EXPERIENCED INDIVIDUALS
risk_exp_final <- fread(file="path to file")

plot_experienced <- 
  ggplot(risk_exp_final,
         aes(x=month, y=risk)) +
  geom_line(aes(x=month, y=risk, color=group2, linetype=group2)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=group2), alpha=0.1) +
  xlab("Months") +
  scale_x_continuous(breaks=seq(0, 48, by=4)) +
  ylab("Cumulative Incidence (%)") + 
  scale_y_continuous(limits=c(0, 0.025),
                     breaks=c(0, 0.005, 0.010, 0.015, 0.020, 0.025),
                     labels=c("0.0%", "0.5%", "1.0%",
                              "1.5%", "2.0%", "2.5%")) +
  theme_bw() +
  ggtitle("b.) ART-experienced individuals") +
  theme(axis.text.y = element_text(color="black", size=8),
        axis.text.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=9), 
        axis.title.x = element_text(color="black", size=9), 
        legend.position = c(0.20, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 10), 
        legend.key=element_blank(),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.2),
        plot.title = element_text(size = 11, face="bold")) +
  scale_linetype_manual(values=c("solid", "dotdash")) +
  scale_color_manual(values=c("indianred4", "steelblue4"))


# LAYOUT OF BOTH PLOTS TOGETHER
plot_final <- grid.arrange(plot_naive, plot_experienced, nrow=2)

