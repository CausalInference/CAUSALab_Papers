############################################################################

# Title: Recommendation-based Physical Activity and Dietary Interventions
#       for Adults Diagnosed with Breast or Prostate Cancer  

############################################################################

# Programmer: Emma McGee

# Date: July 1, 2024

# Purpose of Program: Create figure 3 (modified target trial results)

# Statistical Analyses: 
#  None

############################################################################

# Load packages
library("readxl")
library("ggplot2")
library("grid")
library("RColorBrewer")

# Check color names
f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
(cols <- f("RdYlBu"))

# Set working directory
main.dir = ""
setwd(main.dir)

############################################################################
######                   JOINT LIFESTYLE STRATEGY                     ######
############################################################################

### BREAST CANCER ###

df_b <- read_excel("figure_3_breast_tt1.xlsx")
df_b$RD <-as.numeric(df_b$RD)
df_b$RD_LL <-as.numeric(df_b$RD_LL)
df_b$RD_UL <-as.numeric(df_b$RD_UL)

### PROSTATE CANCER ###

df_p <- read_excel("figure_3_prostate_tt1.xlsx")
df_p$RD <-as.numeric(df_p$RD)
df_p$RD_LL <-as.numeric(df_p$RD_LL)
df_p$RD_UL <-as.numeric(df_p$RD_UL)

### COMBINE BREAST AND PROSTATE ###

df_all <- rbind(df_b, df_p)

### CREATE PLOT ###
plot_p_rd <- ggplot(data=df_all, aes(y=Analysis, x=RD, xmin=RD_LL, xmax=RD_UL, color=Type)) +
  geom_point(size=1.8, aes(shape = Type)) + 
  geom_errorbarh(height=.17, linewidth=0.8) +
  scale_x_continuous(breaks=c(-25, -20, -15, -10,  -5, 0, 5, 10), 
                     limits=c(-26, 6)) +
  scale_y_reverse(breaks=c(df_all$Analysis), labels=df_all$Name) +
  scale_color_manual(breaks=c('Initial specification', 
                              'Modifications to eligibility criteria', 
                              'Modifications to strategies', 
                              'Modifications to both eligibility criteria and strategies',
                              'Modifications to outcome',
                              'Modifications to adjustment covariates'),
                    values=c("#8C510A","#D73027","#FDAE61","#ABDDA4","#003C30","#35978F")) +
  scale_shape_manual(breaks=c('Initial specification', 
                              'Modifications to eligibility criteria', 
                              'Modifications to strategies', 
                              'Modifications to both eligibility criteria and strategies',
                              'Modifications to outcome',
                              'Modifications to adjustment covariates'),
                     values=c(0,17,16,1,18,15)) +
  labs(title='', x='20-year risk difference (%)', y = '') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=0.4) +
  facet_grid(~ Cancer) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        legend.title = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "none",
        legend.text = element_text(size = 12))
plot_p_rd
warnings()

ggsave("/figure_3_tt1.png", width = 40, height = 20, units = "cm")



############################################################################
######                       ALCOHOL STRATEGY                         ######
############################################################################

### BREAST CANCER ###

df_b <- read_excel("figure_3_breast_tt2.xlsx")
df_b$RD <-as.numeric(df_b$RD)
df_b$RD_LL <-as.numeric(df_b$RD_LL)
df_b$RD_UL <-as.numeric(df_b$RD_UL)

### PROSTATE CANCER ###

df_p <- read_excel("figure_3_prostate_tt2.xlsx")
df_p$RD <-as.numeric(df_p$RD)
df_p$RD_LL <-as.numeric(df_p$RD_LL)
df_p$RD_UL <-as.numeric(df_p$RD_UL)

### COMBINE BREAST AND PROSTATE ###

df_all <- rbind(df_b, df_p)

#### CREATE PLOT ###
plot_p_rd <- ggplot(data=df_all, aes(y=Analysis, x=RD, xmin=RD_LL, xmax=RD_UL, color=Type)) +
  geom_point(size=1.8, aes(shape = Type)) + 
  geom_errorbarh(height=.17, linewidth=0.8) +
  scale_x_continuous(breaks=c(-25, -20, -15, -10,  -5, 0, 5, 10), 
                     limits=c(-6.5, 11)) +
  scale_y_reverse(breaks=c(df_all$Analysis), labels=df_all$Name) +
  scale_color_manual(breaks=c('Initial specification', 
                              'Modifications to eligibility criteria', 
                              'Modifications to strategies', 
                              'Modifications to both eligibility criteria and strategies',
                              'Modifications to outcome',
                              'Modifications to adjustment covariates'),
                     values=c("#8C510A","#D73027","#FDAE61","#ABDDA4","#003C30","#35978F")) +
  scale_shape_manual(breaks=c('Initial specification', 
                              'Modifications to eligibility criteria', 
                              'Modifications to strategies', 
                              'Modifications to both eligibility criteria and strategies',
                              'Modifications to outcome',
                              'Modifications to adjustment covariates'),
                     values=c(0,17,16,1,18,15)) +
  labs(title='', x='20-year risk difference (%)', y = '') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=0.4) +
  facet_grid(~ Cancer) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        legend.title = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "none",
        legend.text = element_text(size = 12))
plot_p_rd
warnings()

ggsave("figure_3_tt2.png", width = 40, height = 20, units = "cm")


############################################################################
######                            LEGEND                              ######
############################################################################


#### CREATE PLOT ###
plot_p_rd <- ggplot(data=df_all, aes(y=Analysis, x=RD, xmin=RD_LL, xmax=RD_UL, color=Type)) +
  geom_point(size=1.8, aes(shape = Type)) + 
  geom_errorbarh(height=.17, linewidth=0.8) +
  scale_x_continuous(breaks=c(-25, -20, -15, -10,  -5, 0, 5, 10), 
                     limits=c(-6, 11)) +
  scale_y_reverse(breaks=c(df_all$Analysis), labels=df_all$Name) +
  scale_color_manual(breaks=c('Initial specification', 
                              'Modifications to eligibility criteria', 
                              'Modifications to strategies', 
                              'Modifications to both eligibility criteria and strategies',
                              'Modifications to outcome',
                              'Modifications to adjustment covariates'),
                     values=c("#8C510A","#D73027","#FDAE61","#ABDDA4","#003C30","#35978F")) +
  scale_shape_manual(breaks=c('Initial specification', 
                              'Modifications to eligibility criteria', 
                              'Modifications to strategies', 
                              'Modifications to both eligibility criteria and strategies',
                              'Modifications to outcome',
                              'Modifications to adjustment covariates'),
                     values=c(0,17,16,1,18,15)) +
  labs(title='', x='20-year risk difference (%)', y = '') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=0.4) +
  facet_grid(~ Cancer) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        legend.title = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "bottom",
        legend.text = element_text(size = 12))
plot_p_rd
warnings()

ggsave("figure_3_legend.png", width = 55, height = 20, units = "cm")
