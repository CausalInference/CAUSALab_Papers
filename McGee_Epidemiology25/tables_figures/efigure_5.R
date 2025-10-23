############################################################################

# Title: Recommendation-based Physical Activity and Dietary Interventions
#       for Adults Diagnosed with Breast or Prostate Cancer  

############################################################################

# Programmer: Emma McGee 

# Date: July 1, 2024

# Purpose of Program: Create efigure 5 (specific components of the joint intervention)

# Statistical Analyses: 
#  None

############################################################################

# Load packages
library("readxl")
library("ggplot2")
library("grid")
main.dir = ""
setwd(main.dir)


############################################################
#### TARGET TRIAL 1 - JOINT INTERVENTION
##########################################################

### Prostate ###
df_p <- read_excel("efigure_5_prostate.xlsx")
df_p <- df_p[is.na(df_p$Cancer) == FALSE,]

# Transform to numeric and format
df_p$RD <- as.numeric(df_p$RD)
df_p$RD_LL <- as.numeric(df_p$RD_LL)
df_p$RD_UL <- as.numeric(df_p$RD_UL)

### Breast ###
df_b <- read_excel("efigure_5_breast.xlsx")
df_b <- df_b[is.na(df_b$Cancer) == FALSE,]

# Transform to numeric and format
df_b$RD <- as.numeric(df_b$RD)
df_b$RD_LL <- as.numeric(df_b$RD_LL)
df_b$RD_UL <- as.numeric(df_b$RD_UL)

### Combine data ###
df <- rbind(df_b, df_p)



### Create plot ###
df_p$Analysis <- seq(1, length(df_p$Analysis)*1, 1)

plot_rd <- ggplot(data=df, aes(y=Analysis, x=RD, xmin=RD_LL, xmax=RD_UL)) +
  geom_point(size=2, shape = 15) + 
  geom_errorbarh(height=0.2, linewidth=0.8) +
  scale_x_continuous(breaks=c(-11, -10, -9 ,-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3),
                     limits=c(-11, 3)) +
  scale_y_reverse(breaks=c(df_p$Analysis), labels=df_p$Name) +
  labs(title='', x='', y = '') +
  #geom_vline(xintercept=4.8, color='pink', linetype='dashed') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=0.4, linewidth=0.8) +
  facet_grid(~ Cancer) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6.5),
        axis.text.y = element_text(size = 6.5),
        axis.title.x = element_text(size = 6.5, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 6.5),
        strip.text.x = element_text(size = 6.5),
        panel.spacing.x = unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 6.5)) +
  labs(x = "20-year risk difference") 
plot_rd
warnings()


# Create plot
ggsave(filename=file.path(main.dir,"efigure_5.png"), device = "png", 
       width = 7.5, height = 4.2, units = "in", 
       plot = plot_rd) 
