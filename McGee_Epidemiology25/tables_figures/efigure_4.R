############################################################################

# Title: Recommendation-based Physical Activity and Dietary Interventions
#       for Adults Diagnosed with Breast or Prostate Cancer  

############################################################################

# Programmer: Emma McGee

# Date: July 1, 2024

# Purpose of Program: Create efigure 4 (means of BMI and fiber over time)

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
df_p <- read_excel("efigure_4_prostate.xlsx")
df_p <- df_p[is.na(df_p$mean_bmi) == FALSE,]

df_p$mean_bmi <- as.numeric(df_p$mean_bmi)
df_p$mean_bmi_ul <- as.numeric(df_p$mean_bmi_ul)
df_p$mean_bmi_ll <- as.numeric(df_p$mean_bmi_ll)
df_p$mean_fiber <- as.numeric(df_p$mean_fiber)
df_p$mean_fiber_ul <- as.numeric(df_p$mean_fiber_ul)
df_p$mean_fiber_ll <- as.numeric(df_p$mean_fiber_ll)


### Breast ###
df_b <- read_excel("efigure_4_breast.xlsx")
df_b <- df_b[is.na(df_b$mean_bmi) == FALSE,]

df_b$mean_bmi <- as.numeric(df_b$mean_bmi)
df_b$mean_bmi_ul <- as.numeric(df_b$mean_bmi_ul)
df_b$mean_bmi_ll <- as.numeric(df_b$mean_bmi_ll)
df_b$mean_fiber <- as.numeric(df_b$mean_fiber)
df_b$mean_fiber_ul <- as.numeric(df_b$mean_fiber_ul)
df_b$mean_fiber_ll <- as.numeric(df_b$mean_fiber_ll)

### Combine data ###
df <- rbind(df_b, df_p)

####################
####     BMI    ####
####################

# Plot
plot_rd <- ggplot(data=df, aes(x=year, y=mean_bmi, ymin=mean_bmi_ll, ymax=mean_bmi_ul, color=strategy, group=strategy, shape=strategy, fill=strategy)) +
  geom_point(size=2) + 
  geom_line(linewidth=0.5) +
  geom_ribbon(alpha=0.2, linewidth = 0) +
  scale_y_continuous(breaks=c(22, 23, 24, 25, 26, 27, 28), limits=c(22, 28)) +
  scale_x_continuous(breaks=c(df$year)) +
  labs(title='', x='', y='') +
  geom_hline(yintercept=25, color='black', linetype='dashed', alpha=0.4, linewidth=0.8) +
  facet_grid(~ Cancer) +
  scale_color_manual(labels=c("Jointly adhere to 7 physical activity and dietary recommendations", "No intervention"),
                     values=c("#922b3e", "#4e82b4"),
                      breaks =c("Jointly adhere to 7 physical activity and dietary recommendations", "No intervention")) +
  scale_shape_manual(labels=c("Jointly adhere to 7 physical activity and dietary recommendations", "No intervention"),
                     values=c(17, 16)) +
  scale_fill_manual(labels=c("Jointly adhere to 7 physical activity and dietary recommendations", "No intervention"),
                    values=c( "#922b3e", "#4e82b4")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=6.5),
        axis.text.y=element_text(size=6.5),
        axis.title.x=element_text(size=6.5, margin=margin(t=10, r=0, b=0, l=0)),
        axis.title.y=element_text(size=6.5),
        strip.text.x=element_text(size=6.5),
        panel.spacing.x=unit(0.5, "lines"),
        panel.border=element_rect(color="black", fill=NA, linewidth=1),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=6.5)) +
  labs(x="Year of follow-up", y=bquote("Body mass index, kg/m"^2)) 
plot_rd

ggsave(filename=file.path(main.dir,"efigure_4_bmi.png"), device = "png", 
       width = 7.5, height = 3.1, units = "in", 
       plot = plot_rd) 


####################
####   FIBER    ####
####################

# Plot
plot_rd <- ggplot(data=df, aes(x=year, y=mean_fiber, ymin=mean_fiber_ll, ymax=mean_fiber_ul, color=strategy, group=strategy, shape=strategy, fill=strategy)) +
  geom_point(size=2, na.rm = TRUE) + 
  geom_line(linewidth=0.5, na.rm = TRUE, data=df[!is.na(df$mean_fiber),]) +
  geom_ribbon(alpha=0.2, linewidth = 0, na.rm = TRUE) +
  scale_y_continuous(breaks=c(18, 20, 22, 24, 26, 28, 30, 32, 34), limits=c(18, 34)) +
  scale_x_continuous(breaks=c(df$year)) +
  labs(title='', x='', y='') +
  geom_hline(yintercept=30, color='black', linetype='dashed', alpha=0.4, linewidth=0.8) +
  facet_grid(~ Cancer) +
  scale_color_manual(labels=c("Jointly adhere to 7 physical activity and dietary recommendations", "No intervention"),
                     values=c("#922b3e", "#4e82b4"),
                     breaks =c("Jointly adhere to 7 physical activity and dietary recommendations", "No intervention")) +
  scale_shape_manual(labels=c("Jointly adhere to 7 physical activity and dietary recommendations", "No intervention"),
                     values=c(17, 16)) +
  scale_fill_manual(labels=c("Jointly adhere to 7 physical activity and dietary recommendations", "No intervention"),
                    values=c( "#922b3e", "#4e82b4")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=6.5),
        axis.text.y=element_text(size=6.5),
        axis.title.x=element_text(size=6.5, margin=margin(t=10, r=0, b=0, l=0)),
        axis.title.y=element_text(size=6.5),
        strip.text.x=element_text(size=6.5),
        panel.spacing.x=unit(0.5, "lines"),
        panel.border=element_rect(color="black", fill=NA, linewidth=1),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=6.5)) +
  labs(x="Year of follow-up", y="Fiber, grams/day") 
plot_rd


ggsave(filename=file.path(main.dir,"efigure_4_fiber.png"), device = "png", 
       width = 7.5, height = 3.1, units = "in", 
       plot = plot_rd)
