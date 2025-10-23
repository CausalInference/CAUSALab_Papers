############################################################################

# Title: Recommendation-based Physical Activity and Dietary Interventions
#       for Adults Diagnosed with Breast or Prostate Cancer  

############################################################################

# Programmer: Emma McGee 

# Date: July 1, 2024

# Purpose of Program: Create figure 2 (adjusted survival curves for initial specification)

# Statistical Analyses: 
#  None

############################################################################

# Load required packages

library(haven)
library(plyr)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(tibble)
library(tidyverse)
library(grid)
library(gridExtra)

############################################################################
######                   JOINT LIFESTYLE STRATEGY                     ######
############################################################################

###############################
###### PROSTATE CANCER ########
###############################

forgraphs <- read_sas("forgraphs_prostate.sas7bdat")
main.dir = ""
setwd(main.dir)

head(forgraphs)
forgraphs_final <- forgraphs

# Check that risks for each intervention are correct
surv<- forgraphs_final[11,]
risk <-1-surv
risk

# Take columns for joint intervention and natural course
frame1<-forgraphs_final[c(2,5)]

# Rename columns
frame2<-rename(frame1, c("No intervention" = "surv0",
                         "Adhere to 7 physical activity and dietary recommendations"= "surv1"))

### Create Figure ###

# Prep data for plotting
time<-c(0,2,4,6,8,10,12,14,16,18,20)
frame3min2 <- cbind(frame2, time)
frame3min3 <- as_tibble(frame3min2)
long <- gather(frame3min3, Intervention, survprob, `No intervention`:`Adhere to 7 physical activity and dietary recommendations`, factor_key=TRUE)
gline = linesGrob(y = c(0.2, 1),x = c(-.01, .01),  gp = gpar(col = "black", lwd = 1.4))

# Plot
plot <- ggplot(data=long, aes(x=factor(time), y=survprob, group=Intervention, colour=Intervention))+
  geom_line(linewidth=0.7) +
  geom_point(aes(shape = Intervention), size=2) +
  scale_y_continuous(breaks=c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits=c(0.40, 1),
                     labels=c("0%", "50%","60%","70%","80%","90%","100%")) +
  scale_x_discrete(expand=c(0.05,0.05)) +
  ylab("Survival probability") +
  xlab("Year of follow-up") +
  scale_color_manual(labels = c("No intervention",
                                "Jointly adhere to 7 physical activity and dietary recommendations"),
                     values = c("#4e82b4", "#922b3e")) +
  scale_shape_manual(labels = c("No intervention",
                                "Jointly adhere to 7 physical activity and dietary recommendations"),
                     values = c(16, 17)) +
  theme_minimal()+
  theme(axis.text = element_text(size=11, colour = "black"), 
        axis.line.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=12, colour = "black"),
        axis.ticks.x = element_line(linewidth = 0.8), 
        axis.ticks.y = element_line(linewidth = 0.8),
        axis.ticks.length = unit(5, "pt"))+
  annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  annotate(geom = "segment", x = -Inf, xend = -Inf, y =  .461, yend = .441,
           linetype = "solid", color = "white", lwd=2)+
  annotation_custom(gline, ymin=.45, ymax=.47, xmin=-Inf, xmax=Inf) +
  annotation_custom(gline, ymin=.43, ymax=.45, xmin=-Inf, xmax=Inf)+
  font("xlab",size=11)+
  font("ylab",size=11)+
  font("legend.text",size=9.5)
plot

ggsave(filename="figure_2_prostate.png", device = "png", width = 5, height = 2.8, units = "in", plot = plot)  


###############################
######  BREAST CANCER  ########
###############################

forgraphs <- read_sas("forgraphs_breast.sas7bdat")
main.dir = ""
setwd(main.dir)

head(forgraphs)
forgraphs_final <- forgraphs

# Check that risks for each intervention are correct
surv<- forgraphs_final[11,]
risk <-1-surv
risk

# Take columns for joint intervention and natural course
frame1<-forgraphs_final[c(2,5)]

# Rename columns
frame2<-rename(frame1, c("No intervention" = "surv0",
                         "Adhere to 7 physical activity and dietary recommendations"= "surv1"))

### Create Figure ###

# Prep data for plotting
frame3min2 <- cbind(frame2, time)
frame3min3 <- as_tibble(frame3min2)
long <- gather(frame3min3, Intervention, survprob, `No intervention`:`Adhere to 7 physical activity and dietary recommendations`, factor_key=TRUE)

# Plot 
plot <- ggplot(data=long, aes(x=factor(time), y=survprob, group=Intervention, colour=Intervention))+
  geom_line(linewidth=0.7) +
  geom_point(aes(shape = Intervention), size=2) +
  scale_y_continuous(breaks=c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits=c(0.40, 1),
                     labels=c("0%", "50%","60%","70%","80%","90%","100%")) +
  scale_x_discrete(expand=c(0.05,0.05)) +
  ylab("Survival probability") +
  xlab("Year of follow-up") +
  scale_color_manual(labels = c("No intervention",
                                "Jointly adhere to 7 physical activity and dietary recommendations"),
                     values = c("#4e82b4", "#922b3e")) +
  scale_shape_manual(labels = c("No intervention",
                                "Jointly adhere to 7 physical activity and dietary recommendations"),
                     values = c(16, 17)) +
  theme_minimal()+
  theme(axis.text = element_text(size=11, colour = "black"), 
        axis.line.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=12, colour = "black"),
        axis.ticks.x = element_line(linewidth = 0.8), 
        axis.ticks.y = element_line(linewidth = 0.8),
        axis.ticks.length = unit(5, "pt"))+
  annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  annotate(geom = "segment", x = -Inf, xend = -Inf, y =  .461, yend = .441,
           linetype = "solid", color = "white", lwd=2)+
  annotation_custom(gline, ymin=.45, ymax=.47, xmin=-Inf, xmax=Inf) +
  annotation_custom(gline, ymin=.43, ymax=.45, xmin=-Inf, xmax=Inf)+
  font("xlab",size=11)+
  font("ylab",size=11)+
  font("legend.text",size=9.5)
plot

ggsave(filename="figure_2_breast.png", device = "png", width = 5, height = 2.8, units = "in", plot = plot)  


############################################################################
######                       ALCOHOL STRATEGY                         ######
############################################################################

###############################
###### PROSTATE CANCER ########
###############################

forgraphs <- read_sas("forgraphs_prostate.sas7bdat")
main.dir = ""
setwd(main.dir)

head(forgraphs)
forgraphs_final <- forgraphs

# Check that risks for each intervention are correct
risk <-1-surv
risk

# Take columns for alcohol intervention and natural course
frame1<-forgraphs_final[c(2,8)]

# Rename columns
frame2<-rename(frame1, c("No intervention" = "surv0" ,
                         "Do not consume alcohol" ="surv10" ))

### Create Figure ###

# Prep data for plotting
frame3min2 <- cbind(frame2, time)
frame3min3 <- as_tibble(frame3min2)
long <- gather(frame3min3, Intervention, survprob, `No intervention`:`Do not consume alcohol`, factor_key=TRUE)
gline = linesGrob(y = c(0.2, 1),x = c(-.01, .01),  gp = gpar(col = "black", lwd = 1.4))

# Plot
plot <- ggplot(data=long, aes(x=factor(time), y=survprob, group=Intervention, colour=Intervention))+
  geom_line(linewidth=0.7) +
  geom_point(aes(shape = Intervention), size=2) +
  scale_y_continuous(breaks=c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits=c(0.40, 1),
                     labels=c("0%", "50%","60%","70%","80%","90%","100%")) +
  scale_x_discrete(expand=c(0.05,0.05)) +
  ylab("Survival probability") +
  xlab("Year of follow-up") +
  scale_color_manual(labels = c("No intervention",
                                "Do not consume alcohol"),
                     values = c("#4e82b4", "Goldenrod 1")) +
  scale_shape_manual(labels = c("No intervention",
                                "Do not consume alcohol"),
                     values = c(16, 17)) +
  theme_minimal()+
  theme(axis.text = element_text(size=11, colour = "black"), 
        axis.line.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=12, colour = "black"),
        axis.ticks.x = element_line(linewidth = 0.8), 
        axis.ticks.y = element_line(linewidth = 0.8),
        axis.ticks.length = unit(5, "pt"))+
  annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  annotate(geom = "segment", x = -Inf, xend = -Inf, y =  .461, yend = .441,
           linetype = "solid", color = "white", lwd=2)+
  annotation_custom(gline, ymin=.45, ymax=.47, xmin=-Inf, xmax=Inf) +
  annotation_custom(gline, ymin=.43, ymax=.45, xmin=-Inf, xmax=Inf)+
  font("xlab",size=11)+
  font("ylab",size=11)+
  font("legend.text",size=9.5)
plot

 
ggsave(filename="figure_2_prostate_alc.png", device = "png", width = 5, height = 2.8, units = "in", plot = plot)  



###############################
######  BREAST CANCER  ########
###############################

forgraphs <- read_sas("forgraphs_breast.sas7bdat")
main.dir = ""
setwd(main.dir)

head(forgraphs)
forgraphs_final <- forgraphs

# Check that risks for each intervention are correct
surv<- forgraphs_final[11,]
risk <-1-surv
risk

# Take columns for joint intervention and natural course
frame1<-forgraphs_final[c(2,8)]

# Rename columns
frame2<-rename(frame1, c("No intervention" = "surv0",
                         "Do not consume alcohol" = "surv10"))

### Create Figure ###

# Prep data for plotting
frame3min2 <- cbind(frame2, time)
frame3min3 <- as_tibble(frame3min2)
long <- gather(frame3min3, Intervention, survprob, `No intervention`:`Do not consume alcohol`, factor_key=TRUE)

# Plot
plot <- ggplot(data=long, aes(x=factor(time), y=survprob, group=Intervention, colour=Intervention))+
  geom_line(linewidth=0.7) +
  geom_point(aes(shape = Intervention), size=2) +
  scale_y_continuous(breaks=c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits=c(0.40, 1),
                     labels=c("0%", "50%","60%","70%","80%","90%","100%")) +
  scale_x_discrete(expand=c(0.05,0.05)) +
  ylab("Survival probability") +
  xlab("Year of follow-up") +
  scale_color_manual(labels = c("No intervention",
                                "Do not consume alcohol"),
                     values = c("#4e82b4", "Goldenrod 1")) +
  scale_shape_manual(labels = c("No intervention",
                                "Do not consume alcohol"),
                     values = c(16, 17)) +
  theme_minimal()+
  theme(axis.text = element_text(size=11, colour = "black"), 
        axis.line.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=12, colour = "black"),
        axis.ticks.x = element_line(linewidth = 0.8), 
        axis.ticks.y = element_line(linewidth = 0.8),
        axis.ticks.length = unit(5, "pt"))+
  annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  annotate(geom = "segment", x = -Inf, xend = -Inf, y =  .461, yend = .441,
           linetype = "solid", color = "white", lwd=2)+
  annotation_custom(gline, ymin=.45, ymax=.47, xmin=-Inf, xmax=Inf) +
  annotation_custom(gline, ymin=.43, ymax=.45, xmin=-Inf, xmax=Inf)+
  font("xlab",size=11)+
  font("ylab",size=11)+
  font("legend.text",size=9.5)
plot

ggsave(filename="figure_2_breast_alc.png", device = "png", width = 5, height = 2.8, units = "in", plot = plot)  



############################################################################
######                            LEGEND                              ######
############################################################################


# Prep data for plotting
frame1<-forgraphs_final[c(2,5,8)]
frame2<-rename(frame1, c("No intervention" = "surv0",
                         "Adhere to 7 physical activity and dietary recommendations" = "surv1",
                         "Do not consume alcohol" = "surv10"))
frame3min2 <- cbind(frame2, time)
frame3min3 <- as_tibble(frame3min2)
long <- gather(frame3min3, Intervention, survprob, `No intervention`:`Do not consume alcohol`, factor_key=TRUE)

# Plot with legend
plot <- ggplot(data=long, aes(x=factor(time), y=survprob, group=Intervention, colour=Intervention))+
  geom_line(linewidth=0.7) +
  geom_point(aes(shape = Intervention), size=2) +
  scale_y_continuous(breaks=c(0.5, 0.6, 0.7, 0.8, 0.9, 1), limits=c(0.45, 1)) +
  scale_x_discrete(expand=c(0.05,0.05)) +
  ylab("Survival probability") +
  xlab("Year of follow-up") +
  scale_color_manual(labels = c("No intervention",
                                "Jointly adhere to 7 physical activity and dietary recommendations",
                                "Do not consume alcohol"),
                     values = c("#4e82b4", "#922b3e", "Goldenrod 1")) +
  scale_shape_manual(labels = c("No intervention",
                                "Jointly adhere to 7 physical activity and dietary recommendations",
                                "Do not consume alcohol"),
                     values = c(16, 17, 17)) +
  theme_minimal()+
  theme(axis.text = element_text(size=11, colour = "black"), 
        axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12, colour = "black"),
        axis.ticks.x = element_line(linewidth = 0.8), 
        axis.ticks.y = element_line(linewidth = 0.8),
        axis.ticks.length = unit(5, "pt"))+
  font("xlab",size=11)+
  font("ylab",size=11)+
  font("legend.text",size=9.5)
plot

ggsave(filename="figure_2_legend.png", device = "png", width = 10, height = 5, units = "in", plot = plot)  
