## Computing bounds on individual counterfactual risk curves and bounds on ATE in STAR*D trial

# load necessary libraries
library(tidyverse)
library(readxl)
library(stringr)
library(gridExtra)
library(viridis)
library(grid)

##load data
stardkm <- read.csv("stard_data.csv")

##rename strata
stardkm$strata <- ifelse(stardkm$strata=="as.factor(strategy)=1", 1, ifelse(
  stardkm$strata=="as.factor(strategy)=2", 2, 3
))

##eliminate unecessary columns
stardkm <- stardkm %>% select(strata, time, n.risk, n.event)


## separate each strategy arm into separate dataset
strat1 <- stardkm %>% filter(strata == 1)
strat2 <- stardkm %>% filter(strata == 2)
strat3 <- stardkm %>% filter(strata == 3)

## add zeros to reflect that risks are estimated at end of interval
strat1$time <- strat1$time + 0.1
zero1 <- as.data.frame(t(c(1,0.0, 698, 0)))
colnames(zero1) <- c("strata", "time", "n.risk", "n.event")
strat1 <- as.data.frame(rbind(zero1, strat1))

strat2$time <- strat2$time + 0.1
zero2 <- as.data.frame(t(c(1,0.0, 492, 0)))
colnames(zero2) <- c("strata", "time", "n.risk", "n.event")
strat2 <- as.data.frame(rbind(zero2, strat2))

strat3$time <- strat3$time + 0.1
zero3 <- as.data.frame(t(c(1,0.0, 525,0)))
colnames(zero3) <- c("strata", "time", "n.risk", "n.event")
strat3 <- as.data.frame(rbind(zero3, strat3))



## create function that indexes number adherent (number at risk at beginning + 
## number events that already occurred in arm)

arm_bounds <- function(dat, events, num_at_risk, n, time){

# create vector indexing number of events in prev 
nevents <- c(0,events)
nevents <- nevents[0:length(events)]


# getting cumulative incidence at each time
nevents.tot <- cumsum(nevents)
events.tot <- cumsum(events)

# adding number at risk (adherent and no event) to number events at each time point
dat$n.comp <- num_at_risk + nevents.tot

# proportion of cumulative events in population adherent
dat$prop.event <- events.tot/dat$n.comp

# get proportion non-adherent out of 971 at each time point (971-proportion adherent (adherent or had event))
dat$n.nonad <- n - dat$n.comp
dat$prop.nonad <- dat$n.nonad/n

# get proportion adherent out of 971 at each time point
dat$prop.ad <- dat$n.comp/n

# calculate lower bound on counterfactual outcome at each time
dat$lb <- dat$prop.event*dat$prop.ad 

#calculate upper bound on counterfactual outcome
dat$ub <- dat$prop.event*dat$prop.ad + dat$prop.nonad

bds <- cbind(time, dat$lb, dat$ub)
return(bds)
}

strat1bds <- as.data.frame(arm_bounds(strat1, strat1$n.event, strat1$n.risk, 971, strat1$time))
colnames(strat1bds) <- c("time", "lb", "ub")
strat1bds$arm <- 1
strat2bds <- as.data.frame(arm_bounds(strat2, strat2$n.event, strat2$n.risk, 971, strat2$time))
colnames(strat2bds) <- c("time", "lb", "ub")
strat2bds$arm <- 2
strat3bds <- as.data.frame(arm_bounds(strat3, strat3$n.event, strat3$n.risk, 971, strat3$time))
colnames(strat3bds) <- c("time", "lb", "ub")
strat3bds$arm <- 3

## create graph of counterfactual outcome bounds over time
# combine data
stratbds <- rbind(strat1bds, strat2bds) #, strat3bds)
stratbds$arm <- as.factor(stratbds$arm)


# create graph
cuminc_plot <- ggplot(data = stratbds) + 
  geom_ribbon(aes(x = time, ymin = lb, ymax = ub, fill = arm, color=arm), alpha = 0.4) + 
  xlab("Month") + ylab(expression(E(Y^a))) + 
  labs(title = "B") + 
  theme(plot.caption = element_text(size = 5, hjust=0), legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), plot.title = element_text(size = 5), 
        axis.title = element_text(size = 5), axis.text = element_text(size=12), 
        legend.key.height = unit(0.1, 'cm'), legend.key.width = unit(0.09, 'cm'), 
        legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0), 
        legend.box.spacing = unit(2, "pt"), plot.margin = margin(0,0,0,0)) +
  scale_fill_viridis(name = "Treatment Arm",
                   labels = c('Sequential monotherapy','Sequential dual therapy', 'Guidelines therapy'),discrete = TRUE) + 
  scale_colour_viridis(name = "Treatment Arm",
                    labels = c('Sequential monotherapy','Sequential dual therapy', 'Guidelines therapy'), discrete = TRUE) + 
  scale_x_continuous(breaks=seq(0,9, by = 1))
  



#get bounds on ACE comparing strat 1 and 2 at 0 months, 3 months, 6 months, 9 months
#strat1
bds_arm1 <- stratbds %>% filter(arm==1)
bds_arm1 <- bds_arm1 %>% filter(time == 0.0 | time == 1.0 | time > 3.0 & time < 3.2 |
                                  time > 5.8 & time < 6.0 | time > 8.7)
bds_arm1$time <- round(bds_arm1$time)
bds_arm1$lb1 <- bds_arm1$lb
bds_arm1$ub1 <- bds_arm1$ub
bds_arm1 <- bds_arm1 %>% select(time, lb1, ub1)

# strat2
bds_arm2 <- stratbds %>% filter(arm==2)
bds_arm2 <- bds_arm2 %>% filter(time == 0.0 | time == 1.0 | time > 3.0 & time < 3.2 |
                                  time > 5.8 & time < 6.0 | time > 8.7)
bds_arm2$time <- round(bds_arm2$time)
bds_arm2$lb2 <- bds_arm2$lb
bds_arm2$ub2 <- bds_arm2$ub
bds_arm2 <- bds_arm2 %>% select(time, lb2, ub2)

# combine and get ACE bounds
ace_bds <- merge(bds_arm1, bds_arm2, by = "time")
ace_bds$ace_lb <- ace_bds$lb1 - ace_bds$ub2
ace_bds$ace_ub <- ace_bds$ub1 - ace_bds$lb2
ace_bds <- ace_bds %>% select(time, ace_lb, ace_ub)

# graph ACE bounds
ace_plot <- ggplot(data = ace_bds) + geom_errorbar(aes(x = time, ymin = ace_lb, ymax = ace_ub)) + 
  labs(title = "A") +
  xlab("Month") + ylab("Risk Difference") + scale_x_continuous(breaks=c(0,3,6,9)) +
  theme(plot.caption = element_text(size = 5, hjust=0), legend.text = element_text(size = 3), 
        legend.title = element_text(size = 5), plot.title = element_text(size = 5), 
        axis.title = element_text(size = 12), axis.text = element_text(size=12)) +
  scale_y_continuous(breaks = c(-1, -0.8, -.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), 
                     limits = c(-1,1)) + 
  scale_x_continuous(breaks=seq(0,9, by = 1))


##generate overall title and caption
tg <- textGrob(label = "Figure. Assumption-free bounds on the average causal effect of sequential monotherapy vs sequential dual non-SSRI therapy, and on
               counterfactual risk of remission under different treatment regimes", gp = gpar(fontsize = 5))
cg <- textGrob(label = "Figure Panel A shows bounds on the average causal effect of sequential monotherapy vs. sequential dual non-SSRI therapy on risk of remission at 0, 1, 3, 6, and
9 months of followup among participants in the STAR*D trial who were initially assigned to citalopram and experienced treatment failure. Figure Panel B shows
bounds on the counterfactual risk of remission under all three regimes of interest (sequential monotherapy, sequential dual non-SSRI therapy, sequential
guidelines-based therapy) through followup in the same participants. All data were drawn from a secondary analysis of the STAR*D trial conducted by
Szmulewicz et al.", , gp = gpar(fontsize = 5))

plotlist <- list(ace_plot, cuminc_plot)
g <- c(list(tg), plotlist)
#version with caption as part of it
#g <- c(list(tg), plotlist, list(cg))
N = length(plotlist)
laym <- rbind(rep(1, N), rep(2:3), rep(4, N))
gridimage <- grid.arrange(grobs = g, layout_matrix = laym, heights = c(1,9.5,2.5))


  









