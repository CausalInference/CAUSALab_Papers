library(here)
library(ggplot2)
library(tidyverse)
library(glue)
library(ggpubr)
library(scales)

plotted_results <- list()

DATE <- "2023-11-27"

FILE_SUFFIX <- glue("{DATE}")

# load plots
outputs <- read_csv(here::here("outputs", glue("outputs_{FILE_SUFFIX}.csv")))

# function to plot figure
plot_by_analysis <- function(outputs, input_outcome_prevalence, input_outcome_treatment_rr, 
                             input_outcome_mediator_rr, input_correlation, 
                            ystub,  xlabel, ylabel, ylimits, exponentiate=FALSE, yticks=waiver(), exclude=c("none"),str_pad_val=25) {
  
  # prepare data
  plot_data <- outputs %>%
    filter(treatment_prevalence == 0.5, covariate_prevalence == 0.1, mediator_prevalence == 0.02, 
           mediator_covariate_rr == 1, outcome_prevalence == input_outcome_prevalence, 
           outcome_treatment_rr == input_outcome_treatment_rr, outcome_mediator_rr == input_outcome_mediator_rr,
           outcome_covariate_rr == 1, correlation == input_correlation) %>%
    select(mediator_treatment_rr, starts_with(ystub)) %>%
    pivot_longer(cols=starts_with(ystub), names_to="Analysis", names_prefix=ystub, values_to="result") %>% 
    filter(!(Analysis %in% exclude)) %>%
    mutate(Analysis=factor(Analysis)) %>% 
    mutate(Analysis=droplevels(recode_factor(Analysis, naive_analysis=str_pad("Infant-level: no GEE", str_pad_val, side="right"), 
                                  gee_exchangeable = str_pad("Infant-level: GEE exchangeable", str_pad_val, side="right"), gee_independence=str_pad("Infant-level: GEE independence", str_pad_val, side="right"),
                                  pregnancy_level=str_pad("Pregnancy-level", str_pad_val,side="right"), restriction_singletons=str_pad("Restriction to singletons", str_pad_val, side="right"))))
  inds <- 1:length(levels(plot_data$Analysis))
  
  # exponentiate if RR                                                                                                                                                                                                     pregnancy_level="Pregnancy-level", restriction_singletons="Restriction to singletons")) 
  if (exponentiate) {
    plot_data <- plot_data %>% mutate(result=exp(result))
  }
  
  # create plot
  plotted_figure <- plot_data %>%
    ggplot(aes(x = mediator_treatment_rr, y=result)) + 
    geom_line(aes(colour=Analysis, linetype=Analysis), alpha=0.5, linewidth=1.5) + 
    geom_point(aes(colour=Analysis, shape=Analysis), size=2) +
    scale_linetype_manual(values=c("solid", "dashed", "dotdash", "dotted", "longdash")[inds]) +
    scale_shape_manual(values=c(3,1,4,2,5)[inds]) +
    theme_minimal(base_size=12) +
    scale_y_continuous(trans="log", breaks=yticks, limits=c(ylimits[1],ylimits[2])) +
    xlab(xlabel) +
    ylab(ylabel)
  
  return(plotted_figure)
}


##############################
# Plot  comparing restriction to singletons to pregnancy-level analyses # 

# create plots ICC 0
for (analysis_pair in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
  plotted_results[[glue("trr{analysis_pair[1]}mrr{analysis_pair[2]}p0.03c0")]] <- plot_by_analysis(outputs, 0.03,analysis_pair[1],analysis_pair[2],0, "mean_estimate_", 
                                                                 xlabel="Treatment-multiples RR", 
                                                                 ylabel="Estimated treatment-outcome RR", 
                                                                 yticks=c(1,1.5,2,2.5,3),
                                                                 ylimits=c(0.95,3),
                                                                 exponentiate = TRUE, 
                                                                 exclude=c("gee_exchangeable", "naive_analysis")) + 
    geom_label(aes(x=2.4, y=2.7), label=list(bquote(atop(RR[AY*"|"*M] == .(analysis_pair[1]), textstyle(RR[MY*"|"*A] == .(analysis_pair[2]))))), parse=TRUE, size=5) + 
    theme(legend.text=element_text(size=14), legend.title = element_text(size = 14),
          axis.title=element_text(size=14), axis.text = element_text(size=10))
}


# create plots ICC 0.5
for (analysis_pair in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
  plotted_results[[glue("trr{analysis_pair[1]}mrr{analysis_pair[2]}p0.03c0.5")]] <- plot_by_analysis(outputs, 0.03,analysis_pair[1],analysis_pair[2],0.5, "mean_estimate_", 
                                                                                                   xlabel="Treatment-multiples RR", 
                                                                                                   ylabel="Estimated treatment-outcome RR", 
                                                                                                   yticks=c(1,1.5,2,2.5,3),
                                                                                                   ylimits=c(0.95,3),
                                                                                                   exponentiate = TRUE, 
                                                                                                   exclude=c("gee_exchangeable", "naive_analysis")) + 
    geom_label(aes(x=2.4, y=2.7), label=list(bquote(atop(RR[AY*"|"*M] == .(analysis_pair[1]), textstyle(RR[MY*"|"*A] == .(analysis_pair[2]))))), parse=TRUE, size=5) +
    theme(legend.text=element_text(size=14), legend.title = element_text(size = 14),
          axis.title=element_text(size=14), axis.text = element_text(size=10))
}
  
  
  

# create plots ICC 1
for (analysis_pair in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
  plotted_results[[glue("trr{analysis_pair[1]}mrr{analysis_pair[2]}p0.03c1")]] <- plot_by_analysis(outputs, 0.03,analysis_pair[1],analysis_pair[2],1, "mean_estimate_", 
                                                                                                   xlabel="Treatment-multiples RR", 
                                                                                                   ylabel="Estimated treatment-outcome RR", 
                                                                                                   yticks=c(1,1.5,2,2.5,3),
                                                                                                   ylimits=c(0.95,3),
                                                                                                   exponentiate = TRUE, 
                                                                                                   exclude=c("gee_exchangeable", "naive_analysis")) + 
    geom_label(aes(x=2.4, y=2.7), label=list(bquote(atop(RR[AY*"|"*M] == .(analysis_pair[1]), textstyle(RR[MY*"|"*A] == .(analysis_pair[2]))))), parse=TRUE, size=5) + 
    theme(legend.text=element_text(size=14), legend.title = element_text(size = 14),
          axis.title=element_text(size=14), axis.text = element_text(size=10))
}


# Output plot as JPEG (FIGURE 2)
jpeg(glue("outputs/rr_dodec_plot_rr.jpeg"), units="mm", width=300, height=400, res=300)
ggarrange(plotted_results$trr1mrr1p0.03c0, plotted_results$trr1mrr1p0.03c0.5,
          plotted_results$trr1mrr1p0.03c1, plotted_results$trr1mrr2p0.03c0,
          plotted_results$trr1mrr2p0.03c0.5,plotted_results$trr1mrr2p0.03c1,
          plotted_results$trr2mrr1p0.03c0, plotted_results$trr2mrr1p0.03c0.5,
          plotted_results$trr2mrr1p0.03c1, plotted_results$trr2mrr2p0.03c0, 
          plotted_results$trr2mrr2p0.03c0.5, plotted_results$trr2mrr2p0.03c1,
          nrow=4,ncol=3,common.legend=TRUE, labels=c("A","B","C", "D", "E", "F", "G","H", "I", "J","K", "L"), hjust=-0.2) 
dev.off()

for (analysis_pair in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
  plotted_results[[glue("trr{analysis_pair[1]}mrr{analysis_pair[2]}p0.10c0")]] <- plot_by_analysis(outputs, 0.10,analysis_pair[1],analysis_pair[2],0, "mean_estimate_", 
                                                                                                  xlabel="Treatment-multiples RR", 
                                                                                                  ylabel="Estimated treatment-outcome RR", 
                                                                                                  yticks=c(1,1.5,2,2.5,3),
                                                                                                  ylimits=c(0.95,3),
                                                                                                  exponentiate = TRUE, 
                                                                                                  exclude=c("gee_exchangeable", "naive_analysis"))
}

jpeg(glue("outputs/rr_quad_plot_rr_c1_out0.1.jpeg"), units="mm", width=200, height=200, res=300)
ggarrange(plotted_results$trr1mrr1p0.10c0, plotted_results$trr1mrr2p0.10c0,
          plotted_results$trr2mrr1p0.10c0, plotted_results$trr2mrr2p0.10c0,
          nrow=2,ncol=2,common.legend=TRUE, labels=c("A","B","C", "D"), hjust=-0.2)
dev.off()


# Size of estimated standard error ----------------------------------------

# plot with mediator-outcome RR of 2 and treatment-outcome RR of 2
plotted_results$both2_empSE <- plot_by_analysis(outputs, 2, 2, "empSE_estimate_", 
                                          xlabel="Treatment-multiples RR", 
                                          ylabel="Empirical standard error", 
                                          ylimits=c(0.05,0.15),
                                          exponentiate = FALSE)

# plot with mediator-outcome RR of 1
plotted_results$treatout2_empSE <- plot_by_analysis(outputs, 2, 1, "empSE_estimate_", 
                                              xlabel="Treatment-multiples RR", 
                                              ylabel="Empirical standard error", 
                                              ylimits=c(0.05,0.15),
                                              exponentiate = FALSE)

# plot with mediator-outcome RR of 2 and treatment-outcome RR of 1
plotted_results$mediatorout2_empSE <- plot_by_analysis(outputs, 1, 2, "empSE_estimate_", 
                                                       xlabel="Treatment-multiples RR", 
                                                       ylabel="Empirical standard error", 
                                                       ylimits=c(0.05,0.15),
                                                       exponentiate = FALSE)

# plot with mediator-outcome RR of 1 and treatment-outcome RR of 1
plotted_results$both1_empSE <- plot_by_analysis(outputs, 1, 1, "empSE_estimate_", 
                                                xlabel="Treatment-multiples RR", 
                                                ylabel="Empirical standard error", 
                                                ylimits=c(0.05,0.15),
                                                exponentiate = FALSE)


jpeg(glue("outputs/combined_plot_empSE_{FILE_SUFFIX}.jpeg"), units="mm", width=250, height=250, res=300)
ggarrange(plotted_results$both2_empSE + theme(axis.title.x=element_blank()), 
          plotted_results$treatout2_empSE + theme(axis.title.x=element_blank(), axis.title.y=element_blank()), 
          plotted_results$mediatorout2_empSE, 
          plotted_results$both1_empSE + theme(axis.title.y=element_blank()), 
          nrow=2, ncol=2, 
          common.legend=TRUE, labels=c("A","B","C","D"),
          hjust=-0.2) 
dev.off()


# Size of estimated standard error ----------------------------------------

# plot with mediator-outcome RR of 2 and treatment-outcome RR of 2
plotted_results$both2_relative_error_SE <- plot_by_analysis(outputs, 2, 2, "relative_error_SE_", 
                                                xlabel="Treatment-multiples RR", 
                                                ylabel="Relative error", 
                                                ylimits=c(90,110),
                                                exponentiate = FALSE)

# plot with mediator-outcome RR of 1
plotted_results$treatout2_relative_error_SE <- plot_by_analysis(outputs, 2, 1, "relative_error_SE_", 
                                                    xlabel="Treatment-multiples RR", 
                                                    ylabel="Relative error", 
                                                    ylimits=c(90,110),
                                                    exponentiate = FALSE)

# plot with mediator-outcome RR of 2 and treatment-outcome RR of 1
plotted_results$mediatorout2_relative_error_SE <- plot_by_analysis(outputs, 1, 2, "relative_error_SE_", 
                                                       xlabel="Treatment-multiples RR", 
                                                       ylabel="Relative error", 
                                                       ylimits=c(90,110),
                                                       exponentiate = FALSE)

# plot with mediator-outcome RR of 1 and treatment-outcome RR of 1
plotted_results$both1_relative_error_SE <- plot_by_analysis(outputs, 1, 1, "relative_error_SE_", 
                                                xlabel="Treatment-multiples RR", 
                                                ylabel="Relative error", 
                                                ylimits=c(90,110),
                                                exponentiate = FALSE)


jpeg(glue("outputs/combined_plot_relative_error_SE_{FILE_SUFFIX}.jpeg"), units="mm", width=250, height=250, res=300)
ggarrange(plotted_results$both2_relative_error_SE + theme(axis.title.x=element_blank()), 
          plotted_results$treatout2_relative_error_SE + theme(axis.title.x=element_blank(), axis.title.y=element_blank()), 
          plotted_results$mediatorout2_relative_error_SE, 
          plotted_results$both1_relative_error_SE + theme(axis.title.y=element_blank()), 
          nrow=2, ncol=2, 
          common.legend=TRUE, labels=c("A","B","C","D"),
          hjust=-0.2) 
dev.off()





quad plot
for (analysis_pair in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
  plotted_results[[glue("trr{analysis_pair[1]}mrr{analysis_pair[2]}p0.03c0")]] <- plot_by_analysis(outputs, 0.03,analysis_pair[1],analysis_pair[2],0, "mean_estimate_", 
                                                                 xlabel="Treatment-multiples RR", 
                                                                 ylabel="Estimated treatment-outcome RR", 
                                                                 yticks=c(1,1.5,2,2.5,3),
                                                                 ylimits=c(0.95,3),
                                                                 exponentiate = TRUE, 
                                                                 exclude=c("gee_exchangeable", "naive_analysis")) + 
    geom_label(aes(x=2.4, y=2.7), label=list(bquote(atop(RR[AY*"|"*M] == .(analysis_pair[1]), textstyle(RR[MY*"|"*A] == .(analysis_pair[2]))))), parse=TRUE, size=5) + 
    theme(legend.text=element_text(size=14), legend.title = element_text(size = 14),
          axis.title=element_text(size=14), axis.text = element_text(size=10))
}



jpeg(glue("outputs/rr_quad_plot_rr_c0.jpeg"), units="mm", width=200, height=200, res=300)
ggarrange(plotted_results$trr1mrr1p0.03c0, plotted_results$trr1mrr2p0.03c0,
          plotted_results$trr2mrr1p0.03c0, plotted_results$trr2mrr2p0.03c0,
          nrow=2,ncol=2,common.legend=TRUE, labels=c("A","B","C", "D"), hjust=-0.2)
dev.off()

  for (analysis_pair in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
    plotted_results[[glue("trr{analysis_pair[1]}mrr{analysis_pair[2]}p0.03c0.5")]] <- plot_by_analysis(outputs, 0.03,analysis_pair[1],analysis_pair[2],0.5, "mean_estimate_", 
                                                                                                     xlabel="Treatment-multiples RR", 
                                                                                                     ylabel="Estimated treatment-outcome RR", 
                                                                                                     yticks=c(1,1.5,2,2.5,3),
                                                                                                     ylimits=c(0.95,3),
                                                                                                     exponentiate = TRUE, 
                                                                                                     exclude=c("gee_exchangeable", "naive_analysis")) + 
      geom_label(aes(x=2.4, y=2.7), label=list(bquote(atop(RR[AY*"|"*M] == .(analysis_pair[1]), textstyle(RR[MY*"|"*A] == .(analysis_pair[2]))))), parse=TRUE, size=5) +
      theme(legend.text=element_text(size=14), legend.title = element_text(size = 14),
            axis.title=element_text(size=14), axis.text = element_text(size=10))
  }
  
  
  
  jpeg(glue("outputs/rr_quad_plot_rr_c0.5.jpeg"), units="mm", width=200, height=200, res=300)
  ggarrange(plotted_results$trr1mrr1p0.03c0.5, plotted_results$trr1mrr2p0.03c0.5,
            plotted_results$trr2mrr1p0.03c0.5, plotted_results$trr2mrr2p0.03c0.5,
            nrow=2,ncol=2,common.legend=TRUE, labels=c("A","B","C", "D"), hjust=-0.2)
  dev.off()


for (analysis_pair in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
  plotted_results[[glue("trr{analysis_pair[1]}mrr{analysis_pair[2]}p0.03c1")]] <- plot_by_analysis(outputs, 0.03,analysis_pair[1],analysis_pair[2],1, "mean_estimate_", 
                                                                                                   xlabel="Treatment-multiples RR", 
                                                                                                   ylabel="Estimated treatment-outcome RR", 
                                                                                                   yticks=c(1,1.5,2,2.5,3),
                                                                                                   ylimits=c(0.95,3),
                                                                                                   exponentiate = TRUE, 
                                                                                                   exclude=c("gee_exchangeable", "naive_analysis")) + 
    geom_label(aes(x=2.4, y=2.7), label=list(bquote(atop(RR[AY*"|"*M] == .(analysis_pair[1]), textstyle(RR[MY*"|"*A] == .(analysis_pair[2]))))), parse=TRUE, size=5) + 
    theme(legend.text=element_text(size=14), legend.title = element_text(size = 14),
          axis.title=element_text(size=14), axis.text = element_text(size=10))
}


jpeg(glue("outputs/rr_quad_plot_rr_c1.jpeg"), units="mm", width=200, height=200, res=300)
ggarrange(plotted_results$trr1mrr1p0.03c1, plotted_results$trr1mrr2p0.03c1,
          plotted_results$trr2mrr1p0.03c1, plotted_results$trr2mrr2p0.03c1,
          nrow=2,ncol=2,common.legend=TRUE, labels=c("A","B","C", "D"), hjust=-0.2)
dev.off()

jpeg(glue("outputs/rr_dodec_plot_rr.jpeg"), units="mm", width=300, height=400, res=300)
ggarrange(plotted_results$trr1mrr1p0.03c0, plotted_results$trr1mrr1p0.03c0.5,
          plotted_results$trr1mrr1p0.03c1, plotted_results$trr1mrr2p0.03c0,
          plotted_results$trr1mrr2p0.03c0.5,plotted_results$trr1mrr2p0.03c1,
          plotted_results$trr2mrr1p0.03c0, plotted_results$trr2mrr1p0.03c0.5,
          plotted_results$trr2mrr1p0.03c1, plotted_results$trr2mrr2p0.03c0, 
          plotted_results$trr2mrr2p0.03c0.5, plotted_results$trr2mrr2p0.03c1,
          nrow=4,ncol=3,common.legend=TRUE, labels=c("A","B","C", "D", "E", "F", "G","H", "I", "J","K", "L"), hjust=-0.2) 
dev.off()

for (analysis_pair in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
  plotted_results[[glue("trr{analysis_pair[1]}mrr{analysis_pair[2]}p0.10c0")]] <- plot_by_analysis(outputs, 0.10,analysis_pair[1],analysis_pair[2],0, "mean_estimate_", 
                                                                                                  xlabel="Treatment-multiples RR", 
                                                                                                  ylabel="Estimated treatment-outcome RR", 
                                                                                                  yticks=c(1,1.5,2,2.5,3),
                                                                                                  ylimits=c(0.95,3),
                                                                                                  exponentiate = TRUE, 
                                                                                                  exclude=c("gee_exchangeable", "naive_analysis"))
}

jpeg(glue("outputs/rr_quad_plot_rr_c1_out0.1.jpeg"), units="mm", width=200, height=200, res=300)
ggarrange(plotted_results$trr1mrr1p0.10c0, plotted_results$trr1mrr2p0.10c0,
          plotted_results$trr2mrr1p0.10c0, plotted_results$trr2mrr2p0.10c0,
          nrow=2,ncol=2,common.legend=TRUE, labels=c("A","B","C", "D"), hjust=-0.2)
dev.off()


# Size of estimated standard error ----------------------------------------

# plot with mediator-outcome RR of 2 and treatment-outcome RR of 2
plotted_results$both2_empSE <- plot_by_analysis(outputs, 2, 2, "empSE_estimate_", 
                                          xlabel="Treatment-multiples RR", 
                                          ylabel="Empirical standard error", 
                                          ylimits=c(0.05,0.15),
                                          exponentiate = FALSE)

# plot with mediator-outcome RR of 1
plotted_results$treatout2_empSE <- plot_by_analysis(outputs, 2, 1, "empSE_estimate_", 
                                              xlabel="Treatment-multiples RR", 
                                              ylabel="Empirical standard error", 
                                              ylimits=c(0.05,0.15),
                                              exponentiate = FALSE)

# plot with mediator-outcome RR of 2 and treatment-outcome RR of 1
plotted_results$mediatorout2_empSE <- plot_by_analysis(outputs, 1, 2, "empSE_estimate_", 
                                                       xlabel="Treatment-multiples RR", 
                                                       ylabel="Empirical standard error", 
                                                       ylimits=c(0.05,0.15),
                                                       exponentiate = FALSE)

# plot with mediator-outcome RR of 1 and treatment-outcome RR of 1
plotted_results$both1_empSE <- plot_by_analysis(outputs, 1, 1, "empSE_estimate_", 
                                                xlabel="Treatment-multiples RR", 
                                                ylabel="Empirical standard error", 
                                                ylimits=c(0.05,0.15),
                                                exponentiate = FALSE)


jpeg(glue("outputs/combined_plot_empSE_{FILE_SUFFIX}.jpeg"), units="mm", width=250, height=250, res=300)
ggarrange(plotted_results$both2_empSE + theme(axis.title.x=element_blank()), 
          plotted_results$treatout2_empSE + theme(axis.title.x=element_blank(), axis.title.y=element_blank()), 
          plotted_results$mediatorout2_empSE, 
          plotted_results$both1_empSE + theme(axis.title.y=element_blank()), 
          nrow=2, ncol=2, 
          common.legend=TRUE, labels=c("A","B","C","D"),
          hjust=-0.2) 
dev.off()


# Size of estimated standard error ----------------------------------------

# plot with mediator-outcome RR of 2 and treatment-outcome RR of 2
plotted_results$both2_relative_error_SE <- plot_by_analysis(outputs, 2, 2, "relative_error_SE_", 
                                                xlabel="Treatment-multiples RR", 
                                                ylabel="Relative error", 
                                                ylimits=c(90,110),
                                                exponentiate = FALSE)

# plot with mediator-outcome RR of 1
plotted_results$treatout2_relative_error_SE <- plot_by_analysis(outputs, 2, 1, "relative_error_SE_", 
                                                    xlabel="Treatment-multiples RR", 
                                                    ylabel="Relative error", 
                                                    ylimits=c(90,110),
                                                    exponentiate = FALSE)

# plot with mediator-outcome RR of 2 and treatment-outcome RR of 1
plotted_results$mediatorout2_relative_error_SE <- plot_by_analysis(outputs, 1, 2, "relative_error_SE_", 
                                                       xlabel="Treatment-multiples RR", 
                                                       ylabel="Relative error", 
                                                       ylimits=c(90,110),
                                                       exponentiate = FALSE)

# plot with mediator-outcome RR of 1 and treatment-outcome RR of 1
plotted_results$both1_relative_error_SE <- plot_by_analysis(outputs, 1, 1, "relative_error_SE_", 
                                                xlabel="Treatment-multiples RR", 
                                                ylabel="Relative error", 
                                                ylimits=c(90,110),
                                                exponentiate = FALSE)


jpeg(glue("outputs/combined_plot_relative_error_SE_{FILE_SUFFIX}.jpeg"), units="mm", width=250, height=250, res=300)
ggarrange(plotted_results$both2_relative_error_SE + theme(axis.title.x=element_blank()), 
          plotted_results$treatout2_relative_error_SE + theme(axis.title.x=element_blank(), axis.title.y=element_blank()), 
          plotted_results$mediatorout2_relative_error_SE, 
          plotted_results$both1_relative_error_SE + theme(axis.title.y=element_blank()), 
          nrow=2, ncol=2, 
          common.legend=TRUE, labels=c("A","B","C","D"),
          hjust=-0.2) 
dev.off()





