#####
#
# This script is a part of a set of analytic programs related to the following manuscript(s):
#  - Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
#  - Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors. 
#
# Author(s): VA-CAUSAL Methods Core
# Version: July 2023.
#
######

plot_function <- function(number.resamples=500, 
                          treat, event, variant,
                          risk_time=(24*7),
                          cde=FALSE, subgroup="none",
                          yaxis=0.006, ybreak=0.001, yaccuracy=0.1){
  # load point estimates

  if(subgroup=="none"){subgroup <- NULL}
  load(paste0(prefix, "/results/output-", treat, "-", event, "-exact-", variant,
              if(!is.null(subgroup)){subgroup},
              if(cde){"-cde"}, 
              "_1to1.Rda"))
  point_est <- output$km$point_est
  
  # load bootstrap estimates
  bs_ests_0 <- pbapply::pblapply(1:number.resamples, function(x) {
    load(paste0(prefix, "/bootstrap/output-", treat, "-", event, "-vaccines-", x, 
                "-", variant,
                ifelse(!is.null(subgroup), subgroup, ""),
                ifelse(cde, "-cde", ""),
                "_1to1.Rda"))
    output$bs.km$risk
  }) %>% as.data.frame
  colnames(bs_ests_0) <- NULL
  
	xbreaker<-ifelse(risk_time<14,1,
				ifelse(risk_time%%14==0,14,7)) #add breaking by 14 days if follow-up/risk time plot uses 14+ days
 
  bs_ests_upper.percentile <- apply(bs_ests_0, 1, quantile, probs=0.975) # percentile based bootstrap
  bs_ests_lower.percentile <- apply(bs_ests_0, 1, quantile, probs=0.025)
  
  point_est <- point_est %>% mutate(group=case_when(treat=="PvM" & group==1 ~ "BNT162b2",
                                                    treat=="PvM" & group==0 ~ "mRNA-1273",
                                                    treat=="MvP" & group==0 ~ "BNT162b2",
                                                    treat=="MvP" & group==1 ~ "mRNA-1273"))
  result <- cbind(point_est, 
                  "lcl"=bs_ests_lower.percentile, 
                  "ucl"=bs_ests_upper.percentile) %>% 
    filter(t<=risk_time)

	write.csv(result,paste0(prefix,"/plots/SD_", treat, "-", event, "-", variant,
              if(!is.null(subgroup)){subgroup},"-fup",risk_time,".csv"))  #writing out source data for plot per request from Nature Microbiology paper
  
  event_plotname <- fcase(event== "covidpos" , "Documented SARS-CoV-2 Infection",
                          event== "covidpossymp" , "Symptomatic Covid-19",
                          event== "covidhosp" , "Covid-19 Hospitalization",
                          event== "covidICU" ,"Covid-19 ICU Admission",
                          event== "COVIDdeath" , "Covid-19 Death",
                          event== "notcoviddeath" , "Non-Covid-19 Death"
  )  
  
  
  #define legend labels based on supplied treat variable 
  legend_treat1 <- fcase(treat=="PvM","BNT162b2", 
                         treat== "MvP","mRNA-1273")  
  legend_treat0 <- fcase(treat=="PvM","mRNA-1273", 
                         treat== "MvP", "BNT162b2")  
  
  
  ggplot(data=result,
         mapping=aes(x=t, y=risk, group=as.factor(group))) +
    geom_line(size=0.5, aes(color=as.factor(group), linetype=as.factor(group))) +
    geom_ribbon(aes(x=t,ymin=lcl,ymax=ucl,fill=as.factor(group)),alpha=0.15,show.legend = FALSE) +
    scale_y_continuous(breaks = seq(0,yaxis,by=ybreak),
                       limits = c(0,yaxis),
                       name = "Cumulative incidence",
                       labels = scales::percent_format(accuracy=yaccuracy),
                       expand=expansion(mult=c(0.01,0.01)))+
    scale_x_continuous(breaks = seq(0,risk_time,by=xbreaker),
                       limits = c(0,risk_time),
                       name = "Days",
                       expand=expansion(mult=c(0.01,0.01))) +
    ggtitle(event_plotname) +
    theme_classic(base_size=9) +
    theme(legend.justification = c(0,0), 
          legend.position = c(0.05, 0.6),
          legend.title=element_blank(),
          axis.title = element_text(size=9),
          plot.title = element_text(size=9)) +
    scale_color_manual(labels=c(legend_treat1,legend_treat0),values=c("#374E55FF","#DF8F44FF"),breaks=c(legend_treat1,legend_treat0)) +
    scale_linetype_manual(labels=c(legend_treat1,legend_treat0),values=c("solid","dashed"),breaks=c(legend_treat1,legend_treat0))+
    scale_fill_manual(labels=c(legend_treat1,legend_treat0),values=c("#374E55FF","#DF8F44FF"),breaks=c(legend_treat1,legend_treat0))
  
}
