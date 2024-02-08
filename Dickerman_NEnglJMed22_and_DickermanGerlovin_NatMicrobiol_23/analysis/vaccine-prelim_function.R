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

prelim_function <- function(event="covidpos", treat="PvM", days_fu=7, risk_time=168, 
                             method=c("exact", "cem"),
                             matching.formula = paste0("group.binary ~ Age_at_index + sex + race + VISN + Caldate_num + urbanicity"), 
                             number.rows=Inf, variant="alpha" 
                             , ...){

  # Load the server/path specifications for the R session
  preamble(...); source(paste0(prefix, "/strip-glm_function.R"))

  # Initialize an r-log
  if(server){
    my_log <<- file(paste0(prefix, "/r-logs/log-vaccine-prelim-", treat, "-", event,"-",method,"-",variant,"_1to1.txt"))
    sink(my_log, append=TRUE)
    sink(my_log, append=TRUE, type="message")
  }

  print(variant)

   #define legend labels based on supplied treat variable 
   legend_treat1 <- fcase(treat=="PvM","BNT162b2", 
                   treat=="MvP","mRNA-1273") 

   legend_treat0 <- fcase(treat=="PvM","mRNA-1273", 
                treat=="MvP","BNT162b2")

   # Perform matching to allow for preliminary checks
     if(method=="exact"){
	set.seed(5) 
      m.out <- match.data(matchit(formula = as.formula(matching.formula),
                                  data=fread(file=paste0(prefix, "/data/cleaned-dat-", treat,"-",event,"-",variant,".csv"), nrows=number.rows) %>% 
                                    filter(!is.na(Age_at_index)), # removed n=2 observations with missing age
                                  method="exact",
                                  verbose=TRUE, #progress bar
                                  k2k=TRUE))[,.(newid, group.binary, outcome, newperiod, subclass, days)]         
      n_1 <- unique(m.out[group.binary==1,][,num_1 := .N, by=.(subclass)][, .(subclass, num_1)])[]
      n_0 <- unique(m.out[group.binary==0][, num_0 := .N, by=.(subclass)][, .(subclass, num_0)])[]
      n <- n_1[n_0, on=.(subclass)][, num := pmin(num_1, num_0)][]
      m.out.num <- n[m.out, on = .(subclass)][]
      m.out <- rbindlist(list(m.out.num[group.binary==0,][, .SD[sample(x=.N, size=num)], by=.(subclass)],
                              m.out.num[group.binary==1,][, .SD[sample(x=.N, size=num)], by=.(subclass)]))
    }
    if(method=="cem"){
	set.seed(5) 

		if(variant %in% c("alpha","delta")) {
      		m.out <- match.data(matchit(formula = as.formula(matching.formula),
                                  data=fread(file=paste0(prefix, "/data/cleaned-dat-", treat,"-",event,"-",variant,".csv"), nrows=number.rows) %>% 
                                    filter(!is.na(Age_at_index)), 
                                  method="cem",
                                  verbose=TRUE, #progress bar
				      cutpoints = list(Age_at_index = 0,	
						     sex = 0,
						     race = 0,
						     VISN = 0,
						     Caldate_num = 0,
						     urbanicity = 0), 
                                  k2k=TRUE))[,.(newid, group.binary, outcome, newperiod, subclass, days)]	 

    	} else if(variant %in% c("delom_booster","omicron_boost")){
      		m.out <- match.data(matchit(formula = as.formula(matching.formula),
                                  data=fread(file=paste0(prefix, "/data/cleaned-dat-", treat,"-",event,"-",variant,".csv"), nrows=number.rows) %>% 
                                    filter(!is.na(Age_at_index)), 
                                  method="cem",
                                  verbose=TRUE, #progress bar
				      cutpoints = list(Age_at_index = 0,	
						     sex = 0,
						     race = 0,
						     VISN = 0,
						     Caldate_num = 0,
						     urbanicity = 0,
							 dose2_month = 0, 
							 ntests = 0), 
                                  k2k=TRUE))[,.(newid, group.binary, outcome, newperiod, subclass, days)]	 

    }
	}

  print("matching done")
  

	# May want to uncomment the bracketed lines to not save the datasets for all analyses
   #if(event=="covidpos"){
		fwrite(m.out[,.(newid, group.binary, outcome, newperiod, subclass, days)], 
		file=paste0(prefix, "/data/matched-dat-",treat,"-",event,"-",method,"-",variant,"_1to1.csv"))
   #}
  print("saving matched dataset done")
  
   # Save descriptives to log for preliminary checks
  print(paste("Matched data, number of treat1 for",treat,"comparison:",table(m.out$group.binary)[2]))
  print(paste("Matched data, number of treat0 for",treat,"comparison:",table(m.out$group.binary)[1]))
  print(paste("Matched data, number of events, total:",table(m.out$outcome[m.out$days<=risk_time])[2]))
  print(paste("Matched data, number of events, in treat1:",table(m.out$outcome[m.out$group.binary==1 & m.out$days<=risk_time])[2]))
  print(paste("Matched data, number of events, in treat0:",table(m.out$outcome[m.out$group.binary==0 & m.out$days<=risk_time])[2]))
  print(paste("Matched data, follow-up time, median:",quantile(pmin(m.out$days,risk_time),0.50))) 
  print(paste("Matched data, follow-up time, 25th percentile:",quantile(pmin(m.out$days,risk_time),0.25))) 
  print(paste("Matched data, follow-up time, 75th percentile:",quantile(pmin(m.out$days,risk_time),0.75))) 
  print(paste("Matched data (treat=1), follow-up time, median:",quantile(pmin(m.out$days[m.out$group.binary==1],risk_time),0.50)))  
  print(paste("Matched data (treat=1), follow-up time, 25th percentile:",quantile(pmin(m.out$days[m.out$group.binary==1],risk_time),0.25))) 
  print(paste("Matched data (treat=1), follow-up time, 75th percentile:",quantile(pmin(m.out$days[m.out$group.binary==1],risk_time),0.75))) 
  print(paste("Matched data (treat=0), follow-up time, median:",quantile(pmin(m.out$days[m.out$group.binary==0],risk_time),0.50))) 
  print(paste("Matched data (treat=0), follow-up time, 25th percentile:",quantile(pmin(m.out$days[m.out$group.binary==0],risk_time),0.25))) 
  print(paste("Matched data (treat=0), follow-up time, 75th percentile:",quantile(pmin(m.out$days[m.out$group.binary==0],risk_time),0.75))) 

    # Check the Kaplan-Meier curves for first 7-10 days (negative control, no expected difference)
    km.data <- survfit(Surv(days, outcome) ~ group.binary, data=m.out %>% 
                         mutate(group.binary = factor(group.binary)) %>% 
                         mutate(group.binary = forcats::fct_relevel(group.binary, "0", after=Inf)))
  
    mypal=ggsci::pal_jama("default",alpha=1)(2)

    neg.control.km.plot <- suppressWarnings(survminer::ggsurvplot(km.data, data=m.out, 
                                                             palette = c(mypal[1], mypal[2]),
                                                             size=1,
                                                             censor=FALSE,
                                                             fun="event",
                                                             linetype=c(1,2),
                                                             font.main=24,
                                                             legend=c(0.12, 0.8),
                                                             legend.title="",
                                                             legend.labs=c(paste0(legend_treat1), paste0(legend_treat0)),
                                                             risk.table=TRUE, 
                                                             risk.table.height = 0.25,
                                                             tables.y.text.col=FALSE,
                                                             fontsize=10,
                                                             xlim=c(0,days_fu),
                                                             break.x.by=1, 
                                                             xlab="Days", 
                                                             ylim=c(0,0.005),
                                                             surv.scale="percent",
                                                             ylab="Cumulative incidence (%)",
                                                             ggtheme=theme_survminer(font.main = c(22, "black"),
                                                                                     font.submain = c(22, "black"),
                                                                                     font.caption = c(22, "black"),
                                                                                     font.legend = c(22,"black"),
                                                                                     font.x = c(22, "black"),
                                                                                     font.y = c(22, "black"),
                                                                                     font.tickslab = c(22, "black")),
                                                             tables.theme=theme_cleantable()))


    #save negative control KM plot (pdf)
    pdf(paste0(prefix,"/plots/figure-neg-control-km-plot_",treat,"-",event,"-",method,"-",variant,"_1to1.pdf"),width=14, height=7)
    print(neg.control.km.plot,newpage=FALSE)
    dev.off()

}
