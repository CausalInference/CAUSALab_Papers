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

vaccine_function <- function(event, treat,  variant = "alpha", days_fu=178, risk_time=168, 
                             method=c("exact", "cem"),
                             matching.formula = paste0("group.binary ~ Age_at_index + sex + race + VISN + Caldate_num + urbanicity"), 
                             number.rows=Inf, bsdat=NULL,
                             seed=NULL, bootstrap=FALSE, negcontrol=TRUE,
                             KM_analysis=TRUE, 
                             subgroup=NULL, ...){
  
  # Run the 'prefix' specification - Note: this requires an explicit sourcing of the preamble.R program prior to running this function
  preamble(...); 
  
  # Load the function that strips away extra information from glm model output - this reduces the load on local memory when using the model object
  source(paste0(prefix, "/strip-glm_function.R"))
  
  # Setup the log for non-bootstrapped analytic calls on the SAS Grid side
  if(server & !bootstrap){
    my_log <<- file(paste0(prefix, "/r-logs/log-vaccine-analysis-", treat, "-", event,"-",method,"-",variant,"-",subgroup,"_1to1.txt"))
    sink(my_log, append=TRUE)
    sink(my_log, append=TRUE, type="message")
  }
  
  print(variant)
  
  # Subset the data for sensitivity analyses - Note: this process is already done to the bootstrapping samples prior to resampling, which is why this logic is not necessary to run 'again'
  if(!bootstrap){
    dat <- fread(file=paste0(prefix, "/data/cleaned-dat-", treat,"-",event,"-",variant,".csv"), nrows=number.rows)[!is.na(Age_at_index)] # removed observations with missing age
    if(!is.null(subgroup)){
      #subgroup <- paste0("-",subgroup)				
      if(subgroup=="white"){dat <- dat[race==0,]}
      if(subgroup=="black"){dat <- dat[race==1,]}
      if(subgroup=="old70"){dat <- dat[age>=70,]}
      if(subgroup=="young70"){dat <- dat[age<70,]}
      if(subgroup=="diff23_67"){dat <- dat[diff23_3grp==1,]}   
      if(subgroup=="diff23_8"){dat <- dat[diff23_3grp==2,]}   
      if(subgroup=="diff23_9"){dat <- dat[diff23_3grp==3,]}   
      if(subgroup=="nopriorcov"){dat <- dat[covidpriordose3==0,]}   
      if(subgroup=="trueboost"){dat <- dat[dose3_booster==1,]}   
      if(subgroup=="prim_pfizer"){dat <- dat[dose2_man=="Pfizer",]}   
      if(subgroup=="prim_moderna"){dat <- dat[dose2_man=="Moderna",]}   
    }
    
    
    # Exact matching process
    if(method=="exact"){
      set.seed(5) # Seed selected in our analyses. Changing the seed or failing to specify one will result in different estimates each time the analysis is performed.
      
      # Matching process
      m.out <- match.data(matchit(formula = as.formula(matching.formula),
                                  data=dat, 
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
    # Coarsened exact matching - Note: in this matching call, the cut-points are specified in a list, thus the list differs between analyses
    if(method=="cem"){
      set.seed(5) # Seed selected in our analyses. Changing the seed or failing to specify one will result in different estimates each time the analysis is performed.
      
      if(variant %in% c("alpha","delta")){
        m.out <- match.data(matchit(formula = as.formula(matching.formula),
                                    data=dat,
                                    method="cem",
                                    verbose=TRUE, #progress bar
                                    cutpoints = list(Age_at_index = 0,	# cutpoints option for exact matching on already-coarsened bins
                                                     sex = 0,
                                                     race = 0,
                                                     VISN = 0,
                                                     Caldate_num = 0,
                                                     urbanicity = 0),
                                    k2k=TRUE))[,.(newid, group.binary, outcome, newperiod, subclass, days)]
      } else if(variant %in% c("delom_booster","omicron_boost")){
        m.out <- match.data(matchit(formula = as.formula(matching.formula),
                                    data=dat,
                                    method="cem",
                                    verbose=TRUE, #progress bar
                                    cutpoints = list(Age_at_index = 0,	# cutpoints option for exact matching on already-coarsened bins
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
  } else if (bootstrap){ # receives dataset from vaccine-bootstrap-resample_function.R
    if(method=="exact"){
      m.out <- match.data(matchit(formula = as.formula(matching.formula),
                                  data=bsdat,
                                  method="exact",
                                  verbose=TRUE, #progress bar
                                  k2k=TRUE))[,.(newid, group.binary, outcome, newperiod, subclass, days)]	
      n_1 <- unique(m.out[group.binary==1,][,num_1 := .N, by=.(subclass)][, .(subclass, num_1)])[]
      n_0 <- unique(m.out[group.binary==0][, num_0 := .N, by=.(subclass)][, .(subclass, num_0)])[]
      n <- n_1[n_0, on=.(subclass)][, num := pmin(num_1, num_0)][]
      m.out.num <- n[m.out, on = .(subclass)][]
      m.out <- rbindlist(list(m.out.num[group.binary==0,][, .SD[sample(x=.N, size=num)], by=.(subclass)],
                              m.out.num[group.binary==1,][, .SD[sample(x=.N, size=num)], by=.(subclass)])) # creates N:N (per subclass) for exact matching
    }
    if(method=="cem"){
      if(variant %in% c("alpha","delta")){
        m.out <- match.data(matchit(formula = as.formula(matching.formula),
                                    data=bsdat,
                                    method="cem",
                                    verbose=TRUE, # progress bar
                                    cutpoints = list(Age_at_index = 0,	# added cutpoints option to ensure exact matching on already-coarsened bins
                                                     sex = 0,
                                                     race = 0,
                                                     VISN = 0,
                                                     Caldate_num = 0,
                                                     urbanicity = 0),
                                    k2k=TRUE))[,.(newid, group.binary, outcome, newperiod, subclass, days)]
      } else if(variant %in% c("delom_booster","omicron_boost")){
        m.out <- match.data(matchit(formula = as.formula(matching.formula),
                                    data=bsdat,
                                    method="cem",
                                    verbose=TRUE, # progress bar
                                    cutpoints = list(Age_at_index = 0,	# added cutpoints option to ensure exact matching on already-coarsened bins  
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
  } 
  print("matching done")
  
  # Output the matched dataset to a CSD
  if(!bootstrap & !is.null(subgroup) & event=="covidpos"){
    fwrite(m.out[,.(newid, group.binary, outcome, newperiod, subclass, days)], 
           file=paste0(prefix, "/data/matched-dat-",treat,"-",event,"-",method,"-",variant,"-",subgroup,"_1to1.csv"))
  }
  print("saving matched dataset done")
  
  # Calculate and print some descriptive statistics about number in each arm, number or events, and follow-up times
  output <- list() 
  
  if(!bootstrap){
    print(paste("Race table:",table(m.out$race)[1]))
    print(paste("Matched data, number of treat1 for",treat,"comparison:",table(m.out$group.binary)[2]))
    print(paste("Matched data, number of treat0 for",treat,"comparison:",table(m.out$group.binary)[1]))
    print(paste("Matched data, number of events, total:",table(m.out$outcome[m.out$days<=risk_time])[2]))
    print(paste("Matched data, number of events, in treat1:",table(m.out$outcome[m.out$group.binary==1 & m.out$days<=risk_time])[2]))
    print(paste("Matched data, number of events, in treat0:",table(m.out$outcome[m.out$group.binary==0 & m.out$days<=risk_time])[2]))
    print(paste("Matched data, follow-up time, median:",quantile(pmin(m.out$days,risk_time),0.50))) 
    print(paste("Matched data, follow-up time, 25th percentile:",quantile(pmin(m.out$days,risk_time),0.25)))
    print(paste("Matched data, follow-up time, 75th percentile:",quantile(pmin(m.out$days,risk_time),0.75))) 
    
    #grab results easily -- this is used for outputting data that can later be used for summary statistics
    output$descriptive <- list() 
    output$descriptive$n <- nrow(m.out) # total number
    output$descriptive$treat1 <- table(m.out$group.binary)[["1"]]
    output$descriptive$treat0 <- table(m.out$group.binary)[["0"]]
    
    output$descriptive$events_total <- output$descriptive$events_treat1 <- output$descriptive$events_treat0 <- list()
    output$descriptive$events_total[as.character(risk_time)] <- sum(m.out$outcome[m.out$days<=risk_time])
    output$descriptive$events_treat1[as.character(risk_time)] <- sum(m.out$outcome[m.out$group.binary==1 & m.out$days<=risk_time])
    output$descriptive$events_treat0[as.character(risk_time)] <- sum(m.out$outcome[m.out$group.binary==0 & m.out$days<=risk_time])
    
    output$descriptive$events_total["10"] <- sum(m.out$outcome[m.out$days<=10])
    output$descriptive$events_treat1["10"] <- sum(m.out$outcome[m.out$group.binary==1 & m.out$days<=10])
    output$descriptive$events_treat0["10"] <- sum(m.out$outcome[m.out$group.binary==0 & m.out$days<=10])
    
    output$descriptive$events_total["7"] <- sum(m.out$outcome[m.out$days<=7])
    output$descriptive$events_treat1["7"] <- sum(m.out$outcome[m.out$group.binary==1 & m.out$days<=7])
    output$descriptive$events_treat0["7"] <- sum(m.out$outcome[m.out$group.binary==0 & m.out$days<=7])
    
    output$descriptive$t_median <- paste("Matched data, follow-up time, median:",quantile(pmin(m.out$days,risk_time),0.50))   
    output$descriptive$t_q1 <- paste("Matched data, follow-up time (days), 25th percentile:",quantile(pmin(m.out$days,risk_time),0.25)) 
    output$descriptive$t_q3 <- paste("Matched data, follow-up time (days), 75th percentile:",quantile(pmin(m.out$days,risk_time),0.75))
    
  }
  
  # Estimate the Kaplan-Meier survival curves
  if(KM_analysis) {			
    km.data <- survfit(Surv(days, outcome) ~ group.binary, data=m.out %>% 
                         mutate(group.binary = factor(group.binary)) %>% 
                         mutate(group.binary = forcats::fct_relevel(group.binary, "0", after=Inf)))
    if(!bootstrap){
      output$km.risk <- 1-summary(km.data, times=risk_time,extend=TRUE)$surv #for treat1, for treat0
      output$km.risk_diff <- output$km.risk[1]-output$km.risk[2] #risk treat1 - risk treat0
      output$km.risk_ratio <- output$km.risk[1]/output$km.risk[2] #risk treat1 / risk treat0
      if(negcontrol){
        output$km.risk_10 <- 1-summary(km.data, times=10,extend=TRUE)$surv #for treat1, for treat0
        output$km.risk_diff_10 <- output$km.risk_10[1]-output$km.risk_10[2] #risk treat1 - risk treat0
        output$km.risk_ratio_10 <- output$km.risk_10[1]/output$km.risk_10[2] #risk treat1 / risk treat0
        output$km.risk_7 <- 1-summary(km.data, times=7,extend=TRUE)$surv #for treat1, for treat0 
        output$km.risk_diff_7 <- output$km.risk_7[1]-output$km.risk_7[2] #risk treat1 - risk treat0 
        output$km.risk_ratio_7 <- output$km.risk_7[1]/output$km.risk_7[2] #risk treat1 / risk treat0 
        
      }    	
      
      # save point estimates of risk at all times t, for purpose of plots and tables
      km.data2 <-summary(km.data,times=0:days_fu,extend=TRUE)
      output$km$point_est  <- data.frame(t=km.data2$time, risk=(1-km.data2$surv), 
                                         group=c(rep(1, times=days_fu+1),
                                                 rep(0, times=days_fu+1)))
      fwrite(output$km$point_est, file=paste0(prefix, "/data/km-point-est-",treat,"-",event,"-",method,"-",variant,subgroup,"_1to1.csv"))
    } else if (bootstrap){
      # extend to all times t, and output risk instead of surv 
      km.data2 <-summary(km.data,times=0:days_fu,extend=TRUE)
      output$bs.km <- data.frame(t=km.data2$time, risk=(1-km.data2$surv), 
                                 group=c(rep(1, times=days_fu+1),
                                         rep(0, times=days_fu+1))) 
    }
  }
  
  # Output the final results to an Rdata file - these objects contain estimates, descriptive statistics, and other information needed when aggregating results
  if(!bootstrap){
    print(output)
    ifelse(!dir.exists(file.path(prefix, "results")), dir.create(file.path(prefix, "results")), FALSE)
    save(output, file=paste0(prefix, "/results/output-",treat,"-",event,"-",method,"-",variant,subgroup,"_1to1.Rda"))
    if(server){sink(); sink(type="message"); closeAllConnections()}
  }
  
  # Save bootstrap output in the bootstrap folder
  if(bootstrap){
    save(output, file=paste0(prefix, "/bootstrap/output-",treat,"-",event,"-vaccines-",seed,"-",variant,subgroup,"_1to1.Rda")) # to save each resample separately
  }
}
