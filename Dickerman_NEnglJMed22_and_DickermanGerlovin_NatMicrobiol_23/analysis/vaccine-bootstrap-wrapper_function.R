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

wrapper_function <- function(quickcheck=FALSE, bsnum, cores, total.resamples=500, 
                             treat, event, variant = "alpha", 
                             matching.formula = paste0("group.binary ~ Age_at_index + sex + race + VISN + Caldate_num + urbanicity"), 
				 risk_time=168, days_fu=178, 
				 subgroup="none", 
                             ...){
  # Load the server/filepath specifications
  preamble(server=TRUE)
  if(subgroup=="none"){subgroup <- NULL} 

  # Set the number of bootstraps to run in the R session
  number.reps <- ceiling(as.numeric(total.resamples)/as.numeric(cores))
  if(quickcheck){number.reps=2}
  
  # Initialize log file
  if(length(subgroup)==0) {
	my_log <- file(paste0(prefix, "/r-logs/",treat,"-",event,"-", as.numeric(bsnum),"-",variant,"_1to1.txt")) 
  } else {
	my_log <- file(paste0(prefix, "/r-logs/",treat,"-",event,"-", as.numeric(bsnum),"-",variant,"-",subgroup,"_1to1.txt")) 
  }
  sink(my_log, append=TRUE)
  sink(my_log, append=TRUE, type="message")
  
  # Load the necessary libraries for each bootstrap process
  suppressPackageStartupMessages({library(tidyverse); library(reshape2); library(survival); library(splines)
    library(speedglm);  library(lubridate); library(data.table)})
  t1 <- Sys.time() #start timing the runs

bs_results <- lapply(1:number.reps,
                                  function(x) {
                                    source(paste0(prefix, "/vaccine-bootstrap-resample_function.R"))
                                    source(paste0(prefix, "/vaccine-analysis_function.R"))
                                    
                                    output <- bootstrap_resample_function(treat=treat, event=event,
										seed=number.reps*(as.numeric(bsnum)-1)+x,
										subgroup=subgroup,variant=variant) %>% 
                                      					vaccine_function(bootstrap=TRUE, bsdat=., subgroup=subgroup,      
												treat=treat, event=event, server=server,
												method="exact",matching.formula=matching.formula,
												risk_time=risk_time, days_fu=days_fu, variant=variant, 
                                                       						seed=number.reps*(as.numeric(bsnum)-1)+x, 
														...) 
                                    # options(warn=oldw)
                                    print(paste0("bootstrapnumber=",number.reps*(as.numeric(bsnum)-1)+x))
                                    return(output)
                                  }
)
t2 <- Sys.time() #stop the timer
print(t2-t1)
sink()
sink(type="message")
closeAllConnections()
}
