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

table_2_function <- function(cde=FALSE, treat, subgroup=NULL, event, variant, risk_time){
  
  combined_results <- bs_results_summary(treat=treat, event=event, variant=variant, cde=cde, 
                                         subgroup=subgroup, risk_time=risk_time)
  
  # sample size and number of events for each intervention
  num <- combined_results %>% dplyr::select(n, bnt.num.events, moderna.num.events)
  # risk
  risk <- risk_point_ests_function(d=combined_results) %>% as.data.frame %>% dplyr::select("BNT162b2","mRNA-1273")
  
  # risk difference and risk ratio
  rd_rr <- rd_rr_function(d=combined_results, treat=treat)
  
  subg <- ifelse(is.null(subgroup),"none",subgroup)
  combined <- cbind(variant=variant, subgroup=subg,event=event, num, risk, rd_rr)
  return(combined)
}

parallel_table_2_function <- function(cde=FALSE, treat, subgroup=NULL, 
                                      variant, risk_time,
                                      outlist=c("covidpos","covidpossymp","covidhosp","covidICU","COVIDdeath","notcoviddeath")){
  
  temp <- lapply(1:length(outlist),
                   function(numx)
                     table_2_function(treat=treat,
                                      event=outlist[numx],
                                      subgroup=subgroup,
                                      cde=cde, 
                                      variant=variant,
                                      risk_time=risk_time))
  
  temp2 <- do.call(rbind,temp)
  
  return(temp2)
}
