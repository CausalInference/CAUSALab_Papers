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

rd_rr_function <- function(d, treat){
    rd <- paste0(d %>% 
                   dplyr::select_at(vars(contains("riskdiff"))) %>% 
                   magrittr::multiply_by(., 1e4) %>% 
                   sprintf("%.1f", .),
                 " (",
                 d$risk_diff_lcl %>% 
                   magrittr::multiply_by(., 1e4) %>%
                   sprintf("%.1f", .),
                 ", ",
                 d$risk_diff_ucl %>% 
                   magrittr::multiply_by(., 1e4) %>% 
                   sprintf("%.1f", .),
                 ")")
    
    rr <- paste0(d %>% 
                   dplyr::select_at(vars(contains("riskratio"))) %>% 
                   sprintf("%.2f", .),
                 " (",
                 d$risk_ratio_lcl %>% 
                   sprintf("%.2f", .),
                 ", ",
                 d$risk_ratio_ucl %>% 
                   sprintf("%.2f", .),
                 ")")
    
    temp <- c(rd, rr)
    
    temp2 <- temp %>% 
      data.frame %>% 
      t() %>% 
      data.frame
    
  colnames(temp2) <- c(paste0(treat, "-RD"), paste0(treat, "-RR"))
  
  return(temp2)
}

