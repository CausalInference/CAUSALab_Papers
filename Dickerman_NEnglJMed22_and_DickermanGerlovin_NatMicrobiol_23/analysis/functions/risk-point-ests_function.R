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

# outputs risk table with point estimates for manuscript
risk_point_ests_function <- function(d){
    ref <- d$ref
    contrast <- d$contrast
    ref_ests <- paste0(d %>% 
                         dplyr::select_at(vars(contains(ref))) %>% 
                         dplyr::select_at(vars(contains("risk_"))) %>% 
                         magrittr::multiply_by(., 1e4) %>% 
                         .[,1] %>% 
                         sprintf("%.1f", .),
                       " (",
                       d %>% 
                         dplyr::select_at(vars(contains(ref))) %>% 
                         dplyr::select_at(vars(contains("lcl"))) %>% 
                         magrittr::multiply_by(., 1e4) %>% 
                         .[,1] %>% 
                         sprintf("%.1f", .),
                       ", ",
                       d %>% 
                         dplyr::select_at(vars(contains(ref))) %>% 
                         dplyr::select_at(vars(contains("ucl"))) %>% 
                         magrittr::multiply_by(., 1e4) %>% 
                         .[,1] %>% 
                         sprintf("%.1f", .),
                       ")")
    
    contrast_ests <- paste0(d %>% 
                              dplyr::select_at(vars(contains(contrast))) %>% 
                              dplyr::select_at(vars(contains("risk_"))) %>% 
                              magrittr::multiply_by(., 1e4) %>% 
                              .[,1] %>% 
                              sprintf("%.1f", .),
                            " (",
                            d %>% 
                              dplyr::select_at(vars(contains(contrast))) %>% 
                              dplyr::select_at(vars(contains("lcl"))) %>% 
                              magrittr::multiply_by(., 1e4) %>% 
                              .[,1] %>% 
                              sprintf("%.1f", .),
                            ", ",
                            d %>% 
                              dplyr::select_at(vars(contains(contrast))) %>% 
                              dplyr::select_at(vars(contains("ucl"))) %>% 
                              magrittr::multiply_by(., 1e4) %>% 
                              .[,1] %>% 
                              sprintf("%.1f", .),
                            ")")
    
    rbind(ref_ests, contrast_ests)
    
    temp <- data.frame(ref_ests, contrast_ests) %>% 
      rename(!!d$ref:="ref_ests", !!d$contrast:="contrast_ests")
    
  return(temp)
}


  