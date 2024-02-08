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

main_estimates <- function(treat, event, variant,
                           risk_time,
                           cde=FALSE, subgroup=NULL) {
  
  ## load the point estimate results file
  load(paste0(prefix, "/results/output-", treat, "-", event, "-exact-", variant,
              if(!is.null(subgroup)){subgroup},
              if(cde){"-cde"}, 
              "_1to1.Rda"))
  
  ## empiric sample estimates
  if(treat=="PvM"){
    descriptive <- data.frame("n"=output$descriptive[["n"]], 
                              "bnt-num-events"=
                                ifelse(is.null(output$descriptive$events_treat1[[as.character(risk_time)]]),
                                       "0",
                                       output$descriptive$events_treat1[[as.character(risk_time)]] ), 
                              "moderna-num-events"=
                                ifelse(is.null(output$descriptive$events_treat0[[as.character(risk_time)]]),
                                       "0",
                                       output$descriptive$events_treat0[[as.character(risk_time)]] )
    )
  }
  if(treat=="MvP"){
    descriptive <- data.frame("n"=output$descriptive[["n"]], 
                              "moderna-num-events"=
                                ifelse(is.null(output$descriptive$events_treat1[[as.character(risk_time)]]),
                                       "0",
                                       output$descriptive$events_treat1[[as.character(risk_time)]] ), 
                              "bnt-num-events"=
                                ifelse(is.null(output$descriptive$events_treat0[[as.character(risk_time)]]),
                                       "0",
                                       output$descriptive$events_treat0[[as.character(risk_time)]] )
    )
  }
  
  tempa <- list()
  
  tempa$point_est <- output$km$point_est
  
  tempa$risk_diff <- output$km$point_est %>% 
    pivot_wider(., names_from=group, values_from=risk) %>% 
    mutate(riskdiff=`1`-`0`) %>% 
    dplyr::select(t, riskdiff)
  
  tempa$risk_ratio <- output$km$point_est %>% 
    pivot_wider(., names_from=group, values_from=risk) %>% 
    mutate(riskratio=`1`/`0`) %>% 
    dplyr::select(t, riskratio)
  
  tempa$risk_ratio$riskratio[tempa$risk_ratio$t==0] <- 1
  
  tempa$descriptive <- descriptive
  
  return(tempa)
  
  
}

bs_results_summary <- function(number.resamples=500, 
                               calc_point_est=TRUE, calc_risk_diff=TRUE, 
                               treat, event, variant,
                               risk_time=(18*7),
                               cde=FALSE, subgroup=NULL){
  
  # Get the primary point estimate results
  mainres <- main_estimates(treat=treat, event=event, variant=variant, risk_time=risk_time, cde=cde, subgroup=subgroup)
  point_est <- mainres$point_est
  risk_diff <- mainres$risk_diff
  risk_ratio <- mainres$risk_ratio
  descriptive <- mainres$descriptive
  rm(mainres)
  
  bs_results <- list()
  ## collect up the bootstrapped results
    bs_combined <- lapply(1:number.resamples, function(x) {
      
      #load an individual dataset, x
      load(paste0(prefix, "/bootstrap/output-", treat, "-", event, "-vaccines-", x, 
                  "-", variant,
                  ifelse(!is.null(subgroup), subgroup, ""),
                  ifelse(cde, "-cde", ""),
                  "_1to1.Rda"))
      
      temp <- list()
      
      # grab the estimates
      temp$ests <- output$bs.km$risk
      
      temp$rd <- output$bs.km %>% 
        pivot_wider(., names_from=group, values_from=risk) %>% 
        mutate(riskdiff=`1`-`0`) %>% 
        .$riskdiff
      
      bs_risk_ratio <- output$bs.km %>% 
        pivot_wider(., names_from=group, values_from=risk) %>% 
        mutate(riskratio=`1`/`0`) %>% 
        .$riskratio
        
      bs_risk_ratio[1] <- 1
      
      temp$rr <- bs_risk_ratio
      
      return(temp)
    })
    
    bs_ests_0 <- lapply(1:number.resamples, function(x){
      bs_combined[[x]]$ests
    } ) %>% data.frame
    
    assertthat::are_equal(length(bs_ests_0[50,] %>% as.numeric), 
                          length(bs_ests_0[50,] %>% as.numeric %>% unique)) # check no issue with unique seed in parallel bootstrap
    assertthat::are_equal(length(bs_ests_0[400,] %>% as.numeric), 
                          length(bs_ests_0[400,] %>% as.numeric %>% unique))
    
    colnames(bs_ests_0) <- NULL
    
    # percentile-based bootstrap bounds at each time point
    bs_ests_upper.percentile <- apply(bs_ests_0, 1, quantile, probs=0.975) 
    bs_ests_lower.percentile <- apply(bs_ests_0, 1, quantile, probs=0.025)
    
    # bring the point estimates and bootstraps together
    result_point_est0 <- cbind(point_est %>% 
                                 mutate(group.binary=group,
                                        group=case_when(treat=="PvM" & group==1 ~ "BNT162b2",
                                                        treat=="PvM" & group==0 ~ "mRNA-1273",
                                                        treat=="MvP" & group==0 ~ "BNT162b2",
                                                        treat=="MvP" & group==1 ~ "mRNA-1273")), 
                               "lcl"=bs_ests_lower.percentile, 
                               "ucl"=bs_ests_upper.percentile
                               ) %>%
      filter(t<=risk_time)
    
    groups <- table(result_point_est0$group) %>% names
    
    ref <- unique(result_point_est0$group[result_point_est0$group.binary==0])
    
    contrast <- unique(result_point_est0$group[result_point_est0$group.binary==1])
    
    result_point_est0$group.binary <- NULL
    
    bs_results$result_point_est <- result_point_est0 %>% 
      pivot_wider(., names_from=group, values_from=c(risk, lcl, ucl), names_sep="_") %>% 
      dplyr::select_at(vars("t", contains(groups[1]), contains(groups[2]))) %>% 
      mutate(ref=ref, contrast=contrast)
    
  ## end point est
    
  ## rd and rr
    bs_rd_0 <- lapply(1:number.resamples, function(x){
      bs_combined[[x]]$rd
    } ) %>% data.frame
    bs_rr_0 <- lapply(1:number.resamples, function(x){
      bs_combined[[x]]$rr
    } ) %>% data.frame
    
    colnames(bs_rd_0) <- NULL
    
    # percentile-based risk-differences at each time-point
    bs_rd_upper.percentile <- apply(bs_rd_0, 1, quantile, probs=0.975)
    bs_rd_lower.percentile <- apply(bs_rd_0, 1, quantile, probs=0.025)
    
    colnames(bs_rr_0) <- NULL
    
    # percentile-based risk-ratio at each time-point
    bs_rr_upper.percentile <- apply(bs_rr_0, 1, quantile, probs=0.975, na.rm=TRUE) # removes early undefined risk ratios
    bs_rr_lower.percentile <- apply(bs_rr_0, 1, quantile, probs=0.025, na.rm=TRUE)
    
    bs_results$result_risk_diff <- cbind(risk_diff, 
                                         # "risk_diff_lcl"=risk_diff$riskdiff - bs_rd_lower.empiric, 
                                         # "risk_diff_ucl"=risk_diff$riskdiff - bs_rd_upper.empiric,
                                         "risk_diff_lcl"=bs_rd_lower.percentile, 
                                         "risk_diff_ucl"=bs_rd_upper.percentile) %>% 
      filter(t<=risk_time) %>% 
      rename(!!(paste0("riskdiff_ref=",ref)) := riskdiff)
    
    bs_results$result_risk_ratio <- cbind(risk_ratio, 
                                          # "risk_ratio_lcl"=risk_ratio$riskratio - bs_rr_lower.empiric, 
                                          # "risk_ratio_ucl"=risk_ratio$riskratio - bs_rr_upper.empiric, 
                                         "risk_ratio_lcl"=bs_rr_lower.percentile,
                                         "risk_ratio_ucl"=bs_rr_upper.percentile) %>% 
      filter(t<=risk_time) %>% 
      rename(!!(paste0("riskratio_ref=",ref)) := riskratio)
    ## end rd and rr
    
    combined_result <- cbind(descriptive,
                             left_join(
      left_join(bs_results$result_risk_diff,
                bs_results$result_risk_ratio, by="t"),
      bs_results$result_point_est,
      by="t") %>% filter(t==risk_time)
    )
    
    return(combined_result)
}
