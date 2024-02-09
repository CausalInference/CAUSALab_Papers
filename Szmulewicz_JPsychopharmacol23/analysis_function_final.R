

analysis_function <- function(.dat, .formula_y_d, .formula_t_d, .formula_t_n=NULL, consecutive=1, bootstrap, bs_num=NULL){
  library(data.table)
  setDT(.dat)
  
  if(bootstrap){
    set.seed(bs_num)
    bs_ids <- table(sample(unique(.dat$SubjectNo), size=length(unique(.dat$SubjectNo)), replace=TRUE)) %>% as.data.frame(., stringsAsFactors=FALSE)
    bs_ids$SubjectNo <- as.numeric(bs_ids$Var1)
    bs_ids$Var1 <- NULL
    .dat <- left_join(bs_ids, .dat, by="SubjectNo")
    print("bootstrap resampling done")
  } else if (!bootstrap){
    Freq=1
  }
  
  tv <- ifelse(is.null(.formula_t_d),FALSE,TRUE) # automatically time-varying analysis if .formula_t_d non-null
  if(tv){
    if(!is.null(.formula_t_n)){
    t_n <- glm(formula = as.formula(.formula_t_n), 
               data = .dat %>% filter(allowed_disc==0 & clinical_visit==1 & cumvisit>consecutive & forb_medic==0),
               weights = .dat %>% filter(allowed_disc==0 & clinical_visit==1 & cumvisit>consecutive & forb_medic==0) %>% .$Freq,
               family = quasibinomial())
    } else if(is.null(.formula_t_n)){.dat$pred_t_n <- 1}
    
    t_d <- glm(formula = as.formula(.formula_t_d),
               data = .dat %>% filter(allowed_disc==0 & clinical_visit==1 & cumvisit>consecutive & forb_medic==0),
               weights = .dat  %>% filter(allowed_disc==0 & clinical_visit==1 & cumvisit>consecutive & forb_medic==0) %>% .$Freq,
               family = quasibinomial())
    .dat$pred_t_n <- predict(t_n, newdata=.dat, type="response")
    .dat$pred_t_d <- predict(t_d, newdata=.dat, type="response")
    
    .dat$w0 <- 1/.dat$pred_t_d
    .dat$s_w0 <- .dat$pred_t_n/.dat$pred_t_d
    
    # weight models
    .dat <- .dat %>% mutate(w0 = case_when (clinical_visit==1 & unscheduled==0 & allowed_disc==0 & cumvisit>consecutive & non_adherence==0 & forb_medic==0 & visit.week>1 ~ w0, #added nonadherence & visit.week here instead ofas extra step
                                            visit.week==0 ~ 1, 
                                            non_adherence==1 ~ 0,
                                            TRUE ~ 1),
                            s_w0 = case_when (clinical_visit==1 & unscheduled==0 & allowed_disc==0 & cumvisit>consecutive & non_adherence==0 & forb_medic==0 & visit.week>1 ~ s_w0,
                                              visit.week==0 ~ 1, 
                                              non_adherence==1 ~ 0,
                                              TRUE ~ 1))
    .dat$w <- with(.dat, ave(w0, SubjectNo, FUN=cumprod))
    .dat$s_w <- with(.dat, ave(s_w0, SubjectNo, FUN=cumprod))
    print("treatment model done")
  } else if(!tv){
    .dat$w <- 1
    .dat$s_w <- 1
  }
  
  m1 <- glm(formula = as.formula(.formula_y_d), ## formula must include "PO==0"
            data = .dat, 
            weights = s_w*Freq,
            family=quasibinomial())
  print("outcome model done")
  
  # dataset for outcome prediction
  data_for_risk_estimates_0a <- .dat%>%filter(visit.week==0) %>% select(-c(visit.week, visit.week2))
  data_for_risk_estimates_0b <- data.frame(SubjectNo=rep(unique(.dat$SubjectNo), each=53),
                                           visit.week=rep(0:52, times=length(unique(.dat$SubjectNo))), stringsAsFactors = FALSE)
  data_for_risk_estimates_0b$visit.week2 <- data_for_risk_estimates_0b$visit.week^2
  data_for_risk_estimates <- left_join(data_for_risk_estimates_0b, data_for_risk_estimates_0a, by="SubjectNo")
  
  data_for_risk_estimates$p_0 <- predict(object = m1, newdata=data_for_risk_estimates %>% mutate(lithium=0), 
                      type="response")
  data_for_risk_estimates$p_1 <- predict(object = m1, newdata=data_for_risk_estimates %>% mutate(lithium=1), 
                      type="response")
  .results <- data_for_risk_estimates %>% group_by(SubjectNo) %>% arrange(visit.week) %>% 
    mutate(risk_0 = 1 - cumprod(p_0),
           risk_1 = 1 - cumprod(p_1)) %>% 
    ungroup() %>% 
    filter(visit.week==52) %>% 
    dplyr::select(risk_0, risk_1) %>% 
    apply(., 2, mean) %>% t() %>% as.data.frame() %>% 
    mutate(rd = risk_1-risk_0,
           rr = risk_1/risk_0)
  return(.results)
  
}

