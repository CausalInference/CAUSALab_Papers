# Generate simulate data
simulate_data <- function(sample_size, params) {
  # simulate data
  sim_data <- tibble(
    pregid = 1:sample_size,
    a = rbinom(sample_size, 1, params$treatment_prevalence),
    l = rbinom(sample_size, 1, params$covariate_prevalence),
    m = rbinom(sample_size, 1, exp(log(params$mediator_prevalence) + log(params$mediator_treatment_rr)*a + log(params$mediator_covariate_rr)*l)),
    mean_y = exp(log(params$outcome_prevalence) + log(params$outcome_treatment_rr)*a + log(params$outcome_mediator_rr)*m + log(params$outcome_covariate_rr)*l),
    y = rbinom(sample_size, 1, mean_y)
  )
  
  # duplicate twins
  sim_data <- sim_data %>% bind_rows(sim_data %>% filter(m == 1))
  
  # separate twins and singletons
  sim_data_singletons <- sim_data %>% filter(m==0)
  sim_data_twins <- sim_data %>% filter(m==1)
  
  # determine outcome in second twin using method of Qaqish
  sim_data_twins <- sim_data_twins %>% 
    group_by(pregid) %>%
    mutate(twin_id = row_number()) %>% 
    ungroup() 
  
  sim_data_first_twin <- sim_data_twins %>% filter(twin_id == 1)
  sim_data_second_twin <- sim_data_twins %>% filter(twin_id == 2)
  
  sim_data_second_twin <- sim_data_second_twin %>% 
    mutate(mean_y = ((1 - params$correlation)*mean_y + params$correlation*y)) %>%
    mutate(y = rbinom(nrow(sim_data_second_twin), 1, mean_y))
  
  sim_data <- bind_rows(sim_data_first_twin, sim_data_second_twin, sim_data_singletons) %>% select(!c(mean_y, twin_id))
  
  return(sim_data)
}

# calculate performance measures for all iterations in one scenario
calculate_performance_measures <- function(sim_estimates) {
  sim_estimates <- sim_estimates %>% summarise(
    across(starts_with("estimate"), ~ mean(.x, na.rm=TRUE), .names="mean_{.col}"), 
    across(starts_with("estimate"), ~ sqrt(var(.x, na.rm=TRUE)), .names="empSE_{.col}"),
    across(starts_with("std.error"), ~ sqrt(mean(.x^2)), .names="avg.Mod.{.col}"),
    across(starts_with("estimate"), ~ sum(as.integer(is.na(.x))), .names="no_missing_{.col}"),
    across(starts_with("risk"), ~ mean(.x, na.rm=TRUE), .names="mean_{.col}"), 
    across(starts_with("risk"), ~ sqrt(var(.x, na.rm=TRUE)), .names="empSE_{.col}"),
  ) 
  
  estimate_names <- str_extract(colnames(sim_estimates)[str_detect(colnames(sim_estimates), "^empSE(?!_risk)")], "(?<=empSE_estimate_).+")
  for (estimate_name in estimate_names) {
    sim_estimates[glue("relative_error_SE_{estimate_name}")] <- 100 * sim_estimates[glue("avg.Mod.std.error_{estimate_name}")] / sim_estimates[glue("empSE_estimate_{estimate_name}")]
  }

  return(sim_estimates)
}


# analyse simulated data
analyse_data <- function(sim_data) {
  
  outputs <- NULL
  
  # naive analysis ----------------------------------------------------------
   
  naive_model <- glm(y ~ a, data=sim_data, family="poisson") 
    outputs <- naive_model %>%
      broom::tidy() %>% 
      filter(term == "a") %>% 
      select(estimate) %>%
      mutate(std.error = coeftest(naive_model, vcov = vcovHC(naive_model, type="HC1"))["a", "Std. Error"]) %>% 
      mutate(method="naive_analysis") %>%
      bind_rows(outputs)
  
  # restriction to singletons -----------------------------------------------
   singleton_model <- glm(y ~ a, data=sim_data %>% filter(m==0), family="poisson") 
    outputs <- singleton_model %>%
      broom::tidy() %>% 
      filter(term == "a") %>% 
      select(estimate) %>%
      mutate(std.error = coeftest(singleton_model, vcov = vcovHC(singleton_model, type="HC1"))["a", "Std. Error"]) %>% 
      mutate(method="restriction_singletons") %>% 
      bind_rows(outputs)
  
  # analysis at the pregnancy-level -----------------------------------------
  pregnancy_level_data <- sim_data %>% 
    group_by(pregid) %>% 
    summarise(y = max(y), a=first(a)) %>% 
    ungroup()
  
  plevel_model <- glm(y ~ a, data=pregnancy_level_data, family="poisson") 
  
  outputs <- plevel_model %>%
    broom::tidy() %>% 
    filter(term == "a") %>% 
    select(estimate) %>%
    mutate(std.error = coeftest(plevel_model, vcov = vcovHC(plevel_model, type="HC1"))["a", "Std. Error"]) %>% 
    mutate(method="pregnancy_level") %>% 
    bind_rows(outputs)
  
  # gee analysis exchangeable -----------------------------------------------
  gee_estimate <- gee(y ~ a, data=sim_data %>% arrange(pregid), id=pregid, family="poisson", corstr="exchangeable") 
  outputs <- tibble_row(estimate = gee_estimate$coefficients[[2]], std.error=sqrt(gee_estimate$robust.variance[2,2])) %>% 
    mutate(method="gee_exchangeable") %>% bind_rows(outputs)
  
  # gee analysis independence -----------------------------------------------
  gee_estimate <- gee(y ~ a, data=sim_data %>% arrange(pregid), id=pregid, family="poisson", corstr="independence") 
  outputs <- tibble_row(estimate = gee_estimate$coefficients[[2]], std.error=sqrt(gee_estimate$robust.variance[2,2])) %>% 
    mutate(method="gee_independence") %>% bind_rows(outputs)
  
  return(outputs)
}