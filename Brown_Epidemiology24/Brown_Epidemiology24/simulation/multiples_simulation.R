library(here)
source(here::here("scripts", "load_globals.R"))

# specify parameters
no_pregnancies <- 10000L
no_iterations <- 1000L

# set random seed
set.seed(13290, "L'Ecuyer")
mc.reset.stream()

# master function to run simulation
run_simulation <- function(no_pregnancies, params) {
  sim_data <- simulate_data(no_pregnancies, params)
  rr_results <- analyse_data(sim_data) %>% 
    pivot_wider(names_from=method, values_from=c(estimate, std.error))
  return(rr_results)
}

# load scenarios
scenarios <- expand_grid(treatment_prevalence=0.5, covariate_prevalence=0.1, 
                         mediator_prevalence=0.02, mediator_treatment_rr=seq(1,10,0.5),
                         mediator_covariate_rr=1, outcome_prevalence=c(0.03, 0.1), 
                         outcome_treatment_rr=c(1,2), outcome_mediator_rr=c(1,2), 
                         outcome_covariate_rr=1, correlation=c(0, 0.5, 1)) %>%
  filter(!(outcome_prevalence == 0.1 & correlation > 0)) %>% 
  mutate(scenario_no = row_number())
  
# practice example
sim_data <- simulate_data(no_pregnancies, scenarios[1,])
# run_simulation(no_pregnancies, scenarios[1,])


# run simulations 
outputs <- NULL
for (i in 1:nrow(scenarios)) {
  print(i)
  outputs <- mclapply(1:no_iterations, function(x) run_simulation(no_pregnancies, scenarios[i,]), 
                      mc.cores=8, mc.set.seed = TRUE) %>% 
    bind_rows() %>% 
    calculate_performance_measures() %>%
    mutate(scenario_no = i) %>%
    bind_rows(outputs) 
}

# join outputs with scenarios
outputs <- inner_join(scenarios, outputs, by="scenario_no")

# write outputs
outputs %>% write_csv(here::here("outputs", glue("outputs_{Sys.Date()}.csv")))



