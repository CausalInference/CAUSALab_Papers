source("scripts/load_globals.R")


# Restrict to single gestation -----------------------------------------------------------------------------------------------------

models <- list()
multiple_stub <- "multiples_new"

# load models
for (comparison in c("_artiui_wcalyear", "_nonuse_wcalyear")) {
  for (outcome in c("malformation", "nicu_infant")) {
    if (comparison %in% c("_artiui_wcalyear")) {
      exposure <- "art_only"
    } else {
      exposure <- "art"
    }

    models[[glue("singleton_crude_{outcome}{comparison}")]] <- readRDS(file=glue("models/singleton_crude_pois_model_{multiple_stub}_{outcome}{comparison}.rds")) %>% 
      get_robust_estimate(exposure) %>% mutate(type = "Crude", comparison=comparison, outcome=outcome)
    models[[glue("singleton_adjusted_{outcome}{comparison}")]] <- readRDS(file=glue("models/singleton_adj_pois_model_{multiple_stub}_{outcome}{comparison}.rds")) %>% 
      get_robust_estimate(exposure) %>% mutate(type = "Adjusted", comparison=comparison, outcome=outcome)
  }
}

# process model outputs
model_outputs <- models %>% bind_rows() %>% 
  mutate(across(c(mean, lower, upper), ~ exp(.))) %>%
  mutate(rr = as.character(glue("{scales::number(mean, accuracy=0.01)} ({scales::number(lower, accuracy=0.01)}-{scales::number(upper, accuracy=0.01)})")))  %>%
  mutate(comparison = case_when(comparison == "_artiui_wcalyear" ~ "IUI", comparison == "_nonuse_wcalyear" ~ "Non-use", .default="Error"),
         outcome = case_when(outcome == "malformation" ~ "Malformation", outcome == "nicu_infant" ~ "NICU admission"))


# save outputs to latex table
output_table <- model_outputs %>% pivot_wider(id_cols=c(comparison, outcome), names_from=c(type), values_from=rr)
print(xtable(as.data.frame(output_table), type = "latex"), file = "output/output_table.tex")