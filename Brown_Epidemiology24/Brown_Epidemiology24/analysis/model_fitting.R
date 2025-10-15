source("scripts/load_globals.R")
source(here::here("R", "process_data.R"))

# specify parameters if running from script
if (sys.nframe() == 0) {
  comparison <- c("_artiui", "_nonuse","_artiui_wcalyear", "_nonuse_wcalyear")[4]
  outcome <- c("malformation","nicu_infant", "male")[1]
  period <- c("preoct15", "postoct15", "")[3]
  exposure <- c("art_only", "art")[2]
}

# specify non-use stub
if (comparison %in% c("_nonuse", "_nonuse_wcalyear")) {
  nonuse_stub <- "_nonuse"
} else {
  nonuse_stub <- ""
}

# load data
multiples <- read_parquet(glue("{JEREMY_DATA}/new{nonuse_stub}_cohort.parquet")) %>% process_data()

if (comparison %in% c("_nonuse", "_nonuse_wcalyear")) {
  multiples <- multiples %>% filter(any_iui == 0)
}

# add new definition of multiples
multiples <- multiples %>% left_join(read_parquet(glue("{JEREMY_DATA}/redefined_multiples.parquet")))
assert_that(all(!is.na(multiples$possible_multiple)))
assert_that(all(!is.na(multiples$probable_multiple)))

# restrict on period
if (period == "preoct15") {
 multiples <- multiples %>% filter(enddate < dmy("01/10/2015")) 
 file_stub <- paste0(comparison, "_preoct15")
} else if (period == "postoct15") {
  multiples <- multiples %>% filter(enddate >= dmy("01/10/2015")) 
  file_stub <- paste0(comparison, "_postoct15")
} else {
  file_stub <- comparison
}


# remove variables not needed for analysis (helps with datadist)
analysis_data <- multiples %>% select(pregid, delivery_year, multiples_new, probable_multiple, possible_multiple, nicu_infant, malformation, malf_excl_chrom, cns, eye, ear, cardio, genit, vas, resp, oral, gastro, urinary, muscu, limb,
                                      art_only, art, age_cat, AGE, obesity, tobacco, alcohol, drug_abuse,
                                      diabetes, hypertension, epilepsy, terat, antidep, anticonv, opioids, antibiotic, SEX) %>%
  mutate(delivery_year = factor(as.character(delivery_year)), male=as.integer(SEX == 1))


gformula_formula <- as.formula(glue("{outcome} ~ {exposure} + AGE + I(AGE^2) + delivery_year + obesity + tobacco + alcohol + drug_abuse + diabetes + hypertension + epilepsy + terat + antidep + anticonv + opioids + antibiotic"))
ps_formula <- as.formula(glue("{exposure} ~ AGE + I(AGE^2) + delivery_year + obesity + tobacco + alcohol + drug_abuse + diabetes + hypertension + epilepsy + terat + antidep + anticonv + opioids + antibiotic"))
selection_formula <- as.formula(glue("selected ~ {exposure} + (AGE + I(AGE^2) + delivery_year)^2 + obesity + tobacco + alcohol + drug_abuse + diabetes + hypertension + epilepsy + terat + antidep + anticonv + opioids + antibiotic"))

# Restriction to singletons ---------------------------------------------------------------------------------------

singleton_analysis <- function(analysis_data, multiple_var, outcome, exposure, file_stub, gformula_formula) {
  # prepare data
  singleton_data <- analysis_data[analysis_data[[multiple_var]] == 0,]
  ddist <- datadist(singleton_data)
  options("datadist" = "ddist")
  
  # fit poisson models
  singleton_crude_pois_model <- glm(as.formula(glue("{outcome} ~ {exposure}")), data=singleton_data, 
                                    family="poisson")
  singleton_adj_pois_model <- glm(as.formula(glue("{outcome} ~ {exposure} + AGE + I(AGE^2)  + delivery_year + obesity + tobacco + alcohol + drug_abuse + diabetes + hypertension + epilepsy + terat + antidep + anticonv + opioids + antibiotic")), data=singleton_data, 
                                  family="poisson")
  saveRDS(singleton_crude_pois_model, file=glue("models/singleton_crude_pois_model_{multiple_var}_{outcome}{file_stub}.rds"))
  saveRDS(singleton_adj_pois_model, file=glue("models/singleton_adj_pois_model_{multiple_var}_{outcome}{file_stub}.rds"))
  
  # estimate risk difference
  singleton_rd <- calc_risk_estimates_delta(singleton_data, exposure, outcome, gformula_formula)
  write_csv(singleton_rd, glue("output/singleton_rd_{multiple_var}_{outcome}{file_stub}.csv"))
}

# overall
singleton_analysis(analysis_data, "probable_multiple", outcome, exposure, file_stub, gformula_formula)
# possible multiple
singleton_analysis(analysis_data, "possible_multiple", outcome, exposure, file_stub, gformula_formula)
# any multiple
singleton_analysis(analysis_data, "multiples_new", outcome, exposure, file_stub, gformula_formula)


# Count at pregnancy-level ----------------------------------------------------------------------------------------

plevel_analysis <- function(analysis_data, outcome, exposure, file_stub, gformula_formula) {
  plevel_data <- analysis_data %>% group_by(pregid) %>% 
    mutate(infant_no = row_number(), across(c(malformation, malf_excl_chrom,
                                              cns, ear, eye, cardio, vas, resp, oral, gastro, genit, urinary, muscu, limb), ~ max(., na.rm=TRUE))) %>% 
    ungroup() %>% filter(infant_no == 1)
  ddist <- datadist(plevel_data)
  options("datadist" = "ddist")
  
  # fit poisson models
  plevel_crude_pois_model <- glm(as.formula(glue("{outcome} ~ {exposure}")), data=plevel_data, 
                                 family="poisson")
  plevel_adj_pois_model <- glm(as.formula(glue("{outcome} ~ {exposure} + AGE + I(AGE^2) + delivery_year + obesity + tobacco + alcohol + drug_abuse + diabetes + hypertension + epilepsy + terat + antidep + anticonv + opioids + antibiotic")), data=plevel_data, 
                               family="poisson")
  saveRDS(plevel_crude_pois_model, file=glue("models/plevel_crude_pois_model_{outcome}{file_stub}.rds"))
  saveRDS(plevel_adj_pois_model, file=glue("models/plevel_adj_pois_model_{outcome}{file_stub}.rds"))
  
  # estimate risk difference
  plevel_rd <- calc_risk_estimates_delta(plevel_data, exposure, outcome, gformula_formula)
  write_csv(plevel_rd, glue("output/plevel_rd_{outcome}{file_stub}.csv"))
}


plevel_analysis(analysis_data, outcome, exposure, file_stub, gformula_formula)

# Restricted to complete data no reweighting----------------------------------------------------------------------------------------

restricted_analysis_no_ipcw <- function(analysis_data, outcome, exposure, file_stub) {
  
  restricted_data <- analysis_data %>% group_by(pregid) %>% mutate(no_infants=n()) %>% 
    ungroup() %>% mutate(selected = as.integer(no_infants > 1 | probable_multiple == 0))
  
  # prepare data
  restricted_data <- restricted_data %>% filter(selected == 1)

  # prepare data
  ddist <- datadist(restricted_data)
  options("datadist" = "ddist")
  
  # fit GEE poisson models
  restricted_no_ipcw_crude_pois_model <- geeglm(as.formula(glue("{outcome} ~ {exposure}")), 
                                        data=restricted_data, id=pregid, family="poisson", corstr="independence")
  restricted_no_ipcw_adj_pois_model <- geeglm(as.formula(glue("{outcome} ~ {exposure} + AGE + I(AGE^2) + delivery_year + obesity + tobacco + alcohol + drug_abuse + diabetes + hypertension + epilepsy + terat + antidep + anticonv + opioids + antibiotic")), 
                                      data=restricted_data, id=pregid, family="poisson", corstr="independence")
  
  saveRDS(restricted_no_ipcw_crude_pois_model, file=glue("models/restricted_no_ipcw_crude_pois_model_{outcome}{file_stub}.rds"))
  saveRDS(restricted_no_ipcw_adj_pois_model, file=glue("models/restricted_no_ipcw_adj_pois_model_{outcome}{file_stub}.rds"))
  
}

restricted_analysis_no_ipcw(analysis_data, outcome, exposure, file_stub)

# Restricted to complete data with reweighting----------------------------------------------------------------------------------------

restricted_analysis <- function(analysis_data, outcome, exposure, file_stub) {
  
  restricted_data <- analysis_data %>% group_by(pregid) %>% mutate(no_infants=n()) %>% 
    ungroup() %>% mutate(selected = as.integer(no_infants > 1 | probable_multiple == 0))
  
  # fit selection model at maternal level amongst twins
  selection_model <- glm(selection_formula, data= restricted_data %>% group_by(pregid)  %>% filter(row_number()==1 & probable_multiple == 1) %>% ungroup(), family="binomial") 
  
  # prepare data
  restricted_data <- restricted_data %>% filter(selected == 1)
  
  # add weights
  restricted_data[["probselect"]] <- predict(selection_model, newdata=restricted_data, type="response")
  restricted_data <- restricted_data %>% mutate(ipcw = (1-probable_multiple)*1 + probable_multiple*(1/(probselect)))
  
  # prepare data
  ddist <- datadist(restricted_data)
  options("datadist" = "ddist")
  
  # fit GEE poisson models
  restricted_crude_pois_model <- geeglm(as.formula(glue("{outcome} ~ {exposure}")), 
                                      data=restricted_data, id=pregid, family="poisson", weights=ipcw, corstr="independence")
  restricted_adj_pois_model <- geeglm(as.formula(glue("{outcome} ~ {exposure} + AGE + I(AGE^2) + delivery_year + obesity + tobacco + alcohol + drug_abuse + diabetes + hypertension + epilepsy + terat + antidep + anticonv + opioids + antibiotic")), 
                                        data=restricted_data, id=pregid, family="poisson", weights=ipcw, corstr="independence")
  
  saveRDS(restricted_crude_pois_model, file=glue("models/restricted_crude_pois_model_{outcome}{file_stub}.rds"))
  saveRDS(restricted_adj_pois_model, file=glue("models/restricted_adj_pois_model_{outcome}{file_stub}.rds"))
  
}

restricted_analysis(analysis_data, outcome, exposure, file_stub)
