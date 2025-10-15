  source("scripts/load_globals.R")
  source(here::here("R", "process_data.R"))
  file_stub <- c("_iuiart", "_nonuse")[2]
  
  # specify file stub for outputting results
  if (file_stub %in% c("_nonuse")) {
    nonuse_stub <- "_nonuse"
  } else {
    nonuse_stub <- ""
  }
  
  # load data
  multiples <- read_parquet(glue("{JEREMY_DATA}/new{nonuse_stub}_cohort.parquet")) %>% process_data()
  
  # remove IUI users if comparison is non-use
  if (file_stub %in% c("_nonuse")) {
    multiples <- multiples %>% filter(any_iui == 0)
  }
  
  # add new definition of multiples
  multiples <- multiples %>% left_join(read_parquet(glue("{JEREMY_DATA}/redefined_multiples.parquet")))
  assert_that(all(!is.na(multiples$possible_multiple)))
  multiples <- multiples %>% mutate(pretty_multiples = if_else(probable_multiple == 1, "Multiple", "Singleton"))
  
  # create formatted exposure variable
  if (file_stub %in% c("_iuiart")) {
    multiples <- multiples %>% mutate(pretty_art = factor(if_else(art_only == 1, "ART", "IUI"), levels=c("IUI", "ART")))
  } else if (file_stub %in% "_nonuse") {
    multiples <- multiples %>% mutate(pretty_art = factor(if_else(art == 1, "ART", "Non-use"), levels=c("Non-use", "ART")))
  } else {
    stop("Unrecognised file_stub")
  }
  
  # create spreadsheet for outputs
  wb <- createWorkbook()
  set_wb_options(wb)
  
  # add no. linked infants
  multiples <- multiples %>% group_by(pregid) %>% mutate(two_linked=as.integer((n() > 1))) %>% ungroup()

  # create baseline table - maternal counts (TABLE 2)
  baseline_table_mat <- multiples %>% filter(inf_num == 1) %>% tbl_summary(by=pretty_art, include = c(age_cat, probable_multiple, possible_multiple, multiples_new, two_linked, tobacco, alcohol, drug_abuse, obesity, terat, antidep, 
                                                                                                      anticonv, opioids, codeine, antibiotic, epilepsy, hypertension, diabetes),
                                                                           missing_text="Missing", statistic = list(all_categorical() ~ "{n} ({scales::number(100 * as.numeric(gsub(',', '', n))/as.numeric(gsub(',', '', N)), accuracy = 0.1)})"),
                                                                           label = list(age_cat ~ "Age", probable_multiple ~ "Multiple birth", possible_multiple ~ "Possible multiple birth",
                                                                                        multiples_new ~ "Any multiple gestation", two_linked ~ "Linked to two infant IDs",
                                                                                        tobacco ~ "Tobacco use", alcohol ~ "Alcohol abuse or dependence",
                                                                                        obesity ~ "Obesity", terat ~ "Definite or suspected teratogen", antidep ~ "Antidepressants",
                                                                                        drug_abuse ~ "Drug abuse or dependence",
                                                                                        anticonv ~ "Anticonvulsants", opioids ~ "Opioids", codeine ~ "Codeine",
                                                                                        antibiotic ~ "Antibiotics", epilepsy ~ "Epilepsy", hypertension ~ "Hypertension",
                                                                                        diabetes ~ "Diabetes"))
  
  # save as word table and spreadsheet
  baseline_table_mat %>% as_flex_table() %>% save_as_docx(path=glue("output/baseline_table_mothers{file_stub}_{Sys.Date()}.docx"))
  addWorksheet(wb, sheetName = "Baseline - counts per mother", gridLines = TRUE)
  writeDataTableAuto(wb, sheet = 1, x = baseline_table_mat %>% as_tibble() %>% add_row(`**Characteristic**` = "Note: counts per mother"))
  
  # Create outcome table by multiple birth/no linked infants for ART group (TABLE 3)
  outcome_table_by_no_infants <- multiples %>% filter(pretty_art == "ART") %>%
    group_by(pregid) %>%  mutate(no_infants=n()) %>% ungroup() %>%
    mutate(group_type = case_when(probable_multiple == 1 & no_infants == 1 ~ "Multiple birth (1 linked infant)",
                             probable_multiple == 1 & no_infants == 2 ~ "Multiple birth (2 linked infants)",
                             probable_multiple == 0 ~ "Singleton")) %>%
    tbl_summary(by=group_type, include = c(malformation, malf_excl_chrom,
                                           cns, ear, eye, cardio, vas, resp, oral, gastro, genit, urinary, muscu, limb, other, nicu_infant),
                missing_text="Missing", statistic = list(all_categorical() ~ "{n} ({scales::number(100 * as.numeric(gsub(',', '', n))/as.numeric(gsub(',', '', N)), accuracy = 0.1)})"),
                label = list(malformation ~ "Any malformation", malf_excl_chrom ~ "Non-chromosomal malformation", 
                             cns ~ "Central nervous system", ear ~ "Ear", eye ~ "Eye", cardio ~ "Cardiac", vas ~ "Other vascular",
                             resp ~ "Respiratory", oral ~ "Orofacial", gastro ~ "Gastrointestinal", genit ~ "Genital", urinary ~ "Urinary",
                             muscu ~ "Muscular", limb ~ "Limb", other ~ "Other", nicu_infant ~ "NICU admission"))
  
  addWorksheet(wb, sheetName = "Outcomes by no infants", gridLines = TRUE)
  writeDataTableAuto(wb, sheet = 7, x = outcome_table_by_no_infants %>% as_tibble() %>% add_row(`**Characteristic**` = "Note: counts per infants"))
  outcome_table_by_no_infants %>% as_flex_table() %>% save_as_docx(path=glue("output/outcomes_by_infant{file_stub}.docx"))
  
  
  # Calculate intra-class correlation among ART group
  res <- NULL
  for (outcome in c("malformation", "malf_excl_chrom",
                    "cns", "ear", "eye", "cardio", "vas", "resp", "oral", "gastro", "genit", "urinary", "muscu", "limb", "other", "nicu_infant")) {
    icc_outcome <- multiples %>% 
      group_by(pregid) %>% filter(pretty_art == "ART") %>%
      mutate(no_infants=n()) %>% filter(no_infants == 2) %>%
      mutate(inf_num = row_number()) %>% select(pregid, inf_num, all_of(outcome)) %>% 
      pivot_wider(names_from=inf_num, values_from=all_of(outcome), id_cols=pregid) %>% ungroup() %>% select(!pregid) %>% icc()
    res <- bind_rows(res, tibble_row(outcome=outcome, icc=icc_outcome$value, conf.low=icc_outcome$lbound, conf.high=icc_outcome$ubound))
  }
  addWorksheet(wb, sheetName = "ICC by outcome", gridLines = TRUE)
  writeDataTableAuto(wb, sheet = 8, x = res)
  saveWorkbook(wb, glue("output/baseline{file_stub}_v2.xlsx"), overwrite = TRUE)
  
                             label = list(age_cat ~ "Age", probable_multiple ~ "Multiple birth", possible_multiple ~ "Possible multiple birth",
                                                                     multiples_new ~ "Any multiple gestation",two_linked ~ "Linked to two infant IDs",
                                                                     tobacco ~ "Tobacco use", alcohol ~ "Alcohol abuse or dependence",
                                                                     obesity ~ "Obesity", terat ~ "Definite or suspected teratogen", antidep ~ "Antidepressants",
                                                                     drug_abuse ~ "Drug abuse or dependence",
                                                                     anticonv ~ "Anticonvulsants", opioids ~ "Opioids", codeine ~ "Codeine",
                                                                     antibiotic ~ "Antibiotics", epilepsy ~ "Epilepsy", hypertension ~ "Hypertension",
                                                                     diabetes ~ "Diabetes")) 
  
  addWorksheet(wb, sheetName = "Baseline - per mom post-10-2015", gridLines = TRUE)
  writeDataTableAuto(wb, sheet = 3, x = baseline_table_mompost2015 %>% as_tibble() %>% add_row(`**Characteristic**` = "Note: counts per mother"))
  
  
  # create outcome table at maternal level
  outcome_table_by_treatment <- multiples %>% 
    group_by(pregid) %>%
    summarise(across(c(malformation, malf_excl_chrom,
                    cns, ear, eye, cardio, vas, resp, oral, gastro, genit, urinary, muscu, limb, other, nicu_infant), ~ max(., na.rm=TRUE)), pretty_art=first(pretty_art)) %>% 
    tbl_summary(by=pretty_art, include = c(malformation, malf_excl_chrom,
                                                                           cns, ear, eye, cardio, vas, resp, oral, gastro, genit, urinary, muscu, limb, other, nicu_infant),
                missing_text="Missing", statistic = list(all_categorical() ~ "{n} ({scales::number(100 * as.numeric(gsub(',', '', n))/as.numeric(gsub(',', '', N)), accuracy = 0.1)})"),
                                             label = list(malformation ~ "Any malformation", malf_excl_chrom ~ "Non-chromosomal malformation", 
                                                          cns ~ "Central nervous system", ear ~ "Ear", eye ~ "Eye", cardio ~ "Cardiac", vas ~ "Other vascular",
                                                          resp ~ "Respiratory", oral ~ "Orofacial", gastro ~ "Gastrointestinal", genit ~ "Genital", urinary ~ "Urinary",
                                                          muscu ~ "Muscular", limb ~ "Limb", other ~ "Other", nicu_infant ~ "NICU admission"))
  addWorksheet(wb, sheetName = "Outcomes by Rx", gridLines = TRUE)
  writeDataTableAuto(wb, sheet = 4, x = outcome_table_by_treatment %>% as_tibble() %>% add_row(`**Characteristic**` = "Note: counts per mother"))
  
  
  # create outcome table at maternal level pre Oct 2015
  outcome_table_by_treatment_pre2015 <- multiples %>% filter(enddate < dmy("01/10/2015")) %>%
    group_by(pregid) %>%
    summarise(across(c(malformation, malf_excl_chrom,
                       cns, ear, eye, cardio, vas, resp, oral, gastro, genit, urinary, muscu, limb, other, nicu_infant), ~ max(., na.rm=TRUE)), pretty_art=first(pretty_art)) %>% 
    tbl_summary(by=pretty_art, include = c(malformation, malf_excl_chrom,
                                           cns, ear, eye, cardio, vas, resp, oral, gastro, genit, urinary, muscu, limb, other, nicu_infant),
                missing_text="Missing", statistic = list(all_categorical() ~ "{n} ({scales::number(100 * as.numeric(gsub(',', '', n))/as.numeric(gsub(',', '', N)), accuracy = 0.1)})"),
                label = list(malformation ~ "Any malformation", malf_excl_chrom ~ "Non-chromosomal malformation", 
                             cns ~ "Central nervous system", ear ~ "Ear", eye ~ "Eye", cardio ~ "Cardiac", vas ~ "Other vascular",
                             resp ~ "Respiratory", oral ~ "Orofacial", gastro ~ "Gastrointestinal", genit ~ "Genital", urinary ~ "Urinary",
                             muscu ~ "Muscular", limb ~ "Limb", other ~ "Other", nicu_infant ~ "NICU admission"))
  addWorksheet(wb, sheetName = "Outcomes by Rx pre-10-15", gridLines = TRUE)
  writeDataTableAuto(wb, sheet = 5, x = outcome_table_by_treatment_pre2015 %>% as_tibble() %>% add_row(`**Characteristic**` = "Note: counts per mother"))
  
  # create outcome table at maternal level post Oct 2015
  outcome_table_by_treatment_post2015 <- multiples %>% filter(enddate >= dmy("01/10/2015")) %>%
    group_by(pregid) %>%
    summarise(across(c(malformation, malf_excl_chrom,
                       cns, ear, eye, cardio, vas, resp, oral, gastro, genit, urinary, muscu, limb, other, nicu_infant), ~ max(., na.rm=TRUE)), pretty_art=first(pretty_art)) %>% 
    tbl_summary(by=pretty_art, include = c(malformation, malf_excl_chrom,
                                           cns, ear, eye, cardio, vas, resp, oral, gastro, genit, urinary, muscu, limb, other, nicu_infant),
                missing_text="Missing", statistic = list(all_categorical() ~ "{n} ({scales::number(100 * as.numeric(gsub(',', '', n))/as.numeric(gsub(',', '', N)), accuracy = 0.1)})"),
                label = list(malformation ~ "Any malformation", malf_excl_chrom ~ "Non-chromosomal malformation", 
                             cns ~ "Central nervous system", ear ~ "Ear", eye ~ "Eye", cardio ~ "Cardiac", vas ~ "Other vascular",
                             resp ~ "Respiratory", oral ~ "Orofacial", gastro ~ "Gastrointestinal", genit ~ "Genital", urinary ~ "Urinary",
                             muscu ~ "Muscular", limb ~ "Limb", other ~ "Other", nicu_infant ~ "NICU admission"))
  addWorksheet(wb, sheetName = "Outcomes by Rx post-10-15", gridLines = TRUE)
  writeDataTableAuto(wb, sheet = 6, x = outcome_table_by_treatment_post2015 %>% as_tibble() %>% add_row(`**Characteristic**` = "Note: counts per mother"))
  
  # Create outcome table by multiple birth/no linked infants for ART group
  outcome_table_by_no_infants <- multiples %>% filter(pretty_art == "ART") %>%
    group_by(pregid) %>%  mutate(no_infants=n()) %>% ungroup() %>%
    mutate(group_type = case_when(probable_multiple == 1 & no_infants == 1 ~ "Multiple birth (1 linked infant)",
                             probable_multiple == 1 & no_infants == 2 ~ "Multiple birth (2 linked infants)",
                             probable_multiple == 0 ~ "Singleton")) %>%
    tbl_summary(by=group_type, include = c(malformation, malf_excl_chrom,
                                           cns, ear, eye, cardio, vas, resp, oral, gastro, genit, urinary, muscu, limb, other, nicu_infant),
                missing_text="Missing", statistic = list(all_categorical() ~ "{n} ({scales::number(100 * as.numeric(gsub(',', '', n))/as.numeric(gsub(',', '', N)), accuracy = 0.1)})"),
                label = list(malformation ~ "Any malformation", malf_excl_chrom ~ "Non-chromosomal malformation", 
                             cns ~ "Central nervous system", ear ~ "Ear", eye ~ "Eye", cardio ~ "Cardiac", vas ~ "Other vascular",
                             resp ~ "Respiratory", oral ~ "Orofacial", gastro ~ "Gastrointestinal", genit ~ "Genital", urinary ~ "Urinary",
                             muscu ~ "Muscular", limb ~ "Limb", other ~ "Other", nicu_infant ~ "NICU admission"))
  
  addWorksheet(wb, sheetName = "Outcomes by no infants", gridLines = TRUE)
  writeDataTableAuto(wb, sheet = 7, x = outcome_table_by_no_infants %>% as_tibble() %>% add_row(`**Characteristic**` = "Note: counts per infants"))
  outcome_table_by_no_infants %>% as_flex_table() %>% save_as_docx(path=glue("output/outcomes_by_infant{file_stub}.docx"))
  
  
  # Calculate intra-class correlation among ART group
  res <- NULL
  for (outcome in c("malformation", "malf_excl_chrom",
                    "cns", "ear", "eye", "cardio", "vas", "resp", "oral", "gastro", "genit", "urinary", "muscu", "limb", "other", "nicu_infant")) {
    icc_outcome <- multiples %>% 
      group_by(pregid) %>% filter(pretty_art == "ART") %>%
      mutate(no_infants=n()) %>% filter(no_infants == 2) %>%
      mutate(inf_num = row_number()) %>% select(pregid, inf_num, all_of(outcome)) %>% 
      pivot_wider(names_from=inf_num, values_from=all_of(outcome), id_cols=pregid) %>% ungroup() %>% select(!pregid) %>% icc()
    res <- bind_rows(res, tibble_row(outcome=outcome, icc=icc_outcome$value, conf.low=icc_outcome$lbound, conf.high=icc_outcome$ubound))
  }
  addWorksheet(wb, sheetName = "ICC by outcome", gridLines = TRUE)
  writeDataTableAuto(wb, sheet = 8, x = res)
  saveWorkbook(wb, glue("output/baseline{file_stub}_v2.xlsx"), overwrite = TRUE)
  
