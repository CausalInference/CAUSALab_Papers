################################################################
# CAUSALAB
# /n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_analyses_clean.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-24
# 1) Purpose: Make functions to clean data for the following g-formula

extract_prebase_base <- function(dta) {
  dta_prebaseline <- dta %>% 
    filter(period == 0) %>% # lag-1 term at period 0 is the pre-baseline value
    mutate(
           pre_act = act_l1,
           pre_totfv = totfv_l1,
           pre_totgrn = totgrn_l1, 
           pre_junk = junk_l1,
           pre_totred = totred_l1,
           pre_totproc = totproc_l1,
           pre_sodajuice = sodajuice_l1,
           pre_totalc = totalc_l1,
           pre_cal = cal_l1,
           pre_hbp = hbp_l1,
           pre_chl = chl_l1,
           pre_dia = dia_l1,
           pre_xcond2 = xcond2_l1,
           pre_xcond3 = xcond3_l1,
           pre_bmi = bmi_l1,
           pre_hlthy = hlthy_l1,
           hlthy_base = hlthy,
           age_base = age,
           base_dia = dia,
           base_hbp = hbp,
           base_chl = chl) %>%
    select(id,
           pre_act, 
           pre_totfv, 
           pre_totgrn, 
           pre_junk,
           pre_totred,
           pre_totproc,
           pre_sodajuice,
           pre_totalc,
           pre_cal,
           pre_hbp,
           pre_chl,
           pre_dia,
           pre_xcond2,
           pre_xcond3,
           pre_bmi,
           pre_hlthy,
           hlthy_base,
           age_base,
           base_dia,
           base_hbp,
           base_chl) 
  
  # merge the baseline information to every row of data
  dta <- merge(dta, dta_prebaseline, by = "id") %>%
    select(id, period, event, censor, dead_other,
           xcond2, xcond3, pca_Y,
           age_base, base_dia, base_hbp, base_chl,fhxmi, fhxpca90, hlthy_base,
           act, totfv, totgrn, junk, totred, totproc, sodajuice, totalc, cal,
           hbp, chl, dia, bmi, smkhx, psa,
           pre_act, pre_totfv, pre_totgrn, pre_junk, pre_totred, pre_totproc, pre_sodajuice,
           pre_totalc, pre_cal, pre_hbp, pre_chl, pre_dia, pre_xcond2, pre_xcond3, 
           pre_bmi, pre_hlthy)%>%
    arrange(id, period) 
  
  return(dta)
}


# log transform energy intake and BMI
trans_var <- function(dta) {
  dta %>%
    mutate(
      cal= log(cal),
      bmi = log(bmi)) %>%
    return()
}


###########################
# make a new function for coarsened time period
coarsen <- function(x){
  return( ifelse(x<2, 2, floor(x / 2) + 1))
}

