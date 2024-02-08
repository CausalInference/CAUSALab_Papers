/* 

This script is a part of a set of analytic programs related to the following manuscript(s):
- Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
- Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors.  

Author(s): VA-CAUSAL Methods Core
Version: July 2023. 

*/

/*Table 2, 3 and supplemental results by overall an subgroup analyses*/

proc iml;

submit / R;
source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/preamble.R") # change appropriately to P drive when storing results there
preamble(server=TRUE) # if using P drive, set this to "FALSE"
source(paste0(prefix,"/functions/bs-results-summary_function.R"))
source(paste0(prefix,"/functions/risk-point-ests_function.R"))
source(paste0(prefix,"/functions/rd_rr_function.R"))
source(paste0(prefix,"/functions/table2-function.R"))

#Aggregate results from all outcomes for a given scenario
delom_tab2 <- parallel_table_2_function(treat="PvM", subgroup=NULL, variant="delom_booster", 
                          risk_time=(16*7))
omicron_tab2 <- parallel_table_2_function(treat="PvM", subgroup=NULL, variant="omicron_boost", 
                          risk_time=(9*7), outlist=c("covidpos"))

#Aggregate results from specific outcomes for a given scenario

## DELTA-OMICRON Combined
# old70 - all
old70_delom <- parallel_table_2_function(treat="PvM", subgroup="old70", variant="delom_booster", 
                          risk_time=(16*7))

# young70 - all, except covid death
young70_delom <- parallel_table_2_function(treat="PvM", subgroup="young70", variant="delom_booster", 
                          risk_time=(16*7),
							outlist=c("covidpos","covidpossymp","covidhosp","covidICU"))

# white - all
white_delom <- parallel_table_2_function(treat="PvM", subgroup="white", variant="delom_booster", 
                          risk_time=(16*7))

# black - covidpos, covidpossymp, covidhosp 
black_delom <- parallel_table_2_function(treat="PvM", subgroup="black", variant="delom_booster", 
                          risk_time=(16*7), 
                          outlist=c("covidpos","covidpossymp","covidhosp"))

# diff23_67 - covidpos, covidpossymp, covidhosp 
diff23_67_delom <- parallel_table_2_function(treat="PvM", subgroup="diff23_67", variant="delom_booster", 
                          risk_time=(16*7), 
                          outlist=c("covidpos","covidpossymp","covidhosp"))

# diff23_8 - all, except covid death 
diff23_8_delom <- parallel_table_2_function(treat="PvM", subgroup="diff23_8", variant="delom_booster", 
                          risk_time=(16*7),
							outlist=c("covidpos","covidpossymp","covidhosp","covidICU"))

# diff23_9 - all, except covid death
diff23_9_delom <- parallel_table_2_function(treat="PvM", subgroup="diff23_9", variant="delom_booster", 
                          risk_time=(16*7),
							outlist=c("covidpos","covidpossymp","covidhosp","covidICU"))

# nopriorcov - all
noprior_delom <- parallel_table_2_function(treat="PvM", subgroup="nopriorcov", variant="delom_booster", 
                          risk_time=(16*7))

# trueboost - all
trueboost_delom <- parallel_table_2_function(treat="PvM", subgroup="trueboost", variant="delom_booster", 
                          risk_time=(16*7))

# prim_pfizer - covidpos
pfizer_delom <- table_2_function(treat="PvM", subgroup="prim_pfizer", variant="delom_booster", 
                          risk_time=(16*7), 
                          event="covidpos")

# prim_moderna - covidpos 
moderna_delom <- parallel_table_2_function(treat="PvM", subgroup="prim_moderna", variant="delom_booster", 
                          risk_time=(16*7), 
                          outlist=c("covidpos"))

final_table2 <- rbind(delom_tab2, young70_delom, old70_delom, white_delom, black_delom, 
						diff23_67_delom, diff23_8_delom, diff23_9_delom,
						noprior_delom, 
						trueboost_delom, pfizer_delom, moderna_delom,
						omicron_tab2)

write.csv(final_table2,
		"[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/AllBootstrappedResults.csv",
		row.names=FALSE)

endsubmit;
run;

*** Everything related to the parts 5 and 6 tables and figures;

/*************************************************************************************************/

libname sasin '[SAS_folder_ORD_Project]/Vax_Standard/analysis/data';

/*proc contents data=sasin.matched_PM_omboost; run;*/
/*proc contents data=sasin.matched_PM_delboost; run;*/
/*proc contents data=sasin.matched_PM_delom; run;*/

*********************************************;
** MATCHED TABLES 	;
*********************************************;

*********************************************;
** DELTA-OMICRON COMBINED PERIOD	;
*********************************************;

*This is used for the comparison of test counts in Table 2;
proc means data=sasin.matched_PM_delom mean median q1 q3 range sum;
	class dose3_man;
	var testcnt_post;
run;

proc iml;

run ExportDataSettoR("sasin.matched_PM_delom","matched_PM_merged");
submit / R;

pacman::p_load(MatchIt,ggQC,cobalt,ggplot2,tidyr,viridis,extrafont,dplyr,survival,survminer,tableone)


names(matched_PM_merged) <- tolower(names(matched_PM_merged))

myVars = c("age_at_index","age_cat","male","female","race","ethnicity",
"urbanicity",
"smoking_status","cond_lung","cond_vasc","cond_hypt","cond_diab","cond_ckd",
"cond_cld","cond_5yr_cancer","immunsupp",
"obesity","cond_dement","cond_sud",
"diff23_3grp","BNT162b2","mRNA1273","ntests","pcp5","flu5",
"dose2_month","covidpriordose3","nopriorcov","same_primary","dose3_booster","tot_fup","testcnt_post")

catVars=c("age_cat","male","female","race","ethnicity",
"urbanicity",
"smoking_status","cond_lung","cond_vasc","cond_hypt","cond_diab","cond_ckd",
"cond_cld","cond_5yr_cancer","immunsupp",
"obesity","cond_dement","cond_sud",
"diff23_3grp","BNT162b2","mRNA1273","ntests","pcp5","flu5",
"dose2_month","covidpriordose3","nopriorcov","same_primary","dose3_booster")


matched_PM_merged <- matched_PM_merged %>% 
					mutate(dose3 = factor(group_binary, levels=c(1,0), labels=c('BNT162b2','mRNA-1273'))) %>%
					mutate(current_smoker = case_when(smoking_status==1 ~ 1,
					TRUE ~ 0)) %>%
					mutate(same_primary = case_when(dose2_man==dose3_man ~ 1, 
					TRUE ~ 0)) %>% # variable denoting whether the primary manufacturer matches dose 3 (01mar2022 HG)
					#mutate(sex = factor(sex, levels=c(1,0), labels=c('Male','Female'))) %>%
					mutate(age_cat = factor(age_cat, levels = c(0,1,2,3,4,5),
											labels = c('18 to 39 years','40 to 49 years','50 to 59 years',
														'60 to 69 years', '70 to 79 years','>=80 years'))) %>%
					mutate(male = case_when(sex==1 ~ 1, TRUE ~ 0)) %>%
					mutate(female = case_when(sex==0 ~ 1, TRUE ~ 0)) %>%
					mutate(black = case_when(race==1 ~ 1, TRUE ~ 0)) %>%
					mutate(white = case_when(race==0 ~ 1, TRUE ~ 0)) %>%
					mutate(race = factor(race, levels=c(1,2,3,0), labels=c('Black','Other','Unknown','White'))) %>%
					mutate(ethnicity = factor(ethnicity, levels=c(1,0,2), labels=c('Hispanic','Not Hispanic','Unknown'))) %>%
					mutate(smoking_status = factor(smoking_status, levels=c(1,2,0), labels=c('Current','Former','Never'))) %>%
					mutate(dose2_month = factor(dose2_month, levels=c(1,2,3,4,5,6,7,8,9,10,11,12),
								labels=c('January','February','March','April','May','June','July','August','September','October','November','December'))) %>%
					mutate(Tests0 = case_when(ntests==0 ~ 1, TRUE ~ 0)) %>%
					mutate(Tests1 = case_when(ntests==1 ~ 1, TRUE ~ 0)) %>%
					mutate(Tests2 = case_when(ntests==2 ~ 1, TRUE ~ 0)) %>%
					mutate(ntests = factor(ntests, levels=c(0,1,2), 
								labels=c('0 tests','1 test','>=2 tests'))) %>%
					mutate(diff23_3grp = factor(diff23_3grp, levels=c(1,2,3), labels=c('6-<8','8-<9','>=9'))) %>%
					mutate(pcp5 = factor(pcp5, levels=c(0,1,2,3), labels=c('1-9','10-19','20-29','>=30'))) %>%
					mutate(flu5 = factor(flu5, levels=c(0,1,2,3), labels=c('0 vaccinations','1 or 2','3 or 4','>=5'))) %>%
					mutate(BNT162b2 = case_when(dose2_man=='Pfizer' ~ 1, TRUE ~ 0)) %>%
					mutate(mRNA1273 = case_when(dose2_man=='Moderna' ~ 1, TRUE ~ 0)) %>%
					mutate(nopriorcov = case_when(covidpriordose3==1 ~ 0, TRUE ~ 1))


table_matched <- CreateTableOne(vars=myVars,
								data=matched_PM_merged,
								factorVars = catVars,
								strata="dose3",
								addOverall=TRUE)

dir <- "[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/"

#save table_matched (csv)
write.csv(print(table_matched,nonnormal = c('age_at_index','tot_fup',"testcnt_post"),
			quote=FALSE,noSpaces=TRUE,printToggle = FALSE),
          file=paste0(dir,"table-matched-big_PvM-exact-delom_booster.csv"))

endsubmit;
quit;


*********************************************;
** OMICRON Booster ERA : Matched table	;
*********************************************;

proc means data=sasin.matched_PM_omboost mean median q1 q3 range sum;
	class dose3_man;
	var testcnt_post;
run;

proc iml;

run ExportDataSettoR("sasin.matched_PM_omboost","matched_PM_merged");
submit / R;

pacman::p_load(MatchIt,ggQC,cobalt,ggplot2,tidyr,viridis,extrafont,dplyr,survival,survminer,tableone)


names(matched_PM_merged) <- tolower(names(matched_PM_merged))

myVars = c("age_at_index","age_cat","male","female","race","ethnicity",
"urbanicity","smoking_status","cond_lung","cond_vasc","cond_hypt","cond_diab","cond_ckd",
"cond_cld","cond_5yr_cancer","immunsupp",
"obesity","cond_dement","cond_sud",
"diff23_3grp","BNT162b2","mRNA1273","ntests","pcp5","flu5",
"dose2_month","covidpriordose3","nopriorcov","same_primary","dose3_booster","tot_fup","testcnt_post")

catVars=c("age_cat","male","female","race","ethnicity",
"urbanicity","smoking_status","cond_lung","cond_vasc","cond_hypt","cond_diab","cond_ckd",
"cond_cld","cond_5yr_cancer","immunsupp",
"obesity","cond_dement","cond_sud",
"diff23_3grp","BNT162b2","mRNA1273","ntests","pcp5","flu5",
"dose2_month","covidpriordose3","nopriorcov","same_primary","dose3_booster")


matched_PM_merged <- matched_PM_merged %>% 
					mutate(current_smoker = case_when(smoking_status==1 ~ 1,
					TRUE ~ 0)) %>%
					mutate(dose3 = factor(group_binary, levels=c(1,0), labels=c('BNT162b2','mRNA-1273'))) %>%
					mutate(same_primary = case_when(dose2_man==dose3_man ~ 1, 
					TRUE ~ 0)) %>% # variable denoting whether the primary manufacturer matches dose 3 (01mar2022 HG)
					#mutate(sex = factor(sex, levels=c(1,0), labels=c('Male','Female'))) %>%
					mutate(age_cat = factor(age_cat, levels = c(0,1,2,3,4,5),
											labels = c('18 to 39 years','40 to 49 years','50 to 59 years',
														'60 to 69 years', '70 to 79 years','>=80 years'))) %>%
					mutate(male = case_when(sex==1 ~ 1, TRUE ~ 0)) %>%
					mutate(female = case_when(sex==0 ~ 1, TRUE ~ 0)) %>%
					mutate(black = case_when(race==1 ~ 1, TRUE ~ 0)) %>%
					mutate(white = case_when(race==0 ~ 1, TRUE ~ 0)) %>%
					mutate(race = factor(race, levels=c(1,2,3,0), labels=c('Black','Other','Unknown','White'))) %>%
					mutate(ethnicity = factor(ethnicity, levels=c(1,0,2), labels=c('Hispanic','Not Hispanic','Unknown'))) %>%
					mutate(smoking_status = factor(smoking_status, levels=c(1,2,0), labels=c('Current','Former','Never'))) %>%
					mutate(dose2_month = factor(dose2_month, levels=c(1,2,3,4,5,6,7,8,9,10,11,12),
								labels=c('January','February','March','April','May','June','July','August','September','October','November','December'))) %>%
					mutate(Tests0 = case_when(ntests==0 ~ 1, TRUE ~ 0)) %>%
					mutate(Tests1 = case_when(ntests==1 ~ 1, TRUE ~ 0)) %>%
					mutate(Tests2 = case_when(ntests==2 ~ 1, TRUE ~ 0)) %>%
					mutate(ntests = factor(ntests, levels=c(0,1,2), 
								labels=c('0 tests','1 test','>=2 tests'))) %>%
					mutate(diff23_3grp = factor(diff23_3grp, levels=c(1,2,3), labels=c('6-<8','8-<9','>=9'))) %>%
					mutate(pcp5 = factor(pcp5, levels=c(0,1,2,3), labels=c('1-9','10-19','20-29','>=30'))) %>%
					mutate(flu5 = factor(flu5, levels=c(0,1,2,3), labels=c('0 vaccinations','1 or 2','3 or 4','>=5'))) %>%
					mutate(BNT162b2 = case_when(dose2_man=='Pfizer' ~ 1, TRUE ~ 0)) %>%
					mutate(mRNA1273 = case_when(dose2_man=='Moderna' ~ 1, TRUE ~ 0)) %>%
					mutate(nopriorcov = case_when(covidpriordose3==1 ~ 0, TRUE ~ 1))


table_matched <- CreateTableOne(vars=myVars,
								data=matched_PM_merged,
								factorVars = catVars,
								strata="dose3",
								addOverall=TRUE)

dir <- "[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/"

#save table_matched (csv)
write.csv(print(table_matched,nonnormal =c('age_at_index','tot_fup',"testcnt_post"), #added 08mar2022
			quote=FALSE,noSpaces=TRUE,printToggle = FALSE),
          file=paste0(dir,"table-matched-big_PvM-exact-omicron_boost.csv"))

endsubmit;
quit;

*********************************************;
** ELIGIBILITY TABLES	;
*********************************************;

*********************************************;
** DELTA-OMICRON COMBINED TIME PERIOD	;
*********************************************;
proc means data=sasin.elig_pm_delom mean median q1 q3 range sum;
	class dose3_man;
	var testcnt_post;
run;

proc iml;

run ExportDataSettoR("sasin.elig_pm_delom","eligibles");
submit / R;

pacman::p_load(MatchIt,ggQC,cobalt,ggplot2,tidyr,viridis,extrafont,dplyr,survival,survminer,tableone)


names(eligibles) <- tolower(names(eligibles))

myVars = c("age_at_index","age_cat","male","female","race","ethnicity",
"urbanicity","smoking_status","cond_lung","cond_vasc","cond_hypt","cond_diab","cond_ckd",
"cond_cld","cond_5yr_cancer","immunsupp",
"obesity","cond_dement","cond_sud",
"diff23_3grp","BNT162b2","mRNA1273","ntests","pcp5","flu5",
"dose2_month","covidpriordose3","nopriorcov","same_primary","dose3_booster","tot_fup","testcnt_post")

catVars=c("age_cat","male","female","race","ethnicity",
"urbanicity","smoking_status","cond_lung","cond_vasc","cond_hypt","cond_diab","cond_ckd",
"cond_cld","cond_5yr_cancer","immunsupp",
"obesity","cond_dement","cond_sud",
"diff23_3grp","BNT162b2","mRNA1273","ntests","pcp5","flu5",
"dose2_month","covidpriordose3","nopriorcov","same_primary","dose3_booster")


eligibles <- eligibles %>% 
					mutate(dose3 = factor(group_binary, levels=c(1,0), labels=c('BNT162b2','mRNA-1273'))) %>%
					mutate(same_primary = case_when(dose2_man==dose3_man ~ 1, 
					TRUE ~ 0)) %>% # variable denoting whether the primary manufacturer matches dose 3 (01mar2022 HG)
					mutate(age_cat = factor(age_cat, levels = c(0,1,2,3,4,5),
											labels = c('18 to 39 years','40 to 49 years','50 to 59 years',
														'60 to 69 years', '70 to 79 years','>=80 years'))) %>%
					mutate(male = case_when(sex==1 ~ 1, TRUE ~ 0)) %>%
					mutate(female = case_when(sex==0 ~ 1, TRUE ~ 0)) %>%
					mutate(race = factor(race, levels=c(1,2,3,0), labels=c('Black','Other','Unknown','White'))) %>%
					mutate(ethnicity = factor(ethnicity, levels=c(1,0,2), labels=c('Hispanic','Not Hispanic','Unknown'))) %>%
					mutate(smoking_status = factor(smoking_status, levels=c(1,2,0), labels=c('Current','Former','Never'))) %>%
					mutate(dose2_month = factor(dose2_month, levels=c(1,2,3,4,5,6,7,8,9,10,11,12),
								labels=c('January','February','March','April','May','June','July','August','September','October','November','December'))) %>%
					mutate(ntests = factor(ntests, levels=c(0,1,2), 
								labels=c('0 tests','1 test','>=2 tests'))) %>%
					mutate(diff23_3grp = factor(diff23_3grp, levels=c(1,2,3), labels=c('6-<8','8-<9','>=9'))) %>%
					mutate(pcp5 = factor(pcp5, levels=c(0,1,2,3), labels=c('1-9','10-19','20-29','>=30'))) %>%
					mutate(flu5 = factor(flu5, levels=c(0,1,2,3), labels=c('0 vaccinations','1 or 2','3 or 4','>=5'))) %>%
					mutate(BNT162b2 = case_when(dose2_man=='Pfizer' ~ 1, TRUE ~ 0)) %>%
					mutate(mRNA1273 = case_when(dose2_man=='Moderna' ~ 1, TRUE ~ 0)) %>%
					mutate(nopriorcov = case_when(covidpriordose3==1 ~ 0, TRUE ~ 1))

table_elig_PM <- CreateTableOne(vars=myVars,
								data=eligibles,
								factorVars = catVars,
								strata="dose3",
								addOverall=TRUE)

dir <- "[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/"

#save table_elig_PM (csv)
write.csv(print(table_elig_PM,nonnormal=c('age_at_index','tot_fup',"testcnt_post"),quote=FALSE,noSpaces=TRUE,printToggle = FALSE),
          file=paste0(dir,"table-elig-big_PvM-exact-delom_booster.csv"))



endsubmit;
quit;

*********************************************;
** OMICRON Booster ERA	;
*********************************************;

proc means data=sasin.elig_pm_omboost mean median q1 q3 range sum;
	class dose3_man;
	var testcnt_post;
run;

proc iml;

run ExportDataSettoR("sasin.elig_pm_omboost","eligibles");
submit / R;

pacman::p_load(MatchIt,ggQC,cobalt,ggplot2,tidyr,viridis,extrafont,dplyr,survival,survminer,tableone)


names(eligibles) <- tolower(names(eligibles))

myVars = c("age_at_index","age_cat","male","female","race","ethnicity",
"urbanicity","smoking_status","cond_lung","cond_vasc","cond_hypt","cond_diab","cond_ckd",
"cond_cld","cond_5yr_cancer","immunsupp",
"obesity","cond_dement","cond_sud",
"diff23_3grp","BNT162b2","mRNA1273","ntests","pcp5","flu5",
"dose2_month","covidpriordose3","nopriorcov","same_primary","dose3_booster","tot_fup","testcnt_post")

catVars=c("age_cat","male","female","race","ethnicity",
"urbanicity","smoking_status","cond_lung","cond_vasc","cond_hypt","cond_diab","cond_ckd",
"cond_cld","cond_5yr_cancer","immunsupp",
"obesity","cond_dement","cond_sud",
"diff23_3grp","BNT162b2","mRNA1273","ntests","pcp5","flu5",
"dose2_month","covidpriordose3","nopriorcov","same_primary","dose3_booster")


eligibles <- eligibles %>% 
					mutate(dose3 = factor(group_binary, levels=c(1,0), labels=c('BNT162b2','mRNA-1273'))) %>%
					mutate(same_primary = case_when(dose2_man==dose3_man ~ 1, 
					TRUE ~ 0)) %>% # variable denoting whether the primary manufacturer matches dose 3 (01mar2022 HG)
					mutate(age_cat = factor(age_cat, levels = c(0,1,2,3,4,5),
											labels = c('18 to 39 years','40 to 49 years','50 to 59 years',
														'60 to 69 years', '70 to 79 years','>=80 years'))) %>%
					mutate(male = case_when(sex==1 ~ 1, TRUE ~ 0)) %>%
					mutate(female = case_when(sex==0 ~ 1, TRUE ~ 0)) %>%
					mutate(race = factor(race, levels=c(1,2,3,0), labels=c('Black','Other','Unknown','White'))) %>%
					mutate(ethnicity = factor(ethnicity, levels=c(1,0,2), labels=c('Hispanic','Not Hispanic','Unknown'))) %>%
					mutate(smoking_status = factor(smoking_status, levels=c(1,2,0), labels=c('Current','Former','Never'))) %>%
					mutate(dose2_month = factor(dose2_month, levels=c(1,2,3,4,5,6,7,8,9,10,11,12),
								labels=c('January','February','March','April','May','June','July','August','September','October','November','December'))) %>%
					mutate(ntests = factor(ntests, levels=c(0,1,2), 
								labels=c('0 tests','1 test','>=2 tests'))) %>%
					mutate(diff23_3grp = factor(diff23_3grp, levels=c(1,2,3), labels=c('<8','8','>=9'))) %>%
					mutate(pcp5 = factor(pcp5, levels=c(0,1,2,3), labels=c('1-9','10-19','20-29','>=30'))) %>%
					mutate(flu5 = factor(flu5, levels=c(0,1,2,3), labels=c('0 vaccinations','1 or 2','3 or 4','>=5'))) %>%
					mutate(BNT162b2 = case_when(dose2_man=='Pfizer' ~ 1, TRUE ~ 0)) %>%
					mutate(mRNA1273 = case_when(dose2_man=='Moderna' ~ 1, TRUE ~ 0)) %>%
					mutate(nopriorcov = case_when(covidpriordose3==1 ~ 0, TRUE ~ 1))

table_elig_PM <- CreateTableOne(vars=myVars,
								data=eligibles,
								factorVars = catVars,
								strata="dose3",
								addOverall=TRUE)

dir <- "[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/"

#save table_elig_PM (csv)
write.csv(print(table_elig_PM,nonnormal=c('age_at_index','tot_fup',"testcnt_post"),quote=FALSE,noSpaces=TRUE,printToggle = FALSE),
          file=paste0(dir,"table-elig-big_PvM-exact-omicron_boost.csv"))



endsubmit;
quit;
