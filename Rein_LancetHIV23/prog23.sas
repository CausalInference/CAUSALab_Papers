

*********************************************************************************************************************************** 

CAUSALab 

File...........: Estimation of cumulative incidence of cardiovascular events in ART-experienced people with HIV 
				 using an ITT analysis and adjusting for baseline confounding and 
				 restricting INSTI initiators to top three most used INSTI regimens

Project........: Integrase strand-transfer inhibitor use and cardiovascular events in adults with HIV: 
				 an emulation of target trials in the HIV-CAUSAL Collaboration 
				 and the Antiretroviral Therapy Cohort Collaboration

Author.........: Sophia Rein

Date Created...: 2nd of May 2023

*------------------------------------------------------------------------------------------------------------

Purpose........: Analysis for paper / CAUSALab transparency initiative

*------------------------------------------------------------------------------------------------------------

***********************************************************************************************************************************; 
data insti; set 'path to data'; run;

/* Variable indicating whether INSTI initiators have taken one of the three most used regimens: 'common_drug_exp' */
data insti; set insti;
if (drugs_v="/ / 3TC ABC / DTG / /" or drugs_v="/ / FTC TDF / EVG / /" or drugs_v="/ / FTC TAF / EVG / /") and elig_switch=1 then common_drug_exp=1;
else common_drug_exp=0;
run;
data insti; set insti;
if elig_switch=1 and has_insti=1 and common_drug_exp ne 1 then elig_switch=0;
run;

/* Include INITIATORS MACRO for sequential trial emulation */
%include "path to macro";

/* before running macro: sort data by id and period */
proc sort data=insti;
by id month2;
run;

/**********************************************************************************************************************************/
/* ITT ANALYSIS IN ART-EXPERIENCED INDIVIDUALS: SENSITIVITY ANALYSIS: Restricting initiators to top three most used INSTI regimens 
/**********************************************************************************************************************************/

%INITIATORS(
datain = insti, 
id = id,
period = month2,
treatment = has_insti, 	
outcome = event_new,
first_followup = 0,
last_followup = 47,
eligible = elig_switch,
include_expansion_time_case = 1,                   
include_followup_time_case = 1, 
followup_time_as_spline = 1, 
outcomeCov = female age age1 tsinceart tsinceart1 
			 minority_ethnic unknown_ethnic aquitaine_cohort athena_cohort bmc_cohort ccc_cohort 
			 coris_cohort icona_cohort primo_cohort sac_cohort seropri_cohort shcs_cohort vacs_cohort 
			 hetero_mode idu_mode other_no_mode 
			 aids cd4_v_ln cd4_v_l1 hcv_cat has_hbv missing_hbv 
			 bmi_elev bmi_missing chol_high chol_missing 
			 has_hypertension missing_hypertension abc_6mo
			 smoking exsmoking missing_smoking 
			 diabetes_new has_ckd missing_ckd,
outcomeClass = ,   
outcomeCov_var = female age age1 tsinceart tsinceart1 
			     minority_ethnic unknown_ethnic aquitaine_cohort athena_cohort bmc_cohort ccc_cohort 
			 	 coris_cohort icona_cohort primo_cohort sac_cohort seropri_cohort shcs_cohort vacs_cohort 
			 	 hetero_mode idu_mode other_no_mode 
			 	 aids cd4_v_ln cd4_v_l1 hcv_cat has_hbv missing_hbv 
			 	 bmi_elev bmi_missing chol_high chol_missing 
			 	 has_hypertension missing_hypertension abc_6mo
			 	 smoking exsmoking missing_smoking 
			 	 diabetes_new has_ckd missing_ckd,
final_analysis = 1,  
calculate_var = 0,
final_method=0, 
use_censor = -1, 
use_weights = 0,
print_option = , 
risk = 1,
sorted = 1, 
bootstrap = 1,  
nboot = 500, 
bootstart = 0 ,
bootend = 500,
seed = 1234,                                 
datalib = work, 
dataname = boot_exp_sens2
);


