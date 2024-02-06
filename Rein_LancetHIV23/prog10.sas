

*********************************************************************************************************************************** 

CAUSALab 

File...........: Estimation of cumulative incidence of cardiovascular events in ART-experienced people with HIV 
				 using a per-protocol analysis with censoring if/when individuals deviate from their initial treatment strategy for >2 months

Project........: Integrase strand-transfer inhibitor use and cardiovascular events in adults with HIV: 
				 an emulation of target trials in the HIV-CAUSAL Collaboration 
				 and the Antiretroviral Therapy Cohort Collaboration

Author.........: Sophia Rein

Date Created...: 2nd of May 2023

*------------------------------------------------------------------------------------------------------------

Purpose........: Analysis for paper / CAUSALab transparency initiative

*------------------------------------------------------------------------------------------------------------

***********************************************************************************************************************************; 
data experienced_pp; set 'path to data'; run;

/* Include INITIATORS MACRO for sequential trial emulation */
%include "path to macro";

/* before running macro: sort data by id and period */
proc sort data=experienced_pp;
by id month2;
run;

/**********************************************************************************************************************************/
/* PER-PROTOCOL ANALYSIS IN ART-EXPERIENCED INDIVIDUALS */
/**********************************************************************************************************************************/

%INITIATORS(
datain = experienced_pp,
id = id,
period = month2,
treatment = has_insti_pp, 	
outcome = event_new,
first_followup = 0,
last_followup = 47,
eligible = elig_switch,
eligible_wts_0=  , 
eligible_wts_1=  ,
include_expansion_time_case = 1,                 
include_followup_time_case = 1, 
followup_time_as_spline = 1,  
outcomeCov = ,    
outcomeClass = , 
outcomeCov_var = , 
use_censor = 1 ,
use_weights= 1, 
use_stabilized_weights = 0 ,
	run_unweighted_analysis = 0,
	run_weighted_analysis = 0,
    run_p99_analysis = 1,
    run_user_limits_analysis = 0,                       
    lower_weight =  ,     
    upper_weight =  ,   
model_switchn = ,
cov_switchn = ,      
class_switchn = ,    
model_switchd = female age age1 tsinceart tsinceart1 
			 minority_ethnic unknown_ethnic 
			 aquitaine_cohort athena_cohort bmc_cohort ccc_cohort coris_cohort icona_cohort 
			 primo_cohort sac_cohort seropri_cohort shcs_cohort vacs_cohort 
			 hetero_mode idu_mode other_no_mode aids cd4_v_ln cd4_v_l1 hcv_cat has_hbv missing_hbv 
			 bmi_elev bmi_missing chol_high chol_missing 
			 has_hypertension missing_hypertension abc_6mo_bl smoking exsmoking missing_smoking
			 diabetes_new has_ckd missing_ckd , 
cov_switchd = female age age1 tsinceart tsinceart1 
			 minority_ethnic unknown_ethnic 
			 aquitaine_cohort athena_cohort bmc_cohort ccc_cohort coris_cohort icona_cohort 
			 primo_cohort sac_cohort seropri_cohort shcs_cohort vacs_cohort 
			 hetero_mode idu_mode other_no_mode aids cd4_v_ln cd4_v_l1 hcv_cat has_hbv missing_hbv 
			 bmi_elev bmi_missing chol_high chol_missing 
			 has_hypertension missing_hypertension abc_6mo_bl smoking exsmoking missing_smoking
			 diabetes_new has_ckd missing_ckd , 
class_switchd = , 
cense		 = , 
model_censen = ,
cov_censen   = ,      
class_censen = ,    
model_censed =  ,
cov_censed   =  ,
class_censed = , 
eligible_cense = ,
final_analysis = 1, 
run_diagnostics = 0 ,
final_method = 0 ,            
calculate_var = 0 , 
print_option = , 
risk = 1,
sorted = 1,
bootstrap = 1,  
nboot = 500,
bootstart = 0,
bootend = 500,
seed = 1234,                                               
datalib = work, 
dataname = boot_exp_pp
);


