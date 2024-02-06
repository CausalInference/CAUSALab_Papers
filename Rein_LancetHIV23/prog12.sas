

*********************************************************************************************************************************** 

CAUSALab 

File...........: Estimation of cumulative incidence of cardiovascular events in ART-naive people with HIV 
				 using an ITT analysis, adjusting for baseline confounding and using a relaxed definition of eligibility:
				 instead of requiring HIV-RNA measurement in the current month, CD4 and HIV-RNA measurement only required in last 3 months

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

data insti; set insti;
elig_init_new = 0 ;
if last_visit_cd4 <= 3 and last_visit_rna <= 3 and cancer_v ne 1 and art_naive_l1 = 1 and art_naive = 0 and rna_v > 50 then elig_init_new = 1;
run;

/* Include INITIATORS MACRO for sequential trial emulation */
%include "path to macro";

/* before running macro: sort data by id and period */
proc sort data=insti;
by id month2;
run;

/**********************************************************************************************************************************/
/* ITT ANALYSIS IN ART-NAIVE INDIVIDUALS: SENSITIVITY ANALYSIS: USING RELAXED DEFINITION OF ELIGIBILITY
/**********************************************************************************************************************************/

%INITIATORS(
datain = insti, 
id = id,
period = month2,
treatment = has_insti, 	
outcome = event_new,
first_followup = 0,
last_followup = 47,
eligible = elig_init_new,
include_expansion_time_case = 1,                   
include_followup_time_case = 1, 
followup_time_as_spline = 1, 
outcomeCov = female age age1 minority_ethnic unknown_ethnic 
			 aquitaine_cohort athena_cohort bmc_cohort ccc_cohort coris_cohort icona_cohort 
			 primo_cohort sac_cohort shcs_cohort vacs_cohort 
			 hetero_mode idu_mode other_no_mode 
			 aids cd4_v_ln cd4_v_l1 rna_v_ln rna_v_l1 hcv_cat has_hbv missing_hbv 
			 bmi_elev bmi_missing chol_high chol_missing 
			 has_hypertension missing_hypertension smoking exsmoking missing_smoking 
			 has_abc diabetes_new has_ckd missing_ckd , 
outcomeClass = ,   
outcomeCov_var = female age age1 minority_ethnic unknown_ethnic 
			 	 aquitaine_cohort athena_cohort bmc_cohort ccc_cohort coris_cohort icona_cohort 
				 primo_cohort sac_cohort shcs_cohort vacs_cohort 
			 	 hetero_mode idu_mode other_no_mode 
			 	 aids cd4_v_ln cd4_v_l1 rna_v_ln rna_v_l1 hcv_cat has_hbv missing_hbv 
			 	 bmi_elev bmi_missing chol_high chol_missing 
			 	 has_hypertension missing_hypertension smoking exsmoking missing_smoking 
			 	 has_abc diabetes_new has_ckd missing_ckd ,
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
dataname = boot_naive_sens1
);

