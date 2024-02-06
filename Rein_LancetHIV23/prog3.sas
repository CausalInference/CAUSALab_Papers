

*********************************************************************************************************************************** 

CAUSALab 

File...........: Estimation of cumulative incidence of cardiovascular events in ART-naive people with HIV 
				 using an ITT analysis and adjusted for time, cohort, age and sex

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

/* Include INITIATORS MACRO for sequential trial emulation */
%include "path to macro";

/* before running macro: sort data by id and period */
proc sort data=insti;
by id month2;
run;

/**********************************************************************************************************************************/
/* ITT ANALYSIS IN ART-NAIVE INDIVIDUALS: MODEL ADJUSTED FOR TIME, COHORT, AGE AND SEX */
/**********************************************************************************************************************************/

%INITIATORS(
datain = insti, 
id = id,
period = month2,
treatment = has_insti, 	
outcome = event_new,
first_followup = 0,
last_followup = 47,
eligible = elig_init,
include_expansion_time_case = 1,                   
include_followup_time_case = 1, 
followup_time_as_spline = 1, 
outcomeCov = female age age1 
			 aquitaine_cohort athena_cohort bmc_cohort ccc_cohort coris_cohort icona_cohort 
			 primo_cohort sac_cohort shcs_cohort vacs_cohort , 
outcomeClass = ,   
outcomeCov_var = female age age1 
			 	 aquitaine_cohort athena_cohort bmc_cohort ccc_cohort coris_cohort icona_cohort 
				 primo_cohort sac_cohort shcs_cohort vacs_cohort ,
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
dataname = boot_naive_itt_2
);

