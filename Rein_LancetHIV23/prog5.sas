

*********************************************************************************************************************************** 

CAUSALab 

File...........: Estimation of cumulative incidence of cardiovascular events in ART-experienced people with HIV 
				 using an ITT analysis and only adjusting for time 

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
/* ITT ANALYSIS IN ART-EXPERIENCED INDIVIDUALS: MODEL ONLY ADJUSTED FOR TIME */
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
outcomeCov = ,
outcomeClass = ,   
outcomeCov_var = ,
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
dataname = boot_exp_itt_3
);


