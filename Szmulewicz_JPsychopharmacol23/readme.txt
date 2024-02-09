/*
This is a README outlining the analytic programs used to estimate the per-protocol effect of lithium treatment in preventing recurrent suicide attempts in the following manuscript:
- Szmulewicz AG, Madenci A, Ferguson R, Liang MH, Lew R, Katz IR, Hernan MA. Estimating the per-protocol effect of lithium on suicidality in a randomized trial of individuals with depression or bipolar disorder. Journal of Psychopharmacology. In press.
  
It is assumed that the input dataset for these programs contains all necessary analytic variables, as described in the accompanying codebook. 

Authors: VA-CAUSAL Methods Core
Version: March 2023. 

*/ 

***********************
ANALYSIS
***********************

"analysis_function_final" is a single script that contains the functions required to
	1) estimate the denominators of the IP weights, 
	2) estimate the numerators of the IP weights, 
	3) run the outcome model to generate point estimates 
	4) perform bootstrap resample,
	5) prepare summary data for tables and figures
	and 6) generate confidence intervals, and Manuscript tables 

Some notes that may help clarifying the analysis_functions:

- The weights need to be estimated in those weeks in which the probability of adhering to the treatment strategies is not known. In the piece of code below, we fit the models only to the weeks in which there were clinical visits (i.e., a patient's non-adherence could only be detected at clinical visits) that were scheduled, in which there were no allowed reasons for discontinuing (e.g., toxicity), it was not the first clinical visit, and there was no concomitant use of a non-allowed medication. 

# weight models
    .dat <- .dat %>% mutate(w0 = case_when (clinical_visit==1 & unscheduled==0 & allowed_disc==0 & cumvisit>consecutive & non_adherence==0 & forb_medic==0 & visit.week>1 ~ w0, #added nonadherence & visit.week here instead ofas extra step
                                            visit.week==0 ~ 1, 
                                            non_adherence==1 ~ 0,
                                            TRUE ~ 1),
                            s_w0 = case_when (clinical_visit==1 & unscheduled==0 & allowed_disc==0 & cumvisit>consecutive & non_adherence==0 & forb_medic==0 & visit.week>1 ~ s_w0,
                                              visit.week==0 ~ 1, 
                                              non_adherence==1 ~ 0,
                                              TRUE ~ 1))


"analysis_script_final" is a single script that contains the functions required to run the following analyses. Please note that each of these analyses require a dataset with a different format:
	1) intention-to-treat analysis
	2) intention-to-treat analysis restricted to those participants with complete baseline covariates
	3) intention-to-treat analysis restricted to those participants with complete baseline covariates and censoring at dropout
	4) per-protocol analysis

