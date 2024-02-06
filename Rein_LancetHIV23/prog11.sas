

*********************************************************************************************************************************** 

CAUSALab 

File...........: Five most frequently used INSTI-based and non-INSTI-based ART regimens at trial baseline (in descending order)

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

/**********************************************************************/
/********* ART-NAIVE INDIVIDUALS
/**********************************************************************/

/* FIVE MOST FREQUENTLY USED INSTI-BASED REGIMENS */
proc freq data=insti noprint;
table drugs_v / out=counts_naive_insti;
where elig_init=1 and has_insti=1;
run;
proc sort data=counts_naive_insti;
by descending count;
run;
proc print data=counts_naive_insti(obs=5);
run;

/* FIVE MOST FREQUENTLY USED NON-INSTI-BASED REGIMENS */
proc freq data=insti noprint;
table drugs_v / out=counts_naive_noninsti;
where elig_init=1 and has_insti=0;
run;
proc sort data=counts_naive_noninsti;
by descending count;
run;
proc print data=counts_naive_noninsti(obs=5);
run;


/**********************************************************************/
/********* ART-EXPERIENCED INDIVIDUALS
/**********************************************************************/

/* FIVE MOST FREQUENTLY USED INSTI-BASED REGIMENS */
proc freq data=insti noprint;
table drugs_v / out=counts_exp_insti;
where elig_switch=1 and has_insti=1;
run;
proc sort data=counts_exp_insti;
by descending count;
run;
proc print data=counts_exp_insti(obs=5);
run;

/* FIVE MOST FREQUENTLY USED NON-INSTI-BASED REGIMENS */
proc freq data=insti noprint;
table drugs_v / out=counts_exp_noninsti;
where elig_switch=1 and has_insti=0;
run;
proc sort data=counts_exp_noninsti;
by descending count;
run;
proc print data=counts_exp_noninsti(obs=5);
run;

