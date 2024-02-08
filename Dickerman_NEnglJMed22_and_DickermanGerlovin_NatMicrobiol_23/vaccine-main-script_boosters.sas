/* 

This script is a part of a set of analytic programs related to the following manuscript(s):
- Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
- Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors.  

Author(s): VA-CAUSAL Methods Core
Version: July 2023. 

*/

/****************************************************************************************/
/*											*/
/*1. create datasets for analysis with all combinations					*/
/*											*/
/****************************************************************************************/

%let event1 = covidpos;
%let event2 = covidpossymp;
%let event3 = covidhosp;
%let event4 = covidICU;
%let event5 = COVIDdeath;
%let event6 = notcoviddeath;
%macro cleandat(treat, variant); 
	%do i=1 %to 6;
		%mygsub(jobname=vacc-dat, 
				mypgm=[SAS_folder_ORD_Project]/Vax_Standard/1-vaccine-prep-data.sas, 
				parms=&&event&i.^&treat.^&variant.)
	%end;
%mend cleandat;
%cleandat(treat=PvM,variant=delom_booster);  ** PRIMARY ANALYSIS;
%cleandat(treat=PvM,variant=omicron_boost);  

/****************************************************************************************/
/*											*/
/*	check number of people in each cleaned dataset - 				*/
/*			Use this to determine max follow-up time per analysis		*/
/*											*/
/****************************************************************************************/

** Combined delta-omicron period (Primary analysis);
proc import 
	datafile="[SAS_folder_ORD_Project]/Vax_Standard/analysis/data/cleaned-dat-PvM-covidpos-delom_booster.csv"
	out=matched_PM1
	dbms=csv
	replace;
	getnames=yes;
run;
proc freq  data=matched_PM1;
table 'group.binary'n;
run;
proc univariate data=matched_PM1;
	var days; 
run;
proc import datafile="[SAS_folder_ORD_Project]/Vax_Standard/analysis/data/cleaned-dat-PvM-COVIDdeath-delom_booster.csv"
	out=matched_PM2
	dbms=csv
	replace;
	getnames=yes;
run;
proc freq  data=matched_PM2;
table 'group.binary'n;
run;
proc univariate data=matched_PM2;
	var days; 
run;

** Omicron-only period;
proc import datafile="[SAS_folder_ORD_Project]/Vax_Standard/analysis/data/cleaned-dat-PvM-covidpos-omicron_boost.csv"
	out=matched_PM1
	dbms=csv
	replace;
	getnames=yes;
run;
proc freq  data=matched_PM1;
table 'group.binary'n;
run;
proc univariate data=matched_PM1;
	var days;
run;
proc import datafile="[SAS_folder_ORD_Project]/Vax_Standard/analysis/data/cleaned-dat-PvM-COVIDdeath-omicron_boost.csv"
	out=matched_PM2
	dbms=csv
	replace;
	getnames=yes;
run;
proc freq  data=matched_PM2;
table 'group.binary'n;
run;
proc univariate data=matched_PM2;
	var days;
run;


/****************************************************************************************/
/*											*/
/*2a. preliminaries 									*/
/*	-negative control KM plots to determine final set of matching factors		*/
/* 	-output matched dataset								*/
/* 	-output event counts and follow-up time stats (for reporting and to determine 
		which analyses to proceed with - must have >10 events)			*/			
/*											*/
/****************************************************************************************/

%let event1 = covidpos;
%let event2 = covidpossymp;
%let event3 = covidhosp ;
%let event4 = covidICU;
%let event5 = COVIDdeath;
%let event6 = notcoviddeath;
%macro prelim(treat, variant);
	%do i=1 %to 6;
	%mygsub(jobname=vaccine_prelim, 
			parms=&&event&i.^&treat.^&variant.,
			mypgm=[SAS_folder_ORD_Project]/Vax_Standard/2a-vaccine-prelim.sas);
%end;
%mend prelim;

%prelim(treat=PvM,variant=delom_booster);
%prelim(treat=PvM,variant=omicron_boost);

/*%prelim(treat=MvP,variant=delom_booster);  *alternative contrast approach;*/

/****************************************************************************************/
/*											*/
/*2. point estimates									*/
/*											*/
/****************************************************************************************/

%let event1 = covidpos;
%let event2 = covidpossymp;
%let event3 = covidhosp ;
%let event4 = covidICU;
%let event5 = COVIDdeath;
%let event6 = notcoviddeath;
%macro point_est(treat, variant);
	%do i=1 %to 6;
	%mygsub(jobname=vaccine_point, 
			parms=&&event&i.^&treat^&variant,
			mypgm=[SAS_folder_ORD_Project]/Vax_Standard/2-vaccine-analysis.sas);
%end;
%mend point_est;

%point_est(treat=PvM, variant=delom_booster);
%point_est(treat=PvM, variant=omicron_boost);

/*%point_est(treat=MvP,variant=delom_booster);  *alternative contrast approach;*/


/****************************************************************************************/
/*											*/
/*  2B. point estimates for subgroup analyses 						*/
/*											*/
/****************************************************************************************/
%let event1 = covidpos;
%let event2 = covidpossymp;
%let event3 = covidhosp ;
%let event4 = covidICU;
%let event5 = COVIDdeath;
%let event6 = notcoviddeath;
%macro point_est_subgroup(treat, variant, subgroup);
	%do i=1 %to 6;
	%mygsub(jobname=vaccine_point, 
			parms=&&event&i.^&treat^&variant^&subgroup,
			mypgm=[SAS_folder_ORD_Project]/Vax_Standard/2-vaccine-analysis-subgroup.sas);
%end;
%mend point_est_subgroup;

%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=old70);
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=young70);
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=black);
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=white);
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=diff23_67); *difference between dose2 and dose3 6-7 months;
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=diff23_8) ;*difference between dose2 and dose3 8 months;
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=diff23_9); *difference between dose2 and dose3 9+ months;
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=nopriorcov); *no history of covid before dose3;
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=trueboost); *dose3 recorded as true booster dose;
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=prim_pfizer); *primary series = pfizer; 
%point_est_subgroup(treat=PvM, variant=delom_booster, subgroup=prim_moderna); *primary series = moderna;


/****************************************************************************************************************************************/
/*																	*/	
/* 3. bootstrap 															*/
/*																	*/
/* PRIOR TO RUNNING BOOTSTRAPPING, CHECK TO SEE THAT THERE ARE AT LEAST FIVE EVENTS IN EACH TREATMENT ARM FOR A GIVEN			*/
/* COMPARISON (TREAT), OUTCOME (EVENT), AND EVALUATION PERIOD (VARIANT)									*/
/* THIS CAN BE DONE BY RUNNING THE FOLLOWING SCRIPT ON THE RApp SERVER: Vax_Standard > 3a-vaccine-SubgroupCounts_BeforeBootstrap.sas	*/
/*																	*/
/*	NOTE: In the section, below, all bootstraps have been commented out that are not feasible to run based on this criteria		*/	
/*																	*/	
/****************************************************************************************************************************************/

/* PRIOR TO RUNNING BOOTSTRAPPING, CHECK TO SEE THAT THERE ARE AT LEAST FIVE EVENTS IN EACH TREATMENT ARM FOR A GIVEN			*/
/* COMPARISON (TREAT), OUTCOME (EVENT), AND EVALUATION PERIOD (VARIANT)									*/
	%mygsub(jobname=vaccine_preboot, 
			mypgm=[SAS_folder_ORD_Project]/Vax_Standard/3a-vaccine-SubgroupCounts_BeforeBootstrap.sas);

%macro bs_mac(cores, event, treat, subgroup, variant);
	%do i=1 %to &cores.;
	%mygsub(jobname=vaccine_bs&i., 
			parms=&event.^&treat.^&i.^&cores.^&subgroup.^&variant.,
			mypgm=[SAS_folder_ORD_Project]/Vax_Standard/3-vaccine-bootstrap.sas);
%end;
%mend bs_mac;

*** COMBINED DELTA-OMICRON (PRIMARY) ANALYSIS;
%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=none, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=none, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=none, variant=delom_booster); 
%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=none, variant=delom_booster); 
%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=none, variant=delom_booster); 
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=none, variant=delom_booster); 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=old70, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=old70, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=old70, variant=delom_booster);
%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=old70, variant=delom_booster); 
%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=old70, variant=delom_booster);  
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=old70, variant=delom_booster); *negative control; 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=young70, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=young70, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=young70, variant=delom_booster);
%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=young70, variant=delom_booster); 
/*%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=young70, variant=delom_booster);  **Not performed, too few cases; */ 
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=young70, variant=delom_booster); *negative control; 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=white, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=white, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=white, variant=delom_booster);
%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=white, variant=delom_booster); 
%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=white, variant=delom_booster);  
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=white, variant=delom_booster); *negative control; 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=black, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=black, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=black, variant=delom_booster);
%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=black, variant=delom_booster); 
/*%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=black, variant=delom_booster);  **Not performed, too few cases; */
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=black, variant=delom_booster); *negative control; 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=diff23_67, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=diff23_67, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=diff23_67, variant=delom_booster);
%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=diff23_67, variant=delom_booster); 
/*%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=diff23_67, variant=delom_booster);  **Not performed, too few cases; */
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=diff23_67, variant=delom_booster); *negative control; 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=diff23_8, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=diff23_8, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=diff23_8, variant=delom_booster);
%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=diff23_8, variant=delom_booster); 
%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=diff23_8, variant=delom_booster);  
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=diff23_8, variant=delom_booster); *negative control; 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=diff23_9, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=diff23_9, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=diff23_9, variant=delom_booster);
/*%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=diff23_9, variant=delom_booster); **Not performed, too few cases; */
%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=diff23_9, variant=delom_booster);  
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=diff23_9, variant=delom_booster); *negative control; 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=nopriorcov, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=nopriorcov, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=nopriorcov, variant=delom_booster);
%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=nopriorcov, variant=delom_booster); 
%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=nopriorcov, variant=delom_booster);  
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=nopriorcov, variant=delom_booster); *negative control; 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=trueboost, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=trueboost, variant=delom_booster); 
%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=trueboost, variant=delom_booster);
%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=trueboost, variant=delom_booster); 
%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=trueboost, variant=delom_booster);  
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=trueboost, variant=delom_booster); *negative control; 

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=prim_pfizer, variant=delom_booster); 
/*%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=prim_pfizer, variant=delom_booster);   **Not performed, too few cases; */
/*%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=prim_pfizer, variant=delom_booster);  **Not performed, too few cases; */
/*%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=prim_pfizer, variant=delom_booster);   **Not performed, too few cases; */
/*%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=prim_pfizer, variant=delom_booster);    **Not performed, too few cases; */
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=prim_pfizer, variant=delom_booster); *negative control;  

%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=prim_moderna, variant=delom_booster); 
%bs_mac(cores=20, event=covidpossymp, treat=PvM, subgroup=prim_moderna, variant=delom_booster);  
/*%bs_mac(cores=20, event=covidhosp, treat=PvM, subgroup=prim_moderna, variant=delom_booster);  **Not performed, too few cases; */
/*%bs_mac(cores=20, event=covidICU, treat=PvM, subgroup=prim_moderna, variant=delom_booster);   **Not performed, too few cases; */
/*%bs_mac(cores=20, event=COVIDdeath, treat=PvM, subgroup=prim_moderna, variant=delom_booster);   */
%bs_mac(cores=20, event=notcoviddeath, treat=PvM, subgroup=prim_moderna, variant=delom_booster); *negative control;  

*** OMICRON-ONLY ANALYSIS;
%bs_mac(cores=20, event=covidpos, treat=PvM, subgroup=none, variant=omicron_boost); 


%gsub_status();

/****************************************************************************************/
/*											*/
/* 4. Aggregate the results	(Tables for Boosters)					*/
/*											*/
/****************************************************************************************/

** Either run this as a job (as below), or navigate to the program and execute it using the SASApp server;
%mygsub(jobname=boosters_prepformanuscript, 
			mypgm=[SAS_folder_ORD_Project]/Vax_Standard/4-Booster_PrepTablesFiguresData.sas);
			

/****************************************************************************************/
/*											*/
/* 5. Tables and Figures for manuscript							*/
/*											*/
/****************************************************************************************/
** Either run this as a job (as below), or navigate to the program and execute it using the RApp server;

* This program currently outputs the full set of bootstrapped results across all subgroups, time-periods, and events;
			*Added outputs for eligible and matched characteristics of all scenarios;
%mygsub(jobname=boosters_tables, 
			mypgm=[SAS_folder_ORD_Project]/Vax_Standard/5-Booster_ManuscriptTables.sas);

** All figures considered in manuscript and appendix;
%mygsub(jobname=boosters_figures, 
			mypgm=[SAS_folder_ORD_Project]/Vax_Standard/5-Booster_ManuscriptFigures.sas);
			

/****************************************************************************************/
/*											*/
/* NEXT STEPS:										*/
/*											*/
/****************************************************************************************/

* Copy all the BoosterManuscript folder to the P drive (Needed for downloading);
%sysexec( cp -a -p [SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/
				/[Pdrive_linux_filepath]/Vax_Standard);


/****************************************************************************************/
/*											*/
/* MOVING FILES BETWEEN SERVERS:							*/
/*											*/
/****************************************************************************************/

* Creating a backup in my folder of the files (output*);
%sysexec( cp -a -r -p [SAS_folder_ORD_Project]/Vax_Standard/analysis/bootstrap/output*
				/data/dart/[ORD_PROJECT]/Programs/[LastName]/analysis/bootstrap);

* Remove old files (output*) from P drive folder;
%sysexec (rm -f /[Pdrive_linux_filepath]/Vax_Standard/analysis/bootstrap/output*);

* Copy all the bootstrap files to P drive folder;
%sysexec( cp -a -p [SAS_folder_ORD_Project]/Vax_Standard/analysis/bootstrap/output*
				/[Pdrive_linux_filepath]/Vax_Standard/analysis/bootstrap);

* Archive all the results files to P drive folder;
%sysexec( mv /[Pdrive_linux_filepath]/boosters/FromGrid/analysis/results/output*
				/[Pdrive_linux_filepath]/Vax_Standard/analysis/results/output_archive);

* Copy all the results files to P drive folder;
%sysexec( cp -a -p [SAS_folder_ORD_Project]/Vax_Standard/analysis/results/output*
				/[Pdrive_linux_filepath]/Vax_Standard/analysis/results);

