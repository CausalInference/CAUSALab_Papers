/* 

This script is a part of a set of analytic programs related to the following manuscript(s):
- Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
- Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors.  

Author(s): VA-CAUSAL Methods Core
Version: July 2023. 

*/

****************************************************;
*1. Read csv files of matched data into SAS
****************************************************;

libname sasin '[SAS_folder_ORD_Project]/Data/booster'; *files for boosters analysis;
libname sasout '[SAS_folder_ORD_Project]/Vax_Standard/analysis/data'; 

***************;
** READ-In the Data    **;
***************;

proc import datafile="[SAS_folder_ORD_Project]/Vax_Standard/analysis/data/matched-dat-PvM-covidpos-exact-delom_booster_1to1.csv"
	out=matched_PM_delom
	dbms=csv
	replace;
	getnames=yes;
run;

proc import datafile="[SAS_folder_ORD_Project]/Vax_Standard/analysis/data/matched-dat-PvM-covidpos-exact-omicron_boost_1to1.csv"
	out=matched_PM_omicron_boost
	dbms=csv
	replace;
	getnames=yes;
run;

proc import datafile="[SAS_folder_ORD_Project]/Vax_Standard/analysis/data/cleaned-dat-PvM-covidpos-delom_booster.csv"
	out=cleaned_PM_delom
	dbms=csv
	replace;
	getnames=yes;
run;

proc import datafile="[SAS_folder_ORD_Project]/Vax_Standard/analysis/data/cleaned-dat-PvM-covidpos-omicron_boost.csv"
	out=cleaned_PM_omicron_boost
	dbms=csv
	replace;
	getnames=yes;
run;


****************************************************;
*2. Sort 
****************************************************;
proc sort data=matched_PM_delom; by newid; run;
proc sort data=cleaned_PM_delom; by newid; run;
proc sort data=sasin.finaldata_delom_death; by newid caldate; run;

proc sort data=matched_PM_omicron_boost; by newid; run;
proc sort data=cleaned_PM_omicron_boost; by newid; run;
proc sort data=sasin.finaldata_om_death; by newid caldate; run;


****************************************************;
*3. Merge with big data
****************************************************;
*We are only concerned with baseline covariates, 
so only keep first row for each person in the big data to reduce size;

options compress=binary; 

********************;
** P v M - Combined Delta-Omicron period **;
********************;

proc sql;
	select count(distinct patienticn), dose3_man
	from sasin.FinalData_delom_death
	where eligible_PM=1
	group by dose3_man;
quit;

data baseline_char_delom;
merge matched_PM_delom (in=a keep=newid days rename=(days=days_matched))
		cleaned_PM_delom (in=b rename=(age_at_index=age_at_index_coarsened)) 
	  sasin.FinalData_delom_death (in=c where=(Caldate=dose3_dt) rename=(tot_fup=tot_fup_orig)); *start-date for boosters data is 20Oct2021;
by newid;
if c; *keep all individuals;

if a then matched_PM=1; 
if b then elig_PM=1;
group_binary = 'group.binary'n;

*age category;
age_cat=.;
if 18<=Age_at_index<40 then age_cat=0;
else if 40<=Age_at_index<50 then age_cat=1;
else if 50<=Age_at_index<60 then age_cat=2;
else if 60<=Age_at_index<70 then age_cat=3;
else if 70<=Age_at_index<80 then age_cat=4;
else if 80<=Age_at_index then age_cat=5;

if days>112 then tot_fup=112; else tot_fup=days; *cap follow-up at the 95%;

monthdiff23=intck('month',dose2_dt,dose3_dt,'C'); *calculate months differences based on exact elapsed time ;

if monthdiff23=. then diff23_3grp = .; 
else if monthdiff23<8 then diff23_3grp = 1; *6-7 months;
else if monthdiff23<9 then diff23_3grp = 2; *8 months;
else diff23_3grp=3; *9+ months;

keep PvM MvP newid patienticn matched_PM elig_PM group_binary
		Age_at_index age_cat sex race ethnicity visn urbanicity 
		smoking_status cond_lung cond_vasc cond_hypt cond_diab cond_ckd cond_livsev
		COND_5yr_CANCER obesity pcp5 flu5 ImmunSupp 
		cond_dement cond_sud cond_cld dose2_month ntests tot_tests
		covidpriordose3 diff23_3grp dose2_man dose3_man dose3_booster dose3_dt
		daydiff23 monthdiff23
		tot_fup tot_fup_orig days days_matched; 

run;

*Grab last date of follow-up (from covidpos) for everyone to get the counts of testing after dose 3;
proc sort data=sasin.finaldata_delom_covidpos out=testdat; by newid caldate; run;
proc sort data=baseline_char_delom; by newid; run;
data testdat2;
	merge testdat (in=a keep=newid testcnt_post pcrcnt_post) 
			baseline_char_delom (in=b keep=newid); 
	by newid;
	if a and b;
	if last.newid;
run;


proc sort data=testdat2; by newid; run;
data sasout.matched_PM_delom sasout.elig_PM_delom checkme;
	merge baseline_char_delom (in=a) testdat2 (in=b);
	by newid;
	if a;

	if matched_PM=1 then output sasout.matched_PM_delom;
	if elig_PM=1 then output sasout.elig_PM_delom;

	output checkme;
run;

proc freq data=checkme;
	tables elig_PM*matched_PM*dose3_man*group_binary*pvm/list missing;
run;



********************;
** P v M - omicron_boost **;
********************;
data baseline_char_omboost;
merge matched_PM_omicron_boost (in=a keep=newid days rename=(days=days_matched))
		cleaned_PM_omicron_boost (in=b rename=(age_at_index=age_at_index_coarsened)) 
	  sasin.finaldata_om_death (in=c where=(Caldate=dose3_dt) rename=(tot_fup=tot_fup_orig)); 

by newid;
if c; *keep all baseline characteristics for everyone;

if a then matched_PM=1;
if b then elig_PM=1;
group_binary = 'group.binary'n;

*age category;
age_cat=.;
if 18<=Age_at_index<40 then age_cat=0;
else if 40<=Age_at_index<50 then age_cat=1;
else if 50<=Age_at_index<60 then age_cat=2;
else if 60<=Age_at_index<70 then age_cat=3;
else if 70<=Age_at_index<80 then age_cat=4;
else if 80<=Age_at_index then age_cat=5;

if days>63 then tot_fup=63; else tot_fup=days; *cap at the 95%;

monthdiff23=intck('month',dose2_dt,dose3_dt,'C'); *calculate months differences based on exact elapsed time;
if monthdiff23=. then diff23_3grp = .; 
else if monthdiff23<8 then diff23_3grp = 1; *6-7 months;
else if monthdiff23<9 then diff23_3grp = 2; *8 months;
else diff23_3grp=3; *9+ months;

keep PvM MvP newid patienticn matched_PM elig_PM group_binary
		Age_at_index age_cat sex race ethnicity visn urbanicity 
		smoking_status cond_lung cond_vasc cond_hypt cond_diab cond_ckd cond_livsev
		COND_5yr_CANCER obesity pcp5 flu5 ImmunSupp 
		cond_dement cond_sud cond_cld dose2_month ntests tot_tests 
		covidpriordose3 diff23_3grp dose2_man dose3_man dose3_booster dose3_dt
		daydiff23 monthdiff23
		tot_fup tot_fup_orig days days_matched; 

run;

*Grab last date of follow-up (from covidpos) for everyone to get the counts of testing after dose 3;
proc sort data=sasin.finaldata_om_covidpos out=testdat_om; by newid caldate; run;
proc sort data=baseline_char_omboost; by newid; run;
data testdat2_om;
	merge testdat_om (in=a keep=newid testcnt_post pcrcnt_post) 
			baseline_char_omboost (in=b keep=newid); 
	by newid;
	if a and b;
	if last.newid;
run;

proc sort data=testdat2_om; by newid; run;
data sasout.matched_PM_omboost sasout.elig_PM_omboost checkme_om;
	merge baseline_char_omboost (in=a) testdat2_om (in=b);
	by newid;
	if a;

	if matched_PM=1 then output sasout.matched_PM_omboost;
	if elig_PM=1 then output sasout.elig_PM_omboost;

	output checkme_om;
run;

proc freq data=checkme_om;
	tables elig_PM*matched_PM*dose3_man*group_binary*PvM/list missing;
run;
