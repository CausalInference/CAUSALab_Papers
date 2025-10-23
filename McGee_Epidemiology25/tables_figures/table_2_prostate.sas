****************************************************************************

Title: Recommendation-based Physical Activity and Dietary Interventions
       for Adults Diagnosed with Breast or Prostate Cancer 

****************************************************************************

Programmer: Emma McGee

Date: March 23, 2023

Purpose of Program: Create descriptive table of baseline characteristics

Statistical Analyses:
  Descriptive statistics

****************************************************************************;

filename hpstools '';
filename channing '';

libname  readfmt  '';

options mautosource sasautos=(sasautos hpstools channing) ;
options fmtsearch=(readfmt);
options ls=88 ps=54 nocenter obs=max spool;

libname gcomp '';

****************************************************************************;
********************     FORMAT DATA FOR TABLE        **********************
****************************************************************************;
 
proc format;

value stagef	1='T1'
				2='T2'
				3='T3'
				4='T4';


value gleasonf	1 = "<7"
				2 = "7"
				3 = ">7";


value psacatf
				1 = "<4"
				2 = "4-10"
				3 = ">=10";


value txf 		1="Radical prostatectomy"
				2="Radiation therapy"
				3="Hormones"
				4="Watchful waiting"
				5="Other"
				6="Missing";
run;


data table1;
set gcomp.bytimes_prostate;

* Define time from diagnosis to baseline;
time_dx_base = firstpostdt - dtdx_prostate;

* Define a new stage variable with T4 as a separate category;
if tnmclin=1 then stage=1;                 	/* Stage 1 (T1 N0 M0) */
else if tnmclin in(2,3,4) then stage=2;			/* Stage 2 (T2notbc N0 M0 [T2 or T2a] or T2b N0 M0 or T2c N0 M0) */
else if tnmclin in(5,6,7) then stage=3;		/* Stage 3+ but M0 (T3notb N0 M0 [T3 or T3a] or T3b N0 M0) */
else if tnmclin in (8) then stage=4;			/* Stage T4 N0 M0 or N1 M0 */
else stage=.;						                    /* Missing	*/
* Note that there is no one of type 1 (T4 N0 M0) in the final analytic dataset;


* Restrict to relevant variables;
keep 
baselineage agedx race fhxmi fhxca bmi smkhx smk
modact vigact act 
totfv totgrn junk totred totproc sodajuice totalc 
aofib
cal
dia hbp chl
dxyear
time_dx_base
tx stage gleason psacat
firstpostdt period irt98 irt00  irt02 irt04 irt06 irt08 irt10 irt12 irt14;

if period ne 0 then delete;

* Format clinical variables;
format stage stagef.;
format gleason gleasonf.;
format psacat psacatf.;
format tx txf.;
run;

* Check categorical variables;
proc freq data=table1;
tables race smkhx smk fhxmi fhxca  psacat stage gleason tx 
dia hbp chl /missing;
run;

* Check distributions of continuous variables to determine if mean(sd) or median(iqr) will be reported;
proc means data=table1 n mean median min max nmiss;
var baselineage agedx time_dx_base bmi modact vigact act
totfv totgrn junk totred totproc sodajuice totalc cal aofib;
run;

ods pdf file="check_data_distributions_table_2_prostate.pdf";
proc univariate data=table1 noprint;
var 
baselineage agedx time_dx_base bmi modact vigact act
totfv totgrn junk totred totproc sodajuice totalc cal aofib;
histogram 
baselineage agedx time_dx_base bmi modact vigact act
totfv totgrn junk totred totproc sodajuice totalc cal aofib;
run;
ods pdf close;

* PA, whole grains, processed food, red meat, processed meat, ssb and alc are right skewed;
* Will present median (IQR) in Table;

****************************************************************************;
************************     CREATE TABLE        ***************************
****************************************************************************;

%table1(data=table1,
        noexp=T,
	    varlist= baselineage race modact vigact act totfv totgrn junk totred totproc sodajuice totalc cal 
		bmi smkhx smk dia hbp chl fhxmi dxyear time_dx_base stage tx psacat gleason,
	    file = table_2_prostate,
	    cat= race fhxmi smkhx smk dia hbp chl,
	    poly= stage tx psacat gleason ,
		mdn= baselineage modact vigact act totfv totgrn junk totred totproc sodajuice totalc cal bmi dxyear time_dx_base,
		pctn=pctn,
		dec=1,
        ageadj=F,
	    header='Table 2 Prostate Cancer');
run;

* There is no missing data on any of these variables;