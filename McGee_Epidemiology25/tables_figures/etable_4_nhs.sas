****************************************************************************

Title: Recommendation-based Physical Activity and Dietary Interventions
       for Adults Diagnosed with Breast or Prostate Cancer 

****************************************************************************

Programmer: Emma McGee

Date: March 22, 2023

Purpose of Program: Create descriptive table of baseline characteristics - restricted to NHS

Statistical Analyses:
  Descriptive statistics

****************************************************************************;

filename nhstools '';
filename local '';
libname library '';
options mautosource sasautos=(local nhstools);   *** path to macros  ***;
options fmtsearch=(readfmt);                     *** path to formats ***;
options nocenter;
options linesize=125 pagesize=78 ;

libname gcomp '';

%include 'readcheck.sas';
%include 'varcheck.sas';
%include 'quintile.sas';

****************************************************************************;


****************************************************************************;
********************     FORMAT DATA FOR TABLE        **********************
****************************************************************************;
 
data bytimes;
	set gcomp.bytimes_nhs;
run;

proc format;

value stagef	1='T1'
				2='T2'
				3='T3';

run;


data table1;
set bytimes;

* Define time from diagnosis to baseline;
time_dx_base = firstpostdt - dtdx_breast;

if time_dx_base > 60 then delete;

* Postmenopausal hormone therapy use (current or past user vs. never user);
pmh_predx_1=0;
if (pmh_predx=2 or pmh_predx=3) then pmh_predx_1=1;

* Restrict to relevant variables;
keep 
baselineage white mnp_predx pmh_predx_1 
agedx fhxmi bmi smkhx smk
modact vigact act 
totfv totgrn junk totred totproc sodajuice totalc 
cal
dia hbp chl
yeardx
time_dx_base
stage
surgery chemotx radiotx hormtx est_status
firstpostdt period irt98 irt00  irt02 irt04 irt06 irt08 irt10 irt12
id
;

if period ne 0 then delete;

* Format clinical variables;
format stage stagef.;

run;

****************************************************************************;
************************     CREATE TABLE        ***************************
****************************************************************************;

%table1(data=table1,
        noexp=T,
	    varlist= baselineage white modact vigact act totfv totgrn junk totred totproc sodajuice totalc cal
		bmi smkhx smk dia hbp chl fhxmi yeardx time_dx_base stage est_status hormtx chemotx radiotx mnp_predx pmh_predx_1,
	    file = etable_4_nhs,
	    cat= white smkhx smk dia hbp chl fhxmi est_status hormtx chemotx radiotx mnp_predx pmh_predx_1,
	    poly= stage ,
		mdn= baselineage modact vigact act totfv totgrn junk totred totproc sodajuice totalc cal bmi yeardx time_dx_base,
		pctn=pctn,
		dec=1,
        ageadj=F,
	    header='eTable 4 Breast - NHS');
run;
