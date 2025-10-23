****************************************************************************

Title: Recommendation-based Physical Activity and Dietary Interventions
       for Adults Diagnosed with Breast or Prostate Cancer  

****************************************************************************

Programmer: Emma McGee

Date: July 1, 2024

Purpose of Program: Calculate proportion adherent to recommendations at baseline

Statistical Analyses:
  Descriptive statistics

****************************************************************************;

filename nhstools '/proj/nhsass/nhsas00/nhstools/sasautos/';
filename local '/usr/local/channing/sasautos/';
libname library '/proj/nhsass/nhsas00/formats/';
options mautosource sasautos=(local nhstools);   *** path to macros  ***;
options fmtsearch=(readfmt);                     *** path to formats ***;
options nocenter;
options linesize=125 pagesize=78 ;

libname gcomp '';

****************************************************************************;

data test;
set gcomp.bytimes_breast;

* Create baseline adherence variables;

*activity met-hours/week;
act_base=0;
if act_b>=7.5 then act_base=1;

*fruit and veg servings/day;
fv_base=0;
if totfv_b>=5 then fv_base=1;

*wholegrain/legume servings/day;
grain_base=0;
if totgrn_b>=3 then grain_base=1;

*junk servings/day;
junk_base=0;
if junk_b<=1 then junk_base=1;

*red meat - recorded in servings/week;
red_base=0;
if totred_b<=3 then red_base=1;

*processed meat - recorded in servings/week;
proc_base=0;
if totproc_b<1 then proc_base=1;

*sugar sweetened drinks;
sug_base=0;
if sodajuice_b=0 then sug_base=1;
  
*alc;
alc_base=0;
if totalc_b=0 then alc_base=1;

*all;
all_base=0;
if act_base=1 and fv_base=1 and grain_base=1 and junk_base=1 and red_base=1 and proc_base=1 and sug_base=1 then all_base=1;

* all diet;
diet_base=0;
if fv_base=1 and grain_base=1 and junk_base=1 and red_base=1 and proc_base=1 and sug_base=1 then diet_base=1;


run;

* Report % adherent at baseline;
ods tagsets.excelxp
  file = "results.xls"
  style=minimal;
proc means mean maxdec = 3;
where period = 0;
var all_base diet_base act_base fv_base grain_base junk_base red_base proc_base sug_base alc_base ;
run;
ods tagsets.excelxp close; 
