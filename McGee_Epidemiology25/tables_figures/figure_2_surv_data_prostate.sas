****************************************************************************

Title: Recommendation-based Physical Activity and Dietary Interventions
       for Adults Diagnosed with Breast or Prostate Cancer  

****************************************************************************

Programmer: Emma McGee

Date: July 3, 2024

Purpose of Program: Create data for survival curves (Figure 2)

Statistical Analyses:
  None

****************************************************************************;

libname results '';

data simsurv;
set results.simsurv_0_100 results.simsurv_101_200 
results.simsurv_201_300 results.simsurv_301_400 results.simsurv_401_500;
run;


data forsurv ;
set simsurv (keep = surv0-surv10 int _sample_);
run;


*** Natural course ***;

data nc (keep = time surv _sample_) ;
set forsurv(where = (int = 0)) ;
array surva{*} surv0-surv10 ;
do time = 0 to 10 ;
   i = time + 1 ;
   surv = surva[i];
   output ;
end;
run;

proc sort data = nc ;
by time _sample_ ;
run;


data nc0 ;
set nc (where = (_sample_ = 0));
keep time surv;
run;


proc univariate data = nc (where = (_sample_ > 0)) noprint ;
var surv ;
output out = ncbounds
       pctlpre = surv0 
       pctlname = _pct025 _pct975 
       pctlpts = 2.5 97.5 ;
by time ;
run;


data ncgraph ;
merge nc0 ncbounds ;
by time ;
run;

*** Intervention 1 - Recommendation-based physical activity and dietary intervention (follow all recommendations) ***;

data int1 ;
set forsurv(where = (int = 1));
array surva{*} surv0-surv10 ;
do time = 0 to 10 ;
   i = time + 1 ;
   surv = surva[i];
   output ;
end;
run;

data int1_0 ;
set int1 (where = (_sample_ = 0));
keep time surv;
run;

proc sort data = int1 ;
by time _sample_;
run;

proc univariate data = int1 (where = (_sample_ > 0)) noprint;
var surv  ;
output out = int1bounds
       pctlpre = surv1  
       pctlname = _pct025 _pct975 
       pctlpts = 2.5 97.5 ;
by time ;       
run;

data int1graph ;
merge int1_0 int1bounds ;
by time ;
run;

*** Intervention 10 - follow the alcohol recommendation ***;

data int10 ;
set forsurv(where = (int = 10));
array surva{*} surv0-surv10 ;
do time = 0 to 10 ;
   i = time + 1 ;
   surv = surva[i];
   output ;
end;
run;

data int10_0 ;
set int10 (where = (_sample_ = 0));
keep time surv;
run;

proc sort data = int10 ;
by time _sample_;
run;

proc univariate data = int10 (where = (_sample_ > 0)) noprint;
var surv  ;
output out = int10bounds
       pctlpre = surv10  
       pctlname = _pct025 _pct975 
       pctlpts = 2.5 97.5 ;
by time ;       
run;

data int10graph ;
merge int10_0 int10bounds ;
by time ;
run;

data allgraphs ;
merge 
ncgraph (rename = (surv = surv0)) 
int1graph (rename = (surv = surv1)) 
int10graph (rename = (surv = surv10));
by time ;
run;

* Save dataset for graphs;
data results.forgraphs_prostate;
set allgraphs;
run;


