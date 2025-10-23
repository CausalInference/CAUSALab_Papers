****************************************************************************

Title: Recommendation-based Physical Activity and Dietary Interventions
       for Adults Diagnosed with Breast or Prostate Cancer  

****************************************************************************

Programmer: Emma McGee

Date: July 1, 2024

Purpose of Program: Run parametric g-formula for the the effect of a recommendation-based 
                    physical activity and dietary intervention on 20-year risk of mortality 
                    Change time-varying covariates which are modeled as categorical to splines
                    Spline knots are placed near the 10th, 50th, and 90th percentiles of the distribution (except for BMI which was already set close to these percentiles at meaningful categories)

Statistical Analyses:
  Parametric g-formula

****************************************************************************;

%include 'gformula4.0.sas';

libname results '';
libname gcomp '';

options linesize=88 pagesize=54;
options mprint;
options notes;

footnote 'Effects of interventions on 20-year risk of death -- Prostate cancer';

      %bootstrap_results(
              bootlib = results ,
              outc = event,
              outctype = binsurv,
              bootname = intervWCRF ,
              check_cov_models = 1 ,
              covmeandata = covmean, 
              observed_surv =  obssurv, 
              combine_survdata = 1,               
	            survdata=simsurv, 
              print_cov_means = 1,
              savecovmean = 0,
              time = period ,
              timepoints = 10,
              ncov = 11,
              numparts = 5,
              samplestart = 0 101 201 301 401,
              sampleend = 100 200 300 400 500,
              numboot = 500,
              numint = 2,
              refint = 0,
              resultsdata = results.results,
              rungraphs = 0
              );

* Check number of events;
proc freq data =  gcomp.bytimes_prostate;
    tables event / out = n_events;
    where event = 1;
run;

* Create formatted results;
data results;
    merge results.results n_events(keep = COUNT);
    format pd pd_llim95 pd_ulim95 rd RD_llim95 RD_ulim95 comma9.1;
    risk_95 = cat(put(pd, comma9.1), ' (', strip(put(pd_llim95, comma9.1)), 
                  ', ', strip(put(pd_ulim95, comma9.1)), ')');
    rd_95 = cat(put(rd, comma9.1), ' (', strip(put(RD_llim95, comma9.1)),
                   ', ', strip(put(RD_ulim95, comma9.1)), ')');
    rr_95 = cat(put(rr, comma9.2), ' (', strip(put(RR_llim95, comma9.2)),
                   ', ', strip(put(RR_ulim95, comma9.2)), ')');
    a_intervened = cat(put(averinterv, comma9.1), ' %');
    c_intervened = cat(put(intervened, comma9.1), ' %');
    if int = 0 then do;
        rd_95 = "0 (reference)";
        rr_95 = "1 (reference)";
        a_intervened = "--";
        c_intervened = "--";
    end;

    retain _count;
    if not missing(count) then _count=count;
    else count =_count;
    drop _count;

run;

* Export results to Excel;
ods tagsets.excelxp
  file = "results.xls"
  style=minimal
  options ( absolute_column_width = '15,5,10,10,10,10,10') ;
  proc print data = results noobs;
       var int2 ssize count risk_95 rd_95 rr_95 ;
  run;
ods tagsets.excelxp close; 