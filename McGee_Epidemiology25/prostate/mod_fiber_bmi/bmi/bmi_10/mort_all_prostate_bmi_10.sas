****************************************************************************

Title: Recommendation-based Physical Activity and Dietary Interventions
       for Adults Diagnosed with Breast or Prostate Cancer  

****************************************************************************

Programmer: Emma McGee 

Date: July 1, 2024

Purpose of Program: Run parametric g-formula for the the effect of a recommendation-based 
                    physical activity and dietary intervention on continuous BMI 10 years after baseline

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
              outc = bmi,
              outctype = conteofu,
              bootname = intervWCRF ,
              check_cov_models = 1 ,
              covmeandata = covmean, 
              observed_surv =  obssurv, 
              combine_survdata = 1,               
	            survdata=simsurv, 
              print_cov_means = 1,
              savecovmean = 0,
              time = period ,
              timepoints = 5,
              ncov = 10,
              numparts = 5,
              samplestart = 0 101 201 301 401,
              sampleend = 100 200 300 400 500,
              numboot = 500,
              numint = 2,
              refint = 0,
              resultsdata = results.results,
              rungraphs = 0
              );

* Create formatted results;
data results;
    set results.results;
    format sbmi sbmi_llim95 sbmi_ulim95 rd RD_llim95 RD_ulim95 comma9.1;
    sbmi_95 = cat(put(sbmi, comma9.1), ' (', strip(put(sbmi_llim95, comma9.1)), 
                  ', ', strip(put(sbmi_ulim95, comma9.1)), ')');
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
run;

* Export results to Excel;
ods tagsets.excelxp
  file = "results.xls"
  style=minimal
  options ( absolute_column_width = '15,5,10,10,10,10,10') ;
  proc print data = results noobs;
       var int2 ssize sbmi_95 rd_95 ;
  run;
ods tagsets.excelxp close; 