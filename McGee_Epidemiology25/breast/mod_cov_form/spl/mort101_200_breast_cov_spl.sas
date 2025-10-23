****************************************************************************

Title: Recommendation-based Physical Activity and Dietary Interventions
       for Adults Diagnosed with Breast or Prostate Cancer  

****************************************************************************

Programmer: Emma McGee

Date: July 1, 2024

Purpose of Program: Run parametric g-formula for the effect of a recommendation-based 
                    physical activity and dietary intervention on 20-year risk of mortality 
                    Change time-varying covariates which are modeled as categorical to splines
                    Spline knots are placed near the 10th, 50th, and 90th percentiles of the distribution (except for BMI which was already set close to these percentiles at meaningful categories)

Statistical Analyses:
  Parametric g-formula

****************************************************************************;

filename nhstools '/proj/nhsass/nhsas00/nhstools/sasautos'; 
filename hpstools '/proj/hpsass/hpsas00/hpstools/sasautos/';
filename channing '/usr/local/channing/sasautos';
libname nhsfmt '/proj/nhsass/nhsas00/formatsv9/';

options fmtsearch=(nhsfmt);    
options linesize=130 pagesize=78;

%include 'gformula4.0.sas';

libname results '';
libname gcomp '';

footnote 'Effects of WCRF interventions on 20-year risk of death -- Breast cancer';

****************************************************************************;

****************************************************************************;
*****************     SPECIFY THE INTERVENTION(S)      *********************
****************************************************************************;


* INTERVENTION 1 - All recommendations;
%let interv1  =
    intno     = 1,
    intlabel  = 'All participants follow all WCRF/AICR recommendations in all intervals',
    intcond   = (xcond = 0),
    nintvar   = 7,
    intvar1   = act,
    inttype1  = 2,
    intmin1   = 7.5,
    intpr1    = 1,
    inttimes1 = 0 1 2 3 4 5 6 7 8 9,
    intvar2   = totfv,
    inttype2  = 2,
    intmin2   = 5,
    intpr2    = 1,
    inttimes2 = 0 1 2 3 4 5 6 7 8 9,
    intvar3   = totgrn,
    inttype3  = 2,
    intmin3   = 3,
    intpr3    = 1,
    inttimes3 = 0 1 2 3 4 5 6 7 8 9,
    intvar4   = junk,
    inttype4  = 2,
    intmax4   = 1,
    intpr4    = 1,
    inttimes4 = 0 1 2 3 4 5 6 7 8 9,
    intvar5   = totred,
    inttype5  = 2,
    intmax5   = 3,
    intpr5    = 1,
    inttimes5 = 0 1 2 3 4 5 6 7 8 9,
    intvar6   = totproc,
    inttype6  = 2,
    intmax6   = 0.99,
    intpr6    = 1,
    inttimes6 = 0 1 2 3 4 5 6 7 8 9,
    intvar7   = sodajuice,
    inttype7  = 1,
    intvalue7 = 0,
    intpr7    = 1,
    inttimes7 = 0 1 2 3 4 5 6 7 8 9;

* INTERVENTION 2 - Alcohol;
%let interv2  =
    intno     = 2,
    intlabel  = 'All participants follow the alcohol recommendation in all intervals',
    intcond   = (xcond = 0),
    nintvar   = 1,
    intvar1   = totalc,
    inttype1  = 1,
    intvalue1 = 0,
    intpr1    = 1,
    inttimes1 = 0 1 2 3 4 5 6 7 8 9;



****************************************************************************;
***********************       FORMAT DATA        ***************************
****************************************************************************;

data bytimes;
    set gcomp.bytimes_breast;
   
    * Defining baseline age;
	baseage = baselineage;
	baseage_sq = baseage **2;
    baseage_cub = baseage **3;

    * Make categorical age variables;
    baseage_1=0; baseage_2=0; baseage_3=0; baseage_4=0; baseage_5=0; baseage_6=0;
    baseage_7=0; baseage_8=0; baseage_9=0;
    if baseage <45 then baseage_1=1;
    if baseage>=45 and baseage<50 then baseage_2=1;
    if baseage>=50 and baseage<55 then baseage_3=1;
    if baseage>=55 and baseage<60 then baseage_4=1;
    if baseage>=60 and baseage<65 then baseage_5=1;
	if baseage>=65 and baseage<70 then baseage_6=1;
	if baseage>=70 and baseage<75 then baseage_7=1;
    if baseage>=75 and baseage<80 then baseage_8=1;
	if baseage>=80 then baseage_9=1;
    
    * Create splines for age (knots at 50, 65, and 75);
    _kd_= (75 - 50)**.666666666666 ;  
    baseage_spl1=max((baseage-50)/_kd_,0)**3+((65-50)*max((baseage-75)/_kd_,0)**3 
    -(75-50)*max((baseage-65)/_kd_,0)**3)/(75-65);

run;

* Create non time-varying baseline variable measured at the last pre-diagnostic questionnaire or baseline questionnaire ;

data bytimes2;
    set gcomp.bytimes_breast;
    lbmi_pre=lbmi_predx;
	bmi_pre=bmi_predx;
    lcal_pre=lcal_predx;

    act_pre=act_predx;
    totfv_pre=totfv_predx;
    totgrn_pre=totgrn_predx;
    junk_pre=junk_predx;
    totred_pre=totred_predx;
    totproc_pre=totproc_predx;
    sodajuice_pre=sodajuice_predx;
    totalc_pre=totalc_predx;

	hlthy_pre=hlthy_predx;
    hlthy=hlthy_b;

    if period=0;

   keep id lbmi_pre bmi_pre lcal_pre 
        act_pre totfv_pre totgrn_pre junk_pre totred_pre totproc_pre sodajuice_pre totalc_pre
        hlthy_pre hlthy;
   run;

proc sort data=bytimes2;
   by id;
   run;
   
proc sort data=bytimes;
   by id;
   run;


* Create additional categorical variables;
data bytimes;
   merge bytimes bytimes2;
   by id;

* BMI categories: <18.5, 18.5-24.9, 25-29.9, >=30 kg/m2; 
* log transformed for modeling;
    lbmi_pre_1=0; lbmi_pre_2=0; lbmi_pre_3=0; lbmi_pre_4=0;
    if lbmi_pre <log(18.5) then lbmi_pre_1=1;
	if lbmi_pre>=log(18.5) and lbmi_pre<log(25) then lbmi_pre_2=1;
	if lbmi_pre>=log(25) and lbmi_pre<log(30) then lbmi_pre_3=1;
	if lbmi_pre>=log(30) then lbmi_pre_4=1;

* Aspirin;
    asp_b_1=0; asp_b_2=0; asp_b_3=0;
    if aspirin_b =1 then asp_b_1=1;
	if aspirin_b =2 then asp_b_2=1;
	if aspirin_b =3 then asp_b_3=1;

* Total energy intake categories;
* log transformed for modeling;
    lcal_pre_1=0; lcal_pre_2=0; lcal_pre_3=0; lcal_pre_4=0;
    if lcal_pre <6.95 then lcal_pre_1=1;
	if lcal_pre>=6.95 and lcal_pre<7.40 then lcal_pre_2=1;
	if lcal_pre>=7.40 and lcal_pre<7.80 then lcal_pre_3=1;
    if lcal_pre>=7.80 then lcal_pre_4=1;

* Physical activity; 
    act_pre_1=0; act_pre_2=0; act_pre_3=0; act_pre_4=0;
    if act_pre <0.5 then act_pre_1=1;
	if act_pre>=0.5 and act_pre<7.5 then act_pre_2=1;
    if act_pre>=7.5 and act_pre<44 then act_pre_3=1;
	if act_pre>=44 then act_pre_4=1;

* Fruit and vegetables;
    totfv_pre_1=0; totfv_pre_2=0; totfv_pre_3=0; totfv_pre_4=0;
    if totfv_pre <2 then totfv_pre_1=1;
	if totfv_pre>=2 and totfv_pre<5 then totfv_pre_2=1;
	if totfv_pre>=5 and totfv_pre<7.9 then totfv_pre_3=1;
	if totfv_pre>=7.9 then totfv_pre_4=1;

* Wholegrains;
    totgrn_pre_1=0; totgrn_pre_2=0; totgrn_pre_3=0; totgrn_pre_4=0;
    if totgrn_pre <0.4 then totgrn_pre_1=1;
	if totgrn_pre>=0.4 and totgrn_pre<1.4 then totgrn_pre_2=1;
	if totgrn_pre>=1.4 and totgrn_pre<3 then totgrn_pre_3=1;
	if totgrn_pre>=3 then totgrn_pre_4=1;

* Processed foods;
    junk_pre_1=0; junk_pre_2=0; junk_pre_3=0; junk_pre_4=0;
    if junk_pre <1 then junk_pre_1=1;
	if junk_pre>=1 and junk_pre<2 then junk_pre_2=1;
	if junk_pre>=2 and junk_pre<4.6 then junk_pre_3=1;
	if junk_pre>=4.6 then junk_pre_4=1;

* Red meat;
    totred_pre_1=0; totred_pre_2=0; totred_pre_3=0; totred_pre_4=0;
    if totred_pre <0.8 then totred_pre_1=1;
	if totred_pre>=0.8 and totred_pre<3 then totred_pre_2=1;
	if totred_pre>=3 and totred_pre<7.5 then totred_pre_3=1;
	if totred_pre>=7.5 then totred_pre_4=1;

* Processed meat;
    totproc_pre_1=0; totproc_pre_2=0; totproc_pre_3=0; totproc_pre_4=0;
    if totproc_pre <0.1 then totproc_pre_1=1;
	if totproc_pre>=0.1 and totproc_pre<1 then totproc_pre_2=1;
	if totproc_pre>=1 and totproc_pre<4 then totproc_pre_3=1;
	if totproc_pre>=4 then totproc_pre_4=1;

* Sugar-sweetened beverages;
    sodajuice_pre_1=0; sodajuice_pre_2=0; sodajuice_pre_3=0; sodajuice_pre_4=0;
    if sodajuice_pre <0.1 then sodajuice_pre_1=1;
	if sodajuice_pre>=0.1 and sodajuice_pre<1 then sodajuice_pre_2=1;
	if sodajuice_pre>=1 and sodajuice_pre<2 then sodajuice_pre_3=1;
	if sodajuice_pre>=2 then sodajuice_pre_4=1;

* Alcohol;
    totalc_pre_1=0; totalc_pre_2=0; totalc_pre_3=0; totalc_pre_4=0;
    if totalc_pre <0.01 then totalc_pre_1=1;
	if totalc_pre>=0.01 and totalc_pre<1 then totalc_pre_2=1;
	if totalc_pre>=1 and totalc_pre<3 then totalc_pre_3=1;
	if totalc_pre>=3 then totalc_pre_4=1;

* Healthy behavior score prior to diagnosis;
   hlthy_pre_1=0; hlthy_pre_2=0; hlthy_pre_3=0; hlthy_pre_4=0;
   if hlthy_pre =0 then hlthy_pre_1=1;     *note: hlthy_pre=0 will be the reference group;
   if hlthy_pre =1 then hlthy_pre_2=1; 
   if hlthy_pre =2 then hlthy_pre_3=1; 
   if hlthy_pre =3 then hlthy_pre_4=1; 

* Healthy behavior score at baseline;
   hlthy_1=0; hlthy_2=0; hlthy_3=0; hlthy_4=0;
   if hlthy =0 then hlthy_1=1;     *note: hlthy=0 will be the reference group;
   if hlthy =1 then hlthy_2=1; 
   if hlthy =2 then hlthy_3=1; 
   if hlthy =3 then hlthy_4=1; 

* Clinical covariates;

* Stage;
   stage_1=0; stage_2=0; stage_3=0;
   if stage=1 then stage_1=1;
   if stage=2 then stage_2=1;
   if stage=3 then stage_3=1; 

* Time between diagnosis and baseline;
   time_dx = firstpostdt - dtdx_breast;
   time_dx_sq = time_dx*time_dx;

* Calendar year of diagnosis;
   yeardx_1=0; yeardx_2=0; yeardx_3=0; yeardx_4=0;
   if yeardx<1995 then yeardx_1 = 1;
   if yeardx>=1995 and yeardx<2000 then yeardx_2 = 1;
   if yeardx>=2000 and yeardx<2005 then yeardx_3 = 1;
   if yeardx>=2005 then yeardx_4 = 1;

   yeardx_sq = yeardx*yeardx;
   yeardx_cub = yeardx*yeardx*yeardx;

* Treatment;
  treat_1 = 0; treat_2 = 0; treat_4 = 0; treat_5 = 0;
  if hormtx = 1 and radiotx = 0 and chemotx = 0 then treat_1 =1;
  if hormtx = 0 and radiotx = 1 and chemotx = 0 then treat_2 =1;
  if hormtx = 0 and radiotx = 0 and chemotx = 1 then treat_2 =1;
  if hormtx = 0 and radiotx = 0 and chemotx = 0 then treat_4 =1;
  if hormtx = 1 and radiotx = 1 and chemotx = 0 then treat_5 =1;
  if hormtx = 1 and radiotx = 0 and chemotx = 1 then treat_5 =1;
  if hormtx = 0 and radiotx = 1 and chemotx = 1 then treat_5 =1;
  if hormtx = 1 and radiotx = 1 and chemotx = 1 then treat_5 =1;

  radio_chemo_1 = 0; radio_chemo_2 = 0; radio_chemo_3 = 0;
  if radiotx = 1 and chemotx = 0 then radio_chemo_1 =1;
  if radiotx = 0 and chemotx = 1 then radio_chemo_1 =1;
  if radiotx = 1 and chemotx = 1 then radio_chemo_2 =1;
  if radiotx = 0 and chemotx = 0 then radio_chemo_3 =1;

* Postmenopausal hormone therapy use (current or past user vs. never user);
   pmh_predx_1=0;
   if (pmh_predx=2 or pmh_predx=3) then pmh_predx_1=1;

   cohort_1 = .;
   if cohort = 2 then cohort_1 =0;
   if cohort = 1 then cohort_1 =1;

run;

proc means data= bytimes min p10 p50 p90 max;
var act totfv totgrn junk totred totproc sodajuice totalc lcal cal lbmi bmi;
run;

****************************************************************************;
*******************     RUN PARAMETRIC G-FORMULA      ***********************
****************************************************************************;

%gformula(
    data = bytimes, 
    id = id, 
    time = period, 
    timepoints = 10,
    timeptype = concat,
    timeknots = 1 2 3 4 5 6 7 8 9,
    outc = event,
    outctype = binsurv,  

    censor = censor,
    compevent = ,
    compevent_cens = 0, 

    wherevars = qx_year,

fixedcov = baseage baseage_spl1
        stage_2 stage_3
        hormtx
        radiotx
        chemotx
        est_status   
        fhxmi
        mnp_predx 
        pmh_predx_1
        act_pre_1 act_pre_2 act_pre_3
        totfv_pre_1 totfv_pre_2 totfv_pre_3
        totgrn_pre_1 totgrn_pre_2 totgrn_pre_3
        junk_pre_1 junk_pre_2 junk_pre_3
        totred_pre_1 totred_pre_2 totred_pre_3
        totproc_pre_1 totproc_pre_2 totproc_pre_3
        sodajuice_pre_1 sodajuice_pre_2 sodajuice_pre_3
        totalc_pre_1 totalc_pre_2 totalc_pre_3
        lcal_pre_1 lcal_pre_2 lcal_pre_3
        lbmi_pre_1 lbmi_pre_2 lbmi_pre_3
        smkhx
        hbp_chl_dia_b
		,

    ncov=12,
	
    cov1 = qx_year, cov1otype = 0, cov1ptype = conbin, cov1inc = 2,  
                    cov1mtype = nocheck,
    cov2  = xcond,  cov2otype   = 2,  cov2ptype  = tsswitch1,
    cov3  = lbmi,   cov3otype   = 3,  cov3ptype  = lag1spl,  cov3knots = 2.9177707 3.2188758 3.4011974,
    cov4  = lcal,   cov4otype   = 3,  cov4ptype  = lag1spl,   cov4knots = 6.95  7.40 7.80,
                     cov4wherem=(not (qx_year in (1988, 1992, 1996, 2000, 2004, 2008, 2012, 2014,
                                                  1989, 1993, 1997, 2001, 2005, 2009, 2013, 2017))), 
    cov5  = totalc,  cov5otype  = 4,  cov5ptype  = lag1spl,   cov5knots = 0.01  1 3, 
                     cov5wherem=(not (qx_year in (1988, 1992, 1996, 2000, 2004, 2008, 2012, 2014,
                                                  1989, 1993, 1997, 2001, 2005, 2009, 2013, 2017))),  
    cov6  = sodajuice, cov6otype  = 4,  cov6ptype  = lag1qdc, 
                     cov6wherem=(not (qx_year in (1988, 1992, 1996, 2000, 2004, 2008, 2012, 2014,
                                                  1989, 1993, 1997, 2001, 2005, 2009, 2013, 2017))), 
    cov7  = junk,   cov7otype   = 3,  cov7ptype  = lag1qdc,
                     cov7wherem=(not (qx_year in (1988, 1992, 1996, 2000, 2004, 2008, 2012, 2014,
                                                  1989, 1993, 1997, 2001, 2005, 2009, 2013, 2017))), 
    cov8  = totproc, cov8otype  = 4,  cov8ptype  = lag1spl,   cov8knots = 0.1  1   4, 
                     cov8wherem=(not (qx_year in (1988, 1992, 1996, 2000, 2004, 2008, 2012, 2014,
                                                  1989, 1993, 1997, 2001, 2005, 2009, 2013, 2017))), 
    cov9  = totfv,  cov9otype   = 3,  cov9ptype  = lag1spl,   cov9knots = 2  5   7.9, 
                     cov9wherem=(not (qx_year in (1988, 1992, 1996, 2000, 2004, 2008, 2012, 2014,
                                                  1989, 1993, 1997, 2001, 2005, 2009, 2013, 2017))),
    cov10  = totred, cov10otype   = 3,   cov10ptype  = lag1spl,   cov10knots = 0.8  3   7.5, 
                     cov10wherem=(not (qx_year in (1988, 1992, 1996, 2000, 2004, 2008, 2012, 2014,
                                                  1989, 1993, 1997, 2001, 2005, 2009, 2013, 2017))),   
    cov11  = totgrn, cov11otype   = 3,  cov11ptype  = lag1spl,   cov11knots = 0.4  1.4   3, 
                     cov11wherem=(not (qx_year in (1988, 1992, 1996, 2000, 2004, 2008, 2012, 2014,
                                                   1989, 1993, 1997, 2001, 2005, 2009, 2013, 2017))), 
    cov12  = act,    cov12otype   = 3,  cov12ptype  = lag1spl,  cov12knots = 0.5 7.5 44,
                     cov12wherem=(not (qx_year in (1990, 2002,
                                                   1993, 1995, 1999, 2003, 2007, 2011, 2015))), 
    
    seed= 2394,
  
    check_cov_models= 1,
    print_cov_means = 0,
    save_raw_covmean = 1,

    savelib = results,
    survdata = simsurv,
    covmeandata = covmean,
    intervname = intervWCRF,
    observed_surv = obssurv,
    betadata = betadata0a,

    nsimul=10000,
    nsamples = 500, 
    sample_start = 101,
    sample_end = 200,

    rungraphs=0, 
    printlogstats = 0,

    numint=2
    );
