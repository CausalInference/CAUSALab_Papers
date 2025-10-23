****************************************************************************

Title: Recommendation-based Physical Activity and Dietary Interventions
       for Adults Diagnosed with Breast or Prostate Cancer  

****************************************************************************

Programmer: Emma McGee

Date: July 1, 2024

Purpose of Program: Run parametric g-formula for the the effect of a recommendation-based 
                    physical activity and dietary intervention on 20-year risk of mortality 
                    Do not excuse individuals from following any of the recommendations if they develop an xcond condition

Statistical Analyses:
  Parametric g-formula

****************************************************************************;

%include 'gformula4.0.sas';

libname results '';
libname gcomp '';

footnote 'Effects of WCRF interventions on 20-year risk of death -- Prostate cancer';

****************************************************************************;
*****************     SPECIFY THE INTERVENTION(S)      *********************
****************************************************************************;

* Specify macros for skipped timepoints;
%macro lcalmacro; 
       lcal=lcal_l1;
%mend ;
%macro totfvmacro; 
       totfv=totfv_l1;
%mend ;
%macro totgrnmacro; 
       totgrn=totgrn_l1;
%mend ;
%macro junkmacro; 
       junk=junk_l1;
%mend ;
%macro totredmacro; 
       totred=totred_l1;
%mend ;
%macro totprocmacro; 
       totproc=totproc_l1;
%mend ;
%macro sodajuicemacro; 
       sodajuice=sodajuice_l1;
%mend ;
%macro totalcmacro; 
       totalc=totalc_l1;
%mend ;


* INTERVENTION 1 - All recommendations;
%let interv1  =
    intno     = 1,
    intlabel  = 'All participants follow all WCRF/AICR recommendations in all intervals',
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
    set gcomp.bytimes_prostate;
   
    * Define baseline age;
	baseage = baselineage;
	baseage_sq = baseage **2;
    baseage_cub = baseage **3;

    * Make categorical baseline age variables;
    baseage_1=0; baseage_2=0; baseage_3=0; baseage_4=0; baseage_5=0; baseage_6=0;
    if baseage <60 then baseage_1=1;
    if baseage>=60 and baseage<65 then baseage_2=1;
	if baseage>=65 and baseage<70 then baseage_3=1;
	if baseage>=70 and baseage<75 then baseage_4=1;
    if baseage>=75 and baseage<80 then baseage_5=1;
	if baseage>=80 then baseage_6=1;

    * Create splines for age (knots at 61, 71, and 80);
    _kd_= (80 - 61)**.666666666666 ;  
    baseage_spl1=max((baseage-61)/_kd_,0)**3+((71-61)*max((baseage-80)/_kd_,0)**3 
    -(80-61)*max((baseage-71)/_kd_,0)**3)/(80-71);

run;

* Create non time-varying baseline variable measured at the last pre-diagnostic questionnaire or baseline questionnaire ;

data bytimes2;
    set gcomp.bytimes_prostate;

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

data bytimes;
   merge bytimes bytimes2;
   by id;

* BMI categories: <20, 20-24.9, 25-29.9, >=30 kg/m2; 
* log transformed for modeling;
    lbmi_pre_1=0; lbmi_pre_2=0; lbmi_pre_3=0; lbmi_pre_4=0;
    if lbmi_pre <log(20) then lbmi_pre_1=1;
	if lbmi_pre>=log(20) and lbmi_pre<log(25) then lbmi_pre_2=1;
	if lbmi_pre>=log(25) and lbmi_pre<log(30) then lbmi_pre_3=1;
	if lbmi_pre>=log(30) then lbmi_pre_4=1;

* Aspirin;
    asp_b_1=0; asp_b_2=0; asp_b_3=0;
    if asp_b =1 then asp_b_1=1;
	if asp_b =2 then asp_b_2=1;
	if asp_b =3 then asp_b_3=1;

* Total energy intake categories: <1800, 1800-2499, >=2500;
* log transformed for modeling;
    lcal_pre_1=0; lcal_pre_2=0; lcal_pre_3=0;
    if lcal_pre <log(1800) then lcal_pre_1=1;
	if lcal_pre>=log(1800) and lcal_pre<log(2500) then lcal_pre_2=1;
	if lcal_pre>=log(2500) then lcal_pre_3=1;

* Physical activity; 
    act_pre_1=0; act_pre_2=0; act_pre_3=0; act_pre_4=0;
    if act_pre <1.5 then act_pre_1=1;
	if act_pre>=1.5 and act_pre<7.5 then act_pre_2=1;
    if act_pre>=7.5 and act_pre<65 then act_pre_3=1;
	if act_pre>=65 then act_pre_4=1;

* Fruits and vegetables;
    totfv_pre_1=0; totfv_pre_2=0; totfv_pre_3=0; totfv_pre_4=0;
    if totfv_pre <2.4 then totfv_pre_1=1;
	if totfv_pre>=2.4 and totfv_pre<5 then totfv_pre_2=1;
	if totfv_pre>=5 and totfv_pre<9 then totfv_pre_3=1;
	if totfv_pre>=9 then totfv_pre_4=1;
  
* Wholegrains;
    totgrn_pre_1=0; totgrn_pre_2=0; totgrn_pre_3=0; totgrn_pre_4=0;
    if totgrn_pre <0.5 then totgrn_pre_1=1;
	if totgrn_pre>=0.5 and totgrn_pre<1.5 then totgrn_pre_2=1;
	if totgrn_pre>=1.5 and totgrn_pre<3 then totgrn_pre_3=1;
	if totgrn_pre>=3 then totgrn_pre_4=1;

* Processed foods;
    junk_pre_1=0; junk_pre_2=0; junk_pre_3=0; junk_pre_4=0;
    if junk_pre <1 then junk_pre_1=1;
	if junk_pre>=1 and junk_pre<2.5 then junk_pre_2=1;
	if junk_pre>=2.5 and junk_pre<5 then junk_pre_3=1;
	if junk_pre>=5 then junk_pre_4=1;

* Red meat;
    totred_pre_1=0; totred_pre_2=0; totred_pre_3=0; totred_pre_4=0;
    if totred_pre <1 then totred_pre_1=1;
	if totred_pre>=1 and totred_pre<3 then totred_pre_2=1;
	if totred_pre>=3 and totred_pre<7 then totred_pre_3=1;
	if totred_pre>=7 then totred_pre_4=1;

* Processed meat;
    totproc_pre_1=0; totproc_pre_2=0; totproc_pre_3=0; totproc_pre_4=0;
    if totproc_pre <0.1 then totproc_pre_1=1;
	if totproc_pre>=0.1 and totproc_pre<1 then totproc_pre_2=1;
	if totproc_pre>=1 and totproc_pre<6 then totproc_pre_3=1;
	if totproc_pre>=6 then totproc_pre_4=1;

* Sugar-sweetened beverages;
    sodajuice_pre_1=0; sodajuice_pre_2=0; sodajuice_pre_3=0; sodajuice_pre_4=0;
    if sodajuice_pre <0.1 then sodajuice_pre_1=1;
	if sodajuice_pre>=0.1 and sodajuice_pre<1 then sodajuice_pre_2=1;
	if sodajuice_pre>=1 and sodajuice_pre<2 then sodajuice_pre_3=1;
	if sodajuice_pre>=2 then sodajuice_pre_4=1;

* Alcohol;
    totalc_pre_1=0; totalc_pre_2=0; totalc_pre_3=0; totalc_pre_4=0;
    if totalc_pre <0.1 then totalc_pre_1=1;
	if totalc_pre>=0.1 and totalc_pre<0.5 then totalc_pre_2=1;
	if totalc_pre>=0.5 and totalc_pre<3 then totalc_pre_3=1;
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

* Primary cancer therapy;
   treat_1=0; treat_2=0; treat_3=0;
   if txb=1 then treat_1=1;
   if txb=2 then treat_2=1; 
   if txb=3 then treat_3=1;

* Stage;
   if stage=1 then stage_1=0;
   if stage in(2,3) then stage_1=1; 

* PSA;
   if psacat=1 then psa_1=0;
   if psacat in(2,3) then psa_1=1;

* Gleason grade;
   gleason_1=0; gleason_2=0; gleason_3=0;
   if gleason=1 then gleason_1=1;
   if gleason=2 then gleason_2=1;
   if gleason=3 then gleason_3=1;

* Time between diagnosis and baseline;
   time_dx = firstpostdt - dtdx_prostate;
   time_dx_sq = time_dx*time_dx;

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

fixedcov = baseage baseage_spl1
        stage_1
        gleason_2 gleason_3
        psa_1
        treat_2 treat_3
        fhxmi
        act_pre_1 act_pre_2 act_pre_3
        totfv_pre_1 totfv_pre_2 totfv_pre_3
        totgrn_pre_1 totgrn_pre_2 totgrn_pre_3
        junk_pre_1 junk_pre_2 junk_pre_3
        totred_pre_1 totred_pre_2 totred_pre_3
        totproc_pre_1 totproc_pre_2 totproc_pre_3
        sodajuice_pre_1 sodajuice_pre_2 sodajuice_pre_3
        totalc_pre_1 totalc_pre_2 totalc_pre_3
        lcal_pre_1 lcal_pre_2
        lbmi_pre_1 lbmi_pre_2 lbmi_pre_3
        smkhx
        hbp_chl_dia_b
		,

    ncov=11,
	
    cov1  = xcond,  cov1otype   = 2,  cov1ptype  = tsswitch1,
    cov2  = lbmi,   cov2otype   = 3,  cov2ptype  = lag1cat,  cov2knots = 2.99573227 3.2188758 3.4011974,
    cov3  = lcal,   cov3otype   = 3,  cov3ptype  = lag1cat,   cov3knots = 7.49554194  7.82404601,
                    cov3wherem=(not (period in (1, 3, 5, 7, 9))), 
                    cov3wherenosim=(period in (1, 3, 5, 7, 9)), 
                    cov3nosimelsemacro=lcalmacro,
    cov4  = totalc,    cov4otype  = 4,  cov4ptype  = lag1spl,   cov4knots = 0.1  0.5  3,
                    cov4wherem=(not (period in (1, 3, 5, 7, 9))), 
                    cov4wherenosim=(period in (1, 3, 5, 7, 9)), 
                    cov4nosimelsemacro=totalcmacro,
    cov5  = sodajuice, cov5otype  = 4,  cov5ptype  = lag1qdc, 
                    cov5wherem=(not (period in (1, 3, 5, 7, 9))),
                    cov5wherenosim=(period in (1, 3, 5, 7, 9)),
                    cov5nosimelsemacro=sodajuicemacro,
    cov6  = junk,   cov6otype   = 3,  cov6ptype  = lag1qdc,
                    cov6wherem=(not (period in (1, 3, 5, 7, 9))),
                    cov6wherenosim=(period in (1, 3, 5, 7, 9)), 
                    cov6nosimelsemacro=junkmacro,
    cov7  = totproc, cov7otype  = 4,  cov7ptype  = lag1spl,   cov7knots = 0.1  1   6,
                    cov7wherem=(not (period in (1, 3, 5, 7, 9))),
                    cov7wherenosim=(period in (1, 3, 5, 7, 9)), 
                    cov7nosimelsemacro=totprocmacro,
    cov8  = totfv,  cov8otype   = 3,  cov8ptype  = lag1spl,   cov8knots = 2.4  5   9, 
                    cov8wherem=(not (period in (1, 3, 5, 7, 9))),
                    cov8wherenosim=(period in (1, 3, 5, 7, 9)), 
                    cov8nosimelsemacro=totfvmacro,
    cov9  = totred, cov9otype   = 3,   cov9ptype  = lag1spl,   cov9knots = 1  3   7,
                    cov9wherem=(not (period in (1, 3, 5, 7, 9))), 
                    cov9wherenosim=(period in (1, 3, 5, 7, 9)), 
                    cov9nosimelsemacro=totredmacro, 
    cov10  = totgrn, cov10otype   = 3,  cov10ptype  = lag1spl,   cov10knots = 0.5  1.5   3,
                    cov10wherem=(not (period in (1, 3, 5, 7, 9))),
                    cov10wherenosim=(period in (1, 3, 5, 7, 9)),
                    cov10nosimelsemacro=totgrnmacro,
    cov11  = act,    cov11otype   = 3,  cov11ptype  = lag1spl,  cov11knots = 1.5 7.5 65,
    
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
    sample_start = 401,
    sample_end = 500,

    rungraphs=0, 
    printlogstats = 0,

    numint=2
    );
