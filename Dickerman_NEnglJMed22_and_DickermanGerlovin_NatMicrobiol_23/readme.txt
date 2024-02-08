/*
This is a README outlining the analytic programs used to evaluate the comparative effectiveness of mRNA-based Covid-19 vaccines in the following manuscripts:
- Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
- Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors.  

It is assumed that the input dataset for these programs contains all necessary analytic variables, as described in the accompanying codebook. Certain programs in this package are flexible enough to accommodate applications with different goals (e.g., comparative effectiveness of primary series, third doses); our examples for their invocations focus on one potential application (comparative effectiveness of third doses) for simplicity.  

Author(s): VA-CAUSAL Methods Core
Version: July 2023. 

*/ 

*********************
ANALYSIS
*********************

"vaccine-main-script.sas" is a single script to
	1) prepare data, 
	2) (a) run negative control analyses and other preliminaries, (b) generate point estimates 
	3) (a) check counts then (b) perform bootstrap resample,
	4) prepare summary data for tables and figures
	and (5) generate confidence intervals, Kaplan-Meier plots, and Manuscript tables 

Across the inner R programs, the following arguments persist:
- "treat" (e.g., PvM, MvP) 
- "event" (e.g., covidpos, covidpossymp, covidhosp, covidICU, COVIDdeath, notcoviddeath)
- "method" takes either "exact" (N:N exact matching) or "cem" (coarsened exact matching)	
- "variant" is used to specify the four different time-period analyses included. This is used for naming outputted logs, reading in cleaned datasets, naming saved matched data, naming negative control KM plot, and assigning default values for the inner scripts. Possible values are "alpha", "delta", "delom_booster", and "omicron_boost"

- "matching formula" is specified based on the analysis (see "variant")
* "alpha" and "delta" analyses:  group.binary ~ Age_at_index + sex + race + VISN + Caldate_num + urbanicity
* "delom_booster" and "omicron_boost" analyses:  group.binary ~ Age_at_index + sex + race + VISN + Caldate_num + urbanicity + dose2_month + ntests

- "risk_time" takes number of days to be considered for computing event counts, follow-up statistics (default 168 days = 24 weeks = approximately 6 months). Different defaults are set in the SAS kick-off programs based on analysis. 	
* "alpha" analyses: 168 days (24*7)
* "delta" analyses: 83 days
* "delom_booster" analyses: 112 days (16*7)
* "omicron_boost" analyses: 63 days (9*7)
		
- "days_fu" takes number of days to be considered as negative outcome control in the preliminaries program, and maximum number of days in the vaccine_analysis function. Different defaults are set in the SAS kick-off programs based on analysis and R function used (prelim_function vs. vaccine_function):
* "alpha" analyses: prelim_function(days_fu=10), vaccine_function(days_fu=178)
* "delta" analyses: prelim_function(days_fu=10), vaccine_function(days_fu=83)
* "delom_booster" analyses: prelim_function(days_fu=7), vaccine_function(days_fu=118)
* "omicron_boost" analyses: prelim_function(days_fu=7), vaccine_function(days_fu=65)



1. macro "cleandat" runs script "1-vaccine-prep-data.sas" (which calls "vacccine-data-prep_function.R") on the server.		
- Function: Prepares separate datasets for each of the events of interest. The data preparation coarsens continuous variables and only includes relevant eligible trials.

After running through the data preparation, we checked the follow-up periods for the earliest possible outcome (min(eof,covidpos)) and latest possible outcome (min(eof,coviddeath)) to set the parameters for risk time and days of follow-up in the subsequent analytic programs.


2a. macro "prelim" runs script "2a-vaccine-prelim.sas" (which calls "vaccine-prelim_function.R") on the server
- Function: Accomplish the following tasks, pre-analysis:
	- generate negative control plots for visual assessment of whether nonparametrically estimated cumulative incidence curves overlap, to determine final set of matching factors
	- output matched dataset for each comparison (e.g., PvM) for later creation of Table 1 and covariate balance plots
	- output event counts and follow-up time stats (for reporting and to determine which analyses to proceed with - must have >10 events) *Note this last step is done only in the log file. 


2b. macro "point_est" runs script "2-vaccine-analysis.sas" (which calls "vaccine-analysis_function.R") on the server
			
- options "KM_analysis" (default) to estimate nonparametrically 
- option "negcontrol" (default) outputs non-parametric 10 day point estimate
- results saved in ~/analysis/results folder as .Rda
- argument "subgroup" is NULL
- outputs counts and follow-up time descriptives into dataset for later use (same values as in the log from 2a-vaccine-prelim)


2c. macro "point_est_subgroup" runs script "2-vaccine-analysis-subgroup.sas" (which calls "vaccine-analysis_function.R") on the server
- same as #2b above except "point_est_subgroup" macro also takes argument "subgroup" for subsetting data for subgroup analyses, naming saved risk estimates at each time t for later plots, naming outputted results

subgroup values:
	* white and black -- based on race
	* old70 and young70 -- groups of individuals 70 years and older (old70) or younger than 70 (young70)
	* diff23_67, diff23_8, diff23_9 -- based on dose separation between dose 2 and 3 (in days, approximating months), 6-7 months (diff23_67: diff23_3grp=1), 8 months (diff23_8: diff23_3grp=2), 9+ months (diff23_9: diff23_3grp=3)
	* nopriorcov -- those with no previous covid positive statuses
	* trueboost -- sensitivity analysis where dose3 identified as a booster dose (dose3_booster=1)
	* prim_pfizer and prim_moderna -- analyses based on the primary series manufacturer to allow for comparison of moderna and pfizer third doses within each primary series (pfizer: prim_pfizer; moderna: prim_moderna)


3. Bootstrapping process
3a) Check the counts for the subgroups prior to running any bootstrapping. 
	* Program used for Boosters manuscript: ~/3a-vaccine-SubgroupCounts_BeforeBootstrap.sas

3b) macro "bs_mac" runs script "3-vaccine-bootstrap.sas" (which calls "wrapper_function.R") on the server
- outputs non-parametric estimates (by default)
- option "quickcheck" runs 2 resamples per core to test the function
- default number of resamples is 500 (can be changed with "total.resamples" option)
- estimates from each bootstrap resample saved separately in ~/analysis/bootstrap as .Rda files
- subgroup = "none" allows the user to run the same bs_mac for the non-subgroup analyses


4. Data post-processing and preparation for Manuscripts
	* Program used for Boosters manuscript: ~/4-Booster_PrepTablesFiguresData.sas
	* Combines the final datasets from pre-step 1 (to get extra baseline variables) and the output dataset for matched and eligible individuals used in the Covid-19 positive and Covid-19 death outcomes. 

5. Programs to combine results across analyses, variations, and bootstrapping files for use in the Manuscript Tables and Figures. Both programs can be initialized through the %mygsub() commands inside of the vaccine-main-script.sas, or can be manually executed from a SAS EG connection to the RApp server.

5-Booster_ManuscriptFigures.sas
- calls plot_function() to create Kaplan-Meier plots with confidence bounds from the bootstrapped and point estimate data. This is done for the overall results and negative control plots.
	- argument "total.number" specifies number of bootstrap samples to pool (default 500)
	- argument "risk_time" specifies time at which you want 95% CI around estimated risk, risk ratio, risk difference
	- argument "yaxis" specifies height of y-axis of KM plot
	- argument "ybreak" specifies breaks for y-axis of KM plot
	- argument "yaccuracy" specifies accuracy for y-axis labels 
	- argument "yname" specifies label for y-axis of KM plot 
- cowplot package used to allow for multi-panel plotting
- results saved in ~/BoosterManuscript folder
- Love (balance) plots created using matched datasets from step 4 and matchIt package

5-Booster_ManuscriptTables.sas
- The purpose of this file is to aggregate the results into table 2 for the manuscript and produce the tables of baseline characteristics by vaccine manufacturer in the matched and eligible populations for both primary and secondary analyses. 
- loads the table2-function program from within the ~/analysis/functions/ folder. This is described in more detail, below, under "Other" section
		

Other: 
preamble.R
	* Defines the prefix folder paths for analyses on the SASGrid server (SERVER=TRUE) vs. P drive (SERVER=FALSE)
	* loads the necessary packages needed in the R environment

strip-glm_function.R
	* Trims model output when running PLR in vaccine-analysis_function.R

table2-function.R
	* defines table_2_function, which loads local preamble (server=FALSE), and some functions from the ~/analysis/functions folder (bs-results-summary_function.R, risk-point-ests_function.R, and rd_rr_function.R). More detail on each, below. 
	* function pulls together the summary point estimates and bootstrapped bounds, as well as, descriptive statistics for a particular treat, event, variant, subgroup, and risk_time specification.
	* also defines the parallel_table_2_function which allows for parallelization of the six outcomes (event) across a given treat, variant, risk_time, and subgroup specification.

bs-results-summary_function.R
	* pulls in the output files (~/results/output...) and bootstraps (~/bootstrap/output) for a given treat, event, variant, risk_time, and subgroup specification
	* defaults to 500 resamples
	* descriptives pulled from output (n, bnt-num-events, moderna-num-events) 
	* confidence intervals based on the 97.5 and 2.5 quantiles

risk-point-ests_function.R
	* parses through dataset to capture point estimates based on particular contrast selected and calculates the risks per 10,000 persons

rd_rr_function.R
	* extracts risk ratios and risk differences from the input dataset
	* called from within table_2_function

Notes:
- logs for each R script saved in ~/r-logs folder
