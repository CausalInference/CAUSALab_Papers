/* 

This script is a part of a set of analytic programs related to the following manuscript(s):
- Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
- Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors.  

Author(s): VA-CAUSAL Methods Core
Version: July 2023. 

*/

proc iml;
	%let event0 = %scan(&SYSPARM.,1,^); %put &event0.;
	event1 = "&event0"; %put event1; *copies to a character string;

	%let treat0 = %scan(&SYSPARM.,2,^); %put &treat0.;
	treat1 = "&treat0"; %put treat1; *copies to a character string;

	bsnum0 = %scan(&SYSPARM.,3,^);
	cores0 = %scan(&SYSPARM.,4,^);

	%let subgroup0 = %scan(&SYSPARM.,5,^); %put &subgroup0.; 
	subgroup1 = "&subgroup0"; %put subgroup1; *copies to a character string; 

	%let variant0 = %scan(&SYSPARM.,6,^); %put &variant0.; 
	variant1 = "&variant0"; %put variant1; 

	if variant1="delom_booster" then do; risktime = 112; daysfu = 118; end; 
	else if variant1="omicron_boost" then do; risktime = 63; daysfu = 65; end;
	else if variant1="alpha" then do; risktime = 168; daysfu = 178; end;
	else if variant1="delta" then do; risktime = 83; daysfu = 83; end;

	%put risktime; %put daysfu;

submit event1 treat1 bsnum0 cores0 subgroup1 variant1 risktime daysfu/ R;
	event <- "&event1"
	treat <- "&treat1"
	bsnum <- "&bsnum0"
	cores <- "&cores0" 
	variant <- "&variant1" 
	subgroup <- "&subgroup1" 
	risk_time <- "&risktime" 
	days_fu <- "&daysfu" 

	source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/preamble.R", echo=TRUE)
	source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/vaccine-analysis_function.R", echo=TRUE) #PvM analyses
  	source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/vaccine-bootstrap-wrapper_function.R", echo=TRUE) #PvM analyses

	if(variant %in% c("delom_booster","omicron_boost")) {
		matchform <- c("Age_at_index + sex + race + VISN + Caldate_num + urbanicity + dose2_month + ntests")
	} else if(variant %in% c("alpha","delta","safety")){
		matchform <- c("Age_at_index + sex + race + VISN + Caldate_num + urbanicity")
	}

	wrapper_function(treat=treat, event=event, number.rows=Inf, subgroup=subgroup, 
						quickcheck=FALSE, bsnum=bsnum, cores=cores, total.resamples=500, 
						matching.formula = paste0("group.binary ~ ",matchform), 
						variant = variant,
						days_fu=as.numeric(days_fu), 
						risk_time=as.numeric(risk_time))
endsubmit;
quit;

