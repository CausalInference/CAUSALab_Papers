/* 

This script is a part of a set of analytic programs related to the following manuscript(s):
- Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
- Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors.  

Author(s): VA-CAUSAL Methods Core
Version: July 2023. 

*/

/*iml to read in each separate dataset for vaccines study*/
libname vboost '[SAS_folder_ORD_Project]/Data/booster'; *files for boosters analysis;
libname vdata '[SAS_folder_ORD_Project]/Data'; *alpha and delta datasets;

proc iml;
	libname vboost '[SAS_folder_ORD_Project]/Data/booster'; 
	libname vdata '[SAS_folder_ORD_Project]/Data';

	* Grab the event of interest;
	%let event0 = %scan(&SYSPARM.,1,^); %put &event0.;
	event1 = "&event0"; %put event1; *copies to a character string;

	* Define the treatment contrast to be used;
	%let treat0 = %scan(&SYSPARM.,2,^); %put &treat0.;
	treat1 = "&treat0"; %put treat1; *copies to a character string;

	* Define the variant/timeframe/analysis;
	%let variant0 = %scan(&SYSPARM.,3,^); %put &variant0.; 
	variant1 = "&variant0"; %put variant1; *copies to a character string; 

	*Pull datasets into process based on event, contrast, and analysis of interest;
	if variant1="omicron_boost" then do;
		run ExportDataSettoR("vboost.finaldata_om_&event0","dat"); *omicron-only boosters analysis;
	end;
	else if variant1="delom_booster" then do;
		run ExportDataSettoR("vboost.finaldata_delom_&event0","dat"); *primary delta-omicron boosters analysis;
	end;
	else if variant1="alpha" then do;
		run ExportDataSettoR("vdata.finaldata_&event0","dat"); *alpha-period dataset for first paper;
	end;
	else if variant1="delta" then do;
		run ExportDataSettoR("vdata.finaldata_&event0._delta","dat"); *delta-period dataset for first paper;
	end;

submit event1 treat1 variant1/ R;
	event <- "&event1"
	treat <- "&treat1"
	variant <- "&variant1"
	source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/preamble.R", echo=TRUE) 
	source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/vaccine-data-prep_function.R", echo=TRUE) 
	if(variant %in% c("delom_booster","omicron_boost")) {
		matchform <- c("Age_at_index + sex + race + VISN + Caldate_num + urbanicity + dose2_month + ntests")
	} else if(variant %in% c("alpha","delta","safety")){
		matchform <- c("Age_at_index + sex + race + VISN + Caldate_num + urbanicity")
	}
	data_prep(dat=dat, event=event, treat=treat, variant=variant, server=TRUE, 
			matching.formula = paste0("group.binary ~ ",matchform), 
			age.bin = 5, date.bin = 5 
	)
endsubmit; quit;

