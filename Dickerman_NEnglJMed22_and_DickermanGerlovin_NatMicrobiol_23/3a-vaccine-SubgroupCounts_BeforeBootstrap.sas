/* 

This script is a part of a set of analytic programs related to the following manuscript(s):
- Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
- Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors.  

Author(s): VA-CAUSAL Methods Core
Version: July 2023. 

*/

proc iml;

submit /R;
	source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/preamble.R", echo=FALSE)
	preamble(server=TRUE)
	counts_summary <- function(treat, event, variant, risk_time=(16*7),subgroup=NULL){
		load(paste0("[SAS_folder_ORD_Project]/Vax_Standard/analysis/results/output-", 
				treat, "-", event, "-exact-", variant,
              	if(!is.null(subgroup)){subgroup},
              	"_1to1.Rda"))
  
		if(variant=="delom_booster") { 
			risk_time <- (16*7) 
		} else if(variant=="omicron_boost") { 
			risk_time <- (9*7) 

		## empiric sample estimates
		if(treat=="PvM"){
		descriptive <- data.frame("n"=output$descriptive[["n"]], 
		                            "bnt-num"=output$descriptive[["treat1"]],
		                            "moderna-num"=output$descriptive[["treat0"]],
		                            "bnt-num-events"=ifelse(is.null(output$descriptive$events_treat1[[as.character(risk_time)]]),"0",output$descriptive$events_treat1[[as.character(risk_time)]] ), 
		                            "moderna-num-events"=ifelse(is.null(output$descriptive$events_treat0[[as.character(risk_time)]]),"0",output$descriptive$events_treat0[[as.character(risk_time)]] )
		                            )
		}

		if(is.null(subgroup)) {
			subname <- "None"
		} else {
			subname <- subgroup
		}

		return(cbind(subname,event,variant,descriptive))
	}

	outcomes <- c("covidpos","covidpossymp","covidhosp","covidICU","COVIDdeath","notcoviddeath")
	subgroups <- c("old70","young70","white","black","diff23_67","diff23_8","diff23_9",
  					"nopriorcov","trueboost","prim_pfizer","prim_moderna")

	across_outs_delom <- function(subg) {
							do.call(rbind,
									lapply(1:length(outcomes),
											function(y) counts_summary(treat="PvM",
																		event=outcomes[y],
																		variant="delom_booster",
																		subgroup=subg)))
									}
	events_by_subgroup_delom <- do.call(rbind,lapply(1:length(subgroups), function(x) across_outs_delom(subg=subgroups[x])))
	delom_main <- across_outs_delom(subg=NULL)

	om_main <- counts_summary(treat="PvM",event="covidpos",variant="omicron_boost") #default is subgroup=NULL

	events_by_subgroup <- rbind(delom_main,events_by_subgroup_delom,om_main) #events_by_subgroup_om #events_by_subgroup_del

	fwrite(events_by_subgroup, 
		file="[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/subgroup_counts_summary.csv");

endsubmit;
run;
