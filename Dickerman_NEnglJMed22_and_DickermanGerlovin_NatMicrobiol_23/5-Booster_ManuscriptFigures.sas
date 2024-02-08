/* 

This script is a part of a set of analytic programs related to the following manuscript(s):
- Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
- Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors.  

Author(s): VA-CAUSAL Methods Core
Version: July 2023. 

*/

/*Figures*/
*NOTE THIS PROGRAM MUST BE EXECUTED FROM MYGSUB OR USING THE RApp ON THE SAS ENTERPRISE GUIDE;

/********************************************************/
/*														*/
/*
** Figure 2 - multi-panel plot of all outcomes
*/
/*														*/
/********************************************************/

proc iml;

submit / R;
source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/preamble.R") #preamble location file
preamble(server=TRUE)
library("cowplot")
source(paste0(prefix,"/functions/plotCumInc_function.R"))

dir <- "[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/" #output directory

a <- plot_function(treat="PvM", event="covidpos", subgroup="none",
                    risk_time=(16*7), yaxis=0.04, ybreak=0.01, yaccuracy=0.5,
                    variant = "delom_booster")

b <- plot_function(treat="PvM", event="covidpossymp", subgroup="none",
                    risk_time=(16*7), yaxis=0.005, ybreak=0.001, yaccuracy=0.01,
                     variant = "delom_booster")

c <- plot_function(treat="PvM", event="covidhosp", subgroup="none",
                    risk_time=(16*7), yaxis=0.005, ybreak=0.001, yaccuracy=0.01,
                    variant = "delom_booster")

d <- plot_function(treat="PvM", event="covidICU", subgroup="none",
                    risk_time=(16*7), yaxis=0.0015, ybreak=0.0005, yaccuracy=0.025,
                    variant = "delom_booster")

e <- plot_function(treat="PvM", event="COVIDdeath", subgroup="none",
                    risk_time=(16*7), yaxis=0.0015, ybreak=0.0005, yaccuracy=0.025,
                    variant = "delom_booster")

plotter2 <- cowplot::plot_grid(a,b,c,d,e,labels="AUTO",label_size=10, ncol=2)

  ggsave(filename=paste0(dir,"Fig2.pdf"),
         width=5.9,height=8.2,
         plot=plotter2)

endsubmit;
run;


/********************************************************/
/*														*/
/*
** Negative Control Plot 1 (Figure S2) - Covid positive symptomatic early
*/
/*														*/
/********************************************************/

proc iml;

submit / R;
source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/preamble.R")
preamble(server=TRUE)
source(paste0(prefix,"/functions/plotCumInc_function.R"))

dir <- "[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/"

neg_control_1 <-  plot_function(treat="PvM", event="covidpossymp", subgroup="none",
                risk_time=(1*7), yaxis=0.005, ybreak=0.001, yaccuracy=0.025,
                variant = "delom_booster")

ggsave(filename=paste0(dir,"ED_Fig2.pdf"),
       width=3.4,height=2.4,
       plot=neg_control_1)

endsubmit;
run;

/********************************************************/
/*														*/
/*
** Negative Control Plot 2 (Figure S3) - Non-Covid Death
*/
/*														*/
/********************************************************/

proc iml;

submit / R;
source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/preamble.R")
preamble(server=TRUE)
source(paste0(prefix,"/functions/plotCumInc_function.R"))

dir <- "[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/"

neg_control_2 <-  plot_function(treat="PvM", event="notcoviddeath", subgroup="none",
                risk_time=(16*7), yaxis=0.01, ybreak=0.002, yaccuracy=0.05,
                variant = "delom_booster")

ggsave(filename=paste0(dir,"ED_Fig3.pdf"),
       width=3.4,height=2.4,
       plot=neg_control_2)

endsubmit;
run;


/********************************************************/
/*														*/
/*
** Omicron-only sensitivity analysis (Covid positive)
*/
/*														*/
/********************************************************/

proc iml;

submit / R;
source("[SAS_folder_ORD_Project]/Vax_Standard/analysis/preamble.R")
preamble(server=TRUE)
source(paste0(prefix,"/functions/plotCumInc_function.R"))

dir <- "[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/"

omicron_covidpos <-  plot_function(treat="PvM", event="covidpos", subgroup="none",
                                risk_time=(9*7), yaxis=0.03, ybreak=0.01, yaccuracy=0.5,
                                variant = "omicron_boost")

ggsave(filename=paste0(dir,"Fig3.pdf"),
       width=3.4,height=2.4,
       plot=omicron_covidpos)

endsubmit;
run;

/********************************************************/
/*														*/
/*
** LOVE PLOTS
*/
/*														*/
/********************************************************/

libname sasin '[SAS_folder_ORD_Project]/Vax_Standard/analysis/data';

/*proc contents data=sasin.matched_PM_omboost; run;*/
/*proc contents data=sasin.matched_PM_delom; run;*/


proc iml;

run ExportDataSettoR("sasin.matched_PM_delom","matched_PM_delom");
run ExportDataSettoR("sasin.matched_PM_omboost","matched_PM_omboost");
submit / R;
library("cowplot")

pacman::p_load(MatchIt,ggQC,cobalt,ggplot2,tidyr,viridis,extrafont,dplyr,survival,survminer,tableone)
names(matched_PM_delom) <- tolower(names(matched_PM_delom))
names(matched_PM_omboost) <- tolower(names(matched_PM_omboost))

new.names <- c(age_at_index = "Age",
				male = "Male", female = "Female",
				white = "Race, White", black = "Race, Black",
			   urbanicity = "Urban residence",
			   dose3_dt = "Date of dose 3",
			   monthdiff23 = "Months between doses 2 and 3",
			   tot_tests = "Number of tests in the year before dose 3",
			   current_smoker ="Current smoker",
			   cond_lung = "Chronic lung disease",
			   cond_vasc = "Cardiovascular disease",
			   cond_hypt = "Hypertension",
			   cond_diab = "Diabetes",
			   cond_ckd = "Chronic kidney disease",
				cond_cld = "Chronic liver disease", 
				cond_5yr_cancer = "Cancer",
				immunsupp = "Immunocompromised state",
				obesity = "Obesity",
				cond_dement = "Dementia", 
				cond_sud = "Substance use disorder") 

#Function to make the variables look nice when it comes to plotting and formatting (apply to the data before running matchit)
prep_loveplot <- function(matched_dat){

	matched_dat2 <- matched_dat %>% 
					mutate(dose3 = factor(group_binary, levels=c(1,0), labels=c('BNT162b2','mRNA-1273'))) %>%
					mutate(current_smoker = case_when(smoking_status==1 ~ 1,
					TRUE ~ 0)) %>%
					mutate(same_primary = case_when(dose2_man==dose3_man ~ 1, 
					TRUE ~ 0)) %>% # variable denoting whether the primary manufacturer matches dose 3 (01mar2022 HG)
					#mutate(sex = factor(sex, levels=c(1,0), labels=c('Male','Female'))) %>%
					mutate(age_cat = factor(age_cat, levels = c(0,1,2,3,4,5),
											labels = c('18 to 39 years','40 to 49 years','50 to 59 years',
														'60 to 69 years', '70 to 79 years','>=80 years'))) %>%
					mutate(male = case_when(sex==1 ~ 1, TRUE ~ 0)) %>%
					mutate(female = case_when(sex==0 ~ 1, TRUE ~ 0)) %>%
					mutate(black = case_when(race==1 ~ 1, TRUE ~ 0)) %>%
					mutate(white = case_when(race==0 ~ 1, TRUE ~ 0)) %>%
					mutate(race = factor(race, levels=c(1,2,3,0), labels=c('Black','Other','Unknown','White'))) %>%
					mutate(ethnicity = factor(ethnicity, levels=c(1,0,2), labels=c('Hispanic','Not Hispanic','Unknown'))) %>%
					mutate(smoking_status = factor(smoking_status, levels=c(1,2,0), labels=c('Current','Former','Never'))) %>%
					mutate(January = case_when(dose2_month==1 ~ 1, TRUE ~ 0)) %>%
					mutate(February = case_when(dose2_month==2 ~ 1, TRUE ~ 0)) %>%
					mutate(March = case_when(dose2_month==3 ~ 1, TRUE ~ 0)) %>%
					mutate(April = case_when(dose2_month==4 ~ 1, TRUE ~ 0)) %>%
					mutate(May = case_when(dose2_month==5 ~ 1, TRUE ~ 0)) %>%
					mutate(June = case_when(dose2_month==6 ~ 1, TRUE ~ 0)) %>%
					mutate(July = case_when(dose2_month==7 ~ 1, TRUE ~ 0)) %>%
					mutate(August = case_when(dose2_month==8 ~ 1, TRUE ~ 0)) %>%
					mutate(September = case_when(dose2_month==9 ~ 1, TRUE ~ 0)) %>%
					mutate(October = case_when(dose2_month==10 ~ 1, TRUE ~ 0)) %>%
					mutate(November = case_when(dose2_month==11 ~ 1, TRUE ~ 0)) %>%
					mutate(December = case_when(dose2_month==12 ~ 1, TRUE ~ 0)) %>%
					mutate(dose2_month = factor(dose2_month, levels=c(1,2,3,4,5,6,7,8,9,10,11,12),
								labels=c('January','February','March','April','May','June','July','August','September','October','November','December'))) %>%
					mutate(Tests0 = case_when(ntests==0 ~ 1, TRUE ~ 0)) %>%
					mutate(Tests1 = case_when(ntests==1 ~ 1, TRUE ~ 0)) %>%
					mutate(Tests2 = case_when(ntests==2 ~ 1, TRUE ~ 0)) %>%
					mutate(ntests = factor(ntests, levels=c(0,1,2), 
								labels=c('0 tests','1 test','>=2 tests'))) %>%
					mutate(diff23_3grp = factor(diff23_3grp, levels=c(1,2,3), labels=c('6-<8','8-<9','>=9'))) %>%
					mutate(pcp5 = factor(pcp5, levels=c(0,1,2,3), labels=c('1-9','10-19','20-29','>=30'))) %>%
					mutate(flu5 = factor(flu5, levels=c(0,1,2,3), labels=c('0 vaccinations','1 or 2','3 or 4','>=5'))) %>%
					mutate(BNT162b2 = case_when(dose2_man=='Pfizer' ~ 1, TRUE ~ 0)) %>%
					mutate(mRNA1273 = case_when(dose2_man=='Moderna' ~ 1, TRUE ~ 0)) %>%
					mutate(nopriorcov = case_when(covidpriordose3==1 ~ 0, TRUE ~ 1))



	#love plot values
	m.out1 <- matchit(formula = as.formula(group_binary ~ age_at_index + male + female + white + black + urbanicity + 
	dose3_dt + monthdiff23 + tot_tests + 
	current_smoker + cond_lung + cond_vasc + cond_hypt + cond_diab + cond_ckd + 
	cond_cld + cond_5yr_cancer + immunsupp + obesity + cond_dement + cond_sud) , 
	                  data=matched_dat2, method = NULL) 

	return(m.out1)

}


dir <- "[SAS_folder_ORD_Project]/Vax_Standard/BoosterManuscript/"

delom_love <- prep_loveplot(matched_dat=matched_PM_delom)
omboost_love <- prep_loveplot(matched_dat=matched_PM_omboost)

fig_sd_delom <- bal.tab(delom_love)$Balance
write.csv(fig_sd_delom, file=paste0(dir,"SD_ED_Fig1A.csv"))

fig_sd_omboost <- bal.tab(omboost_love)$Balance
write.csv(fig_sd_omboost, file=paste0(dir,"SD_ED_Fig1B.csv"))


  ### PLOT ###
fig_love_delom <-  love.plot(delom_love, 
            binary="raw",
            stars="std",
            drop.distance = TRUE,
            stats=c("mean.diffs"),
            #var.order="unadjusted",
            thresholds=c(m=0.1),
            limits=list(m=c(-0.12,0.12)),
            var.names=new.names,
            position="none",
			title="Matched Population from Primary Analysis during a Period \n Spanning Delta- and Omicron-Variant Predominance",
            themes = list(m = theme(text = element_text(size=9,colour="black"), #family="Arial",
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                    panel.border = element_blank(),axis.line = element_line(colour = "black"))))

fig_love_omboost <-  love.plot(omboost_love, 
            binary="raw",
            stars="std",
            drop.distance = TRUE,
            stats=c("mean.diffs"),
            #var.order="unadjusted",
            thresholds=c(m=0.1),
            limits=list(m=c(-0.12,0.12)),
            var.names=new.names,
            position="none",
			title="Matched Population from Secondary Analysis \n during a Period of Omicron-Variant Predominance",
            themes = list(m = theme(text = element_text(size=9,colour="black"), #family="Arial",
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                    panel.border = element_blank(),axis.line = element_line(colour = "black"))))

plotter <- cowplot::plot_grid(fig_love_delom,fig_love_omboost,labels="AUTO",label_size=10, ncol=1)



#save covariate balance (love) plot (pdf)
ggsave(plotter, 
       file=paste0(dir,"ED_Fig1.pdf"),
       width=5.9, height=8.2)


endsubmit;
quit;
