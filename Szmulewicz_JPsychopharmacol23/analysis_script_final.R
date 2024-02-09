#################################################################################################################
# CAUSALab 
# FILE: analysis_script_final
# PROJECT: Estimating the per-protocol effect of lithium on suicidality in a randomized trial of individuals with
#          major depression or bipolar disorder
# AUTHORS: Alejandro Szmulewicz and Arin Madenci
# DATE: 23 February 2023
# MANUSCRIPT: Estimating the per-protocol effect of lithium on suicidality in a randomized trial of individuals with
#             major depression or bipolar disorder, Journal of Psychopharmalogy, 2023. 
#################################################################################################################


# The code below reproduces the tables with main findings in the manuscript: Table 2 and eTable 5
# Each row of those tables represents a separate analysis with its own data cleaning file.
# but analyses require the same function file (analysis_function_final)

######################################################
## ROW 1: Intention-to-treat analysis
######################################################
# data cleaning
source(".../01-intention-to-treat analysis/DataCleaning_itt.R")
# functions
source(".../analysis_function_final.R")

formula_y_d = "PO==0 ~ visit.week + visit.week2 + lithium + I(visit.week*lithium) + I(visit.week2*lithium)"

analysis_function(.dat=finaldata, .formula_y_d=formula_y_d, .formula_t_d=NULL, .formula_t_n=NULL, bootstrap=FALSE)

bs_ests <- pbapply::pblapply(1:2000, function(x){
  analysis_function(.dat=finaldata, .formula_y_d=formula_y_d,bootstrap=TRUE,.formula_t_d=NULL, .formula_t_n=NULL, bs_num=x)
})

apply(do.call(rbind, bs_ests), 2, quantile, probs=c(0.025, 0.975))


##########################################################
## ROW 2: Intention-to-treat analysis restricted to those
##        with complete baseline covariates
#########################################################
# data cleaning
source(".../DataCleaning_ittcompletebasecovs.R")
# functions
source(".../analysis_function_final.R")

formula_y_d = "PO==0 ~ visit.week + visit.week2 + lithium + I(visit.week*lithium) + I(visit.week2*lithium)"

analysis_function(.dat=finaldata, .formula_y_d=formula_y_d, .formula_t_d=NULL, .formula_t_n=NULL, bootstrap=FALSE)

bs_ests <- pbapply::pblapply(1:2000, function(x){
  analysis_function(.dat=finaldata, .formula_y_d=formula_y_d,bootstrap=TRUE,.formula_t_d=NULL, .formula_t_n=NULL, bs_num=x)
})

apply(do.call(rbind, bs_ests), 2, quantile, probs=c(0.025, 0.975))

##########################################################
## ROW 3: Restricted intention-to-treat analysis with
##        censoring at dropout
#########################################################
# data cleaning
source(".../DataCleaning_censoringatdropout.R")
# functions
source(".../analysis_function_final.R")

# 3.1. Unadjusted
formula_y_d = "PO==0 ~ visit.week + visit.week2 + lithium + I(visit.week*lithium) + I(visit.week2*lithium)"

analysis_function(.dat=finaldata, .formula_y_d=formula_y_d, .formula_t_d=NULL, .formula_t_n=NULL, bootstrap=FALSE)

bs_ests <- pbapply::pblapply(1:2000, function(x){
  analysis_function(.dat=finaldata, .formula_y_d=formula_y_d,bootstrap=TRUE,.formula_t_d=NULL, .formula_t_n=NULL, bs_num=x)
})

apply(do.call(rbind, bs_ests), 2, quantile, probs=c(0.025, 0.975))

# 3.2. Baseline adjusted:
formula_y_d = "PO==0 ~ visit.week + visit.week2 + lithium + I(visit.week*lithium) + I(visit.week2*lithium) +
                          SEX + age + I(age*age) + White +
                         as.factor(prior_suicide_cat) + has_bipolar_disorder + columbia_ideation_bas + PTSD +
                         as.factor(baseline_dep_cat)"

analysis_function(.dat=finaldata, .formula_y_d=formula_y_d, .formula_t_d=NULL, .formula_t_n=NULL, bootstrap=FALSE)

bs_ests <- pbapply::pblapply(1:2000, function(x){
  analysis_function(.dat=finaldata, .formula_y_d=formula_y_d,bootstrap=TRUE,.formula_t_d=NULL, .formula_t_n=NULL, bs_num=x)
})

apply(do.call(rbind, bs_ests), 2, quantile, probs=c(0.025, 0.975))

##########################################################
## ROW 4: Per-protocol analysis
#########################################################
# data cleaning
options(warn=-1)
source(".../DataCleaning_ipw.R")
# functions
source(".../analysis_function_final.R")

# 4.1. Unadjusted
formula_y_d = "PO==0 ~ visit.week + visit.week2 + lithium"
analysis_function(.dat=finaldata, .formula_y_d=formula_y_d, .formula_t_d=NULL, .formula_t_n=NULL, bootstrap=FALSE)
bs_ests <- pbapply::pblapply(1:2000, function(x){
  analysis_function(.dat=finaldata, .formula_y_d=formula_y_d, .formula_t_d=NULL, .formula_t_n=NULL, bootstrap=TRUE, bs_num=x)
})

apply(do.call(rbind, bs_ests), 2, quantile, probs=c(0.025, 0.975))

# 4.2. Baseline adjusted
analysis_function(.dat=finaldata, 
                  .formula_y_d="PO==0 ~ visit.week + visit.week2 + lithium + White + education + as.factor(prior_suicide_cat) + age+I(age*age)+
                         has_bipolar_disorder + columbia_ideation_bas + PTSD + PERS + 
                         SEX+marital+subs_abuse + 
                         I(lithium*visit.week) + I(lithium*visit.week2) +
                         as.factor(baseline_dep_cat2)", 
                  .formula_t_d=NULL, .formula_t_n=NULL, bootstrap=FALSE)

bs_ests <- pbapply::pblapply(1:2000, function(x){
  analysis_function(.dat=finaldata, 
                    .formula_y_d="PO==0 ~ visit.week + visit.week2 + lithium + White + education + as.factor(prior_suicide_cat) + age+I(age*age)+
                         has_bipolar_disorder + columbia_ideation_bas + PTSD + PERS + 
                         SEX+marital+subs_abuse + 
                         I(lithium*visit.week) + I(lithium*visit.week2) +
                         as.factor(baseline_dep_cat2)", 
                    .formula_t_d=NULL, .formula_t_n=NULL, bootstrap=TRUE, bs_num=x)
})

apply(do.call(rbind, bs_ests), 2, quantile, probs=c(0.025, 0.975))

# 4.3. Baseline and post-baseline adjustment
formula_y_d = "PO==0 ~visit.week + visit.week2 + lithium + 
                        I(lithium*visit.week) + I(lithium*visit.week2) +
                        as.factor(prior_suicide_cat) + has_bipolar_disorder + columbia_ideation_bas  + PTSD"

formula_t_d <- "as.factor(non_adherence==0) ~ lithium*(visit.week + visit.week2 +  White + education + as.factor(prior_suicide_cat) + 
                   has_bipolar_disorder + columbia_ideation_bas + PTSD + PERS +age + I(age*age) +
                   SEX + marital + subs_abuse + lag_antipsy + er.visit +
                   baseline_dep_cat2 +
                   change_ideation_lagged_cat + depression_change_lagged + 
                   I(depression_change_lagged*depression_change_lagged) +
                   side.eff_bin)"

formula_t_n <- "as.factor(non_adherence==0) ~ lithium*(visit.week + visit.week2 +  #interaction with lithium instead of splitting dataset
                                                                  as.factor(prior_suicide_cat) + has_bipolar_disorder + columbia_ideation_bas + PTSD)"

analysis_function(.dat=finaldata, .formula_y_d=formula_y_d, .formula_t_d=formula_t_d, .formula_t_n=formula_t_n, bootstrap=FALSE)
bs_ests <- pbapply::pblapply(1:2000, function(x){
  analysis_function(.dat=finaldata, .formula_y_d=formula_y_d, .formula_t_d=formula_t_d, .formula_t_n=formula_t_n, bootstrap=TRUE, bs_num=x)
})

apply(do.call(rbind, bs_ests), 2, quantile, probs=c(0.025, 0.975))


