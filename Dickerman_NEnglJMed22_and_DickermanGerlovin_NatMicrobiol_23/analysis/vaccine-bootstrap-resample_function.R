#####
#
# This script is a part of a set of analytic programs related to the following manuscript(s):
#  - Dickerman BA, Gerlovin H, Madenci AL, Kurgansky KE, Ferolito BR, Figueroa Muniz MJ, et al. Comparative effectiveness of BNT162b2 and mRNA-1273 vaccines in U.S.veterans. N Engl J Med. 2022;386(2):105-15.
#  - Dickerman BA*, Gerlovin H*, Madenci AL, Figueroa Muniz MJ, Wise JK, Adhikari N, et al. Comparative effectiveness of third doses of mRNA-based COVID-19 vaccines in US veterans. Nat Microbiol. 2023 Jan;8(1):55-63. *Co-first authors. 
#
# Author(s): VA-CAUSAL Methods Core
# Version: July 2023.
#
######

# bootstrap resample


bootstrap_resample_function <- function(treat, event, variant="alpha", 
                                        subgroup=NULL, number.rows=Inf,		
                                        seed, N=length(unique(dat$newid))) {

  set.seed(seed)
  ifelse(!dir.exists(file.path(prefix, "bootstrap")), dir.create(file.path(prefix, "bootstrap")), FALSE)
  
  dat <- fread(file=paste0(prefix, "/data/cleaned-dat-", treat,"-",event,"-",variant,".csv"), nrows=number.rows)[!is.na(Age_at_index)] 

  if(!is.null(subgroup)){
    #subgroup <- paste0("-",subgroup)				
    if(subgroup=="white"){dat <- dat[race==0,]}
    if(subgroup=="black"){dat <- dat[race==1,]}
    if(subgroup=="old70"){dat <- dat[age>=70,]}
    if(subgroup=="young70"){dat <- dat[age<70,]}
      if(subgroup=="diff23_67"){dat <- dat[diff23_3grp==1,]}  
      if(subgroup=="diff23_8"){dat <- dat[diff23_3grp==2,]}  
      if(subgroup=="diff23_9"){dat <- dat[diff23_3grp==3,]}  
      if(subgroup=="nopriorcov"){dat <- dat[covidpriordose3==0,]}  
      if(subgroup=="trueboost"){dat <- dat[dose3_booster==1,]}  
      if(subgroup=="prim_pfizer"){dat <- dat[dose2_man=="Pfizer",]}  
      if(subgroup=="prim_moderna"){dat <- dat[dose2_man=="Moderna",]}  
  }
  bs.ids <- sample(unique(as.numeric(dat$newid)), size=N, replace = TRUE)
  bs.ids.df <- left_join(data.frame(newid=unique(bs.ids)), 
                          data.frame(table(bs.ids)) %>% mutate(newid=as.numeric(as.character(bs.ids))) %>% dplyr::select(-bs.ids),
                          by="newid"
  )
  
  bs.dat <- left_join(bs.ids.df, 
                   dat, by="newid")
  
  bs.dat$bsid <- 1:nrow(bs.dat)
  bs.dat <- setDT(bs.dat)[rep(seq(.N), Freq)]  
  # options(warn=oldw)
  return(bs.dat)
}

