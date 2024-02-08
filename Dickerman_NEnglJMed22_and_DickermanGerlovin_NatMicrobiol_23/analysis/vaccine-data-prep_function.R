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

data_prep <- function(dat, treat, event,
                      matching.formula = paste0("group.binary ~ Age_at_index + sex + race + VISN + Caldate_num + urbanicity"), 
                      age.bin = 5, date.bin = 5, variant="alpha", ...){ 
  # Load server/filepath specifications
  preamble(...)
  
  # Initialize log
  if(server){
    my_log <<- file(paste0(prefix, "/r-logs/log-",treat,"-",event,"-",variant,".txt"))
    sink(my_log, append=TRUE)
    sink(my_log, append=TRUE, type="message")
  }

  print(variant)
  print(matching.formula)

  # Keep only the variables needed in subgroup (stratified) analyses and for matching purposes 
  keep.vars = unique(c("newid","outcome","newperiod","group.binary","eligible","days","race","urbanicity",
						"smoking_status","pcp5","flu5","age",
						"diff23_3grp","dose3_man","dose2_man","covidpriordose3","dose3_booster",
                       	all.vars(as.formula(matching.formula)))) 

  compare <- fcase(treat=="PvM","keepme_PM",
					treat=="MvP","keepme_PM") 

  elig <- fcase(treat=="PvM","eligible_PM", 
					treat== "MvP","eligible_PM") 

  print(treat)
  print(compare); print(names(dat))
  
  dat <- setDT(dat)[(get(compare)==1),]; print(paste0("subset by ", compare)) #selects only the comparison of interest by the keepme_XY variable
  
  dat[, `:=` (newperiod=seq(.N), outcome=max(get(event))), 
      by=(newid)
  ][,days:=(max(newperiod)-1), by=(newid)  ## by convention: -1 because newperiod otherwise ranges 1:179 (instead of 0:178)
  ]; print("reorganized after subset")

  # Create variables for matching on Age (Age_at_index) and subsetting by age (age)
  dat$age <- dat$Age_at_index 	#for later subgroup analysis restriction 
  dat$Age_at_index <- as.numeric(cut(dat$Age_at_index, breaks=seq(min(dat$Age_at_index)-1, max(dat$Age_at_index)+1, by=age.bin), include.lowest = TRUE)) 

  # Create numeric version of groups for matching on calendar date
  dat$Caldate_num <- as.numeric(as.Date(dat$Caldate, format="%d%b%Y"))
  dat$Caldate_num <-  as.numeric(cut(dat$Caldate_num, breaks=seq(min(dat$Caldate_num)-1, max(dat$Caldate_num)+1, by=date.bin), include.lowest = FALSE)); print("coarsening")

  # Rename variables in the dataset for the contrast of interest
  setnames(dat, c(treat, elig), c("group.binary", "eligible"))

  # Write-out the dataset with eligible individuals by treatment and event type
  fwrite(dat[eligible==1,.SD, .SDcols=keep.vars], file=paste0(prefix, "/data/cleaned-dat-",treat,"-",event,"-",variant,".csv")) 

  # For datasets on the p-drive, reduce to the variables needed to save space
  if(!server){return(dat[eligible==1,.SD, .SDcols=keep.vars])}

  # Output a log showing completed data prep
  if(server){sink(); sink(type="message"); closeAllConnections()}
}
