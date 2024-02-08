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

# preamble
preamble <- function(server=TRUE){
  if(server){
    prefix <<- "[SAS_folder_ORD_Project]/Vax_Standard/analysis" #location of files on SAS Grid
  } else {
    prefix <<- "[P_drive_ORD_Project]/Vax_Standard/analysis"  
  } 

  server<<- server
  lapply(c("survival","survminer","tidyverse", "data.table", "MatchIt"), function(x) suppressPackageStartupMessages(library(x, character.only=TRUE, quietly=TRUE)))
  #print("loaded packages")
}
