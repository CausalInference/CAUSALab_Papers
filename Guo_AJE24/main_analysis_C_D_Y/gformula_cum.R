################################################################
# /n/hpnh/Users/fyguo/R_proj/code_Review/main_analysis_C_D_Y/gformula_tot_1115.R
# Programer: Fuyu Guo
# Last updated: 2023-Jun-24
# This code is aimed to generate a function to calculate the cumulative sum of histvars
# I modified the code from cumavg.
cum <- function (pool, histvars, time_name, t, id_name, below_zero_indicator = TRUE) 
{
  if (t == 0 & !below_zero_indicator) {
    lapply(histvars, FUN = function(histvar) {
      pool[get(time_name) == t, `:=`((paste("cum_", 
                                            histvar, sep = "")), 
                                     as.double(pool[get(time_name) == t][[histvar]]))]
    })
  }
  else {
    current_ids <- unique(pool[get(time_name) == t][[id_name]])
    colnam <- colnames(pool)
    if (!(paste("cum_", "_", histvars[1], sep = "") %in% 
          colnam)) {
      id_factor <- is.factor(pool[[id_name]])
      if (id_factor) {
        lapply(histvars, FUN = function(histvar) {
          pool[get(time_name) == t, `:=`((paste("cum_", 
                                                histvar, sep = "")), 
                                         as.double(tapply(pool[get(id_name) %in% current_ids & get(time_name) <= t][[histvar]], 
                                                          droplevels(pool[get(id_name) %in% current_ids & get(time_name) <= t][[id_name]]), 
                                                          FUN = sum)))]})
      }
      else {
        lapply(histvars, FUN = function(histvar) {
          pool[get(time_name) == t, `:=`((paste("cum_", 
                                                histvar, sep = "")), as.double(tapply(pool[get(id_name) %in% 
                                                                                             current_ids & get(time_name) <= t][[histvar]], 
                                                                                      pool[get(id_name) %in% current_ids & get(time_name) <= 
                                                                                             t][[id_name]], FUN = sum)))]
        })
      }
    }
    else {
      for (histvar in histvars) {
        pool[get(time_name) == t, `:=`((paste("cum_", 
                                              histvar, sep = "")), as.double((pool[get(id_name) %in% 
                                                                                     current_ids & get(time_name) == (t - 1)][[paste("cum_", 
                                                                                                                                     histvar, sep = "")]]  + pool[get(id_name) %in% 
                                                                                                                                                                                 current_ids & get(time_name) == t][[histvar]])))]
      }
    }
  }
}

#################################################
# this function is to add a function like tsswitch in SAS Marcro
# when the histvars switch from 0 to 1
# then tsint_histvars = cumulative sum of histvars * histvars
# cumulative sum of histvars is the time since histvars switched from 0 to 1


tsswitch <- function(pool, histvars, time_name, t, id_name) {
  for (histvar in histvars){
    pool[get(time_name) == t,
         `:=`(paste("ts_int", histvar, sep = "_"),
              as.double(pool[get(time_name) == t][[paste("cum", histvar, sep = "_")]]*pool[get(time_name) == t][[histvar]] - 1))]
    
    # this line is to make sure that the minimum ts_int is 0.
    pool[get(paste("ts_int", histvar, sep = "_")) < 0,
         `:=`(paste("ts_int", histvar, sep = "_"),0)]
  }
}



