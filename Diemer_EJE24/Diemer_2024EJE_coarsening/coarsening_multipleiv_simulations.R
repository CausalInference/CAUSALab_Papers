###Simulations using multiple proposed instruments
##created by: Lizzie Diemer
##created on: 04/19/2022
#######

# loading key library
library(parallel)
library(tidyverse)
library(haven)
library(kableExtra)

#setting up cluster
##setting up parallelization
num_cores <- detectCores() - 1
my_cluster <- parallel::makeCluster(num_cores, setup_timeout = 0.5)

##Running everything in clusterEvalQ so that is included on each cluster machine
clusterEvalQ(my_cluster,{
  ##set random seed
  set.seed(587643)
  ## loading necessary libraries
  library(tidyverse)
  
  ##instrumental inequalities
  ##function that returns maximum value of the instrumental inequalities for a single IV
  ##Inputs:
  ##data=dataset
  ##y = character string naming outcome
  ##x= character string naming exposure
  ##instrument = character string naming proposed IV
  ##weight = optional argument naming weight to be used. If none is provided, defaults to 1.
  #########################
  run_instrumental_inequalities_singlejointiv <- function(data, y, x,
                                                          instrument, weight = NULL){
    # # Set variable names
    data$Y <- data[[y]]
    data$X <- data[[x]]
    data$IV <- data[[instrument]]
    if(is.null(weight)){data$weight <- 1}else{data$weight <- data[[weight]]}
    
    # Count the frequencies of the different values for IV
    IV_counts <- data %>% 
      group_by(IV) %>% 
      tally(., wt = weight) %>% 
      rename(n_IV=n)
    
    # Count the frequencies of the different combinations for IV, X and Y
    grouped_counts <- data %>% 
      group_by(IV, X, Y) %>% 
      tally(., wt = weight)
    
    # Merge counts and fill na with 0
    grouped_counts <- merge(grouped_counts, IV_counts, by='IV')
    grouped_counts$n_IV[is.na(grouped_counts$n_IV)] <- 0
    
    # Calculate proportions
    grouped_counts$p = grouped_counts$n / grouped_counts$n_IV
    
    # Calculate the max per X, Y (over IV), sum per X (over Y), calculate max
    # Equation [3] Pearl, J. (1995)
    
    pearl <- grouped_counts %>%
      group_by(X, Y) %>%
      summarise(max_p = max(p)) %>%
      group_by(X) %>%
      summarise(sum_maxp=sum(max_p))
    
    max_ineq <- max(pearl$sum_maxp)
    
    return(max_ineq)
    
  }
  
  
  ##Function to apply instrumental inequalities across multiple proposed instruments
  ##Inputs:
  ##data=dataset
  ##outcome = character string naming outcome
  ##exposure= character string naming exposure
  ##IV = character vector naming proposed IVs
  ## input_weight: the name of an externally calculated weight variable (character)
  ## generate_weight: options to generate unstabilized IP weight for instrument based on input covariates
  ## weight_covs: character vector naming covariates to be used in generation in IP weights
  ## count_output: options to limit additional output of function. Options are 
  ## "unweighted" (calculates numbers of contributing rows without weights), "weighted"
  ## (calculates number of individuals contributing in weighted pseudopopulation), or 
  ## "full" (gives both weighted and unweighted number of individuals in levels of joint instrument).
  #########################
  instrumental_inequalities_multiv <- function(datasetname, IV, exposure, outcome, input_weight = NULL, 
                                               generate_weights = FALSE, weight_covs = NULL, count_output= "full"){
    
    k <- length(IV)
    
    # create list of sets of instruments
    mylist <- lapply(seq_along(IV), function(i) combn(IV, i, FUN = list))
    
    mylist <- flatten(mylist)
    
    
    # create summary table
    summarydat <- matrix(nrow=length(mylist), ncol=1)
    colnames(summarydat) <- c('BP Inequalities Max Value')
    
    rownames(summarydat) <- c(mylist)
    
    #get names alone for covariates for weights generated
    covnames <- weight_covs[!(str_detect(weight_covs, "\\*|\\^"))]
    
    # create variable and run results for each possible combination of proposed IVs
    for (i in 1:length(mylist)){
      dat <- datasetname %>% select(IV, everything())  
      dat$a <- dat[[exposure]] 
      dat$y <- dat[[outcome]]
      if(is.null(input_weight)){dat$input_weight <- 1}else{dat$input_weight <- dat[[input_weight]]}
      dat <- dat%>% select(IV, a, y, input_weight, covnames)
      dat <- dat %>% drop_na()
      n_uniq_Y <- length(unique(dat$Y))
      
      #create new joint variable
      IVT = mylist[i]
      dat <- unite_(dat,"jointIV", flatten(IVT), remove = FALSE)
      
      
      #generate weights for jointIV if generated weights turned on
      if(generate_weights==TRUE && is.null(weight_covs)==FALSE && length(unique(dat$jointIV)) > 2){
        denom <- reformulate(weight_covs, response = "jointIV")
        
        denom_obj <- multinom(denom, dat, trace=FALSE)
        p_denom <- as.data.frame(predict(denom_obj, type="probs"))
        for (g in 1:length(sort(unique(dat$jointIV)))){
          dat$denom[with(dat, jointIV == sort(unique(dat$jointIV))[g])] <- p_denom[dat$jointIV==sort(unique(dat$jointIV))[g],g]
        }
        dat$gen_weight <- 1/dat$denom
      }else if(generate_weights==TRUE && is.null(weight_covs)==FALSE && length(unique(dat$jointIV))<3){
        denom <- reformulate(weight_covs, response = "as.numeric(jointIV)")
        denom_obj <- glm(denom, dat, family="binomial")
        dat$denom <- predict(denom_obj, type="response")
        dat$gen_weight <- 1/dat$denom
      }else if(generate_weights==TRUE && is.null(weight_covs)==TRUE){
        dat$gen_weight <- 1
        print("No covariates supplied, IP weights set to 1")
      } else{
        dat$gen_weight <- 1
      }
      
      #multiply input weight and generated weight together
      dat$weight <- dat$input_weight*dat$gen_weight
      
      #running instrumental inequalities function
      combo <- run_instrumental_inequalities_singlejointiv(data=dat,
                                                           y="y", x="a", 
                                                           instrument="jointIV", weight="weight")
      
      ineq <- combo
      
      #print IV inequalities held or no
      
      summarydat[i, 1] <- ineq
      
      
      
      print(paste0(i," of ",length(mylist), " combinations complete"))
    }
    
    
    #make copy of max values of BP inequalities
    max_bp_ineqs <- unname(summarydat[, 1])
    rm(summarydat)
    
    
    
    
    resultslist <- list(mylist, k, max_bp_ineqs)
    return(resultslist)
  }
  
  ##Data generating function
  ##parameters: n = number individuals in sample, nz = number of proposed instruments,
  ##z_prob = vector of probabilities of zs (for generation), b_prev_a = baseline 
  ##prevalence of exposure (beta 0 for a), b_ua = beta for effect of u on a, 
  ##b_za = vector of betas for effect of zs on exposure, b_zy = vector of betas for 
  ##effects of z on y, b_uy = effect of u on y, b_prev_y = baseline prev of outcome 
  ##(beta 0 for y).
  
  
  datagen <- function(n=1000, nz=1, z_size = rep(1, 1), z_prob=rep(0.5,1), b_prev_a0 = 2.43, b_ua=0.1, 
                      b_za0=rep(0.022,1), b_a1 =  0.58, 
                      b_za1 = rep(0.58, 1), 
                      a_amount_var, 
                      y_formula){
    ##library so that can use pipe
    library(magrittr)
    
    
    ##creating inverse logit function
    inv.logit <- function(x) {
      return(exp(x)/(1+exp(x)))
    }
    
    # Create nz instruments
    for (i in 1:nz) {assign(paste0("z", i), rbinom(n, size = z_size[i], prob = z_prob[i]))}
    
    # U, A and Y
    ##u is standard normal
    u <- rnorm(n = n, mean = 0, sd = 1)
    ##creating formula for a - p(a) = b_prev_a + b_ua*u + z1 +z2+ ..+z-nz
    a_ever_formula <- paste0("b_prev_a0 + b_ua*u + ", paste0(b_za0, "*",
                                                             paste0("z",1:nz), collapse= " + "))
    ##evaluating function to run a - a is random binomial of size 1, with probability inv.logit(formula above)
    a_ever <- rbinom(n, size = 1, prob = inv.logit(eval(parse(text = a_ever_formula))))
    ##a amount formula
    a_amount_formula <- paste0("b_a1 + b_ua*u + ", paste0(b_za1, "*", 
                                                          paste0("z", 1:nz), collapse = " + "))
    ##evaluating function to run a_amount - log-normal distribution
    a_amount <- rlnorm(n = n, meanlog = eval(parse(text = a_amount_formula)), sdlog = a_amount_var)
    ##generating true value of a
    a <- ifelse(a_ever==0, 0, a_amount)
    ##because y value varies across simulations, introduce as parameter
    ##evaluating formula to get y vector
    y <- rbinom(n, size = 1, prob = inv.logit(eval(parse(text = y_formula))))
    #binding zs together - use mget because it evaluates pasted names as objects, 
    ##not just names. do.call creates function call from function and list containing
    ##arguments
    test <- do.call(cbind, mget(paste0("z", 1:nz)))
    ##cbinding the rest of this
    test <- cbind(test, u, a, y) %>% data.frame
    
    return(list(data=test,
                params=as.list(match.call())))
  }
  
  ## function to run all other functions within lapply    
  sim_sum <- function(i){
    ds <- datagen(n=param_grid$n[i], nz=param_grid$nz[i], 
                  z_size = param_grid$z_size[i], 
                  z_prob=eval(parse(text = param_grid$z_prob[i])), 
                  b_prev_a0 = param_grid$b_prev_a0[i], b_ua=param_grid$b_ua[i], 
                  b_za0=rep(param_grid$b_za0[i],param_grid$nz[i]), 
                  b_a1 =  param_grid$b_a1[i], 
                  b_za1 = eval(parse(text = param_grid$b_za1[i])),
                  a_amount_var= param_grid$a_amount_var[i], 
                  y_formula = param_grid$y_formula[i])
    
    
    ymean <- mean(ds[[1]]$y)
    amean <- mean(ds[[1]]$a)
    zmean <- mean(ds[[1]]$z1)
    ##individuals within drinking categories
    ds[[1]]$a_cat <- ifelse(ds[[1]]$a == 0, 1, 
                            ifelse(ds[[1]]$a>0 & ds[[1]]$a<1, 2,
                                   ifelse(ds[[1]]$a>=1 & ds[[1]]$a<=2, 3, 
                                          ifelse(ds[[1]]$a >2 & ds[[1]]$a <=4,
                                                 4, 5))))
    freq <- ds[[1]] %>%
      group_by(a_cat) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n))
    freq <- freq$freq
    ##strength of IV
    ivst_lm <- lm(a ~ z1, data = ds[[1]])
    iv_st <- ivst_lm$coefficients[[2]]
    
    ##dichotomized a
    ds[[1]]$a_bin <- ifelse(ds[[1]]$a==0, 0, 1)
    
    ##deciles of a
    ds[[1]]$a_decile <- ntile(x=ds[[1]]$a, 10)
    
    ## median split
    ds[[1]]$a_median <- ntile(x=ds[[1]]$a, 2)
    ##instrumental inequalities max value
    ##binary
    iv <- sapply(1:param_grid$nz[i], function(x){paste0("z",x)})
    test <- instrumental_inequalities_multiv(datasetname = ds[[1]], IV = iv,
                                             exposure = "a_bin", outcome = "y")
    checklength <- rapply(test[1], length, how = "list")
    ineqlengths <-as.data.frame(cbind(test[[3]], unlist(checklength)))
    colnames(ineqlengths) <- c("ineq", "niv")
    maxviol <- ineqlengths %>% group_by(niv) %>% summarise(max = max(ineq))
    
    ##add bit where 1 if violated, 0 if not
    maxviol$ineq_viol <- ifelse(maxviol$max > 1, 1, 0)
    
    ##now need way to transform this into arbitrarily long vector named for # of instru
    
    testvec <- t(maxviol$ineq_viol)
    niv_names <- sapply(maxviol$niv, function(x){paste0("ineqviol_niv",x)})
    colnames(testvec) <- niv_names
    
    ##same for categorical a
    test <- instrumental_inequalities_multiv(datasetname = ds[[1]], IV = iv,
                                             exposure = "a_cat", outcome = "y")
    checklength <- rapply(test[1], length, how = "list")
    ineqlengths <-as.data.frame(cbind(test[[3]], unlist(checklength)))
    colnames(ineqlengths) <- c("ineq", "niv")
    cat_viol <- ineqlengths %>% group_by(niv) %>% summarise(max = max(ineq))
    
    ##add bit where 1 if violated, 0 if not
    cat_viol$ineq_viol <- ifelse(cat_viol$max > 1, 1, 0)
    
    ##now need way to transform this into arbitrarily long vector named for # of instru
    
    cat_vec <- t(cat_viol$ineq_viol)
    niv_names <- sapply(cat_viol$niv, function(x){paste0("ineqviol_niv",x)})
    colnames(cat_vec) <- niv_names
    
    ##and same but with deciles
    test <- instrumental_inequalities_multiv(datasetname = ds[[1]], IV = iv,
                                             exposure = "a_decile", outcome = "y")
    checklength <- rapply(test[1], length, how = "list")
    ineqlengths <-as.data.frame(cbind(test[[3]], unlist(checklength)))
    colnames(ineqlengths) <- c("ineq", "niv")
    maxviol <- ineqlengths %>% group_by(niv) %>% summarise(max = max(ineq))
    
    ##add bit where 1 if violated, 0 if not
    maxviol$ineq_viol <- ifelse(maxviol$max > 1, 1, 0)
    
    ##now need way to transform this into arbitrarily long vector named for # of instru
    
    dec_vec <- t(maxviol$ineq_viol)
    niv_names <- sapply(maxviol$niv, function(x){paste0("ineqviol_niv",x)})
    colnames(dec_vec) <- niv_names
    
    ##and same but with median split
    test <- instrumental_inequalities_multiv(datasetname = ds[[1]], IV = iv,
                                             exposure = "a_median", outcome = "y")
    checklength <- rapply(test[1], length, how = "list")
    ineqlengths <-as.data.frame(cbind(test[[3]], unlist(checklength)))
    colnames(ineqlengths) <- c("ineq", "niv")
    maxviol <- ineqlengths %>% group_by(niv) %>% summarise(max = max(ineq))
    
    ##add bit where 1 if violated, 0 if not
    maxviol$ineq_viol <- ifelse(maxviol$max > 1, 1, 0)
    
    ##now need way to transform this into arbitrarily long vector named for # of instru
    
    med_vec <- t(maxviol$ineq_viol)
    niv_names <- sapply(maxviol$niv, function(x){paste0("ineqviol_niv",x)})
    colnames(med_vec) <- niv_names
    
    return(c(amean, ymean, testvec, cat_vec, dec_vec, med_vec))
  }
})
# setting parameters that vary over simulations
beta2 <- c(seq(from = 1, to = 1.2, by = 0.1), seq(from=1.5, to = 7.9, by = 0.5))
yformula <- sapply(beta2, FUN = function(x){paste0("-1.56 + log(",x,")*a + 0.1*u")})
yformula2 <- sapply(beta2, FUN = function(x){paste0("-1.56 + 0.02765*a + log(", x, ")*(a*a) + 0.1*u")})
yformula <- c(yformula, yformula2)
param_grid <- expand.grid(n = c(1000, 10000, 100000),  nz = 4, z_size = c(1,2),
                          b_za1 = c("c(0.0108, 0.02089, 0.02595, 0.01372)", "rep(0.005, 4)"),
                          y_formula = yformula,
                          stringsAsFactors = FALSE)

#setting parameters fixed over simulations
param_grid$b_prev_a0 <- 2.43
param_grid$b_ua <- 0.1
param_grid$b_za0 <- 0.022
param_grid$b_a1 <- 0.58
param_grid$a_amount_var <- 0.451
param_grid$z_prob <- "c(0.718, 0.601, 0.597, 0.508)"

##number replications
n_reps = 500
clusterExport(my_cluster, var = c("n_reps", "param_grid"))

sim_res <-parLapply(my_cluster, seq_len(nrow(param_grid)), function(i){
  sim <- replicate(n_reps,sim_sum(i))
  simsum_ymean <- rowSums(sim)[2]/n_reps
  simsum_amean <- rowSums(sim)[1]/n_reps
  
  
  ##assigning values for binary ineq
  numbers <- 1:param_grid$nz[i]
  for (j in numbers){
    assign(paste0("binviol_niv",j), (rowSums(sim)[2+j])/n_reps)
  }
  
  ##assigning values for categorical ineq
  for (j in numbers){
    assign(paste0("catviol_niv",j), (rowSums(sim)[2+j + max(numbers)])/n_reps)
  }
  
  ##assigning values for decile ineq
  for (j in numbers){
    assign(paste0("decviol_niv",j), (rowSums(sim)[2+j + 2*max(numbers)])/n_reps)
  }
  
  ##assigning values for decile ineq
  for (j in numbers){
    assign(paste0("medviol_niv",j), (rowSums(sim)[2+j + 3*max(numbers)])/n_reps)
  }
  
  
  return(c(yprop_mean = simsum_ymean, amean_total = simsum_amean, 
           binviol_niv1 = binviol_niv1, binviol_niv2 = binviol_niv2, 
           binviol_niv3 = binviol_niv3, binviol_niv4 = binviol_niv4,
           catviol_niv1 = catviol_niv1, catviol_niv2 = catviol_niv2,
           catviol_niv3 = catviol_niv3, catviol_niv4 = catviol_niv4,
           decviol_niv1 = decviol_niv1, decviol_niv2 = decviol_niv2,
           decviol_niv3 = decviol_niv3, decviol_niv4 = decviol_niv4,
           medviol_niv1 = medviol_niv1, medviol_niv2 = medviol_niv2,
           medviol_niv3 = medviol_niv3, medviol_niv4 = medviol_niv4
  ))}) %>% do.call(rbind, .) %>%
  as.data.frame()

sim_res_complete <- cbind(param_grid, sim_res)

stopCluster(my_cluster)

save(sim_res_complete, file="./results/multiv_results.Rda")