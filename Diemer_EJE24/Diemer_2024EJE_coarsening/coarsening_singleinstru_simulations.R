##Dataset creation for coarsening simulations (single instrument case)
##created by: Lizzie Diemer
##created 02/08/2022
##############################


# loading key library
library(parallel)
library(tidyverse)

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
  
  ##instrumental inequalities, attributable to Kelly Guo
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
  
  
  ##Data generating function
  ##parameters: n = number individuals in sample, nz = number of proposed instruments,
  ##z_prob = vector of probabilities of zs (for generation), b_prev_a = baseline 
  ##prevalence of exposure (beta 0 for a), b_ua = beta for effect of u on a, 
  ##b_za = vector of betas for effect of zs on exposure, b_zy = vector of betas for 
  ##effects of z on y, b_uy = effect of u on y, b_prev_y = baseline prev of outcome 
  ##(beta 0 for y).
  
  
  datagen <- function(n=1000, nz=1, z_prob=rep(0.5,1), b_prev_a0 = 2.43, b_ua=0.1, 
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
    for (i in 1:nz) {assign(paste0("z", i), rbinom(n, size = 1, prob = z_prob[i]))}
    
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
  
  
  
  ##combined function to check means & prevs
  sim_sum <- function(i){
    ds <- datagen(n=param_grid$n[i], nz=param_grid$nz[i], z_prob=rep(param_grid$z_prob[i],1), 
                  b_prev_a0 = param_grid$b_prev_a0[i], b_ua = param_grid$b_ua[i], 
                  b_za0=rep(param_grid$b_za0[i],1), b_a1 =  param_grid$b_a1[i], 
                  b_za1 = rep(param_grid$b_za1[i], 1), 
                  a_amount_var = param_grid$a_amount_var[i], 
                  y_formula = param_grid$y_formula[i])
    ##proportions of binary vars
    ymean <- mean(ds[[1]]$y)
    amean <- mean(ds[[1]]$a)
    zmean <- mean(ds[[1]]$z1)
    ##individuals within drinking categories
    ds[[1]]$a_cat <- ifelse(ds[[1]]$a == 0, 1, 
                            ifelse(ds[[1]]$a>0 & ds[[1]]$a<1, 2,
                                   ifelse(ds[[1]]$a>=1 & ds[[1]]$a<=2, 3, 
                                          ifelse(ds[[1]]$a >2 & ds[[1]]$a <=4, 4, 5))))
    
    ##strength of IV
    ivst_lm <- lm(a ~ z1, data = ds[[1]])
    iv_st <- ivst_lm$coefficients[[2]]
    
    ##dichotomized a
    ds[[1]]$a_bin <- ifelse(ds[[1]]$a==0, 0, 1)
    
    ##deciles of a
    ds[[1]]$a_decile <- ntile(x=ds[[1]]$a, 10)
    
    
    ##instrumental inequalities max value
    ##binary
    ineq_bin <- run_instrumental_inequalities_singlejointiv(data=ds[[1]], y = "y", 
                                                            x = "a_bin", 
                                                            instrument = "z1")
    
    ##categories
    ineq_cat <- run_instrumental_inequalities_singlejointiv(data = ds[[1]], y = "y", 
                                                            x = "a_cat", 
                                                            instrument = "z1")
    
    ## deciles
    ineq_dec <- run_instrumental_inequalities_singlejointiv(data=ds[[1]], y="y", 
                                                            x = "a_decile", 
                                                            instrument = "z1")
    
    
    return(c(amean, ymean, ineq_bin, ineq_cat, ineq_dec))
    
  }
})

##setting up varied parameter table
beta2 <- c(seq(from = 1, to = 1.2, by = 0.1), seq(from=1.5, to = 7.9, by = 0.5))
yformula <- sapply(beta2, FUN = function(x){paste0("-1.56 + log(",x,")*a + 0.1*u")})
yformula2 <- sapply(beta2, FUN = function(x){paste0("-1.56 + 0.02765*a + log(", x, ")*(a*a) + 0.1*u")})
yformula <- c(yformula, yformula2)
param_grid <- expand.grid(n = c(1000),  b_za0 = c(0.022, 0.002), 
                          y_formula = yformula, stringsAsFactors = FALSE)
param_grid$b_za1 <- ifelse(param_grid$b_za0==0.022, .89, .005)

#setting parameters that are fixed across simulations
param_grid$nz <- 1
param_grid$z_prob <- 0.5
param_grid$b_prev_a0 <- 2.43
param_grid$b_ua <- 0.1
param_grid$b_a1 <- 0.58
param_grid$a_amount_var <- 0.451

##number replications
n_reps = 500
clusterExport(my_cluster, var = c("n_reps", "param_grid"))

##running simulation
sim_res <-parLapply(my_cluster, seq_len(nrow(param_grid)), function(i){
  sim <- replicate(n_reps,sim_sum(i))
  simsum_ymean <- rowSums(sim)[1]/n_reps
  simsum_amean <- rowSums(sim)[2]/n_reps
  ineq_bin_full <- sim[3,]
  ineq_cat_full <- sim[4,]
  ineq_dec_full <- sim[5,]
  
  return(c(yprop_mean = simsum_ymean, amean_total = simsum_amean,
           ineq_bin_full, ineq_cat_full,ineq_dec_full
           
  ))}) %>% do.call(rbind, .) %>%
  as.data.frame()

sim_res_complete <- cbind(param_grid, sim_res)

stopCluster(my_cluster)

save(sim_res_complete, file="./results/results_4ivs.Rda")  