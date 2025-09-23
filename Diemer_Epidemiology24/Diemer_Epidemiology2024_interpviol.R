##########################################################################
##Simulations to evaluate the effects of the instrumental inequalities for
##varying combinations of proposed joint instruments
##Author: Elizabeth Diemer
##Published in: DOI: 10.1097/EDE.0000000000001751
##########################################################################
##Load necessary packages
if(!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)
if(!require(data.table)) install.packages("tidyverse")
library(data.table)
if(!require(gridExtra)) install.packages("gridExtra")
library(gridExtra)
# Function Definitions
#Data Generating Function
##Data generating function
##parameters: n = number individuals in sample, nz = number of proposed
##instruments, z_prob = vector of probabilities of zs (for generation),
##b_prev_a = baseline prevalence of exposure (beta 0 for a),
##b_ua = beta for effect of u on a, b_za = vector of betas for effect of zs on
##exposure, b_zy = vector of betas for effects of z on y, b_uy = effect of u on
##y, b_prev_y = baseline prev of outcome
##(beta 0 for y).
datagen <- function(n=1000, nz=1, z_prob=rep(0.5,1), b_prev_a0 = 2.43,
                    b_ua=0.1, b_za0=rep(0.022,1), y_formula){
  ##library so that can use pipe
  library(magrittr)
  ##creating inverse logit function
  inv.logit <- function(x) {
    return(exp(x)/(1+exp(x)))
  }
  # Create nz instruments
  for (i in 1:nz) {assign(paste0("z", i), rbinom(n, size = 1,
                                                 prob = z_prob[i]))}
  # U, A and Y
  ##u is standard normal
  u <- rnorm(n = n, mean = 0, sd = 1)
  ##creating formula for a - p(a) = b_prev_a + b_ua*u + z1 +z2+ ..+z-nz
  a_ever_formula <- paste0("b_prev_a0 + b_ua*u + ", paste0(b_za0, "*",
                                                           paste0("z",1:nz),
                                                           collapse= " + "))
  ##evaluating function to run a - a is random binomial of size 1, with
  #probability inv.logit(formula above)
  a <- rbinom(n, size = 1, prob = inv.logit(eval(parse(text = a_ever_formula)
  )))
  ##because y value varies across simulations, introduce as parameter
  ##evaluating formula to get y vector
  y <- rbinom(n, size = 1, prob = inv.logit(eval(parse(text = y_formula))))
  #binding zs together - use mget because it evaluates pasted names as objects,
  ##not just names. do.call creates function call from function and list
  ##containing arguments
  test <- do.call(cbind, mget(paste0("z", 1:nz)))
  test <- cbind(test, u, a, y) %>% data.frame
  return(list(data=test,
              params=as.list(match.call())))
}
## internal function to get maximum value of the instrumental inequalities for
##single given joint instrument
## NOTE: ALL VARIABLES CAN BE MULTICATEGORICAL BUT CANNOT BE CONTINUOUS
## ARGUMENTS
## data: names of dataset (data.frame)
## instrument: name of instrument variable (character)
## x: name of exposure variable (character)
## y: name of outcome variable (character)
run_instrumental_inequalities_singlejointiv <- function(data, y, x,
                                                        instrument,
                                                        weight = NULL){
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
    summarise(max_p = max(p),.groups="keep") %>%
    group_by(X) %>%
    summarise(sum_maxp=sum(max_p))
  max_ineq <- max(pearl$sum_maxp)
  return(max_ineq)
}
## Function applying the instrumental inequalities for a given exposure-outcome
## pair across all combinations of multiple proposed instruments.
##
## ARGUMENTS:
## datasetname: the dataset to be used (data.frame)
## IV: a character vector containing the names of the variables proposed as
##instruments (character vector)
## exposure: the name of the exposure of interest (character)
## outcome: the name of the outcome of interest (character)
## single_weight: the name of an externally calculated weight variable
## (character)
## count_output: options to limit additional output of function. Options are
## "unweighted" (calculates numbers of contributing rows without weights),
## "weighted" (calculates number of individuals contributing in weighted
## pseudopopulation), or "full" (gives both weighted and unweighted number of
## individuals in levels of joint instrument).
##
## The function outputs a list of 4 results - the first is a summary table of the
## findings, and the second, third, and fourth are information necessary for the
## creation of the bfi graph visualizations of the results.
## The function will remove any rows with missing data.
instrumental_inequalities_multiv <- function(datasetname, IV, exposure, outcome,
                                             input_weight = NULL,
                                             generate_weights = FALSE,
                                             weight_covs = NULL,
                                             count_output= "full"){
  k <- length(IV)
  # create list of sets of instruments
  mylist <- lapply(seq_along(IV), function(i) combn(IV, i, FUN = list))
  mylist <- flatten(mylist)
  # check GRS and add to list
  datasetname$AlleleScore <- apply(datasetname[IV], 1, sum)
  mylist <- c(mylist, "AlleleScore")
  # create summary table
  summarydat <- matrix(nrow=length(mylist), ncol=2)
  colnames(summarydat) <- c('max_value_bp_ineqs', 'bp_ineqs_hold')
  rownames(summarydat) <- c(mylist)
  #get names alone for covariates for weights generated
  covnames <- weight_covs[!(str_detect(weight_covs, "\\*|\\^"))]
  # create variable and run results for each possible combination of proposed IVs
  for (i in 1:length(mylist)){
    dat <- datasetname %>% select(all_of(IV), everything())
    dat$a <- dat[[exposure]]
    dat$y <- dat[[outcome]]
    if(is.null(input_weight)){dat$input_weight <- 1
    }else{dat$input_weight <- dat[[input_weight]]}
    dat <- dat%>% select(all_of(IV), AlleleScore, a, y, input_weight,
                         all_of(covnames))
    dat <- dat %>% drop_na()
    n_uniq_Y <- length(unique(dat$Y))
    #create new joint variable
    IVT = mylist[i]
    #dat <- unite_(dat,"jointIV", flatten(IVT), remove = FALSE)
    dat <- unite(dat,"jointIV", unlist(IVT), remove = FALSE)
    #generate weights for jointIV if generated weights turned on
    if(generate_weights==TRUE && is.null(weight_covs)==FALSE &&
       length(unique(dat$jointIV)) > 2){
      denom <- reformulate(weight_covs, response = "jointIV")
      denom_obj <- multinom(denom, dat, trace=FALSE)
      p_denom <- as.data.frame(predict(denom_obj, type="probs"))
      for (g in 1:length(sort(unique(dat$jointIV)))){
        dat$denom[with(dat, jointIV == sort(unique(dat$jointIV))[g])] <-
          p_denom[dat$jointIV==sort(unique(dat$jointIV))[g],g]
      }
      dat$gen_weight <- 1/dat$denom
    }else if(generate_weights==TRUE && is.null(weight_covs)==FALSE &&
             length(unique(dat$jointIV))<3){
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
                                                         instrument="jointIV",
                                                         weight="weight")
    ineq <- combo
    #print IV inequalities held or no
    summarydat[i, 2] <- if (ineq<=1) {1} else {0}
    summarydat[i, 1] <- ineq
  }
  return(summarydat)
}
##Function to create instrumental inequality plots
##
##ARGUMENTS:
##instru: list containing names of all possible combinations of variables -
##resultslist[[2]] from instrumental_inequalities_multiv function (list)
##k: number of variable jointly proposed as instruments - resultslist[[3]]
##from instrumental_inequalities_multiv function (integer)
##ineqs: vector of maximum value of the instrumental inequalities for
##each combination of variables - resultslist[[4]] from
##instrumental_inequalities_multiv function (vector)
##title: optional title of plot (character)
##
##Required arguments are supplied by instrumental_inequalities_multiv function
##as second, third, and fourth objects on output list of results. The required
##inputs from the instrumental_inequalities_multiv function should be double
##bracketed.
plot_instrumental_inequalities <- function(instru,k,ineqs, title){
  ##determining title
  if(missing(title)){title<-NA}
  ##generating nodes dataset
  nodes <- data.frame(id=1:length(flatten(instru)))
  nodes$label1 <- lapply(1:length(flatten(instru)), function(i){
    flatten(instru)[[i]]})
  ##generate y position - aligning along y spots
  nodes$y <- 1
  for (j in 1:k){nodes<-nodes %>% mutate(y=ifelse(label1==instru[[j]], j, y))}
  nodes$y <- as.numeric(nodes$y)
  nodes$y <- nodes$y*10
  ##label of group - creates their x axis coordinates
  nodes$x <- flatten(lapply(1:length(instru), function(i){rep(
    length(instru[1:i]), length(flatten(instru[i])))}))
  nodes$x <- unlist(nodes$x)
  nodes$x <- nodes$x*10
  ##edges generation
  edges <- data.frame(id=1:length(unique(nodes$y)))
  edges$fromy <- unique(nodes$y)
  edges$fromx <- tapply(nodes$x, nodes$y, min)
  edges$toy <- edges$fromy
  edges$tox <- tapply(nodes$x, nodes$y, max)
  vertedges <- data.frame(id=1:length(unique(nodes$x)))
  vertedges$fromx <- unique(nodes$x)
  vertedges$tox <- vertedges$fromx
  vertedges$fromy <- tapply(nodes$y, nodes$x, min)
  vertedges$toy <- tapply(nodes$y, nodes$x, max)
  edges <- rbind(edges, vertedges)
  ## size and color of nodes
  printsumnum <- function(i){print(ineqs[i])}
  nodes <- nodes %>% rowwise %>% mutate(colorfactor=printsumnum(x/10))
  nodes$ii <- ifelse(nodes$colorfactor<=1, NA, cut(as.numeric(nodes$colorfactor),
                                                   breaks=seq(1, 2, len=100),
                                                   include.lowest=TRUE))
  nodes$color <- ifelse(nodes$colorfactor<=1, colorRampPalette("grey99")(1),
                        colorRampPalette(c("gray60", "gray22"))(99)[nodes$ii])
  ##node labels
  nodes$label1 <- ifelse(nodes$id<=k, nodes$label1,
                         ifelse(nodes$id==length(flatten(instru)),
                                nodes$label1, ""))
  ##generating plot itself
  layout(matrix(1:2,nrow=1),widths=c(0.8,0.2))
  par(mar=c(5.1,4.1,4.1,1.0))
  plot(c(-2,max(nodes$x)+5),c(0,max(nodes$y)+5),type = 'n',
       axes = F,xlab = '', ylab = '', main = title)
  segments(edges$fromx, edges$fromy, x1= edges$tox, y1= edges$toy)
  points(nodes$x, nodes$y, pch=21, cex=3.5, bg=nodes$color)
  text(nodes$x[1:length(flatten(instru))-1],
       nodes$y[1:length(flatten(instru))-1],
       nodes$label1[1:length(flatten(instru))-1], pos=1, offset=1)
  text(nodes$x[length(flatten(instru))]+2, nodes$y[length(flatten(instru))],
       nodes$label1[length(flatten(instru))], pos=3, offset=1.5)
  legend_image <- as.raster(matrix(
    colorRampPalette(c('gray22', 'gray60'))(99), ncol=1))
  par(mar=c(5.1,1.0,4.1,2.1))
  plot(c(0,2),c(0,2),type = 'n', axes = F,xlab = '', ylab = '', main = 'Legend')
  text(x=1.5, y = c(.5, seq(1,2,by=.25)), labels = c(0, 1, seq(1.25,2,by=.25)))
  rasterImage(legend_image, 0, 1, 1,2)
  rect(0.025,.5,1,1, col='grey99', border='black')
}
#set seed for simulations
set.seed(22026)
#Generate Study 1
y_formula <-"-1.56 + log(3)*a +0.1*u + 1.5*a*u + log(7.2)*z3"
study1 <- datagen(n=1000000, nz=3, z_prob=c(0.8,0.8, 0.5), b_prev_a0 = -1.63,
                  b_ua=0.1, b_za0=c(log(1.1),log(1.05),log(3.96)),
                  y_formula=y_formula)
# Generate Study 2
y_formula <-"-1.56 + log(3)*a +0.1*u + -1.5*a*u + log(2.7)*z1 + log(2.7)*z2"
study2 <- datagen(n=1000000, nz=3, z_prob=c(0.8,0.8, 0.5), b_prev_a0 = -1.63,
                  b_ua=0.1, b_za0=c(log(1.1),log(1.05),log(3.96)),
                  y_formula=y_formula)
# Apply the instrumental inequalities to Study 1
study1_ineqs <- instrumental_inequalities_multiv(datasetname=study1[[1]],
                                                 IV=c("z1","z2","z3"),
                                                 exposure="a", outcome="y", input_weight = NULL,
                                                 generate_weights = FALSE, weight_covs = NULL,
                                                 count_output= "full")
# plot inequalities
# create list of sets of instruments
mylist <- lapply(seq_along(c("z1","z2","z3")),
                 function(i) combn(c("z1","z2","z3"), i, FUN = list))
instru <- flatten(mylist)
study1_fig <- plot_instrumental_inequalities(instru=instru, k=3,
                                             ineqs = study1_ineqs,
                                             title= "Study 1")
# Apply the instrumental inequalities to Study 2
study2_ineqs <- instrumental_inequalities_multiv(datasetname=study2[[1]],
                                                 IV=c("z1","z2","z3"),
                                                 exposure="a", outcome="y", input_weight = NULL,
                                                 generate_weights = FALSE, weight_covs = NULL,
                                                 count_output= "full")
study2_fig <- plot_instrumental_inequalities(instru=instru, k=3,
                                             ineqs = study2_ineqs,
                                             title= "Study 2")
##Add both plots to single figure
study1_fig <- plot_instrumental_inequalities(instru=instru, k=3,
                                             ineqs = study1_ineqs,
                                             title= "Study 1")
study2_fig <- plot_instrumental_inequalities(instru=instru, k=3,
                                             ineqs = study2_ineqs,
                                             title= "Study 2")
